use malachite::num::arithmetic::traits::{Abs, Lcm};
use malachite::num::basic::traits::{One, Zero};
use malachite::{Natural, Rational};
use mendeleev::{ALL_ELEMENTS, Element};
use std::cmp::{max, min, Ordering};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::iter::zip;
use std::mem;


#[test]
fn main2() {
    //let mut equ = Equation::from_latex(r"[Cr(N_{2}H_{4}CO)_{6}]_{4}[Cr(CN)_{6}]_{3}+KMnO_{4}+H_{2}SO_{4} \longrightarrow K_{2}Cr_{2}O_{7}+MnSO_{4}+CO_{2}+KNO_{3}+K_{2}SO_{4}+H_{2}O").unwrap();
    let mut equ = Equation::from_latex(r"K \longrightarrow K^+").unwrap();
    equ.solve().unwrap();
    println!("{}", equ.solution_str().unwrap());
}



/// Errors that can occur during solving system of linear equations
#[derive(Copy, Clone, Debug, Eq, Hash, PartialEq)]
pub enum Errors {
    /// Entered equation is invalid
    InvalidEquation,
    /// Solution was calculated, but is invalid
    InvalidSolution,
    /// Entered element is invalid
    InvalidElement,

    /// Matrix has wrong dimensions (rows and columns)
    WrongMatrixDimensions,
    /// There is no solution to the system of linear equations
    NoSystemSolution,
}
impl Display for Errors {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Errors::WrongMatrixDimensions => write!(f, "Wrong matrix dimensions"),
            Errors::NoSystemSolution => write!(f, "No solution"),
            Errors::InvalidElement => write!(f, "Invalid element"),
            Errors::InvalidEquation => write!(f, "Invalid equation"),
            Errors::InvalidSolution => write!(f, "Invalid solution"),
        }
    }
}
impl Error for Errors {}





/// A struct that represents a chemical equation (e.g. 2H2 + O2 -> 2H2O)
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Equation {
    /// String from which the equation was parsed
    original_str: String,
    /// A vector of reactants
    reactants: Vec<Compound>,
    /// A vector of products
    products: Vec<Compound>,
    /// A vector of solutions for reactants (stoichiometric coefficients)
    solutions_reactants: Option<Vec<i64>>,
    /// A vector of solutions for products (stoichiometric coefficients)
    solutions_products: Option<Vec<i64>>,
}
impl Equation {
    pub fn from_latex(input: &str) -> Result<Self, Errors> {
        let original_str = input;
        if input.split(r"\longrightarrow").count() != 2 { return Err(Errors::InvalidEquation); }
        let reactants_str = input.split(r"\longrightarrow").next().unwrap();
        let products_str = input.split(r"\longrightarrow").last().unwrap();

        let sanitize_equation = |string: &str| -> String {

            // remove any unused character
            let string = string.chars().filter(|&c| {
                c.is_ascii_lowercase() ||
                c.is_ascii_uppercase() ||
                c.is_ascii_digit() ||
                ['_', '^', '{', '}', '[', ']', '(', ')', '+', '-'].contains(&c)
            }).collect::<String>();

            // remove invalid digits in a string
            let mut digits_locs = HashSet::new();
            let mut inside_parentheses = false;
            let mut prev_char = ' ';
            for (i, c) in string.chars().enumerate() {
                if c == '{' {
                    inside_parentheses = true;
                } else if c == '}' {
                    inside_parentheses = false;
                } else if c.is_ascii_digit() && !inside_parentheses && !['^', '_'].contains(&prev_char) {
                    digits_locs.insert(i);
                }
                prev_char = c;
            }
            let string = string.chars().enumerate().filter(|(i, _)| !digits_locs.contains(i)).map(|(_, c)| c).collect::<String>();

            string
        };

        let process_side = |string: &str| -> Result<Vec<Compound>, Errors> {
            let mut compounds = Vec::new();

            let mut plus_locs = Vec::new();
            let mut prev_char = ' ';
            let mut inside_parentheses = false;
            for (i, c) in string.chars().enumerate() {
                if c == '{' {
                    inside_parentheses = true;
                } else if c == '}' {
                    inside_parentheses = false;
                } else if c == '+' && !inside_parentheses && prev_char != '^' {
                    plus_locs.push(i);
                }
                prev_char = c;
            }
            plus_locs.push(string.len());

            for plus_loc in 0..plus_locs.len() {
                let compound_str = string.chars().enumerate().filter(|(i, _)| {
                    if plus_loc == 0 {
                        *i < plus_locs[plus_loc]
                    } else {
                        *i > plus_locs[plus_loc - 1] && *i < plus_locs[plus_loc]
                    }
                }).map(|(_, c)| c).collect::<String>();

                if !compound_str.starts_with('e') {  // if it does start with e, it is an electron, we will add it later after calculations
                    compounds.push(Compound::from_latex(&compound_str)?);
                }
            }

            Ok(compounds)
        };

        let reactants = process_side(&sanitize_equation(reactants_str))?;
        let products = process_side(&sanitize_equation(products_str))?;

        Ok(Self {
            original_str: String::from(original_str),
            reactants,
            products,
            solutions_reactants: None,
            solutions_products: None,
        })
    }

    /// Solves the equation
    /// # Returns
    /// * `Ok` - if the equation was solved successfully
    /// * `Err` - if the equation was not solved successfully
    pub fn solve(&mut self) -> Result<(), Errors> {
        // get all the elements used in the equation (each element corresponds to one equation) (number of rows)
        let mut elements = HashSet::new();
        for compound in self.reactants.iter().chain(self.products.iter()) {
            for (element, _) in compound.elements.iter() {
                elements.insert(element);
            }
        }

        // number of coefficients in the matrix (number of columns)
        let coeff_count = self.reactants.len() + self.products.len() + 1;

        // construct matrix of coefficients
        let mut matrix = vec![vec![0; coeff_count]; elements.len() + 1];

        // set last equation in the matrix (last stoichiometric coefficient to 1, even if it isn't, it will be fixed later)
        matrix[elements.len()][coeff_count - 1] = 1;
        matrix[elements.len()][coeff_count - 2] = 1;

        // set other equations in the matrix
        for (i, element) in elements.iter().enumerate() {
            let mut col = 0;
            for compound in self.reactants.iter() {
                for (e, q) in compound.elements.iter() {
                    if e == *element {
                        matrix[i][col] = *q;
                        break;
                    }
                }
                col += 1;
            }
            for compound in self.products.iter() {
                for (e, q) in compound.elements.iter() {
                    if e == *element {
                        matrix[i][col] = -*q;
                        break;
                    }
                }
                col += 1;
            }
        }

        let solutions = solve_equations(&matrix)?;
        if solutions.iter().any(|x| *x <= 0) { return Err(Errors::InvalidSolution); }

        let (reactants_solutions, products_solutions) = solutions.split_at(self.reactants.len());
        let mut reactants_solutions = reactants_solutions.to_vec();
        let mut products_solutions = products_solutions.to_vec();

        // check if solutions are correct
        let mut reactants_element_counts = HashMap::new();
        let mut reactants_electrons = 0;
        for (i, coeff) in reactants_solutions.iter().enumerate() {
            for (element, q) in self.reactants[i].elements.iter() {
                *reactants_element_counts.entry(element).or_insert(0) += q * coeff;
            }
            reactants_electrons += self.reactants[i].electrons() * coeff;
        }
        let mut products_element_counts = HashMap::new();
        let mut products_electrons = 0;
        for (i, coeff) in products_solutions.iter().enumerate() {
            for (element, q) in self.products[i].elements.iter() {
                *products_element_counts.entry(element).or_insert(0) += q * coeff;
            }
            products_electrons += self.products[i].electrons() * coeff;
        }
        if reactants_element_counts != products_element_counts { return Err(Errors::InvalidSolution); }

        // check if electrons are balanced
        match reactants_electrons.cmp(&products_electrons) {
            Ordering::Less => { reactants_solutions.push(products_electrons - reactants_electrons) },
            Ordering::Greater => { products_solutions.push(reactants_electrons - products_electrons) },
            Ordering::Equal => {},
        }

        self.solutions_reactants = Some(reactants_solutions);
        self.solutions_products = Some(products_solutions);

        Ok(())
    }

    /// Returns the original string from which the equation was parsed
    /// # Returns
    /// * `&str` - original string
    /// # Example
    /// ```
    /// use stoichio::Equation;
    ///
    /// let equation_str = r"H_2 + O_2 \longrightarrow H_2O";
    /// let equation = Equation::from_latex(equation_str).unwrap();
    ///
    /// assert_eq!(equation.original_str(), equation_str);
    /// ```
    pub fn original_str(&self) -> &str {
        &self.original_str
    }

    /// Returns the vector of reactants
    /// # Returns
    /// * `&Vec<Compound>` - vector of reactants
    /// # Example
    /// ```
    /// use stoichio::{Compound, Equation};
    ///
    /// let equation_str = r"H_2 + O_2 \longrightarrow H_2O";
    /// let equation = Equation::from_latex(equation_str).unwrap();
    ///
    /// let expected = vec![
    ///     Compound::from_latex("H_2").unwrap(),
    ///     Compound::from_latex("O_2").unwrap(),
    /// ];
    ///
    /// assert_eq!(equation.reactants(), &expected);
    /// ```
    pub fn reactants(&self) -> &Vec<Compound> {
        &self.reactants
    }

    /// Returns the vector of products
    /// # Returns
    /// * `&Vec<Compound>` - vector of products
    /// # Example
    /// ```
    /// use stoichio::{Compound, Equation};
    ///
    /// let equation_str = r"H_2 + O_2 \longrightarrow H_2O";
    /// let equation = Equation::from_latex(equation_str).unwrap();
    ///
    /// let expected = vec![Compound::from_latex("H_2O").unwrap()];
    ///
    /// assert_eq!(equation.products(), &expected);
    /// ```
    pub fn products(&self) -> &Vec<Compound> {
        &self.products
    }

    /// Returns the vector of solutions for reactants (stoichiometric coefficients)
    /// It might include extra solution (stoichiometric coefficient for electrons if needed)
    /// # Returns
    /// * `Option<&Vec<i64>>` - vector of solutions for reactants
    /// # Example
    /// ```
    /// use stoichio::Equation;
    ///
    /// let equation_str = r"H_2 + O_2 \longrightarrow H_2O";
    /// let mut equation = Equation::from_latex(equation_str).unwrap();
    /// equation.solve().unwrap();
    ///
    /// let solutions = equation.solution_reactants().unwrap();
    /// assert_eq!(solutions, &[2, 1]);
    pub fn solution_reactants(&self) -> Option<&Vec<i64>> {
        self.solutions_reactants.as_ref()
    }

    /// Returns the vector of solutions for products (stoichiometric coefficients)
    /// It might include extra solution (stoichiometric coefficient for electrons if needed)
    /// # Returns
    /// * `Option<&Vec<i64>>` - vector of solutions for products
    /// # Example
    /// ```
    /// use stoichio::Equation;
    ///
    /// let equation_str = r"H_2 + O_2 \longrightarrow H_2O";
    /// let mut equation = Equation::from_latex(equation_str).unwrap();
    /// equation.solve().unwrap();
    ///
    /// let solutions = equation.solution_products().unwrap();
    /// assert_eq!(solutions, &[2]);
    /// ```
    pub fn solution_products(&self) -> Option<&Vec<i64>> {
        self.solutions_products.as_ref()
    }

    /// Returns the solution of the equation as a string
    /// # Returns
    /// * `Option<String>` - solution of the equation as a string
    /// # Example
    /// ```
    /// use stoichio::Equation;
    ///
    /// let equation_str = r"H_{2} + O_{2} \longrightarrow H_{2}O";
    /// let mut equation = Equation::from_latex(equation_str).unwrap();
    /// equation.solve().unwrap();
    ///
    /// assert_eq!(equation.solution_str().unwrap(), "2H_{2} + O_{2} \\longrightarrow 2H_{2}O");
    /// ```
    pub fn solution_str(&self) -> Option<String> {
        // check if solutions are available
        if self.solutions_reactants.is_none() || self.solutions_products.is_none() { return None; }

        // get references to solutions
        let sols_reacts = self.solutions_reactants.as_ref().unwrap();
        let sols_prods = self.solutions_products.as_ref().unwrap();

        // generate strings for reactants and products
        let mut reactants_str = String::new();
        let mut products_str = String::new();

        for (i, (reactant, quantity)) in zip(self.reactants.iter(), sols_reacts.iter()).enumerate() {
            if i != 0 { reactants_str.push_str(" + "); }
            if *quantity != 1 {
                reactants_str.push_str(&quantity.to_string());
            }
            reactants_str.push_str(&reactant.original_str);
        }

        for (i, (product, quantity)) in zip(self.products.iter(), sols_prods.iter()).enumerate() {
            if i != 0 { products_str.push_str(" + "); }
            if *quantity != 1 {
                products_str.push_str(&quantity.to_string());
            }
            products_str.push_str(&product.original_str);
        }

        // strings generated don't include electrons, we will add them now (if needed)
        if self.reactants.len() < sols_reacts.len() {
            reactants_str.push_str(&format!(" + {}e^{{-}}", if sols_reacts[sols_reacts.len() - 1] > 1 { sols_reacts[sols_reacts.len() - 1].to_string() } else { "".to_string() }));
        } else if self.products.len() < sols_prods.len() {
            products_str.push_str(&format!(" + {}e^{{-}}", if sols_prods[sols_prods.len() - 1] > 1 { sols_prods[sols_prods.len() - 1].to_string() } else { "".to_string() }));
        }

        Some(format!("{} \\longrightarrow {}", reactants_str, products_str))
    }
}

/// A struct that represents a chemical compound (e.g. H2O, NaCl, ...)
/// # Example
/// ```
/// use stoichio::Compound;
///
/// let compound_str = r"H_{2}O";
/// let compound = Compound::from_latex(compound_str).unwrap();
///
/// assert_eq!(compound.original_str(), compound_str);
/// ```
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Compound {
    /// String from which the compound was parsed
    original_str: String,
    /// HashMap of elements and their quantities
    elements: HashMap<Element, i64>,
    /// Offset of electrons in the compound
    electron_offset: i64,
}
impl Compound {
    /// Create new compound from LaTeX string
    /// # Arguments
    /// * `input` - LaTeX string
    /// # Returns
    /// * `Ok` - compound
    /// * `Err` - error that occurred during parsing
    /// # Example
    /// ```
    /// use stoichio::Compound;
    /// use mendeleev::Element;
    /// use std::collections::HashMap;
    ///
    /// let compound_str = r"C_{2}H_{5}OH";
    /// let compound = Compound::from_latex(compound_str).unwrap();
    ///
    /// let expected_elements = HashMap::from([(Element::C, 2), (Element::H, 6), (Element::O, 1)]);
    ///
    /// assert_eq!(compound.elements(), &expected_elements);
    /// assert_eq!(compound.electrons(), 26);
    /// assert_eq!(compound.electron_offset(), 0);
    /// assert_eq!(compound.original_str(), compound_str);
    /// ```
    pub fn from_latex(input: &str) -> Result<Self, Errors> {
        let original_str = input.to_string();
        let input = input.replace('[', "(");
        let input = input.replace(']', ")");

        if (input.chars().filter(|&c| c == '{').count() != input.chars().filter(|&c| c == '}').count()) ||
            (input.chars().filter(|&c| c == '(').count() != input.chars().filter(|&c| c == ')').count()) {

            return Err(Errors::InvalidEquation);
        }

        #[derive(Copy, Clone, Debug)]
        enum Token {
            Element(Element),
            Number(i64),
            Superscript,
            Subscript,
            Quantity(i64),
            Electrons(i64),
            OpenParentheses,
            CloseParentheses,
        }
        let mut letters: Vec<char> = input.chars().collect();
        let mut tokens: Vec<Token> = Vec::new();

        while !letters.is_empty() {
            match &letters[0] {
                '(' => {
                    tokens.push(Token::OpenParentheses);
                    letters.remove(0);
                },
                ')' => {
                    tokens.push(Token::CloseParentheses);
                    letters.remove(0);
                },
                ('A'..='Z') => {
                    match letters.get(1) {
                        Some(next_letter) => {
                            if next_letter.is_ascii_lowercase() {
                                let pattern = format!("{}{}", letters[0], next_letter);
                                tokens.push(Token::Element(*ALL_ELEMENTS.iter().find(|e| e.symbol() == pattern).ok_or(Errors::InvalidElement)?));
                                letters.remove(0);
                                letters.remove(0);
                            } else {
                                let pattern = letters[0].to_string();
                                tokens.push(Token::Element(*ALL_ELEMENTS.iter().find(|e| e.symbol() == pattern).ok_or(Errors::InvalidElement)?));
                                letters.remove(0);
                            }
                        },
                        None => {
                            let pattern = letters[0].to_string();
                            tokens.push(Token::Element(*ALL_ELEMENTS.iter().find(|e| e.symbol() == pattern).ok_or(Errors::InvalidEquation)?));
                            letters.remove(0);
                        },
                    }
                },
                '{' => {
                    match letters.iter().find(|&&c| c == '}') {
                        Some(_) => {
                            let mut quantity = String::new();
                            letters.remove(0);
                            while letters[0] != '}' {
                                quantity.push(letters[0]);
                                letters.remove(0);
                            }
                            letters.remove(0);
                            match quantity.parse::<i64>() {
                                Ok(num) => tokens.push(Token::Number(num)),
                                Err(_) => {
                                    if !quantity.ends_with('+') && !quantity.ends_with('-') {
                                        return Err(Errors::InvalidEquation);
                                    }
                                    if quantity.chars().count() == 1 {
                                        quantity.insert(0, '1');
                                    }
                                    let last_char = quantity.pop().unwrap();
                                    quantity.insert(0, last_char);
                                    tokens.push(Token::Number(quantity.parse::<i64>().map_err(|_| Errors::InvalidEquation)?));
                                },
                            }
                        },
                        None => return Err(Errors::InvalidEquation),
                    }
                },
                '^' => {
                    tokens.push(Token::Superscript);
                    letters.remove(0);
                },
                '_' => {
                    tokens.push(Token::Subscript);
                    letters.remove(0);
                },
                '+' => {
                    tokens.push(Token::Number(1));
                    letters.remove(0);
                },
                '-' => {
                    tokens.push(Token::Number(-1));
                    letters.remove(0);
                },
                ('0'..='9') => {
                    tokens.push(Token::Number(letters[0].to_digit(10).unwrap() as i64));
                    letters.remove(0);
                }
                _ => return Err(Errors::InvalidEquation),
            }
        }

        while let Some(pos) = tokens.iter().position(|t| matches!(t, Token::Superscript)) {
            let next_pos = tokens.get(pos + 1).copied();
            match next_pos {
                Some(Token::Number(num)) => {
                    tokens.remove(pos + 1);
                    tokens.remove(pos);
                    tokens.insert(pos, Token::Electrons(num));
                },
                Some(_) | None => return Err(Errors::InvalidEquation),
            }
        }

        while let Some(pos) = tokens.iter().position(|t| matches!(t, Token::Subscript)) {
            let next_pos = tokens.get(pos + 1).copied();
            match next_pos {
                Some(Token::Number(num)) => {
                    tokens.remove(pos + 1);
                    tokens.remove(pos);
                    tokens.insert(pos, Token::Quantity(num));
                },
                Some(_) | None => return Err(Errors::InvalidEquation),
            }
        }

        for token in &tokens {
            match token {
                Token::Element(_) => {},
                Token::Quantity(_) => {},
                Token::Electrons(_) => {},
                Token::OpenParentheses => {},
                Token::CloseParentheses => {},
                _ => return Err(Errors::InvalidEquation),
            }
        }

        let mut added = 0;
        for i in 0..tokens.len() {
            if let Token::Element(_) = tokens[i + added] {
                let next_token = tokens.get(i + added + 1).copied();
                match next_token {
                    Some(Token::Quantity(_)) => {},
                    Some(_) => {
                        tokens.insert(i + added + 1, Token::Quantity(1));
                        added += 1;
                    },
                    None => {
                        tokens.push(Token::Quantity(1));
                    },
                }
            }
        }

        let mut changes = true;
        while changes {
            changes = false;
            for i in 0..(tokens.len() - 1) {
                if let Token::Electrons(_) = tokens[i] {
                    if let Token::Quantity(_) = tokens[i + 1] {
                        (tokens[i], tokens[i + 1]) = (tokens[i + 1], tokens[i]);
                        changes = true;
                    }
                }
            }
        }

        let mut parentheses_stack = Vec::new();
        let mut deleted = 0;
        let mut skip = 0;
        for i in 0..tokens.len() {
            if skip > 0 {
                skip -= 1;
                continue;
            }
            match tokens[i - deleted] {
                Token::OpenParentheses => parentheses_stack.push(i - deleted),
                Token::CloseParentheses => {
                    let open_parentheses_pos = parentheses_stack.pop().ok_or(Errors::InvalidEquation)?;
                    let next_token = tokens.get(i - deleted + 1).copied();
                    match next_token {
                        Some(Token::Quantity(num)) => {
                            for inner_token_pos in (open_parentheses_pos + 1)..(i - deleted) {
                                let inner_token = tokens[inner_token_pos];
                                match inner_token {
                                    Token::Quantity(q) => {
                                        tokens[inner_token_pos] = Token::Quantity(q * num);
                                    },
                                    Token::Electrons(e) => {
                                        tokens[inner_token_pos] = Token::Electrons(e * num);
                                    },
                                    _ => {},
                                }
                            }
                            tokens.remove(i - deleted);
                            tokens.remove(i - deleted);
                            tokens.remove(open_parentheses_pos);
                            deleted += 3;
                            skip = 1;
                        },
                        Some(Token::Electrons(_)) => {
                            (tokens[i - deleted], tokens[i - deleted + 1]) = (tokens[i - deleted + 1], tokens[i - deleted]);
                            parentheses_stack.push(open_parentheses_pos);
                        },
                        Some(_) | None => {
                            tokens.remove(i - deleted);
                            tokens.remove(open_parentheses_pos);
                            deleted += 2;
                        },
                    }
                },
                _ => {},
            }
        }

        if !parentheses_stack.is_empty() { return Err(Errors::InvalidEquation); }

        for token in &tokens {
            match token {
                Token::Element(_) => {},
                Token::Quantity(_) => {},
                Token::Electrons(_) => {},
                _ => return Err(Errors::InvalidEquation),
            }
        }

        let mut electron_offset = 0;
        while let Some(pos) = tokens.iter().position(|t| matches!(t, Token::Electrons(_))) {
            match tokens[pos] {
                Token::Electrons(num) => {
                    electron_offset += num;
                    tokens.remove(pos);
                },
                _ => unreachable!(),
            }
        }

        // count elements
        let mut elements = HashMap::new();
        let mut skip = 0;
        for i in 0..tokens.len() {
            if skip > 0 {
                skip -= 1;
                continue;
            }
            if let Token::Element(e) = tokens[i] {
                let next_token = tokens.get(i + 1).copied();
                match next_token {
                    Some(Token::Quantity(q)) => {
                        *elements.entry(e).or_insert(0) += q;
                        skip = 1;
                    },
                    Some(Token::Element(_)) => {},
                    Some(_) => return Err(Errors::InvalidEquation),
                    None => {
                        *elements.entry(e).or_insert(0) += 1;
                    },
                }
            } else {
                return Err(Errors::InvalidEquation);
            }
        }

        Ok(Self {
            original_str,
            elements,
            electron_offset,
        })
    }

    /// Returns the original string from which the compound was parsed
    /// # Returns
    /// * `&str` - original string
    /// # Example
    /// ```
    /// use stoichio::Compound;
    ///
    /// let compound_str = r"H_{2}O";
    /// let compound = Compound::from_latex(compound_str).unwrap();
    ///
    /// assert_eq!(compound.original_str(), compound_str);
    /// ```
    pub fn original_str(&self) -> &str {
        &self.original_str
    }

    /// Returns the HashMap of elements and their quantities in the compound
    /// For example, in the compound H_{2}O the HashMap will be {H: 2, O: 1}
    /// # Returns
    /// * `&HashMap<Element, i64>` - HashMap of elements and their quantities
    /// # Example
    /// ```
    /// use stoichio::Compound;
    /// use mendeleev::Element;
    /// use std::collections::HashMap;
    ///
    /// let compound_str = r"H_{2}O";
    /// let compound = Compound::from_latex(compound_str).unwrap();
    ///
    /// let expected = HashMap::from([(Element::H, 2), (Element::O, 1)]);
    /// assert_eq!(compound.elements(), &expected);
    /// ```
    pub fn elements(&self) -> &HashMap<Element, i64> {
        &self.elements
    }

    /// Returns the number of electrons in the compound
    /// This takes into account the electron offset
    /// For example, in the compound SO_4^{2-} the number of electrons is 50 (the compound has 2 electrons more than it should have)
    /// # Returns
    /// * `i64` - number of electrons
    /// # Example
    /// ```
    /// use stoichio::Compound;
    ///
    /// let compound_str = r"SO_{4}^{2-}";
    /// let compound = Compound::from_latex(compound_str).unwrap();
    ///
    /// assert_eq!(compound.electrons(), 50);
    /// ```
    pub fn electrons(&self) -> i64 {
        self.elements.iter().map(|(e, q)| i64::from(e.atomic_number()) * q).sum::<i64>() - self.electron_offset
    }

    /// Returns the offset of electrons in the compound
    /// For example, in the compound K^{+} the offset is 1 (the compound has 1 electron less than it should have)
    /// # Returns
    /// * `i64` - offset of electrons
    /// # Example
    /// ```
    /// use stoichio::Compound;
    ///
    /// let compound_str = r"K^{+}";
    /// let compound = Compound::from_latex(compound_str).unwrap();
    ///
    /// assert_eq!(compound.electron_offset(), 1);
    /// ```
    pub fn electron_offset(&self) -> i64 {
        self.electron_offset
    }
}





/// Solves system of linear equations.
/// # Arguments
/// * `matrix` - matrix of coefficients of equations, last column is used to store solutions
/// # Returns
/// * `Ok` - vector of solutions
/// * `Err` - error that occurred during solving equations
/// # Example
/// ```
/// use stoichio::solve_equations;
///
/// // we will solve this system of linear equations:
/// // 2x + y - z = 8
/// // -3x - y + 2z = -11
/// // -2x + y + 2z = -3
///
/// // the matrix of coefficients is:
/// //  2   1  -1   8
/// // -3  -1   2 -11
/// // -2   1   2  -3
/// let matrix = vec![
///     vec![2, 1, -1, 8],
///     vec![-3, -1, 2, -11],
///     vec![-2, 1, 2, -3],
/// ];
///
/// // we will get vector of solutions:
/// // [2, 3, -1]
/// // this means that x = 2, y = 3, z = -1
///
/// let solutions = solve_equations(&matrix).unwrap();
/// assert_eq!(solutions, vec![2, 3, -1]);
/// ```
pub fn solve_equations(matrix: &Vec<Vec<i64>>) -> Result<Vec<i64>, Errors> {
    // compute dimensions (m x n)
    // m - number of rows
    // n - number of columns
    let mut m = matrix.len();  // rows
    if m == 0 { return Err(Errors::WrongMatrixDimensions); }

    let mut n = matrix[0].len();  // cols
    if n <= 1 { return Err(Errors::WrongMatrixDimensions); }
    for row in matrix.iter().skip(1) {
        if row.len() != n { return Err(Errors::WrongMatrixDimensions); }  // loop through other rows and check if they have the same number of columns as the first one
    }
    n -= 1;  // last column is used only to store solutions so we don't count it as a column

    // if there are more columns than rows, then there are more unknowns than equations, therefore there is no solution
    if m < n { return Err(Errors::NoSystemSolution); }

    // if there are more rows than columns, then there are more equations than unknowns which usually means that the system is overdetermined and there is no solution
    // but that we will check later

    // we create a new matrix that uses Rational numbers instead of i64
    let mut matrix = matrix
        .iter()
        .map(|row| row
            .iter()
            .map(|&x| Rational::from(x))
            .collect()
        )
        .collect::<Vec<Vec<Rational>>>();

    // perform gaussian elimination on the matrix
    gaussian_elimination(&mut matrix, m, n);

    // check if there is no solution
    // there is no solution if there is a column of zeros in the matrix
    for column in 0..n {
        if matrix.iter().all(|row| row[column] == Rational::ZERO) {
            return Err(Errors::NoSystemSolution);
        }
    }

    // remove zero rows (if there are any)
    // any row after the n rows must contain only zeros (because of gaussian elimination)
    // if there are rows in which all coefficients are zero, but the solution is not zero, then there is no solution
    while m > n {
        if matrix[m - 1].iter().any(|x| x != &Rational::ZERO) { return Err(Errors::NoSystemSolution); }
        m -= 1;
        matrix.pop();
    }

    // m and n should be the same by now

    // reduce matrix to reduced row echelon form from which we can just read the solutions
    reduce_row_echelon(&mut matrix, n);

    // store solutions in a vector
    let mut solutions: Vec<Rational> = matrix.iter().map(|row| row[n].clone()).collect();

    // multiply solutions by the least common multiple of denominators to get integer solutions
    let mut lcm = Natural::ONE;
    for solu in solutions.iter() {
        lcm = lcm.lcm(solu.denominator_ref());
    }
    for element in &mut solutions {
        *element = &*element * &Rational::from(&lcm);
    }

    Ok(solutions.iter().map(|x| i64::try_from(x).unwrap()).collect())
}

/// Performs Gaussian elimination on augmented matrices (additional column stores solutions)
/// # Arguments
/// * `matrix` - matrix of coefficients of equations, last column is used to store solutions
/// * `m` - number of rows
/// * `n` - number of columns (without the last column that stores solutions)
/// # Example
/// ```
/// use stoichio::gaussian_elimination;
/// use malachite::Rational;
/// use std::str::FromStr;
///
/// let mut matrix = vec![
///     vec![Rational::from(2), Rational::from(1), Rational::from(-1), Rational::from(8)],
///     vec![Rational::from(-3), Rational::from(-1), Rational::from(2), Rational::from(-11)],
///     vec![Rational::from(-2), Rational::from(1), Rational::from(2), Rational::from(-3)],
/// ];
///
/// gaussian_elimination(&mut matrix, 3, 3);
///
/// assert_eq!(matrix, vec![
///     vec![Rational::from(-3), Rational::from(-1), Rational::from(2), Rational::from(-11)],
///     vec![Rational::from(0), Rational::from_str("5/3").unwrap(), Rational::from_str("2/3").unwrap(), Rational::from_str("13/3").unwrap()],
///     vec![Rational::from(0), Rational::from(0), Rational::from_str("1/5").unwrap(), Rational::from_str("-1/5").unwrap()],
/// ]);
/// ```
pub fn gaussian_elimination(matrix: &mut [Vec<Rational>], m: usize, n: usize) {
    let mut row = 0;
    let mut col = 0;
    while row < m && col < n {
        let mut i_max = row;
        for (i, row_n) in matrix.iter().enumerate().skip(row + 1) {
            if (&row_n[col]).abs() > (&matrix[row][col]).abs() {
                i_max = i;
            }
        }

        if matrix[i_max][col] == 0 {
            col += 1;
        } else {
            swap_rows(row, i_max, matrix);
            for i in (row + 1)..m {
                let f = &matrix[i][col] / &matrix[row][col];
                matrix[i][col] = Rational::ZERO;
                for j in (col + 1)..(n + 1) {
                    let div_amount = &f * &matrix[row][j];
                    matrix[i][j] -= div_amount;
                }
            }

            row += 1;
            col += 1;
        }
    }
}

/// Swaps two rows in a matrix
/// Used in Gaussian elimination
/// # Arguments
/// * `r1` - index of the first row
/// * `r2` - index of the second row
/// * `matrix` - matrix of coefficients of equations, last column is used to store solutions
/// # Example
/// ```
/// use stoichio::swap_rows;
///
/// let mut matrix = vec![
///     vec![1, 2, 3],
///     vec![4, 5, 6],
///     vec![7, 8, 9],
/// ];
///
/// swap_rows(0, 2, &mut matrix);
///
/// assert_eq!(matrix, vec![
///     vec![7, 8, 9],
///     vec![4, 5, 6],
///     vec![1, 2, 3],
/// ]);
/// ```
#[inline(always)]
pub fn swap_rows<T>(r1: usize, r2: usize, matrix: &mut [Vec<T>]) {
    if r1 != r2 {
        let bigger_r = max(r1, r2);
        let smaller_r = min(r1, r2);
        let (top, bot) = matrix.split_at_mut(bigger_r);  // splits before bigger_r so index 0 in bot will be bigger_r
        mem::swap(&mut top[smaller_r], &mut bot[0])
    }
}

/// Reduces matrix in row echelon form to reduced row echelon form
/// # Arguments
/// * `matrix` - matrix of coefficients of equations, last column is used to store solutions
/// * `n` - number of rows and columns (should be the same by now) (without the last column that stores solutions)
/// # Example
/// ```
/// use stoichio::reduce_row_echelon;
/// use malachite::Rational;
/// use std::str::FromStr;
///
/// let mut matrix = vec![
///     vec![Rational::from(-3), Rational::from(-1), Rational::from(2), Rational::from(-11)],
///     vec![Rational::from(0), Rational::from_str("5/3").unwrap(), Rational::from_str("2/3").unwrap(), Rational::from_str("13/3").unwrap()],
///     vec![Rational::from(0), Rational::from(0), Rational::from_str("1/5").unwrap(), Rational::from_str("-1/5").unwrap()],
/// ];
///
/// reduce_row_echelon(&mut matrix, 3);
///
/// assert_eq!(matrix, vec![
///     vec![Rational::from(1), Rational::from(0), Rational::from(0), Rational::from(2)],
///     vec![Rational::from(0), Rational::from(1), Rational::from(0), Rational::from(3)],
///     vec![Rational::from(0), Rational::from(0), Rational::from(1), Rational::from(-1)],
/// ]);
/// ```
pub fn reduce_row_echelon(matrix: &mut [Vec<Rational>], n: usize) {
    for r in (1..n).rev() {
        for r_temp in 0..r {
            // only process cells in a column above leading coefficient of the current row and the result column
            // because other coefficients are already zero (we are processing matrix from bottom to top)
            let factor = &matrix[r_temp][r] / &matrix[r][r];
            matrix[r_temp][r] = Rational::ZERO;
            let div_amount = &factor * &matrix[r][n];
            matrix[r_temp][n] -= div_amount;
        }
    }

    for (i, row) in matrix.iter_mut().enumerate() {
        let factor = Rational::ONE / &row[i];
        row[i] = Rational::ONE;
        row[n] *= factor;
    }
}
