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





/// Horizontal arrows in LaTeX
pub const LATEX_ARROWS: [&str; 23] = [
    r"\leftarrow",
    r"\rightarrow",
    r"\longleftarrow",
    r"\longrightarrow",
    r"\Leftarrow",
    r"\Rightarrow",
    r"\Longleftarrow",
    r"\Longrightarrow",
    r"\leftrightarrow",
    r"\longleftrightarrow",
    r"\Leftrightarrow",
    r"\Longleftrightarrow",
    r"\leftrightarrows",
    r"\leftharpoonup",
    r"\leftharpoondown",
    r"\rightharpoonup",
    r"\rightharpoondown",
    r"\leftrightharpoons",
    r"\rightleftharpoons",
    r"\hookleftarrow",
    r"\hookrightarrow",
    r"\mapsto",
    r"\longmapsto",
];





/// Errors that can occur during solving system of linear equations
#[derive(Copy, Clone, Debug, Eq, Hash, PartialEq)]
pub enum StoichioError {
    /// Entered equation is invalid
    InvalidEquation,
    /// Solution was calculated, but is invalid
    InvalidSolution,
    /// Entered element is invalid
    InvalidElement,
    /// There should be exactly one arrow in the equation
    InvalidArrowCount,
    /// Entered reaction environment is invalid
    InvalidEnvironment,

    /// Matrix has wrong dimensions (rows and columns)
    WrongMatrixDimensions,
    /// There is no solution to the system of linear equations
    NoSystemSolution,
}
impl Display for StoichioError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            StoichioError::WrongMatrixDimensions => write!(f, "Wrong matrix dimensions"),
            StoichioError::NoSystemSolution => write!(f, "No solution"),
            StoichioError::InvalidElement => write!(f, "Invalid element"),
            StoichioError::InvalidEquation => write!(f, "Invalid equation"),
            StoichioError::InvalidSolution => write!(f, "Invalid solution"),
            StoichioError::InvalidArrowCount => write!(f, "There should be exactly one arrow in the equation"),
            StoichioError::InvalidEnvironment => write!(f, "Entered reaction environment is invalid"),
        }
    }
}
impl Error for StoichioError {}

/// Stores the information about the environment in which the reaction is happening
#[derive(Copy, Clone, Debug, Eq, Hash, PartialEq)]
pub enum EnvironmentType {
    /// Acidic environment
    Acidic,
    /// Basic environment
    Basic,
    /// Neutral environment
    Neutral
}





/// A struct that represents a chemical equation (e.g. 2H2 + O2 -> 2H2O)
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Equation {
    /// String from which the equation was parsed
    original_str: String,
    /// Arrow type used in the equation
    arrow_type: String,
    /// Environment in which the reaction is happening
    environment_type: EnvironmentType,
    /// Stores whether the equation is a half-reaction (the electrons are included in the equation)
    half_reaction: bool,
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
    /// Create new equation from LaTeX string
    /// The equation should contain exactly one arrow
    /// If you need to specify in which environment the reaction is happening, add (acid.) for acidic environment and (base.) for basic environment to the end of equation
    /// If you want to have electrons in the equation (for half-reactions), use e^- somewhere in the equation
    /// Note that electrons should not be used in the equation if it is not a half-reaction, so if you specify environment, you should not use electrons
    /// # Arguments
    /// * `input` - LaTeX string
    /// # Returns
    /// * `Ok` - equation
    /// * `Err` - error that occurred during parsing
    /// # Example
    /// ```
    /// use stoichio::{Compound, Equation};
    ///
    /// let equation_str = r"H_{2} + O_{2} \longrightarrow H_{2}O";
    /// let equation = Equation::from_latex(equation_str).unwrap();
    ///
    /// let expected_reactants = vec![
    ///     Compound::from_latex("H_{2}").unwrap(),
    ///     Compound::from_latex("O_{2}").unwrap(),
    /// ];
    /// let expected_products = vec![Compound::from_latex("H_{2}O").unwrap()];
    ///
    /// assert_eq!(equation.original_str(), equation_str);
    /// assert_eq!(equation.reactants(), &expected_reactants);
    /// assert_eq!(equation.products(), &expected_products);
    /// ```
    pub fn from_latex(input: &str) -> Result<Self, StoichioError> {
        // store original string
        let original_str = input;

        // get input as String
        let input = input.to_string();

        // check input for arrows
        let mut latex_arrow_counts: [usize; 23] = [0; 23];
        for (i, arrow_type) in LATEX_ARROWS.iter().enumerate() {
            latex_arrow_counts[i] += input.split(arrow_type).count() - 1;
        }
        latex_arrow_counts[8] -= latex_arrow_counts[12];  // at this point, latex_arrow_counts[12] is the number of \leftrightarrows, we need to subtract it from the number of \leftrightarrow since it is a subset of \leftrightarrow
        if latex_arrow_counts.iter().sum::<usize>() != 1 { return Err(StoichioError::InvalidArrowCount); }
        let arrow_type = LATEX_ARROWS[latex_arrow_counts.iter().position(|&x| x == 1).unwrap()];

        // check input for environment (acidic or basic)
        let base = input.contains("(base.)");
        let acid = input.contains("(acid.)");
        if base && acid { return Err(StoichioError::InvalidEnvironment); }
        let environment_type = if base { EnvironmentType::Basic } else if acid { EnvironmentType::Acidic } else { EnvironmentType::Neutral };
        let input = input.replace("(base.)", "");
        let input = input.replace("(acid.)", "");

        // split input string into reactants and products
        let reactants_str = input.split(arrow_type).next().unwrap();
        let products_str = input.split(arrow_type).last().unwrap();

        // flag that stores whether the equation is a half-reaction (the electrons are included in the equation)
        // it is set to true if electrons are detected in the equation
        let mut half_reaction = false;

        // cleans up one side of the equation
        let sanitize_side = |string: &str| -> String {

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

        // processes one side of the equation
        let mut process_side = |string: &str| -> Result<Vec<Compound>, StoichioError> {
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
                } else {
                    half_reaction = true;
                }
            }

            Ok(compounds)
        };

        // get vectors of reactants and products
        let reactants = process_side(&sanitize_side(reactants_str))?;
        let products = process_side(&sanitize_side(products_str))?;

        Ok(Self {
            original_str: String::from(original_str),
            arrow_type: String::from(arrow_type),
            environment_type,
            half_reaction,
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
    pub fn solve(&mut self) -> Result<(), StoichioError> {
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
        if solutions.iter().any(|x| *x <= 0) { return Err(StoichioError::InvalidSolution); }

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
        if reactants_element_counts != products_element_counts { return Err(StoichioError::InvalidSolution); }

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

    /// Returns the arrow type used in the equation
    /// # Returns
    /// * `&str` - arrow type
    /// # Example
    /// ```
    /// use stoichio::Equation;
    /// use stoichio::LATEX_ARROWS;
    ///
    /// let equation_str = r"H_2 + O_2 \longrightarrow H_2O";
    /// let equation = Equation::from_latex(equation_str).unwrap();
    ///
    /// let arrow_type = equation.arrow_type();
    ///
    /// assert_eq!(arrow_type, r"\longrightarrow");
    /// assert!(LATEX_ARROWS.contains(&arrow_type));
    /// ```
    pub fn arrow_type(&self) -> &str {
        &self.arrow_type
    }

    /// Returns the environment in which the reaction is happening
    /// # Returns
    /// * `EnvironmentType` - environment in which the reaction is happening
    /// # Example
    /// ```
    /// use stoichio::{Equation, EnvironmentType};
    ///
    /// let equation_str = r"H_2 + O_2 \longrightarrow H_2O (acid.)";
    /// let equation = Equation::from_latex(equation_str).unwrap();
    ///
    /// assert_eq!(equation.environment_type(), EnvironmentType::Acidic);
    ///
    /// let equation_str = r"H_2 + O_2 \longrightarrow H_2O (base.)";
    /// let equation = Equation::from_latex(equation_str).unwrap();
    ///
    /// assert_eq!(equation.environment_type(), EnvironmentType::Basic);
    ///
    /// let equation_str = r"H_2 + O_2 \longrightarrow H_2O";
    /// let equation = Equation::from_latex(equation_str).unwrap();
    ///
    /// assert_eq!(equation.environment_type(), EnvironmentType::Neutral);
    /// ```
    pub fn environment_type(&self) -> EnvironmentType {
        self.environment_type
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

        Some(format!("{} {} {}", reactants_str, self.arrow_type, products_str))
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
    pub fn from_latex(input: &str) -> Result<Self, StoichioError> {
        let original_str = input.to_string();
        let input = input.replace('[', "(");
        let input = input.replace(']', ")");

        if (input.chars().filter(|&c| c == '{').count() != input.chars().filter(|&c| c == '}').count()) ||
            (input.chars().filter(|&c| c == '(').count() != input.chars().filter(|&c| c == ')').count()) {

            return Err(StoichioError::InvalidEquation);
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
                                tokens.push(Token::Element(*ALL_ELEMENTS.iter().find(|e| e.symbol() == pattern).ok_or(StoichioError::InvalidElement)?));
                                letters.remove(0);
                                letters.remove(0);
                            } else {
                                let pattern = letters[0].to_string();
                                tokens.push(Token::Element(*ALL_ELEMENTS.iter().find(|e| e.symbol() == pattern).ok_or(StoichioError::InvalidElement)?));
                                letters.remove(0);
                            }
                        },
                        None => {
                            let pattern = letters[0].to_string();
                            tokens.push(Token::Element(*ALL_ELEMENTS.iter().find(|e| e.symbol() == pattern).ok_or(StoichioError::InvalidEquation)?));
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
                                        return Err(StoichioError::InvalidEquation);
                                    }
                                    if quantity.chars().count() == 1 {
                                        quantity.insert(0, '1');
                                    }
                                    let last_char = quantity.pop().unwrap();
                                    quantity.insert(0, last_char);
                                    tokens.push(Token::Number(quantity.parse::<i64>().map_err(|_| StoichioError::InvalidEquation)?));
                                },
                            }
                        },
                        None => return Err(StoichioError::InvalidEquation),
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
                _ => return Err(StoichioError::InvalidEquation),
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
                Some(_) | None => return Err(StoichioError::InvalidEquation),
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
                Some(_) | None => return Err(StoichioError::InvalidEquation),
            }
        }

        for token in &tokens {
            match token {
                Token::Element(_) => {},
                Token::Quantity(_) => {},
                Token::Electrons(_) => {},
                Token::OpenParentheses => {},
                Token::CloseParentheses => {},
                _ => return Err(StoichioError::InvalidEquation),
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
                    let open_parentheses_pos = parentheses_stack.pop().ok_or(StoichioError::InvalidEquation)?;
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

        if !parentheses_stack.is_empty() { return Err(StoichioError::InvalidEquation); }

        for token in &tokens {
            match token {
                Token::Element(_) => {},
                Token::Quantity(_) => {},
                Token::Electrons(_) => {},
                _ => return Err(StoichioError::InvalidEquation),
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
                    Some(_) => return Err(StoichioError::InvalidEquation),
                    None => {
                        *elements.entry(e).or_insert(0) += 1;
                    },
                }
            } else {
                return Err(StoichioError::InvalidEquation);
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
pub fn solve_equations(matrix: &Vec<Vec<i64>>) -> Result<Vec<i64>, StoichioError> {
    // compute dimensions (m x n)
    // m - number of rows
    // n - number of columns
    let mut m = matrix.len();  // rows
    if m == 0 { return Err(StoichioError::WrongMatrixDimensions); }

    let mut n = matrix[0].len();  // cols
    if n <= 1 { return Err(StoichioError::WrongMatrixDimensions); }
    for row in matrix.iter().skip(1) {
        if row.len() != n { return Err(StoichioError::WrongMatrixDimensions); }  // loop through other rows and check if they have the same number of columns as the first one
    }
    n -= 1;  // last column is used only to store solutions so we don't count it as a column

    // if there are more columns than rows, then there are more unknowns than equations, therefore there is no solution
    if m < n { return Err(StoichioError::NoSystemSolution); }

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
            return Err(StoichioError::NoSystemSolution);
        }
    }

    // remove zero rows (if there are any)
    // any row after the n rows must contain only zeros (because of gaussian elimination)
    // if there are rows in which all coefficients are zero, but the solution is not zero, then there is no solution
    while m > n {
        if matrix[m - 1].iter().any(|x| x != &Rational::ZERO) { return Err(StoichioError::NoSystemSolution); }
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





#[cfg(test)]
mod tests {
    use super::*;

    fn test_equation(equation: &str, solved_equation: &str) {
        let mut eq = Equation::from_latex(equation).unwrap();
        eq.solve().unwrap();
        let solution = eq.solution_str().unwrap();

        assert_eq!(solution, solved_equation);
    }


    #[test]
    fn arrow_types_latex() {
        for arrow in LATEX_ARROWS {
            let eq_str = format!("H_2 + O_2 {} H_2O", arrow);
            let mut eq = Equation::from_latex(&eq_str).unwrap();
            eq.solve().unwrap();

            assert_eq!(eq.arrow_type(), arrow);
            assert_eq!(eq.solution_str().unwrap(), format!("2H_2 + O_2 {} 2H_2O", arrow));
        }
    }

    #[test]
    fn eq1() {
        test_equation(r"H_{2} + O_{2} \longrightarrow H_{2}O", r"2H_{2} + O_{2} \longrightarrow 2H_{2}O");
    }

    #[test]
    fn eq2() {
        test_equation(r"K \longrightarrow K^+", r"K \longrightarrow K^+ + e^{-}");
    }

    #[test]
    fn eq3() {
        test_equation(r"[Cr(N_{2}H_{4}CO)_{6}]_{4}[Cr(CN)_{6}]_{3}+KMnO_{4}+H_{2}SO_{4} \longrightarrow K_{2}Cr_{2}O_{7}+MnSO_{4}+CO_{2}+KNO_{3}+K_{2}SO_{4}+H_{2}O", r"10[Cr(N_{2}H_{4}CO)_{6}]_{4}[Cr(CN)_{6}]_{3} + 1176KMnO_{4} + 1399H_{2}SO_{4} \longrightarrow 35K_{2}Cr_{2}O_{7} + 1176MnSO_{4} + 420CO_{2} + 660KNO_{3} + 223K_{2}SO_{4} + 1879H_{2}O");
    }

    #[test]
    fn eq4() {
        test_equation(r"P_4O_{10} + H_2O \longrightarrow H_3PO_4", r"P_4O_{10} + 6H_2O \longrightarrow 4H_3PO_4");
    }

    #[test]
    fn eq5() {
        test_equation(r"P_4O_{10} + H_2O \longrightarrow H_3PO_4", r"P_4O_{10} + 6H_2O \longrightarrow 4H_3PO_4");
    }

    #[test]
    fn eq6() {
        test_equation(r"CO_2 + H_2O \longrightarrow C_6H_{12}O_6 + O_2", r"6CO_2 + 6H_2O \longrightarrow C_6H_{12}O_6 + 6O_2");
    }

    #[test]
    fn eq7() {
        test_equation(r"SiCl_4 + H_2O \longrightarrow H_4SiO_4 + HCl", r"SiCl_4 + 4H_2O \longrightarrow H_4SiO_4 + 4HCl");
    }

    #[test]
    fn eq8() {
        test_equation(r"Al + HCl \longrightarrow AlCl_3 + H_2", r"2Al + 6HCl \longrightarrow 2AlCl_3 + 3H_2");
    }

    #[test]
    fn eq9() {
        test_equation(r"Na_2CO_3 + HCl \longrightarrow NaCl + H_2O + CO_2", r"Na_2CO_3 + 2HCl \longrightarrow 2NaCl + H_2O + CO_2");
    }

    #[test]
    fn eq10() {
        test_equation(r"C_7H_6O_2 + O_2 \longrightarrow CO_2 + H_2O", r"2C_7H_6O_2 + 15O_2 \longrightarrow 14CO_2 + 6H_2O");
    }

    #[test]
    fn eq11() {
        test_equation(r"Fe_2(SO_4)_3 + KOH \longrightarrow K_2SO_4 + Fe(OH)_3", r"Fe_2(SO_4)_3 + 6KOH \longrightarrow 3K_2SO_4 + 2Fe(OH)_3");
    }

    #[test]
    fn eq12() {
        test_equation(r"Ca_3(PO_4)_2 + SiO_2 \longrightarrow P_4O_{10} + 3CaSiO_3", r"2Ca_3(PO_4)_2 + 6SiO_2 \longrightarrow P_4O_{10} + 6CaSiO_3");
    }

    #[test]
    fn eq13() {
        test_equation(r"KClO_3 \longrightarrow KClO_4 + KCl", r"4KClO_3 \longrightarrow 3KClO_4 + KCl");
    }

    #[test]
    fn eq14() {
        test_equation(r"Al_2(SO_4)_3 + Ca(OH)_2 \longrightarrow Al(OH)_3 + CaSO_4", r"Al_2(SO_4)_3 + 3Ca(OH)_2 \longrightarrow 2Al(OH)_3 + 3CaSO_4");
    }

    #[test]
    fn eq15() {
        test_equation(r"H_2SO_4 + HI \longrightarrow H_2S + I_2 + H_2O", r"H_2SO_4 + 8HI \longrightarrow H_2S + 4I_2 + 4H_2O");
    }

    #[test]
    fn eq16() {
        test_equation(r"C_2H_6 + O_2 \longrightarrow CO_2 + H_2O", r"2C_2H_6 + 7O_2 \longrightarrow 4CO_2 + 6H_2O");
    }

    #[test]
    fn eq17() {
        test_equation(r"NaN_3 \longrightarrow Na + N_2", r"2NaN_3 \longrightarrow 2Na + 3N_2");
    }

    #[test]
    fn eq18() {
        test_equation(r"Na + Fe_2O_3 \longrightarrow Na_2O + Fe", r"6Na + Fe_2O_3 \longrightarrow 3Na_2O + 2Fe");
    }

    #[test]
    fn eq19() {
        test_equation(r"Mg + N_2 \longrightarrow Mg_3N_2", r"3Mg + N_2 \longrightarrow Mg_3N_2");
    }

    #[test]
    fn eq20() {
        test_equation(r"Na + NH_3 \longrightarrow NaNH_2 + H_2", r"2Na + 2NH_3 \longrightarrow 2NaNH_2 + H_2");
    }

    #[test]
    fn eq21() {
        test_equation(r"Na_2O + CO_2 + H_2O \longrightarrow NaHCO_3", r"Na_2O + 2CO_2 + H_2O \longrightarrow 2NaHCO_3");
    }

    #[test]
    fn eq22() {
        test_equation(r"P_4S_3 + O_2 \longrightarrow P_4O_6 + SO_2", r"P_4S_3 + 6O_2 \longrightarrow P_4O_6 + 3SO_2");
    }

    #[test]
    fn eq23() {
        test_equation(r"Na_3PO_4 + CaCl_2 \longrightarrow Ca_3(PO_4)_2 + NaCl", r"2Na_3PO_4 + 3CaCl_2 \longrightarrow Ca_3(PO_4)_2 + 6NaCl");
    }

    #[test]
    fn eq24() {
        test_equation(r"C_8H_{18} + O_2 \longrightarrow CO_2 + H_2O", r"2C_8H_{18} + 25O_2 \longrightarrow 16CO_2 + 18H_2O");
    }

    #[test]
    fn eq25() {
        test_equation(r"C_2H_6O + O_2 \longrightarrow CO_2 + H_2O", r"C_2H_6O + 3O_2 \longrightarrow 2CO_2 + 3H_2O");
    }

    #[test]
    fn eq26() {
        test_equation(r"Pb(NO_3)_2 + KI \longrightarrow PbI_2 + KNO_3", r"Pb(NO_3)_2 + 2KI \longrightarrow PbI_2 + 2KNO_3");
    }

    #[test]
    fn eq27() {
        test_equation(r"N_2O_5 \longrightarrow NO_2 + O_2", r"2N_2O_5 \longrightarrow 4NO_2 + O_2");
    }

    #[test]
    fn eq28() {
        test_equation(r"KClO_3 \longrightarrow KCl + O_2", r"2KClO_3 \longrightarrow 2KCl + 3O_2");
    }

    #[test]
    fn eq29() {
        test_equation(r"CO + O_2 \longrightarrow CO_2", r"2CO + O_2 \longrightarrow 2CO_2");
    }

    #[test]
    fn eq30() {
        test_equation(r"C_{57}H_{110}O_6 + O_2 \longrightarrow CO_2 + H_2O", r"2C_{57}H_{110}O_6 + 163O_2 \longrightarrow 114CO_2 + 110H_2O");
    }

    #[test]
    fn eq31() {
        test_equation(r"K_4[Fe(SCN)_6] + K_2Cr_2O_7 + H_2SO_4 \longrightarrow Fe_2(SO_4)_3 + Cr_2(SO_4)_3 + CO_2 + H_2O + K_2SO_4 + KNO_3", r"6K_4[Fe(SCN)_6] + 97K_2Cr_2O_7 + 355H_2SO_4 \longrightarrow 3Fe_2(SO_4)_3 + 97Cr_2(SO_4)_3 + 36CO_2 + 355H_2O + 91K_2SO_4 + 36KNO_3");
    }

    #[test]
    fn eq32() {
        test_equation(r"Al + H_2SO_4 \longrightarrow Al_2(SO_4)_3 + H_2", r"2Al + 3H_2SO_4 \longrightarrow Al_2(SO_4)_3 + 3H_2");
    }

    #[test]
    fn eq33() {
        test_equation(r"C_7H_{10}N + O_2 \longrightarrow CO_2 + H_2O + NO_2", r"2C_7H_{10}N + 21O_2 \longrightarrow 14CO_2 + 10H_2O + 2NO_2");
    }

    #[test]
    fn eq34() {
        test_equation(r"Al(OH)_3 + H_2SO_4 \longrightarrow Al_2(SO_4)_3 + H_2O", r"2Al(OH)_3 + 3H_2SO_4 \longrightarrow Al_2(SO_4)_3 + 6H_2O");
    }

    #[test]
    fn eq35() {
        test_equation(r"BaO + Al \longrightarrow BaAl_4 + Al_2O_3", r"3BaO + 14Al \longrightarrow 3BaAl_4 + Al_2O_3");
    }

    #[test]
    fn eq36() {
        test_equation(r"AgN_3 \longrightarrow N_2 + Ag", r"2AgN_3 \longrightarrow 3N_2 + 2Ag");
    }

    #[test]
    fn eq37() {
        test_equation(r"Pt + HNO_3 + HCl \longrightarrow H_2PtCl_6 + NO_2 + H_2O", r"Pt + 4HNO_3 + 6HCl \longrightarrow H_2PtCl_6 + 4NO_2 + 4H_2O");
    }

    #[test]
    fn eq38() {
        test_equation(r"LuCl_3 + Ca \longrightarrow Lu + CaCl_2", r"2LuCl_3 + 3Ca \longrightarrow 2Lu + 3CaCl_2");
    }

    #[test]
    fn eq39() {
        test_equation(r"XeF_6 + H_2O \longrightarrow XeO_3 + HF", r"XeF_6 + 3H_2O \longrightarrow XeO_3 + 6HF");
    }

    #[test]
    fn eq40() {
        test_equation(r"Ba_2XeO_6 + H_2SO_4 \longrightarrow BaSO_4 + H_2O + XeO_4", r"Ba_2XeO_6 + 2H_2SO_4 \longrightarrow 2BaSO_4 + 2H_2O + XeO_4");
    }

    #[test]
    fn eq41() {
        test_equation(r"P_4O_6 + H_2O \longrightarrow H_3PO_3", r"P_4O_6 + 6H_2O \longrightarrow 4H_3PO_3");
    }

    #[test]
    fn eq42() {
        test_equation(r"C_6H_{14} + O_2 \longrightarrow CO_2 + H_2O", r"2C_6H_{14} + 19O_2 \longrightarrow 12CO_2 + 14H_2O");
    }

    #[test]
    fn eq43() {
        test_equation(r"MoS_2 + O_2 \longrightarrow MoO_3 + SO_2", r"2MoS_2 + 7O_2 \longrightarrow 2MoO_3 + 4SO_2");
    }

    #[test]
    fn eq44() {
        test_equation(r"K_2MnF_6 + SbF_5 \longrightarrow KSbF_6 + MnF_3 + F_2", r"2K_2MnF_6 + 4SbF_5 \longrightarrow 4KSbF_6 + 2MnF_3 + F_2");
    }

    #[test]
    fn eq45() {
        test_equation(r"S + HNO_3 \longrightarrow H_2SO_4 + NO_2 + H_2O", r"S + 6HNO_3 \longrightarrow H_2SO_4 + 6NO_2 + 2H_2O");
    }

    #[test]
    fn eq46() {
        test_equation(r"Cu + HNO_3 \longrightarrow Cu(NO_3)_2 + NO + H_2O", r"3Cu + 8HNO_3 \longrightarrow 3Cu(NO_3)_2 + 2NO + 4H_2O");
    }

    #[test]
    fn eq47() {
        test_equation(r"CuS + HNO_3 \longrightarrow CuSO_4 + NO_2 + H_2O", r"CuS + 8HNO_3 \longrightarrow CuSO_4 + 8NO_2 + 4H_2O");
    }

    #[test]
    fn eq48() {
        test_equation(r"Cu_2S + HNO_3 \longrightarrow Cu(NO_3)_2 + CuSO_4 + NO_2 + H_2O", r"Cu_2S + 12HNO_3 \longrightarrow Cu(NO_3)_2 + CuSO_4 + 10NO_2 + 6H_2O");
    }

    #[test]
    fn eq49() {
        test_equation(r"NaBr + NaBrO_3 + H_2SO_4 \longrightarrow Br_2 + Na_2SO_4 + H_2O", r"5NaBr + NaBrO_3 + 3H_2SO_4 \longrightarrow 3Br_2 + 3Na_2SO_4 + 3H_2O");
    }

    #[test]
    fn eq50() {
        test_equation(r"KNO_3 + C_{12}H_{22}O_{11} \longrightarrow N_2 + CO_2 + H_2O + K_2CO_3", r"48KNO_3 + 5C_{12}H_{22}O_{11} \longrightarrow 24N_2 + 36CO_2 + 55H_2O + 24K_2CO_3");
    }
}