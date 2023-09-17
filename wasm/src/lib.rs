use wasm_bindgen::prelude::*;
use stoichio::Equation;


// Export a `equation_io` function from Rust to JavaScript.
#[wasm_bindgen]
/// Process input equation and return output
pub fn equation_io(equation: &str) -> String {
    let mut equation = match Equation::from_latex(equation) {
        Ok(equation) => equation,
        Err(err) => return format!("0{}", err),
    };

    // first char is 1 if success, 0 if error
    match equation.solve() {
        Ok(_) => format!("1{}", equation.solution_str().unwrap()),
        Err(err) => format!("0{}", err),
    }
}
