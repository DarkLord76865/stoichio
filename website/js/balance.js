function balance() {
    const input = document.getElementById('input-equation');
    const output_equation = document.getElementById('equation-2');

    let result_text = equation_io(input.value);
    let result = "";
    if (result_text[0] === "1") {
        result += "$$ " + result_text.substring(1) + " $$";
    } else {
        result += result_text.substring(1);
    }

    output_equation.textContent = result;
    MathJax.typeset([output_equation]);
}
