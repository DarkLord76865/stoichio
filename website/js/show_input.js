const input = document.getElementById('input-equation');
const equation = document.getElementById('equation-1');

input.addEventListener('input', function() {
    let text_value = input.value;
    if (text_value === "") {
        text_value = "-";
    }
    equation.textContent = "$$" + text_value + "$$";
    MathJax.typeset([equation]);
});