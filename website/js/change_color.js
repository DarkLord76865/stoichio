// Get references to the button and the target elements
const change_color_button = document.getElementById("change_color_button");
const change_color_button_next_color = document.getElementById("change_color_button_next_color");
const header = document.getElementById("header_id");
const body = document.getElementById("body_id");
const balance_button = document.getElementById("balance-button");
const input_equation = document.getElementById("input-equation");
const description = document.getElementById("description");
const eq1 = document.getElementById("equation-1");
const eq2 = document.getElementById("equation-2");
const input_label = document.getElementById("in-label");
const output_label = document.getElementById("out-label");
const info_label = document.getElementById("info-label");
const lines = document.querySelectorAll(".horizontal-line, .vertical-line");

// Define an array of colors to cycle through, first color is for the header, second color is for the body, third is for the text
const colors = [
    ["#007aff", "#f2f8ff", "#000000"],
    ["#7a00ff", "#eeeeff", "#000000"],
    ["#007aff", "#1D2226", "#ffffff"],
    ["#7a00ff", "#262632", "#ffffff"],
];

// Initialize a color index (color at index 0 is assumed to be the starting color) it is set to the -1 because on the initializing call to set_colors() it will be incremented to 0
let current_color_index = -1;
// Try to get the color index from local storage (if there isn't one, current_color_index will stay at 0)
const stored_index = localStorage.getItem("color_index");
if (stored_index) {
    current_color_index = parseInt(stored_index) - 1;  // Subtract 1 because it will be incremented to 0 on the first call to set_colors()
}

function set_colors() {
    // Increment the color index, looping back to the start if necessary
    current_color_index = (current_color_index + 1) % colors.length;

    // Store the color index in local storage
    localStorage.setItem("color_index", current_color_index);

    // Change the background color
    header.style.backgroundColor = colors[current_color_index][0];
    body.style.backgroundColor = colors[current_color_index][1];
    change_color_button_next_color.style.backgroundColor = colors[(current_color_index + 1) % colors.length][0];
    description.style.color = colors[current_color_index][2];
    eq1.style.color = colors[current_color_index][2];
    eq2.style.color = colors[current_color_index][2];
    input_label.style.color = colors[current_color_index][2];
    output_label.style.color = colors[current_color_index][2];
    balance_button.style.backgroundColor = colors[current_color_index][0];
    input_equation.style.backgroundColor = colors[current_color_index][1];
    input_equation.style.color = colors[current_color_index][2];
    
    if (colors[current_color_index][2] === "#ffffff") {
        balance_button.style.boxShadow = "0 0 10px 0 rgba(255, 255, 255, 0.5)";
        input_equation.style.boxShadow = "0 0 10px 0 rgba(255, 255, 255, 0.5), inset 0 0 10px 0 rgba(255, 255, 255, 0.2)";
        info_label.style.color = "rgb(180, 180, 180)";
        info_label.style.backgroundColor = "rgba(155, 155, 155, 0.1)";
        lines.forEach(line => {
            line.style.backgroundColor = "rgba(155, 155, 155, 0.1)";
        });
    } else {
        balance_button.style.boxShadow = "0 0 10px 0 rgba(0, 0, 0, 0.5)";
        input_equation.style.boxShadow = "0 0 10px 0 rgba(0, 0, 0, 0.5), inset 0 0 10px 0 rgba(0, 0, 0, 0.2)";
        info_label.style.color = "rgb(75, 75, 75)";
        info_label.style.backgroundColor = "rgba(100, 100, 100, 0.1)";
        lines.forEach(line => {
            line.style.backgroundColor = "rgba(100, 100, 100, 0.1)";
        });
    }
}

set_colors();

// Add a click event listener to the button
change_color_button.addEventListener("click", set_colors);

// Add a focus event listener to the input element
input_equation.addEventListener("focus", function() {
    // Change the CSS property when the element gains focus
    if (colors[current_color_index][2] === "#ffffff") {
        input_equation.style.boxShadow = "0 0 10px 0 rgba(255, 255, 255, 0.75), inset 0 0 10px 0 rgba(255, 255, 255, 0.3)";
    } else {
        input_equation.style.boxShadow = "0 0 10px 0 rgba(0, 0, 0, 0.75), inset 0 0 10px 0 rgba(0, 0, 0, 0.3)";
    }
});

// Add a blur event listener to the input element to reset the style when it loses focus
input_equation.addEventListener("blur", function() {
    // Reset the CSS property when the element loses focus
    if (colors[current_color_index][2] === "#ffffff") {
        input_equation.style.boxShadow = "0 0 10px 0 rgba(255, 255, 255, 0.5), inset 0 0 10px 0 rgba(255, 255, 255, 0.2)";
    } else {
        input_equation.style.boxShadow = "0 0 10px 0 rgba(0, 0, 0, 0.5), inset 0 0 10px 0 rgba(0, 0, 0, 0.2)";
    }
});


// Check if the media query condition is met
if (window.matchMedia("(hover: hover) and (pointer: fine)").matches) {
    balance_button.addEventListener("mouseenter", function() {
        if (colors[current_color_index][2] === "#ffffff") {
            balance_button.style.boxShadow = "0 0 10px 0 rgba(255, 255, 255, 0.9)";
        } else {
            balance_button.style.boxShadow = "0 0 10px 0 rgba(0, 0, 0, 0.9)";
        }
    });
    
    balance_button.addEventListener("mouseleave", function() {
        if (colors[current_color_index][2] === "#ffffff") {
            balance_button.style.boxShadow = "0 0 10px 0 rgba(255, 255, 255, 0.5)";
        } else {
            balance_button.style.boxShadow = "0 0 10px 0 rgba(0, 0, 0, 0.5)";
        }
    });
} else {
    balance_button.addEventListener("mousedown", function() {
        if (colors[current_color_index][2] === "#ffffff") {
            balance_button.style.boxShadow = "0 0 10px 0 rgba(255, 255, 255, 0.9)";
        } else {
            balance_button.style.boxShadow = "0 0 10px 0 rgba(0, 0, 0, 0.9)";
        }
    });
    
    balance_button.addEventListener("mouseup", function() {
        if (colors[current_color_index][2] === "#ffffff") {
            balance_button.style.boxShadow = "0 0 10px 0 rgba(255, 255, 255, 0.5)";
        } else {
            balance_button.style.boxShadow = "0 0 10px 0 rgba(0, 0, 0, 0.5)";
        }
    });
}
