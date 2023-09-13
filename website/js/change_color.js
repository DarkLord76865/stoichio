// Get references to the button and the target elements
const change_color_button = document.getElementById("change_color_button");
const change_color_button_next_color = document.getElementById("change_color_button_next_color");
const header = document.getElementById("header_id");
const body = document.getElementById("body_id");

// Define an array of colors to cycle through, first color is for the header, second color is for the body
const colors = [
    ["#007aff", "#f2f8ff"],
    ["#7a00ff", "#eeeeff"]
];

// Initialize a color index (color at index 0 is assumed to be the starting color)
let current_color_index = 0;
// Try to get the color index from local storage (if there isn't one, current_color_index will stay at 0)
const stored_index = localStorage.getItem("color_index");
if (stored_index) {
    current_color_index = parseInt(stored_index);
}

// Set initial color
header.style.backgroundColor = colors[current_color_index][0];
body.style.backgroundColor = colors[current_color_index][1];
change_color_button_next_color.style.backgroundColor = colors[(current_color_index + 1) % colors.length][0];

// Add a click event listener to the button
change_color_button.addEventListener("click", () => {
    // Increment the color index, looping back to the start if necessary
    current_color_index = (current_color_index + 1) % colors.length;

    // Store the color index in local storage
    localStorage.setItem("color_index", current_color_index);

    // Change the background color
    header.style.backgroundColor = colors[current_color_index][0];
    body.style.backgroundColor = colors[current_color_index][1];
    change_color_button_next_color.style.backgroundColor = colors[(current_color_index + 1) % colors.length][0];
});
