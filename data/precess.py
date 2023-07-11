def process_input(input_lines):
    excluded_lines = set()
    output_lines = []
    num_vertices = set()
    num_edges = 0
    for line in input_lines:
        v1, v2, operator = line.strip().split('\t')

        # Check if the reverse order is already in the excluded lines
        reverse_line = f"{v2}\t{v1}"
        if reverse_line in excluded_lines or f"{v1}\t{v2}" in excluded_lines:
            continue

        # Process the current line
        print(f"Processing: {line}")

        # Add the current line to the excluded lines and output lines
        # Increment vertex count
        num_vertices.add(v1)
        num_vertices.add(v2)
        num_edges += 1

        excluded_lines.add(f"{v1}\t{v2}")
        output_lines.append(line)

    print("Processing complete.")
    print(str(len(num_vertices)))
    print(str(num_edges))
    return output_lines


# Example usage
filename = "temp.txt"
outputname = 'temp1.txt'

with open(filename, "r") as file:
    input_lines = file.readlines()

output_lines = process_input(input_lines)

# Save the non-repeated lines back to the original file
with open(outputname, "w") as file:
    file.writelines(output_lines)

    
