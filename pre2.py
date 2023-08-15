import collections

def process_input(input_lines):
    excluded_lines = set()
    output_lines = []
    vertex_counter = collections.Counter()
    edge_set = set()
    increment_vertices = False # Use this flag to check if vertices need to be incremented

    for line in input_lines:
        v1, v2, operator = line.strip().split()
        v1, v2 = int(v1), int(v2)  # Ensure that the vertices are integers

        # Check if the reverse order is already in the excluded lines
        if (v2, v1) in edge_set or (v1, v2) in edge_set:
            continue
        if v1 == v2:
            continue

        # Add the current line to the excluded lines and vertex counter
        vertex_counter.update([v1, v2])
        edge_set.add((v1, v2))

        if int(operator) > 0:
            output_lines.append(f"{v1}\t{v2}\t{1}\n")
        else:
            output_lines.append(f"{v1}\t{v2}\t{-1}\n")

    # Check if vertex 0 exists, if yes, increment all vertices by 1
    if 0 in vertex_counter:
        increment_vertices = True

    if increment_vertices:
        # Update all vertices in the output lines
        updated_output_lines = []
        for line in output_lines:
            v1, v2, operator = line.strip().split()
            updated_output_lines.append(f"{int(v1) + 1}\t{int(v2) + 1}\t{operator}\n")
        output_lines = updated_output_lines
        edge_set = {(x + 1, y + 1) for x, y in edge_set}

    # Calculate max_vertex and num_edges
    max_vertex = max(vertex_counter) + (1 if increment_vertices else 0)
    num_edges = len(edge_set)

    # Add max_vertex and num_edges at the beginning of the output
    output_lines.insert(0, f'{max_vertex} {num_edges}\n')

    return output_lines


# Example usage
filename = "data/sl.txt"
outputname = 'expf/sl/sl.txt'

with open(filename, "r") as file:
    input_lines = file.readlines()

output_lines = process_input(input_lines)

# Save the non-repeated lines back to the original file
with open(outputname, "w") as file:
    file.writelines(output_lines)
