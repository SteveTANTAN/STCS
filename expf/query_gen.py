import random
from collections import defaultdict
from math import ceil

# Read file and construct graph
try:
    with open('sl/sl.txt', 'r') as file:
        lines = file.read().splitlines()
        if len(lines) < 2:
            print("File has insufficient data.")
            exit()

        num_vertex, num_edge = map(int, lines[0].split())
        graph = defaultdict(int)
        for line in lines[1:]:
            from_vertex, to_vertex, _ = map(int, line.split())
            graph[from_vertex] += 1
            graph[to_vertex] += 1
except FileNotFoundError:
    print("File not found.")
    exit()
except Exception as e:
    print(f"An error occurred: {e}")
    exit()

# Sort vertices by degree
sorted_vertices = sorted(graph.items(), key=lambda x: x[1], reverse=True)

# Partition into 5 equal-sized buckets
bucket_size = ceil(len(sorted_vertices) / 5)

for i in range(5):
    # Extract vertices for this bucket
    start_idx = i * bucket_size
    end_idx = min((i + 1) * bucket_size, len(sorted_vertices))
    bucket = sorted_vertices[start_idx:end_idx]

    # Randomly select 100 vertices from this bucket
    selected_count = min(100, len(bucket))
    selected_vertices = random.sample(bucket, selected_count)
    selected_vertices = [vertex for vertex, degree in selected_vertices]

    # Write to file
    with open(f'sl/Q/sl_{(i + 1) * 20}.txt', 'w') as file:
        for vertex in selected_vertices:
            file.write(str(vertex) + '\n')
