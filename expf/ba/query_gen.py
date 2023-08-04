import random
from collections import defaultdict

# 读取文件并构造图
with open('ba.txt', 'r') as file:
    lines = file.read().splitlines()
    num_vertex, num_edge = map(int, lines[0].split())
    graph = defaultdict(int)
    for line in lines[1:]:
        from_vertex, to_vertex, _ = map(int, line.split())
        graph[from_vertex] += 1
        graph[to_vertex] += 1

# 将顶点按度数排序
sorted_vertices = sorted(graph.items(), key=lambda x: x[1], reverse=True)

# 计算每个阶段的顶点数量
percentiles = [int(len(sorted_vertices)*perc/100) for perc in [20, 40, 60, 80, 100]]

# 为每个阶段选择100个随机顶点并写入文件
for i, perc in enumerate(percentiles):
    selected_vertices = random.sample(sorted_vertices[:perc], 100)
    selected_vertices = [vertex for vertex, degree in selected_vertices]
    with open(f'ba_{(i+1)*20}.txt', 'w') as file:
        for vertex in selected_vertices:
            file.write(str(vertex) + '\n')
