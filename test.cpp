#include <iostream>
#include <vector>
#include <queue>
#include <limits>

using namespace std;

struct EdgePair {
    int v1, v2;
};

int findSmallestDistance(const vector<EdgePair>& all_edge_pairs, const vector<vector<int>>& vec, int start_vertex, int end_vertex) {
    vector<bool> visited(vec.size(), false);
    vector<int> distance(vec.size(), numeric_limits<int>::max());

    queue<int> q;
    q.push(start_vertex);
    visited[start_vertex] = true;
    distance[start_vertex] = 0;

    while (!q.empty()) {
        int current_vertex = q.front();
        q.pop();

        for (int neighbor : vec[current_vertex]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                distance[neighbor] = distance[current_vertex] + 1;
                q.push(neighbor);

                if (neighbor == end_vertex) {
                    return distance[neighbor];
                }
            }
        }
    }

    return -1; // If v1 and v2 are not connected
}

int main() {
    // Assuming you have already read the graph data from the input file
    int num_vertices, num_edges;
    cin >> num_vertices >> num_edges;

    vector<EdgePair> all_edge_pairs(num_edges);
    vector<vector<int>> vec(num_vertices);

    // Read the edge pairs and populate the adjacency list
    for (int i = 0; i < num_edges; i++) {
        int v1, v2;
        cin >> v1 >> v2;
        all_edge_pairs[i].v1 = v1;
        all_edge_pairs[i].v2 = v2;
        vec[v1].push_back(v2);
        vec[v2].push_back(v1);
    }

    int v1, v2;
    cin >> v1 >> v2;

    int smallest_distance = findSmallestDistance(all_edge_pairs, vec, v1, v2);
    cout << "Smallest distance between v1 and v2: " << smallest_distance << endl;

    return 0;
}
