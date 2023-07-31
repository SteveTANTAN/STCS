#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <set>

using namespace std;

struct Edge {
    int from;
    int to;
    int signed_value;

    Edge(int f, int t, int s) : from(f), to(t), signed_value(s) {}
};

class Graph {
public:
    Graph(const string& filename) {
        ifstream infile(filename);
        if (!infile) {
            cerr << "Error: Unable to open the file " << filename << endl;
            exit(1);
        }

        infile >> num_vertices;

        int num_edges;
        infile >> num_edges;

        int from, to, signed_value;
        for (int i = 0; i < num_edges; ++i) {
            infile >> from >> to >> signed_value;
            adjList[from].emplace_back(from, to, signed_value);
        }

        infile.close();
    }

    const vector<Edge>& getNeighbors(int vertex) const {
        return adjList.at(vertex);
    }

    const unordered_map<int, vector<Edge>>& getAdjList() const {
        return adjList;
    }

private:
    int num_vertices;
    unordered_map<int, vector<Edge>> adjList;
};

vector<Edge> findInducedSubgraph(const Graph& graph, const set<int>& vertices) {
    vector<Edge> induced_subgraph;

    for (int vertex : vertices) {
        const auto& adj_vertices = graph.getAdjList();
        if (adj_vertices.find(vertex) != adj_vertices.end()) {
            for (const Edge& edge : adj_vertices.at(vertex)) {
                if (vertices.find(edge.to) != vertices.end()) {
                    induced_subgraph.push_back(edge);
                }
            }
        }
    }

    return induced_subgraph;
}

void writeEdgesToFile(const string& filename, const vector<Edge>& edges) {
    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Error: Unable to open the file " << filename << endl;
        exit(1);
    }

    for (const Edge& edge : edges) {
        outfile << edge.from << " " << edge.to << " " << edge.signed_value << endl;
    }

    outfile.close();
}

int main() {
    // Replace "temp1.txt" with the correct filename containing the graph information
    Graph graph("temp1.txt");

    set<int> vertices_set = {66, 82, 599, 818, 860, 936, 1012, 1029, 2442, 2443, 2458, 7777};
    
    vector<Edge> induced_subgraph = findInducedSubgraph(graph, vertices_set);

    writeEdgesToFile("output.txt", induced_subgraph);

    return 0;
}
