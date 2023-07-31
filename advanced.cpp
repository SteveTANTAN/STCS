#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<string>
#include<vector>
#include<queue>
#include<fstream>
//#include<Windows.h>
#include<cstring>
#include <set>
#include <limits>
#include<sys/time.h>

using namespace std;
const int MAX_V = 1160000;
const int MAX_E = 3000000;
int e_num, v_num , nhop_e_num;

class Edge {
public:
	int v1;
	int v2;
	int sign;
	Edge() :v1(-1), v2(-1), sign(0) {}
};


class Triangle {
public:
	int edge1;
	int edge2;
	int edge3;
	int v1;
	int v2;
	int v3;
	int is_balanced;
	int is_broken;
	Triangle() :edge1(-1), edge2(-1), edge3(-1), v1(-1), v2(-1), v3(-1), is_balanced(0), is_broken(0) {}
};

class Group {
public:
	vector<int> edges;
	int is_delete;
	int typ_edge;
	Group():is_delete(0),typ_edge(-1){}
};

typedef struct {
	Edge all_edge_pairs[MAX_E];
	vector<int> vec[MAX_V];
	vector<int> adj_edge[MAX_V];
	//vector<int> edge_to_edge[MAX_E];
	vector<int> in_which_triangle[MAX_E];
	vector<int> is_linked;
	int book[MAX_V];
	vector<int> is_booked;
	int is_delete_e[MAX_V];
	int is_delete_vec[MAX_E];
	int link[MAX_E];
	int temp_delete_e[MAX_E];
	int break_unb[MAX_E];
	int support[MAX_E];
	int edge_num;


	int unbalance_num = 0;
	vector<Triangle> Triangles;
	int size_of_truss;
	vector<int> followers[MAX_E + 1];
	vector<int> two_dimension;
	int grouped[MAX_E];
	int candidate[MAX_E];
	vector<int> candidates;
	vector<Group> groups;
	vector<int> id;

	int diameter;
	vector<int> path;
} Graph;




clock_t start, finish;
int start_vertex, k, global_hop;


// modify from here
int hop_num[MAX_V];
int visited[MAX_V];

string filename, outname;



Graph* build_graph() {
	Graph *g = new Graph();
	memset(g->is_delete_e, 0, sizeof(g->is_delete_e));
	memset(g->is_delete_vec, 0, sizeof(g->is_delete_vec));
	memset(g->support, 0, sizeof(g->support));
	memset(g->book, 0, sizeof(g->book));
	memset(g->break_unb, 0, sizeof(g->break_unb));
	memset(g->temp_delete_e, 0, sizeof(g->temp_delete_e));
	memset(g->link, 0, sizeof(g->link));
	memset(g->grouped, -1, sizeof(g->grouped));
	memset(g->candidate, 0, sizeof(g->candidate));
	g->two_dimension.resize(MAX_E);
	ifstream inputFile{ filename };
    
    //ifstream inputFile("data/temp1.txt");
    if (!inputFile) {
        cerr << "Input file could not be opened!" << endl;
        exit(EXIT_FAILURE);
    }

    //inputFile >> vertexnum >> edgenum;

	inputFile >> v_num >> e_num;
	for (int i = 0; i < e_num; i++)
	{
		int v1, v2, sign;
		inputFile >> v1 >> v2 >> sign;
		g->all_edge_pairs[i].v1 = v1;
		g->all_edge_pairs[i].v2 = v2;
		g->all_edge_pairs[i].sign = sign;
		g->vec[v1].push_back(v2);
		g->vec[v2].push_back(v1);
		g->adj_edge[v1].push_back(i);
		g->adj_edge[v2].push_back(i);
	}
	for (int i = 0; i < MAX_E; i++) {
		g->followers[MAX_E].push_back(0);
	}
	g->diameter = e_num;
	g->edge_num = e_num;

	return g;
}
bool if_query_inside(Graph* g_ori) {
	for (auto i : g_ori->adj_edge[start_vertex]) {
		if (!g_ori->is_delete_e[i]){
			return true;
		}
	}
	return false;
}

Graph* GetKtrusswith_Nhops(int n, int k, Graph* g_ori) {
	

	Graph *g = new Graph();
	int counts = 0;
	//////////   No subgraph build Method /////////////////////////
	// *g = *g_ori;
	// g->size_of_truss = e_num;
	//////////////////////////////////////////////////////////////


	///////////// build subgraph ///////////////////////////
	memset(g->is_delete_e, 0, sizeof(g->is_delete_e));
	memset(g->is_delete_vec, 0, sizeof(g->is_delete_vec));
	memset(g->support, 0, sizeof(g->support));
	memset(g->book, 0, sizeof(g->book));
	memset(g->break_unb, 0, sizeof(g->break_unb));
	memset(g->temp_delete_e, 0, sizeof(g->temp_delete_e));
	memset(g->link, 0, sizeof(g->link));
	memset(g->grouped, -1, sizeof(g->grouped));
	memset(g->candidate, 0, sizeof(g->candidate));
	g->two_dimension.resize(MAX_E);

	nhop_e_num = 0;

	for (int i = 0; i < g_ori->edge_num; i++)
	{
		int v1 = g_ori->all_edge_pairs[i].v1;
		int v2 = g_ori->all_edge_pairs[i].v2;
		if (hop_num[v1] <= n && hop_num[v2] <= n
		&& hop_num[v1] != -1 && hop_num[v2] != -1) {
			g->all_edge_pairs[nhop_e_num].v1 = v1;
			g->all_edge_pairs[nhop_e_num].v2 = v2;
			g->all_edge_pairs[nhop_e_num].sign = g_ori->all_edge_pairs[i].sign;
			g->vec[v1].push_back(v2);
			g->vec[v2].push_back(v1);
			g->adj_edge[v1].push_back(nhop_e_num);
			g->adj_edge[v2].push_back(nhop_e_num);
			// g->is_delete_e[nhop_e_num] = 0;
			nhop_e_num ++;
		}
	}
	for (int i = 0; i < MAX_E; i++) {
		g->followers[i].push_back(0);
	}
	g->diameter = nhop_e_num;
	g->edge_num = nhop_e_num;
	g->size_of_truss = nhop_e_num;

	for (int i = 0; i < g->edge_num; i++) {
		// If the vertices of the edge are not within n hops of the source, skip this iteration
		// if (hop_num[all_edge_pairs[i].v1] > n || hop_num[all_edge_pairs[i].v2] > n) 
    	// 	continue;
		Triangle temp_triangle;
		// Assigning values to each edge
		temp_triangle.v1 = g->all_edge_pairs[i].v1;
		temp_triangle.v2 = g->all_edge_pairs[i].v2;
		

		temp_triangle.edge1 = i;
		for (int j = 0; j < g->vec[g->all_edge_pairs[i].v1].size(); j++) {
			int v = g->vec[g->all_edge_pairs[i].v1][j];
			// If the vertex v is not within n hops of the source, skip this iteration
			// if (hop_num[v] > n)
			// 	continue;
			
			counts++;
			g->book[v] = 1;
			g->two_dimension[v] = j;
			g->is_booked.push_back(v);
		}
		for (int j = 0; j < g->vec[g->all_edge_pairs[i].v2].size(); j++) {
			int v = g->vec[g->all_edge_pairs[i].v2][j];
			int edg2 = g->adj_edge[g->all_edge_pairs[i].v2][j];
			// If the vertex v is not within n hops of the source, skip this iteration
			// if (hop_num[v] > n)
			// 	continue;
			counts++;
			if (g->book[v]) {
				int edg1 = g->adj_edge[g->all_edge_pairs[i].v1][g->two_dimension[v]]; //?edg1
				if (edg1 > i && edg2 > i) {
					temp_triangle.edge2 = edg1;
					temp_triangle.edge3 = edg2;
					temp_triangle.v3 = v;
					
					if ((g->all_edge_pairs[edg1].sign + g->all_edge_pairs[edg2].sign + g->all_edge_pairs[i].sign) == -1) {
						temp_triangle.is_balanced = 1;
					}
					else if ((g->all_edge_pairs[edg1].sign + g->all_edge_pairs[edg2].sign + g->all_edge_pairs[i].sign) == 3) {
						temp_triangle.is_balanced = 1;
					}
				
					g->Triangles.push_back(temp_triangle);
					g->in_which_triangle[i].push_back(g->Triangles.size() - 1);
					g->in_which_triangle[edg1].push_back(g->Triangles.size() - 1);
					g->in_which_triangle[edg2].push_back(g->Triangles.size() - 1);
				}
			}
		}
		for (int j = 0; j < g->is_booked.size(); j++) {
			g->book[g->is_booked[j]] = 0;
		}
		g->is_booked.clear();
	}
	for (int i = 0; i < g->Triangles.size(); i++) {
		if (!g->Triangles[i].is_balanced)
			g->unbalance_num++;
		else {
			g->support[g->Triangles[i].edge1]++;
			g->support[g->Triangles[i].edge2]++;
			g->support[g->Triangles[i].edge3]++;
		}
	}
	cout << "counts:" << counts << endl;
	//unbalance_num = Triangles.size() - num_of_balance;
	cout << "orangial:" << g->unbalance_num << endl;
	queue<int> q;
	// 找出所有不满足 support的边
	for (int i = 0; i < g->edge_num; i++) {
		if (!g->is_delete_e[i] ) {
			if (g->support[i] < k - 2|| hop_num[g->all_edge_pairs[i].v1] > n || hop_num[g->all_edge_pairs[i].v2] > n
		|| hop_num[g->all_edge_pairs[i].v1] == -1 || hop_num[g->all_edge_pairs[i].v2] == -1) {			
			// if (g->support[i] < k - 2) {
				g->is_delete_e[i] = 1;
				g->size_of_truss--;
				g->temp_delete_e[i] = 1;
				q.push(i);
			}
		}
	}
	while (!q.empty()) {
		int sub = q.front();
		q.pop();
		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
		for (int i = 0; i < g->in_which_triangle[sub].size(); i++) {
			if (!g->Triangles[g->in_which_triangle[sub][i]].is_broken) {
				if (g->Triangles[g->in_which_triangle[sub][i]].is_balanced) {
					// 删除临边
					g->support[g->Triangles[g->in_which_triangle[sub][i]].edge1]--;
					g->support[g->Triangles[g->in_which_triangle[sub][i]].edge2]--;
					g->support[g->Triangles[g->in_which_triangle[sub][i]].edge3]--;
				}
				else 
					g->unbalance_num--;
				// 删除一条边后check 他的临边
				if (!g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge1]) {
					if (g->support[g->Triangles[g->in_which_triangle[sub][i]].edge1] < k - 2)
					{
						q.push(g->Triangles[g->in_which_triangle[sub][i]].edge1);
						g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge1] = 1;
						g->temp_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge1] = 1;
						g->size_of_truss--;
					}
				}
				if (!g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge2]) {
					if (g->support[g->Triangles[g->in_which_triangle[sub][i]].edge2] < k - 2)
					{
						q.push(g->Triangles[g->in_which_triangle[sub][i]].edge2);
						g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge2] = 1;
						g->temp_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge2] = 1;
						g->size_of_truss--;
					}
				}
				if (!g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge3]) {
					if (g->support[g->Triangles[g->in_which_triangle[sub][i]].edge3] < k - 2)
					{
						q.push(g->Triangles[g->in_which_triangle[sub][i]].edge3);
						g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge3] = 1;
						g->temp_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge3] = 1;
						g->size_of_truss--;
					}
				}
				g->Triangles[g->in_which_triangle[sub][i]].is_broken = 1;
			}
		}
	}
	cout << "size of KTruss" << g->size_of_truss << endl;
	cout << "truss unb num:" << g->unbalance_num << endl;

	return g;
	// return g->size_of_truss > 0;
	/*
	if (size_of_truss > 0){
		return true;
	}
	else {
		return false;
	}*/



}







vector<int> findLongestPath(Graph *g) {
	// little optimization by replace v_num by v_num-1
	//int dist[MAX_V][MAX_V] = {MAX_E + 1};
	// cout << "in\n";
	int pathLength = -1;
	vector<int> path;
	
	for (int i = 1; i < v_num - 1; i++) {
		memset(visited, -1, sizeof(visited));
		queue<pair<int,pair<int, vector<int>>>> q;
		vector<int> curr_path;
		curr_path.push_back(i);
		q.push(make_pair(i, make_pair(0, curr_path)));
		
		while (!q.empty()) {
			int size = q.size();
			for (size; size > 0; size--) {
				auto curr_pair = q.front();
				q.pop();
				// cout << "current node is " << curr_pair.first << "\n";
				// cout << "current path size is " << curr_pair.second.first << "\n";
				// cout << "current path is :\n";
				// for (auto v : curr_pair.second.second) {
				// 	cout << v << "->";
				// }
				// cout << "\n";
				if (visited[curr_pair.first] <= curr_pair.second.first && visited[curr_pair.first] != -1) continue;

				visited[curr_pair.first] = curr_pair.second.first;
				if (pathLength < curr_pair.second.first) {
					path = curr_pair.second.second;
					pathLength = curr_pair.second.first;
				}
				//cout << "neighour number is " << g->vec[4].size() << "\n";
				for (int idx = 0; idx < g->adj_edge[curr_pair.first].size(); idx++) {

					if (g->is_delete_e[g->adj_edge[curr_pair.first][idx]]) {
						continue;
					}
					int v;
					if (g->all_edge_pairs[g->adj_edge[curr_pair.first][idx]].v1 == curr_pair.first) {
						v = g->all_edge_pairs[g->adj_edge[curr_pair.first][idx]].v2;
					} else {
						v = g->all_edge_pairs[g->adj_edge[curr_pair.first][idx]].v1;
					}
					// cout << "BFS to " << g->vec[curr_pair.first][idx] << "\n";
					// cout << " newPath is :\n";
					vector<int> nextPath = curr_pair.second.second;
					nextPath.push_back(v);
					// for (auto v : nextPath) {
					// 	cout << v << "->";
					// }
					// cout << "\n";

					q.push(make_pair(v, make_pair(curr_pair.second.first + 1, nextPath)));
					
				}
			}
			
		}
		
	}

	// cout << "================= current longest path is:\n";
	// for (auto v : path) {
	// 	cout << v << "->";
	// }

	//scout << "\n";
	g->diameter = pathLength;
	g->path = path;
	
	return path;
}
vector<int> findLongestDistanceFromStartVertex(Graph *g) {
    vector<int> distances(v_num, -1);
    vector<vector<int>> paths(v_num);
    vector<bool> visited(v_num, false);
    queue<int> q;

    distances[start_vertex] = 0;
    paths[start_vertex].push_back(start_vertex);
    visited[start_vertex] = true;
    q.push(start_vertex);

    int max_distance = 0;
    vector<int> longest_path;
	vector<int> longest_list;

    while (!q.empty()) {
        int curr_vertex = q.front();
        q.pop();

        if(distances[curr_vertex] > max_distance){
            max_distance = distances[curr_vertex];
            longest_path = paths[curr_vertex];
        }

        for (int idx = 0; idx < g->adj_edge[curr_vertex].size(); idx++) {
            if (g->is_delete_e[g->adj_edge[curr_vertex][idx]]) {
                continue;
            }

            int v;
            if (g->all_edge_pairs[g->adj_edge[curr_vertex][idx]].v1 == curr_vertex) {
                v = g->all_edge_pairs[g->adj_edge[curr_vertex][idx]].v2;
            } else {
                v = g->all_edge_pairs[g->adj_edge[curr_vertex][idx]].v1;
            }

            if (!visited[v]) {
                q.push(v);
                visited[v] = true;
                distances[v] = distances[curr_vertex] + 1;
				if (distances[v] > global_hop) longest_list.push_back(v);
                paths[v] = paths[curr_vertex];
                paths[v].push_back(v);
            }
        }
    }

    g->diameter = max_distance;
    g->path = longest_path;

    return longest_list;
}


void print_result(Graph* g){
		// Assuming vertices_set is your set of vertices
	// outname = "baseline_" + filename;
	ofstream outfile(outname);
	// ifstream outfile{ "baseline_result.txt" };
    if (!outfile) {
        cerr << "Error: Unable to open the file " << outname << endl;
        exit(1);
    }
	set<int> vertices_set;
	cout <<  " ====result graph: \n";
	findLongestPath(g);

	// Collect all unique vertices from all_edge_pairs into the set
	for (int i = 0; i < g->edge_num; i++) {
		if (!g->is_delete_e[i]) {
			int v1 = g->all_edge_pairs[i].v1;
			int v2 = g->all_edge_pairs[i].v2;
			int sig = g->all_edge_pairs[i].sign;
			vertices_set.insert(v1);
			vertices_set.insert(v2);
        	outfile << v1 << " " << v2 << " " << sig << endl;
		}
	}
    outfile.close();

	// Print all vertices in the set
	for (const auto& vertex : vertices_set) {
		cout << vertex << " ";
	}
	cout << endl;
	cout<<"final diameter is: "<< g->diameter << endl;
	cout << "final size of KTruss: " << g->size_of_truss << endl;
	cout << "final truss unbalance num:" << g->unbalance_num << endl;
}

bool delete_on_edge(int edge_num, Graph*newG, int update_dia) {


	queue<int> q;


	if (!newG->is_delete_e[edge_num] ) {
		newG->is_delete_e[edge_num] = 1;
		newG->size_of_truss--;
		newG->temp_delete_e[edge_num] = 1;
		q.push(edge_num);
	}


	while (!q.empty()) {
		int sub = q.front();
		q.pop();
		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
		for (int i = 0; i < newG->in_which_triangle[sub].size(); i++) {
			if (!newG->Triangles[newG->in_which_triangle[sub][i]].is_broken) {
				if (newG->Triangles[newG->in_which_triangle[sub][i]].is_balanced) {
					// 删除临边
					newG->support[newG->Triangles[newG->in_which_triangle[sub][i]].edge1]--;
					newG->support[newG->Triangles[newG->in_which_triangle[sub][i]].edge2]--;
					newG->support[newG->Triangles[newG->in_which_triangle[sub][i]].edge3]--;
				}
				else 
					newG->unbalance_num--;
				// 删除一条边后check 他的临边
				if (!newG->is_delete_e[newG->Triangles[newG->in_which_triangle[sub][i]].edge1]) {
					if (newG->support[newG->Triangles[newG->in_which_triangle[sub][i]].edge1] < k - 2)
					{
						q.push(newG->Triangles[newG->in_which_triangle[sub][i]].edge1);
						newG->is_delete_e[newG->Triangles[newG->in_which_triangle[sub][i]].edge1] = 1;
						newG->temp_delete_e[newG->Triangles[newG->in_which_triangle[sub][i]].edge1] = 1;
						newG->size_of_truss--;
					}
				}
				if (!newG->is_delete_e[newG->Triangles[newG->in_which_triangle[sub][i]].edge2]) {
					if (newG->support[newG->Triangles[newG->in_which_triangle[sub][i]].edge2] < k - 2)
					{
						q.push(newG->Triangles[newG->in_which_triangle[sub][i]].edge2);
						newG->is_delete_e[newG->Triangles[newG->in_which_triangle[sub][i]].edge2] = 1;
						newG->temp_delete_e[newG->Triangles[newG->in_which_triangle[sub][i]].edge2] = 1;
						newG->size_of_truss--;
					}
				}
				if (!newG->is_delete_e[newG->Triangles[newG->in_which_triangle[sub][i]].edge3]) {
					if (newG->support[newG->Triangles[newG->in_which_triangle[sub][i]].edge3] < k - 2)
					{
						q.push(newG->Triangles[newG->in_which_triangle[sub][i]].edge3);
						newG->is_delete_e[newG->Triangles[newG->in_which_triangle[sub][i]].edge3] = 1;
						newG->temp_delete_e[newG->Triangles[newG->in_which_triangle[sub][i]].edge3] = 1;
						newG->size_of_truss--;
					}
				}
				newG->Triangles[newG->in_which_triangle[sub][i]].is_broken = 1;
			}
		}
	}

	//cout << newG->size_of_truss <<"abc\n";
	if (newG->size_of_truss > 0 && if_query_inside(newG)) {

		return true;
	} 
	cout << "delete on edge fail, newG->size_of_truss > 0 && if_query_inside(newG)\n";
	return false;
}

bool delete_on_node(int A, Graph*newG) {

	int edge_num = -1;
	queue<int> q;

	for (int i = 0; i < newG->adj_edge[A].size(); ++i) {
        int edge_index = newG->adj_edge[A][i];
        // Check if this edge index appears in the list for v2
		if (!newG->is_delete_e[edge_index] ) {
			if (!delete_on_edge(edge_index,newG, false)) return false;
		}
    }
	return true;

}


// Assuming you have the Graph and related classes defined here...

// Function to check if Graph H1 is a subgraph of Graph H2.
bool isSubgraph(Graph* H1, Graph* H2) {
    // Check if the number of vertices in H1 is greater than H2.
    if (H1->size_of_truss > H2->size_of_truss) {
        return false;
    }

    // Check if all edges of H1 are also present in H2.
    for (int i = 0; i < H1->edge_num; i++) {
        if (!H1->is_delete_e[i] && H2->is_delete_e[i]) {
			return false;
			break;
		}
    }
    return true;
}



bool removeEdgeFromLongestPath(Graph* g) {
	/*
	vector<int> path = findLongestPath(g);
	cout << "Longest Path is :\n";
	for (auto v : path) {
		cout << v << " -> ";
	}
	cout << "\n";
	*/
	//TODO
	//Graph* tempGraph = new Graph();
	//*tempGraph = *g;
	int curr_diameter = g->diameter;
	if (curr_diameter < 1) return false;
	queue<Graph*> Queue;
	Graph *best = new Graph();
	*best = *g;
	Queue.push(best);
	int que_in,que_out = 0;
	while (!Queue.empty()) {
		
		int size = Queue.size();
		Graph *left = Queue.front();
		Queue.pop();
		que_out++;
		
		// cout << "mid1 queue \n";
		Graph *right = new Graph();
		// cout << "mid2 queue\n";
		*right = *left;
		if (delete_on_node(left->path[right->path.size() - 1], left)) {
			bool best = false;
			if ((left->diameter == g->diameter && g->size_of_truss < left->size_of_truss ) ||
			left->diameter < g->diameter) {
				*g = *left;
				best = true;
			}
			if (left->diameter > 2*global_hop) {
			// if (if_query_inside(left) && left->size_of_truss >= g->size_of_truss) {
				Queue.push(left);
				que_in++;
			} else {
				if (!best) delete(left);
			}
		} else {
			delete(left);
		}
		if (delete_on_node(right->path[right->path.size() - 1], right)) {
			bool best = false;
			if ((right->diameter == g->diameter && g->size_of_truss < right->size_of_truss ) ||
			right->diameter < g->diameter) {
				*g = *right;
				best = true;
			}

			if (right->diameter > 2*global_hop) {
				// if (right->diameter >= g->diameter && )
				Queue.push(right);
				que_in++;

			} else {
				if (!best) delete(right);
			}
		} else {
			delete(right);
		}
		queue<Graph*> filteredQueue;
	// Check each subgraph in Queue.
		while (!Queue.empty()) {
			Graph* currentGraph = Queue.front();
			Queue.pop();

			// Check if currentGraph is a subgraph of any other subgraph in filteredQueue.
			bool isSubgraphOfOther = false;
			queue<Graph*> tempQueue = filteredQueue; // Copy filteredQueue to a temporary queue.

			while (!tempQueue.empty()) {
				Graph* graphInFiltered = tempQueue.front();
				tempQueue.pop();

				if (isSubgraph(currentGraph, graphInFiltered)) {
					isSubgraphOfOther = true;
					// cout << "prunning queue\n";
					break;
				}
			}

			// If currentGraph is not a subgraph of any other subgraph in filteredQueue, keep it.
			if (!isSubgraphOfOther) {
				filteredQueue.push(currentGraph);
			} else {
				// If currentGraph is a subgraph of some other subgraph, you may want to free its memory.
				delete currentGraph;
			}
		}

		// Now, filteredQueue contains only the subgraphs that are not subgraphs of any other subgraph.
		// You can continue using filteredQueue or copy its contents back to the original Queue.

		// Clean up the original Queue if needed.
		while (!Queue.empty()) {
			delete Queue.front();
			Queue.pop();
		}

		// Copy the filtered subgraphs back to the original Queue.
		Queue = filteredQueue;

		    // Remember to free memory of the Graph objects properly.
		// while (!Queue.empty()) {
		// 	delete Queue.front();
		// 	Queue.pop();
		// }

	}

	cout<< "que_in" << que_in<<"\n";
	cout<< "que_out" << que_out<<"\n";
	return true;
}


bool removeNegativeTriangle(Graph* g) {
	int curr_diameter = g->diameter;
	if (curr_diameter < 1) return false;



    for (auto unb : g->Triangles) {
        if (!unb.is_balanced && !unb.is_broken) {
            if (hop_num[g->all_edge_pairs[unb.edge1].v1] = hop_num[g->all_edge_pairs[unb.edge1].v2]) {
                if (!delete_on_edge(unb.edge1, g,true)) return false;
            } else if (hop_num[g->all_edge_pairs[unb.edge2].v1] = hop_num[g->all_edge_pairs[unb.edge2].v2]) {
                if (!delete_on_edge(unb.edge1, g,true)) return false;
            } else {
                if (!delete_on_edge(unb.edge1, g,true)) return false;
            }
        }
    }
	return true;
}



void GetmaximumKtruss(Graph *g) {
	g->size_of_truss = e_num;
	int counts = 0;
	// for (int i = 0; i < e_num; i++) {
	// 	if (hop_num[all_edge_pairs[i].v1] <= n || hop_num[all_edge_pairs[i].v2] <= n) size_of_truss++;
	// }
	for (int i = 0; i < e_num; i++) {
		// If the vertices of the edge are not within n hops of the source, skip this iteration
		// if (hop_num[all_edge_pairs[i].v1] > n || hop_num[all_edge_pairs[i].v2] > n) 
    	// 	continue;
		Triangle temp_triangle;
		// Assigning values to each edge
		temp_triangle.v1 = g->all_edge_pairs[i].v1;
		temp_triangle.v2 = g->all_edge_pairs[i].v2;
		

		temp_triangle.edge1 = i;
		for (int j = 0; j < g->vec[g->all_edge_pairs[i].v1].size(); j++) {
			int v = g->vec[g->all_edge_pairs[i].v1][j];
			// If the vertex v is not within n hops of the source, skip this iteration
			// if (hop_num[v] > n)
			// 	continue;
			
			counts++;
			g->book[v] = 1;
			g->two_dimension[v] = j;
			g->is_booked.push_back(v);
		}
		for (int j = 0; j < g->vec[g->all_edge_pairs[i].v2].size(); j++) {
			int v = g->vec[g->all_edge_pairs[i].v2][j];
			int edg2 = g->adj_edge[g->all_edge_pairs[i].v2][j];
			// If the vertex v is not within n hops of the source, skip this iteration
			// if (hop_num[v] > n)
			// 	continue;
			counts++;
			if (g->book[v]) {
				int edg1 = g->adj_edge[g->all_edge_pairs[i].v1][g->two_dimension[v]]; //?edg1
				if (edg1 > i && edg2 > i) {
					temp_triangle.edge2 = edg1;
					temp_triangle.edge3 = edg2;
					temp_triangle.v3 = v;
					
					if ((g->all_edge_pairs[edg1].sign + g->all_edge_pairs[edg2].sign + g->all_edge_pairs[i].sign) == -1) {
						temp_triangle.is_balanced = 1;
					}
					else if ((g->all_edge_pairs[edg1].sign + g->all_edge_pairs[edg2].sign + g->all_edge_pairs[i].sign) == 3) {
						temp_triangle.is_balanced = 1;
					}
				
					g->Triangles.push_back(temp_triangle);
					g->in_which_triangle[i].push_back(g->Triangles.size() - 1);
					g->in_which_triangle[edg1].push_back(g->Triangles.size() - 1);
					g->in_which_triangle[edg2].push_back(g->Triangles.size() - 1);
				}
			}
		}
		for (int j = 0; j < g->is_booked.size(); j++) {
			g->book[g->is_booked[j]] = 0;
		}
		g->is_booked.clear();
	}
	for (int i = 0; i < g->Triangles.size(); i++) {
		if (!g->Triangles[i].is_balanced)
			g->unbalance_num++;
		else {
			g->support[g->Triangles[i].edge1]++;
			g->support[g->Triangles[i].edge2]++;
			g->support[g->Triangles[i].edge3]++;
		}
	}
	cout << "counts:" << counts << endl;
	//unbalance_num = Triangles.size() - num_of_balance;
	cout << "orangial:" << g->unbalance_num << endl;
	queue<int> q;
	// 找出所有不满足 support的边
	for (int i = 0; i < e_num; i++) {
		if (!g->is_delete_e[i] ) {
			if (g->support[i] < k - 2) {
				g->is_delete_e[i] = 1;
				g->size_of_truss--;
				g->temp_delete_e[i] = 1;
				q.push(i);
			}
		}
	}
	while (!q.empty()) {
		int sub = q.front();
		q.pop();
		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
		for (int i = 0; i < g->in_which_triangle[sub].size(); i++) {
			if (!g->Triangles[g->in_which_triangle[sub][i]].is_broken) {
				if (g->Triangles[g->in_which_triangle[sub][i]].is_balanced) {
					// 删除临边
					g->support[g->Triangles[g->in_which_triangle[sub][i]].edge1]--;
					g->support[g->Triangles[g->in_which_triangle[sub][i]].edge2]--;
					g->support[g->Triangles[g->in_which_triangle[sub][i]].edge3]--;
				}
				else 
					g->unbalance_num--;
				// 删除一条边后check 他的临边
				if (!g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge1]) {
					if (g->support[g->Triangles[g->in_which_triangle[sub][i]].edge1] < k - 2)
					{
						q.push(g->Triangles[g->in_which_triangle[sub][i]].edge1);
						g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge1] = 1;
						g->temp_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge1] = 1;
						g->size_of_truss--;
					}
				}
				if (!g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge2]) {
					if (g->support[g->Triangles[g->in_which_triangle[sub][i]].edge2] < k - 2)
					{
						q.push(g->Triangles[g->in_which_triangle[sub][i]].edge2);
						g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge2] = 1;
						g->temp_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge2] = 1;
						g->size_of_truss--;
					}
				}
				if (!g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge3]) {
					if (g->support[g->Triangles[g->in_which_triangle[sub][i]].edge3] < k - 2)
					{
						q.push(g->Triangles[g->in_which_triangle[sub][i]].edge3);
						g->is_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge3] = 1;
						g->temp_delete_e[g->Triangles[g->in_which_triangle[sub][i]].edge3] = 1;
						g->size_of_truss--;
					}
				}
				g->Triangles[g->in_which_triangle[sub][i]].is_broken = 1;
			}
		}
	}
	cout << "size of KTruss" << g->size_of_truss << endl;
	cout << "truss unb num:" << g->unbalance_num << endl;






}




bool GetKtruss(int src, int k, Graph* g) {

	// check if the original graph has k-truss including src

	Graph *new_g = new Graph();
	*new_g = * g;
	GetmaximumKtruss(new_g);
	if (!if_query_inside(new_g)) {
		cout << "Cannot obtain any result \n";
		exit(1);
	}
	// print_result(new_g);
	delete(new_g);




	// record all nodes with k hops, hop is initially set to 2
	memset(hop_num, -1, sizeof(hop_num));
	queue<pair<int,int>> q;
	int hop = 2;
	int max_hop = -1;
	q.push(make_pair(src,0));
	hop_num[src] = 0;
	
	while (!q.empty()) {
		int curr_v= q.front().first;
		int curr_hop = q.front().second;
		q.pop();
		for (int i = 0; i < g->vec[curr_v].size(); i++) {
			int v = g->vec[curr_v][i];
			if (hop_num[v] == -1 || hop_num[v] > curr_hop + 1) {
				hop_num[v] = curr_hop + 1;
				if (max_hop < curr_hop + 1) max_hop = curr_hop + 1;
				q.push(make_pair(v, curr_hop + 1));
			}
		}
	}
	cout << "max hop is " << max_hop << "\n";

	// add one hop at a time
	
	bool isSuccess = false;

	int curr_hop = 1;

	while (curr_hop >= 1 && curr_hop <= max_hop) {
		cout << "--------when hop is " << curr_hop << "\n";

		Graph* g_hop =  GetKtrusswith_Nhops(curr_hop, k, g);
		
		if (g_hop->size_of_truss > 0 && if_query_inside(g_hop)) {

			global_hop = curr_hop;
			Graph *newG = new Graph();
			*newG = *g_hop;
			cout << "===================Deleting NegativeTriangle ==================\n";
			if (removeNegativeTriangle(newG)) {
				vector<int> longest_list = findLongestDistanceFromStartVertex(newG);
				if (newG->diameter <= (curr_hop)) {
					*g = *newG;
					return true;
				} else {
					cout<<"---start to delete node \n"<<endl;
					for (auto i : longest_list) {
						delete_on_node(i,newG);
					}
					vector<int> longest_list = findLongestDistanceFromStartVertex(newG);
					if (newG->diameter <= (curr_hop)) {
						*g = *newG;
						return true;
					} 
				}
            }
		} 
		curr_hop += 1;
		delete(g_hop);
	}
	cout<<"---calculation fail! \n"<<endl;
	return false;
}








int main() {


    filename = "data/temp1";
	outname = filename + "_solution.txt";
	filename += ".txt";
	
	Graph* g = build_graph();
	cout << "graph build" << endl;


    set<int> C;
    
    cout << "Enter the start vertex: ";
    cin >> start_vertex;
	cout << "Enter the K: ";
    cin >> k;
	struct timeval start, end;
    double timeuse;
	bool result;
    gettimeofday(&start, NULL);
	result = GetKtruss(start_vertex, k, g);
    gettimeofday(&end, NULL);
    timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
   	cout << "baslinerun time: " << timeuse << '\n' << endl;
	if (result) print_result(g);

	return 0;
}
