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
#include <algorithm>

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
bool delete_on_radius(Graph *g_hop);


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
    *g = *g_ori;
	int counts = 0;

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


}







vector<int> findLongestPath(Graph *g) {
	// little optimization by replace v_num by v_num-1
	//int dist[MAX_V][MAX_V] = {MAX_E + 1};
	// cout << "in\n";
	int pathLength = -1;
	vector<int> path;
	
	for (int i = 1; i < v_num - 1; i++) {
		if (g->is_delete_vec[i]) continue; // skip the vertices marked as deleted

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
// vector<int> findLongestDistanceFromStartVertex(Graph *g) {
//     vector<int> distances(MAX_V, -1);
//     vector<vector<int>> paths(MAX_V);
//     vector<bool> visited(MAX_V, false);
//     queue<int> q;

//     distances[start_vertex] = 0;
//     paths[start_vertex].push_back(start_vertex);
//     visited[start_vertex] = true;
//     q.push(start_vertex);

//     int max_distance = 0;
//     vector<int> longest_path;
// 	vector<int> longest_list;

//     while (!q.empty()) {
//         int curr_vertex = q.front();
//         q.pop();

//         if(distances[curr_vertex] > max_distance){
//             max_distance = distances[curr_vertex];
//             longest_path = paths[curr_vertex];
//         }

//         for (int idx = 0; idx < g->adj_edge[curr_vertex].size(); idx++) {
//             if (g->is_delete_e[g->adj_edge[curr_vertex][idx]]) {
//                 continue;
//             }

//             int v;
//             if (g->all_edge_pairs[g->adj_edge[curr_vertex][idx]].v1 == curr_vertex) {
//                 v = g->all_edge_pairs[g->adj_edge[curr_vertex][idx]].v2;
//             } else {
//                 v = g->all_edge_pairs[g->adj_edge[curr_vertex][idx]].v1;
//             }

//             if (!visited[v]) {
//                 q.push(v);
//                 visited[v] = true;
//                 distances[v] = distances[curr_vertex] + 1;
// 				if (distances[v] > global_hop) longest_list.push_back(v);
//                 paths[v] = paths[curr_vertex];
//                 paths[v].push_back(v);
//             }
//         }
//     }

//     g->diameter = max_distance;
//     g->path = longest_path;

//     return longest_list;
// }
vector<int> findLongestDistanceFromStartVertex(Graph *g) {


    int max_distance = 0;
    vector<int> longest_path;
	vector<int> longest_list;
	memset(hop_num, -1, sizeof(hop_num));


	queue<pair<int,int>> q;
	int hop = 2;
	int max_hop = -1;
	q.push(make_pair(start_vertex,0));
	hop_num[start_vertex] = 0;
	
	while (!q.empty()) {
		int curr_v = q.front().first;
		int curr_hop = q.front().second;
		q.pop();
		for (int i = 0; i < g->adj_edge[curr_v].size(); i++) {
			int edgeIndex = g->adj_edge[curr_v][i];
			if (g->is_delete_e[edgeIndex] || g->temp_delete_e[edgeIndex] ) continue;

			int v = (g->all_edge_pairs[edgeIndex].v1 == curr_v)
					? g->all_edge_pairs[edgeIndex].v2
					: g->all_edge_pairs[edgeIndex].v1;

			if (hop_num[v] == -1 || hop_num[v] > curr_hop + 1) {
				hop_num[v] = curr_hop + 1;
				if (max_hop < curr_hop + 1) max_hop = curr_hop + 1;
				if (hop_num[v] > global_hop) longest_list.push_back(v);
				q.push(make_pair(v, curr_hop + 1));
			}
		}
	}
	
    g->diameter = max_hop;

	// cout<<max_hop<<endl;

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
		if (update_dia) findLongestDistanceFromStartVertex(newG);

		return true;
	} 
	cout << "delete on edge fail, newG->size_of_truss > 0 && if_query_inside(newG)\n";
	return false;
}

bool delete_on_node(int A, Graph*newG) {

	int edge_num = -1;
	queue<int> q;
	if(newG->is_delete_vec[A]) return true;
	for (int i = 0; i < newG->adj_edge[A].size(); ++i) {
        int edge_index = newG->adj_edge[A][i];
        // Check if this edge index appears in the list for v2
		if (!newG->is_delete_e[edge_index] ) {
			if (!delete_on_edge(edge_index,newG, false)) return false;
		}
    }
	newG->is_delete_vec[A] = 1;
	if (newG->size_of_truss > 0) {
		// findLongestPath(newG);
		//cout << "delete on edge ended\n";
		return true;
	} 
	return false;

}
bool delete_on_radius(Graph *g_hop) {
	if (g_hop->size_of_truss == 0) return false;
	vector<int> longest_list = findLongestDistanceFromStartVertex(g_hop);

	while (g_hop->diameter > (global_hop)){
		for (auto i : longest_list) {
			if ( !delete_on_node(i,g_hop)) return false;
		}
		longest_list = findLongestDistanceFromStartVertex(g_hop);
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





bool quickremoveNegativeTriangle(Graph* g) {
	int curr_diameter = g->diameter;
	if (curr_diameter < 1) return false;
	if (g->unbalance_num <= 0) return true;

	if (curr_diameter < 1) return false;

    for (auto unb : g->Triangles) {

		if (!unb.is_balanced && !unb.is_broken){
        // if (!unb.is_balanced && !unb.is_broken) {
			bool runned = false;
			// cout<<"check1\n";

            if (!runned && hop_num[g->all_edge_pairs[unb.edge1].v1] == hop_num[g->all_edge_pairs[unb.edge1].v2]) {
				runned = true;
				// cout<<"check2\n";
				Graph *temp_g = new Graph();
				*temp_g = *g;
				if (!delete_on_edge(unb.edge1, g,false)) {
					*g = *temp_g;
					if (!delete_on_edge(unb.edge2, g,false)) {
						*g = *temp_g;
						if (!delete_on_edge(unb.edge3, g,false)) {
							*g = *temp_g;
							delete(temp_g);
							return false;
						}
					}

				};
				delete(temp_g);
                // delete_on_edge(unb.edge1, g,false);

            } 
			
			else if (!runned && hop_num[g->all_edge_pairs[unb.edge2].v1] == hop_num[g->all_edge_pairs[unb.edge2].v2]) {
				runned = true;
				// cout<<"check3\n";

                // delete_on_edge(unb.edge2, g,false);
				Graph *temp_g = new Graph();
				*temp_g = *g;
				if (!delete_on_edge(unb.edge2, g,false)) {
					*g = *temp_g;
					if (!delete_on_edge(unb.edge1, g,false)) {
						*g = *temp_g;
						if (!delete_on_edge(unb.edge3, g,false)) {
							*g = *temp_g;
							delete(temp_g);
							return false;
						}
					}

				};
				delete(temp_g);
            } 

			else if (!runned && hop_num[g->all_edge_pairs[unb.edge3].v1] == hop_num[g->all_edge_pairs[unb.edge3].v2]) {

				runned = true;
				// cout<<"check4\n";
                // delete_on_edge(unb.edge3, g,false);
				Graph *temp_g = new Graph();
				*temp_g = *g;
				if (!delete_on_edge(unb.edge3, g,false)) {
					*g = *temp_g;
					if (!delete_on_edge(unb.edge1, g,false)) {
						*g = *temp_g;
						if (!delete_on_edge(unb.edge2, g,false)) {
							*g = *temp_g;
							delete(temp_g);
							return false;
						}
					}

				};
				delete(temp_g);

            }
			// cout<<"check=5\n";


        }
    }

	
	return delete_on_radius(g);;
}



bool removeNegativeTriangle(Graph* g) {
	int curr_diameter = g->diameter;
	if (curr_diameter < 1) return false;

	queue<Graph*> Queue;
	Graph *best_g = new Graph();
	*best_g = *g;
	// g->diameter ++;
	Queue.push(best_g);
	int que_in,que_out = 0;
    quickremoveNegativeTriangle(g);
	while (!Queue.empty()) {
		
		int size = Queue.size();
		Graph *graph_first = Queue.front();
		Queue.pop();
		que_out++;

		for (auto unb : graph_first->Triangles) {
			if (!unb.is_balanced && !unb.is_broken) {

				// cout << "test 1111\n";
				Graph *left = new Graph();
				*left = *graph_first;

				Graph *middle = new Graph();
				*middle = *graph_first;

				Graph *right = new Graph();
				*right = *graph_first;
				if (delete_on_edge(unb.edge1, left,true)) {
					// bool best = false;
					if ((left->diameter == g->diameter && left->unbalance_num == 0
					&& g->size_of_truss < left->size_of_truss) ||
					left->diameter < g->diameter && left->unbalance_num == 0) {
						*g = *left;
						// best = true;
					}
					if (left->diameter >= global_hop) {
					// if (left->unbalance_num != 0) {
						Queue.push(left);
						que_in++;
					} else {
						delete(left);
						// if (!best) delete(left);
					}
				} else {
					delete(left);
				}
				if (delete_on_edge(unb.edge2, middle,true)) {
					// bool best = false;
					if ((middle->diameter == g->diameter && middle->unbalance_num == 0
					&& g->size_of_truss < middle->size_of_truss ) ||
					middle->diameter < g->diameter && middle->unbalance_num == 0)  {
						*g = *middle;
						// best = true;
					}

					if (middle->diameter >= global_hop) {
					// if (middle->unbalance_num != 0) {
						Queue.push(middle);
						que_in++;

					} else {
						delete(middle);
						// if (!best) delete(middle);
					}
				} else {
					delete(middle);
				}
				if (delete_on_edge(unb.edge3, right,true)) {
					// bool best = false;
					if ((right->diameter == g->diameter && right->unbalance_num == 0
					&& g->size_of_truss < right->size_of_truss ) ||
					right->diameter < g->diameter && right->unbalance_num == 0) {
						*g = *right;
						// best = true;
					}

					if (right->diameter >= global_hop) {
					// if (right->unbalance_num != 0) {
						Queue.push(right);
						que_in++;

					} else {
						delete(right);
						
						// if (!best) delete(right);
					}
				} else {
					delete(right);
				}

                break;


			}
		}
		delete(graph_first);

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
	}
    if (g->unbalance_num >0) return false;
	return true;
}
vector<int> quick_delete(Graph* g, Triangle unb) {
	findLongestDistanceFromStartVertex(g);
	vector<int> deleted_edge;
	if (hop_num[g->all_edge_pairs[unb.edge1].v1] == hop_num[g->all_edge_pairs[unb.edge1].v2]) {
		// if(check_follower(unb.edge1,g).size() == 1){
		deleted_edge.push_back(unb.edge1);
		// }
	} 
	if (hop_num[g->all_edge_pairs[unb.edge2].v1] == hop_num[g->all_edge_pairs[unb.edge2].v2]) {
		deleted_edge.push_back(unb.edge2);

		// if(check_follower(unb.edge2,g).size() == 1){
		// 	deleted_edge.push_back(unb.edge2);
		// }

	} 
	if (hop_num[g->all_edge_pairs[unb.edge3].v1] == hop_num[g->all_edge_pairs[unb.edge3].v2]) {
		deleted_edge.push_back(unb.edge3);

		// if(check_follower(unb.edge3,g).size() == 1){
		// 	deleted_edge.push_back(unb.edge3);
		// }
	}
	return deleted_edge;

}
// void candidate_lemma(Graph* g) {
	
// 	for (int i = 0; i < g->Triangles.size(); i++) {
// 		if (!g->Triangles[i].is_broken) {
// 			if (!g->Triangles[i].is_balanced) {

//                 // vector<int> quick_result = quick_delete(g, g->Triangles[i]);
//                 vector<int> quick_result;
//                 quick_result.push_back(g->Triangles[i].edge1);
//                 quick_result.push_back(g->Triangles[i].edge2);
//                 quick_result.push_back(g->Triangles[i].edge3);
//                 for (auto edge: quick_result) {
//                     if (!g->candidate[edge])
//                     {
//                         g->candidates.push_back(edge);
//                         g->candidate[edge] = 1;
//                     }
//                 }
// 			}
// 		}
// 	}
// 	int num = 0;
// 	for (int i = 0; i < g->candidates.size(); i++) {
// 		Group temp_group;
// 		if (g->support[g->candidates[i]] == k - 2) {
// 			if (g->grouped[g->candidates[i]] == -1) {
// 				queue<int> q;
// 				q.push(g->candidates[i]);
// 				g->grouped[g->candidates[i]] = g->groups.size();
// 				temp_group.edges.push_back(g->candidates[i]);
// 				while (!q.empty()) {
// 					int sub = q.front();
// 					q.pop();
// 					for (int j = 0; j < g->in_which_triangle[sub].size(); j++) {
// 						int temp = g->in_which_triangle[sub][j];
// 						if (!g->Triangles[temp].is_broken) {
// 							if (g->Triangles[temp].is_balanced) {
// 								if (g->support[g->Triangles[temp].edge1] == k - 2 && g->grouped[g->Triangles[temp].edge1] == -1) {
// 									q.push(g->Triangles[temp].edge1);
// 									temp_group.edges.push_back(g->Triangles[temp].edge1);
// 									g->grouped[g->Triangles[temp].edge1] = g->groups.size();

// 								}
// 								if (g->support[g->Triangles[temp].edge2] == k - 2 && g->grouped[g->Triangles[temp].edge2] == -1) {
// 									q.push(g->Triangles[temp].edge2);
// 									temp_group.edges.push_back(g->Triangles[temp].edge2);
// 									g->grouped[g->Triangles[temp].edge2] = g->groups.size();

// 								}
// 								if (g->support[g->Triangles[temp].edge3] == k - 2 && g->grouped[g->Triangles[temp].edge3] == -1) {
// 									q.push(g->Triangles[temp].edge3);
// 									temp_group.edges.push_back(g->Triangles[temp].edge3);
// 									g->grouped[g->Triangles[temp].edge3] = g->groups.size();

// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 		else if (g->support[g->candidates[i]] > k - 2)
// 		{
// 			g->grouped[g->candidates[i]] = g->groups.size();
// 			temp_group.edges.push_back(g->candidates[i]);
// 		}

// 		if (temp_group.edges.size() > 0)
// 		{
// 			int min = MAX_E + 1;
// 			for (int j = 0; j < temp_group.edges.size(); j++)
// 			{
// 				if (temp_group.edges[j] < min)
// 					min = temp_group.edges[j];
// 			}
// 			temp_group.typ_edge = min;
// 			g->groups.push_back(temp_group);
// 		}
// 	}

// }
void candidate_lemma(Graph* g, int version) {
    memset(g->break_unb, 0, sizeof(g->break_unb));
	memset(g->temp_delete_e, 0, sizeof(g->temp_delete_e));
	memset(g->link, 0, sizeof(g->link));
	memset(g->grouped, -1, sizeof(g->grouped));
	memset(g->candidate, 0, sizeof(g->candidate));
    g->groups.clear();
    g->candidates.clear();
	for (int i = 0; i < g->Triangles.size(); i++) {
		if (!g->Triangles[i].is_broken) {
			if (!g->Triangles[i].is_balanced) {
                if (version == 2) {
                    if (!g->candidate[g->Triangles[i].edge1])
                    {
                        g->candidates.push_back(g->Triangles[i].edge1);
                        g->candidate[g->Triangles[i].edge1] = 1;
                    }
                    if (!g->candidate[g->Triangles[i].edge2])
                    {
                        g->candidates.push_back(g->Triangles[i].edge2);
                        g->candidate[g->Triangles[i].edge2] = 1;
                    }
                    if (!g->candidate[g->Triangles[i].edge3])
                    {
                        g->candidates.push_back(g->Triangles[i].edge3);
                        g->candidate[g->Triangles[i].edge3] = 1;
                    }
                } else {
                    for (auto edge: quick_delete(g, g->Triangles[i])) {
                        if (!g->candidate[edge])
                        {
                            g->candidates.push_back(edge);
                            g->candidate[edge] = 1;
                        }
                    }
                }

			}
		}
	}
	int num = 0;
	for (int i = 0; i < g->candidates.size(); i++) {
		Group temp_group;
		if (g->support[g->candidates[i]] == k - 2) {
			if (g->grouped[g->candidates[i]] == -1) {
				queue<int> q;
				q.push(g->candidates[i]);
				g->grouped[g->candidates[i]] = g->groups.size();
				temp_group.edges.push_back(g->candidates[i]);
				while (!q.empty()) {
					int sub = q.front();
					q.pop();
					for (int j = 0; j < g->in_which_triangle[sub].size(); j++) {
						int temp = g->in_which_triangle[sub][j];
						if (!g->Triangles[temp].is_broken) {
							if (g->Triangles[temp].is_balanced) {
								if (g->support[g->Triangles[temp].edge1] == k - 2 && g->grouped[g->Triangles[temp].edge1] == -1) {
									q.push(g->Triangles[temp].edge1);
									temp_group.edges.push_back(g->Triangles[temp].edge1);
									g->grouped[g->Triangles[temp].edge1] = g->groups.size();

								}
								if (g->support[g->Triangles[temp].edge2] == k - 2 && g->grouped[g->Triangles[temp].edge2] == -1) {
									q.push(g->Triangles[temp].edge2);
									temp_group.edges.push_back(g->Triangles[temp].edge2);
									g->grouped[g->Triangles[temp].edge2] = g->groups.size();

								}
								if (g->support[g->Triangles[temp].edge3] == k - 2 && g->grouped[g->Triangles[temp].edge3] == -1) {
									q.push(g->Triangles[temp].edge3);
									temp_group.edges.push_back(g->Triangles[temp].edge3);
									g->grouped[g->Triangles[temp].edge3] = g->groups.size();

								}
							}
						}
					}
				}
			}
		}
		else if (g->support[g->candidates[i]] > k - 2)
		{
			g->grouped[g->candidates[i]] = g->groups.size();
			temp_group.edges.push_back(g->candidates[i]);
		}

		if (temp_group.edges.size() > 0)
		{
			int min = MAX_E + 1;
			for (int j = 0; j < temp_group.edges.size(); j++)
			{
				if (temp_group.edges[j] < min)
					min = temp_group.edges[j];
			}
			temp_group.typ_edge = min;
			g->groups.push_back(temp_group);
		}
	}

}
int update_follower(Graph *g, int i, int min) {

    g->break_unb[i] = 0;
    g->followers[i].clear();
    if (!g->groups[i].is_delete) {
        queue<int> q;
        for (int j = 0; j < g->groups[i].edges.size(); j++) {
            if (!g->is_delete_e[g->groups[i].edges[j]]) {
                q.push(g->groups[i].edges[j]);
                g->temp_delete_e[g->groups[i].edges[j]] = 1;
                g->followers[i].push_back(g->groups[i].edges[j]);
            }	
        }
        while (!q.empty()) {
            //counts++;
            int sub = q.front();
            q.pop();
            if (sub < g->groups[i].typ_edge)
                g->groups[i].typ_edge = sub;
            for (int j = 0; j < g->in_which_triangle[sub].size(); j++) {
                int temp = g->in_which_triangle[sub][j];
                if (!g->Triangles[temp].is_broken) {
                    if (g->Triangles[temp].is_balanced)
                    {
                        g->link[g->Triangles[temp].edge1]++;
                        g->link[g->Triangles[temp].edge2]++;
                        g->link[g->Triangles[temp].edge3]++;
                    }
                    else
                        g->break_unb[i]++;
                    if (g->link[g->Triangles[temp].edge1] == 1)
                        g->is_linked.push_back(g->Triangles[temp].edge1);
                    if (g->link[g->Triangles[temp].edge2] == 1)
                        g->is_linked.push_back(g->Triangles[temp].edge2);
                    if (g->link[g->Triangles[temp].edge3] == 1)
                        g->is_linked.push_back(g->Triangles[temp].edge3);

                    if (g->support[g->Triangles[temp].edge1] - g->link[g->Triangles[temp].edge1] < k - 2 && !g->temp_delete_e[g->Triangles[temp].edge1])
                    {
                        g->temp_delete_e[g->Triangles[temp].edge1] = 1;
                        q.push(g->Triangles[temp].edge1);
                        g->followers[i].push_back(g->Triangles[temp].edge1);
                    }


                    if (g->support[g->Triangles[temp].edge2] - g->link[g->Triangles[temp].edge2] < k - 2 && !g->temp_delete_e[g->Triangles[temp].edge2])
                    {
                        g->temp_delete_e[g->Triangles[temp].edge2] = 1;
                        q.push(g->Triangles[temp].edge2);
                        g->followers[i].push_back(g->Triangles[temp].edge2);
                    }


                    if (g->support[g->Triangles[temp].edge3] - g->link[g->Triangles[temp].edge3] < k - 2 && !g->temp_delete_e[g->Triangles[temp].edge3])
                    {
                        g->temp_delete_e[g->Triangles[temp].edge3] = 1;
                        q.push(g->Triangles[temp].edge3);
                        g->followers[i].push_back(g->Triangles[temp].edge3);
                    }

                    // for (auto A:  findLongestDistanceFromStartVertex(g) ) {
                    //     for (int i_o = 0; i_o < g->adj_edge[A].size(); i_o++) {
                    //         int edge_index_ = g->adj_edge[A][i_o];
                    //         // Check if this edge index appears in the list for v2
                    //         if (!g->temp_delete_e[edge_index_] ) {
                    //             g->temp_delete_e[edge_index_] = 1;
                    //             q.push(edge_index_);
                    //             g->followers[i].push_back(edge_index_);
                    //         }
                    //     }
                    // }

                }
            }
        }
        //if (break_unb[i] == 0)
            //groups[i].is_delete = 1;
        bool query_inside = false;
        for (auto edge_d : g->adj_edge[start_vertex]) {
            if (!g->is_delete_e[edge_d] && !g->temp_delete_e[edge_d]){
                query_inside = true;
                break;
            }
        }
        if (true) {
            if (g->break_unb[i]>0)
            {
                if (g->followers[i].size() < g->followers[min].size())
                    min = i;
                else if (g->followers[i].size() == g->followers[min].size()) {
                    if (g->groups[i].typ_edge < g->groups[min].typ_edge)
                        min = i;
                }
            }
        } else {
            g->groups[i].is_delete = true;
            g->break_unb[i] = 0;
            g->followers[i].clear();


        }
        for (int j = 0; j < g->is_linked.size(); j++)
        {
            g->link[g->is_linked[j]] = 0;
            g->temp_delete_e[g->is_linked[j]] = 0;
        }
        g->is_linked.clear();
    }
    return min;

}


void Get_result(Graph *g) {

    int counts = 0;
    // int min = MAX_E;
    // for (int i = 0; i < g->groups.size(); i++){
    //     min = update_follower(g, i, min);
    // }
	while (g->unbalance_num > 0) {
        counts++;
        int min = MAX_E;
        candidate_lemma(g, 1);
        for (int i = 0; i < g->groups.size(); i++){
            min = update_follower(g, i, min);
        }
        // cout<<g->groups.size()<<"   test\n";
		if (min != MAX_E)
		{
            cout<< min<<endl;
            g->is_delete_e[min] = 1;
            g->temp_delete_e[min] = 1;
			for (int i = 0; i < g->followers[min].size(); i++)
			{
                int del_edge = g->followers[min][i];
				g->is_delete_e[del_edge] = 1;
				g->temp_delete_e[del_edge] = 1;

				//support[followers[min][i]] = 0;
			}
			g->size_of_truss -= g->followers[min].size();
			for (int i = 0; i < g->followers[min].size(); i++) {
				for (int j = 0; j < g->in_which_triangle[g->followers[min][i]].size(); j++) {
					int temp = g->in_which_triangle[g->followers[min][i]][j];
					if (!g->Triangles[temp].is_broken) {
						if (g->Triangles[temp].is_balanced) {
							if (!g->is_delete_e[g->Triangles[temp].edge1])
								--g->support[g->Triangles[temp].edge1];
							if (!g->is_delete_e[g->Triangles[temp].edge2])
								--g->support[g->Triangles[temp].edge2];
							if (!g->is_delete_e[g->Triangles[temp].edge3])
								--g->support[g->Triangles[temp].edge3];
						}
						else
							g->unbalance_num--;
						g->Triangles[temp].is_broken = 1;
					}
				}
			}
			g->groups[min].is_delete = 1;
			if (counts % 100 == 0) {
				cout << "sizeof truss" << g->size_of_truss << endl;
				cout << "unb num:" << g->unbalance_num << endl;
				cout << counts << " " << g->groups[min].typ_edge << endl;
			}
            delete_on_radius(g);


            // int cur_min = min;
            // vector<int> change_set;
            // for (int i = 0; i < g->followers[cur_min].size(); i++)
			// {
            //     int group_id = g->grouped[g->followers[cur_min][i]];
            //     if (group_id > -1){
            //         // if (group_id != cur_min) {
            //         bool t1 = false;
            //         for (auto edge_inside: g->groups[group_id].edges) {
            //             if (!g->is_delete_e[edge_inside]) {
            //                 change_set.push_back(group_id);
            //                 t1 = true;
            //                 break;
            //             }
            //         }
            //         if (!t1) g->groups[group_id].is_delete = 1;                        
            //         // }
            //     }
			// 	//support[followers[min][i]] = 0;
			// }
            // std::sort(change_set.begin(), change_set.end());
            // change_set.erase(std::unique(change_set.begin(), change_set.end()), change_set.end());

            // for (auto id : change_set)//对changeSet中节点的相关信息进行更新
            // {
            //     if (!g->groups[id].is_delete) update_follower(g, id, min);
            // }

            // min = MAX_E;
            // for (int i = 0; i < g->groups.size(); i++){
            //     if (g->break_unb[i]>0 && (!g->groups[i].is_delete))
            //     {
            //         if (g->followers[i].size() < g->followers[min].size())
            //             min = i;
            //         else if (g->followers[i].size() == g->followers[min].size()) {
            //             if (g->groups[i].typ_edge < g->groups[min].typ_edge)
            //                 min = i;
            //         }
            //     }
            // }
			//out << "sizeof truss" << size_of_truss << endl;
			//out << "unb num:" << unbalance_num << endl;
			//out << counts << " " << groups[min].typ_edge << endl;
			
		} 
        // cout<<"test2\n";
	
	}
	cout << "delete " << counts << " times" << endl;
	
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


	GetmaximumKtruss(g);
	if (!if_query_inside(g)) {
		cout << "Cannot obtain any result \n";
		return false;

	}


	// record all nodes with k hops, hop is initially set to 2
	memset(hop_num, -1, sizeof(hop_num));
	queue<pair<int,int>> q;
	int hop = 2;
	int max_hop = -1;
	q.push(make_pair(src,0));
	hop_num[src] = 0;
	
	while (!q.empty()) {
		int curr_v = q.front().first;
		int curr_hop = q.front().second;
		q.pop();
		for (int i = 0; i < g->adj_edge[curr_v].size(); i++) {
			int edgeIndex = g->adj_edge[curr_v][i];
			int v = (g->all_edge_pairs[edgeIndex].v1 == curr_v)
					? g->all_edge_pairs[edgeIndex].v2
					: g->all_edge_pairs[edgeIndex].v1;
			if (g->is_delete_e[edgeIndex]) continue;

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
		Graph *g_hop = new Graph();
		g_hop =  GetKtrusswith_Nhops(curr_hop, k, g);
		
		if (g_hop->size_of_truss > 0 && if_query_inside(g_hop)) {

			global_hop = curr_hop;
			cout << "===================Deleting NegativeTriangle ==================\n";
            if(g_hop->unbalance_num > 0)  {
                candidate_lemma(g_hop, 1);
                Get_result(g_hop);
			    // if (!removeNegativeTriangle(g_hop))	goto hop_add;


            }
            if (!if_query_inside(g_hop)) {
                cout << "Cannot obtain any result \n";
                // return false;

	        }
			
            vector<int> longest_list = findLongestDistanceFromStartVertex(g_hop);

			cout<<"---start to delete node \n"<<endl;
			while (g_hop->diameter > (curr_hop)) {
				for (auto i : longest_list) {
					if ( !delete_on_node(i,g_hop)) goto hop_add;
				}
				longest_list = findLongestDistanceFromStartVertex(g_hop);
			}
            cout<<"check3\n";

			*g = *g_hop;
			return true;
		} 
		hop_add:
		curr_hop += 1;
		delete(g_hop);
	}
	cout<<"---calculation fail! \n"<<endl;
	return false;
}






int main(int argc, char** argv) {
	if (argc < 4) {
		printf("Not enough args: ./ad start_vertex k filename\n");
		return 0;
	}

    start_vertex = atoi(argv[1]);
	k = atoi(argv[2]);
    filename = argv[3];

	outname = filename + "_ad_solution.txt";
	filename += ".txt";
	
	Graph* g = build_graph();
	cout << "graph build" << endl;
    cout << "filename: " << filename <<endl;
    cout << "start vertex: " << start_vertex<<endl;
    cout << "k: " << k <<endl;
    cout << "|V|: " << v_num <<endl;
    cout << "|E|: " << e_num <<endl;

	struct timeval start, end;
    double timeuse;
	bool result;
    gettimeofday(&start, NULL);
	result = GetKtruss(start_vertex, k, g);
    gettimeofday(&end, NULL);
    timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
   	cout << "advanced-exact time: " << timeuse << '\n' << endl;
	if (result) print_result(g);

	return 0;
}