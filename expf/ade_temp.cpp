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
	// vector<int> is_linked;
	int book[MAX_V];
	vector<int> is_booked;
	int is_delete_e[MAX_E];
	int is_delete_vec[MAX_V];
	// int link[MAX_E];
	int temp_delete_e[MAX_E];
	int break_unb[MAX_E];
	int support[MAX_E];
	int edge_num;


	int unbalance_num = 0;
	vector<Triangle> Triangles;
	int size_of_truss;
	// vector<int> followers[MAX_E + 1];
	vector<int> two_dimension;
	// int grouped[MAX_E];
	// int candidate[MAX_E];
	// vector<int> candidates;
	// vector<Group> groups;
	// vector<int> id;

	int diameter;
	vector<int> path;
} Graph;




clock_t start, finish;
int start_vertex, k, global_hop;


// modify from here
int hop_num[MAX_V];
int node_deleted[MAX_V];
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
	// memset(g->link, 0, sizeof(g->link));
	// memset(g->grouped, -1, sizeof(g->grouped));
	// memset(g->candidate, 0, sizeof(g->candidate));
	g->two_dimension.resize(MAX_E);
	ifstream inputFile{ filename };
    
    //ifstream inputFile("data/temp1.txt");
    if (!inputFile) {
        cerr << "Input file could not be opened!" << endl;
        exit(EXIT_FAILURE);
    }

    //inputFile >> vertexnum >> edgenum;

	inputFile >> v_num >> e_num;
	int max_v_number = 0;

	for (int i = 0; i < e_num; i++)
	{
		int v1, v2, sign;
		inputFile >> v1 >> v2 >> sign;
		g->all_edge_pairs[i].v1 = v1;
		g->all_edge_pairs[i].v2 = v2;
		max_v_number = max(max(max_v_number,v1),v2);

		g->all_edge_pairs[i].sign = sign;
		g->vec[v1].push_back(v2);
		g->vec[v2].push_back(v1);
		g->adj_edge[v1].push_back(i);
		g->adj_edge[v2].push_back(i);
	}
	// for (int i = 0; i < MAX_E; i++) {
	// 	g->followers[MAX_E].push_back(0);
	// }
	g->diameter = e_num;
	g->edge_num = e_num;
	v_num = max_v_number+1;


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
	memset(g->is_delete_vec, 1, sizeof(g->is_delete_vec));
	memset(g->support, 0, sizeof(g->support));
	memset(g->book, 0, sizeof(g->book));
	memset(g->break_unb, 0, sizeof(g->break_unb));
	memset(g->temp_delete_e, 0, sizeof(g->temp_delete_e));
	// memset(g->link, 0, sizeof(g->link));
	// memset(g->grouped, -1, sizeof(g->grouped));
	// memset(g->candidate, 0, sizeof(g->candidate));
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
			g->is_delete_vec[v1] = 0;
			g->is_delete_vec[v2] = 0;
			// g->is_delete_e[nhop_e_num] = 0;
			
			nhop_e_num ++;
		} 
	}
	// for (int i = 0; i < MAX_E; i++) {
	// 	g->followers[i].push_back(0);
	// }
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
	// cout << "counts:" << counts << endl;
	//unbalance_num = Triangles.size() - num_of_balance;
	// cout << "orangial:" << g->unbalance_num << endl;
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
	// cout << "size of KTruss: " << g->size_of_truss << endl;
	// cout << "truss unb num: " << g->unbalance_num << endl;

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


// void findLongestPath(Graph *g) {
// 	int pathLength = -1;

// 	for (int i = 1; i < v_num - 1; i++) {
// 		if (g->is_delete_vec[i]) continue; // skip the vertices marked as deleted

// 		memset(visited, -1, sizeof(visited));
// 		queue<pair<int,pair<int, int>>> q;
// 		vector<vector<int>> paths(v_num);
// 		paths[i].push_back(i);
// 		q.emplace(i, make_pair(0, i));

// 		while (!q.empty()) {
// 			auto curr_pair = q.front();
// 			q.pop();

// 			if (visited[curr_pair.first] < curr_pair.second.first && visited[curr_pair.first] != -1) continue;

// 			visited[curr_pair.first] = curr_pair.second.first;

// 			if (pathLength < curr_pair.second.first) {
// 				g->paths.clear();
// 				g->paths.push_back(paths[curr_pair.first]);
// 				pathLength = curr_pair.second.first;
// 			} else if (pathLength == curr_pair.second.first) {
// 				g->paths.push_back(paths[curr_pair.first]);
// 			}

// 			for (int idx = 0; idx < g->adj_edge[curr_pair.first].size(); idx++) {
// 				if (g->is_delete_e[g->adj_edge[curr_pair.first][idx]]) {
// 					continue;
// 				}
// 				int v;
// 				if (g->all_edge_pairs[g->adj_edge[curr_pair.first][idx]].v1 == curr_pair.first) {
// 					v = g->all_edge_pairs[g->adj_edge[curr_pair.first][idx]].v2;
// 				} else {
// 					v = g->all_edge_pairs[g->adj_edge[curr_pair.first][idx]].v1;
// 				}

// 				// Skip the vertices marked as deleted
// 				if (g->is_delete_vec[v]) continue;

// 				paths[v] = paths[curr_pair.first];
// 				paths[v].push_back(v);
				
// 				q.emplace(v, make_pair(curr_pair.second.first + 1, v));
// 			}
// 		}
// 	}
// 	g->diameter = pathLength;
// }




vector<int> findLongestPath(Graph *g) {
	// little optimization by replace v_num by v_num-1
	//int dist[MAX_V][MAX_V] = {MAX_E + 1};
	// cout << "in\n";
	int pathLength = -1;
	vector<int> path;
	
	for (int i = 1; i < v_num-1; i++) {
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
vector<int> findLongestDistanceFromStartVertex(Graph *g) {
    vector<int> distances(MAX_V, -1);
    vector<vector<int>> paths(MAX_V);
    vector<bool> visited(MAX_V, false);
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
	findLongestPath(g);
	set<int> vertices_set;
	cout <<  " ====result graph: \n";
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
	cout << "final truss unbalance num: " << g->unbalance_num << endl;
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
	// cout << "delete on edge fail\n";
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
        if (!unb.is_balanced && !unb.is_broken) {
            if (hop_num[g->all_edge_pairs[unb.edge1].v1] = hop_num[g->all_edge_pairs[unb.edge1].v2]) {
                if (!delete_on_edge(unb.edge1, g,false)) return false;
            } else if (hop_num[g->all_edge_pairs[unb.edge2].v1] = hop_num[g->all_edge_pairs[unb.edge2].v2]) {
                if (!delete_on_edge(unb.edge1, g,false)) return false;
            } else {
                if (!delete_on_edge(unb.edge1, g,false)) return false;
            }
        }
    }
	return true;
}


bool removeNegativeTriangle(Graph* g) {
	int curr_diameter = g->diameter;
	if (curr_diameter < 1) return false;
	if (g->unbalance_num <= 0) return true;

	queue<Graph*> Queue;
	Graph *best_g = new Graph();
	*best_g = *g;
	// g->diameter ++;
	Queue.push(best_g);
	int que_in,que_out = 0;
	Graph *temp_g = new Graph();
	*temp_g = *g;
    if (!quickremoveNegativeTriangle(temp_g)) {
		cout << "quicl delete fail\n";
		*temp_g = *g;
		bool check = false;
		for (auto unb : temp_g->Triangles) {
			if (!unb.is_balanced && !unb.is_broken) {
				if (!check && delete_on_edge(unb.edge1, temp_g,false)) check = true;
				if (!check && delete_on_edge(unb.edge2, temp_g,false)) check = true;
				if (!check && delete_on_edge(unb.edge2, temp_g,false)) check = true;
				break;
			}
		}
		if (!check) {
			cout<< "cannot delete \n";
			delete(temp_g);
			return false;
		}

	}
	*g = *temp_g;
	delete(temp_g);
	findLongestDistanceFromStartVertex(g);

	while (!Queue.empty()) {
		
		int size = Queue.size();
		Graph *graph_first = Queue.front();
		Queue.pop();
		queue<Graph*> filteredQueue;

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
						filteredQueue.push(left);
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
						filteredQueue.push(middle);
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
						filteredQueue.push(right);
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

		queue<Graph*> TempQueue;
		while (!filteredQueue.empty()) {
			Graph* currentGraph = filteredQueue.front();
			filteredQueue.pop();

			// Check if currentGraph is a subgraph of any other subgraph in filteredQueue.
			bool isSubgraphOfOther = false;
			queue<Graph*> tempQueue; // Copy filteredQueue to a temporary queue.

			while (!Queue.empty()) {
				Graph* graphInFiltered = Queue.front();
				Queue.pop();

				if (!isSubgraphOfOther && isSubgraph(currentGraph, graphInFiltered)) {
					isSubgraphOfOther = true;
					// cout << "prunning queue\n";
					delete (currentGraph);
				}
				tempQueue.push(graphInFiltered);
			}

			// If currentGraph is not a subgraph of any other subgraph in filteredQueue, keep it.
			if (!isSubgraphOfOther) {
				Queue.push(currentGraph);
			} 
		}

		// Now, filteredQueue contains only the subgraphs that are not subgraphs of any other subgraph.
		// You can continue using filteredQueue or copy its contents back to the original Queue.

		// Clean up the original Queue if needed.
		while (!filteredQueue.empty()) {
			delete filteredQueue.front();
			filteredQueue.pop();
		}

	}
    if (g->unbalance_num >0) return false;

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
	// cout << "counts:" << counts << endl;
	// //unbalance_num = Triangles.size() - num_of_balance;
	// cout << "orangial:" << g->unbalance_num << endl;
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
	// cout << "size of KTruss" << g->size_of_truss << endl;
	// cout << "truss unb num:" << g->unbalance_num << endl;






}




bool GetKtruss(int src, int k, Graph* g) {

	// check if the original graph has k-truss including src


	GetmaximumKtruss(g);
	if (!if_query_inside(g)) {
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
	// cout << "max hop is " << max_hop << "\n";

	// add one hop at a time
	
	bool isSuccess = false;

	int curr_hop = 1;

	while (curr_hop >= 1 && curr_hop <= max_hop) {
		// cout << "--------when hop is " << curr_hop << "\n";
		Graph *g_hop = new Graph();
		g_hop =  GetKtrusswith_Nhops(curr_hop, k, g);
		
		if (g_hop->size_of_truss > 0 && if_query_inside(g_hop)) {

			global_hop = curr_hop;
			cout << "===================Deleting NegativeTriangle ==================\n";
            
			if (!removeNegativeTriangle(g_hop))	goto hop_add;
            vector<int> longest_list = findLongestDistanceFromStartVertex(g_hop);

			while (g_hop->diameter > (curr_hop)) {
				cout<<"---start to delete diameter \n"<<endl;

				for (auto i : longest_list) {
					if ( !delete_on_node(i,g_hop)) goto hop_add;
				}
				longest_list = findLongestDistanceFromStartVertex(g_hop);
			}
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





double cc_vertex(Graph* g, int v1)  {


    vector<bool> visited(MAX_V, false);
    int count = 0;
    queue<int> q;

    visited[v1] = true; // 将 v1 标记为已访问
    q.push(v1); // 将 v1 放入队列

    while (!q.empty()) { // 当队列不为空时
        int v = q.front(); // 取出队列的第一个元素
        q.pop(); // 从队列中删除第一个元素
        count++; // 访问的 vertex 数量增加
        // 遍历当前 vertex 所有的相邻 vertex
        for (int i = 0; i < g->adj_edge[v].size(); i++) {
            int edge_id = g->adj_edge[v][i];  // 相邻的边的 ID
            int v2 = g->all_edge_pairs[edge_id].v1;
            if (v2 == v) v2 = g->all_edge_pairs[edge_id].v2;


            if (!visited[v2] && !g->is_delete_e[edge_id]) {  // 如果相邻的 vertex 没有被访问，并且连接的边没有被删除
                visited[v2] = true; // 标记相邻 vertex 为已访问
                q.push(v2); // 将相邻 vertex 放入队列
            }
        }
    }
    return count;


}

double calcualte_vertex (Graph* g) {
    int total = 0;
	set<int> vertices_set;


	// Collect all unique vertices from all_edge_pairs into the set
	for (int i = 0; i < g->edge_num; i++) {
		if (!g->is_delete_e[i]) {
			int v1 = g->all_edge_pairs[i].v1;
			int v2 = g->all_edge_pairs[i].v2;
			// int sig = g->all_edge_pairs[i].sign;
			vertices_set.insert(v1);
			vertices_set.insert(v2);
        	// outfile << v1 << " " << v2 << " " << sig << endl;
		}
	}
    return vertices_set.size();
}

double calcualte_edge (Graph* g_ori) {
    int total = 0;

    for (int i = 0; i < g_ori->edge_num; i++) {
        if (!g_ori->is_delete_e[i]){
            total++;
        }
    }
    
    return total;
}




int main(int argc, char** argv)  {
	if (argc < 6) {
		printf("Not enough args: ./base  k datafilename queryfilename i_th answer_floder\n");
		return 0;
	}

    // start_vertex = atoi(argv[1]);
	k = atoi(argv[1]);
    filename = argv[2];
    string query_file = argv[3];
	int query_number = atoi(argv[4]) - 1;

	// outname = strcat(strcat(argv[3], "_f/"), strcat(argv[4], ".txt"));
	// outname = sprintf(outname, "%s_f/%s.txt", argv[3], argv[4]);
	filename += ".txt";
	query_file += ".txt";
	



	Graph* original = build_graph();
	// cout << "graph build" << endl;
    // cout << "filename: " << filename <<endl;
    // cout << "start vertex: " << start_vertex<<endl;
    // cout << "k: " << k <<endl;
    // cout << "|V|: " << v_num <<endl;
    // cout << "|E|: " << e_num <<endl;

	struct timeval start, end;
    double timeuse;
    double timetotal = 0.00;
    double total_diameter = 0.00;
    double total_size_of_truss = 0.00;
    double total_unbalance_num=0.00;
    double total_percentage = 0.00;
    double total_density = 0.00;
    double attempt_time = 0.00;

    vector<int> query_data;
    // 打开文件
    ifstream file(query_file); // 这里替换为你的文件名
    if (!file) {
        cerr << "Unable to open file\n";
        // exit(1); // 退出程序，因为文件打开失败
    }

    // 从文件中读取数据
    string line;
    while (getline(file, line)) {
        query_data.push_back(stoi(line));
    }
    file.close();


    // for (auto query: query_data) {

	cout << "query_data: NO"<<query_number<<"; vertex_id: "<<query_data[query_number] << endl;
	start_vertex = query_data[query_number];
	Graph *g = (Graph*) malloc(sizeof(Graph));
	*g = *original;
	gettimeofday(&start, NULL);
	bool result = GetKtruss(start_vertex, k, g);
	gettimeofday(&end, NULL);
	timetotal = timetotal + ((end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0);

	if (result) {
		cout << "calculation success!"<<endl;
		attempt_time ++;
		total_diameter += g->diameter;
		total_size_of_truss += g->size_of_truss;
		total_unbalance_num += g->unbalance_num;
		// int calcualte_vertex, calcualte_edge;
		double vertexCount  = calcualte_vertex(g);
		double edgeCount  = calcualte_edge(g);
		total_density = total_density + ((2*edgeCount) / ((vertexCount) * (vertexCount-1)));

		Graph *temp = (Graph*) malloc(sizeof(Graph));
		*temp = *original;
		GetmaximumKtruss(temp);
		double large_graph = cc_vertex(temp, start_vertex);

		total_percentage += (vertexCount/large_graph);
		free(temp);
		
	}
	free(g);
    

    // string directory = string(argv[5]);
    // std::filesystem::create_directories(directory);
    // outname = directory + "/"+string(argv[4])+".txt";
    outname = string(argv[5]);
	std::ofstream outfile(outname, std::ios::app);

	// ofstream outfile(outname);
	if (!outfile) {
        cerr << "Error: Unable to open the file " << outname << endl;
		// ofstream newFile(outname);
    }

	outfile <<  "====advanced_exact result: \n";
	// Collect all unique vertices from all_edge_pairs into the set
    outfile << filename <<" " << query_file  << endl;
	outfile << "query_data_No: "<<query_number<<"; vertex_id: "<<query_data[query_number] << endl;
    outfile <<"Attempt: " << attempt_time  << endl;
    outfile <<"time: " << timetotal << endl;
    outfile <<"diameter: " << total_diameter << endl;
    outfile <<"total_size_of_truss: " << total_size_of_truss << endl;
    outfile <<"total_unbalance_num: " << total_unbalance_num << endl;
    outfile <<"total_percentage: " << total_percentage << endl;
    outfile <<"total_density: " << total_density << endl;
    outfile.close();
	return 0;
}
