#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<string>
#include<vector>
#include<queue>
#include<fstream>
//#include<Windows.h>
#include<cstring>
#include<set>
#include<limits>
#include<sys/time.h>
#include<algorithm>


using namespace std;
const int MAX_V = 516000;
const int MAX_E = 1000000;
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
	int id;

	bool is_balanced;
	bool is_broken;
	Triangle() :edge1(-1), edge2(-1), edge3(-1),  id(-1), is_balanced(false), is_broken(false) {}
};
vector<Triangle> Triangles;

Edge all_edge_pairs[MAX_E];
// vector<int> vec[MAX_V];
vector<int> adj_edge[MAX_V];
//vector<int> edge_to_edge[MAX_E];
vector<int> in_which_triangle[MAX_E];
typedef struct {

	vector<bool> is_delete_e{vector<bool>(MAX_E,false)};
	vector<bool> is_delete_vec{vector<bool>(MAX_V,false)};
	vector<bool> Triangle_balance{vector<bool>(MAX_E,false)};
	vector<bool> Triangle_broken{vector<bool>(MAX_E,false)};
	int support[MAX_E];
	int edge_num;
	int unbalance_num = 0;
	// vector<Triangle> Triangles;
	int size_of_truss;

	
	int diameter;
	int point1;
	int point2;
	// vector<int> path;
	// vector<int> path;
} Graph;




clock_t start, finish;
int start_vertex, k, global_hop;


// modify from here
int hop_num[MAX_V];
int node_deleted[MAX_V];
int visited[MAX_V];

string filename, outname;

int first_round = 0;
vector<Graph*> graphPtrs;

bool delete_on_radius(Graph *g_hop);

Graph* build_graph() {
	Graph *g = new Graph();
	// memset(g->is_delete_e, false, sizeof(g->is_delete_e));
	// memset(g->is_delete_vec, false, sizeof(g->is_delete_vec));
	memset(g->support, 0, sizeof(g->support));
	// memset(book, 0, sizeof(book));
	// memset(g->break_unb, 0, sizeof(g->break_unb));
	// memset(g->temp_delete_e, 0, sizeof(g->temp_delete_e));
	// memset(g->link, 0, sizeof(g->link));
	// memset(g->grouped, -1, sizeof(g->grouped));
	// memset(g->candidate, 0, sizeof(g->candidate));
	// g->two_dimension.resize(MAX_E);
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
		all_edge_pairs[i].v1 = v1;
		all_edge_pairs[i].v2 = v2;
		g->is_delete_e[i] = false;
		g->is_delete_vec[v1] = false;
		g->is_delete_vec[v2] = false;
		max_v_number = max(max(max_v_number,v1),v2);

		all_edge_pairs[i].sign = sign;
		// g->vec[v1].push_back(v2);
		// g->vec[v2].push_back(v1);
		adj_edge[v1].push_back(i);
		adj_edge[v2].push_back(i);
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
	for (auto i : adj_edge[start_vertex]) {
		if (!g_ori->is_delete_e[i]){
			return true;
		}
	}
	return false;
}


Graph* GetKtrusswith_Nhops(int n, int k, Graph* g_ori) {
	

	Graph *g = new Graph();
	*g = *g_ori;
	// int counts = 0;
	// int two_dimension[MAX_E];
	// int book[MAX_V];
	// vector<int> is_booked;

	// //////////   No subgraph build Method /////////////////////////
	// // *g = *g_ori;
	// // g->size_of_truss = e_num;
	// //////////////////////////////////////////////////////////////


	// ///////////// build subgraph ///////////////////////////
	// // memset(g->is_delete_e, 0, sizeof(g->is_delete_e));
	// // memset(g->is_delete_vec, 1, sizeof(g->is_delete_vec));
	// memset(g->support, 0, sizeof(g->support));

	// Triangles.clear();
	// for (int i = 0; i < g_ori->edge_num; i++)
	// {
	// 	int v1 = all_edge_pairs[i].v1;
	// 	int v2 = all_edge_pairs[i].v2;
	// 	if (hop_num[v1] <= n && hop_num[v2] <= n
	// 	&& hop_num[v1] != -1 && hop_num[v2] != -1) {
			
	// 		// all_edge_pairs[nhop_e_num].v1 = v1;
	// 		// all_edge_pairs[nhop_e_num].v2 = v2;
	// 		// all_edge_pairs[nhop_e_num].sign = g_ori->all_edge_pairs[i].sign;
	// 		// g->vec[v1].push_back(v2);
	// 		// g->vec[v2].push_back(v1);
	// 		// adj_edge[v1].push_back(nhop_e_num);
	// 		// adj_edge[v2].push_back(nhop_e_num);
	// 		g->is_delete_vec[v1] = 0;
	// 		g->is_delete_vec[v2] = 0;
	// 		// g->is_delete_e[nhop_e_num] = 0;
			
	// 		// nhop_e_num ++;
	// 	} 
	// }
	// // for (int i = 0; i < MAX_E; i++) {
	// // 	g->followers[i].push_back(0);
	// // }
	// g->diameter = g_ori->edge_num;
	// g->edge_num = g_ori->edge_num;
	// g->size_of_truss = g_ori->edge_num;

	// for (int i = 0; i < g->edge_num; i++) {

	// 	Triangle temp_triangle;
	// 	// Assigning values to each edge
	// 	// temp_triangle.v1 = all_edge_pairs[i].v1;
	// 	// temp_triangle.v2 = all_edge_pairs[i].v2;
		

	// 	temp_triangle.edge1 = i;

	// 	for (int j = 0; j < adj_edge[all_edge_pairs[i].v1].size(); j++) {
	// 		int edgeIndex = adj_edge[all_edge_pairs[i].v1][j];
	// 		int v = (all_edge_pairs[edgeIndex].v1 == all_edge_pairs[i].v1)
	// 				? all_edge_pairs[edgeIndex].v2
	// 				: all_edge_pairs[edgeIndex].v1;

	// 		// If the vertex v is not within n hops of the source, skip this iteration
	// 		// if (hop_num[v] > n)
	// 		// 	continue;

	// 		counts++;
	// 		book[v] = 1;
	// 		two_dimension[v] = j;
	// 		is_booked.push_back(v);
	// 	}
	// 	for (int j = 0; j < adj_edge[all_edge_pairs[i].v2].size(); j++) {
	// 		int edg2 = adj_edge[all_edge_pairs[i].v2][j];
	// 		int v = (all_edge_pairs[edg2].v1 == all_edge_pairs[i].v2)
	// 				? all_edge_pairs[edg2].v2
	// 				: all_edge_pairs[edg2].v1;

	// 		// If the vertex v is not within n hops of the source, skip this iteration
	// 		// if (hop_num[v] > n)
	// 		// 	continue;

	// 		counts++;
	// 		if (book[v]) {
	// 			int edg1 = adj_edge[all_edge_pairs[i].v1][two_dimension[v]]; //?edg1
	// 			if (edg1 > i && edg2 > i) {
	// 				temp_triangle.edge2 = edg1;
	// 				temp_triangle.edge3 = edg2;
	// 				// temp_triangle.v3 = v;
	// 				temp_triangle.id = Triangles.size();
					
	// 				if ((all_edge_pairs[edg1].sign + all_edge_pairs[edg2].sign + all_edge_pairs[i].sign) == -1) {
	// 					temp_triangle.is_balanced = 1;
	// 					g->Triangle_balance[temp_triangle.id] = true;
	// 				}
	// 				else if ((all_edge_pairs[edg1].sign + all_edge_pairs[edg2].sign + all_edge_pairs[i].sign) == 3) {
	// 					temp_triangle.is_balanced = 1;
	// 					g->Triangle_balance[temp_triangle.id] = true;

	// 				}
	// 				Triangles.push_back(temp_triangle);
	// 				in_which_triangle[i].push_back(Triangles.size() - 1);
	// 				in_which_triangle[edg1].push_back(Triangles.size() - 1);
	// 				in_which_triangle[edg2].push_back(Triangles.size() - 1);
	// 			}
	// 		}
	// 	}

	// 	for (int j = 0; j < is_booked.size(); j++) {
	// 		book[is_booked[j]] = 0;
	// 	}
	// 	is_booked.clear();
	// }
	// for (int i = 0; i < Triangles.size(); i++) {
	// 	Triangles[i].id = i;
	// 	g->Triangle_balance[i] = Triangles[i].is_balanced;
	// 	g->Triangle_broken[i] = Triangles[i].is_broken;

	// }
	// for (int i = 0; i < Triangles.size(); i++) {
	// 	// cout<<"size 1: "<< i << "size2:  "<< Triangles[i].id<<endl;
	// 	if (!Triangles[i].is_balanced)
	// 		g->unbalance_num++;
	// 	else {
	// 		g->support[Triangles[i].edge1]++;
	// 		g->support[Triangles[i].edge2]++;
	// 		g->support[Triangles[i].edge3]++;
	// 	}
	// 	if (Triangles[i].is_balanced != g->Triangle_balance[Triangles[i].id]) {
	// 		cout << "error\n";
	// 	}

	// }
	// cout << "counts:" << counts << endl;
	// unbalance_num = Triangles.size() - num_of_balance;
	queue<int> q;
	// 找出所有不满足 support的边
	// 找出所有不满足 support的边
	for (int i = 0; i < g->edge_num; i++) {
		if (!g->is_delete_e[i] ) {
			if (hop_num[all_edge_pairs[i].v1] > n || hop_num[all_edge_pairs[i].v2] > n
		|| hop_num[all_edge_pairs[i].v1] == -1 || hop_num[all_edge_pairs[i].v2] == -1) {			
			// if (g->support[i] < k - 2) {
				g->is_delete_e[i] = 1;
				g->size_of_truss--;
				// g->temp_delete_e[i] = 1;
				q.push(i);
			}
		}
	}
	while (!q.empty()) {
		int sub = q.front();
		q.pop();
		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
		for (int i = 0; i < in_which_triangle[sub].size(); i++) {
			if (!g->Triangle_broken[Triangles[in_which_triangle[sub][i]].id]){
				if (g->Triangle_balance[Triangles[in_which_triangle[sub][i]].id]) {
					// 删除临边
					g->support[Triangles[in_which_triangle[sub][i]].edge1]--;
					g->support[Triangles[in_which_triangle[sub][i]].edge2]--;
					g->support[Triangles[in_which_triangle[sub][i]].edge3]--;
				}
				else 
					g->unbalance_num--;
				// 删除一条边后check 他的临边
				if (!g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1]) {
					if (g->support[Triangles[in_which_triangle[sub][i]].edge1] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge1);
						g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						// g->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						g->size_of_truss--;
					}
				}
				if (!g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2]) {
					if (g->support[Triangles[in_which_triangle[sub][i]].edge2] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge2);
						g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						// g->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						g->size_of_truss--;
					}
				}
				if (!g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3]) {
					if (g->support[Triangles[in_which_triangle[sub][i]].edge3] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge3);
						g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						// g->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						g->size_of_truss--;
					}
				}

				g->Triangle_broken[Triangles[in_which_triangle[sub][i]].id] = true;
			}
		}
	}
	cout << "size of KTruss: " << g->size_of_truss << endl;
	cout << "truss unb num: " << g->unbalance_num << endl;

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

// 			for (int idx = 0; idx < adj_edge[curr_pair.first].size(); idx++) {
// 				if (g->is_delete_e[adj_edge[curr_pair.first][idx]]) {
// 					continue;
// 				}
// 				int v;
// 				if (all_edge_pairs[adj_edge[curr_pair.first][idx]].v1 == curr_pair.first) {
// 					v = all_edge_pairs[adj_edge[curr_pair.first][idx]].v2;
// 				} else {
// 					v = all_edge_pairs[adj_edge[curr_pair.first][idx]].v1;
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







void findLongestPath(Graph *g) {
	// little optimization by replace v_num by v_num-1
	//int dist[MAX_V][MAX_V] = {MAX_E + 1};
	// cout << "in\n";
	int pathLength = -1;
	vector<int> path;
	
	for (int i = 1; i < v_num; i++) {
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
				if (visited[curr_pair.first] <= curr_pair.second.first && visited[curr_pair.first] != -1) continue;

				visited[curr_pair.first] = curr_pair.second.first;
				if (pathLength < curr_pair.second.first) {
					path = curr_pair.second.second;
					pathLength = curr_pair.second.first;

				}
				//cout << "neighour number is " << g->vec[4].size() << "\n";
				for (int idx = 0; idx < adj_edge[curr_pair.first].size(); idx++) {

					if (g->is_delete_e[adj_edge[curr_pair.first][idx]]) {
						continue;
					}
					int v;
					if (all_edge_pairs[adj_edge[curr_pair.first][idx]].v1 == curr_pair.first) {
						v = all_edge_pairs[adj_edge[curr_pair.first][idx]].v2;
					} else {
						v = all_edge_pairs[adj_edge[curr_pair.first][idx]].v1;
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
		if (pathLength >= 2*global_hop) {
			g->diameter = pathLength;
			// g->path = path;
			g->point1 = path[0];
			g->point2 = path[path.size()-1];
			return;
		}
		
	}

	result:
	g->diameter = pathLength;
	// g->path = path;
	g->point1 = path[0];
	g->point2 = path[path.size()-1];
	return;

	// return path;
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

        for (int idx = 0; idx < adj_edge[curr_vertex].size(); idx++) {
            if (g->is_delete_e[adj_edge[curr_vertex][idx]]) {
                continue;
            }

            int v;
            if (all_edge_pairs[adj_edge[curr_vertex][idx]].v1 == curr_vertex) {
                v = all_edge_pairs[adj_edge[curr_vertex][idx]].v2;
            } else {
                v = all_edge_pairs[adj_edge[curr_vertex][idx]].v1;
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
	//g->point1 = longest_path[0];
	//g->point2 = longest_path[longest_path.size()-1];	

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
			int v1 = all_edge_pairs[i].v1;
			int v2 = all_edge_pairs[i].v2;
			int sig = all_edge_pairs[i].sign;
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
		// cout<<"11\n";
		// newG->temp_delete_e[edge_num] = 1;
		q.push(edge_num);
	}

	// cout<<"11\n";

	while (!q.empty()) {
		int sub = q.front();
		q.pop();
		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
		for (int i = 0; i < in_which_triangle[sub].size(); i++) {
			if (!newG->Triangle_broken[Triangles[in_which_triangle[sub][i]].id]){
				if (newG->Triangle_balance[Triangles[in_which_triangle[sub][i]].id]) {
					// 删除临边
					newG->support[Triangles[in_which_triangle[sub][i]].edge1]--;
					newG->support[Triangles[in_which_triangle[sub][i]].edge2]--;
					newG->support[Triangles[in_which_triangle[sub][i]].edge3]--;
				}
				else 
					newG->unbalance_num--;
				// 删除一条边后check 他的临边
				if (!newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1]) {
					if (newG->support[Triangles[in_which_triangle[sub][i]].edge1] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge1);
						newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						// newG->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						newG->size_of_truss--;
					}
				}
				if (!newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2]) {
					if (newG->support[Triangles[in_which_triangle[sub][i]].edge2] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge2);
						newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						// newG->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						newG->size_of_truss--;
					}
				}
				if (!newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3]) {
					if (newG->support[Triangles[in_which_triangle[sub][i]].edge3] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge3);
						newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						// newG->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						newG->size_of_truss--;
					}
				}
				newG->Triangle_broken[Triangles[in_which_triangle[sub][i]].id] = true;
			}
		}
	}
	// cout<<"12\n";


	// if (newG->unbalance_num == 0 && newG->size_of_truss > 0){
	// 	return delete_on_radius(newG);

	// }
	// cout<<"13\n";

	if (newG->size_of_truss > 0 && if_query_inside(newG)) {
		//if (update_dia) findLongestDistanceFromStartVertex(newG);
		if (update_dia) findLongestPath(newG);
		return true;
	} 

	return false;
}

bool delete_on_node(int A, Graph*newG) {

	int edge_num = -1;
	queue<int> q;
	if(newG->is_delete_vec[A]) return true;
	for (int i = 0; i < adj_edge[A].size(); ++i) {
        int edge_index = adj_edge[A][i];
        // Check if this edge index appears in the list for v2
		if (!newG->is_delete_e[edge_index] ) {
			if (!delete_on_edge(edge_index,newG, false)) return false;
		}
    }
	newG->is_delete_vec[A] = 1;
	if (newG->size_of_truss > 0) {
		findLongestPath(newG);
		//cout << "delete on edge ended\n";
		return true;
	} 
	return false;

}
bool delete_on_radius(Graph *g_hop) {
	vector<int> longest_list = findLongestDistanceFromStartVertex(g_hop);
	
	while (g_hop->diameter > (global_hop)) {
		for (auto i : longest_list) {
			if ( !delete_on_node(i,g_hop)) return false;
		}
		findLongestDistanceFromStartVertex(g_hop);
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
	int que_in = 0;
	int que_out = 0;
	if (curr_diameter < 1) return false;
	if (curr_diameter == 1) {
		if(first_round != global_hop) {
			graphPtrs.push_back(g);
		} 
		first_round = global_hop;
		return true;
	}
	queue<Graph*> Queue;
	vector<Graph*> graph_list;
	Graph *best = new Graph();
	*best = *g;
	Queue.push(best);

	
	while (!Queue.empty()) {
		
		int size = Queue.size();
		Graph *left = Queue.front();
		queue<Graph*> filteredQueue;
		Queue.pop();
		que_out++;
		
		// cout << "mid1 queue \n";
		Graph *right = new Graph();
		*right = *left;
		// cout<<"1\n";
		if (delete_on_node(left->point1, left)) {
			bool best = false;
			if ((left->diameter == g->diameter && g->size_of_truss < left->size_of_truss ) ||
			left->diameter < g->diameter) {
				// Graph *swap_temp = g;
				// delete(swap_temp);
				*g = *left;
				best = true;
			}

			if (left->diameter > global_hop) {
			// if (if_query_inside(left) && left->size_of_truss >= g->size_of_truss) {
				que_in++;
				filteredQueue.push(left);
			} else {
				graph_list.push_back(left);
			}
		} else {
			delete(left);
		}
		// cout<<"2\n";

		if (delete_on_node(right->point2, right)) {
			bool best = false;
			if ((right->diameter == g->diameter && g->size_of_truss < right->size_of_truss ) ||
			right->diameter < g->diameter) {
				// Graph *swap_temp = g;
				// delete(swap_temp);

				*g = *right;
				best = true;

			}

			if (right->diameter > global_hop) {
				// if (right->diameter >= g->diameter && )
				que_in++;
				
				filteredQueue.push(right);

			} else {
				graph_list.push_back(right);

			}
		} else {
			delete(right);
		}
		// cout<<"3\n";

	// Check each subgraph in Queue.
		// queue<Graph*> TempQueue;
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

				tempQueue.push(currentGraph);
			} 
			Queue = move(tempQueue);

		}

		// Now, filteredQueue contains only the subgraphs that are not subgraphs of any other subgraph.
		// You can continue using filteredQueue or copy its contents back to the original Queue.

		// Clean up the original Queue if needed.
		while (!filteredQueue.empty()) {
			Graph *temp = filteredQueue.front();
			filteredQueue.pop();
			delete(temp);
		}
		// cout <<"que_in"<<que_in<<endl;
		// cout <<"que_out"<<que_out+Queue.size()<<endl;

		// Copy the filtered subgraphs back to the original Queue.

	}

	if(first_round != global_hop) {
		graphPtrs = graph_list;


	} else {
		for (Graph* graph : graph_list) {
			delete graph;  // 释放动态分配的内存
		}
		graph_list.clear(); 
	}

	first_round = global_hop;

	return true;
}



bool quickremoveNegativeTriangle(Graph* g) {
	int curr_diameter = g->diameter;
	if (curr_diameter < 1) return false;
	if (g->unbalance_num <= 0) return true;

	if (curr_diameter < 1) return false;

    for (auto unb : Triangles) {

		if (!g->Triangle_balance[unb.id] && !g->Triangle_broken[unb.id]){
        // if (!unb.is_balanced && !unb.is_broken) {
			bool runned = false;
			// cout<<"check1\n";

            if (!runned && hop_num[all_edge_pairs[unb.edge1].v1] == hop_num[all_edge_pairs[unb.edge1].v2]) {
				runned = true;
				// cout<<"check2\n";
				if (!delete_on_edge(unb.edge1, g,false)) return false;
                // delete_on_edge(unb.edge1, g,false);

            } 
			
			else if (!runned && hop_num[all_edge_pairs[unb.edge2].v1] == hop_num[all_edge_pairs[unb.edge2].v2]) {
				runned = true;
				// cout<<"check3\n";

                // delete_on_edge(unb.edge2, g,false);
				if (!delete_on_edge(unb.edge2, g,false)) return false;
            } 

			else if (!runned && hop_num[all_edge_pairs[unb.edge3].v1] == hop_num[all_edge_pairs[unb.edge3].v2]) {

				runned = true;
				// cout<<"check4\n";
                // delete_on_edge(unb.edge3, g,false);
                if (!delete_on_edge(unb.edge3, g,false)) return false;

            }
			// cout<<"check=5\n";


        }
    }

	
	return delete_on_radius(g);;
}



bool removeNegativeTriangle(Graph* g) {
	int curr_diameter = g->diameter;
	if (curr_diameter < 1) return false;
	if (g->unbalance_num <= 0) return true;
	cout<<"===remove unbalanced==\n"<<endl;	


	queue<Graph*> Queue;
	Graph *best_g = new Graph();
	*best_g = *g;
	// g->diameter ++;
	Queue.push(best_g);
	int que_in,que_out = 0;
	int count  = 0;
	Graph *temp_g = new Graph();
	*temp_g = *g;
	// cout<<"1\n";
    if (!quickremoveNegativeTriangle(temp_g)) {
		cout << "quicl delete fail\n";
		*temp_g = *g;
		bool check = false;
		// cout<<"1.1\n";

		for (auto unb : Triangles) {
			// if (!unb.is_balanced && !unb.is_broken) {
			if (!temp_g->Triangle_balance[unb.id] && !temp_g->Triangle_broken[unb.id]){
			
				if (!check && delete_on_edge(unb.edge1, temp_g,false)) check = true;
				if (!check && delete_on_edge(unb.edge2, temp_g,false)) check = true;
				if (!check && delete_on_edge(unb.edge2, temp_g,false)) check = true;
				break;
			}
		}
		if (!check) {
			cout<< "cannot delete \n";
			delete(temp_g);
			temp_g = NULL;

			return false;
		}

	}
	// cout<<"2\n";

	*g = *temp_g;
	delete(temp_g);
	temp_g = NULL;

	findLongestPath(g);
	// cout<<"3\n";

	while (!Queue.empty()) {
		
		int size = Queue.size();
		Graph *first = Queue.front();
		Queue.pop();
		queue<Graph*> filteredQueue;

		que_out++;

		for (auto unb : Triangles) {
			// if (!unb.is_balanced && !unb.is_broken) {
			if (!first->Triangle_balance[unb.id] && !first->Triangle_broken[unb.id]){


				Graph *left = new Graph();
				*left = *first;
				Graph *middle = new Graph();
				*middle = *first;

				Graph *right = new Graph();
				*right = *first;
				int max_value = 0;
				Graph *max_graph;
				// cout<<"4\n";

				if (delete_on_edge(unb.edge1, left,true) && (left->size_of_truss- left->unbalance_num)> g->size_of_truss) {
					// bool best = false;
					// cout<<"4.1\n";
					// cout<< "checkkkk11\n";

					if ((left->diameter == g->diameter && left->unbalance_num == 0
					&& g->size_of_truss < left->size_of_truss) ||
					left->diameter < g->diameter && left->unbalance_num == 0) {
						*g = *left;
						// best = true;
						delete(left);
						left = NULL;


					} else {
						if (left->diameter >= global_hop && left->unbalance_num != 0) {
						// if (left->unbalance_num != 0) {
							// if (left->size_of_truss > max_value) 
							filteredQueue.push(left);
							que_in++;
						} else {
							delete(left);
							left = NULL;

							// if (!best) delete(left);
						}
					}
					// cout<<"4.2\n";


				} else {
					delete(left);
					left = NULL;

				}
				// cout<<"5\n";


				if (delete_on_edge(unb.edge2, middle,true) && (middle->size_of_truss-middle->unbalance_num) > g->size_of_truss) {
					// bool best = false;
					if ((middle->diameter == g->diameter && middle->unbalance_num == 0
					&& g->size_of_truss < middle->size_of_truss ) ||
					middle->diameter < g->diameter && middle->unbalance_num == 0)  {
						*g = *middle;
						delete(middle);
						middle = NULL;


						// best = true;
					} else {
						if (middle->diameter >= global_hop && middle->unbalance_num != 0) {
						// if (middle->unbalance_num != 0) {
							filteredQueue.push(middle);
							que_in++;

						} else {
							delete(middle);
							middle = NULL;

							// if (!best) delete(middle);
						}
					}


				} else {
					delete(middle);
					middle = NULL;

				}
				// cout<< "checkkkk3\n";

				if (delete_on_edge(unb.edge3, right,true) && (right->size_of_truss-right->unbalance_num) > g->size_of_truss) {
					// bool best = false;
					if ((right->diameter == g->diameter && right->unbalance_num == 0
					&& g->size_of_truss < right->size_of_truss ) ||
					right->diameter < g->diameter && right->unbalance_num == 0) {
						*g = *right;
						delete(right);
						right = NULL;

						// best = true;
					} else {
						if (right->diameter >= global_hop && right->unbalance_num != 0) {
						// if (right->unbalance_num != 0) {
							filteredQueue.push(right);
							que_in++;
						} else {
							delete(right);
							right = NULL;

							// if (!best) delete(right);
						}
					}


				} else {
					delete(right);
					right = NULL;

				}
				// cout<< "checkkkk4\n";
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

						tempQueue.push(currentGraph);
					} 
					Queue = move(tempQueue);

				}

				// Now, filteredQueue contains only the subgraphs that are not subgraphs of any other subgraph.
				// You can continue using filteredQueue or copy its contents back to the original Queue.

				// Clean up the original Queue if needed.
				while (!filteredQueue.empty()) {
					Graph *temp = filteredQueue.front();
					filteredQueue.pop();
					delete(temp);
				}
                break;


			}
		}
		// if (count %10) {
		// 	cout<< "count size: "<< count<< endl;
		// 	cout<< "queue size: "<< Queue.size()<< endl;
		// }
		delete(first);
		first = NULL;

	}
    if (g->unbalance_num >0) return false;
	return true;

}


// bool removeNegativeTriangle(Graph* g) {
//     int curr_diameter = g->diameter;
//     if (curr_diameter < 1) return false;
//     if (g->unbalance_num <= 0) return true;
//     cout << "===remove unbalanced==\n" << endl;

//     std::priority_queue<Graph*, std::vector<Graph*>, std::function<bool(Graph*, Graph*)>> Queue(
//         [](Graph* a, Graph* b) -> bool {
//             return a->size_of_truss < b->size_of_truss; // Largest size_of_truss first
//         }
//     );

//     Graph* best_g = new Graph();
//     *best_g = *g;
//     Queue.push(best_g);
//     int que_in, que_out = 0;
//     Graph* temp_g = new Graph();
//     *temp_g = *g;

//     if (!quickremoveNegativeTriangle(temp_g)) {
//         cout << "quick delete fail\n";
//         *temp_g = *g;
//         bool check = false;
//         for (auto unb : Triangles) {
//             if (!unb.is_balanced && !unb.is_broken) {
//                 if (!check && delete_on_edge(unb.edge1, temp_g, false)) check = true;
//                 if (!check && delete_on_edge(unb.edge2, temp_g, false)) check = true;
//                 if (!check && delete_on_edge(unb.edge3, temp_g, false)) check = true;
//                 break;
//             }
//         }
//         if (!check) {
//             cout << "cannot delete \n";
//             delete (temp_g);
//             return false;
//         }
//     }
//     *g = *temp_g;
//     delete (temp_g);
//     findLongestPath(g);

//     while (!Queue.empty()) {
//         Graph* first = Queue.top(); // Use top instead of front
//         Queue.pop();

//         // If this graph has unbalance_num == 0 and its size_of_truss is largest, use it.
//         if (first->unbalance_num == 0) {
//             *g = *first;
//             return true;
//         }

//         std::priority_queue<Graph*, std::vector<Graph*>, std::function<bool(Graph*, Graph*)>> filteredQueue(
//             [](Graph* a, Graph* b) -> bool {
//                 return a->size_of_truss < b->size_of_truss;
//             }
//         );

//         que_out++;
//         // ... same code as before		
// 		for (auto unb : Triangles) {
// 			if (!unb.is_balanced && !unb.is_broken) {

// 				Graph *left = new Graph();
// 				*left = *first;
// 				Graph *middle = new Graph();
// 				*middle = *first;

// 				Graph *right = new Graph();
// 				*right = *first;
// 				int max_value = 0;
// 				Graph *max_graph;
// 				if (delete_on_edge(unb.edge1, left,true) && (left->size_of_truss- left->unbalance_num)> g->size_of_truss) {
// 					// bool best = false;
// 					if ((left->diameter == g->diameter && left->unbalance_num == 0
// 					&& g->size_of_truss < left->size_of_truss) ||
// 					left->diameter < g->diameter && left->unbalance_num == 0) {
// 						*g = *left;
// 						// best = true;
// 						delete(left);

// 					} else {
// 						if (left->diameter >= global_hop && left->unbalance_num != 0) {
// 						// if (left->unbalance_num != 0) {
// 							// if (left->size_of_truss > max_value) 
// 							filteredQueue.push(left);
// 							que_in++;
// 						} else {
// 							delete(left);
// 							// if (!best) delete(left);
// 						}
// 					}

// 				} else {
// 					delete(left);
// 				}
// 				if (delete_on_edge(unb.edge2, middle,true) && (middle->size_of_truss-middle->unbalance_num) > g->size_of_truss) {
// 					// bool best = false;
// 					if ((middle->diameter == g->diameter && middle->unbalance_num == 0
// 					&& g->size_of_truss < middle->size_of_truss ) ||
// 					middle->diameter < g->diameter && middle->unbalance_num == 0)  {
// 						*g = *middle;
// 						delete(middle);

// 						// best = true;
// 					} else {
// 						if (middle->diameter >= global_hop && middle->unbalance_num != 0) {
// 						// if (middle->unbalance_num != 0) {
// 							filteredQueue.push(middle);
// 							que_in++;

// 						} else {
// 							delete(middle);
// 							// if (!best) delete(middle);
// 						}
// 					}


// 				} else {
// 					delete(middle);
// 				}
// 				if (delete_on_edge(unb.edge3, right,true) && (right->size_of_truss-right->unbalance_num) > g->size_of_truss) {
// 					// bool best = false;
// 					if ((right->diameter == g->diameter && right->unbalance_num == 0
// 					&& g->size_of_truss < right->size_of_truss ) ||
// 					right->diameter < g->diameter && right->unbalance_num == 0) {
// 						*g = *right;
// 						delete(right);
// 						// best = true;
// 					} else {
// 						if (right->diameter >= global_hop && right->unbalance_num != 0) {
// 						// if (right->unbalance_num != 0) {
// 							filteredQueue.push(right);
// 							que_in++;
// 						} else {
// 							delete(right);
// 							// if (!best) delete(right);
// 						}
// 					}


// 				} else {
// 					delete(right);
// 				}

//                 break;


// 			}
// 		} 

//         while (!filteredQueue.empty()) {
//             Graph* currentGraph = filteredQueue.top();
//             filteredQueue.pop();

//             // ... same code as before ...

//             std::priority_queue<Graph*, std::vector<Graph*>, std::function<bool(Graph*, Graph*)>> tempQueue(
//                 [](Graph* a, Graph* b) -> bool {
//                     return a->size_of_truss < b->size_of_truss;
//                 }
//             );
// 			bool isSubgraphOfOther = false;
// 			// queue<Graph*> tempQueue; // Copy filteredQueue to a temporary queue.

// 			while (!Queue.empty()) {
// 				Graph* graphInFiltered = filteredQueue.top();
// 				Queue.pop();

// 				if (!isSubgraphOfOther && isSubgraph(currentGraph, graphInFiltered)) {
// 					isSubgraphOfOther = true;
// 					// cout << "prunning queue\n";
// 					delete (currentGraph);
// 				}
// 				tempQueue.push(graphInFiltered);
// 			}

// 			// If currentGraph is not a subgraph of any other subgraph in filteredQueue, keep it.
// 			if (!isSubgraphOfOther) {
// 				tempQueue.push(currentGraph);
// 			} 


//             while (!tempQueue.empty()) {
//                 Queue.push(tempQueue.top());
//                 tempQueue.pop();
//             }
//         }
//     }

//     // If no suitable graph found
//     if (g->unbalance_num > 0) return false;
//     return true;
// }


void GetmaximumKtruss(Graph *g) {
	g->size_of_truss = e_num;
	int counts = 0;
	int book[MAX_V];
	vector<int> is_booked;
	int two_dimension[MAX_E];
	memset(g->support, 0, sizeof(g->support));
	Triangles.clear();


	for (int i = 0; i < e_num; i++) {

		Triangle temp_triangle;


		temp_triangle.edge1 = i;

		for (int j = 0; j < adj_edge[all_edge_pairs[i].v1].size(); j++) {
			int edgeIndex = adj_edge[all_edge_pairs[i].v1][j];
			int v = (all_edge_pairs[edgeIndex].v1 == all_edge_pairs[i].v1)
					? all_edge_pairs[edgeIndex].v2
					: all_edge_pairs[edgeIndex].v1;

			// If the vertex v is not within n hops of the source, skip this iteration
			// if (hop_num[v] > n)
			// 	continue;

			counts++;
			book[v] = 1;
			two_dimension[v] = j;
			is_booked.push_back(v);
		}
		for (int j = 0; j < adj_edge[all_edge_pairs[i].v2].size(); j++) {
			int edg2 = adj_edge[all_edge_pairs[i].v2][j];
			int v = (all_edge_pairs[edg2].v1 == all_edge_pairs[i].v2)
					? all_edge_pairs[edg2].v2
					: all_edge_pairs[edg2].v1;

			// If the vertex v is not within n hops of the source, skip this iteration
			// if (hop_num[v] > n)
			// 	continue;

			counts++;
			if (book[v]) {
				int edg1 = adj_edge[all_edge_pairs[i].v1][two_dimension[v]]; //?edg1
				if (edg1 > i && edg2 > i) {
					temp_triangle.edge2 = edg1;
					temp_triangle.edge3 = edg2;
					// temp_triangle.v3 = v;
					// temp_triangle.id = Triangles.size();
					
					if ((all_edge_pairs[edg1].sign + all_edge_pairs[edg2].sign + all_edge_pairs[i].sign) == -1) {
						temp_triangle.is_balanced = 1;
						// g->Triangle_balance[temp_triangle.id] = true;
					}
					else if ((all_edge_pairs[edg1].sign + all_edge_pairs[edg2].sign + all_edge_pairs[i].sign) == 3) {
						temp_triangle.is_balanced = 1;
						// g->Triangle_balance[temp_triangle.id] = true;

					}
				
					Triangles.push_back(temp_triangle);
					in_which_triangle[i].push_back(Triangles.size() - 1);
					in_which_triangle[edg1].push_back(Triangles.size() - 1);
					in_which_triangle[edg2].push_back(Triangles.size() - 1);
				}
			}
		}

		for (int j = 0; j < is_booked.size(); j++) {
			book[is_booked[j]] = 0;
		}
		is_booked.clear();
	}

	for (int i = 0; i < Triangles.size(); i++) {
		Triangles[i].id = i;
		g->Triangle_balance[i] = Triangles[i].is_balanced;
		g->Triangle_broken[i] = Triangles[i].is_broken;

	}

	for (int i = 0; i < Triangles.size(); i++) {
		if (!g->Triangle_balance[Triangles[i].id])
			g->unbalance_num++;
		else {
			g->support[Triangles[i].edge1]++;
			g->support[Triangles[i].edge2]++;
			g->support[Triangles[i].edge3]++;
		}
		if (Triangles[i].is_balanced != g->Triangle_balance[Triangles[i].id]) {
			cout << "error\n";
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
				// g->temp_delete_e[i] = 1;
				q.push(i);
			}
		}
	}

	while (!q.empty()) {
		int sub = q.front();
		q.pop();
		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
		for (int i = 0; i < in_which_triangle[sub].size(); i++) {

			if (!g->Triangle_broken[in_which_triangle[sub][i]]){
			// if(!g->Triangle_broken[Triangles[in_which_triangle[sub][i]].id]){
				
				if (g->Triangle_balance[Triangles[in_which_triangle[sub][i]].id]) {
					// 删除临边
					g->support[Triangles[in_which_triangle[sub][i]].edge1]--;
					g->support[Triangles[in_which_triangle[sub][i]].edge2]--;
					g->support[Triangles[in_which_triangle[sub][i]].edge3]--;
				}
				else 
					g->unbalance_num--;
				// 删除一条边后check 他的临边

				if (!g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1]) {
					if (g->support[Triangles[in_which_triangle[sub][i]].edge1] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge1);
						g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						// g->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						g->size_of_truss--;
					}
				}

				if (!g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2]) {
					if (g->support[Triangles[in_which_triangle[sub][i]].edge2] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge2);
						g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						// g->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						g->size_of_truss--;
					}
				}

				if (!g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3]) {
					if (g->support[Triangles[in_which_triangle[sub][i]].edge3] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge3);
						g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						// g->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						g->size_of_truss--;
					}
				}

				g->Triangle_broken[Triangles[in_which_triangle[sub][i]].id] = true;

			}

		}
	}

	// cout << "size of KTruss" << g->size_of_truss << endl;
	// cout << "truss unb num:" << g->unbalance_num << endl;






}


auto compare = [](const Graph* a, const Graph* b) {
	// cout << "123\n";
	// cout << a->diameter << "\n";
    if (a->diameter == b->diameter) {
		// cout << "======================================\n";
		// cout << a->diameter << "\n";
        return a->size_of_truss > b->size_of_truss;
    } 
    return a->diameter < b->diameter;
};

bool GetKtruss(int src, int k, Graph* g) {

	// check if the original graph has k-truss including src

	cout<<"check11\n";

	GetmaximumKtruss(g);
	if (!if_query_inside(g)) {
		return false;
	}

	cout<<"check12\n";

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
		for (int i = 0; i < adj_edge[curr_v].size(); i++) {
			int edgeIndex = adj_edge[curr_v][i];
			int v = (all_edge_pairs[edgeIndex].v1 == curr_v)
					? all_edge_pairs[edgeIndex].v2
					: all_edge_pairs[edgeIndex].v1;

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
	cout<<"check13\n";

	while (curr_hop >= 1 && curr_hop <= max_hop) {
		// cout << "--------when hop is " << curr_hop << "\n";
		Graph *g_hop = new Graph();
		g_hop =  GetKtrusswith_Nhops(curr_hop, k, g);
		
		if (g_hop->size_of_truss > 0 && if_query_inside(g_hop)) {

			global_hop = curr_hop;

			cout << "===================Deleting diameter ==================\n";
			bool if_nega_change = false;
			findLongestPath(g_hop);
			bool negative = false;	
			// print_result(g_hop);
			cout<<"check14\n";

			if (removeEdgeFromLongestPath(g_hop) && g_hop->diameter <= global_hop) {
				// cout<<global_hop <<"  test\n";
				cout<<"check15\n";

				while (!graphPtrs.empty()) {
					cout << "candidate:"<<graphPtrs.size() << "\n";
					sort(graphPtrs.begin(), graphPtrs.end(), compare);
					// cout<<"test1\n";
					// cout <<"dia1:"<< g_hop->diameter << "\n";
					// cout <<"dia2:"<< graphPtrs[0]->diameter << "\n";
					// for (Graph* graph : graphPtrs) {
					// 	cout <<"size:"<< graph->size_of_truss << "\n";

					// }
					cout<<"check17\n";

					if (graphPtrs[0]->unbalance_num > 0) {
						// cout <<"dia1:"<< graphPtrs[0]->size_of_truss << "\n";
						// cout <<"dia1:"<< graphPtrs[3]->size_of_truss << "\n";
						
						if (!removeNegativeTriangle(graphPtrs[0]) || graphPtrs[0]->diameter > global_hop) {
							// cout<<"test2\n";
							if (!graphPtrs.empty()) {
								delete graphPtrs[0];  // 释放动态分配的内存
								graphPtrs.erase(graphPtrs.begin());  // 删除vector中的第一个元素
							}
						}
						// cout <<"dia1:"<< graphPtrs[0]->size_of_truss << "\n";


					} else {
						*g = *graphPtrs[0];
						return true;
					}
				}
				cout<<"check16\n";

				
			}

		} 		
		hop_count:
		for (Graph* graph : graphPtrs) {
			delete graph;  // 释放动态分配的内存
		}
		graphPtrs.clear(); 
		curr_hop += 1;
		delete(g_hop);
	}
	cout<<"---calculation fail! \n"<<endl;
	return false;
}






int main(int argc, char** argv) {
	if (argc < 4) {
		printf("Not enough args: ./ade start_vertex k filename\n");
		return 0;
	}

    start_vertex = atoi(argv[1]);
	k = atoi(argv[2]);
    filename = argv[3];

	outname = filename + "_ade_solution.txt";
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
   	cout << "base_line time: " << timeuse << '\n' << endl;
	if (result) print_result(g);

	return 0;
}
