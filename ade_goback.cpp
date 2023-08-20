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
	vector<bool> Triangle_balance{vector<bool>(5000000,false)};
	vector<bool> Triangle_broken{vector<bool>(5000000,false)};
	int support[MAX_E];
	int edge_num;
	int unbalance_num = 0;
	// vector<Triangle> Triangles;
	int size_of_truss;

	int diameter;

} Graph;




clock_t start, finish;
int start_vertex, k, global_hop;


// modify from here
int hop_num[MAX_V];
int node_deleted[MAX_V];
int visited[MAX_V];

string filename, outname;

bool delete_on_radius(Graph *g_hop);
vector<int> findLongestDistanceFromStartVertex(Graph *g);

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
	findLongestDistanceFromStartVertex(g); 

	queue<int> q;
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
		for (int i = 0; i < adj_edge[curr_v].size(); i++) {
			int edgeIndex = adj_edge[curr_v][i];
			if (g->is_delete_e[edgeIndex]) continue;

			int v = (all_edge_pairs[edgeIndex].v1 == curr_v)
					? all_edge_pairs[edgeIndex].v2
					: all_edge_pairs[edgeIndex].v1;

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
		
	}


	g->diameter = pathLength;
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
	delete_on_radius(newG);
	if (newG->size_of_truss > 0 && if_query_inside(newG)) {
		if (update_dia) findLongestDistanceFromStartVertex(newG);
		return true;
	} 

	return false;
}
vector<int> check_follower(int edge_num, Graph*newG){
	vector<int> delCost_vertices;
    //因为这个计算delCost不能破坏原本的值，所以就用了copy
    vector<bool> is_delete_e_copy = newG->is_delete_e;
	int support_copy[MAX_E];
	copy(newG->support, newG->support + MAX_E, support_copy);
	vector<bool> Triangle_broken_copy = newG->Triangle_broken;


	queue<int> q;


	if (!is_delete_e_copy[edge_num] ) {
		is_delete_e_copy[edge_num] = 1;
		q.push(edge_num);
		delCost_vertices.push_back(edge_num);
	}

	// cout<<"11\n";

	while (!q.empty()) {
		int sub = q.front();
		q.pop();
		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
		for (int i = 0; i < in_which_triangle[sub].size(); i++) {
			if (!Triangle_broken_copy[Triangles[in_which_triangle[sub][i]].id]){
				if (newG->Triangle_balance[Triangles[in_which_triangle[sub][i]].id]) {
					// 删除临边
					support_copy[Triangles[in_which_triangle[sub][i]].edge1]--;
					support_copy[Triangles[in_which_triangle[sub][i]].edge2]--;
					support_copy[Triangles[in_which_triangle[sub][i]].edge3]--;
				}
				else 
				// 删除一条边后check 他的临边
				if (!is_delete_e_copy[Triangles[in_which_triangle[sub][i]].edge1]) {
					if (support_copy[Triangles[in_which_triangle[sub][i]].edge1] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge1);
						delCost_vertices.push_back(Triangles[in_which_triangle[sub][i]].edge1);

						is_delete_e_copy[Triangles[in_which_triangle[sub][i]].edge1] = 1;
					}
				}
				if (!is_delete_e_copy[Triangles[in_which_triangle[sub][i]].edge2]) {
					if (support_copy[Triangles[in_which_triangle[sub][i]].edge2] < k - 2)
					{
						delCost_vertices.push_back(Triangles[in_which_triangle[sub][i]].edge2);
						q.push(Triangles[in_which_triangle[sub][i]].edge2);
						is_delete_e_copy[Triangles[in_which_triangle[sub][i]].edge2] = 1;
					}
				}
				if (!is_delete_e_copy[Triangles[in_which_triangle[sub][i]].edge3]) {
					if (support_copy[Triangles[in_which_triangle[sub][i]].edge3] < k - 2)
					{
						delCost_vertices.push_back(Triangles[in_which_triangle[sub][i]].edge3);
						q.push(Triangles[in_which_triangle[sub][i]].edge3);
						is_delete_e_copy[Triangles[in_which_triangle[sub][i]].edge3] = 1;
					}
				}
				Triangle_broken_copy[Triangles[in_which_triangle[sub][i]].id] = true;
			}
		}
	}
	is_delete_e_copy.clear();
	Triangle_broken_copy.clear();

	return delCost_vertices;
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

    for (auto unb : Triangles) {

		if (!g->Triangle_balance[unb.id] && !g->Triangle_broken[unb.id]){
        // if (!unb.is_balanced && !unb.is_broken) {
			bool runned = false;
			// cout<<"check1\n";

            if (!runned && hop_num[all_edge_pairs[unb.edge1].v1] == hop_num[all_edge_pairs[unb.edge1].v2]) {
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
			
			else if (!runned && hop_num[all_edge_pairs[unb.edge2].v1] == hop_num[all_edge_pairs[unb.edge2].v2]) {
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

			else if (!runned && hop_num[all_edge_pairs[unb.edge3].v1] == hop_num[all_edge_pairs[unb.edge3].v2]) {

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

auto compare = [](const Graph* a, const Graph* b) {
	// cout << "123\n";
	// cout << a->diameter << "\n";
    if (a->diameter == b->diameter) {
		// cout << "======================================\n";
		// cout << a->diameter << "\n";
        return a->size_of_truss > b->size_of_truss;
    } 
    return a->diameter > b->diameter;
};
int quick_delete(Graph* g, Triangle unb) {
	if (hop_num[all_edge_pairs[unb.edge1].v1] == hop_num[all_edge_pairs[unb.edge1].v2]) {
		return unb.edge1;
	} else if (hop_num[all_edge_pairs[unb.edge2].v1] == hop_num[all_edge_pairs[unb.edge2].v2]) {
			return unb.edge2;
	} else return unb.edge3;

}




bool removeNegativeTriangle(Graph* g) {
	int curr_diameter = g->diameter;
	if (curr_diameter < 1) return false;
	cout<<"test1\n";
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
		cout<< "cannot delete \n";
		delete(temp_g);
		temp_g = NULL;
		return false;
	}
	// cout<<"2\n";

	*g = *temp_g;
	delete(temp_g);
	temp_g = NULL;

	findLongestDistanceFromStartVertex(g);
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
				int check = quick_delete(first, unb);
				vector<int> v1 = check_follower(check,first);
				cout<< "size1:" <<v1.size()<<endl;
			
				if (v1.size() == 1) {
					Graph *left = new Graph();
					*left = *first;
					// cout<<"4.1\n";

					if (delete_on_edge(check, left,false) && delete_on_radius(left)) {
						// bool best = false;
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
					
				} else {
					Graph *left = new Graph();
					*left = *first;
					Graph *middle = new Graph();
					*middle = *first;

					Graph *right = new Graph();
					*right = *first;
					int max_value = 0;
					Graph *max_graph;
					int quick_d = quick_delete(left, unb);
					// cout<<"4\n";
					if (delete_on_edge(unb.edge1, left,false) && delete_on_radius(left)) {
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


					if (delete_on_edge(unb.edge2, middle,false)&&delete_on_radius(middle)) {
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

					if (delete_on_edge(unb.edge3, right,false)&&delete_on_radius(right))  {
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
							currentGraph = NULL;

						}
						tempQueue.push(graphInFiltered);
					}

					// If currentGraph is not a subgraph of any other subgraph in filteredQueue, keep it.
					if (!isSubgraphOfOther) {
						tempQueue.push(currentGraph);
						count++;
						// cout<< "count size: "<< count<< endl;

						// Queue.push(currentGraph);
					} 
					while (!tempQueue.empty()) {
						Queue.push(tempQueue.front());
						tempQueue.pop();
					}
				}
                break;


			}
		}
		// if (count %10) {

		// }
		delete(first);
		first = NULL;

	}
	cout<< "q-in size: "<< que_in << endl;
	cout<< "queue out: "<< que_out + Queue.size()<< endl;
    if (g->unbalance_num >0) return false;
	return true;

}





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
	// cout<<"1\n";

	while (!q.empty()) {
		int sub = q.front();
		q.pop();
		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
		for (int i = 0; i < in_which_triangle[sub].size(); i++) {
			// cout<<"1.0\n";

			if (!g->Triangle_broken[in_which_triangle[sub][i]]){
			// if(!g->Triangle_broken[Triangles[in_which_triangle[sub][i]].id]){
				// cout<<"1.1\n";
				
				if (g->Triangle_balance[Triangles[in_which_triangle[sub][i]].id]) {
					// 删除临边
					g->support[Triangles[in_which_triangle[sub][i]].edge1]--;
					g->support[Triangles[in_which_triangle[sub][i]].edge2]--;
					g->support[Triangles[in_which_triangle[sub][i]].edge3]--;
				}
				else 
					g->unbalance_num--;
				// 删除一条边后check 他的临边
				// cout<<"1.2\n";

				if (!g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1]) {
					if (g->support[Triangles[in_which_triangle[sub][i]].edge1] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge1);
						g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						// g->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						g->size_of_truss--;
					}
				}
				// cout<<"1.3\n";

				if (!g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2]) {
					if (g->support[Triangles[in_which_triangle[sub][i]].edge2] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge2);
						g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						// g->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						g->size_of_truss--;
					}
				}
				// cout<<"1.4\n";

				if (!g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3]) {
					if (g->support[Triangles[in_which_triangle[sub][i]].edge3] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge3);
						g->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						// g->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						g->size_of_truss--;
					}
				}
				// cout<<"1.5\n";

				g->Triangle_broken[Triangles[in_which_triangle[sub][i]].id] = true;
				// cout<<"1.6\n";

			}
			// cout<<"1.7\n";

		}
	}
	// cout<<"2\n";

	// cout << "size of KTruss" << g->size_of_truss << endl;
	// cout << "truss unb num:" << g->unbalance_num << endl;






}


// bool delete_on_edge(int edge_num, Graph*newG, int update_dia) {


// 	queue<int> q;


// 	if (!newG->is_delete_e[edge_num] ) {
// 		newG->is_delete_e[edge_num] = 1;
// 		goBack1.push_back(edge_num);
// 		newG->size_of_truss--;
// 		del_count++;
// 		// cout<<"11\n";
// 		// newG->temp_delete_e[edge_num] = 1;
// 		q.push(edge_num);
// 	}

// 	// cout<<"11\n";

// 	while (!q.empty()) {
// 		int sub = q.front();
// 		q.pop();
// 		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
// 		for (int i = 0; i < in_which_triangle[sub].size(); i++) {
// 			if (!newG->Triangle_broken[Triangles[in_which_triangle[sub][i]].id]){
// 				if (newG->Triangle_balance[Triangles[in_which_triangle[sub][i]].id]) {
// 					// 删除临边
// 					newG->support[Triangles[in_which_triangle[sub][i]].edge1]--;
// 					newG->support[Triangles[in_which_triangle[sub][i]].edge2]--;
// 					newG->support[Triangles[in_which_triangle[sub][i]].edge3]--;
// 				}
// 				else 
// 					newG->unbalance_num--;
// 				// 删除一条边后check 他的临边
// 				if (!newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1]) {
// 					if (newG->support[Triangles[in_which_triangle[sub][i]].edge1] < k - 2)
// 					{
// 						q.push(Triangles[in_which_triangle[sub][i]].edge1);
// 						newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
// 						// newG->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
// 						newG->size_of_truss--;
// 					}
// 				}
// 				if (!newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2]) {
// 					if (newG->support[Triangles[in_which_triangle[sub][i]].edge2] < k - 2)
// 					{
// 						q.push(Triangles[in_which_triangle[sub][i]].edge2);
// 						newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
// 						// newG->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
// 						newG->size_of_truss--;
// 					}
// 				}
// 				if (!newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3]) {
// 					if (newG->support[Triangles[in_which_triangle[sub][i]].edge3] < k - 2)
// 					{
// 						q.push(Triangles[in_which_triangle[sub][i]].edge3);
// 						newG->is_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
// 						// newG->temp_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
// 						newG->size_of_truss--;
// 					}
// 				}
// 				newG->Triangle_broken[Triangles[in_which_triangle[sub][i]].id] = true;
// 			}
// 		}
// 	}
// 	// cout<<"12\n";


// 	// if (newG->unbalance_num == 0 && newG->size_of_truss > 0){
// 	// 	return delete_on_radius(newG);

// 	// }
// 	// cout<<"13\n";

// 	if (newG->size_of_truss > 0 && if_query_inside(newG)) {
// 		if (update_dia) findLongestDistanceFromStartVertex(newG);
// 		return true;
// 	} 

// 	return false;
// }

int exact(Graph*g,int maxSize) {
	int process_id = -11;//当前层需要处理的negtriangle的编号
	Triangle temp;
	for (auto unb : Triangles) {
		// if (!unb.is_balanced && !unb.is_broken) {
		if (!g->Triangle_balance[unb.id] && !g->Triangle_broken[unb.id]){
			process_id = unb.id;
			temp = unb;
			break;
		}
	}
	if(process_id == -11)
		return g->size_of_truss;


	int id1 = temp.edge1;//这个负三角形的三个顶点的id
	Graph *temp_g = new Graph();
	*temp_g = *g; 
	// del(id1, k, del_count, goBack1, goBack2);
	delete_on_edge(id1, g, false);
	maxSize = max(maxSize, exact(g, maxSize));
	*g = *temp_g;
	// delete(temp);


	int id2 = temp.edge2;
	// Graph *temp = new Graph();
	// *temp = *g; 
	// del(id1, k, del_count, goBack1, goBack2);
	delete_on_edge(id2, g, false);
	maxSize = max(maxSize, exact(g, maxSize));
	*g = *temp_g;
	// delete(temp);



	int id3 = temp.edge3;
	// Graph *temp = new Graph();
	// *temp = *g; 
	// // del(id1, k, del_count, goBack1, goBack2);
	delete_on_edge(id3, g, false);
	maxSize = max(maxSize, exact(g, maxSize));
	*g = *temp_g;
	delete(temp_g);

	return maxSize;

}

int exact_removeNegativeTriangle(Graph* g, int maxSize) {
    int process_id = -11; // Current triangle to process
    Triangle temp;

    for (auto unb : Triangles) {
        if (!g->Triangle_balance[unb.id] && !g->Triangle_broken[unb.id]) {
            process_id = unb.id;
            temp = unb;
            break;
        }
    }

    if (process_id == -11) {
        return g->size_of_truss;
    }

    // Store the current state of the graph before trying each edge removal.
    Graph *temp_g = new Graph();
    *temp_g = *g; 

    for(int i = 0; i < 3; i++) {
        // Select which edge to delete based on the iteration
        int edgeToDelete;
        switch(i) {
            case 0:
                edgeToDelete = temp.edge1;
                break;
            case 1:
                edgeToDelete = temp.edge2;
                break;
            case 2:
                edgeToDelete = temp.edge3;
                break;
        }
		int result = 0;
        if(delete_on_edge(edgeToDelete, g, false) && delete_on_radius(g)) {
			// PRUNING rules:
			// We'll only proceed with this graph state if the current state is not worse than the previously best
			if (g->diameter >= global_hop) {
				result = exact_removeNegativeTriangle(g, maxSize);
				
			}
			// We'll only update maxSize if the result is strictly better than the previous maxSize
		}
		if (result > maxSize) {
			maxSize = result;
		}
        *g = *temp_g; // Reset the graph to the stored state
    }

    delete(temp_g);

    return maxSize;
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
		int curr_v = q.front().first;
		int curr_hop = q.front().second;
		q.pop();
		for (int i = 0; i < adj_edge[curr_v].size(); i++) {
			int edgeIndex = adj_edge[curr_v][i];
			int v = (all_edge_pairs[edgeIndex].v1 == curr_v)
					? all_edge_pairs[edgeIndex].v2
					: all_edge_pairs[edgeIndex].v1;
			if (g->is_delete_e[edgeIndex]) continue;

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
		cout << "--------when hop is " << curr_hop << "\n";
		Graph *g_hop = new Graph();
		g_hop =  GetKtrusswith_Nhops(curr_hop, k, g);
		int maxSize = 0;
		if (g_hop->size_of_truss > 0 && if_query_inside(g_hop)) {

			global_hop = curr_hop;
			cout << "===================Deleting NegativeTriangle ==================\n";
            cout<<"check1\n";
			// if (!removeNegativeTriangle(g_hop))	goto hop_add;
			// maxSize = exact(g_hop,0);
            maxSize = exact_removeNegativeTriangle(g_hop,0);
			cout<<"check2\n";
			cout << "maxSize : " << maxSize << endl;
            // vector<int> longest_list = findLongestDistanceFromStartVertex(g_hop);

			// cout<<"---start to delete node \n"<<endl;
			// while (g_hop->diameter > (curr_hop)) {
			// 	for (auto i : longest_list) {
			// 		if ( !delete_on_node(i,g_hop)) goto hop_add;
			// 	}
			// 	findLongestDistanceFromStartVertex(g_hop);
			// }
            // cout<<"check3\n";

			// *g = *g_hop;
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
   	cout << "advanced time: " << timeuse << '\n' << endl;
	// if (result) print_result(g);

	return 0;
}
