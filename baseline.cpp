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
using namespace std;
const int MAX_V = 1160000;
const int MAX_E = 3000000;
int e_num, v_num;

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
	int is_delete_e[MAX_E];
	int link[MAX_E];
	int temp_delete_e[MAX_E];
	int break_unb[MAX_E];
	int support[MAX_E];


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
} Graph;

Graph *g;


clock_t start, finish;
int k;
// Edge all_edge_pairs[MAX_E];
// vector<int> vec[MAX_V];
// vector<int> adj_edge[MAX_V];
// //vector<int> edge_to_edge[MAX_E];
// vector<int> in_which_triangle[MAX_E];
// vector<int> is_linked;
// int book[MAX_V];
// vector<int> is_booked;
// int is_delete_e[MAX_E];
// int link[MAX_E];
// int temp_delete_e[MAX_E];
// int break_unb[MAX_E];
// int support[MAX_E];

// string filename;
// int unbalance_num = 0;
// vector<Triangle> Triangles;
// int k, size_of_truss;
// vector<int> followers[MAX_E + 1];
// vector<int> two_dimension;
// int grouped[MAX_E];
// int candidate[MAX_E];
// vector<int> candidates;
// vector<Group> groups;
// vector<int> id;

// modify from here
int hop_num[MAX_V];
int visited[MAX_V];

string filename;



void build_graph() {
	g = new Graph();
	memset(g->is_delete_e, 0, sizeof(g->is_delete_e));
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
}


bool GetKtrusswith_Nhops(int n, int k) {
	memset(g->is_delete_e, 0, sizeof(g->is_delete_e));
	memset(g->support, 0, sizeof(g->support));
	memset(g->book, 0, sizeof(g->book));
	memset(g->break_unb, 0, sizeof(g->break_unb));
	memset(g->temp_delete_e, 0, sizeof(g->temp_delete_e));
	memset(g->link, 0, sizeof(g->link));
	memset(g->grouped, -1, sizeof(g->grouped));
	memset(g->candidate, 0, sizeof(g->candidate));
	g->two_dimension.resize(MAX_E);
	
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
			if (g->support[i] < k - 2|| hop_num[g->all_edge_pairs[i].v1] > n || hop_num[g->all_edge_pairs[i].v2] > n
		|| hop_num[g->all_edge_pairs[i].v1] == -1 || hop_num[g->all_edge_pairs[i].v2] == -1) {
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

	return g->size_of_truss > 0;
	/*
	if (size_of_truss > 0){
		return true;
	}
	else {
		return false;
	}*/



}







vector<int> findLongestPath() {
	// little optimization by replace v_num by v_num-1
	//int dist[MAX_V][MAX_V] = {MAX_E + 1};
	int pathLength = -1;
	vector<int> path;
	
	for (int i = 1; i < v_num - 1; i++) {
		memset(visited, 0, sizeof(visited));
		queue<pair<int,pair<int, vector<int>>>> q;
		vector<int> curr_path;
		curr_path.push_back(i);
		q.push(make_pair(i, make_pair(0, curr_path)));
		
		while (!q.empty()) {
			int size = q.size();
			for (size; size > 0; size--) {
				auto curr_pair = q.front();
				q.pop();
				if (visited[curr_pair.first]) continue;

				visited[curr_pair.first] = true;
				if (pathLength < curr_pair.second.first) {
					path = curr_pair.second.second;
					pathLength = curr_pair.second.first;
				}
				
				for (int idx = 0; idx < g->vec[curr_pair.first].size(); idx++) {
					int v = g->vec[curr_pair.first][idx];
					
					vector<int> nextPath = curr_pair.second.second;
					nextPath.push_back(v);
					q.push(make_pair(v, make_pair(curr_pair.second.first + 1, nextPath)));
					
				}
			}
			
		}
		break;
	}
	return path;
}

bool removeEdgeFromLongestPath() {
	vector<int> path = findLongestPath();
	cout << "Longest Path is :\n";
	for (auto v : path) {
		cout << v << " -> ";
	}
	cout << "\n";

	//TODO
	//remove edge from any side of the longest path

	return false;
}
void delete_on_edge(int A, int B){
	int edge_num = -1;
	for (int i = 0; i < g->adj_edge[A].size(); ++i) {
        int edge_index = g->adj_edge[A][i];
        // Check if this edge index appears in the list for v2
        for (int j = 0; j < g->adj_edge[B].size(); ++j) {
            if (g->adj_edge[B][j] == edge_index) {
                // If it does, return this edge index
                edge_num = edge_index;
            }
        }
    }

	queue<int> q;
	if (!g->is_delete_e[edge_num] ) {
		g->is_delete_e[edge_num] = 1;
		g->size_of_truss--;
		g->temp_delete_e[edge_num] = 1;
		q.push(edge_num);
	
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


}


bool removeNegativeTriangle() {
	return false;
}
void print_result(){
		// Assuming vertices_set is your set of vertices
	set<int> vertices_set;
	cout <<  " ====result graph: \n";
	// Collect all unique vertices from all_edge_pairs into the set
	for (int i = 0; i < e_num; i++) {
		if (!g->is_delete_e[i]) {
			vertices_set.insert(g->all_edge_pairs[i].v1);
			vertices_set.insert(g->all_edge_pairs[i].v2);
		}
	}

	// Print all vertices in the set
	for (const auto& vertex : vertices_set) {
		cout << vertex << " ";
	}
	cout << endl;

}

void GetKtruss(int src, int k) {

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

		if (GetKtrusswith_Nhops(curr_hop, k)) {
			print_result();
			return;
			// calculate result for current diameter
			if (!removeNegativeTriangle()) {
				curr_hop += 1;
				continue;
			}
			int curr_diameter = 2 * curr_hop;
			//assume current graph in g and result store in res1;
			while (removeEdgeFromLongestPath()) {
				curr_diameter -= 1;
				if (removeNegativeTriangle()) {
					// update res
					curr_diameter -= 1;
				} else {
					// if exist valid res, return res;
					break;
				}
			}
			 
			// if (removeEdgeFromLongestPath() && removeNegativeTriangle()) return; 
		} 		
		curr_hop += 1;
	}
}





void candidate_lemma() {
	
	for (int i = 0; i < g->Triangles.size(); i++) {
		if (!g->Triangles[i].is_broken) {
			if (!g->Triangles[i].is_balanced) {
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

void Get_result() {
	int counts = 0;
	while (g->unbalance_num > 0) {
		counts++;
		int min = MAX_E;
		for (int i = 0; i < g->groups.size(); i++)
		{
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

						}
					}
				}
				//if (break_unb[i] == 0)
					//groups[i].is_delete = 1;
				if (g->break_unb[i]>0)
				{
					if (g->followers[i].size() < g->followers[min].size())
						min = i;
					else if (g->followers[i].size() == g->followers[min].size()) {
						if (g->groups[i].typ_edge < g->groups[min].typ_edge)
							min = i;
					}
				}

				for (int j = 0; j < g->is_linked.size(); j++)
				{
					g->link[g->is_linked[j]] = 0;
					g->temp_delete_e[g->is_linked[j]] = 0;
				}
				g->is_linked.clear();
			}
			
			
		}

		if (min != MAX_E)
		{
			for (int i = 0; i < g->followers[min].size(); i++)
			{
				g->is_delete_e[g->followers[min][i]] = 1;
				g->temp_delete_e[g->followers[min][i]] = 1;
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
			//out << "sizeof truss" << size_of_truss << endl;
			//out << "unb num:" << unbalance_num << endl;
			//out << counts << " " << groups[min].typ_edge << endl;
			
		}
	
	}
	cout << "delete " << counts << " times" << endl;
	
}


int findSmallestDistance(int start_vertex, int end_vertex) {
    vector<bool> visited(MAX_V, false);
    vector<int> distance(MAX_V, numeric_limits<int>::max());

    queue<int> q;
    q.push(start_vertex);
    visited[start_vertex] = true;
    distance[start_vertex] = 0;

    while (!q.empty()) {
        int current_vertex = q.front();
        q.pop();

        for (int neighbor : g->vec[current_vertex]) {
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

// vector<int> getInducedSubgraph(int start_vertex, int k) {
//     vector<bool> visited(MAX_V, false);
//     vector<int> distance(MAX_V, numeric_limits<int>::max());
//     vector<int> subgraph;

//     queue<int> q;
//     q.push(start_vertex);
//     visited[start_vertex] = true;
//     distance[start_vertex] = 0;
//     subgraph.push_back(start_vertex);

//     while (!q.empty()) {
//         int current_vertex = q.front();
//         q.pop();

//         if (distance[current_vertex] == k) {
//             break; // Stop exploring further at distance k
//         }

//         for (int neighbor : vec[current_vertex]) {
//             if (!visited[neighbor]) {
//                 visited[neighbor] = true;
//                 distance[neighbor] = distance[current_vertex] + 1;
//                 q.push(neighbor);
//                 subgraph.push_back(neighbor);
//             }
//         }
//     }

//     return subgraph;
// }


vector<int> getInducedSubgraph(const set<int>& vertices) {
    vector<bool> visited(MAX_V, false);
    vector<int> induced_subgraph;

    for (int vertex : vertices) {
        visited[vertex] = true;
        induced_subgraph.push_back(vertex);
    }

    return induced_subgraph;
}

// int getDistance(int start_vertex, const set<int>& vertices) {
//     vector<bool> visited(MAX_V, false);
//     vector<int> distance(MAX_V, numeric_limits<int>::max());

//     queue<int> q;
//     q.push(start_vertex);
//     visited[start_vertex] = true;
//     distance[start_vertex] = 0;

//     while (!q.empty()) {
//         int current_vertex = q.front();
//         q.pop();

//         for (int neighbor : vec[current_vertex]) {
//             if (!visited[neighbor]) {
//                 visited[neighbor] = true;
//                 distance[neighbor] = distance[current_vertex] + 1;
//                 q.push(neighbor);

//                 if (vertices.count(neighbor) > 0) {
//                     return distance[neighbor];
//                 }
//             }
//         }
//     }

//     return -1;
// }

int findLargestSmallestDistance(int start_vertex, const set<int>& end_vertices) {
    int largest_smallest_distance = -1;

    for (int end_vertex : end_vertices) {
        vector<bool> visited(MAX_V, false);
        vector<int> distance(MAX_V, numeric_limits<int>::max());

        queue<int> q;
        q.push(start_vertex);
        visited[start_vertex] = true;
        distance[start_vertex] = 0;

        while (!q.empty()) {
            int current_vertex = q.front();
            q.pop();

            for (int neighbor : g->vec[current_vertex]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    distance[neighbor] = distance[current_vertex] + 1;
                    q.push(neighbor);

                    if (neighbor == end_vertex) {
                        if (distance[neighbor] > largest_smallest_distance) {
                            largest_smallest_distance = distance[neighbor];
                        }
                        break;
                    }
                }
            }
        }
    }

    return largest_smallest_distance;
}


// int getDiameter() {
//     int diameter = 0;

//     for (int v = 0; v < MAX_V; v++) {
//         for (int u = v + 1; u < MAX_V; u++) {
//             int distance = getDistance(v, { u });
//             if (distance > diameter) {
//                 diameter = distance;
//             }
//         }
//     }

//     return diameter;
// }



int main() {
	// string outname;
	// scanf("%d", &k);
	// cin >> filename;
	// outname = filename + "5_solution.txt";
	// filename += ".txt";
	// build_graph();
	// cout << "graph build" << endl;
	// start = clock();
	// GetKtruss();
	// finish = clock();
	// cout << "Timecost" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	// cout << "Truss Ready" << endl;
	// candidate_lemma();
	// finish = clock();
	// cout << "Timecost" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	// cout << "candidates ready???" << endl;
	// Get_result();
	// finish = clock();
	// cout << "size of truss" << size_of_truss << endl;
	// cout << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	// ofstream out;
	// out.open(outname, ios::out | ios::app);
	// if (out.fail())
	// 	cout << "error" << endl;
	// out << "size: " << size_of_truss << endl;
	// out << "Timecost" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	// out.close();

	string outname;
	// scanf("%d", &k);
	// cin >> filename;
    filename = "data/test1";
	outname = filename + "5_solution.txt";
	filename += ".txt";
	build_graph();
	cout << "graph build" << endl;

    // int start_vertex, k;
    // cout << "Enter the start vertex: ";
    // cin >> start_vertex;
    // cout << "Enter the distance threshold (k): ";
    // cin >> k;


    // vector<int> subgraph = getInducedSubgraph(start_vertex, k);
    // cout << "Vertices in the induced subgraph: ";
    // for (int vertex : subgraph) {
    //     cout << vertex << " ";
    // }
    // cout << endl;

    set<int> C;
    int start_vertex, k;
    cout << "Enter the start vertex: ";
    cin >> start_vertex;
	cout << "Enter the K: ";
    cin >> k;
	GetKtruss(start_vertex, k);

	/*
    C.insert(start_vertex); // Starting with vertex 0
    int d = 0;

    //while (d <= getDiameter()) {
    while (d <=2) {
        int Nc = C.size();

        for (int v : C) {
            for (int t : vec[v]) {
                if (C.count(t) > 0) {
                    continue;
                }

                int distance = findLargestSmallestDistance(t, C);
                if (distance <= d) {
                    C.insert(t);
                }
            }
        }

        if (C.size() == Nc) {
            d++;
        }
    }

    vector<int> induced_subgraph = getInducedSubgraph(C);

    // Display the induced subgraph
    cout << "Induced Subgraph: ";
    for (int vertex : induced_subgraph) {
        cout << vertex << " ";
    }
    cout << endl;
	*/
	// for (int i = 1; i < MAX_V; i++) {
	// 	if (hop_num[i] < 0) break;
	// 	cout << "hop number for " << i << " is " << hop_num[i] << "\n";
	// }

	return 0;
}