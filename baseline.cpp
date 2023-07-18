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

clock_t start, finish;
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

string filename;
int unbalance_num = 0;
vector<Triangle> Triangles;
int k, size_of_truss;
vector<int> followers[MAX_E + 1];
vector<int> two_dimension;
int grouped[MAX_E];
int candidate[MAX_E];
vector<int> candidates;
vector<Group> groups;
vector<int> id;

// modify from here
int hop_num[MAX_V];



void build_graph() {
	memset(is_delete_e, 0, sizeof(is_delete_e));
	memset(support, 0, sizeof(support));
	memset(book, 0, sizeof(book));
	memset(break_unb, 0, sizeof(break_unb));
	memset(temp_delete_e, 0, sizeof(temp_delete_e));
	memset(link, 0, sizeof(link));
	memset(grouped, -1, sizeof(grouped));
	memset(candidate, 0, sizeof(candidate));
	two_dimension.resize(MAX_E);
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
		all_edge_pairs[i].v1 = v1;
		all_edge_pairs[i].v2 = v2;
		all_edge_pairs[i].sign = sign;
		vec[v1].push_back(v2);
		vec[v2].push_back(v1);
		adj_edge[v1].push_back(i);
		adj_edge[v2].push_back(i);
	}
	for (int i = 0; i < MAX_E; i++) {
		followers[MAX_E].push_back(0);
	}
}
/*
void GetKtruss() {
	
	size_of_truss = e_num;
	int counts = 0;
	for (int i = 0; i < e_num; i++) {
		//if (i % 10000 == 0)
			//cout << i << endl;
		Triangle temp_triangle;
		// 对于每条边进行赋值
		temp_triangle.v1 = all_edge_pairs[i].v1;
		temp_triangle.v2 = all_edge_pairs[i].v2;
		temp_triangle.edge1 = i;
		for (int j = 0; j < vec[all_edge_pairs[i].v1].size(); j++) {
			counts++;
			int v = vec[all_edge_pairs[i].v1][j];
			//int edg1 = adj_edge[all_edge_pairs[i].v1][j];
			book[v] = 1;
			two_dimension[v] = j;
			is_booked.push_back(v);
		}
		for (int j = 0; j < vec[all_edge_pairs[i].v2].size(); j++) {
			counts++;
			int v = vec[all_edge_pairs[i].v2][j];
			int edg2 = adj_edge[all_edge_pairs[i].v2][j];
			if (book[v]) {
				int edg1 = adj_edge[all_edge_pairs[i].v1][two_dimension[v]];//?edg1
				//support[i]++;
				if (edg1 > i&&edg2 > i) {
					temp_triangle.edge2 = edg1;
					temp_triangle.edge3 = edg2;
					temp_triangle.v3 = v;
					
					if ((all_edge_pairs[edg1].sign + all_edge_pairs[edg2].sign + all_edge_pairs[i].sign) == -1) {
						temp_triangle.is_balanced = 1;
					}
					else if ((all_edge_pairs[edg1].sign + all_edge_pairs[edg2].sign + all_edge_pairs[i].sign) == 3) {
						temp_triangle.is_balanced = 1;
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
		if (!Triangles[i].is_balanced)
			unbalance_num++;
		else {
			support[Triangles[i].edge1]++;
			support[Triangles[i].edge2]++;
			support[Triangles[i].edge3]++;
		}
	}
	cout << "counts:" << counts << endl;
	//unbalance_num = Triangles.size() - num_of_balance;
	cout << "orangial:" << unbalance_num << endl;
	queue<int> q;
	// 找出所有不满足 support的边
	for (int i = 0; i < e_num; i++) {
		if (!is_delete_e[i]) {
			if (support[i] < k - 2) {
				is_delete_e[i] = 1;
				size_of_truss--;
				temp_delete_e[i] = 1;
				q.push(i);
			}
		}
	}
	while (!q.empty()) {
		int sub = q.front();
		q.pop();
		//in_which_triangle[sub][i]].edge1 表示包含 边SUB的第 i 个三角形的 三边
		for (int i = 0; i < in_which_triangle[sub].size(); i++) {
			if (!Triangles[in_which_triangle[sub][i]].is_broken) {
				if (Triangles[in_which_triangle[sub][i]].is_balanced) {
					// 删除临边
					support[Triangles[in_which_triangle[sub][i]].edge1]--;
					support[Triangles[in_which_triangle[sub][i]].edge2]--;
					support[Triangles[in_which_triangle[sub][i]].edge3]--;
				}
				else 
					unbalance_num--;
				// 删除一条边后check 他的临边
				if (!is_delete_e[Triangles[in_which_triangle[sub][i]].edge1]) {
					if (support[Triangles[in_which_triangle[sub][i]].edge1] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge1);
						is_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						temp_delete_e[Triangles[in_which_triangle[sub][i]].edge1] = 1;
						size_of_truss--;
					}
				}
				if (!is_delete_e[Triangles[in_which_triangle[sub][i]].edge2]) {
					if (support[Triangles[in_which_triangle[sub][i]].edge2] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge2);
						is_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						temp_delete_e[Triangles[in_which_triangle[sub][i]].edge2] = 1;
						size_of_truss--;
					}
				}
				if (!is_delete_e[Triangles[in_which_triangle[sub][i]].edge3]) {
					if (support[Triangles[in_which_triangle[sub][i]].edge3] < k - 2)
					{
						q.push(Triangles[in_which_triangle[sub][i]].edge3);
						is_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						temp_delete_e[Triangles[in_which_triangle[sub][i]].edge3] = 1;
						size_of_truss--;
					}
				}
				Triangles[in_which_triangle[sub][i]].is_broken = 1;
			}
		}
	}
	cout << "size of KTruss" << size_of_truss << endl;
	cout << "truss unb num:" << unbalance_num << endl;



}
*/

void GetKtruss(int src) {

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
		for (int i = 0; i < vec[curr_v].size(); i++) {
			int v = vec[curr_v][i];
			if (hop_num[v] == -1 || hop_num[v] > curr_hop + 1) {
				hop_num[v] = curr_hop + 1;
				
				//if (curr_hop == hop - 1) continue;
				if (max_hop < curr_hop + 1) max_hop = curr_hop + 1;
				q.push(make_pair(v, curr_hop + 1));
			}
		}
	}
	cout << "max hop is " << max_hop << "\n";
	/*
	for (int i = 0; i < MAX_V; i++) {
		if (hop_num[i] == 3) cout << i << "\n";
	}*/


	// add one hop at a time
	
	bool isSuccess = false;
	//int lo = 1;
	//int hi = max_hop;
	//hop = -1;
	int curr_hop = 1;
	while (curr_hop >= 1 && curr_hop <= max_hop) {
		
		if (getKtrussWithHops(hop)) {
			// should return success
			if (removeEdgeFromLongestPath() && removeNegativeTriangle()) return; 
		} 		
		curr_hop += 1;
	}
}

bool getKtrussWithHops(int hop) {
	return false;
}

bool removeEdgeFromLongestPath() {
	return false;
}

bool removeNegativeTriangle() {
	return false;
}


void candidate_lemma() {
	
	for (int i = 0; i < Triangles.size(); i++) {
		if (!Triangles[i].is_broken) {
			if (!Triangles[i].is_balanced) {
				if (!candidate[Triangles[i].edge1])
				{
					candidates.push_back(Triangles[i].edge1);
					candidate[Triangles[i].edge1] = 1;
				}
				if (!candidate[Triangles[i].edge2])
				{
					candidates.push_back(Triangles[i].edge2);
					candidate[Triangles[i].edge2] = 1;
				}
				if (!candidate[Triangles[i].edge3])
				{
					candidates.push_back(Triangles[i].edge3);
					candidate[Triangles[i].edge3] = 1;
				}
			}
		}
	}
	int num = 0;
	for (int i = 0; i < candidates.size(); i++) {
		Group temp_group;
		if (support[candidates[i]] == k - 2) {
			if (grouped[candidates[i]] == -1) {
				queue<int> q;
				q.push(candidates[i]);
				grouped[candidates[i]] = groups.size();
				temp_group.edges.push_back(candidates[i]);
				while (!q.empty()) {
					int sub = q.front();
					q.pop();
					for (int j = 0; j < in_which_triangle[sub].size(); j++) {
						int temp = in_which_triangle[sub][j];
						if (!Triangles[temp].is_broken) {
							if (Triangles[temp].is_balanced) {
								if (support[Triangles[temp].edge1] == k - 2 && grouped[Triangles[temp].edge1] == -1) {
									q.push(Triangles[temp].edge1);
									temp_group.edges.push_back(Triangles[temp].edge1);
									grouped[Triangles[temp].edge1] = groups.size();

								}
								if (support[Triangles[temp].edge2] == k - 2 && grouped[Triangles[temp].edge2] == -1) {
									q.push(Triangles[temp].edge2);
									temp_group.edges.push_back(Triangles[temp].edge2);
									grouped[Triangles[temp].edge2] = groups.size();

								}
								if (support[Triangles[temp].edge3] == k - 2 && grouped[Triangles[temp].edge3] == -1) {
									q.push(Triangles[temp].edge3);
									temp_group.edges.push_back(Triangles[temp].edge3);
									grouped[Triangles[temp].edge3] = groups.size();

								}
							}
						}
					}
				}
			}
		}
		else if (support[candidates[i]] > k - 2)
		{
			grouped[candidates[i]] = groups.size();
			temp_group.edges.push_back(candidates[i]);
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
			groups.push_back(temp_group);
		}
	}

}

void Get_result() {
	int counts = 0;
	while (unbalance_num > 0) {
		counts++;
		int min = MAX_E;
		for (int i = 0; i < groups.size(); i++)
		{
			break_unb[i] = 0;
			followers[i].clear();
			if (!groups[i].is_delete) {
				queue<int> q;
				for (int j = 0; j < groups[i].edges.size(); j++) {
					if (!is_delete_e[groups[i].edges[j]]) {
						q.push(groups[i].edges[j]);
						temp_delete_e[groups[i].edges[j]] = 1;
						followers[i].push_back(groups[i].edges[j]);
					}	
				}
				while (!q.empty()) {
					//counts++;
					int sub = q.front();
					q.pop();
					if (sub < groups[i].typ_edge)
						groups[i].typ_edge = sub;
					for (int j = 0; j < in_which_triangle[sub].size(); j++) {
						int temp = in_which_triangle[sub][j];
						if (!Triangles[temp].is_broken) {
							if (Triangles[temp].is_balanced)
							{
								link[Triangles[temp].edge1]++;
								link[Triangles[temp].edge2]++;
								link[Triangles[temp].edge3]++;
							}
							else
								break_unb[i]++;
							if (link[Triangles[temp].edge1] == 1)
								is_linked.push_back(Triangles[temp].edge1);
							if (link[Triangles[temp].edge2] == 1)
								is_linked.push_back(Triangles[temp].edge2);
							if (link[Triangles[temp].edge3] == 1)
								is_linked.push_back(Triangles[temp].edge3);

							if (support[Triangles[temp].edge1] - link[Triangles[temp].edge1] < k - 2 && !temp_delete_e[Triangles[temp].edge1])
							{
								temp_delete_e[Triangles[temp].edge1] = 1;
								q.push(Triangles[temp].edge1);
								followers[i].push_back(Triangles[temp].edge1);
							}


							if (support[Triangles[temp].edge2] - link[Triangles[temp].edge2] < k - 2 && !temp_delete_e[Triangles[temp].edge2])
							{
								temp_delete_e[Triangles[temp].edge2] = 1;
								q.push(Triangles[temp].edge2);
								followers[i].push_back(Triangles[temp].edge2);
							}


							if (support[Triangles[temp].edge3] - link[Triangles[temp].edge3] < k - 2 && !temp_delete_e[Triangles[temp].edge3])
							{
								temp_delete_e[Triangles[temp].edge3] = 1;
								q.push(Triangles[temp].edge3);
								followers[i].push_back(Triangles[temp].edge3);
							}

						}
					}
				}
				//if (break_unb[i] == 0)
					//groups[i].is_delete = 1;
				if (break_unb[i]>0)
				{
					if (followers[i].size() < followers[min].size())
						min = i;
					else if (followers[i].size() == followers[min].size()) {
						if (groups[i].typ_edge < groups[min].typ_edge)
							min = i;
					}
				}

				for (int j = 0; j < is_linked.size(); j++)
				{
					link[is_linked[j]] = 0;
					temp_delete_e[is_linked[j]] = 0;
				}
				is_linked.clear();
			}
			
			
		}

		if (min != MAX_E)
		{
			for (int i = 0; i < followers[min].size(); i++)
			{
				is_delete_e[followers[min][i]] = 1;
				temp_delete_e[followers[min][i]] = 1;
				//support[followers[min][i]] = 0;
			}
			size_of_truss -= followers[min].size();
			for (int i = 0; i < followers[min].size(); i++) {
				for (int j = 0; j < in_which_triangle[followers[min][i]].size(); j++) {
					int temp = in_which_triangle[followers[min][i]][j];
					if (!Triangles[temp].is_broken) {
						if (Triangles[temp].is_balanced) {
							if (!is_delete_e[Triangles[temp].edge1])
								--support[Triangles[temp].edge1];
							if (!is_delete_e[Triangles[temp].edge2])
								--support[Triangles[temp].edge2];
							if (!is_delete_e[Triangles[temp].edge3])
								--support[Triangles[temp].edge3];
						}
						else
							unbalance_num--;
						Triangles[temp].is_broken = 1;
					}
				}
			}
			groups[min].is_delete = 1;
			if (counts % 100 == 0) {
				cout << "sizeof truss" << size_of_truss << endl;
				cout << "unb num:" << unbalance_num << endl;
				cout << counts << " " << groups[min].typ_edge << endl;
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

            for (int neighbor : vec[current_vertex]) {
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
	GetKtruss(start_vertex);

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
	for (int i = 1; i < MAX_V; i++) {
		if (hop_num[i] < 0) break;
		cout << "hop number for " << i << " is " << hop_num[i] << "\n";
	}

	return 0;
}