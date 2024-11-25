#ifndef GRAPH_H_
#define GRAPH_H_

#include "utilities/Defines.h"
#include "utilities/Utility.h"
#include "utilities/Timer.h"
#include "utilities/readFile.h"
#include "utilities/tools.h"
#include "utilities/Time.h"
#include "utilities/log.h"
#include "utilities/ThreadPool.h"

using lint = long long;
#define min(a,b) ((a)<(b)?(a):(b))

const lint LINT_INF = (lint)2000000000*(lint)2000000000;
typedef struct edge{
	ui src;
	ui end;
	ui to;
}__attribute__((aligned(4))) edge;

class Graph {
public:
	std::string dir; //input graph directory
	ui n,original_n; //#nodes of the graph
	ui m,original_m; //#edges of the graph
	ui maxDeg; //max degree
	ui number;
	ui max_id;  //used for reordering vertex id
	ui pa; //pseudoarboricity
	ui max_pa; //maximum pseudoarboricity
	ui app_pa; //approximate pseudoarboricity
	ui tau; //the threshold of pseudoarboricity
	ui min_wei; //weight of densest subraph
	ui left, right; 
	ui thread_number;

	ui *original_pstart;
	ui *original_pend;
	ui *original_edgeList;

	ui *pstart; //start positions of neighbors of vertices in the array "edges"
	ui *pend; //end positions of neighbors of vertices in the array "edges"
	ui *edgeList; //concatenation of neighbors of all vertices
	ui *eid; // id of each edge
	ui *deg; // the out degree of each vertex
	ui *weigh; //the influence of each vertex
	bool *vis; // flag of original graph
	bool *cur_vis; // flag of integer denset subgraph
	edge *orien_edges;
	bool isFinish; // if density of densest less than tau, finish execution


	//pseudoarboricity
	ui *p;
	ui test_value;

	//flow 
	ui *dist;
	ui *cur;


	std::vector<ui> res;
	std::unordered_map<ui,ui> map_id;
	std::unordered_map<ui,ui> ori_map;
	std::unordered_map<ui,ui> final_map;
	std::unordered_set<ui> isContain; // used in recording integer_densest_sub

	std::ofstream fout;
	Clock graphClock_;

	std::mutex visited_mutex;
    std::unordered_map<int, bool> visited;

public:
	Graph();
	Graph(const char *_dir) ;
	~Graph() ;
	bool get_finished(){
		return isFinish;
	}
	void partition_by_edge(ui k,const char *graph_name);
	void write_graph_binary(const char *graph_name);
	void read_graph_binary(const char * graph_name) ;
	void read_graph(const char *input_file) ;
	void read_graph_without_vm(const char *input_file) ;
	ui get_node_id(ui id);

	ui get_pseudoarboricity(){
		return pa;
	}
	void set_tau(ui _tau){
		tau = _tau;
		// printf("tau: %u\n",tau);
	}
	void set_thrnmb(ui _tau){
		thread_number = _tau;
	}
	ui get_tau(){
		return tau;
	}
	ui get_max_pseudoarb(){
		return max_pa;
	}
	void influential_community();
	void recursive_function(ui left, ui right, ui ver_n, ui edg_m, ui *pstart, ui *pend, ui *edges, ui *cur_eid);
	void compute_mid_pseudoarb(ui mid, bool is_reduced = false);
	void load_vertex_weight(const char *input_file);
	void extract_subgraph_from_original();
	void generate_non_unique_weights(std::string ss);
	void generate_unique_weights(std::string ss);
	void max_pseudoarb_community();
	void partition_graph(ui left_, ui right_, ui mid, Graph* g_left, Graph* g_right);
	void transfer_graphinfo_parallel(Graph *g);
	void transfer_graphinfo_parallel_plus(Graph *g);
	void transfer_graphinfo_parallel_right(Graph *g);
	void transfer_graphinfo(Graph *g);

	void malloc_glocal_info(ui n, ui m);
	int getRunningThreadCount();


public:
	ui core_decomposition(ui *peel_sequence, ui *core, const ui n, const ui *pstart, const ui *pend, const ui *edges, const ui *eid);
	ui core_decomposition_parallel(ui *peel_sequence, ui *core);
	void scan(ui *core, ui level, ui* curr, int *curr_tail);
	void process_sublevel(ui *curr, int curr_tail, ui *core, ui level, ui *next, int *next_tail);
	void process_sublevel_plus(ui *curr, int curr_tail, ui *core, ui level, ui *next, int *next_tail);

	void greedy_and_reduce(ui *peel_sequence, ui *core, ui mid = 0);
	void greedy_and_reduce_plus(ui *peel_sequence, ui *core, ui mid = 0);
	void greedy_and_reduce_plus_plus(ui *peel_sequence, ui *core, ui mid = 0);
	void exact_pseudoarboricity_computation();
	void find_maximum_weight_subgraph_peeling();
	void find_maximum_weight_subgraph_2approximate();
	void find_maximum_weight_subgraph_2approximate_binary();
	ui get_max_d();
	

	bool ReTest();
	bool DinicBFS(std::unordered_set<ui>& nodes_to_decrease);
	bool DinicDFS(ui x);
	void find_integer_densest_sub(std::unordered_set<ui>& isContainVer);
	bool isPowerOfTwo(ui x){
		return (x > 0) && ((x & (x - 1)) == 0);
	}
	void construct_subgraph(std::vector<ui> &vertexSet, ui pre_count);

};

#endif /* GRAPH_H_ */