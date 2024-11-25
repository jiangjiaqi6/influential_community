#ifndef GRAPH_H_
#define GRAPH_H_

#include "utilities/Defines.h"
#include "utilities/Utility.h"
#include "utilities/Timer.h"
#include "utilities/readFile.h"
#include "utilities/tools.h"
#include "utilities/Time.h"
#include "utilities/log.h"

using lint = long long;
#define min(a,b) ((a)<(b)?(a):(b))

const lint LINT_INF = (lint)2000000000*(lint)2000000000;
typedef struct edge{
	ui src;
	ui end;
	ui to;
}edge;

class Graph {
private:
	std::string dir; //input graph directory
	ui n,original_n; //#nodes of the graph
	ui m,original_m, update_m; //#edges of the graph
	ui maxDeg; //max degree
	ui number;
	ui max_id;
	ui pa; //pseudoarboricity
	ui max_pa; //maximum pseudoarboricity
	ui app_pa; //approximate pseudoarboricity
	ui tau; //the threshold of pseudoarboricity
	ui min_wei; //weight of densest subraph

	ui *original_pstart;
	ui *original_pend;
	ui *original_edgeList;

	ui *pstart; //start positions of neighbors of vertices in the array "edges"
	ui *pend; //end positions of neighbors of vertices in the array "edges"
	ui *edgeList; //concatenation of neighbors of all vertices
	ui *eid; // id of each edge
	ui *deg; // the out degree of each vertex
	ui *weigh; //the influence of each vertex
	ui *sub_weigh; //the influence of each vertex in the induced subgraph
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

	//output file
	// std::vector<std::pair<ui,double>> out_res;  
	std::vector<std::pair<ui,ui>> out_res;

public:
	Graph(const char *_dir) ;
	~Graph() ;
	bool get_finished(){ return isFinish;}
	void set_tau(ui _tau){tau = _tau;}
	void read_graph_binary(const char * graph_name) ;
	void read_graph(const char *input_file) ;
	void read_graph_without_vm(const char *input_file) ;
	void load_vertex_weight(const char *input_file);
	void generate_non_unique_weights(std::string ss);
	void generate_unique_weights(std::string ss);
	void generate_pagerank_unique_weights(std::string ss);
	void extract_subgraph_from_original();

	ui get_node_id(ui id); //used in reordering vertex id when reading original graph
	
	// integer densest subgraph model
	void influential_community();
	void influential_community_divide_and_conquer();
	void div_con_framework();
	void recursive_function(ui left, ui right, ui ver_n, ui edg_m, ui *pstart, ui *pend, ui *edges, ui *cur_eid);
	void non_recursive_function();
	ui get_pseudoarboricity(){ return pa;}
	
	

	// kmax-core model
	void influential_community_based_core();
	void greedy_and_reduce_kmax_core(ui *peel_sequence, ui *core);
	void maxweight_core_computation(ui max_core);
	void dfs_deletion(ui u, ui max_core, ui del_v, ui del_m);


	void write_output_weight_density(char * name, ui flag);
	void write_output_prune_effect(char * name, ui flag);

private:
	ui core_decomposition(ui *peel_sequence, ui *core, const ui n, const ui *pstart, const ui *pend, const ui *edges, const ui *eid);
	void greedy_and_reduce(ui *peel_sequence, ui *core, ui mid = 0);
	void greedy_and_reduce_plus(ui *peel_sequence, ui *core, ui input_mid = 0);

	void exact_pseudoarboricity_computation();
	void find_maximum_weight_subgrah_peeling();
	void find_maximum_weight_subgrah_2approximate();
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