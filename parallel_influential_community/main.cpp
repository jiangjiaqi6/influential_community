#include "Graph.h"

void print_usage() ;



void process_in_threadpool_plus(Graph *graph, int thread_num){
	ThreadPool pool(thread_num);
	std::mutex queue_mutex;
	// graph->extract_subgraph_from_original();
	ui *peel_sequence = new ui[graph->n]();
	ui *core = new ui[graph->n]();
	ui max_core = graph->core_decomposition_parallel(peel_sequence,core);
	// graph->greedy_and_reduce_plus(peel_sequence, core);
	graph->greedy_and_reduce_plus_plus(peel_sequence, core);
	graph->exact_pseudoarboricity_computation();
    graph->left = graph->tau;
    graph->right = graph->max_pa;
	
    {
		if(graph->edgeList != nullptr) delete[] graph->edgeList;
		graph->edgeList = new ui[2*graph->original_m];
		if(graph->pstart != nullptr) delete[] graph->pstart;
		graph->pstart = new ui[graph->original_n]();
		if(graph->pend != nullptr) delete[] graph->pend;
		graph->pend = new ui[graph->original_n]();
		if(graph->deg != nullptr) delete[] graph->deg;
		graph->deg = new ui[graph->original_n]();
    }

	


    std::condition_variable condition;
    std::atomic<int> active_tasks(0);


	graph->ori_map.clear();  // in greedy_and_reduce function
	ui mid = (graph->left + graph->right) / 2;
	// graph->greedy_and_reduce_plus(peel_sequence,core,mid);
	graph->greedy_and_reduce_plus_plus(peel_sequence,core,mid);
	delete[] peel_sequence;
	delete[] core;
    graph->p = new ui[graph->n]();
	graph->cur_vis = new bool[graph->n]();
    graph->test_value = mid-1;
	graph->find_maximum_weight_subgraph_2approximate_binary();
    delete []graph->p;
    delete []graph->cur_vis;

    std::function<void(Graph*)> process_task = [&](Graph* g) {
        ui left_ = g->left, right_ = g->right, mid = (left_ + right_) / 2;
        if (left_ > right_) {
            delete g;
        } else {
            g->compute_mid_pseudoarb(mid);
            Graph* g_left = new Graph();
            Graph* g_right = new Graph();
            g->partition_graph(left_, right_, mid, g_left, g_right);
			delete g;

			if(g_left->left <= g_left->right){
				pool.enqueue( g_left->n,process_task, g_left);
				active_tasks+=1;
			}
			else
				delete g_left;
			if(g_right->left <= g_right->right){
				pool.enqueue(g_right->n,process_task, g_right);
				active_tasks+=1;
			}
			else
				delete g_right;			            
        }
        active_tasks--;
        condition.notify_one();
    };

	Graph* g_left = new Graph();
	Graph* g_right = new Graph();
	graph->partition_graph(graph->left, graph->right, mid, g_left, g_right);
	pool.enqueue(g_left->n,process_task, g_left);
	active_tasks+=1;
	pool.enqueue(g_right->n,process_task, g_right);
	active_tasks+=1;

	delete graph;
	// graph = nullptr;

    std::unique_lock<std::mutex> lock(queue_mutex);
    condition.wait(lock, [&active_tasks] { return active_tasks.load() == 0; });
}

void filter_from_global_graph(Graph* global_g, ui pa_, Graph* local_g, ui* peel_sequence, ui* core, ui* min_wei_thread){
	
	ui idx = 0, n = global_g->n;
	std::unordered_map<ui,ui> gv_to_lv;  //global vertex info to local vertex info
	bool *is_in = new bool[n]();
	std::vector<ui> l_weigh;
	ui tid = omp_get_thread_num();
	for(int i = n-1; i >=0; i--){
		ui now = peel_sequence[i];
		if(global_g->weigh[now] <= min_wei_thread[tid]) continue;
		if(core[now] < pa_) break;
		gv_to_lv[now] = idx++;
		is_in[now] = true;
		l_weigh.push_back(global_g->weigh[now]);
	}
	ui local_v = idx, local_m = 0;
	ui *degree = new ui[local_v]();
	for(ui i = 0; i < n; i++) {
		if(is_in[i]){
			degree[gv_to_lv[i]] = 0;
			for(ui j = global_g->pstart[i]; j < global_g->pend[i]; j++) if(is_in[global_g->edgeList[j]]) ++ degree[gv_to_lv[i]];
			local_m += degree[gv_to_lv[i]];	
		}	
	}
	local_m /= 2;
	local_g->malloc_glocal_info(local_v,local_m);
	memcpy(local_g->weigh,l_weigh.data(),l_weigh.size()*sizeof(ui));
	local_g->pstart[0] = 0;
	for(ui i = 1; i < local_v; i++)
	{
		local_g->pstart[i] = local_g->pstart[i-1] + degree[i-1];
		local_g->pend[i-1] = local_g->pstart[i];
	}
	local_g->pend[local_v-1] = local_g->pstart[local_v-1] + degree[local_v-1];

	idx = 0;
	#pragma omp parallel for
	for(ui i = 0; i < n; i++) {
		if(is_in[i]){
			ui u = gv_to_lv[i];
			local_g->pend[u] = local_g->pstart[u];
			for(ui j = global_g->pstart[i]; j < global_g->pend[i]; j++){
				ui g_v = global_g->edgeList[j];
				if(is_in[g_v]){
					ui v = gv_to_lv[g_v];
					local_g->edgeList[local_g->pend[u]++] = v;
				}
			} 
		}	
	}
	#pragma omp parallel for
	for(ui i = 0; i < local_v; i++){	
		std::sort(local_g->edgeList+local_g->pstart[i],local_g->edgeList+local_g->pstart[i]+degree[i]);
	}
	for(ui i = 0; i < local_v; i++) {
		ui u = i;
		local_g->ori_map[i] = i;
		for(ui j = local_g->pstart[i]; j < local_g->pend[i]; j++){
			ui v = local_g->edgeList[j];
			if(v > u){
				local_g->eid[j] = idx;
				local_g->eid[local_g->pstart[v]++] = idx++;
			}
		} 	
	}
	
	#pragma omp parallel for
	for(ui i = 0; i < local_v; i++)
	{
		local_g->pstart[i] = local_g->pend[i] - degree[i];
	}
	memset(is_in,false,sizeof(bool)*n);
	for(int i = local_v-1; i >= 0; i--){
		ui u = i;
		for(ui j = local_g->pstart[u]; j < local_g->pend[u]; j++){
			ui v = local_g->edgeList[j];
			if(is_in[v])continue;
			local_g->deg[u]++;
			local_g->orien_edges[local_g->eid[j]] = (edge){.src = u, .end = v, .to = v};
		}
		is_in[u] = true;
	}
	delete[] degree;
	delete[] is_in;
	for(ui i = 0; i < local_v; i++){
		ui u = i;
		for(ui j = local_g->pstart[u]; j < local_g->pend[u]; j++){
			edge& ne = local_g->orien_edges[local_g->eid[j]];
			if(ne.to != u) continue;
			ui tar = ne.src == u ? ne.end : ne.src;
			if(local_g->deg[u] <= local_g->deg[tar] - 2){
				ne.to = tar;
				local_g->deg[u]++;
				local_g->deg[tar]--;
			}
		}
	}
}

void compute_task(Graph* graph, ui pa_, ui* peel_sequence, ui* core, ui* min_wei_thread){
	Graph *local_g = new Graph();
	ui tid = omp_get_thread_num();
	filter_from_global_graph(graph,pa_,local_g,peel_sequence,core,min_wei_thread);
	local_g->compute_mid_pseudoarb(pa_,false);
	min_wei_thread[tid] = local_g->min_wei;
	delete local_g;
}

void process_in_parallel_block(Graph *graph, int thread_num){
	graph->extract_subgraph_from_original();
	graph->ori_map.clear();  // in greedy_and_reduce function
	ui *peel_sequence = new ui[graph->n]();
	ui *core = new ui[graph->n]();
	ui *peel_seq = new ui[graph->n]();
	ui *core_ = new ui[graph->n]();

	// ui max_core = graph->core_decomposition(peel_sequence, core, graph->n, graph->pstart, graph->pend, graph->edgeList, graph->eid);  
	ui max_core = graph->core_decomposition_parallel(peel_sequence,core);
	memcpy(peel_seq,peel_sequence,sizeof(ui)*graph->n);
	memcpy(core_,core,sizeof(ui)*graph->n);
	graph->greedy_and_reduce_plus(peel_sequence, core);
	delete[] peel_sequence;
	delete[] core;
	graph->exact_pseudoarboricity_computation();
	graph->n = graph->original_n, graph->m = graph->original_m;
	delete[] graph->pstart;
	delete[] graph->pend;
	delete[] graph->edgeList;
	delete[] graph->orien_edges;
	graph->orien_edges = nullptr;

	graph->pstart = graph->original_pstart;
	graph->original_pstart = nullptr;
	graph->pend = graph->original_pend;
	graph->original_pend = nullptr;
	graph->edgeList = graph->original_edgeList;
	graph->original_edgeList = nullptr;
	ui max_pa = graph->max_pa, tau = graph->tau;
	
	ui *min_wei_thread = new ui[thread_num]();

	#pragma omp parallel for
	for(int i = max_pa-1; i >= tau; i--){
		compute_task(graph,i,peel_seq,core_,min_wei_thread);
	}

	delete graph;
	delete[] peel_seq;
	delete[] core_;
}

int main(int argc, char *argv[]) {
#ifndef NDEBUG
	printf("!!! You may want to define NDEBUG in utilities/Defines.h to get better performance!\n");
#endif
	// if(argc < 2) {
	// 	print_usage();
	// 	return 0;
	// }
#ifndef NDEBUG
	printf("**** DensestSubgraph (Debug) build at %s %s ***\n", __TIME__, __DATE__);
#else
	printf("**** DensestSubgraph (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif
	Clock allClock("All");
    log_info(allClock.Start());
	FILE* dataset = fopen("./dataset.txt", "r");
	char graph_name[100];
	int thres = 0, thread_num = 0;
	while (fscanf(dataset, "%s %d %d", &graph_name, &thres, &thread_num) == 3)
	{
		Timer timer;
		Graph *graph = new Graph(graph_name);
		graph->read_graph_binary(graph_name);
		// graph->read_graph_without_vm(graph_name);
		log_info(allClock.Count("graph : %s ,thread_num: %u",graph_name,thread_num));
		// graph->write_graph_binary(graph_name);
		omp_set_num_threads(thread_num);
		graph->set_tau(thres);
		graph->set_thrnmb(thread_num);
		graph->load_vertex_weight(graph_name);
		#ifdef TP
		process_in_threadpool_plus(graph,thread_num);
		#else
		process_in_parallel_block(graph,thread_num);
		#endif

		size_t peakRSS = Utility::getPeakRSS();
		log_info(allClock.Count("Finish\n"));
			
	}
	fclose(dataset);

	return 0;
}

void print_usage() {
	printf("Usage: [1]exe [2]graph-dir [3 optional] \"output\"\n");
}