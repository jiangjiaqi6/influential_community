#include "Graph.h"

using namespace std;
bool TESTING = true;
namespace TEST
{
	unsigned p0, iteration_time, phase_number;
}

Graph::Graph(const char *_dir):graphClock_("graph process") 
{
	log_info(graphClock_.Start());
	// dir = string(_dir);
	// ui index_s = 0, index_e = 0;
	// for(ui i = 0; i < dir.size(); i++){
	// 	if(dir[i] == '\\')
	// 		index_s = i+1;
	// 	else if(dir[i] == '.')
	// 		index_e = i;
	// }
	// string ss = dir.substr(index_s);

	// string ss(_dir);
	// log_info(graphClock_.Count("graph file name: %s", _dir));
	// ss = "result/out-" + ss;
	// fout.open(ss);

    // dir = "/home/jjq/research/ktruss/ktruss-data/";
	dir = "/mnt/hdd0/jjq/graphData/";

	n = 0, original_n = 0;
	m = 0, original_m = 0;
	max_id = 0;
	maxDeg = 0;
	number = 1;
	pa = 0;
	max_pa = 0;
	app_pa = 0;
    left = right = 0;
	thread_number = 0;
    min_wei = INT32_MAX;

	// global info
	original_pstart = nullptr;
	original_pend = nullptr;
	original_edgeList = nullptr;
	vis = nullptr;
	weigh = nullptr;
	

	// subgraph info
	pstart = nullptr;
	pend = nullptr;
	edgeList = nullptr;
	orien_edges = nullptr;
	eid = nullptr;
	deg = nullptr;

	res.clear();
	isFinish = false;
}

Graph::Graph():graphClock_("graph process") 
{
	graphClock_.Start();
	dir = "/mnt/hdd0/jjq/graphData/";

	n = 0, original_n = 0;
	m = 0, original_m = 0;
	max_id = 0;
	maxDeg = 0;
	number = 1;
	pa = 0;
	max_pa = 0;
	app_pa = 0;
    left = right = 0;
	thread_number = 0;
    min_wei = INT32_MAX;

	// global info
	original_pstart = nullptr;
	original_pend = nullptr;
	original_edgeList = nullptr;
	vis = nullptr;
	weigh = nullptr;
	

	// subgraph info
	pstart = nullptr;
	pend = nullptr;
	edgeList = nullptr;
	orien_edges = nullptr;
	eid = nullptr;
	deg = nullptr;

	res.clear();
	isFinish = false;
}

Graph::~Graph() {
	// fout.close();
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(pend != nullptr) {
		delete[] pend;
		pend = nullptr;
	}
	if(edgeList != nullptr) {
		delete[] edgeList;
		edgeList = nullptr;
	}
	if(original_pstart != nullptr) {
		delete[] original_pstart;
		original_pstart = nullptr;
	}
	if(original_pend != nullptr) {
		delete[] original_pend;
		original_pend = nullptr;
	}
	if(original_edgeList != nullptr) {
		delete[] original_edgeList;
		original_edgeList = nullptr;
	}
    if(eid != nullptr) {
		delete[] eid;
		eid = nullptr;
	}
	if(vis != nullptr) {
		delete[] vis;
		vis = nullptr;
	}
	if(weigh != nullptr) {
		delete[] weigh;
		weigh = nullptr;
	}
	if(deg != nullptr) {
		delete[] deg;
		deg = nullptr;
	}
    if(orien_edges != nullptr) {
		delete[] orien_edges;
		orien_edges = nullptr;
	}
}

void Graph::partition_by_edge(ui percentage, const char *graph_name){
	string ss = dir;
	int i = 0, j = 0, len = strlen(graph_name);
	bool flag = true;
	while(i < len){
		ss += graph_name[i];
		if(graph_name[i] == '/'){
			flag = false;
			break;
		}		
		i++;
	}
	if(flag)
		ss = dir;
	cout<<ss<<endl;
	ss += "scalability/";
	if (0 != access(ss.c_str(), 0)){
        int isCreate = mkdir(ss.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
    }
    ss = ss + to_string(percentage) +"/";
	if (0 != access(ss.c_str(), 0)){
        int isCreate = mkdir(ss.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
    }

	bool *is_in = new bool[original_m]();
	bool *is_vis = new bool[original_n]();
	double stride = 10.0 / percentage;
	
	

	ui sub_edge = (original_m/10) * percentage;
	for(ui i = 0; i < sub_edge; i++){
		ui id = round(i*stride);
		is_in[id] = true;
	}
	ui c = 0;
	unordered_map<ui,ui> verMap;
	vector<pair<ui,ui>> allEdge;
	unordered_set<ui> verSet;

	for(ui i = 0; i < n; i++){
		for(ui j = original_pstart[i]; j < original_pend[i]; j++){
			if(original_edgeList[j] < i) continue;
			if(is_in[c]){
				allEdge.push_back({i,original_edgeList[j]});
				is_vis[i] = true;
				is_vis[original_edgeList[j]] = true;
			}
			c++;
		}
	}
	c = 0;
	for(ui i = 0; i < original_n; i++){
		if(is_vis[i]){
			verMap[i] = c++;
		}
	}

	assert(sub_edge == allEdge.size());
	ui sub_ver = c;
	log_info(graphClock_.Count("sub_v: %u, sub_e: %u\n",sub_ver,sub_edge));
	for(ui i = 0; i < allEdge.size(); i++){
		pair<ui,ui> tmp = allEdge[i];
		allEdge[i] = {verMap[tmp.first],verMap[tmp.second]};
	}
	ui *degree = new ui[sub_ver]();
	for(ui i = 0; i < allEdge.size(); i++){
		pair<ui,ui> tmp = allEdge[i];
		degree[tmp.first]++;
		degree[tmp.second]++;
	}
	pstart = new ui[sub_ver]();
	pend = new ui[sub_ver]();
	edgeList = new ui[2*sub_edge]();
	for(ui i = 1; i < sub_ver; i++){
		pstart[i] = pend[i-1] = pstart[i-1] + degree[i-1];
	}
	pend[sub_ver-1] = pstart[sub_ver-1] + degree[sub_ver-1];
	
	for(ui i = 0; i < allEdge.size(); i++){
		pair<ui,ui> tmp = allEdge[i];
		edgeList[pstart[tmp.first]++] = tmp.second;
		edgeList[pstart[tmp.second]++] = tmp.first;
	}

	for(ui i= 0; i < sub_ver; i++)
		pstart[i]  = pend[i] - degree[i];

	// for(ui i = 0; i < sub_ver; i++){
	// 	printf("\nu: %u\n",i);
	// 	for(ui j = pstart[i]; j < pend[i]; j++){
	// 		printf("%u ",edgeList[j]);
	// 	}
	// }

	FILE *f = Utility::open_file((ss + string("b_degree.bin")).c_str(), "w");
	ui num = sizeof(ui);
	fwrite(&num,sizeof(ui),1,f);
	fwrite(&sub_ver, sizeof(ui), 1, f);
	fwrite(&sub_edge, sizeof(ui), 1, f);
	
	
	fwrite(degree, sizeof(ui), sub_ver, f);

	fclose(f);
	f = Utility::open_file((ss + string("b_adj.bin")).c_str(), "w");
	fwrite(edgeList, sizeof(ui), 2*sub_edge, f);
	fclose(f);

	delete[] degree;
	delete[] is_in;
}

void Graph::write_graph_binary(const char *graph_name){
	printf("# Start writing graph, files \"b_degree.bin\" and \"b_adj.bin\"\n");
    string ss = dir;
	int i = 0, j = 0, len = strlen(graph_name);
	bool flag = true;
	while(i < len){
		ss += graph_name[i];
		if(graph_name[i] == '/'){
			flag = false;
			break;
		}		
		i++;
	}
	if(flag)
		ss = dir;
	cout<<ss<<endl;

	
	FILE *f = Utility::open_file((ss + string("b_degree.bin")).c_str(), "w");
	ui num = sizeof(ui);
	fwrite(&num,sizeof(ui),1,f);
	fwrite(&original_n, sizeof(ui), 1, f);
	fwrite(&original_m, sizeof(ui), 1, f);
	ui *degree = new ui[n]();
	for(ui i= 0; i < n; i++)
		degree[i] = original_pend[i]-original_pstart[i];
	fwrite(degree, sizeof(ui), n, f);

	fclose(f);
	f = Utility::open_file((ss + string("b_adj.bin")).c_str(), "w");
	fwrite(original_edgeList, sizeof(ui), 2*original_m, f);
	fclose(f);

	delete[] degree;
}


void Graph::read_graph_binary(const char * graph_name) {
	printf("# Start reading graph, Require files \"b_degree.bin\" and \"b_adj.bin\"\n");
	string ss = dir;
	int i = 0, j = 0, len = strlen(graph_name);
	bool flag = true;
	while(i < len){
		ss += graph_name[i];
		if(graph_name[i] == '/'){
			flag = false;
			break;
		}		
		i++;
	}
	if(flag)
		ss = dir;
	FILE *f = Utility::open_file((ss + string("b_degree.bin")).c_str(), "rb");

	ui tt;
	fread(&tt, sizeof(ui), 1, f);
	if(tt != sizeof(ui)) {
		printf("sizeof unsigned int is different: b_degree.bin(%u), machine(%lu)\n", tt, sizeof(ui));
		return ;
	}

	fread(&n, sizeof(ui), 1, f);
	fread(&m, sizeof(ui), 1, f);
	original_n = n, original_m = m;

	log_info(graphClock_.Count("*\tn = %s; m = %s (undirected), graph: %s", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m).c_str(),graph_name));

	ui *degree = new ui[n];
	fread(degree, sizeof(ui), n, f);

	fclose(f);

#ifndef NDEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	assert(sum == m*2);
#endif

	f = Utility::open_file((ss + string("b_adj.bin")).c_str(), "rb");

	if(original_pstart != nullptr) delete[] original_pstart;
	original_pstart = new ui[n]();
	if(original_pend != nullptr) delete[] original_pend;
	original_pend = new ui[n];
	if(original_edgeList != nullptr) delete[] original_edgeList;
	original_edgeList = new ui[2*m];
	if(deg != nullptr) delete[] deg;
	deg = new ui[original_n]();
	if(vis != nullptr) delete[] vis;
	vis = new bool[original_n]();

	original_pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) fread(original_edgeList+original_pstart[i], sizeof(ui), degree[i], f);
		if(i < n-1)
			original_pstart[i+1] = original_pend[i] = original_pstart[i] + degree[i];
		else
			original_pend[i] = original_pstart[i] + degree[i];
	}
	// for(ui i = 0; i < n; i++)
	// 	for(ui j = original_pstart[i]; j < original_pend[i]; j++)
	// 		printf("u: %u, v: %u\n",i,original_edgeList[j]);

	fclose(f);
	delete[] degree;
}



ui Graph::get_node_id(ui id){
	if( map_id.find(id) == map_id.end() ) {
		map_id[id] = max_id;
		// org_id.push_back( id );   //org_id[max_id] = id; renturn original vertex id;
		return max_id++;
	}
	return map_id[id];
}


void Graph::read_graph(const char *input_file) {
	printf("# Start reading graph %s\n", input_file);
	FILE *f = Utility::open_file(input_file, "r");

	fscanf(f, "%u%u", &n, &m);
	++ n; 
	m *= 2;
	printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	vector<pair<ui,ui> > vp;
	for(ui i = 0;i < m/2;i ++) {
		ui a, b;
		fscanf(f, "%u%u", &a, &b);
		vp.pb(mp(a,b));
		vp.pb(mp(b,a));
	}
	sort(vp.begin(), vp.end());

	if(pstart != nullptr) delete[] pstart;
	pstart = new ui[n+3];
	if(pend != nullptr) delete[] pend;
	pend = new ui[n+3];
	if(edgeList != nullptr) delete[] edgeList;
	edgeList = new ui[m+4*n];

	pstart[0] = 0;
	ui idx = 0;
	for(ui i = 0;i < n;i ++) {
		pend[i] = pstart[i];
		while(idx < vp.size()&&vp[idx].first == i) edgeList[pend[i] ++] = vp[idx ++].second;
		pstart[i+1] = pend[i];
	}

	fclose(f);

}

//Renumber all vertices in the graph with get_node_id function, starting from 0.
void Graph::read_graph_without_vm(const char *graph_name) {
	FILE *file;
	// char input_file[100] = "/home/jjq/research/ktruss/ktruss-data/";
	// strcat(input_file, graph_name);
    string input_file = dir;
    input_file.append(graph_name);
	file = fopen(input_file.c_str(),"r");
	char line[200];
	vector<pair<ui,ui> > vp;
	while(fgets(line,200,file)){
		if(line[0] == '%' || line[0] == '#')continue;
		ui i = 0;
        ui a = 0;
		while (line[i] < '0' || line[i]>'9') i++;
		while (line[i] >= '0' && line[i] <= '9') a = a * 10 + line[i] - '0', i++;
		ui u = get_node_id(a);
		// u = a;  //vertex id from original graph file
        ui b = 0;
		while (line[i] < '0' || line[i]>'9') i++;
		while (line[i] >= '0' && line[i] <= '9') b = b * 10 + line[i] - '0', i++;
		ui v = get_node_id(b);
		// v = b; //vertex id from original graph file
        if(u == v)continue;
        vp.pb(mp(u,v));
        vp.pb(mp(v,u));
    }
	sort(vp.begin(),vp.end());
    vp.erase(unique(vp.begin(), vp.end()), vp.end()); //Edge deduplication
	n = max_id;
    m = vp.size()/2;

	original_n = n, original_m = m;

	if(original_pstart != nullptr) delete[] original_pstart;
	original_pstart = new ui[n];
	if(original_pend != nullptr) delete[] original_pend;
	original_pend = new ui[n];
	if(original_edgeList != nullptr) delete[] original_edgeList;
	original_edgeList = new ui[m*2];
    if(deg != nullptr) delete[] deg;
	deg = new ui[original_n]();
	if(vis != nullptr) delete[] vis;
	vis = new bool[n]();

	original_pstart[0] = 0;
	ui idx = 0;
	for(ui i = 0;i < n;i ++) {
		original_pend[i] = original_pstart[i];
		while(idx < vp.size()&&vp[idx].first == i) original_edgeList[original_pend[i] ++] = vp[idx ++].second;
		if(i < n-1)
			original_pstart[i+1] = original_pend[i];
		maxDeg = (original_pend[i]-original_pstart[i]) > maxDeg ? (original_pend[i]-original_pstart[i]) : maxDeg;
        deg[i] = original_pend[i]-original_pstart[i];
	}
	// for(ui i = 0; i < n; i++)
	// 	for(ui j = original_pstart[i]; j < original_pend[i]; j++)
	// 		printf("u: %u, v: %u\n",i,original_edgeList[j]);
	log_info(graphClock_.Count("*\tn = %s; m = %s ; maxDeg = %s (undirected)", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m).c_str(),Utility::integer_to_string(maxDeg).c_str()));
}

void Graph::generate_non_unique_weights(string ss){
//weight: [1,maxDeg]

	srand(time(0));
	/* Seed */ 
	random_device rd; 
	/* Random number generator */ 
	default_random_engine generator(rd()); 
	/* Distribution on which to apply the generator */ 
	std::uniform_int_distribution<long long unsigned> distribution(0,0xFFFFFFFFFFFFFFFF);
	ofstream fwrite(ss);
	if (fwrite.is_open()) {
		for(int i = 0; i < n; i++){
			ui wei_u = distribution(generator) % maxDeg + 1;
			fwrite << wei_u;
			fwrite << endl;  
			weigh[i] = wei_u;	
		}
		fwrite.close();
	} else {
		log_debug(graphClock_.Count("Failed to open file for writing"));
	}
}

// Function to generate unique weights for each vertex in graph
void Graph::generate_unique_weights(string ss){
    if (n > m) {
        std::cerr << "Number of vertices cannot be greater than the number of edges." << std::endl;
        return;
    }

    // Initialize a vector with sequential numbers from 1 to num_edges
    std::vector<ui> all_weights(m);
    for (int i = 0; i < m; ++i) {
        all_weights[i] = i + 1;
    }
    // Use a random device to seed the random number generator
    std::random_device rd;
    std::mt19937 g(rd());
    // Shuffle the weights to ensure they are randomly distributed
    std::shuffle(all_weights.begin(), all_weights.end(), g);
    // Assign the first num_vertices weights to the vertices
	ofstream fwrite(ss);
	if (fwrite.is_open()) {
		for(int i = 0; i < n; i++){
			fwrite << all_weights[i];
			fwrite << endl;  
			weigh[i] = all_weights[i];
		}
		fwrite.close();
	} else {
		log_debug(graphClock_.Count("Failed to open file for writing."));
	}
}

void Graph::load_vertex_weight(const char *graph_name){
	// string ss = "/home/jjq/research/ktruss/ktruss-data/";
    string ss = dir;
	int i = 0, len = strlen(graph_name);
	while(i < len && graph_name[i] != '.'){
		ss += graph_name[i++];
	}
	ss += "-weight.txt";
	// std::cout<<ss<<endl;
	ifstream input(ss);
	weigh = new ui[n]();
	
	if(input.good()){  //read prepared file 
		int i = 0;
    	assert(input.is_open());
		string line;
		while (!input.eof()) {
			getline(input, line);
			ui x = 0;
			auto it = line.begin();
			if (*it == '%' || *it == '#' || it == line.end()) {
				continue;
			}
			while (it != line.end() && (*it < '0' || *it > '9')) it++;
			while (it != line.end() && (*it >= '0' && *it <= '9')) x = x * 10 + (*it++) - '0';
			weigh[i++] = x;
		}
		// readFile rf(ss.c_str());
		// while(!rf.empty()){
		// 	ui u = rf.getUInt();
		// 	weigh[i++] = u;
		// }
		assert(i == n);
	}
	else{   //generate random number into file
	    generate_unique_weights(ss);
	}
}


void Graph::influential_community(){
	ori_map.clear();  // in greedy_and_reduce function
	ui *peel_sequence = new ui[n]();
	ui *core = new ui[n]();
	ui max_core = core_decomposition(peel_sequence, core, n, pstart, pend, edgeList, eid);
	greedy_and_reduce(peel_sequence, core);
	exact_pseudoarboricity_computation();
	delete[] peel_sequence;
	delete[] core;
}


void Graph::max_pseudoarb_community(){
	extract_subgraph_from_original();
	ori_map.clear();  // in greedy_and_reduce function
	ui *peel_sequence = new ui[n]();
	ui *core = new ui[n]();
	// ui max_core = core_decomposition(peel_sequence, core, n, pstart, pend, edgeList, eid);
	ui max_core = core_decomposition_parallel(peel_sequence,core);
	// greedy_and_reduce(peel_sequence, core);
	greedy_and_reduce_plus(peel_sequence, core);
    delete[] peel_sequence;
	delete[] core;
	exact_pseudoarboricity_computation();
    left = tau;
    right = max_pa;
}


void Graph::partition_graph(ui left_, ui right_, ui mid, Graph* g_left, Graph* g_right){
    // assert(min_wei != INT32_MAX);
    g_left->left = INT32_MAX, g_left->right = 0;
    if(left_ <= mid-1){
        g_left->transfer_graphinfo_parallel_plus(this);
        g_left->left = left_, g_left->right = mid-1;
    }
    g_right->left = INT32_MAX, g_right->right = 0;
    if(mid+1 <= right_){
        for(int i = 0; i < original_n; i++) vis[i] = true;

        // for(auto it = isContain.begin(); it != isContain.end(); it++)
        //     vis[ori_map[*it]] = false;
        // g_right->transfer_graphinfo_parallel_plus(this);

		for(auto it = isContain.begin(); it != isContain.end(); it++)
            vis[*it] = false;
		g_right->transfer_graphinfo_parallel_right(this);

        g_right->left = mid+1, g_right->right = right_;
    }
}

void Graph::compute_mid_pseudoarb(ui mid, bool is_reduced){
	if(!is_reduced){
		ui *peel_sequence = new ui[original_n]();
		ui *core = new ui[original_n]();
		ori_map.clear();  // in greedy_and_reduce function
		// ui max_core = core_decomposition(peel_sequence, core, n, pstart, pend, edgeList, eid);
		ui max_core = core_decomposition_parallel(peel_sequence, core);
		// greedy_and_reduce(peel_sequence,core,mid);

		#ifdef TP
		greedy_and_reduce_plus_plus(peel_sequence,core,mid);
		#else
		greedy_and_reduce_plus(peel_sequence,core,mid);
		#endif

		delete[] peel_sequence;
		delete[] core;
	}
    p = new ui[n]();
	cur_vis = new bool[n]();
    test_value = mid-1;
	find_maximum_weight_subgraph_2approximate_binary();
    // find_maximum_weight_subgraph_2approximate();
    delete []p;
    delete []cur_vis;
}


// core decomposition:
//		"peel_sequence" stores the vertices in degenerarcy order
//		"core[u]" stores the core value of vertex u
//		"n, pstart, pend, edges" are used to represent the graph
// return the maximum core value
ui Graph::core_decomposition(ui *peel_sequence, ui *core, const ui n, const ui *pstart, const ui *pend, const ui *edgeList, const ui *eid) {
#ifndef NDEBUG
	printf("Start core decomposition\n");
#endif
	Timer timer;
	ui *degree = new ui[n];
	for(ui i = 0;i < n;i ++) degree[i] = pend[i]-pstart[i];

	ui *rid = new ui[n];
	ui *id = peel_sequence;
	memset(id, 0, sizeof(ui)*n);
	for(ui i = 0;i < n;i ++) ++ id[degree[i]];
	for(ui i = 1;i < n;i ++) id[i] += id[i-1];

	for(ui i = 0;i < n;i ++) rid[i] = -- id[degree[i]];
	for(ui i = 0;i < n;i ++) id[rid[i]] = i;

	ui *degree_start = new ui[n+1];
	for(ui i = 0, j = 0;i <= n;i ++) {
		while(j < n&&degree[id[j]] < i) ++ j;
		degree_start[i] = j;
	}

	ui max_core = 0;
	for(ui i = 0;i < n;i ++) {
		ui u = id[i];
		assert(degree_start[degree[u]] == i);
		if(degree[u] > max_core) max_core = degree[u];
		core[u] = max_core;

		++ degree_start[degree[u]];
		if(degree[u] == 0) continue;

		degree_start[degree[u]-1] = degree_start[degree[u]];
		for(ui j = pstart[u];j < pend[u];j ++) if(rid[edgeList[j]] > i) {
			ui v = edgeList[j];
			ui pos1 = degree_start[degree[v]], pos2 = rid[v];
			std::swap(id[pos1], id[pos2]);
			rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
			++ degree_start[degree[v]];
			-- degree[v];

			/* non optimal method, corresponding to greedy_and_reduce function */
			// deg[u]++;   
			// orien_edges[eid[j]].to = v;
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) {
		ui cnt = 0;
		ui u = peel_sequence[i];
		for(ui j = pstart[u];j < pend[u];j ++) if(rid[edgeList[j]] > i) ++ cnt;
		assert(cnt == degree[u]);
	}
#endif

	delete[] degree;
	delete[] degree_start;
	delete[] rid;
	 
	log_info(graphClock_.Count("core decomposition time: %s s, max_core: %u", Utility::integer_to_string(timer.elapsed()).c_str(),max_core));

	return max_core;
}

void Graph::scan(ui *core, ui level, ui* curr, int* curr_tail){
	const ui BUFFER_SIZE_BYTES = 2048; //size of cache line
	const ui BUFFER_SIZE = BUFFER_SIZE_BYTES / sizeof(ui);
	ui buff[BUFFER_SIZE];
	ui index = 0;
	#pragma omp for schedule(static)
	for(ui i = 0; i < n; i++){
		if(core[i] == level){
			buff[index++] = i;
			if (index >= BUFFER_SIZE) {
                ui tempIdx = __sync_fetch_and_add(curr_tail, BUFFER_SIZE);
                for (ui j = 0; j < BUFFER_SIZE; j++) {
                    curr[tempIdx+j] = buff[j];
                }
                index = 0;
            }
		}
	}
	if (index > 0) {
        ui tempIdx = __sync_fetch_and_add(curr_tail, index);
        for (ui j = 0; j < index; j++)
            curr[tempIdx+j] = buff[j];
    }
	#pragma omp barrier
}

void Graph::process_sublevel(ui *curr, int curr_tail, ui *core, ui level, ui *next, int *next_tail){
	const ui BUFFER_SIZE_BYTES = 2048;
	const ui BUFFER_SIZE = BUFFER_SIZE_BYTES / sizeof(ui);
	ui buff[BUFFER_SIZE];
	ui index = 0;
	#pragma omp for schedule(static)
	for(ui i = 0; i < curr_tail; i++){
		ui v = curr[i];
		for(ui j = pstart[v]; j < pend[v]; j++){
			ui u = edgeList[j];
			ui deg_u = core[u];
			if (deg_u > level) {
                int du =  __sync_fetch_and_sub(&core[u], 1);
                if (du == (level+1)) {
                    buff[index++] = u;
                    if (index >= BUFFER_SIZE) {
                        ui tempIdx = __sync_fetch_and_add(next_tail, BUFFER_SIZE);
                        for(ui bufIdx = 0; bufIdx < BUFFER_SIZE; bufIdx++)
                            next [tempIdx + bufIdx] = buff[bufIdx];
                        index = 0;
                    }
                }
            }
		}
	}
	if (index > 0) {
        ui tempIdx =  __sync_fetch_and_add(next_tail, index);;
        for (ui bufIdx = 0; bufIdx < index; bufIdx++)
            next [tempIdx + bufIdx] = buff[bufIdx];
    }
#pragma omp barrier

#pragma omp for schedule(static)
    for (ui i=0; i<*next_tail; i++) {
        ui u = next[i];
        if (core[u] != level)
            core[u] = level;
    }
#pragma omp barrier
}



void Graph::process_sublevel_plus(ui *curr, int curr_tail, ui *core, ui level, ui *next, int *next_tail){
	const ui BUFFER_SIZE_BYTES = 2048;
	const ui BUFFER_SIZE = BUFFER_SIZE_BYTES / sizeof(ui);
	ui buff[BUFFER_SIZE];
	ui index = 0;
	#pragma omp for schedule(static)
	for(ui i = 0; i < curr_tail; i++){
		ui v = curr[i];
		for(ui j = original_pstart[v]; j < original_pend[v]; j++){
			ui u = original_edgeList[j];
			ui deg_u = core[u];
			if (deg_u > level) {
                int du =  __sync_fetch_and_sub(&core[u], 1);
                if (du == (level+1)) {
                    buff[index++] = u;
                    if (index >= BUFFER_SIZE) {
                        ui tempIdx = __sync_fetch_and_add(next_tail, BUFFER_SIZE);
                        for(ui bufIdx = 0; bufIdx < BUFFER_SIZE; bufIdx++)
                            next [tempIdx + bufIdx] = buff[bufIdx];
                        index = 0;
                    }
                }
            }
		}
	}
	if (index > 0) {
        ui tempIdx =  __sync_fetch_and_add(next_tail, index);;
        for (ui bufIdx = 0; bufIdx < index; bufIdx++)
            next [tempIdx + bufIdx] = buff[bufIdx];
    }
#pragma omp barrier

#pragma omp for schedule(static)
    for (ui i=0; i<*next_tail; i++) {
        ui u = next[i];
        if (core[u] != level)
            core[u] = level;
    }
#pragma omp barrier
}




ui Graph::core_decomposition_parallel(ui *peel_sequence, ui *core){
	n = original_n, m = original_m;
	ui *curr = new ui[n]();
	assert(curr != NULL);
	ui *next = new ui[n]();
	assert(next != NULL);
	int curr_tail = 0, next_tail = 0;
	int thrnmb = getRunningThreadCount();
	#pragma omp parallel
	{
		ui tid = omp_get_thread_num();
		int to_do = n;
		ui level = 0;
		#pragma omp for schedule(static)
		for(ui i = 0; i < n; i++){
			#ifdef TP
			core[i] = original_pend[i] - original_pstart[i];
			#else
			core[i] = pend[i] - pstart[i];
			#endif
		}
		while(to_do > 0){
			scan(core,level,curr,&curr_tail);
			while(curr_tail > 0){
				to_do -= curr_tail;
				#ifdef TP
				process_sublevel_plus(curr,curr_tail,core,level,next,&next_tail);
				#else
				process_sublevel(curr,curr_tail,core,level,next,&next_tail);
				#endif
				
				if(tid == 0){
					ui *tmp_curr = curr;
					curr = next;
					next = tmp_curr;
					curr_tail = next_tail;
					next_tail = 0;
				}
				#pragma omp barrier
			}
			level = level + 1;
			#pragma omp barrier
		}
	}
	delete[] curr;
	delete[] next;
	for(int i = 0; i < n; i++)
        peel_sequence[i] = i;
    sort(peel_sequence,peel_sequence+n,[&](const int& a, const int& b)->bool{
		if(core[a]<core[b])
			return true;
		else
			return false;
	});
	ui max_core = core[peel_sequence[n-1]];
	#ifndef NDEBUG
	log_info(graphClock_.Count("core decomposition max_core: %u",max_core));
	#endif
	return max_core;
}

void Graph::greedy_and_reduce(ui *peel_sequence, ui *core, ui input_mid) {
	ui *rid = new ui[n]();
	ui *degree = new ui[n]();
	ui best_n = n, best_m = m;
	ui best_idx = 0, current_m = m;
	if(input_mid == 0){
		for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
		for(ui i = 0;i < n;i ++) {
			if((lint)current_m*(lint)best_n > (lint)best_m*(lint)(n-i)) {
				best_n = n-i;
				best_m = current_m;
				best_idx = i;
				// printf("best_idx: %d, best_n: %d, best_m: %d\n",best_idx,best_n,best_m);
			}
			for(ui j = pstart[peel_sequence[i]];j < pend[peel_sequence[i]];j ++) if(rid[edgeList[j]] > i) -- current_m;
		}
		double density = 1.0*best_m/best_n;
		app_pa = ceil(density);
		// log_info(graphClock_.Count("approximate pa: %u, real core: %u",app_pa,core[peel_sequence[n-best_n]]));
	}
	else{
		// for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
		// for(ui i = 0;i < n;i ++) {
		// 	if((lint)current_m*(lint)best_n > (lint)best_m*(lint)(n-i)) {
		// 		best_n = n-i;
		// 		best_m = current_m;
		// 		best_idx = i;
		// 		double density = 1.0*best_m/best_n;
		// 		if(ceil(density) >= input_mid)
		// 			break;
		// 		// printf("best_idx: %d, best_n: %d, best_m: %d\n",best_idx,best_n,best_m);
		// 	}
		// 	for(ui j = pstart[peel_sequence[i]];j < pend[peel_sequence[i]];j ++) if(rid[edgeList[j]] > i) -- current_m;
		// }
		app_pa = input_mid;
		// log_info(graphClock_.Count("approximate pa: %u, real core: %u",app_pa,core[peel_sequence[n-best_n]]));
	}
	
	ui start_idx = 0;
	for(int i = n-1; i >=0; i--){
		ui now = peel_sequence[i];
		if(core[now] < app_pa) break;
		start_idx = i;
		for(ui j = pstart[now]; j < pend[now]; j++){
			edge& ne = orien_edges[eid[j]];
			if (ne.to != now) continue;
			unsigned tar = (ne.src == now) ? ne.end : ne.src;
			if (core[tar] < app_pa) continue;
			if (deg[now] <= deg[tar] - 2)
			{
				ne.to = tar;
				deg[now]++;
				deg[tar]--;
			}
		}
	}
    // log_info(graphClock_.Count("start_idx %u, best_idx %u",start_idx,best_idx));
	// assert(start_idx <= best_idx);

	for(int i = start_idx; i < n; i++){
		peel_sequence[i-start_idx] = peel_sequence[i];
	} 
	for(int i = 0; i < n; i++) rid[i] = n;
	n -= start_idx;
	vector<ui> ids;
	for(int i = 0; i < n; i++) {
		ids.push_back(peel_sequence[i]);
	}
	sort(ids.begin(),ids.end());
	for(ui i = 0;i < n;i ++) rid[ids[i]] = i;   // correspond to vertex id in the reduced graph
	for(ui i = 0;i < n;i ++) {
		// printf("i: %u, ids[i]: %u\n",i,ids[i]);
		ori_map[i] = ids[i];
		peel_sequence[i] = rid[peel_sequence[i]];
		core[i] = core[ids[i]];
	}
	for(ui i = 0;i < res.size();i ++) res[i] = rid[res[i]];	

	// get the degree of vertices in the reduced graph
	m = 0;
	for(ui i = 0;i < n;i ++) {
		degree[i] = 0;
		for(ui j = pstart[ids[i]];j < pend[ids[i]];j ++) if(rid[edgeList[j]] < n) ++ degree[i];
		m += degree[i];		
	}
	m /= 2;

	ui *tmp_edgeList = new ui[2*m]();
	ui *tmp_eid = new ui[2*m]();
	edge *tmp_orien_edges = new edge[m];
	ui *tmp_pstart = new ui[n]();
	ui *tmp_pend = new ui[n]();
	ui *tmp_deg = new ui[n]();

	tmp_pstart[0] = 0;
	tmp_pend[0] = 0;
	for(int i = 1; i < n; i++){
		tmp_pstart[i] = tmp_pstart[i-1]+degree[i-1];
		tmp_pend[i] = tmp_pstart[i];
	}

	ui id = 0;
	for(int i = 0; i < n; i++){
		tmp_pend[i] = tmp_pstart[i];
		ui old_start = pstart[ids[i]], old_end = pend[ids[i]];
		for(ui j = old_start;j < old_end;j ++) {
			ui v = rid[edgeList[j]];
			if(v < n && v > i) {
				if(orien_edges[eid[j]].src == ori_map[i])
					tmp_orien_edges[id].src = i, tmp_orien_edges[id].end = v;
				else
					tmp_orien_edges[id].src = v, tmp_orien_edges[id].end = i;
				if(orien_edges[eid[j]].to == ori_map[i])
					tmp_orien_edges[id].to = i;
				else
					tmp_orien_edges[id].to = v;
				tmp_eid[tmp_pend[i]] = id;
				tmp_edgeList[tmp_pend[i]++] = v;
				tmp_eid[tmp_pstart[v]] = id++;
				tmp_edgeList[tmp_pstart[v]++] = i;
			}
		}
	}
	for(int i = 0; i < n; i++){
		tmp_pstart[i] = tmp_pend[i]-degree[i];
		tmp_deg[i] = deg[ori_map[i]];
	}
	#ifndef NDEBUG
	log_info(graphClock_.Count("reduced subgraph v: %u, m: %u",n,m));	
	#endif
	
	// for(int i = 0; i < n; i++){
	// 	ui u = i;
	// 	for(int j = tmp_pstart[u]; j < tmp_pend[u]; j++){
	// 		ui v = tmp_edgeList[j];
	// 		if(i > v) continue;
	// 		ui id = tmp_eid[j];
	// 		// printf("u: %u, v: %u, deg[u]: %u, eid: %u, ori->u: %u, ori->v: %u, ori->to: %u\n",ori_map[u],ori_map[v],deg[ori_map[u]],id,orien_edges[id].src,orien_edges[id].end,orien_edges[id].to);
	// 	}
	// }

	delete []rid;
	delete []degree;

	delete []edgeList;
	delete []eid;
	delete []pstart;
	delete []pend; 
	delete []deg;
	delete []orien_edges;
		
	edgeList = tmp_edgeList;
	tmp_edgeList = nullptr;
	eid = tmp_eid;
	tmp_eid = nullptr;
	pstart = tmp_pstart;
	tmp_pstart = nullptr;
	pend = tmp_pend;
	tmp_pend = nullptr;
	deg = tmp_deg;
	tmp_deg = nullptr;
	orien_edges = tmp_orien_edges;
	tmp_orien_edges = nullptr;
}



void Graph::greedy_and_reduce_plus(ui *peel_sequence, ui *core, ui input_mid) {
	ui *rid = new ui[n]();
	ui *degree = new ui[n]();
	ui best_n = n, best_m = m;
	ui best_idx = 0, current_m = m;
	if(input_mid == 0){
		for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
		for(ui i = 0;i < n;i ++) {
			if((lint)current_m*(lint)best_n > (lint)best_m*(lint)(n-i)) {
				best_n = n-i;
				best_m = current_m;
				best_idx = i;
			}
			for(ui j = pstart[peel_sequence[i]];j < pend[peel_sequence[i]];j ++) if(rid[edgeList[j]] > i) -- current_m;
		}
		double density = 1.0*best_m/best_n;
		app_pa = ceil(density);
		// log_info(graphClock_.Count("approximate pa: %u, real core: %u",app_pa,core[peel_sequence[n-best_n]]));
	}
	else{
		app_pa = input_mid;
	}
	
	ui start_idx = 0;
	for(int i = n-1; i >=0; i--){
		ui now = peel_sequence[i];
		if(core[now] < app_pa) break;
		start_idx = i;
	}

	for(int i = start_idx; i < n; i++)
		peel_sequence[i-start_idx] = peel_sequence[i];
	 
	for(int i = 0; i < n; i++) rid[i] = n;
	n -= start_idx;
	vector<ui> ids; // represent original graph vertex id
	for(int i = 0; i < n; i++) {
		ids.push_back(peel_sequence[i]);
	}
	sort(ids.begin(),ids.end());
	for(ui i = 0;i < n;i ++) rid[ids[i]] = i;   // correspond to vertex id in the reduced graph
	for(ui i = 0;i < n;i ++) {
		ori_map[i] = ids[i];
		peel_sequence[i] = rid[peel_sequence[i]];
		core[i] = core[ids[i]];
	}

	// get the degree of vertices in the reduced graph
	m = 0;
	for(ui i = 0;i < n;i ++) {
		degree[i] = 0;
		for(ui j = pstart[ids[i]];j < pend[ids[i]];j ++) if(rid[edgeList[j]] < n) ++ degree[i];
		m += degree[i];		
	}
	m /= 2;

	ui *tmp_edgeList = new ui[2*m]();
	ui *tmp_eid = new ui[2*m]();
	edge *tmp_orien_edges = new edge[m];
	ui *tmp_pstart = new ui[n]();
	ui *tmp_pend = new ui[n]();
	ui *tmp_deg = new ui[n]();

	tmp_pstart[0] = 0;
	tmp_pend[0] = 0;
	for(int i = 1; i < n; i++){
		tmp_pstart[i] = tmp_pstart[i-1]+degree[i-1];
		tmp_pend[i] = tmp_pstart[i];
	}

	ui id = 0;
	for(int i = 0; i < n; i++){
		tmp_pend[i] = tmp_pstart[i];
		ui old_start = pstart[ids[i]], old_end = pend[ids[i]];
		for(ui j = old_start;j < old_end;j ++) {
			ui v = rid[edgeList[j]];
			if(v < n && v > i) {
				tmp_eid[tmp_pend[i]] = id;
				tmp_edgeList[tmp_pend[i]++] = v;
				tmp_eid[tmp_pstart[v]] = id++;
				tmp_edgeList[tmp_pstart[v]++] = i;
			}
		}
	}
	for(int i = 0; i < n; i++){
		tmp_pstart[i] = tmp_pend[i]-degree[i];
	}
	#ifndef NDEBUG
	log_info(graphClock_.Count("reduced subgraph v: %u, m: %u",n,m));
	#endif

	bool *is_in = new bool[n]();

	for(ui i = 0; i < n; i++){
		ui u = peel_sequence[i];
		for(ui j = tmp_pstart[u]; j < tmp_pend[u]; j++){
			ui v = tmp_edgeList[j];
			if(is_in[v]) continue;
			tmp_deg[u]++;
			tmp_orien_edges[tmp_eid[j]] = (edge){.src = u, .end = v, .to = v};
		}
		is_in[u] = true;
	}
	delete[] is_in;

	for(int i = n-1; i >=0; i--){
		ui now = peel_sequence[i];
		for(ui j = tmp_pstart[now]; j < tmp_pend[now]; j++){
			edge& ne = tmp_orien_edges[tmp_eid[j]];
			if (ne.to != now) continue;
			ui tar = (ne.src == now) ? ne.end : ne.src;
			if (tmp_deg[now] <= tmp_deg[tar] - 2)
			{
				ne.to = tar;
				tmp_deg[now]++;
				tmp_deg[tar]--;
			}
		}
	}

	delete []rid;
	delete []degree;
	delete []edgeList;
	delete []pstart;
	delete []pend; 
	delete []deg;
	if(orien_edges != nullptr)
		delete []orien_edges;
	if(eid != nullptr)
		delete []eid;
		
	edgeList = tmp_edgeList;
	tmp_edgeList = nullptr;
	eid = tmp_eid;
	tmp_eid = nullptr;
	pstart = tmp_pstart;
	tmp_pstart = nullptr;
	pend = tmp_pend;
	tmp_pend = nullptr;
	deg = tmp_deg;
	tmp_deg = nullptr;
	orien_edges = tmp_orien_edges;
	tmp_orien_edges = nullptr;

}

// generate reduced subgraph based on original_graph info
void Graph::greedy_and_reduce_plus_plus(ui *peel_sequence, ui *core, ui input_mid) {
	n = original_n, m = original_m;
	ui *rid = new ui[n]();
	ui *degree = new ui[n]();
	ui best_n = n, best_m = m;
	ui best_idx = 0, current_m = m;
	if(input_mid == 0){
		for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
		for(ui i = 0;i < n;i ++) {
			if((lint)current_m*(lint)best_n > (lint)best_m*(lint)(n-i)) {
				best_n = n-i;
				best_m = current_m;
				best_idx = i;
			}
			for(ui j = original_pstart[peel_sequence[i]];j < original_pend[peel_sequence[i]];j ++) if(rid[original_edgeList[j]] > i) -- current_m;
		}
		double density = 1.0*best_m/best_n;
		app_pa = ceil(density);
		// log_info(graphClock_.Count("approximate pa: %u, real core: %u",app_pa,core[peel_sequence[n-best_n]]));
	}
	else{
		app_pa = input_mid;
	}
	
	ui start_idx = 0;
	for(int i = n-1; i >=0; i--){
		ui now = peel_sequence[i];
		if(core[now] < app_pa) break;
		start_idx = i;
	}

	for(int i = start_idx; i < n; i++)
		peel_sequence[i-start_idx] = peel_sequence[i];
	 
	for(int i = 0; i < n; i++) rid[i] = n;
	n -= start_idx;
	vector<ui> ids; // represent original graph vertex id
	for(int i = 0; i < n; i++) {
		ids.push_back(peel_sequence[i]);
	}
	sort(ids.begin(),ids.end());
	for(ui i = 0;i < n;i ++) rid[ids[i]] = i;   // correspond to vertex id in the reduced graph
	for(ui i = 0;i < n;i ++) {
		ori_map[i] = ids[i];
		peel_sequence[i] = rid[peel_sequence[i]];
		core[i] = core[ids[i]];
	}

	// get the degree of vertices in the reduced graph
	m = 0;
	for(ui i = 0;i < n;i ++) {
		degree[i] = 0;
		for(ui j = original_pstart[ids[i]];j < original_pend[ids[i]];j ++) if(rid[original_edgeList[j]] < n) ++ degree[i];
		m += degree[i];		
	}
	m /= 2;

	ui *tmp_edgeList = new ui[2*m]();
	ui *tmp_eid = new ui[2*m]();
	edge *tmp_orien_edges = new edge[m];
	ui *tmp_pstart = new ui[n]();
	ui *tmp_pend = new ui[n]();
	ui *tmp_deg = new ui[n]();

	tmp_pstart[0] = 0;
	tmp_pend[0] = 0;
	for(int i = 1; i < n; i++){
		tmp_pstart[i] = tmp_pstart[i-1]+degree[i-1];
		tmp_pend[i] = tmp_pstart[i];
	}

	ui id = 0;
	for(int i = 0; i < n; i++){
		tmp_pend[i] = tmp_pstart[i];
		ui old_start = original_pstart[ids[i]], old_end = original_pend[ids[i]];
		for(ui j = old_start;j < old_end;j ++) {
			ui v = rid[original_edgeList[j]];
			if(v < n && v > i) {
				tmp_eid[tmp_pend[i]] = id;
				tmp_edgeList[tmp_pend[i]++] = v;
				tmp_eid[tmp_pstart[v]] = id++;
				tmp_edgeList[tmp_pstart[v]++] = i;
			}
		}
	}
	for(int i = 0; i < n; i++){
		tmp_pstart[i] = tmp_pend[i]-degree[i];
	}
	#ifndef NDEBUG
	log_info(graphClock_.Count("reduced subgraph v: %u, m: %u",n,m));
	#endif

	bool *is_in = new bool[n]();

	for(ui i = 0; i < n; i++){
		ui u = peel_sequence[i];
		for(ui j = tmp_pstart[u]; j < tmp_pend[u]; j++){
			ui v = tmp_edgeList[j];
			if(is_in[v]) continue;
			tmp_deg[u]++;
			tmp_orien_edges[tmp_eid[j]] = (edge){.src = u, .end = v, .to = v};
		}
		is_in[u] = true;
	}
	delete[] is_in;

	for(int i = n-1; i >=0; i--){
		ui now = peel_sequence[i];
		for(ui j = tmp_pstart[now]; j < tmp_pend[now]; j++){
			edge& ne = tmp_orien_edges[tmp_eid[j]];
			if (ne.to != now) continue;
			ui tar = (ne.src == now) ? ne.end : ne.src;
			if (tmp_deg[now] <= tmp_deg[tar] - 2)
			{
				ne.to = tar;
				tmp_deg[now]++;
				tmp_deg[tar]--;
			}
		}
	}

	delete []rid;
	delete []degree;

	// delete []edgeList;
	// delete []pstart;
	// delete []pend; 
	// delete []deg;

	if(orien_edges != nullptr)
		delete []orien_edges;
	if(eid != nullptr)
		delete []eid;
	if(deg != nullptr)
		delete []deg;
		
	edgeList = tmp_edgeList;
	tmp_edgeList = nullptr;
	eid = tmp_eid;
	tmp_eid = nullptr;
	pstart = tmp_pstart;
	tmp_pstart = nullptr;
	pend = tmp_pend;
	tmp_pend = nullptr;
	deg = tmp_deg;
	tmp_deg = nullptr;
	orien_edges = tmp_orien_edges;
	tmp_orien_edges = nullptr;
}


void Graph::malloc_glocal_info(ui n_, ui m_){
	if(edgeList == nullptr)
		edgeList = new ui[2*m_]();
	if(eid == nullptr)
		eid = new ui[2*m_]();
	if(orien_edges == nullptr)
		orien_edges = new edge[m_];
	if(pstart == nullptr)
		pstart = new ui[n_]();
	if(pend == nullptr)
		pend = new ui[n_]();
	if(deg == nullptr)
		deg = new ui[n_]();
	if(weigh == nullptr)
		weigh = new ui[n_]();
	if(vis == nullptr)
		vis = new bool[n_]();
	m = original_m = m_, n = original_n = n_;
	#ifndef NDEBUG
	log_info(graphClock_.Count("reduced subgraph v: %u, m: %u",n,m));	
	#endif
}

int Graph::getRunningThreadCount(){
    int threadCount = 0;
    std::ifstream file("/proc/self/status");
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("Threads:") == 0) {
            threadCount = std::stoi(line.substr(9));
            break;
        }
    }
    return threadCount;
}

void Graph::transfer_graphinfo_parallel(Graph *g){
	ui *mp_arr = new ui[g->original_n]();
	std::vector<ui> local_count(omp_get_max_threads(),0);
	// Step 1: Count the number of nodes to be assigned by each thread
	#pragma omp parallel
	{
		ui tid = omp_get_thread_num();
		#pragma omp for schedule(static)
		for(ui i = 0; i < g->original_n; i++){
			if(!g->vis[i])
				local_count[tid]++;
		}
	}
	// Step 2: Compute prefix sums to determine the starting index for each thread
    std::vector<ui> prefix_sums(omp_get_max_threads(), 0);
    for (int i = 1; i < omp_get_max_threads(); i++) {
        prefix_sums[i] = prefix_sums[i - 1] + local_count[i - 1];
    }
    original_n = n = prefix_sums[omp_get_max_threads() - 1] + local_count[omp_get_max_threads() - 1];
	// Step 3: Assign numbers based on the prefix sums
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        ui start_index = prefix_sums[tid];
        #pragma omp for schedule(static)
        for (ui i = 0; i < g->original_n; i++) {
            if (!g->vis[i]) {
                mp_arr[i] = start_index++;
            }
        }
    }

	if(pstart != nullptr) delete[] pstart;
	pstart = new ui[original_n]();
	if(pend != nullptr) delete[] pend;
	pend = new ui[original_n]();
    if(deg != nullptr) delete[] deg;
	deg = new ui[original_n]();
    if(weigh != nullptr) delete[] weigh;
	weigh = new ui[original_n]();
    if(vis != nullptr) delete[] vis;
	vis = new bool[original_n]();

    pstart[0] = 0;
    #pragma omp parallel for schedule(static)
    for (ui i = 0; i < g->original_n; i++) {
        if (!g->vis[i]) {
            ui count = 0;
            for (ui j = g->original_pstart[i]; j < g->original_pend[i]; j++) {
                if (!g->vis[g->original_edgeList[j]]) {
                    count++;
                }
            }
            deg[mp_arr[i]] = count;
        }
    }

    // Step 3: Compute prefix sums for pstart
    for (ui i = 1; i < original_n; i++) {
        pstart[i] = pstart[i - 1] + deg[i - 1];
    }
	original_m = m = (pstart[original_n - 1] + deg[original_n - 1])/2;
    // Step 4: Fill pend and curedge in parallel

	if(edgeList != nullptr) delete[] edgeList;
	edgeList = new ui[(pstart[original_n - 1] + deg[original_n - 1])];

    #pragma omp parallel for schedule(dynamic)
    for (ui i = 0; i < g->original_n; i++) {
        if (!g->vis[i]) {
            ui vc = mp_arr[i];
            weigh[vc] = g->weigh[i];
            pend[vc] = pstart[vc];
            for (ui j = g->original_pstart[i]; j < g->original_pend[i]; j++) {
                if (!g->vis[g->original_edgeList[j]]) {
                    edgeList[pend[vc]] = mp_arr[g->original_edgeList[j]];
                    pend[vc]++;
                }
            }
        }
    }


    if(original_pstart != nullptr) delete[] original_pstart;
	original_pstart = new ui[original_n]();
	if(original_pend != nullptr) delete[] original_pend;
	original_pend = new ui[original_n]();
    if(original_edgeList != nullptr) delete[] original_edgeList;
	original_edgeList = new ui[2*original_m];
    memcpy(original_pstart,pstart,original_n*sizeof(ui));
    memcpy(original_pend,pend,original_n*sizeof(ui));
    memcpy(original_edgeList,edgeList,2*original_m*sizeof(ui));
	#ifndef NDEBUG
    log_info(graphClock_.Count("partitioned graph v: %u, m: %u",n,m));
	#endif
	delete[] mp_arr;
}



void Graph::transfer_graphinfo_parallel_plus(Graph *g){
	ui *mp_arr = new ui[g->original_n]();
	std::vector<ui> local_count(omp_get_max_threads(),0);
	// Step 1: Count the number of nodes to be assigned by each thread
	#pragma omp parallel
	{
		ui tid = omp_get_thread_num();
		#pragma omp for schedule(static)
		for(ui i = 0; i < g->original_n; i++){
			if(!g->vis[i])
				local_count[tid]++;
		}
	}
	// Step 2: Compute prefix sums to determine the starting index for each thread
    std::vector<ui> prefix_sums(omp_get_max_threads(), 0);
    for (int i = 1; i < omp_get_max_threads(); i++) {
        prefix_sums[i] = prefix_sums[i - 1] + local_count[i - 1];
    }
    original_n = n = prefix_sums[omp_get_max_threads() - 1] + local_count[omp_get_max_threads() - 1];
	// Step 3: Assign numbers based on the prefix sums
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        ui start_index = prefix_sums[tid];
        #pragma omp for schedule(static)
        for (ui i = 0; i < g->original_n; i++) {
            if (!g->vis[i]) {
                mp_arr[i] = start_index++;
            }
        }
    }

	// if(pstart != nullptr) delete[] pstart;
	// pstart = new ui[original_n]();
	// if(pend != nullptr) delete[] pend;
	// pend = new ui[original_n]();


    if(deg != nullptr) delete[] deg;
	deg = new ui[original_n]();
	if(original_pstart != nullptr) delete[] original_pstart;
	original_pstart = new ui[original_n]();
	if(original_pend != nullptr) delete[] original_pend;
	original_pend = new ui[original_n]();
    if(weigh != nullptr) delete[] weigh;
	weigh = new ui[original_n]();
    if(vis != nullptr) delete[] vis;
	vis = new bool[original_n]();

    original_pstart[0] = 0;
    #pragma omp parallel for schedule(static)
    for (ui i = 0; i < g->original_n; i++) {
        if (!g->vis[i]) {
            ui count = 0;
            for (ui j = g->original_pstart[i]; j < g->original_pend[i]; j++) {
                if (!g->vis[g->original_edgeList[j]]) {
                    count++;
                }
            }
            deg[mp_arr[i]] = count;
        }
    }

    // Step 3: Compute prefix sums for pstart
    for (ui i = 1; i < original_n; i++) {
        original_pstart[i] = original_pstart[i - 1] + deg[i - 1];
    }
	original_m = m = (original_pstart[original_n - 1] + deg[original_n - 1])/2;

	if(original_edgeList != nullptr) delete[] original_edgeList;
	original_edgeList = new ui[2*original_m];
    // Step 4: Fill pend and curedge in parallel

    #pragma omp parallel for schedule(dynamic)
    for (ui i = 0; i < g->original_n; i++) {
        if (!g->vis[i]) {
            ui vc = mp_arr[i];
            weigh[vc] = g->weigh[i];
            original_pend[vc] = original_pstart[vc];
            for (ui j = g->original_pstart[i]; j < g->original_pend[i]; j++) {
                if (!g->vis[g->original_edgeList[j]]) {
                    original_edgeList[original_pend[vc]] = mp_arr[g->original_edgeList[j]];
                    original_pend[vc]++;
                }
            }
        }
    }


    // if(original_pstart != nullptr) delete[] original_pstart;
	// original_pstart = new ui[original_n]();
	// if(original_pend != nullptr) delete[] original_pend;
	// original_pend = new ui[original_n]();
    // if(original_edgeList != nullptr) delete[] original_edgeList;
	// original_edgeList = new ui[2*original_m];
    // memcpy(original_pstart,pstart,original_n*sizeof(ui));
    // memcpy(original_pend,pend,original_n*sizeof(ui));
    // memcpy(original_edgeList,edgeList,2*original_m*sizeof(ui));
	#ifndef NDEBUG
    log_info(graphClock_.Count("partitioned graph v: %u, m: %u",n,m));
	#endif
	delete[] mp_arr;
}




void Graph::transfer_graphinfo_parallel_right(Graph *g){
	ui *mp_arr = new ui[g->original_n]();
	std::vector<ui> local_count(omp_get_max_threads(),0);
	// Step 1: Count the number of nodes to be assigned by each thread
	#pragma omp parallel
	{
		ui tid = omp_get_thread_num();
		#pragma omp for schedule(static)
		for(ui i = 0; i < g->n; i++){
			if(!g->vis[i])
				local_count[tid]++;
		}
	}
	// Step 2: Compute prefix sums to determine the starting index for each thread
    std::vector<ui> prefix_sums(omp_get_max_threads(), 0);
    for (int i = 1; i < omp_get_max_threads(); i++) {
        prefix_sums[i] = prefix_sums[i - 1] + local_count[i - 1];
    }
    original_n = n = prefix_sums[omp_get_max_threads() - 1] + local_count[omp_get_max_threads() - 1];
	// Step 3: Assign numbers based on the prefix sums
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        ui start_index = prefix_sums[tid];
        #pragma omp for schedule(static)
        for (ui i = 0; i < g->n; i++) {
            if (!g->vis[i]) {
                mp_arr[i] = start_index++;
            }
        }
    }



    if(deg != nullptr) delete[] deg;
	deg = new ui[original_n]();
	if(original_pstart != nullptr) delete[] original_pstart;
	original_pstart = new ui[original_n]();
	if(original_pend != nullptr) delete[] original_pend;
	original_pend = new ui[original_n]();
    if(weigh != nullptr) delete[] weigh;
	weigh = new ui[original_n]();
    if(vis != nullptr) delete[] vis;
	vis = new bool[original_n]();

    original_pstart[0] = 0;
    #pragma omp parallel for schedule(static)
    for (ui i = 0; i < g->n; i++) {
        if (!g->vis[i]) {
            ui count = 0;
            for (ui j = g->pstart[i]; j < g->pend[i]; j++) {
                if (!g->vis[g->edgeList[j]]) {
                    count++;
                }
            }
            deg[mp_arr[i]] = count;
        }
    }

    // Step 3: Compute prefix sums for pstart
    for (ui i = 1; i < original_n; i++) {
        original_pstart[i] = original_pstart[i - 1] + deg[i - 1];
    }
	original_m = m = (original_pstart[original_n - 1] + deg[original_n - 1])/2;

	if(original_edgeList != nullptr) delete[] original_edgeList;
	original_edgeList = new ui[2*original_m];
    // Step 4: Fill pend and curedge in parallel

    #pragma omp parallel for schedule(dynamic)
    for (ui i = 0; i < g->n; i++) {
        if (!g->vis[i]) {
            ui vc = mp_arr[i];
            weigh[vc] = g->weigh[g->ori_map[i]];
            original_pend[vc] = original_pstart[vc];
            for (ui j = g->pstart[i]; j < g->pend[i]; j++) {
                if (!g->vis[g->edgeList[j]]) {
                    original_edgeList[original_pend[vc]] = mp_arr[g->edgeList[j]];
                    original_pend[vc]++;
                }
            }
        }
    }

	#ifndef NDEBUG
    log_info(graphClock_.Count("partitioned graph v: %u, m: %u",n,m));
	#endif
	delete[] mp_arr;
}



void Graph::transfer_graphinfo(Graph *g){
    unordered_map<ui,ui> mp_;
    ui vc = 0;
    for(ui i = 0; i < g->original_n; i++)
        if(!g->vis[i]){
            n++;
            mp_[i] = vc++;
        }
    
	if(pstart != nullptr) delete[] pstart;
	pstart = new ui[original_n]();
	if(pend != nullptr) delete[] pend;
	pend = new ui[original_n]();
    if(deg != nullptr) delete[] deg;
	deg = new ui[original_n]();
    if(weigh != nullptr) delete[] weigh;
	weigh = new ui[original_n]();
    if(vis != nullptr) delete[] vis;
	vis = new bool[original_n]();

    std::vector<ui> curedge;
    pstart[0] = 0;
    vc = 0;
    for(ui i = 0; i < g->original_n; i++){
		if(!g->vis[i]) {
            assert(vc == mp_[i]);
            weigh[vc] = g->weigh[i];
            pend[vc] = pstart[vc];
			for(ui j = g->original_pstart[i]; j < g->original_pend[i]; j++){
				assert(g->original_edgeList[j] < g->original_n);
				if(!g->vis[g->original_edgeList[j]]){
                    curedge.push_back(mp_[g->original_edgeList[j]]);
					pend[vc]++;
				}
			}
            if(vc < original_n-1)
			    pstart[vc+1] = pend[vc];
            deg[vc] = pend[vc]-pstart[vc];
            vc++;
		}
	}
    assert(vc == original_n);
    if(edgeList != nullptr) delete[] edgeList;
	edgeList = new ui[curedge.size()];
    memcpy(edgeList,curedge.data(),curedge.size()*sizeof(ui));
    m = original_m = curedge.size()/2;

    // if(orien_edges != nullptr) delete[] orien_edges;
	// orien_edges = new edge[m];
	// if(eid != nullptr) delete[] eid;
	// eid = new ui[m*2]();

	// ui idx = 0;
	// for(ui i = 0; i < n; i++){
	// 	for(ui j = pstart[i]; j < pend[i]; j++){
	// 		ui v = edgeList[j];
	// 		if(v < i) continue;
	// 		eid[j] = idx;
	// 		eid[pstart[v]++] = idx;
	// 		orien_edges[idx].src = i;
	// 		orien_edges[idx].end = v;
	// 		orien_edges[idx].to = v;
	// 		idx++;
	// 	}
	// }

	// for(int i = 0; i < n; i++){
	// 	pstart[i] = pend[i] - deg[i];
	// 	deg[i] = 0;
	// }

    if(original_pstart != nullptr) delete[] original_pstart;
	original_pstart = new ui[original_n]();
	if(original_pend != nullptr) delete[] original_pend;
	original_pend = new ui[original_n]();
    if(original_edgeList != nullptr) delete[] original_edgeList;
	original_edgeList = new ui[2*original_m];
    memcpy(original_pstart,pstart,original_n*sizeof(ui));
    memcpy(original_pend,pend,original_n*sizeof(ui));
    memcpy(original_edgeList,edgeList,2*original_m*sizeof(ui));
	#ifndef NDEBUG
    log_info(graphClock_.Count("partitioned graph v: %u, m: %u",n,m));
	#endif

}

void Graph::extract_subgraph_from_original(){
	if(edgeList != nullptr) delete[] edgeList;
	edgeList = new ui[2*original_m];
	if(pstart != nullptr) delete[] pstart;
	pstart = new ui[original_n]();
	if(pend != nullptr) delete[] pend;
	pend = new ui[original_n]();
    if(deg != nullptr) delete[] deg;
	deg = new ui[original_n]();
	
    ui vc = 0, ec = 0;
    if(min_wei != INT32_MAX){
        for(ui i = 0; i < original_n; i++){
            pend[i] = pstart[i];
            if(!vis[i]) {
                vc++;
                for(ui j = original_pstart[i]; j < original_pend[i]; j++){
                    assert(original_edgeList[j] < original_n);
                    if(!vis[original_edgeList[j]]){
                        edgeList[ec++] = original_edgeList[j];
                        pend[i]++;
                    }
                }
            }
            if(i < original_n-1)
                pstart[i+1] = pend[i];
            deg[i] = pend[i]-pstart[i];
        }
        n = original_n;
	    m = ec/2;
    }
    else{
        memcpy(edgeList,original_edgeList,2*original_m*sizeof(ui));
        memcpy(pstart,original_pstart,original_n*sizeof(ui));
        memcpy(pend,original_pend,original_n*sizeof(ui));
        n = original_n;
        m = original_m;
    }

	#pragma omp parallel for
    for(ui i = 0; i < n; i++)
        deg[i] = pend[i] - pstart[i];

	#ifndef NDEBUG
	log_info(graphClock_.Count("Update original graph to current graph"));
	#endif

	// if(orien_edges != nullptr) delete[] orien_edges;
	// orien_edges = new edge[m];
	// if(eid != nullptr) delete[] eid;
	// eid = new ui[m*2]();

	// ui idx = 0;
	// // #pragma omp parallel for schedule(dynamic)
	// for(ui i = 0; i < n; i++){
	// 	for(ui j = pstart[i]; j < pend[i]; j++){
	// 		ui v = edgeList[j];
	// 		if(v < i) continue;
	// 		// #pragma omp critical
	// 		{
	// 			eid[j] = idx;
	// 			eid[pstart[v]] = idx;
	// 			pstart[v]++;
	// 			orien_edges[idx] = (edge){.src = i, .end = v, .to = v};
	// 			// orien_edges[idx].src = i;
	// 			// orien_edges[idx].end = v;
	// 			// orien_edges[idx].to = v;
	// 			idx++;
	// 		}
	// 	}
	// }

	// #pragma omp parallel for
	// for(int i = 0; i < n; i++){
	// 	pstart[i] = pend[i] - deg[i];
    //     // printf("start: %u, end: %u, deg: %u   ",pstart[i],pend[i],deg[i]);
	// 	deg[i] = 0;
	// }
	// log_info(graphClock_.Count("Update original graph to current graph"));
}



void Graph::exact_pseudoarboricity_computation(){
	p = new ui[n]();
	cur_vis = new bool[n]();
	if(pa == 0){
		test_value = pa;
		ui p_upper = get_max_d() + 1, p_lower = app_pa;
		ui ori_upper = p_upper;
		while (p_upper > p_lower)
		{
			test_value = (p_upper + p_lower) / 2;
			if (ReTest())
				p_upper = test_value;
			else
				p_lower = test_value + 1;
		}
		pa = max_pa = p_lower;
		#ifndef NDEBUG
        log_info(graphClock_.Count("p_upper: %u, p_lower: %u, max_pa: %u",ori_upper,p_lower,max_pa));
		#endif
	}

    // if(pa == 0){
	// 	test_value = app_pa-1;
    //     unsigned p_upper = get_max_d() + 1, p_lower = app_pa;
	// 	while (!ReTest())
	// 	{
    //         unordered_set<ui> isContain_;
    //         find_integer_densest_sub(isContain_);
    //         ui m_count = 0, cur_min_wei = INT32_MAX;
    //         for(auto ii : isContain_){
    //             cur_min_wei = min(weigh[ori_map[ii]],cur_min_wei);
    //             for(int j = pstart[ii]; j < pend[ii]; j++){
    //                 ui v = edgeList[j];
    //                 if(v < ii) continue;
    //                 if(isContain_.find(v) != isContain_.end() && !cur_vis[v]) m_count++;
    //             }
    //         }
    //         printf("ver: %u, edge: %u, integer density: %f, test_value: %u\n",isContain_.size(),m_count,1.0*m_count/isContain_.size(),test_value);
    //         test_value++;
	// 	}
    //     log_info(graphClock_.Count("p_upper: %u, p_lower: %u, max_pa: %u",get_max_d() + 1,p_lower,test_value));
	// }
	
	// find subgraph with maximum integer density value
	test_value = max_pa-1;
	// find_maximum_weight_subgraph_peeling();  //based peeling method
	// find_maximum_weight_subgraph_2approximate();  //based 2-approximate method
	find_maximum_weight_subgraph_2approximate_binary();  //based 2-approximate and binary
	delete[] p;
	delete[] cur_vis;
}



void Graph::construct_subgraph(std::vector<ui> &vertexSet, ui pre_count){
	ui cur_count = vertexSet.size();
	// printf("cur: %u, pre: %u\n",cur_count,pre_count);
	for(ui i = pre_count; i < cur_count; i++){
		ui u = vertexSet[i];
		cur_vis[u] = false;
		for(int j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(cur_vis[v]) continue;
			ui id = eid[j];
			edge e = orien_edges[id];
			ui from = e.to == e.end ? e.src : e.end;
			deg[from]++;
			// printf("u: %u, v: %u, from: %u, deg: %u\n",u,v,from,deg[from]);
		}
	}

}


void Graph::find_maximum_weight_subgraph_peeling(){
	unordered_set<ui> isContain;
	ui weight = 0;
	ui current_n = n;
	while(1){
		if(ReTest()) break;
		isContain.clear();
		find_integer_densest_sub(isContain);
		min_wei = INT32_MAX;
		for(auto u : isContain){
			// printf("v:%u, w:%u ",u,weigh[ori_map[u]]);
			min_wei = min(min_wei,weigh[ori_map[u]]);
		}
		// printf("remaining vertex in current satisfied graph: %u, isContain.size: %u\n",current_n,isContain.size());

		// ui m_count = 0;
		// for(int i = 0; i < n; i++){
		// 	if(isContain.find(i) != isContain.end() && !cur_vis[i]){
		// 		for(int j = pstart[i]; j < pend[i]; j++){
		// 			ui v = edgeList[j];
		// 			if(v < i) continue;
		// 			if(isContain.find(v) != isContain.end() && !cur_vis[v]) m_count++;
		// 		}
		// 	}
		// }
		// printf("min_wei: %u, v_c: %u, m_c: %u, den: %f\n",min_wei,isContain.size(),m_count,1.0*m_count/isContain.size());

		
		for(int i = 0; i < n; i++){
			if(!cur_vis[i] && weigh[ori_map[i]] <= min_wei){
				for(int j = pstart[i]; j < pend[i]; j++){
					ui v = edgeList[j];
					if(cur_vis[v]) continue;
					ui id = eid[j];
					edge e = orien_edges[id];
					if(e.to == i){
						ui from = e.to == e.end ? e.src : e.end;
						assert(deg[from]>0);
						deg[from]--;
					}
				}
				deg[i] = 0;
				cur_vis[i] = true;
				current_n--;
			}
		}

		//delete vertex with weight no greater than min_wei from original graph
		for(int i = 0; i < original_n; i++)
			if(!vis[i] && weigh[i] <= min_wei)
				vis[i] = true;			
	}
	#ifndef NDEBUG
	printf("pseudoarboricity: %u, maxweight: %u, subgraph.v: %u\n",pa,min_wei,isContain.size());	
	#endif
}
void Graph::find_maximum_weight_subgraph_2approximate_binary(){
	isContain.clear();
	ui weight = 0;
	min_wei = INT32_MAX;
	ui current_n = n;
	if(!ReTest()){
		find_integer_densest_sub(isContain);
		for(int i = 0; i < n; i++) cur_vis[i] = true;  //search for optimal subgraph
		priority_queue<pair<ui,ui>> pq;
		for(auto u : isContain){
			pq.push({weigh[ori_map[u]],u});	
            // printf("u: %u, w: %u, deg[v]: %u\n",u,weigh[ori_map[u]],deg[u]);
		}
		vector<ui> cur_verset;
		ui cur_count = 0, pre_count = 0;
		bool flag = false;
		memset(deg,0,sizeof(ui)*n);
		while (!pq.empty()) {
        	std::pair<ui, ui> top = pq.top();
			pq.pop();
			cur_verset.push_back(top.second);
			cur_count++;
			if(isPowerOfTwo(cur_count) || pq.empty()){
				construct_subgraph(cur_verset,pre_count);
				if(!ReTest()) 
				{
					flag = true;
					break;
				}
				pre_count = cur_count;
			}
			if(flag)
				break;
		}
		if(flag){
			// we adopt binary search to delete vertices
			ui left = pre_count+1, right = cur_count-1, mid = 0, last_mid = 0;
			bool inside = false;
			// log_debug(graphClock_.Count("left: %u, right: %u",left,right));	
			while(left <= right){
				if(!inside){
					mid = (left + right) / 2;
					// delete vertices
					for(ui i = mid; i <= right; i++){
						ui u = cur_verset[i];
						cur_vis[u] = true;	
						for(int j = pstart[u]; j < pend[u]; j++){
							ui v = edgeList[j];
							if(cur_vis[v]) continue;
							ui id = eid[j];
							edge e = orien_edges[id];
							if(e.to == u){
								ui from = e.to == e.end ? e.src : e.end;
								assert(deg[from]>0);
								deg[from]--;
								// printf("update v: %u, deg: %u\n",from,deg[from]);
							}
						}
						deg[u] = 0;
					}
				}	
				// log_debug(graphClock_.Count("left: %u, right: %u, mid: %u",left,right,mid));		
				if(ReTest()){// not exist
					left = mid+1;
					last_mid = (left + right) / 2;
					// restore vertices
					for(ui i = mid; i < last_mid; i++){
						ui u = cur_verset[i];
						cur_vis[u] = false;	
						for(int j = pstart[u]; j < pend[u]; j++){
							ui v = edgeList[j];
							if(cur_vis[v]) continue;
							ui id = eid[j];
							edge e = orien_edges[id];
							if(e.to == u){
								ui from = e.to == e.end ? e.src : e.end;
								assert(deg[from]<=n);
								deg[from]++;
								// printf("update v: %u, deg: %u\n",from,deg[from]);
							}
							else{
								deg[u]++;
							}
						}
					}
					inside = true;
					mid = last_mid;
				}
				else{
					right = mid-1;
					inside = false;
				}
				// unordered_set<ui> isContain_;
				// find_integer_densest_sub(isContain_);
				// ui m_count = 0, cur_min_wei = INT32_MAX;
				// for(auto ii : isContain_){
				// 	cur_min_wei = min(weigh[ori_map[ii]],cur_min_wei);
				// 	for(int j = pstart[ii]; j < pend[ii]; j++){
				// 		ui v = edgeList[j];
				// 		if(v < ii) continue;
				// 		if(isContain_.find(v) != isContain_.end() && !cur_vis[v]) m_count++;
				// 	}
				// }
				// ui u = cur_verset[mid];
				// min_wei = weigh[ori_map[u]];
				
				// if(isContain_.size() != 0)
				// 	log_debug(graphClock_.Count("min_wei: %u, cur_min_wei: %u, isContain.size(): %u, m_c: %u, den: %f, right: %u",weigh[ori_map[u]],cur_min_wei,isContain_.size(),m_count,1.0*m_count/isContain_.size(),right));
			}
			ui u = cur_verset[right];
			min_wei = weigh[ori_map[u]];
			#ifndef NDEBUG
			log_info(graphClock_.Count("pseudoarboricity: %u, maxweight: %u\n",test_value+1,min_wei));
			#endif
			//delete vertex with weight no greater than min_wei from original graph
			for(int i = 0; i < original_n; i++)
				if(!vis[i] && weigh[i] <= min_wei)
					vis[i] = true;
		}
		else{
			log_debug(graphClock_.Count("not find integer densest subgraph in 2-approximation method!"));
		}
	}
    else{
        log_debug(graphClock_.Count("not find integer densest subgraph!"));
    }
}

void Graph::find_maximum_weight_subgraph_2approximate()
{
	isContain.clear();
	ui weight = 0;
	min_wei = INT32_MAX;
	ui current_n = n;
	
	if(!ReTest()){
		find_integer_densest_sub(isContain);
		for(int i = 0; i < n; i++) cur_vis[i] = true;  //search for optimal subgraph
		priority_queue<pair<ui,ui>> pq;
		for(auto u : isContain){
			pq.push({weigh[ori_map[u]],u});	
            // printf("u: %u, w: %u, deg[v]: %u\n",u,weigh[ori_map[u]],deg[u]);
		}
		vector<ui> cur_verset;
		ui cur_count = 0, pre_count = 0;
		bool flag = false;
		memset(deg,0,sizeof(ui)*n);
		while (!pq.empty()) {
        	std::pair<ui, ui> top = pq.top();
			pq.pop();
			cur_verset.push_back(top.second);
			cur_count++;
			if(isPowerOfTwo(cur_count) || pq.empty()){
				construct_subgraph(cur_verset,pre_count);
				if(!ReTest()) 
				{
					flag = true;
					break;
				}
				pre_count = cur_count;
			}
			if(flag)
				break;
		}
		if(flag){
			for(ui i = cur_count-1; i >= pre_count; i--){
				ui u = cur_verset[i];	
                if(cur_vis[u])continue;	
				cur_vis[u] = true;
				for(int j = pstart[u]; j < pend[u]; j++){
					ui v = edgeList[j];
					if(cur_vis[v]) continue;
					ui id = eid[j];
					edge e = orien_edges[id];
					if(e.to == u){
						ui from = e.to == e.end ? e.src : e.end;
						assert(deg[from]>0);
						deg[from]--;
                        // printf("update v: %u, deg: %u\n",from,deg[from]);
					}
				}

				deg[u] = 0;
				if(ReTest()) //here is to be optimized
				{
					cur_vis[u] = false;
					min_wei = weigh[ori_map[u]];
					#ifndef NDEBUG
					log_info(graphClock_.Count("pseudoarboricity: %u, maxweight: %u\n",test_value+1,min_wei));	
					#endif
					break;
				}
				else{
					unordered_set<ui> isContain_;
					find_integer_densest_sub(isContain_);
					ui m_count = 0, cur_min_wei = INT32_MAX;
                    for(auto ii : isContain_){
                        cur_min_wei = min(weigh[ori_map[ii]],cur_min_wei);
                        for(int j = pstart[ii]; j < pend[ii]; j++){
								ui v = edgeList[j];
								if(v < ii) continue;
								if(isContain_.find(v) != isContain_.end() && !cur_vis[v]) m_count++;
							}
                    }
                    for(ui ii = i-1; ii >= pre_count; ii--)
                        if(weigh[ori_map[cur_verset[ii]]] < cur_min_wei)
                            cur_vis[cur_verset[ii]] = true;
					#ifndef NDEBUG
					log_info(graphClock_.Count("min_wei: %u, cur_min_wei: %u, isContain.size(): %u, m_c: %u, den: %f",weigh[ori_map[u]],cur_min_wei,isContain_.size(),m_count,1.0*m_count/isContain_.size()));
					#endif
				}
				
			}
		}
		else{
			log_debug(graphClock_.Count("not find!"));
		}
		//delete vertex with weight no greater than min_wei from original graph
		for(int i = 0; i < original_n; i++)
			if(!vis[i] && weigh[i] <= min_wei)
				vis[i] = true;
	}
    else{
        // for(ui i =0; i < n; i++)
        //     printf("v: %u, deg: %u\n",i,deg[i]);
        log_debug(graphClock_.Count("not find integer densest subgraph!"));
    }
}


void Graph::find_integer_densest_sub(std::unordered_set<ui>& isContain)
{
	queue<ui> Q;
	int pre_num = 0;
	for (unsigned i = 0; i < n; i++)
	{
		if (!cur_vis[i] && deg[i] > test_value){
			Q.push(i);
			isContain.insert(i);
		}
	}
	// printf("prenum: %u\n",isContain.size());
	while(!Q.empty()){
		int x = Q.front();
		Q.pop();
		for (unsigned i = pstart[x]; i < pend[x]; i++)
		{
			ui v = edgeList[i];
			ui id = eid[i];
			edge& ne = orien_edges[id];
			if (ne.to == x || cur_vis[v]) continue;
			ui to = ne.to;
			if(isContain.find(to) != isContain.end()) continue;
			Q.push(to);
			isContain.insert(to);
		}
	}
	// log_info(graphClock_.Count("densest subgraph vertex: %u",isContain.size()));
}
bool Graph::ReTest()
{
	unordered_set<ui> nodes_to_decrease;
	dist = new ui[n]();
	cur = new ui[n]();
	for (ui i = 0; i < n; i++)
	{
		if (!cur_vis[i] && deg[i] > test_value) //vertex with outdegree greater than test_value
		{
			nodes_to_decrease.insert(i);
		}
	}
	while (DinicBFS(nodes_to_decrease))
	{
		memcpy(cur, pstart, n * sizeof(ui));
		vector<ui> to_erase;
		for (auto i : nodes_to_decrease)
		{
			while (deg[i] > test_value)
			{
				p[i] = INF + 1;
				if (!DinicDFS(i))
					break;
			}
			if (deg[i] <= test_value)
				to_erase.push_back(i);
		}
		for (auto i : to_erase)
			nodes_to_decrease.erase(i);
	}

	delete[] dist;
	delete[] cur;
	return nodes_to_decrease.empty();
}

bool Graph::DinicBFS(unordered_set<ui>& nodes_to_decrease)
{
	if (TESTING) TEST::phase_number++;
	if (nodes_to_decrease.empty()) return false;
	bool can_reach_low_indegree = false;
	queue<ui> Q;
	for (unsigned i = 0; i < n; i++)
		dist[i] = INF;
	for (auto i : nodes_to_decrease)
	{
		Q.push(i);
		dist[i] = 0;
	}
	while (!Q.empty())
	{
		unsigned x = Q.front();
		Q.pop();
		for (ui i = pstart[x]; i < pend[x]; i++)
		{
			ui v = edgeList[i];
			if (cur_vis[v]) continue;
			edge& ne = orien_edges[eid[i]];
			if (ne.to == x) continue;
			ui to = ne.to; 
			if (dist[to] != INF) continue;
			dist[to] = dist[x] + 1;
			if (deg[to] < test_value) can_reach_low_indegree = true;
			Q.push(to);
		}
	}
	return can_reach_low_indegree;
}

bool Graph::DinicDFS(ui x)
{
	// if(test_value == 6 && p[x] < INF)
	// 	printf("x: %u, d: %u, from: %u\t",x,deg[x],orien_edges[p[x]].to == orien_edges[p[x]].src ? orien_edges[p[x]].end : orien_edges[p[x]].src);
	if (deg[x] < test_value)
	{
		ui from = orien_edges[p[x]].to == orien_edges[p[x]].src ? orien_edges[p[x]].end : orien_edges[p[x]].src;
		deg[orien_edges[p[x]].to]++;
		assert(deg[from]>0);
		deg[from]--;
		orien_edges[p[x]].to = from;
		return true;
	}
	for (ui& i = cur[x]; i < pend[x]; i++)
	{
		edge& ne = orien_edges[eid[i]];
		if (ne.to == x) continue;
		ui to = ne.to;
		if (dist[to] != dist[x] + 1 || cur_vis[to]) continue;
		p[to] = eid[i];
		if (DinicDFS(to))
		{
			if (p[x] >= INF) return true;
			deg[orien_edges[p[x]].to]++;
			deg[x]--;
			orien_edges[p[x]].to = x;
			return true;
		}
	}
	return false;
}

// get maximum outdegree
unsigned Graph::get_max_d()
{
	unsigned maxd = 0;
	for (unsigned i = 0; i < n; i++)
		maxd = max(maxd, deg[i]);
	return maxd;
}


