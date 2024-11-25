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
	dir = "/mnt/hdd0/jjq/graphData/";  // used for reading binary file
	// string ss(_dir);
	// log_info(graphClock_.Count("graph file name: %s", _dir)); //such as "graph file name: tech-p2p/tech-p2p.edges"
	
	// used for recording output information
	// ss = "result/out-" + ss;
	// fout.open(ss);  

	n = 0, original_n = 0;
	m = 0, original_m = 0;
	max_id = 0;
	maxDeg = 0;
	number = 1;  
	pa = 0;
	max_pa = 0;
	app_pa = 0;

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
	sub_weigh = nullptr;  //the weight in the subgraph
	deg = nullptr;

	eid = nullptr;

	res.clear();  //?
	//terminate condition
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
	if(sub_weigh != nullptr) {
		delete[] sub_weigh;
		sub_weigh = nullptr;
	}
	if(deg != nullptr) {
		delete[] deg;
		deg = nullptr;
	}
	if(orien_edges != nullptr) {
		delete[] orien_edges;
		orien_edges = nullptr;
	}
	if(eid != nullptr) {
		delete[] eid;
		eid = nullptr;
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
		delete[] edgeList;
		edgeList = nullptr;
	}
	if(vis != nullptr) {
		delete[] vis;
		vis = nullptr;
	}
	if(weigh != nullptr) {
		delete[] weigh;
		weigh = nullptr;
	}
}


void Graph::read_graph_binary(const char * graph_name)  // input graph name such as "tech-p2p/tech-p2p.edges"
{
	string ss = dir;   // the common path of the graph data such as "/mnt/hdd0/jjq/graphData/"
	int i = 0, j = 0, len = strlen(graph_name);
	bool flag = true;
	// while(i < len){
	// 	ss += graph_name[i];
	// 	if(graph_name[i] == '/'){
	// 		flag = false;
	// 		break;
	// 	}		
	// 	i++;
	// }
	// if(flag)
	// 	ss = dir;
	
	string path(graph_name);
	size_t pos = path.find_last_of('/');
	string result;
    // extract the str from the last '/'
    if (pos != std::string::npos) {
        result = path.substr(0, pos + 1);
        // std::cout << "correct : " << result << std::endl;
    } else {
        std::cout << "error: not find '/'" << std::endl;
    }
	ss += result;
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

#ifdef DEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	assert(sum == 2*m);
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

	/* delete isolation vertex */

	// unordered_map<ui,ui> tmp_mp;
	// ui vc = 0;
	// for(ui i = 0; i < n;i ++){
	// 	if(degree[i] != 0){
	// 		tmp_mp[i] = vc++;
	// 	}
	// } 

	// log_info(graphClock_.Count("*\tn = %s; m = %s (undirected), graph: %s", Utility::integer_to_string(vc).c_str(), Utility::integer_to_string(m).c_str(),graph_name));

	// for(ui i = 0; i < n; i++){
	// 	for(ui j = original_pstart[i]; j < original_pend[i]; j++){
	// 		original_edgeList[j] = tmp_mp[original_edgeList[j]];
	// 	}
	// 	if(degree[i] == 0) continue;
	// 	original_pstart[tmp_mp[i]] = original_pstart[i];
	// 	original_pend[tmp_mp[i]] = original_pend[i];
	// }

	// original_n = n = vc;


	fclose(f);
	delete[] degree;

}

ui Graph::get_node_id(ui id){
	if( map_id.find(id) == map_id.end() ) {
		map_id[id] = max_id;
		// org_id.push_back( id );
		return max_id++;
	}
	return map_id[id];
}


void Graph::read_graph(const char *input_file) 
{
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
	string ss(input_file);
	while(fgets(line,200,file)){
		if(line[0] == '%' || line[0] == '#')continue;
		ui i = 0;
        ui a = 0;
		while (line[i] < '0' || line[i]>'9') i++;
		while (line[i] >= '0' && line[i] <= '9') a = a * 10 + line[i] - '0', i++;
		ui u = get_node_id(a);
		// u = a;
        ui b = 0;
		while (line[i] < '0' || line[i]>'9') i++;
		while (line[i] >= '0' && line[i] <= '9') b = b * 10 + line[i] - '0', i++;
		ui v = get_node_id(b);
		// v = b;
        if(u == v)continue;
        vp.pb(mp(u,v));
        vp.pb(mp(v,u));
    }
	sort(vp.begin(),vp.end());
    vp.erase(unique(vp.begin(), vp.end()), vp.end());
	n = max_id;
    m = vp.size()/2;

	original_n = n, original_m = m;

	if(original_pstart != nullptr) delete[] original_pstart;
	original_pstart = new ui[n];
	if(original_pend != nullptr) delete[] original_pend;
	original_pend = new ui[n];
	if(original_edgeList != nullptr) delete[] original_edgeList;
	original_edgeList = new ui[m*2];
	if(vis != nullptr) delete[] vis;
	vis = new bool[n]();

	original_pstart[0] = 0;
	ui idx = 0;
	for(ui i = 0;i < n;i ++) {
		original_pend[i] = original_pstart[i];
		while(idx < vp.size()&&vp[idx].first == i) original_edgeList[original_pend[i] ++] = vp[idx ++].second;
		if(i < n-1)
			original_pstart[i+1] = original_pend[i];
		// deg[i] = pend[i]-pstart[i];
		maxDeg = (original_pend[i]-original_pstart[i]) > maxDeg ? (original_pend[i]-original_pstart[i]) : maxDeg;
	}
	

	// idx = 0;
	// for(ui i = 0; i < n; i++){
	// 	for(ui j = pstart[i]; j < pend[i]; j++){
	// 		ui v = edgeList[j];
	// 		if(v < i) continue;
	// 		eid[j] = idx;
	// 		eid[pstart[v]++] = idx;
	// 		orien_edges[idx].src = i;
	// 		orien_edges[idx].end = v;
	// 		orien_edges[idx].to = INT32_MAX;
	// 		idx++;
	// 	}
	// }
	// for(ui i = 0; i < n; i++){
	// 	pstart[i] = pend[i] - deg[i];
	// 	maxDeg = max(maxDeg,deg[i]);
	// 	deg[i] = 0;
	// }
		
	// for(ui i = 0; i < n; i++){
	// 	for(ui j = pstart[i]; j < pend[i]; j++){
	// 		ui v = edgeList[j];
	// 		printf("u: %u, v: %u, eid: %u, orie_u: %u, orie_v: %u, to: %u\n",i,v,eid[j],orien_edges[eid[j]].src,orien_edges[eid[j]].end,orien_edges[eid[j]].to);
	// 	}
	// }
	log_info(graphClock_.Count("*\tn = %s; m = %s (undirected)", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m).c_str()));
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


//generate unique weights for each node with a pagerank method
void Graph::generate_pagerank_unique_weights(string ss){
	if (n > m) {
        std::cerr << "Number of vertices cannot be greater than the number of edges." << std::endl;
        return;
    }
	// PageRank algorithm parameters
	double d = 0.85;      // damping factor
	double epsilon = 1e-6; // convergence threshold
	ui max_iter = 500;    // maximum number of iterations

	vector<double> pagerank(n,1.0/n);
	vector<double> new_pagerank(n,0.0);
	vector<ui> degree(n,0);
	ui zero_deg = 0;
	for(ui i = 0; i < n; i++){
		degree[i] = original_pend[i] - original_pstart[i];
	}

	double max_pr;

	for (ui iter = 0; iter < max_iter; ++iter) {
		max_pr = 0.0;
        // Calculate new PageRank values
        for (ui u = 0; u < n; ++u) {
            double sum = 0.0;
			for(ui j = original_pstart[u]; j < original_pend[u]; j++){
				ui v = original_edgeList[j];
				sum += pagerank[v] / degree[v];
			}
            new_pagerank[u] = (1 - d) / n + d * sum;
        }

        // Check for convergence
        double diff = 0.0;
        for (int u = 0; u < n; ++u) {
            diff += fabs(new_pagerank[u] - pagerank[u]);
            pagerank[u] = new_pagerank[u]; // Update the pagerank value
			max_pr = max(max_pr,pagerank[u]);
        }

        if (diff < epsilon) {
            cout << "Converged after " << iter + 1 << " iterations." << endl;
            break;
        }
    }
	ui double_scale = 1;

	while(max_pr<1){
		double_scale*=10;
		max_pr *= 10;
	}

	// Scale and convert to integer PageRank values
    vector<pair<ui,ui>> node_pagerank(n);
	ui basescale = 1, tmp_m = m;
	while(tmp_m){
		basescale *= 10;
		tmp_m /= 10;
	}
	basescale *= double_scale;
	basescale /= 10;
    double scale_factor = basescale*1.0;  // Scaling factor to amplify differences and convert to integers

	ui max_value_pagerank = 0;
	
    for (ui u = 0; u < n; ++u) {
        ui value_pagerank = static_cast<ui>(pagerank[u] * basescale);
		max_value_pagerank = max(max_value_pagerank,value_pagerank);
		node_pagerank[u]={u, value_pagerank}; 
    }
	max_value_pagerank = max(max_value_pagerank,m)+1;
	ui *count = new ui[m]();
	for (ui u = 0; u < n; ++u) {
		count[node_pagerank[u].second]++;
    }
	
	queue<ui> que;
	for(ui i = 1; i < m; i++){
		if(count[i] == 0)
			que.push(i);
	}

    // Ensure all pagerank values are unique by adjusting any duplicates
    sort(node_pagerank.begin(), node_pagerank.end(), [](const pair<int, int>& a, const pair<int, int>& b) {
        return a.second < b.second;  // Sort by pagerank in descending order
    });
	
    // Ensure all pagerank values are unique by adjusting any duplicates
	ui last;
    for (ui i = 0; i < n; ++i) {
		if(count[node_pagerank[i].second] == 1) continue;
		else{
			// printf("value % u , repeat times: %u!\n",node_pagerank[i].second,count[node_pagerank[i].second]);
			ui j = node_pagerank[i].second;
			i++;
			for(;i < n ; i++){
				if(node_pagerank[i].second != j) {
					i--;break;
				}
				node_pagerank[i].second = que.front();
				que.pop();
				if(que.empty())printf("queue empty error!");
			}
		}
    }

	// for (ui i = 0; i < n; ++i) {
	// 	if(i == 0) {
	// 		last = node_pagerank[i].second;
	// 		printf("i: % u , node_pagerank[i].second: %u!\n",i,node_pagerank[i].second);
	// 		continue;
	// 	}
	// 	if(node_pagerank[i].second != last){
	// 		last = node_pagerank[i].second;
	// 		printf("i: % u , node_pagerank[i].second: %u!\n",i,node_pagerank[i].second);
	// 	}
    // }

	ui *all_weight = new ui[n]();
	for(ui i = 0; i < n; i++){
		all_weight[node_pagerank[i].first] = node_pagerank[i].second;
	}

	ofstream fwrite(ss);
	if (fwrite.is_open()) {
		for(int i = 0; i < n; i++){
			fwrite << all_weight[i];
			fwrite << endl;  
			weigh[i] = all_weight[i];
		}
		fwrite.close();
	} else {
		log_debug(graphClock_.Count("Failed to open file for writing."));
	}
	fwrite.close(); 	

	delete []all_weight;
	delete []count;
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
	string ss = dir;
	int i = 0, len = strlen(graph_name);
	while(i < len && graph_name[i] != '.'){
		ss += graph_name[i++];
	}
	// ss += "-weight1.txt"; //pagerank
	ss += "-weight.txt";
	// std::cout<<ss<<endl;
	ifstream input(ss);
	weigh = new ui[n]();
	
	if(input.good()){  //read prepared file 
		int i = 0;
		// string str_input_file(input_file);
		// ifstream input(str_input_file);
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
			if(i == n) break;
		}
		assert(i == n);
	}
	else{   //generate random number into file
	    generate_unique_weights(ss);
		// generate_pagerank_unique_weights(ss);
	}
}

//influential community based on kmax-core
void Graph::influential_community_based_core(){
	ori_map.clear();  // in greedy_and_reduce function
	ui *peel_sequence = new ui[n]();
	ui *core = new ui[n]();
	ui max_core = core_decomposition(peel_sequence, core, n, pstart, pend, edgeList, eid);
	greedy_and_reduce_kmax_core(peel_sequence, core);
	delete[] peel_sequence;
	delete[] core;
	maxweight_core_computation(max_core);
	if(max_core <= tau+1)
		isFinish = true;
}

void Graph::maxweight_core_computation(ui max_core){
	cur_vis = new bool[n]();
	// priority_queue<pair<ui,ui>> pq;
	priority_queue<pair<ui,ui>, vector<pair<ui,ui>>, std::greater<pair<ui,ui>>> pq;
	for(int i = 0; i < n; i++){
		cur_vis[i] = true;  //search for optimal subgraph
		pq.push({weigh[ori_map[i]],i});	
	} 
	ui kmax_core_n = n, kmax_core_m = m;
	ui del_v = 0, del_m = 0;
	ui tmp_max_wei = 0, tmp_n = 0, tmp_m = 0;
	while(!pq.empty()){
		ui u = pq.top().second;
		ui w = pq.top().first;
		pq.pop();
		// printf("u: %u, weight: %u\n",u,w);
		if(!cur_vis[u]) continue;
		tmp_n = kmax_core_n, tmp_m = kmax_core_m, tmp_max_wei = weigh[ori_map[u]];
		cur_vis[u] = false;
		del_v = 1, del_m = 0;
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(!cur_vis[v]) continue;
			deg[v]--;
			del_m++;
			if(deg[v] < max_core){
				dfs_deletion(v,max_core,del_v,del_m);
			}
		}
		kmax_core_n -= del_v;
		kmax_core_m -= del_m;
	}
	double maxden = tmp_m*1.0/tmp_n;
	out_res.push_back({tmp_max_wei,maxden});
	// log_info(graphClock_.Count("max_wei: %u, kmax_core_n: %u, kmax_core_m: %u, density: %lf",tmp_max_wei,tmp_n,tmp_m,tmp_m*1.0/tmp_n));
	for(int i = 0; i < original_n; i++)
		if(!vis[i] && weigh[i] <= tmp_max_wei)
			vis[i] = true;
	delete []cur_vis;
}

void Graph::dfs_deletion(ui u, ui max_core, ui del_v, ui del_m){
	cur_vis[u] = false;
	del_v++;
	for(ui j = pstart[u]; j < pend[u]; j++){
		ui v = edgeList[j];
		if(!cur_vis[v]) continue;
		deg[v]--;
		del_m++;
		if(deg[v] < max_core){
			dfs_deletion(v,max_core,del_v,del_m);
		}
	}
}

void Graph::greedy_and_reduce_kmax_core(ui *peel_sequence, ui *core){

	ui max_core = core[peel_sequence[n-1]];
	ui n_count = 0, cur_min_wei = INT32_MAX;
	bool *is_in = new bool[n]();
	ui *rid = new ui[n]();
	for(ui i = 0; i < n; i++) rid[i] = n;


	for(ui u = 0; u < n; u++){
		if(core[u] == max_core){
			rid[u] = n_count;
			ori_map[n_count++] = u;
			is_in[u] = true;
		}
	}

	
	ui *tmp_pstart = new ui[n_count]();
	ui *tmp_pend = new ui[n_count]();
	ui *tmp_deg = new ui[n_count]();

	ui m_count = 0;
	for(auto it = ori_map.begin(); it != ori_map.end(); it++){
		ui u = it->second;
		for(int j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(!is_in[v]) continue;
			tmp_deg[it->first]++;
			m_count++;
		}
	}
	ui *tmp_edgeList = new ui[m_count]();
	tmp_pstart[0] = 0;
	for(int i = 1; i < n_count; i++){
		tmp_pend[i-1] = tmp_pstart[i-1]+tmp_deg[i-1];
		tmp_pstart[i] = tmp_pend[i-1];
	}
	tmp_pend[n_count-1] = tmp_pstart[n_count-1]+tmp_deg[n_count-1];

	ui id = 0;
	for(int i = 0; i < n_count; i++){
		ui old_start = pstart[ori_map[i]], old_end = pend[ori_map[i]];
		tmp_pend[i] = tmp_pstart[i];
		for(ui j = old_start;j < old_end;j ++) {
			ui global_v = edgeList[j];
			if(!is_in[global_v]) continue;
			ui v = rid[global_v];
			assert(v<n_count);
			if(v > i) {
				tmp_edgeList[tmp_pend[i]++] = v;
				tmp_edgeList[tmp_pstart[v]++] = i;
			}
		}
	}
	n = n_count, m = m_count/2;
	for(int i = 0; i < n; i++)		tmp_pstart[i] = tmp_pend[i]-tmp_deg[i];

	
	// for(ui i = 0; i < n; i++){
	// 	printf("u: %u, deg[u]: %u\n",i,tmp_deg[i]);
	// 	for(ui j = tmp_pstart[i]; j < tmp_pend[i]; j++){
	// 		printf("%u ",tmp_edgeList[j]);
	// 	}
	// 	printf("\n");
	// }

	delete []is_in;
	delete []rid;

	delete []edgeList;
	delete []pstart;
	delete []pend; 
	delete []deg;
		
	edgeList = tmp_edgeList;
	tmp_edgeList = nullptr;
	pstart = tmp_pstart;
	tmp_pstart = nullptr;
	pend = tmp_pend;
	tmp_pend = nullptr;
	deg = tmp_deg;
	tmp_deg = nullptr;
}


void Graph::write_output_weight_density(char * graph_name, ui para){
	string ss = "";  
	if(para == 0){
		ss += "/home/jjq/research/NetworkFlow/DensestResearch/Influential_Community/exp/core-";
	} 
	else
		ss += "/home/jjq/research/NetworkFlow/DensestResearch/Influential_Community/exp/pa-";
	
	int i = 0, j = 0, len = strlen(graph_name);
	
	while(i < len){
		if(graph_name[i] == '/'){
			break;
		}	
		ss += graph_name[i];	
		i++;
	}
	ss += ".csv";
	// cout<<ss<<endl;
	std::ofstream file(ss);
    if (file.is_open()) {
		if(para == 0)
        	file << "weightCore,denCore\n";
		else
			file << "weightPa,denPa\n";

        for (ui i = 0; i < out_res.size(); ++i) {
            file << out_res[i].first << "," << out_res[i].second << "\n";
        }

        file.close();
		// log_debug(graphClock_.Count("successfully saved in file"));	
    } else {
        std::cerr << "write csv error" << std::endl;
    }
}


void Graph::write_output_prune_effect(char * graph_name, ui para){
	string ss = "";  
	if(para == 0){
		ss += "/home/jjq/research/NetworkFlow/DensestResearch/Influential_Community/exp/cs-";
	} 
	else
		ss += "/home/jjq/research/NetworkFlow/DensestResearch/Influential_Community/exp/dac-";
	
	int i = 0, j = 0, len = strlen(graph_name);
	
	while(i < len){
		if(graph_name[i] == '/'){
			break;
		}	
		ss += graph_name[i];	
		i++;
	}
	ss += ".csv";
	// cout<<ss<<endl;
	std::ofstream file(ss);
    if (file.is_open()) {
        file << "scEdge,graphEdge\n";
        for (ui i = 0; i < out_res.size(); ++i) {
            file << out_res[i].first << "," << out_res[i].second << "\n";
        }
        file.close();
		// log_debug(graphClock_.Count("successfully saved in file"));	
    } else {
        std::cerr << "write csv error" << std::endl;
    }
}


//influential community based on integer densest subgraph
void Graph::influential_community(){
	ori_map.clear();  // in greedy_and_reduce function
	ui *peel_sequence = new ui[n]();
	ui *core = new ui[n]();
	ui max_core = core_decomposition(peel_sequence, core, n, pstart, pend, edgeList, eid);
	// greedy_and_reduce(peel_sequence, core);
	greedy_and_reduce_plus(peel_sequence, core);
	exact_pseudoarboricity_computation();
	delete[] peel_sequence;
	delete[] core;
}

//By default, each vertex has different weights
void Graph::recursive_function(ui left, ui right, ui ver_n, ui edg_m, 
ui *cur_pstart, ui *cur_pend, ui *cur_edgeList, ui *cur_eid)
{
	if(left > right)
		return;
	
	ui mid = (left+right)/2;
	#ifdef DEBUG
	log_info(graphClock_.Count("left: %u, mid: %u, right: %u",left,mid,right));
	#endif	
	ui *peel_sequence = new ui[original_n]();
	ui *core = new ui[original_n]();

	extract_subgraph_from_original();
	ori_map.clear();  // in greedy_and_reduce function
	ui max_core = core_decomposition(peel_sequence, core, n, pstart, pend, edgeList, eid);
	// greedy_and_reduce(peel_sequence,core,mid);
	greedy_and_reduce_plus(peel_sequence, core, mid);
	delete[] peel_sequence;
	delete[] core;
	

	memset(p,0,sizeof(ui)*original_n);
	memset(cur_vis,0,sizeof(bool)*original_n);
	
	assert(mid!=0);
	test_value = mid-1;

	// find_maximum_weight_subgrah_peeling();  //based peeling method
	// find_maximum_weight_subgrah_2approximate();
	find_maximum_weight_subgraph_2approximate_binary();
	bool *backup_vis = new bool[original_n]();
	for(int i = 0; i < original_n; i++) backup_vis[i] = true;
	for(auto it = isContain.begin(); it != isContain.end(); it++)
		backup_vis[ori_map[*it]] = false;
	// assert(min_wei != INT32_MAX);
	recursive_function(left,mid-1,original_n,original_m,original_pstart,original_pend,original_edgeList,eid); //left -- higher weight
	
	for(int i = 0; i < original_n; i++) vis[i] = backup_vis[i];
	recursive_function(mid+1,right,original_n,original_m,original_pstart,original_pend,original_edgeList,eid); //right -- higher density

	delete[] backup_vis;
}

void Graph::div_con_framework(){
	extract_subgraph_from_original();
	ori_map.clear();  // in greedy_and_reduce function
	ui *peel_sequence = new ui[n]();
	ui *core = new ui[n]();
	ui max_core = core_decomposition(peel_sequence, core, n, pstart, pend, edgeList, eid);
	// greedy_and_reduce(peel_sequence, core);
	greedy_and_reduce_plus(peel_sequence, core);
	exact_pseudoarboricity_computation();

	delete[] peel_sequence;
	delete[] core;
	p = new ui[original_n]();
	cur_vis = new bool[original_n]();
	recursive_function(tau,max_pa-1,original_n,original_m,original_pstart,original_pend,original_edgeList,eid);
	delete[] p;
	delete[] cur_vis;
}

void Graph::influential_community_divide_and_conquer(){
	ori_map.clear();  // in greedy_and_reduce function
	ui *peel_sequence = new ui[n]();
	ui *core = new ui[n]();
	ui max_core = core_decomposition(peel_sequence, core, n, pstart, pend, edgeList, eid);
	greedy_and_reduce(peel_sequence, core);
	exact_pseudoarboricity_computation();

	//recursive function
	ui left = tau;
	ui right = pa+1;
	ui mid = (left+right)/2;
	// 选择一个 mid/2 -core 子图S1，并构造成有向图
	// 跑整数值网络流算法得到一个整数值稠密度>= mid 的子图S2
	// 把子图S2顶点按权重排序，采用从大到小遍历插入,找整数值稠密度>= mid 的子图 S3
	// 得到 skyline 上点 (cur_maxweight, mid)
	// skyline 上左边的点构造子图S4，权重大于 cur_maxweight
	// 把子图S4传进 recursive function
	// skyline 上右边的点构造子图S5，基于S2和S3找到整数值稠密度 > mid的子图
	// 把子图S4传进 recursive function
}


// private member functions
// core decomposition:
//		"peel_sequence" stores the vertices in degenerarcy order
//		"core[u]" stores the core value of vertex u
//		"n, pstart, pend, edges" are used to represent the graph
// return the maximum core value
ui Graph::core_decomposition(ui *peel_sequence, ui *core, const ui n, const ui *pstart, const ui *pend, const ui *edgeList, const ui *eid) {
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

	ui max_core = 0, min_core = 0;
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
			/*if use following code, then correspoding to the greedy_and_reduce function*/
			// deg[u]++;
			// orien_edges[eid[j]].to = v;
		}
	}

#ifdef DEBUG
	for(ui i = 0;i < n;i ++) {
		ui cnt = 0;
		ui u = peel_sequence[i];
		for(ui j = pstart[u];j < pend[u];j ++) if(rid[edgeList[j]] > i) ++ cnt;
		assert(cnt == degree[u]);
	}
#endif
	for(ui i = 0; i < n; i++)
	{
		if(core[id[i]]!=0){
			min_core = core[id[i]];
			break;
		}
			
	}

	delete[] degree;
	delete[] degree_start;
	delete[] rid;

	#ifdef DEBUG
	log_info(graphClock_.Count("core decomposition time: %s s, max_core: %u, min_core: %u", Utility::integer_to_string(timer.elapsed()).c_str(),max_core,min_core));
	#endif
	return max_core;
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
		if(degree[i] == 0) printf("error\n");
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
	#ifdef DEBUG
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
			}
			for(ui j = pstart[peel_sequence[i]];j < pend[peel_sequence[i]];j ++) if(rid[edgeList[j]] > i) -- current_m;
		}
		double density = 1.0*best_m/best_n;
		app_pa = ceil(density);
	}
	else{
		app_pa = input_mid;
	}
	
	#ifdef DEBUG
	// log_info(graphClock_.Count("approximate pa: %u",app_pa));
	#endif

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

	for(int i = start_idx; i < n; i++){
		peel_sequence[i-start_idx] = peel_sequence[i];
	} 
	for(int i = 0; i < n; i++) rid[i] = n;
	n -= start_idx;
	// ids.clear(); ids.resize(n);
	vector<ui> ids;
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
		// printf("i: %u, w: %u\n",i,weigh[ori_map[i]]);
	}
	#ifdef DEBUG
	// log_info(graphClock_.Count("reduced subgraph v: %u, m: %u",n,m));
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



void Graph::extract_subgraph_from_original(){
	if(edgeList != nullptr) delete[] edgeList;
	edgeList = new ui[2*original_m];
	if(pstart != nullptr) delete[] pstart;
	pstart = new ui[original_n]();
	if(pend != nullptr) delete[] pend;
	pend = new ui[original_n]();
	if(deg != nullptr) delete[] deg;
	deg = new ui[original_n]();
	
	
	ui e_c = 0, v_c = 0;
	for(ui i = 0; i < original_n; i++){
		pend[i] = pstart[i];
		if(!vis[i]) {
			v_c++;
			for(ui j = original_pstart[i]; j < original_pend[i]; j++){
				assert(original_edgeList[j] < original_n);
				if(!vis[original_edgeList[j]]){
					edgeList[e_c++] = original_edgeList[j];
					pend[i]++;
				}
			}
		}
		if(i < original_n-1)
			pstart[i+1] = pend[i];
		deg[i] = pend[i]-pstart[i];
	}

	n = original_n;
	m = e_c/2;
	update_m = m;

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

	#ifdef DEBUG
	log_info(graphClock_.Count("extracted subgraph ver: %u, edge: %u",v_c,m));
	#endif
}



void Graph::exact_pseudoarboricity_computation(){
	p = new ui[n]();
	cur_vis = new bool[n]();
	if(pa == 0){
		test_value = pa;
		unsigned p_upper = get_max_d() + 1, p_lower = app_pa;
		// log_info(graphClock_.Count("p_upper: %u, p_lower: %u",p_upper,p_lower));
		while (p_upper > p_lower)
		{
			test_value = (p_upper + p_lower) / 2;
			if (ReTest())
				p_upper = test_value;
			else
				p_lower = test_value + 1;
		}
		pa = p_lower;
		max_pa = p_lower;
	}
	// if(pa == 0){
	// 	test_value = app_pa-1;
    //     unsigned p_upper = get_max_d() + 1, p_lower = app_pa;
	// 	log_info(graphClock_.Count("p_upper: %u, p_lower: %u",get_max_d() + 1,p_lower));
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
    //         test_value = m_count/isContain_.size()+1;
	// 	}
    //     log_info(graphClock_.Count("p_upper: %u, p_lower: %u, max_pa: %u",get_max_d() + 1,p_lower,test_value));
	// }
	
	// find subgraph with maximum integer density value
	test_value = pa-1;
	// find_maximum_weight_subgrah_peeling();  //based peeling method
	// find_maximum_weight_subgrah_2approximate();  //based 2-approximate method
	find_maximum_weight_subgraph_2approximate_binary();  //based doubling and binary search method

	pa -= 1;
	if(pa <= tau) isFinish = true;
	delete[] p;
	delete[] cur_vis;
}

void Graph::find_maximum_weight_subgraph_2approximate_binary(){
	isContain.clear();
	ui weight = 0;
	min_wei = INT32_MAX;
	ui current_n = n, current_m = 0, m_count = 0, cur_min_wei = INT32_MAX;
	if(!ReTest()){
		find_integer_densest_sub(isContain);

		// m_count = 0, cur_min_wei = INT32_MAX;
		// for(auto ii : isContain){
		// 	cur_min_wei = min(weigh[ori_map[ii]],cur_min_wei);
		// 	for(int j = pstart[ii]; j < pend[ii]; j++){
		// 		ui v = edgeList[j];
		// 		if(v < ii) continue;
		// 		if(isContain.find(v) != isContain.end() && !cur_vis[v]) m_count++;
		// 	}
		// }
		// min_wei = cur_min_wei;
		// current_m = m_count;

		for(int i = 0; i < n; i++) cur_vis[i] = true;  //search for optimal subgraph
		priority_queue<pair<ui,ui>> pq;
		for(auto u : isContain){
			pq.push({weigh[ori_map[u]],u});	
            // printf("u: %u, w: %u, deg[v]: %u\n",u,weigh[ori_map[u]],deg[u]);
		}
		vector<ui> cur_verset;
		ui cur_count = 0, pre_count = 0;
		bool flag = false;
		double maxden = 0;
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
				// if(!inside){
				// 	unordered_set<ui> isContain_;
				// 	find_integer_densest_sub(isContain_);
				// 	m_count = 0, cur_min_wei = INT32_MAX;
				// 	for(auto ii : isContain_){
				// 		cur_min_wei = min(weigh[ori_map[ii]],cur_min_wei);
				// 		for(int j = pstart[ii]; j < pend[ii]; j++){
				// 			ui v = edgeList[j];
				// 			if(v < ii) continue;
				// 			if(isContain_.find(v) != isContain_.end() && !cur_vis[v]) m_count++;
				// 		}
				// 	}
				// 	min_wei = cur_min_wei;
				// 	maxden = 1.0*m_count/isContain_.size();
				// 	current_m = m_count;
				// 	// if(isContain_.size() != 0)
				// 	// 	log_debug(graphClock_.Count("min_wei: %u, v_c: %u, m_c: %u, den: %f, left: %u, right: %u",cur_min_wei,isContain_.size(),m_count,1.0*m_count/isContain_.size(),left,right));
				// }
			}
			ui u = cur_verset[right];
			min_wei = weigh[ori_map[u]];

			//write output file
			// out_res.push_back({min_wei,maxden});  // record influence and density of the skyline community
			// out_res.push_back({current_m,update_m});  // record influence and density of the skyline community

			#ifdef DEBUG
			log_info(graphClock_.Count("pseudoarboricity: %u, maxweight: %u, density: %lf, current_m: %u, update_m: %u\n",test_value+1,min_wei,maxden,current_m,update_m));
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
        log_debug(graphClock_.Count("pseudoarboricity: %u, threshold: %u , min_wei: %u ,not find integer densest subgraph!",test_value,tau,min_wei));
    }
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
		}
	}

}


void Graph::find_maximum_weight_subgrah_peeling(){
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
	
	printf("pseudoarboricity: %u, maxweight: %u, subgraph.v: %u\n",pa,min_wei,isContain.size());	
}

void Graph::find_maximum_weight_subgrah_2approximate()
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
					}
				}

				deg[u] = 0;
				if(ReTest()) //here is to be optimized
				{
					cur_vis[u] = false;
					min_wei = weigh[ori_map[u]];
					log_info(graphClock_.Count("pseudoarboricity: %u, maxweight: %u\n",test_value+1,min_wei));	
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
					log_info(graphClock_.Count("min_wei: %u, cur_min_wei: %u, isContain.size(): %u, m_c: %u, den: %f",weigh[ori_map[u]],cur_min_wei,isContain_.size(),m_count,1.0*m_count/isContain_.size()));
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
	// printf("afternum: %u\n",isContain.size());
}
bool Graph::ReTest()
{
	// if(test_value > max_pa)
	// 	printf("[error] test_value: %u, max_pa: %u\n",test_value,max_pa);
	unordered_set<ui> nodes_to_decrease;
	dist = new ui[n]();
	cur = new ui[n]();
	for (unsigned i = 0; i < n; i++)
	{
		if (!cur_vis[i] && deg[i] > test_value) //vertex with outdegree greater than test_value
		{
			nodes_to_decrease.insert(i);
		}
	}
	// printf("nodes_to_decrease size: %u\n",nodes_to_decrease.size());
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
	// for(auto u : nodes_to_decrease)
	// 	printf("%u ",u);
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