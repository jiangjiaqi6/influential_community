#include "Graph.h"

void print_usage() ;
size_t getPeakRSS(Clock &allClock);

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

	
	FILE* dataset = fopen("./dataset.txt", "r");
	
	char graph_name[100];
	int thres = 0;
	while (fscanf(dataset, "%s %d", &graph_name, &thres) == 2)
	{
		Clock allClock("All");
    	log_info(allClock.Start());
		Timer timer;
		Graph *graph = new Graph(graph_name);
		// graph->read_graph_without_vm(graph_name);
		graph->read_graph_binary(graph_name);
		log_info(allClock.Count("thres: %u ,read graph time: %s s",thres,Utility::integer_to_string(timer.elapsed()).c_str()));
		graph->set_tau(thres);
		graph->load_vertex_weight(graph_name);

		//baseline: climb stair method
		// while(!graph->get_finished()){
		// 	graph->extract_subgraph_from_original();
		// 	graph->influential_community();
		// }


		//advance: divide and conquer
		graph->div_con_framework();


		getPeakRSS(allClock);

		delete graph;
		log_info(allClock.Count("Finish\n"));
		// printf("total time: %s s\n",Utility::integer_to_string(timer.elapsed()).c_str());
	}

	return 0;
}

void print_usage() {
	printf("Usage: [1]exe [2]graph-dir [3 optional] \"output\"\n");
}

size_t getPeakRSS(Clock &allClock) {
	std::ifstream file("/proc/self/status");
	if (!file.is_open()) {
		std::cerr << "Error: Unable to open /proc/self/status" << std::endl;
		return 0;
	}
	size_t peakRSS = 0;
	for (std::string line; std::getline(file, line);) {
		if (line.rfind("VmPeak:", 0) == 0) {
			peakRSS = std::stod(line.substr(7));
			log_info(allClock.Count("[MEMORY] VmPeak: %f MB", std::stod(line.substr(7)) / 1024.0));
			// std::cout << " [MEMORY] VmPeak: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
		}
		// if (line.rfind("VmSize:", 0) == 0) {
		// std::cout << " [MEMORY] VmSize: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
		// }
		if (line.rfind("VmHWM:", 0) == 0) {
			// std::cout << " [MEMORY] VmHWM: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
			log_info(allClock.Count("[MEMORY] VmHWM: %f MB", std::stod(line.substr(6)) / 1024.0));
		}
		// if (line.rfind("VmRSS:", 0) == 0) {
		// std::cout << " [MEMORY] VmRSS: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
		// }
	}

	return 0; // VmHWM is reported in kB, convert to bytes
}