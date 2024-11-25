#pragma once
#ifndef DEFINES_H_
#define DEFINES_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <random>
#include <assert.h> //stardand header file for c.
#include <omp.h>


#include <sys/types.h>
#include <sys/mman.h>  //memory management
#include <sys/stat.h>
#include <sys/time.h>  //gettimeofday function
#include <fcntl.h>  //IO in linux or unix : open and close function
#include <unistd.h>  //IO in linux or unix : read and write function


#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <chrono>
#include <iomanip>
#include <cstdarg>

#include <fcntl.h> 
#include <sys/socket.h>  
#include <netinet/in.h>  
#include <arpa/inet.h>  
#include <functional>
#include <tuple>

#include <cstdarg>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <functional>
#include <future>
#include <stdexcept>
#include <type_traits>  // 包含 std::invoke_result


// #define NDEBUG // must precede cassert to disable assert.

using ui = unsigned int;

#define pb push_back
#define mp make_pair

const unsigned INF = INT32_MAX;

enum GraphStore { uncompressed, byte_compressed, nibble_compressed };
enum GraphOrientation { original_graph, degree_oriented };

#endif /* DEFINES_H_ */