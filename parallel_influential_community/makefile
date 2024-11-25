# ********     Author: Jiaqi Jiang    ******

# -DTP

all: pseudoarboricity

pseudoarboricity: main.cpp Graph.cpp log.cpp 
	g++ -std=c++17 main.cpp Graph.cpp log.cpp -o main -O3 -g -w -lpthread -fopenmp -DTP
		

clean:
	rm -rf *o 