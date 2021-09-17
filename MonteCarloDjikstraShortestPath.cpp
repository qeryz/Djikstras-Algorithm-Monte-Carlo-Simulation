// Djikstra's Shortest Path using a Monte Carlo Simulation
// Author: Marcos Padilla
// Date: 03/21/2021
// The purpose of the Monte Carlo simulation implementation in this assignment is to use a pseudo-random
// number generator to randomly fill a graph with random edges and their cost, up to a specified density,
// all while calculating the average shortest path cost after running multiple Monte Carlo simulations.

#include <iostream>
#include <climits>
#include <vector>
#include <list>
#include <random>
#include <ctime>
#include <string>
#include <chrono>
using namespace std;

const int INFINIT = INT_MAX;
int numIterations = 1000;	// Number of iterations of Monte Carlo simulations, change to your liking

class Graph {
	private:
	vector<int> adjMatrix;	// Adjacency matrix of the graph
	int numVertices;		// Number of vertices in the graph
	int numEdges;			// Number of edges in the graph

	public:
	Graph(int numVertices);		// Parameterized Graph constructor
	int V();						// Returns the number of vertices in the graph
	int E();						// Returns the number of edges in the graph
	bool adjacent (int x, int y);	// Tests whether there is an edge from node x to node y
	list<int> neighbors(int x);			// Lists all nodes y such that there is an edge from x to y
	int get_edge_value(int x, int y);	// Returns the value associated to the edge (x,y)
	void set_edge_value(int x, int y, int value); // Sets the value associated to the edge (x,y) to value
	void toString();						// Prints the graph with edge values
};

  // Initialize the matrix to zero
Graph::Graph(int numVertices) {
	this->numVertices = numVertices;
	this->numEdges = 0;

	// Creation of n x n adjacency matrix, initialized edges to INFINIT
	adjMatrix.resize(numVertices * numVertices, INFINIT);
	for (int i = 0; i < numVertices; ++i){
		adjMatrix[i * numVertices + i] = 0; // For node(i,j) where i = j, set edge value as 0
	}
}

int Graph::V(){
	return numVertices; // Returns the number of vertices in the graph
}

int Graph::E(){
	return numEdges; // Returns the number of edges in the graph
}

  // Return true if x and y are neighbors, false if not
bool Graph::adjacent (int x, int y){
  	return (get_edge_value(x, y) != INFINIT && get_edge_value(x, y) != 0);
}

list<int> Graph::neighbors(int x){
	list<int> n; // list of neighbors y of x, if any
	for (int y = 0; y < numVertices; y++){
		if (adjacent(x,y) == true){
			n.push_back(y);
		}
	}
	return n;
}

int Graph::get_edge_value(int x, int y){
	return (adjMatrix[x * numVertices + y]); // Return edge value
}

void Graph::set_edge_value(int x, int y, int value){
	if (adjacent(x,y) == false){ // If there is not an existing edge, increment number of edges
  		++numEdges;
	}
	adjMatrix[x * numVertices + y] = adjMatrix[y * numVertices + x] = value; // add new edge value
}

  // Print and format the matrix
void Graph::toString(){
	cout << "Graph: " << endl;
	cout << "    ";
	for (int j = 0; j < numVertices; j++){
		cout << j << " ";
	}
	cout << endl;
	for (int i = 0; i < numVertices; i++) {
    	cout << i << " : ";
    	for (int j = 0; j < 10; j++){
    		if (adjacent(i, j) == true || adjMatrix[i * numVertices + j] == 0){
        		cout << adjMatrix[i * numVertices + j] << " ";
    		}
    		else
    			cout << "-" << " ";
    	}
    	for (int j = 10; j < numVertices; ++j){
    		if (adjacent(i, j) == true || adjMatrix[i * numVertices + j] == 0){
        		cout << " " << adjMatrix[i * numVertices + j] << " ";
    		}
    		else
    			cout << " " << "-" << " ";
    	}
      cout << "\n";
    }
}

class MonteCarlo{
	public:
	MonteCarlo();
	Graph simulation(double density, int minDistRange, int maxDistRange, int vertices);

};

MonteCarlo::MonteCarlo(){

}

// Graph MonteCarlo::simulation ; generates a Monte Carlo simulated graph for a set density
// The following, randEdge, randomly generates an edge distance value between 1 - 10
// The following, randNum, randomly generates a number in the range 1-10. Such that
// any number generated that is <= density * 10 (for instance, density = 0.2 * 10 = 2), will
// set the edge value, randEdge, for node(i,j).
Graph MonteCarlo::simulation(double density, int minDistRange, int maxDistRange, int vertices){
	int maxEdges = ( (vertices * (vertices - 1))/2 ) * density; //Max edges = |V|(|V|- 1)/2 * density

	uniform_int_distribution<unsigned> randEdge(minDistRange, maxDistRange - 1);
	uniform_int_distribution<unsigned> randNum(1, 10);
	unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
	default_random_engine e(seed);

	Graph g(vertices);
	// cout << "NEW GRAPH G CREATED WITH EDGES OF G = " << g.E() << endl;

	while(g.E() != maxEdges){ // Loop until number of edges reaches desired maxEdges needed for desired density
		for (int i = 0; i < g.V(); ++i){
			for (int j = i + 1; j < g.V(); ++j){

				// if randNum probability is less than or equal to density, and if edges dont surpass maxEdges,
				// and if an edge value node(i,j) has not been already set, set an edge value for node (i, j)
				if (randNum(e) <= density * 10 && g.E() < maxEdges && g.get_edge_value(i, j) == INFINIT)
					g.set_edge_value(i, j, randEdge(e));
			}
		}
	}
	return g;
}

class Djikstra{
	public:
	int shortestPath(Graph g, int src); // Finds shortest path from source to rest of vertices in graph g
	int minDistance(Graph g, int dist[], bool pathSet[]); // Helper function for shortestPath function, returns
	// index of the vertex with the shortest distance to be included in the shortest path tree set, pathSet[]
};

int Djikstra::minDistance(Graph g, int dist[], bool pathSet[]){
    // Initialize min value
    int min = INFINIT;
    int min_index;

    for (int v = 0; v < g.V(); v++) {
        if (pathSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;
    }
    return min_index;
}

int Djikstra::shortestPath(Graph g, int src){

	static int count = 0; 	// Counter for total running simulations
	string srcToString = to_string(src); // Converts source int to string, explained later.
	string path[g.V()];		// String array containing the shortest path from source to node i
	int dist[g.V()]; 		// The output array. dist[i] will hold the shortest distance from source to node i
	bool pathSet[g.V()]; // pathSet[i] will be true if node i has been processed in shortest path

    // Initialize all distances as INFINIT, pathSet as false
    for (int i = 0; i < g.V(); i++){
        dist[i] = INFINIT;
        pathSet[i] = false;
        path[i] = srcToString; // first node in shortest path is always the source
    }

    // Distance of source vertex from itself is 0
    dist[src] = 0;

    // Find the shortest path for all the vertices
    for (int vert = 0; vert < g.V() - 1; vert++) {

        // Choose the vertex with the minimum distance from the set of vertices not yet processed.
        int u = minDistance(g, dist, pathSet);

        // Update the chosen vertex as processed
        pathSet[u] = true;

        // Then update the edge value of the adjacent nodes of the chosen vertex
        for (int v = 0; v < g.V(); v++) {
        	int edge = g.get_edge_value(u,v); // Get edge value of node(u,v)

            // Update dist[v] if v is not in pathSet, and if there exists an edge from u to v, and if the total
        	// weight of path from source to v through u is less than the current value of dist[v],
        	// update dist[v] and update the string array, path[v], which displays the current shortest path

            if (!pathSet[v] && edge != INFINIT && dist[u] != INFINIT && dist[u] + edge < dist[v]){
            	path[v] = path[u] + " -> " + to_string(v);
                dist[v] = dist[u] + edge;
            }
        }
    }

    // Print the graph followed by the distance array
    int maxEdges = ((g.V() * (g.V() - 1)) / 2) ; // Max number of edges
    double density = static_cast<double>(g.E())/(maxEdges)*100;	// Calculate real density reached
    int reachableVert = 0, sumPath = 0, avgShortestPath = 0;

    cout << "----------------------------------------------------" << endl;
    cout << "DJIKSTRA'S SHORTEST PATH: MONTECARLO SIMULATION #" << ++count << endl;
    cout << "----------------------------------------------------" << endl;
    cout << "Density: " << density << "%" << endl;
    cout << "Number of edges: " << g.E() << endl;
    cout << "Number of vertices: " << g.V() << endl << endl;


    // Print the current vertex, the shortest path distances, followed by the shortest path itself
    //printf("\nVertex \t Shortest Distance from %d \t Shortest Path \n", src);
    for (int i = 0; i < g.V(); i++){
    	if (dist[i] != INFINIT){
    		printf("Shortest Path (%d to %d) = %d", src, i, dist[i]);
    		cout << ", Path: " << path[i] << endl;
    		reachableVert++; 			// increment vertices that are reachable
    		sumPath += dist[i]; 		// sum shortest distances for each reachable node
    	}
    	else{
    		printf("Shortest Path (%d, %d) = UNREACHABLE", src, i);
			cout << " " << "-" << endl;
    	}
    }
    avgShortestPath = sumPath / reachableVert; // Calculates average shortest path
    cout << "\nAverage Shortest Path for current simulation #" << count << 
		" of density " << density/100 << " is: " << avgShortestPath << endl;
    cout << endl;

	// Check counter if it is to be reset after total number of simulations for the current density
	// signaling a new set of simulations for the next density.
	if ( count == numIterations)
		count = 0;

	// Return average shortest path for current simulation	
	return avgShortestPath;
    // Print graph g
    // g.toString();

}


// Main driver function
int main(){

	int vertices = 50;			// Vertices is set to 50, change to your liking
	int src = 0; 				// Source vertex, change to your liking
	double avgeShortest2 = 0;	// Average shortest path after multiple Monte Carlo simulations of density 0.2
	double avgeShortest4 = 0;	// Average shortest path after multiple Monte Carlo simulations of density 0.4	
	Graph g(vertices);			// Initialize Graph g
	MonteCarlo montecarlo; 		// Constructor for MonteCarlo class
	Djikstra djikstras;			// Constructor for Djikstra class

	// Simulations of Density: 20%
	for (int i = 0; i < numIterations; ++i){
		g = montecarlo.simulation(0.2, 1, 10, vertices);	// Runs a MonteCarlo simulated graph for given density
		avgeShortest2 += djikstras.shortestPath(g, src);	// Returns shortest path using Djikstra's Algorithm
	}
	// End of simulations

	// Simulations of Density: 40%
	for (int i = 0; i < numIterations; ++i){
		g = montecarlo.simulation(0.4, 1, 10, vertices);	// Runs a MonteCarlo simulated graph for given density
		avgeShortest4 += djikstras.shortestPath(g, src);	// Returns shortest path using Djikstra's Algorithm
	}
	// End of simulations

	avgeShortest2 = avgeShortest2 / numIterations;
	cout << "The total average shortest path for " << numIterations << " Monte Carlo simulations of graphs of density 0.2 is: "
		<< avgeShortest2 << endl;

	avgeShortest4 = avgeShortest4 / numIterations;
	cout << "The total average shortest path for " << numIterations << " Monte Carlo simulations of graphs of density 0.4 is: "
		<< avgeShortest4 << endl;

	return 0; // End of program.
}
