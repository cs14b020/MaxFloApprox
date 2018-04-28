#include <unistd.h>
#include "kernels.h"
#include "CycleTimer.h"
#include <cstdlib>
#include <cuda.h>
#include <bits/stdc++.h>
#include <tuple>
#define IDX(i, j, n) ((i) * (n) + (j))

using namespace std;

typedef vector< tuple<int, int, int> > tuple_list;

struct Graph1 {
    int n;
    int m;
    int *capacities;
    tuple_list edges;

};

struct New_graph {
    Graph* graph;
    int k;
    int* removed_set;
};

struct Var_helper {
    int size;
    float sum;
    float sqsum;
};

Graph1* scan() {
    ifstream myfile("input.txt");
    int n, m, a, b, c;
    myfile >> n >> m;
    Graph1* graph = (Graph1 *)malloc(sizeof(Graph1));
    graph->n = n;
    graph->m = m;
    tuple_list edges;
    int *capacities = (int *)calloc((n * n), sizeof(int));
    while(m--){
        myfile >> a >> b >> c;
        capacities[IDX(a,b,n)] = c;
        edges.push_back( tuple<int, int, int>(c, a, b));
    }
    myfile.close();
    graph->capacities = capacities;
    graph->edges = edges;
    return graph;
}

Graph1* scan_format(int no, int id) {
    string filename = "sample_graphs/g_" + to_string(no) + "_" + to_string(id) + ".txt";
    // filename = "small.txt";
    ifstream myfile(filename);
    int i,j, n, m, a, b, c, edge_count=0;
    string s1,s2,s3;
    string line;
    for(i=0;i<7;i++){
        getline(myfile, line);
    }
    myfile>> s1 >> s2 >> n >> m;

    Graph1* graph = (Graph1 *)malloc(sizeof(Graph1));
    graph->n = n;
    tuple_list edges;
    int *capacities = (int *)calloc((n * n), sizeof(int));
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            capacities[IDX(i,j,n)] = 0;
    while(m--){
        myfile >> s1 >> a >> b >> c;
        a--;
        b--;
        c++;
        if(capacities[IDX(a,b,n)] == 0){
            capacities[IDX(a,b,n)] = c;
            edges.push_back( tuple<int, int, int>(c, a, b));
            edge_count++;
        }
    }
    myfile.close();
    graph->m = edge_count;
    graph->capacities = capacities;
    graph->edges = edges;
    return graph;
}

int BFS(Graph *g, int *flowMatrix, int *parents, int *pathCapacities, int s, int t) {
    memset(parents, -1, (g->n * sizeof(int)));
    memset(pathCapacities, 0, (g->n * sizeof(int)));
    parents[s] = s;
    pathCapacities[s] = std::numeric_limits<int>::max();
    std::queue<int> bfsQueue;
    bfsQueue.push(s);
    while (!bfsQueue.empty()) {
        int u = bfsQueue.front();
        bfsQueue.pop();
        for (int v = 0; v < g->n; v++) {
            if (u == v) continue;
            int residual = g->capacities[IDX(u, v, g->n)] - flowMatrix[IDX(u, v, g->n)];
            if ((residual > 0) && (parents[v] == -1)) {
                parents[v] = u;
                pathCapacities[v] = std::min(pathCapacities[u], residual);
                if (v != t) {
                    bfsQueue.push(v);
                } else {
                    int result = pathCapacities[t];
                    return result;
                }
            }
        }
    }
    return 0;
}

// Edmonds-Karp algorithm to find max s-t flow
Flow *edKarpSeq(Graph *g, int s, int t) {
    int flow = 0;
    int *flowMatrix = (int *)calloc((g->n * g->n), sizeof(int));
    int *parents = (int *)malloc(g->n * sizeof(int));
    int *pathCapacities = (int *)calloc(g->n, sizeof(int));
    while (true) {
        int tempCapacity = BFS(g, flowMatrix, parents, pathCapacities, s, t);
        if (tempCapacity == 0) {
            break;
        }
        flow += tempCapacity;
        int v = t;
        // backtrack
        while (v != s) {
            int u = parents[v];
            flowMatrix[IDX(u, v, g->n)] += tempCapacity;
            flowMatrix[IDX(v, u, g->n)] -= tempCapacity;
            v = u;
        }
    }
    Flow *result = (Flow *)malloc(sizeof(Flow));
    result->maxFlow = flow;
    result->finalEdgeFlows = flowMatrix;
    free(parents);
    free(pathCapacities);
    return result;
}

set<int> into_node(Graph* g, int v) {
    int i, n = g->n;
    set<int> ret;
    for(i=0;i<n;i++) {
        if((g->capacities)[IDX(i,v,n)] != 0){
            ret.insert(i);
        }
    }
    return ret;
}

set<int> out_node(Graph* g, int v) {
    int i, n = g->n;
    set<int> ret;
    for(i=0;i<n;i++) {
        if((g->capacities)[IDX(v,i,n)] != 0){
            ret.insert(i);
        }
    }
    return ret;
}

set<int> into_node_intersect(Graph* g, int v, set<int> removal_set) {
    int i, n = g->n;
    set<int> ret;
    for(i=0;i<n;i++) {
        if(removal_set.find(i) != removal_set.end()){
            if((g->capacities)[IDX(v,i,n)] != 0){
                ret.insert(i);
            }
        }
    }
    return ret;
}

set<int> out_node_intersect(Graph* g, int v, set<int> removal_set) {
    int i, n = g->n;
    set<int> ret;
    for(i=0;i<n;i++) {
        if(removal_set.find(i) != removal_set.end()){
            if((g->capacities)[IDX(i,v,n)] != 0){
                ret.insert(i);
            }
        }
    }
    return ret;
}
set<int> into_nodes_set(Graph* g, set<int> removal_set) {
    set<int> ret;
    for(auto x : removal_set) {
        set<int> temp = into_node(g, x);
        ret.insert(temp.begin(), temp.end());
    }
    for(auto x : removal_set) {
        if(ret.find(x) != ret.end())
            ret.erase(x);
    }
    return ret;
}

set<int> out_nodes_set(Graph* g, set<int> removal_set) {
    set<int> ret;
    for(auto x : removal_set) {
        set<int> temp = out_node(g, x);
        ret.insert(temp.begin(), temp.end());
    }
    for(auto x : removal_set) {
        if(ret.find(x) != ret.end())
            ret.erase(x);
    }
    return ret;
}

float into_node_intersect_sum(Graph* g, int v, set<int> removal_set){
    int i, n = g->n;
    float ret=0;
    for(i=0;i<n;i++) {
        if(removal_set.find(i) != removal_set.end()){
            int edge = (g->capacities)[IDX(v,i,n)];
            if(edge != 0){
                ret += edge;
            }
        }
    }
    return ret;
}

float out_node_intersect_sum(Graph* g, int v, set<int> removal_set){
    int i, n = g->n;
    float ret=0;
    for(i=0;i<n;i++) {
        if(removal_set.find(i) != removal_set.end()){
            int edge = (g->capacities)[IDX(v,i,n)];
            if(edge != 0){
                ret += edge;
            }
        }
    }
    return ret;
}

float minimum_edge(Graph* g, set<int> nodes){
    int n = g->n;
    float ret = INT_MAX;
    for(auto x:nodes){
        for(auto y:nodes){
            if(x != y){
                float temp = g->capacities[IDX(x, y, n)];
                if(temp != 0)
                    ret = min(temp, ret);
            }   
        }
    }
    return ret;
}

void del_part_graph(Graph* g, set<int> removal_set, int* connectivity) {
    int i,n = g->n;
    std::vector<tuple<int, int, int> > edges;
    set<int> into_nodes = into_nodes_set(g, removal_set) ;
    set<int> out_nodes = out_nodes_set(g, removal_set);
    float minedge = minimum_edge(g, removal_set);
    int k = removal_set.size();
   for(auto x: into_nodes){
        for(auto y: out_nodes){
            set<int> inner_into = into_node_intersect(g, x, removal_set);
            set<int> inner_out = out_node_intersect(g, y, removal_set);
            float replace = minedge;
            int flag = 0;
            for(auto a: inner_into){
                for(auto b: inner_out){
                    if((a == b) || connectivity[IDX(a,b,n)]){
                        (g->capacities)[IDX(x,y,n)] = replace;
                        flag = 1;
                        break;
                    }
                }
                if(flag) {
                    break;
                }
            }
        }
    }

    // Disconnecting the removaal set from the rest graph
    for(auto x:removal_set){
        for(i=0;i<n;i++){
            (g->capacities)[IDX(i,x,n)] = 0;
            (g->capacities)[IDX(x,i,n)] = 0;
        }
    }
}

New_graph* findRemovedSet(Graph1* g, tuple_list edges, set<int> imp, int* connectivity, float MAXVAR) {
    int i,m,n,a,b;
    m = edges.size();
    n = g->n;;
    float c;
    bool debugFlag = false;
    sort(edges.begin(), edges.end());
    if(debugFlag){
        for(i=0;i<m;i++)
            cout<<get<0>(edges[i])<<" "<<get<1>(edges[i])<<" "<<get<2>(edges[i])<<endl;
    }

    Var_helper* var_helper = (Var_helper*)malloc(sizeof(Var_helper));
    var_helper->size = 0;
    var_helper->sum = 0;
    var_helper->sqsum = 0;
    set<int> removal_set;
    set<int> removed_set;
    Graph* g1 = (Graph *)malloc(sizeof(Graph));
    g1->n = g->n;
    g1->capacities = (int *)malloc(n*n*sizeof(int));
    memcpy(g1->capacities, g->capacities, n*n*sizeof(int));
    for(i=0;i<m;i++){
        a = get<1>(edges[i]);
        b = get<2>(edges[i]);
        c = get<0>(edges[i]);
        if((removed_set.find(a) == removed_set.end()) && (removed_set.find(b) == removed_set.end())){
            if((imp.find(a) == imp.end()) && (imp.find(b) == imp.end())) {
                float sum = var_helper->sum;
                float sqsum = var_helper->sqsum;
                int k = var_helper->size;
                k++;
                sum += c;
                sqsum += (c*c);
                float var = (sqsum/k) - (sum/k)*(sum/k);
                if(debugFlag){
                    cout<<i<<" :"<<m<<": "<<c<<" --- "<<var<<" ----> "<<k<<endl;
                }
                if(var < MAXVAR){
                    var_helper->sum = sum;
                    var_helper->sqsum = sqsum;
                    var_helper->size = k;
                    removal_set.insert(a);
                    removal_set.insert(b);
                    if(debugFlag){
                        cout<<"inserted "<<a<<" "<<b<<endl;
                    }
                   
                }else {
                    if(removal_set.size() > 2) { 
                        removed_set.insert(removal_set.begin(), removal_set.end());
                        if(debugFlag){
                            for(auto x:removal_set){
                                cout<< x ;
                                cout<< " ";
                            }
                            cout<<endl;
                        }
                        del_part_graph(g1, removal_set, connectivity);
                    }
                    removal_set.clear();
                    i--;
                    var_helper->sum = 0;
                    var_helper->sqsum = 0;
                    var_helper->size = 0;
                }
                
            }
        }
    }
    removed_set.insert(removal_set.begin(), removal_set.end());
    if(debugFlag){
        for(auto x:removal_set){
            cout<< x ;
            cout<< " ";
        }
        cout<<endl;
    }
    del_part_graph(g1, removal_set, connectivity);
    removal_set.clear();
    int k = removed_set.size();
    New_graph* ret = (New_graph*)malloc(sizeof(New_graph));
    ret->graph = g1;
    ret->removed_set = (int*) malloc(k*sizeof(int));
    ret->k = k;
    i=0;
    for(auto x:removed_set){
        (ret->removed_set)[i++] = x;
    }
    return ret;
}

Graph* converted_graph(Graph* g, vector<int> consi_vert){
    Graph* ret = (Graph*)malloc(sizeof(Graph));
    int k = consi_vert.size();
    int n = g->n;
    ret->n = k;
    ret->capacities = (int*)malloc(k*k*sizeof(int));
    int i = 0,j=0;
    for(i=0;i<k;i++){
        for(j=0;j<k;j++){
            (ret->capacities)[IDX(i,j,k)] = (g->capacities)[IDX(consi_vert[i],consi_vert[j],n)];
        }   
    }
    return ret;
}

int main(int argc, char** argv) {
    if(argc != 4){
        cout<< "Command Format"<<endl;
        cout<<"./main N id MAXVAR"<<endl;
        cout<< "N -> No. of vertices(100, 1000)"<<endl;
        cout<< "id -> Index among graphs(0 to 49)"<<endl;
        return;
    }
    srand (time(NULL));
    int no = atoi(argv[1]);
    int id = atoi(argv[2]);
    float MAXVAR = stof(argv[3]);
	Graph1 *graph = scan_format(no,id);
	int i,n = graph->n;
    double start, finalTime, origTime, changeTime, preTime;
    int origFlow,changeFlow;
    
    //********** Get connectivity of vertices here****************
	int *W, *con;
	con = (int*)malloc(n*n*sizeof(int));
	cudaMalloc(&W, n*n*sizeof(int));
	cudaMemcpy(W, graph->capacities, n*n*sizeof(int), cudaMemcpyHostToDevice);
    start = CycleTimer::currentSeconds();
    for(int k=0; k<n; k++) {
        parallel_floyd_warshall<<<1, n*n>>>(n, k, W);
        cudaThreadSynchronize();
    }
    finalTime = CycleTimer::currentSeconds() - start;
    // printf("connectivity took %f\n", finalTime);
    cudaMemcpy(con, W, n*n*sizeof(int), cudaMemcpyDeviceToHost);
    // *********** Connectivity part ends here *********************



    int s = rand()%n;
    int t= s;
    while(t == s){
        t = rand()%n;
    }
    set<int> imp;
    imp.insert(s);
    imp.insert(t);

   
    tuple_list edges = graph->edges;
    start = CycleTimer::currentSeconds();
    New_graph* ret = findRemovedSet(graph, edges, imp, con, MAXVAR);
    preTime = CycleTimer::currentSeconds() - start;
    int k = ret->k;
    int* removed_set = ret->removed_set;
    Graph* g1 = ret->graph;

    set<int> consi_vert;
    for(i=0;i<n;i++)
        consi_vert.insert(i);

    for(i=0;i<k;i++)
        consi_vert.erase(removed_set[i]);

    vector<int> v(consi_vert.begin(), consi_vert.end());
    sort(v.begin(),v.end());
    Graph* conv_graph = converted_graph(g1, v);

	Graph *g = (Graph *)malloc(sizeof(Graph));
    g->n = graph->n;
    g->capacities = graph->capacities;
    Flow *result;


    start = CycleTimer::currentSeconds();
    result = edKarpSeq(g, s, t);
    origTime = CycleTimer::currentSeconds() - start;
    origFlow = result->maxFlow;
    // ***** Application Ends here **********************************
    int news = distance(consi_vert.begin(), consi_vert.find(s));
    int newt = distance(consi_vert.begin(), consi_vert.find(t));
    start = CycleTimer::currentSeconds();
    result = edKarpSeq(conv_graph, news, newt);
    changeTime = CycleTimer::currentSeconds() - start;
    changeFlow = result->maxFlow;
    printf("%d %d %f %f %d %f %f %d \n", no,id, MAXVAR, origTime, origFlow, preTime+finalTime, changeTime, changeFlow); 
    consi_vert.clear();
    
    return 0;
}