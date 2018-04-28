#include "graphs.h"

__global__ void parallel_floyd_warshall(int n, int k, int* W);

Flow *pushRelabelLockFreeGPU(Graph *g, int s, int t);