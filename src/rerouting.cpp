#include "rerouting.h"

#include <bits/stdc++.h>

using namespace std;

static const int MAX_EDGE_QUANTITY = 40;
static const int MAX_FREQ_QUANTITY = 32;
static const int MAX_INT = 1e9;

// Parameters
const int ITERATIONS = 30;
const unsigned int SEED = 42;
const int MAX_CONSIDERED_SERVICES = 25;

// Variables
int N, M, W, K;
vector<Edge> edges;
vector<Service> services;
// Current answer state
vector<Route> currentRoutes;
// Frequencies for current paths that were reserved, i's bit of number 1, if
// this frequency is reserved, 0 otherwise
vector<long long> curReservedFrequenciesMask;
// Frequencies for current initials paths that were reserved, i's bit of number
// 1, if this frequency is reserved, 0 otherwise
vector<long long> initReservedFrequenciesMask;
// deleted edges
unordered_set<int> deletedEdges;
// initBelongingService[i][j] = k, if the frequency j at the edge i belongs to the
// service k
int initBelongingService[MAX_EDGE_QUANTITY + 1][MAX_FREQ_QUANTITY + 1];
// Adjacency list
vector<vector<int>> adj;
// Current request number
int currentRequest = 0;
// random engine
default_random_engine randomEngine(SEED);
// Floyd-Warshall algorithm
vector<vector<vector<int>>> dist;
vector<vector<vector<int>>> currentDist;


/**
 * Called once before edge removal requests for initialization
 *
 * @param n quantity of nodes
 * @param m quantity of edges
 * @param w quantity of frequency
 * @param k quantity of services
 * @param initEdges edges
 * @param initServices services
 */
void init(int n, int m, int w, int k, vector<Edge> initEdges,
          vector<Service> initServices) {
    // cerr << "init" << endl;
    N = n;
    M = m;
    W = w;
    K = k;
    currentRequest = 0;
    edges = initEdges;
    services = initServices;
    deletedEdges.clear();
    curReservedFrequenciesMask = vector<long long>(initEdges.size() + 1, 0ll);
    initReservedFrequenciesMask = vector<long long>(initEdges.size() + 1, 0ll);
    currentRoutes.resize(k);
    memset(initBelongingService, 0, sizeof initBelongingService);
    for (const auto &[id, s, t, freq, p] : initServices) {
        currentRoutes[id - 1].service_id = id;
        currentRoutes[id - 1].w = freq;
        currentRoutes[id - 1].p = p;
        for (int edgeId : currentRoutes[id - 1].p) {
            // Reserved the frequency for the edge
            initReservedFrequenciesMask[edgeId] |= (1ll << freq);
            // Reserved the frequency for the edge
            curReservedFrequenciesMask[edgeId] |= (1ll << freq);
            // Indicate freq at edge belongs to the service
            initBelongingService[edgeId][freq] = id;
        }
    }

    if (adj.empty()) {
        adj.resize(N + 1);
        for (const auto &edge : edges) {
            adj[edge.u].push_back(edge.id);
            adj[edge.v].push_back(edge.id);
        }
    }

    dist.resize(W + 1, vector<vector<int>>(N + 1, vector<int>(N + 1, MAX_INT)));
    for(const auto &edge : edges) {
        for(int freq = 1; freq <= W; freq++) {
            dist[freq][edge.u][edge.v] = 1;
            dist[freq][edge.v][edge.u] = 1;
        }
    }
    for(const auto& route : currentRoutes) {
        for(int i = 0; i < route.p.size(); i++) {
            int edgeId = route.p[i];
            Edge &edge = edges[edgeId - 1];
            dist[route.w][edge.u][edge.v] = MAX_INT;
            dist[route.w][edge.v][edge.u] = MAX_INT;
        }
    }

}

/**
 * Find path for service
 *
 * @param service
 * @param reservedFrequenciesMask
 * @param freq
 * @return
 */
vector<int> findPath(const Service &service,
                     const vector<long long> &reservedFrequenciesMask,
                     int freq) {
    int start = service.s;
    int end = service.t;
    vector<int> edgesPath;
    vector<int> bookNode(N + 1, 0);
    vector<int> bookEdge(N + 1, -1);
    queue<int> q;
    q.push(start);
    bookNode[start] = 1;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int edgeId : adj[u]) {
            auto& edge = edges[edgeId - 1];
            int v = (edge.v ^ edge.u) ^ u;
            if (bookNode[v]) continue;
            // Checking if an edge has been removed
            if (deletedEdges.count(edge.id) == 1) continue;
            // Checking whether the frequency at edge is reserved
            if (reservedFrequenciesMask[edge.id] >> freq & 1ll) continue;
            // Checking whether edge was reserved in the initial paths, and if
            // it was reserved, whether it belongs to the current route
            if ((initReservedFrequenciesMask[edge.id] >> freq & 1ll) &&
                initBelongingService[edge.id][freq] != service.id)
                continue;
            bookNode[v] = 1;
            bookEdge[v] = edge.id - 1;
            if (v == end) break;
            q.push(v);
        }
    }
    if (bookEdge[end] == -1) return {};
    int u = end;
    while (u != start) {
        edgesPath.push_back(bookEdge[u] + 1);
        u = (edges[bookEdge[u]].v ^ edges[bookEdge[u]].u) ^ u;
    }
    reverse(edgesPath.begin(), edgesPath.end());
    return edgesPath;
}

/**
 * [Deprecated] 
 * Find multiple paths for a service, prioritizing shorter and diverse paths.
 *
 * @param service The service for which to find paths.
 * @param reservedFrequenciesMask The current frequency reservation mask.
 * @param freq The frequency to use for the paths.
 * @return A vector of paths, sorted by length and diversity.
 */
vector<vector<int>> findMultiPath(const Service &service,
                                  const vector<long long> &reservedFrequenciesMask,
                                  int freq) {
    int start = service.s;
    int end = service.t;
    vector<vector<int>> paths;
    vector<int> bookNode(N + 1, 0);
    vector<vector<int>> bookEdges(N + 1, vector<int>(N+1, -1)); // 记录到达每个节点的边的信息
    queue<pair<int, vector<int>>> q; // 使用pair存储节点和到达该节点的路径
    q.push({start, {}});
    bookNode[start] = 1;
    int shortestPathLength = INT_MAX;
    while (!q.empty()) {
        int u = q.front().first;
        vector<int> currentPath = q.front().second;
        q.pop();

        if (u == end) {
            if (currentPath.size() < shortestPathLength) {
                shortestPathLength = currentPath.size();
                paths.clear();
                paths.push_back(currentPath);
            } else if(currentPath.size() < shortestPathLength + 3) {
                 paths.push_back(currentPath);
            }
            continue;
        }
        
        if (currentPath.size() > shortestPathLength + 2) continue; // 剪枝，超过一定长度不再继续搜索

        for (int edgeId : adj[u]) {
            Edge &edge = edges[edgeId - 1];
            // Checking if an edge has been removed
            if (deletedEdges.count(edge.id) == 1) continue;
            // Checking whether the frequency at edge is reserved
            if (reservedFrequenciesMask[edge.id] >> freq & 1ll) continue;
            // Checking whether edge was reserved in the initial paths, and if
            // it was reserved, whether it belongs to the current route
            if ((initReservedFrequenciesMask[edge.id] >> freq & 1ll) &&
                initBelongingService[edge.id][freq] != service.id)
                continue;
            int v = (edge.v ^ edge.u) ^ u;
            if (bookEdges[u][v] != -1 && bookEdges[u][v] <= currentPath.size() + 1) continue; // 避免重复访问和循环
           
            bookEdges[u][v] = currentPath.size() + 1;
            bookEdges[v][u] = currentPath.size() + 1;
            vector<int> newPath = currentPath;
            newPath.push_back(edgeId);
            q.push({v, newPath});
        }
    }
    
    // Sort paths by diversity (prioritize paths with fewer common edges)
    if (paths.size() > 1) {
        sort(paths.begin() + 1, paths.end(), [&](const vector<int>& a, const vector<int>& b) {
            int commonEdgesAB = 0;
            for (int edgeIdA : a) {
                for (int edgeIdB : paths[0]) { // Compare with the shortest path
                    if (edgeIdA == edgeIdB) {
                        commonEdgesAB++;
                        break;
                    }
                }
            }
            int commonEdgesBB = 0;
            for (int edgeIdA : b) {
                for (int edgeIdB : paths[0]) { // Compare with the shortest path
                    if (edgeIdA == edgeIdB) {
                        commonEdgesBB++;
                        break;
                    }
                }
            }
            return commonEdgesAB < commonEdgesBB;
        });
    }

    return paths;
}


/**
 * Evaluate the impact of a potential path on future flexibility.
 *
 * @param serviceId The ID of the service for which the path is being evaluated.
 * @param path The potential path (sequence of edge IDs).
 * @param freq The frequency being considered for this path.
 * @param reservedFrequenciesMask The current frequency reservation mask.
 * @return A score representing the potential impact of this path on future flexibility.
 */
double evaluatePathImpact(int serviceId, const vector<int>& path, int freq, const vector<long long>& reservedFrequenciesMask, const vector<int>& normalService, vector<int> idleService, int idleIndex) {
    double score1 = 0;
    unordered_set<int> blockedEdges(path.begin(), path.end());
    // Count how many services could potentially be blocked by this reservation in the future
    // Only consider services that haven't been routed yet
    int countIdleService = idleService.size();
    // shuffle(idleService.begin(), idleService.end(), randomEngine);
    // int consideredIdleService = min(MAX_CONSIDERED_SERVICES, countIdleService);
    int consideredIdleService = countIdleService;
    for (int idx = 0; idx < consideredIdleService; ++idx) {
        int otherServiceId = idleService[idx] + 1;
        if (otherServiceId == serviceId) {
            continue;
        }
        if (currentDist[freq][services[otherServiceId - 1].s][services[otherServiceId - 1].t] == MAX_INT) {
            continue;
        }
        // Check if this reservation would block any alternative paths for the other service
        vector<int> otherPath = findPath(services[otherServiceId - 1], reservedFrequenciesMask, freq);
        if (!otherPath.empty()) {
            int overlap = 0;
            // Check if the paths share any edges
            for (int edgeId : otherPath) {
                if (blockedEdges.count(edgeId) > 0) {
                    overlap++;
                }
            }
            score1 -= 1.0 * overlap; // * otherPath.size() / currentDist[freq][services[otherServiceId - 1].s][services[otherServiceId - 1].t];
        }
    }

    // Count how many services could potentially be blocked by this reservation in the future
    double score2 = 0;
    // shuffle(normalService.begin(), normalService.end(), randomEngine);
    int consideredServices = min(MAX_CONSIDERED_SERVICES, (int)normalService.size());
    for (int i = 0; i < consideredServices; ++i) {
        int otherServiceId = normalService[i] + 1;
        if (currentDist[freq][services[otherServiceId - 1].s][services[otherServiceId - 1].t] == MAX_INT) {
            continue;
        }
        // Check if this reservation would block any alternative paths for the other service
        vector<int> otherPath = findPath(services[otherServiceId - 1], reservedFrequenciesMask, freq);
        if (!otherPath.empty()) {
            int overlap = 0;
            // Check if the paths share any edges
            for (int edgeId : otherPath) {
                if (blockedEdges.count(edgeId) > 0) {
                    overlap++;
                }
            }
            score2 -= 1.0 * overlap; // * otherPath.size() / currentDist[freq][services[otherServiceId - 1].s][services[otherServiceId - 1].t];
        }
    }
    int score1Count = consideredIdleService;
    int score2Count = consideredServices;
    return (score1) + (score2);
}

/**
 * Edge deletion request. After the request, it will no longer be possible to
 * use this edge and all edges deleted before. if id turns out to be equal to 0,
 * then you need to return to the original state and ignore all previous
 * deletion requests.
 *
 * @param deletedEdgeId Id of deleted Edge. If id equal 0 need to reset all
 * edges removed before this.
 * @return All laid service paths taking into account remote edges
 */
vector<Route> request(int deletedEdgeId) {
    if (deletedEdgeId == 0) {  // Reset to the original state
        init(N, M, W, K, edges, services);
        return currentRoutes;
    }
    currentRequest++;
    // 1. 释放删除的边上的频率
    // Remember all the deleted edges so as not to use them
    deletedEdges.insert(deletedEdgeId);
    // Clear all failed path by deleted edge deletedEdgeId
    for (Route &route : currentRoutes)
        for (int i = 0; i < route.p.size(); i++)
            if (deletedEdgeId == route.p[i]) {
                // Release the frequency for the edge
                for (int edgeId : route.p) {
                    curReservedFrequenciesMask[edgeId] ^= (1ll << route.w);
                    // update floyd matrix
                    Edge &edge = edges[edgeId - 1];
                    dist[route.w][edge.u][edge.v] = 1;
                    dist[route.w][edge.v][edge.u] = 1;
                }
                route.w = 0;
                route.p.clear();
                break;
            }
    for(int freq = 1; freq <= W; freq++) {
        dist[freq][edges[deletedEdgeId - 1].u][edges[deletedEdgeId - 1].v] = MAX_INT;
        dist[freq][edges[deletedEdgeId - 1].v][edges[deletedEdgeId - 1].u] = MAX_INT;
    }
    // 2. 使用Floyd算法把所有的s-t路径找出来
    // 先计算矩阵
    currentDist = dist;
    for(int k = 1; k <= N; k++) {
        for(int i = 1; i <= N; i++) {
            for(int j = 1; j <= N; j++) {
                for(int freq = 1; freq <= W; freq++) {
                    currentDist[freq][i][j] = min(currentDist[freq][i][j], currentDist[freq][i][k] + currentDist[freq][k][j]);
                }
            }
        }
    }
    // 3. 重新寻找路径
    vector<Route> bestRoutes = currentRoutes;
    vector<long long> bestReservedFrequenciesMask = curReservedFrequenciesMask;
    int minFailedServices = INT_MAX;
    double bestRoutesScore = 0;
    vector<int> normalService;
    vector<int> idleService;
    for(int i = 0; i < K; ++i) {
        if (currentRoutes[i].p.empty()) {
            idleService.push_back(i);
        } else {
            normalService.push_back(i);
        }
    }
    sort(normalService.begin(), normalService.end(), [&](int a, int b) {
        return currentRoutes[a].p.size() > currentRoutes[b].p.size();
    });

    int iterNum = min(ITERATIONS - (currentRequest - 1) * 3, 20);
    for (int iter = 0; iter < iterNum; ++iter) {
        // Prepare the temporary variables
        vector<Route> tempRoutes = currentRoutes;
        vector<long long> tempReservedFrequenciesMask = curReservedFrequenciesMask;
        int failedServices = 0;
        double routesScore = 0;
        if (iter > 0) {
            shuffle(idleService.begin(), idleService.end(), randomEngine);
        }
        // Try to find a path for each service
        int countIdleService = idleService.size();
        for (int idleIndex = 0; idleIndex < countIdleService; ++idleIndex) {
            int serviceId = idleService[idleIndex];
            Route &route = tempRoutes[serviceId];
            Route bestRoute;
            double bestImpactScore = -1e9;
            for (int freq = 1; freq <= W; freq++) {
                if (currentDist[freq][services[route.service_id - 1].s][services[route.service_id - 1].t] == MAX_INT) {
                    continue;
                }
                vector<int> edgesPath = findPath(services[route.service_id - 1], tempReservedFrequenciesMask, freq);
                if (!edgesPath.empty()) {
                    double impactScore = evaluatePathImpact(route.service_id, edgesPath, freq, tempReservedFrequenciesMask, normalService, idleService, idleIndex);
                    if (impactScore > bestImpactScore || (impactScore == bestImpactScore && edgesPath.size() < bestRoute.p.size())) {
                        bestImpactScore = impactScore;
                        bestRoute = {route.service_id, freq, edgesPath};
                    }
                }
            }
            // If no path was found, mark the service as failed
            if (bestRoute.p.empty()) {
                failedServices++;
            } else {
                route = bestRoute;
                routesScore += bestImpactScore;
                // Reserve the frequency for the edge
                for (int edgeId : route.p) {
                    tempReservedFrequenciesMask[edgeId] |= (1ll << route.w);
                }
            }
        }
        // Update the best result
        if (failedServices < minFailedServices || (failedServices == minFailedServices && routesScore > bestRoutesScore)) {
            minFailedServices = failedServices;
            bestRoutes = tempRoutes;
            bestReservedFrequenciesMask = tempReservedFrequenciesMask;
            bestRoutesScore = routesScore;
        }
        // Update idleService
        // vector<int> successService;
        // vector<int> newIdleService;
        // for (const auto& serviceId : idleService) {
        //     if (tempRoutes[serviceId].p.empty()) {
        //         newIdleService.push_back(serviceId);
        //     } else {
        //         successService.push_back(serviceId);
        //     }
        // }
        // cerr << "iter: " << iter << " failedServices: " << failedServices << " successService: " << successService.size() << " newIdleService: " << newIdleService.size() << endl;
        // shuffle(successService.begin(), successService.end(), randomEngine);
        // shuffle(newIdleService.begin(), newIdleService.end(), randomEngine);
        // for (const auto& serviceId : successService) {
        //     newIdleService.push_back(serviceId);
        // }
        // idleService = newIdleService;
    }
    
    // 4. 更新当前的路径
    for(int i = 0; i < K; ++i) {
        if (currentRoutes[i].p.empty() && !bestRoutes[i].p.empty()) {
            for(int edgeId : bestRoutes[i].p) {
                Edge &edge = edges[edgeId - 1];
                dist[bestRoutes[i].w][edge.u][edge.v] = MAX_INT;
                dist[bestRoutes[i].w][edge.v][edge.u] = MAX_INT;
            }
        }
    }
    currentRoutes = bestRoutes;
    curReservedFrequenciesMask = bestReservedFrequenciesMask;

    return currentRoutes;
}