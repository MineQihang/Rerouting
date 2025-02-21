#include <vector>

struct Service {
    int id, s, t, w;
    std::vector<int> p;
};

struct Route {
    int service_id, w;
    std::vector<int> p;
};

struct Edge {
    int id, u, v;
};

void init(int N, int M, int W, int K,
          std::vector<Edge> E,
          std::vector<Service> D);

std::vector<Route> request(int r);