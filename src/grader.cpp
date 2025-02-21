#include <vector>
#include <iostream>
#include <cmath>
#include <ctime>
#include "rerouting.h"

using namespace std;

int main() {
    auto start = clock();
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    int n, m, w, k;
    cin >> n >> m >> w >> k;
    vector<Edge> edges(m);
    int id, u, v;
    int x, y;
    for (int i = 0; i < n; ++i) {
        cin >> id >> x >> y;
    }
    for (int i = 0; i < m; ++i) {
        cin >> id >> u >> v;
        edges[id - 1] = {id, u, v};
    }
    vector<Service> services(k);
    int pLen;
    for (int i = 0; i < k; ++i) {
        cin >> id;
        --id;
        services[id].id = id + 1;
        cin >> services[id].s >> services[id].t >> services[id].w >> pLen;
        services[id].p.resize(pLen);
        for (int j = 0; j < pLen; ++j) {
            cin >> services[id].p[j];
        }
    }
    init(n, m, w, k, edges, services);
    int q;
    cin >> q;
    int r;
    
    double total_score = 3 * 1e7;
    int request_count = 0;
    int subsequence_count = 0;
    
    while (q--) {
        cin >> r;
        
        if (r == 0) {
            subsequence_count++;
            request_count = 0;
        } else {
            request_count++;
        }
        
        vector<Route> answer = request(r);
        answer.resize(k);
        int unsuccessful_services = 0;
        for (const auto& route : answer) {
            cout << route.service_id << ' ' << route.w << ' ' << route.p.size() << ' ';
            for (int edge_id : route.p) {
                cout << edge_id << ' ';
            }
            cout << '\n';
            if (route.p.empty()) {
                unsuccessful_services++;
            }
        }
        cout.flush();

        if (r != 0) {
            double C_q = floor(pow(4, 10 - request_count) * 100 * unsuccessful_services / k);
            total_score -= C_q;
            // cerr << "Request " << r << " score: " << (int)C_q << endl;
        } else {
        //   cerr << "Request " << r << " score: " << 0 << endl;
        }
    }
    cerr << "Total score: " << (int)max(1.0, total_score) << endl;
    auto end = clock();
    cerr << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
}