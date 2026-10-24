// Exact contracted helper-sheet connectivity. No assertions or sampling.
// Mode bfs: recursive core subsets, then adjacency-list breadth-first search.
// Mode dsu: numeric Gosper core masks, then disjoint-set union.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
using namespace std;
constexpr int Q = 2048;
struct Record { int source, label, cost; };
struct Core { uint32_t mask = 0; int weight = 999; };
struct Pair { int a, b, charge, weight; };
int add(int a, int b) { return ((((a >> 4) + (b >> 4)) & 127) << 4) | ((a ^ b) & 15); }
int neg(int a) { return (((-(a >> 4)) & 127) << 4) | (a & 15); }
void need(bool condition, const string &message) { if (!condition) throw runtime_error(message); }
uint64_t mix(uint64_t h, uint64_t x) { return (h ^ x) * UINT64_C(1099511628211); }

int main(int argc, char **argv) {
    need(argc == 3, "usage: sheets INPUT bfs|dsu");
    const string mode = argv[2];
    need(mode == "bfs" || mode == "dsu", "unknown algorithm");
    ifstream in(argv[1]);
    int nc, nh;
    in >> nc >> nh;
    need(nc == 28 && nh == 51, "exact inherited inventory sizes");
    vector<Record> core(nc), helper(nh);
    for (auto &r : core) in >> r.source >> r.label >> r.cost;
    for (auto &r : helper) in >> r.source >> r.label >> r.cost;
    need(bool(in), "complete input");
    uint32_t weight3 = 0, weight4 = 0;
    for (int i = 0; i < nc; ++i) {
        need(core[i].cost >= 2 && core[i].cost <= 4, "core literal clock costs");
        if (core[i].cost == 3) weight3 |= 1u << i;
        if (core[i].cost == 4) weight4 |= 1u << i;
    }
    auto weight = [&](uint32_t mask) {
        return 2 * __builtin_popcount(mask) + __builtin_popcount(mask & weight3)
             + 2 * __builtin_popcount(mask & weight4);
    };
    vector<Pair> pairs;
    for (int i = 0; i < nh; ++i) for (int j = i + 1; j < nh; ++j)
        pairs.push_back({i, j, add(helper[i].label, helper[j].label), helper[i].cost + helper[j].cost});
    const int P = pairs.size();
    vector<tuple<int,int,int>> johnson;
    for (int i = 0; i < P; ++i) for (int j = i + 1; j < P; ++j) {
        const Pair a = pairs[i], b = pairs[j];
        int common = -1;
        if (a.a == b.a || a.a == b.b) common = a.a;
        if (a.b == b.a || a.b == b.b) common = a.b;
        if (common >= 0) johnson.emplace_back(i, j, a.weight + b.weight - 2 * helper[common].cost);
    }
    need(P == 1275 && johnson.size() == 62475, "complete helper Johnson graph");
    cout << "[\n";
    bool comma = false;
    for (int grade = 2; grade <= 4; ++grade) {
        vector<Core> representatives(Q);
        uint64_t core_states = 0;
        auto insert = [&](uint32_t mask, int charge) {
            ++core_states;
            int w = weight(mask);
            auto &old = representatives[charge];
            if (make_pair(w, mask) < make_pair(old.weight, old.mask)) old = {mask, w};
        };
        if (mode == "bfs") {
            auto enumerate = [&](auto &&self, int next, int left, uint32_t mask, int charge) -> void {
                if (left == 0) { insert(mask, charge); return; }
                for (int i = next; i <= nc - left; ++i)
                    self(self, i + 1, left - 1, mask | (1u << i), add(charge, core[i].label));
            };
            enumerate(enumerate, 0, grade, 0, 0);
        } else {
            uint32_t mask = (1u << grade) - 1;
            while (mask < (1u << nc)) {
                int charge = 0;
                for (int i = nc - 1; i >= 0; --i) if (mask >> i & 1) charge = add(charge, core[i].label);
                insert(mask, charge);
                uint32_t low = mask & -mask, next = mask + low;
                mask = next | (((next ^ mask) / low) >> 2);
            }
        }
        const uint64_t expected_states[] = {0, 0, 378, 3276, 20475};
        const int expected_images[] = {0, 0, 351, 1646, 2048};
        int image_size = count_if(representatives.begin(), representatives.end(), [](Core c) { return c.weight != 999; });
        need(core_states == expected_states[grade] && image_size == expected_images[grade], "complete core grade universe");
        vector<uint8_t> work(Q * Q, 255), support(Q * Q, 255);
        for (int a = 0; a < Q; ++a) if (representatives[a].weight != 999)
            for (int b = 0; b < Q; ++b) if (representatives[b].weight != 999) {
                uint32_t changed = representatives[a].mask ^ representatives[b].mask;
                int slots = __builtin_popcount(changed);
                if (slots <= 6) { work[a * Q + b] = weight(changed); support[a * Q + b] = slots; }
            }
        for (int sign : {1, -1}) {
            uint64_t vertices = 0, edges = 0, digest = UINT64_C(1469598103934665603);
            int min_vertices = P, max_vertices = 0, max_clock = 0, max_support = 0;
            for (int total = 0; total < Q; ++total) {
                vector<int> charges(P, -1), parent(P), degree(P);
                vector<vector<int>> adjacency;
                if (mode == "bfs") adjacency.resize(P);
                iota(parent.begin(), parent.end(), 0);
                auto find = [&](int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
                int n = 0, start = -1, components = 0;
                for (int h = 0; h < P; ++h) {
                    int q = add(total, neg(pairs[h].charge));
                    if (sign == -1) q = neg(q);
                    if (representatives[q].weight != 999) {
                        charges[h] = q; start = h; ++n; ++components;
                        need(add(sign == 1 ? q : neg(q), pairs[h].charge) == total, "both prescribed charge coordinates");
                    }
                }
                need(n > 0, "all total charges attained");
                uint64_t local_edges = 0;
                for (auto [a, b, helper_work] : johnson) {
                    if (charges[a] < 0 || charges[b] < 0) continue;
                    int index = charges[a] * Q + charges[b];
                    int w = work[index] + helper_work;
                    if (w > 33) continue;
                    int slots = support[index] + 2;
                    need(slots <= 8, "actual simultaneous slot bound");
                    ++local_edges; ++degree[a]; ++degree[b];
                    max_clock = max(max_clock, w); max_support = max(max_support, slots);
                    digest = mix(digest, uint64_t(total) * P * P + uint64_t(a) * P + b);
                    if (mode == "bfs") { adjacency[a].push_back(b); adjacency[b].push_back(a); }
                    else { int ra = find(a), rb = find(b); if (ra != rb) { parent[ra] = rb; --components; } }
                }
                if (mode == "bfs") {
                    vector<bool> seen(P);
                    queue<int> queue;
                    queue.push(start); seen[start] = true;
                    int reached = 0;
                    while (!queue.empty()) { int a = queue.front(); queue.pop(); ++reached;
                        for (int b : adjacency[a]) if (!seen[b]) { seen[b] = true; queue.push(b); }
                    }
                    need(reached == n, "complete helper-sheet graph connected by BFS");
                } else need(components == 1, "complete helper-sheet graph connected by DSU");
                digest = mix(mix(digest, n), local_edges);
                vertices += n; edges += local_edges;
                min_vertices = min(min_vertices, n); max_vertices = max(max_vertices, n);
            }
            if (comma) cout << ",\n";
            comma = true;
            cout << "  {\"grade\":" << grade << ",\"sign\":" << sign
                 << ",\"core_states\":" << core_states << ",\"core_charge_values\":" << image_size
                 << ",\"total_charge_fibres\":2048,\"sheet_vertices\":" << vertices << ",\"sheet_edges\":" << edges
                 << ",\"minimum_sheets\":" << min_vertices << ",\"maximum_sheets\":" << max_vertices
                 << ",\"maximum_edge_clocks\":" << max_clock << ",\"maximum_edge_slots\":" << max_support
                 << ",\"edge_digest\":\"" << digest << "\",\"all_connected\":true}";
        }
    }
    cout << "\n]\n";
}
