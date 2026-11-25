// Count each unoriented induced cycle once, fixing its least label and first direction.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using Mask = std::uint64_t;
using Rows = std::vector<Mask>;
Mask bit(int v) { return Mask{1} << v; }

void add(Rows &graph, int a, int b) {
    graph[a] |= bit(b);
    graph[b] |= bit(a);
}

void extend(const Rows &g, int root, int first, int last, int vertices, Mask free,
            std::vector<long long> &counts) {
    Mask closers = g[last] & g[root] & free & ~(bit(first + 1) - 1);
    counts[vertices+1] += __builtin_popcountll(closers);
    Mask options = g[last] & free & ~g[root];
    while (options) {
        int v = __builtin_ctzll(options);
        options &= options-1;
        extend(g,root,first,v,vertices+1,free & ~g[last] & ~bit(v),counts);
    }
}

std::vector<long long> census(const Rows &g) {
    int n = g.size();
    std::vector<long long> counts(n+2);
    for (int root=0; root<n; ++root) {
        Mask above = (bit(n)-1) & ~(bit(root+1)-1);
        for (int first=root+1; first<n; ++first) {
            if (!(g[root] & bit(first))) continue;
            Mask second_options = g[first] & above & ~g[root];
            while (second_options) {
                int second = __builtin_ctzll(second_options);
                second_options &= second_options-1;
                extend(g,root,first,second,3,above & ~g[first] & ~bit(first) & ~bit(second),counts);
            }
        }
    }
    return counts;
}

int main() {
    Rows p(10), h(50);
    for (int j=0; j<5; ++j) {
        add(p,j,(j+1)%5);
        add(p,5+j,5+(j+2)%5);
        add(p,j,5+j);
        for (int i=0; i<5; ++i) {
            add(h,5*i+j,5*i+(j+1)%5);
            add(h,25+5*i+j,25+5*i+(j+2)%5);
            for (int k=0; k<5; ++k) add(h,5*i+j,25+5*k+(i*k+j)%5);
        }
    }
    for (const Rows &g : {p,h}) {
        auto counts = census(g);
        std::cout << "order " << g.size() << ", induced length and count\n";
        for (int i=4; i<=int(g.size()); ++i)
            if (counts[i]) std::cout << i << " " << counts[i] << "\n";
    }
}
