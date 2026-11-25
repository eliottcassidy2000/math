// Count disjoint induced 13-cycles through the known Hoffman–Singleton cycle orbit.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using Mask = std::uint64_t;
using Rows = std::array<Mask, 50>;
Mask bit(int x) { return Mask{1} << x; }

Rows hoffman_singleton() {
    Rows graph{};
    auto join = [&](int a, int b) {
        graph[a] |= bit(b);
        graph[b] |= bit(a);
    };
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            join(5*i+j, 5*i+(j+1)%5);
            join(25+5*i+j, 25+5*i+(j+2)%5);
            for (int k = 0; k < 5; ++k) join(5*i+j, 25+5*k+(i*k+j)%5);
        }
    }
    return graph;
}

void extend(const Rows& g, int root, int first, int last, int size,
            Mask used, Mask free, std::vector<Mask>& cycles) {
    if (size == 12) {
        Mask close = free & g[last] & g[root] & ~(bit(first+1)-1);
        while (close) {
            Mask next = close & -close;
            cycles.push_back(used | next);
            close ^= next;
        }
        return;
    }
    Mask options = free & g[last] & ~g[root];
    while (options) {
        int v = __builtin_ctzll(options);
        options &= options-1;
        extend(g, root, first, v, size+1, used | bit(v), free & ~g[last] & ~bit(v), cycles);
    }
}

int main() {
    Rows g = hoffman_singleton();
    std::vector<Mask> cycles;
    for (int r = 0; r < 50; ++r) {
        Mask above = (bit(50)-1) & ~(bit(r+1)-1);
        for (int f = r+1; f < 50; ++f) {
            if (!(g[r] & bit(f))) continue;
            Mask seconds = above & g[f] & ~g[r];
            while (seconds) {
                int v = __builtin_ctzll(seconds);
                seconds &= seconds-1;
                extend(g,r,f,v,3,bit(r)|bit(f)|bit(v),above & ~g[f] & ~bit(v) & ~bit(f),cycles);
            }
        }
    }
    Mask fixed = 0;
    for (int v : {0,1,2,3,28,8,7,6,44,15,16,17,30}) fixed |= bit(v);
    std::array<int,14> intersections{};
    std::vector<Mask> disjoint;
    for (Mask c : cycles) {
        int overlap = __builtin_popcountll(c & fixed);
        ++intersections[overlap];
        if (!overlap) disjoint.push_back(c);
    }
    long long extensions = 0, dead = 0;
    int least = 1000000, most = 0;
    for (Mask c : disjoint) {
        int compatible = 0;
        for (Mask d : disjoint) if (!(c & d)) ++compatible;
        if (!compatible) ++dead;
        if (compatible < least) least = compatible;
        if (compatible > most) most = compatible;
        extensions += compatible;
    }
    if (cycles.size()!=126000 || intersections[13]!=1 || extensions%2) return 1;
    std::cout << "intersection counts";
    for (int count : intersections) std::cout << " " << count;
    std::cout << "\nunordered pairs " << cycles.size()*disjoint.size()/2
              << "\nunordered triples " << cycles.size()*extensions/6
              << "\nthird choices over each pair " << least << ".." << most
              << " dead pairs per fixed cycle " << dead << "\n";
}
