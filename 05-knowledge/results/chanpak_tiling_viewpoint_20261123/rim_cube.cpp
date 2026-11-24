// Enumerate one antipodal chord assignment and its physical rim edits.
// All requests in this replay have order at most eight.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

using Rows = std::vector<int>;
using Heights = std::set<int>;

int paths(const Rows& rows) {
    const int n = rows.size(), limit = 1 << n;
    std::array<std::array<int,8>,256> ways{};
    for (int i = 0; i < n; ++i) ways[1 << i][i] = 1;
    for (int used = 1; used < limit; ++used) {
        for (int bits = used; bits; bits &= bits-1) {
            int u = __builtin_ctz(bits), count = ways[used][u];
            if (!count) continue;
            for (int next = rows[u] & ~used; next; next &= next-1) {
                int v = __builtin_ctz(next);
                ways[used | (1 << v)][v] += count;
            }
        }
    }
    int total = 0;
    for (int count : ways[limit-1]) total += count;
    return total;
}

void toggle(Rows& rows, int i) {
    int j = (i + 1) % rows.size();
    rows[i] ^= 1 << j;
    rows[j] ^= 1 << i;
}

Rows closure(int n, int chord_bits) {
    int m = n/2;
    Rows rows(n);
    for (int i=0; i<n; ++i) rows[i] |= 1 << ((i+1)%n);
    int bit = 0;
    for (int distance=2; distance<=m; ++distance) {
        for (int i=0; i<m; ++i, ++bit) {
            int a=i, b=(i+distance)%n;
            if ((chord_bits >> bit) & 1) std::swap(a,b);
            rows[a] |= 1 << b;
            rows[(b+m)%n] |= 1 << ((a+m)%n);
        }
    }
    return rows;
}

void show(const Heights& list) {
    for (int x : list) std::cout << ' ' << x;
    std::cout << '\n';
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const int n = std::atoi(argv[1]);
    if (n != 4 && n != 6 && n != 8) return 2;
    int m=n/2;
    Heights unit, balanced;
    for (int chords=0; chords < (1 << (m*(m-1))); ++chords) {
        Rows base = closure(n,chords);
        int h=paths(base);
        for (int edge=0; edge<n; ++edge) {
            toggle(base,edge);
            int k=paths(base);
            if (std::abs(h-k)==2) unit.insert(std::min(h,k));
            toggle(base,edge);
        }
        for (int side=0; side<(1 << m); ++side) {
            Rows center=base;
            for (int i=0; i<m; ++i) toggle(center,(side >> i & 1) ? i : i+m);
            int first=paths(center);
            for (int i=0; i<m; ++i) {
                toggle(center,i);
                int middle=paths(center);
                if (std::abs(first-middle)==2) {
                    toggle(center,i+m);
                    int last=paths(center);
                    if (last-first==2*(middle-first)) balanced.insert(std::min(first,last));
                    toggle(center,i+m);
                }
                toggle(center,i);
            }
        }
    }
    show(unit);
    show(balanced);
}
