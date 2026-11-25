// Classify girth-five p-sheet covers with a lifted p-cycle on the branches.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

template<int p>
class RotatingCover {
    using Perm = std::array<unsigned char, p>;
    std::vector<Perm> permutations, inverses;
    std::map<std::vector<int>, int> class_rep;
    std::array<const Perm*, p> edge{};
    std::array<std::uint64_t, (p+1)/2> tried{};
    std::uint64_t passes = 0;

    static std::vector<int> cycle_type(const Perm& perm) {
        unsigned seen = 0;
        std::vector<int> lengths;
        for (int start=0; start<p; ++start) {
            if ((seen >> start) & 1) continue;
            int x=start, length=0;
            do {
                seen |= 1U << x;
                x=perm[x];
                ++length;
            } while (x != start);
            lengths.push_back(length);
        }
        std::sort(lengths.begin(), lengths.end());
        return lengths;
    }

    bool partial_routes_ok() const {
        // Translation invariance makes the source branch 0. Reversing all routes
        // gives inverse permutations, so it suffices to test half the destinations.
        for (int x=0; x<p; ++x) {
            for (int dest=1; dest<=(p-1)/2; ++dest) {
                unsigned seen = edge[dest] ? 1U << (*edge[dest])[x] : 0;
                for (int via=1; via<p; ++via) {
                    int step=(dest+p-via)%p;
                    if (!edge[via] || !edge[step]) continue;
                    unsigned bit=1U << (*edge[step])[(*edge[via])[x]];
                    if (seen & bit) return false;
                    seen |= bit;
                }
            }
        }
        return true;
    }

    void set(int offset, int index) {
        edge[offset]=&permutations[index];
        edge[p-offset]=&inverses[index];
        ++tried[offset];
        if (partial_routes_ok()) {
            if (offset==(p-1)/2) ++passes;
            else for (int next=0; next<int(permutations.size()); ++next)
                set(offset+1, next);
        }
        edge[offset]=edge[p-offset]=nullptr;
    }

public:
    std::uint64_t run() {
        Perm a{};
        for (unsigned char i=0; i<p; ++i) a[i]=i;
        do {
            class_rep.emplace(cycle_type(a), int(permutations.size()));
            permutations.push_back(a);
            Perm inv{};
            for (int i=0; i<p; ++i) inv[a[i]]=i;
            inverses.push_back(inv);
        } while (std::next_permutation(a.begin(), a.end()));
        for (const auto& entry : class_rep) set(1, entry.second);
        std::cout<<"p="<<p<<" tried";
        for (int i=1; i<=(p-1)/2; ++i) std::cout<<" "<<tried[i];
        std::cout<<" full tuples passing "<<passes<<"\n";
        return passes;
    }
};

int main() {
    auto three=RotatingCover<3>().run();
    auto five=RotatingCover<5>().run();
    auto seven=RotatingCover<7>().run();
    return three==0 && five==11 && seven==0 ? 0 : 1;
}
