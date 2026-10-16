// Exact point-envelope event sweep. All labels and clocks come from the
// authenticated Python native graph; no Collatz oracle is used here.
#include <cstdint>
#include <vector>
#include <cstdlib>

extern "C" int point_norms(std::int64_t X, std::int64_t ncols,
                           const std::int64_t* label,
                           const std::int64_t* value,
                           std::int64_t* out) {
    std::vector<std::int64_t> row(ncols, 0);
    for (std::int64_t n=X+1; n<=4*X; ++n) row[label[n]] += value[n];
    std::int64_t norm=0;
    for (auto v: row) norm += std::llabs(v);
    for (std::int64_t N=X; N<2*X; ++N) {
        out[N-X]=norm;
        if (N+1==2*X) break;
        // The fourth entering term and the departing term cancel by two
        // initial halving edges. Check their signed values and endpoints.
        if (label[4*N+4]!=label[N+1] || value[4*N+4]!=value[N+1]) return 1;
        std::int64_t ids[3], change[3]; int used=0;
        for (int j=1; j<=3; ++j) {
            const auto e=label[4*N+j], v=value[4*N+j];
            int a=0; while (a<used && ids[a]!=e) ++a;
            if (a==used) { ids[used]=e; change[used]=0; ++used; }
            change[a]+=v;
        }
        std::int64_t variation=0;
        for (int a=0; a<used; ++a) {
            variation+=std::llabs(change[a]);
            auto &v=row[ids[a]];
            norm-=std::llabs(v); v+=change[a]; norm+=std::llabs(v);
        }
        if (variation>4) return 2;
    }
    return 0;
}
