#include "lrc_rectangle_common.hpp"

#include <chrono>

using namespace lrc_rectangle;

namespace {

std::vector<u32> sorted_candidates() {
    std::vector<OrderedRepair> universe;
    universe.reserve(REPAIR_COUNT);
    u32 mask=(u32{1}<<8)-1;
    const u32 limit=u32{1}<<30;
    while (mask!=0 && mask<limit) {
        universe.push_back({splitmix64(static_cast<u64>(mask)^ORDER_SEED),mask});
        const u32 next=next_mask(mask);
        if (next<=mask) break;
        mask=next;
    }
    require(universe.size()==REPAIR_COUNT,"repair universe changed");
    std::sort(universe.begin(),universe.end(),[](const auto& a,const auto& b) {
        return a.key!=b.key ? a.key<b.key : a.mask<b.mask;
    });
    std::vector<u32> result;
    result.reserve(CANDIDATES);
    for (std::size_t i=0;i<CANDIDATES;++i) result.push_back(universe[i].mask);
    return result;
}

struct Arc { i64 left=0; i64 right=0; i64 prefix=0; };

struct PrimitivePair {
    int q=0;
    int r=0;
    i64 g=0;
    i64 u=0;
    i64 v=0;
    i64 grid=0;
    i64 safe_ticks=0;
    std::vector<Arc> arcs;
};

i64 mod_tick(i128 value,i64 modulus) {
    value%=modulus;
    if (value<0) value+=modulus;
    return static_cast<i64>(value);
}

bool primitive_safe_midpoint(i64 speed,i64 grid,i64 left,i64 right) {
    i128 residue=static_cast<i128>(speed)*(left+right);
    residue%=static_cast<i128>(2)*grid;
    if (residue<0) residue+=static_cast<i128>(2)*grid;
    return static_cast<i128>(7)*residue>=grid &&
           static_cast<i128>(7)*residue<=static_cast<i128>(13)*grid;
}

PrimitivePair build_primitive(int q,int r) {
    require(Q0<=q && q<=Q1 && R0<=r && r<=R1 && q<r,"pair outside rectangle");
    PrimitivePair p;
    p.q=q;
    p.r=r;
    p.g=std::gcd(q,r);
    p.u=q/p.g;
    p.v=r/p.g;
    const i128 grid=static_cast<i128>(14)*p.u*p.v;
    require(grid<=std::numeric_limits<i64>::max(),"primitive grid overflow");
    p.grid=static_cast<i64>(grid);
    std::vector<i64> walls={0,p.grid};
    walls.reserve(static_cast<std::size_t>(2*(p.u+p.v)+2));
    for (i64 i=0;i<p.u;++i) {
        walls.push_back(mod_tick(static_cast<i128>(p.v)*(14*i-1),p.grid));
        walls.push_back(mod_tick(static_cast<i128>(p.v)*(14*i+1),p.grid));
    }
    for (i64 i=0;i<p.v;++i) {
        walls.push_back(mod_tick(static_cast<i128>(p.u)*(14*i-1),p.grid));
        walls.push_back(mod_tick(static_cast<i128>(p.u)*(14*i+1),p.grid));
    }
    std::sort(walls.begin(),walls.end());
    walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
    require(walls.front()==0 && walls.back()==p.grid,"primitive wall endpoints lost");
    i64 running=0;
    for (std::size_t i=0;i+1<walls.size();++i) {
        const i64 left=walls[i],right=walls[i+1];
        const bool safe=primitive_safe_midpoint(p.u,p.grid,left,right) &&
                        primitive_safe_midpoint(p.v,p.grid,left,right);
        if (!safe) continue;
        if (!p.arcs.empty() && p.arcs.back().right==left) {
            p.arcs.back().right=right;
            running+=right-left;
        } else {
            p.arcs.push_back({left,right,running});
            running+=right-left;
        }
    }
    p.safe_ticks=running;
    require(!p.arcs.empty() && p.safe_ticks>0,"primitive safe set empty");
    return p;
}

i128 safe_prefix(const PrimitivePair& p,i64 remainder) {
    require(0<=remainder && remainder<D,"primitive remainder invalid");
    const i128 y=static_cast<i128>(p.grid)*remainder;
    std::size_t low=0,high=p.arcs.size();
    while (low<high) {
        const std::size_t mid=low+(high-low)/2;
        if (static_cast<i128>(p.arcs[mid].right)*D<=y) low=mid+1;
        else high=mid;
    }
    if (low==p.arcs.size()) return static_cast<i128>(p.safe_ticks)*D;
    const Arc& arc=p.arcs[low];
    i128 answer=static_cast<i128>(arc.prefix)*D;
    const i128 left=static_cast<i128>(arc.left)*D;
    if (y>left) answer+=std::min(y,static_cast<i128>(arc.right)*D)-left;
    return answer;
}

i128 primitive_integral(const PrimitivePair& p,i128 z) {
    require(z>=0,"negative primitive coordinate");
    const i128 whole=z/D;
    const i64 remainder=static_cast<i64>(z%D);
    return whole*p.safe_ticks*D+safe_prefix(p,remainder);
}

struct Mass { i128 numerator=0; i128 denominator=1; i128 margin=0; };

Mass exact_mass(const PrimitivePair& p,const Geometry& geometry) {
    i128 numerator=0;
    for (const auto& [left,right]:geometry.components) {
        numerator+=primitive_integral(p,static_cast<i128>(p.g)*right)-
                   primitive_integral(p,static_cast<i128>(p.g)*left);
    }
    const i128 denominator=static_cast<i128>(p.grid)*D*p.g;
    return {numerator,denominator,static_cast<i128>(63)*numerator-4*denominator};
}

bool smaller_gap(i128 margin_a,i128 den_a,i128 margin_b,i128 den_b) {
    return product_less_exact(margin_a,den_b,margin_b,den_a);
}

} // namespace

int main() {
    try {
        const auto started=std::chrono::steady_clock::now();
        check_u256_products();
        const auto pairs=rectangle_pairs();
        const auto cells=build_pool_cells();
        const auto candidates=sorted_candidates();
        std::vector<Geometry> geometry;
        geometry.reserve(CANDIDATES);
        for (u32 repair:candidates) geometry.push_back(build_geometry(repair,cells));

        Fnv64 pair_fnv,candidate8_fnv,candidate16_fnv,matrix_fnv;
        for (auto [q,r]:pairs) { pair_fnv.word(q); pair_fnv.word(r); }
        for (std::size_t i=0;i<candidates.size();++i) {
            if (i<HOSTILE_BUDGET) candidate8_fnv.word(candidates[i]);
            candidate16_fnv.word(candidates[i]);
        }

        std::vector<unsigned char> common(CANDIDATES,1);
        std::vector<i128> best_margin(CANDIDATES,0),best_den(CANDIDATES,1);
        std::vector<int> best_q(CANDIDATES,0),best_r(CANDIDATES,0);
        std::vector<unsigned char> best_set(CANDIDATES,0);
        u64 nonnegative=0;
        u64 equalities=0;
        u64 matrix_cells=0;
        for (auto [q,r]:pairs) {
            const PrimitivePair primitive=build_primitive(q,r);
            for (std::size_t index=0;index<CANDIDATES;++index) {
                const Mass mass=exact_mass(primitive,geometry[index]);
                const bool active=mass.margin>=0;
                matrix_fnv.byte(active?1:0);
                ++matrix_cells;
                if (active) ++nonnegative;
                if (mass.margin==0) ++equalities;
                if (common[index]) {
                    if (!active) {
                        common[index]=0;
                    } else if (!best_set[index] ||
                               smaller_gap(mass.margin,mass.denominator,
                                           best_margin[index],best_den[index])) {
                        best_set[index]=1;
                        best_margin[index]=mass.margin;
                        best_den[index]=mass.denominator;
                        best_q[index]=q;
                        best_r[index]=r;
                    }
                }
            }
        }
        require(matrix_cells==PAIR_COUNT*CANDIDATES,"matrix cell census changed");

        std::vector<u32> hostile_deck,deck;
        std::vector<std::size_t> deck_indices;
        Fnv64 hostile_fnv,deck_fnv;
        for (std::size_t i=0;i<CANDIDATES;++i) if (common[i]) {
            if (i<HOSTILE_BUDGET) {
                hostile_deck.push_back(candidates[i]);
                hostile_fnv.word(candidates[i]);
            }
            deck.push_back(candidates[i]);
            deck_indices.push_back(i);
            deck_fnv.word(candidates[i]);
        }
        require(hostile_deck.size()==2572,"hostile common count changed");
        require(hex64(hostile_fnv.value())=="44a01e1ab114723e","hostile deck FNV changed");
        require(deck.size()==5257,"theorem common count changed");
        require(hex64(deck_fnv.value())=="60f329212844f8ac","theorem deck FNV changed");

        bool global_set=false;
        std::size_t weakest_index=0;
        for (std::size_t index:deck_indices) {
            require(best_margin[index]>0,"common repair has nonpositive weakest margin");
            if (!global_set || smaller_gap(best_margin[index],best_den[index],
                                           best_margin[weakest_index],best_den[weakest_index])) {
                global_set=true;
                weakest_index=index;
            }
        }
        const PrimitivePair weakest_pair=build_primitive(best_q[weakest_index],best_r[weakest_index]);
        const Mass weakest_mass=exact_mass(weakest_pair,geometry[weakest_index]);
        require(weakest_mass.margin==best_margin[weakest_index] &&
                weakest_mass.denominator==best_den[weakest_index],
                "weakest mass replay changed");

        const auto bodies=all_bodies();
        const BodyScan hostile_scan=parallel_body_scan(hostile_deck,bodies);
        const BodyScan theorem_scan=parallel_body_scan(deck,bodies);
        require(hostile_scan.failures==7,"hostile body failure count changed");
        require(theorem_scan.failures==0,"theorem deck misses a body");
        require(theorem_scan.bodies==BODY_COUNT,"theorem body universe incomplete");

        const double seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();
        std::cout << "LRC14_2D_RECTANGLE_PRIMARY_V1\n"
                  << "RECT Q " << Q0 << ".." << Q1 << " R " << R0 << ".." << R1
                  << " PAIRS " << pairs.size() << " PAIR_FNV " << hex64(pair_fnv.value()) << "\n"
                  << "CANDIDATES_8192_FNV " << hex64(candidate8_fnv.value())
                  << " CANDIDATES_16384_FNV " << hex64(candidate16_fnv.value()) << "\n"
                  << "MATRIX_CELLS " << matrix_cells << " NONNEGATIVE " << nonnegative
                  << " EQUALITIES " << equalities << " SIGN_BYTE_FNV "
                  << hex64(matrix_fnv.value()) << "\n"
                  << "HOSTILE_COMMON " << hostile_deck.size() << " FNV "
                  << hex64(hostile_fnv.value()) << " BODY_FAILURES " << hostile_scan.failures
                  << " FIRST_FAILURE 0x" << std::hex << hostile_scan.first_failure << std::dec
                  << " {" << labels(hostile_scan.first_failure) << "} CHECKS "
                  << hostile_scan.checks << " MAX_CHECKS " << hostile_scan.max_checks << "\n"
                  << "THEOREM_COMMON " << deck.size() << " FNV " << hex64(deck_fnv.value())
                  << " BODY_FAILURES " << theorem_scan.failures << " CHECKS "
                  << theorem_scan.checks << " MAX_CHECKS " << theorem_scan.max_checks
                  << " WORST_BODY 0x" << std::hex << theorem_scan.worst_body << std::dec
                  << " {" << labels(theorem_scan.worst_body) << "} WITNESS {"
                  << labels(theorem_scan.worst_repair) << "}\n"
                  << "WEAKEST_NORMALIZED_GAP_PAIR " << best_q[weakest_index] << ','
                  << best_r[weakest_index] << " REPAIR 0x" << std::hex
                  << geometry[weakest_index].repair << std::dec << " {"
                  << labels(geometry[weakest_index].repair) << "} MASS "
                  << fraction(weakest_mass.numerator,weakest_mass.denominator)
                  << " GAP " << fraction(weakest_mass.margin,
                                          static_cast<i128>(63)*weakest_mass.denominator)
                  << " RAW_MARGIN " << decimal(weakest_mass.margin) << "\n"
                  << "SECONDS " << std::fixed << std::setprecision(6) << seconds << "\n"
                  << "VERDICT PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << "\n";
        return 1;
    }
}
