#include "lrc_rectangle_common.hpp"

#include <chrono>
#include <functional>

using namespace lrc_rectangle;

namespace {

bool ordered_less(const OrderedRepair& a,const OrderedRepair& b) {
    return a.key!=b.key ? a.key<b.key : a.mask<b.mask;
}

std::vector<u32> bounded_heap_candidates() {
    struct Less {
        bool operator()(const OrderedRepair& a,const OrderedRepair& b) const {
            return ordered_less(a,b);
        }
    };
    std::priority_queue<OrderedRepair,std::vector<OrderedRepair>,Less> heap;
    u32 mask=(u32{1}<<8)-1;
    const u32 limit=u32{1}<<30;
    u64 count=0;
    while (mask!=0 && mask<limit) {
        const OrderedRepair item{splitmix64(static_cast<u64>(mask)^ORDER_SEED),mask};
        if (heap.size()<CANDIDATES) heap.push(item);
        else if (ordered_less(item,heap.top())) { heap.pop(); heap.push(item); }
        ++count;
        const u32 next=next_mask(mask);
        if (next<=mask) break;
        mask=next;
    }
    require(count==REPAIR_COUNT,"repair universe changed");
    std::vector<OrderedRepair> selected;
    selected.reserve(CANDIDATES);
    while (!heap.empty()) { selected.push_back(heap.top()); heap.pop(); }
    std::sort(selected.begin(),selected.end(),ordered_less);
    std::vector<u32> result;
    result.reserve(CANDIDATES);
    for (const auto& item:selected) result.push_back(item.mask);
    return result;
}

i128 lcm128(i128 a,i128 b) {
    require(a>0 && b>0,"nonpositive lcm input");
    return (a/gcd128(a,b))*b;
}

struct LiteralArc { i128 left=0; i128 right=0; i128 prefix=0; };

struct LiteralPair {
    int q=0;
    int r=0;
    i128 grid=1;
    i128 pool_scale=1;
    std::vector<LiteralArc> arcs;
};

std::vector<std::pair<i128,i128>> speed_intervals(int speed,i128 grid) {
    const i128 denominator=static_cast<i128>(14)*speed;
    require(grid%denominator==0,"literal grid not speed-divisible");
    const i128 quantum=grid/denominator;
    std::vector<std::pair<i128,i128>> result;
    result.reserve(speed);
    for (int tooth=0;tooth<speed;++tooth) {
        result.push_back({static_cast<i128>(14*tooth+1)*quantum,
                          static_cast<i128>(14*tooth+13)*quantum});
    }
    return result;
}

LiteralPair build_literal_pair(int q,int r) {
    require(Q0<=q && q<=Q1 && R0<=r && r<=R1 && q<r,"literal pair outside rectangle");
    LiteralPair p;
    p.q=q;
    p.r=r;
    p.grid=lcm128(lcm128(D,static_cast<i128>(14)*q),static_cast<i128>(14)*r);
    require(p.grid%D==0,"literal grid lost pool denominator");
    p.pool_scale=p.grid/D;
    const auto q_arcs=speed_intervals(q,p.grid);
    const auto r_arcs=speed_intervals(r,p.grid);
    std::size_t i=0,j=0;
    i128 running=0;
    while (i<q_arcs.size() && j<r_arcs.size()) {
        const i128 left=std::max(q_arcs[i].first,r_arcs[j].first);
        const i128 right=std::min(q_arcs[i].second,r_arcs[j].second);
        if (left<right) {
            if (!p.arcs.empty() && p.arcs.back().right==left) {
                p.arcs.back().right=right;
                running+=right-left;
            } else {
                p.arcs.push_back({left,right,running});
                running+=right-left;
            }
        }
        if (q_arcs[i].second<r_arcs[j].second) ++i;
        else if (r_arcs[j].second<q_arcs[i].second) ++j;
        else { ++i; ++j; }
    }
    require(!p.arcs.empty() && running>0,"literal pair safe set empty");
    return p;
}

i128 literal_prefix(const LiteralPair& p,i128 tick) {
    require(0<=tick && tick<=p.grid,"literal prefix outside circle");
    std::size_t low=0,high=p.arcs.size();
    while (low<high) {
        const std::size_t mid=low+(high-low)/2;
        if (p.arcs[mid].right<=tick) low=mid+1;
        else high=mid;
    }
    if (low==p.arcs.size()) {
        const auto& arc=p.arcs.back();
        return arc.prefix+(arc.right-arc.left);
    }
    const LiteralArc& arc=p.arcs[low];
    i128 answer=arc.prefix;
    if (tick>arc.left) answer+=std::min(tick,arc.right)-arc.left;
    return answer;
}

struct Mass { i128 numerator=0; i128 denominator=1; i128 margin=0; };

Mass literal_mass(const LiteralPair& p,const Geometry& geometry) {
    i128 numerator=0;
    for (const auto& [left,right]:geometry.components) {
        const i128 a=static_cast<i128>(left)*p.pool_scale;
        const i128 b=static_cast<i128>(right)*p.pool_scale;
        numerator+=literal_prefix(p,b)-literal_prefix(p,a);
    }
    return {numerator,p.grid,static_cast<i128>(63)*numerator-4*p.grid};
}

bool smaller_gap(i128 margin_a,i128 den_a,i128 margin_b,i128 den_b) {
    return product_less_exact(margin_a,den_b,margin_b,den_a);
}

BodyScan recursive_body_scan(const std::vector<u32>& deck) {
    BodyScan report;
    std::function<void(int,int,u32)> visit=[&](int next,int need,u32 body) {
        if (need==0) {
            u64 used=0;
            u32 witness=0;
            for (u32 repair:deck) {
                ++used;
                if ((repair&body)==0) { witness=repair; break; }
            }
            ++report.bodies;
            report.checks+=used;
            if (witness==0) {
                ++report.failures;
                report.first_failure=std::min(report.first_failure,body);
            }
            if (used>report.max_checks ||
                (used==report.max_checks && body<report.worst_body)) {
                report.max_checks=used;
                report.worst_body=body;
                report.worst_repair=witness;
            }
            return;
        }
        for (int bit=next;bit<=30-need;++bit)
            visit(bit+1,need-1,body|(u32{1}<<bit));
    };
    visit(0,9,0);
    require(report.bodies==BODY_COUNT,"recursive body census changed");
    return report;
}

} // namespace

int main() {
    try {
        const auto started=std::chrono::steady_clock::now();
        check_u256_products();
        const auto pairs=rectangle_pairs();
        const auto cells=build_pool_cells();
        const auto candidates=bounded_heap_candidates();
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
        i128 largest_grid=0;
        int largest_q=0,largest_r=0;
        for (auto [q,r]:pairs) {
            const LiteralPair literal=build_literal_pair(q,r);
            if (literal.grid>largest_grid) {
                largest_grid=literal.grid;
                largest_q=q;
                largest_r=r;
            }
            for (std::size_t index=0;index<CANDIDATES;++index) {
                const Mass mass=literal_mass(literal,geometry[index]);
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
        require(matrix_cells==PAIR_COUNT*CANDIDATES,"literal matrix cell census changed");

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
        require(hostile_deck.size()==2572,"literal hostile common count changed");
        require(hex64(hostile_fnv.value())=="44a01e1ab114723e","literal hostile FNV changed");
        require(deck.size()==5257,"literal theorem common count changed");
        require(hex64(deck_fnv.value())=="60f329212844f8ac","literal theorem FNV changed");

        bool global_set=false;
        std::size_t weakest_index=0;
        for (std::size_t index:deck_indices) {
            require(best_margin[index]>0,"literal common repair has nonpositive margin");
            if (!global_set || smaller_gap(best_margin[index],best_den[index],
                                           best_margin[weakest_index],best_den[weakest_index])) {
                global_set=true;
                weakest_index=index;
            }
        }
        const LiteralPair weakest_pair=build_literal_pair(best_q[weakest_index],best_r[weakest_index]);
        const Mass weakest_mass=literal_mass(weakest_pair,geometry[weakest_index]);
        require(weakest_mass.margin==best_margin[weakest_index] &&
                weakest_mass.denominator==best_den[weakest_index],
                "literal weakest mass replay changed");

        const BodyScan hostile_scan=recursive_body_scan(hostile_deck);
        const BodyScan theorem_scan=recursive_body_scan(deck);
        require(hostile_scan.failures==7,"literal hostile body failures changed");
        require(theorem_scan.failures==0,"literal theorem deck misses a body");

        const double seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();
        std::cout << "LRC14_2D_RECTANGLE_LITERAL_V1\n"
                  << "RECT Q " << Q0 << ".." << Q1 << " R " << R0 << ".." << R1
                  << " PAIRS " << pairs.size() << " PAIR_FNV " << hex64(pair_fnv.value()) << "\n"
                  << "CANDIDATES_8192_FNV " << hex64(candidate8_fnv.value())
                  << " CANDIDATES_16384_FNV " << hex64(candidate16_fnv.value()) << "\n"
                  << "MATRIX_CELLS " << matrix_cells << " NONNEGATIVE " << nonnegative
                  << " EQUALITIES " << equalities << " SIGN_BYTE_FNV "
                  << hex64(matrix_fnv.value()) << "\n"
                  << "LARGEST_LITERAL_GRID " << decimal(largest_grid) << " AT "
                  << largest_q << ',' << largest_r << "\n"
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
