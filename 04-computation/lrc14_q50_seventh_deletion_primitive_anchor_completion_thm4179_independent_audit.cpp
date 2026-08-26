#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr int NEWCOMER = 50;
constexpr std::array<int, 3> ORIGINAL_ANCHOR = {120, 126, 143};
constexpr u64 FNV1A64_OFFSET = UINT64_C(0xcbf29ce484222325);
constexpr u64 FNV1A64_PRIME = UINT64_C(0x100000001b3);

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void force_binary_stdout() {
#ifdef _WIN32
    require(_setmode(_fileno(stdout), _O_BINARY) != -1,
            "failed to put stdout in binary mode");
#endif
}

class Fnv1a64 {
  public:
    void add_byte(std::uint8_t byte) {
        value_ ^= byte;
        value_ *= FNV1A64_PRIME;
    }

    void add_u64_le(u64 word) {
        for (int shift = 0; shift < 64; shift += 8) {
            add_byte(static_cast<std::uint8_t>((word >> shift) & UINT64_C(0xff)));
        }
    }

    void add_i64_le(i64 word) {
        add_u64_le(static_cast<u64>(word));
    }

    u64 value() const { return value_; }

  private:
    u64 value_ = FNV1A64_OFFSET;
};

std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex << std::nouppercase << std::setw(16) << std::setfill('0')
        << value;
    return out.str();
}

std::string triple_text(const std::array<int, 3>& triple) {
    return std::to_string(triple[0]) + "," + std::to_string(triple[1]) +
           "," + std::to_string(triple[2]);
}

std::string labels_of_local_mask(u32 mask, const std::vector<int>& labels) {
    std::ostringstream out;
    bool first = true;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if ((mask & (u32{1} << i)) == 0) continue;
        if (!first) out << ',';
        first = false;
        out << labels[i];
    }
    return out.str();
}

i64 exact_lcm(i64 left, i64 right) {
    const i64 divisor = std::gcd(left, right);
    const i128 value = static_cast<i128>(left / divisor) * right;
    require(value <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(value);
}

u64 binomial(int n, int k) {
    require(n >= 0 && k >= 0 && k <= n, "invalid binomial arguments");
    k = std::min(k, n - k);
    u64 result = 1;
    for (int i = 1; i <= k; ++i) {
        result = result * static_cast<u64>(n - k + i) / static_cast<u64>(i);
    }
    return result;
}

u32 next_combination(u32 mask) {
    const u32 low_bit = mask & (~mask + 1u);
    const u32 ripple = mask + low_bit;
    if (ripple == 0) return 0;
    return ripple | (((mask ^ ripple) >> 2) / low_bit);
}

struct Atom {
    i64 left;
    i64 right;
    u32 failed_pool;
    bool failed_newcomer;
};

bool midpoint_is_safe(int speed, i64 left, i64 right, i64 denominator) {
    // The midpoint is (left+right)/(2*denominator).  If residue/(2D) is
    // its speed-times fractional part, safety is 1/14 <= residue/(2D)
    // <= 13/14.  No midpoint is a wall, but non-strict comparisons mirror
    // the definition exactly.
    const i64 period = 2 * denominator;
    const i64 residue = static_cast<i64>(
        (static_cast<i128>(speed) * (left + right)) % period);
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <=
               static_cast<i128>(13) * denominator;
}

struct Geometry {
    i64 denominator = 0;
    std::vector<i64> walls;
    std::vector<Atom> atoms;
    u64 atom_hash = 0;
    std::array<u64, 31> safe_hist{};
    std::array<u64, 31> failed_newcomer_hist{};
};

Geometry build_joint_geometry() {
    Geometry geometry;
    geometry.denominator = 1;
    for (int speed : POOL) {
        geometry.denominator = exact_lcm(geometry.denominator, 14LL * speed);
    }
    geometry.denominator =
        exact_lcm(geometry.denominator, 14LL * NEWCOMER);

    geometry.walls = {0, geometry.denominator};
    auto append_walls = [&](int speed) {
        const i64 unit = geometry.denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            geometry.walls.push_back((14LL * tooth + 1) * unit);
            geometry.walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : POOL) append_walls(speed);
    append_walls(NEWCOMER);
    std::sort(geometry.walls.begin(), geometry.walls.end());
    geometry.walls.erase(
        std::unique(geometry.walls.begin(), geometry.walls.end()),
        geometry.walls.end());

    require(geometry.walls.size() >= 2, "joint wall set is empty");
    i64 total_length = 0;
    Fnv1a64 atom_hash;
    geometry.atoms.reserve(geometry.walls.size() - 1);
    for (std::size_t index = 0; index + 1 < geometry.walls.size(); ++index) {
        const i64 left = geometry.walls[index];
        const i64 right = geometry.walls[index + 1];
        require(left < right, "joint walls are not strictly increasing");
        u32 failed_pool = 0;
        for (int pool_index = 0; pool_index < 30; ++pool_index) {
            if (!midpoint_is_safe(POOL[pool_index], left, right,
                                  geometry.denominator)) {
                failed_pool |= u32{1} << pool_index;
            }
        }
        const bool failed_newcomer = !midpoint_is_safe(
            NEWCOMER, left, right, geometry.denominator);
        geometry.atoms.push_back(
            {left, right, failed_pool, failed_newcomer});
        total_length += right - left;

        atom_hash.add_u64_le(static_cast<u64>(index));
        atom_hash.add_i64_le(left);
        atom_hash.add_i64_le(right);
        atom_hash.add_u64_le(failed_pool);
        atom_hash.add_u64_le(failed_newcomer ? 1 : 0);
        const int failure_size = std::popcount(failed_pool);
        if (failed_newcomer) {
            ++geometry.failed_newcomer_hist[failure_size];
        } else {
            ++geometry.safe_hist[failure_size];
        }
    }
    require(total_length == geometry.denominator,
            "joint atom lengths do not sum to one circle");
    geometry.atom_hash = atom_hash.value();
    return geometry;
}

bool is_divisor_complete(const std::array<int, 3>& anchor) {
    for (int divisor = 2; divisor <= 14; ++divisor) {
        bool covered = false;
        for (int label : anchor) covered = covered || label % divisor == 0;
        if (!covered) return false;
    }
    return true;
}

int anchor_gcd(const std::array<int, 3>& anchor) {
    return std::gcd(std::gcd(anchor[0], anchor[1]), anchor[2]);
}

std::vector<std::array<int, 3>> enumerate_anchors() {
    std::vector<std::array<int, 3>> anchors;
    for (int i = 0; i < 30; ++i) {
        for (int j = i + 1; j < 30; ++j) {
            for (int k = j + 1; k < 30; ++k) {
                const std::array<int, 3> anchor =
                    {POOL[i], POOL[j], POOL[k]};
                if (anchor_gcd(anchor) == 1 && is_divisor_complete(anchor)) {
                    anchors.push_back(anchor);
                }
            }
        }
    }
    return anchors;
}

u32 pool_mask_of_anchor(const std::array<int, 3>& anchor) {
    u32 mask = 0;
    for (int label : anchor) {
        const auto found = std::find(POOL.begin(), POOL.end(), label);
        require(found != POOL.end(), "anchor label is outside the pool");
        const int index = static_cast<int>(found - POOL.begin());
        mask |= u32{1} << index;
    }
    require(std::popcount(mask) == 3, "anchor does not have three labels");
    return mask;
}

u32 local_to_global_mask(u32 local_mask,
                         const std::vector<int>& local_to_pool) {
    u32 global_mask = 0;
    for (int local = 0; local < 27; ++local) {
        if ((local_mask & (u32{1} << local)) != 0) {
            global_mask |= u32{1} << local_to_pool[local];
        }
    }
    return global_mask;
}

struct Hypergraph {
    int arity = 0;
    u64 candidate_count = 0;
    u64 equality_count = 0;
    std::vector<u32> edges;
    i64 minimum_edge_delta = std::numeric_limits<i64>::max();
    i64 maximum_nonedge_delta = std::numeric_limits<i64>::min();
    u64 edge_hash = 0;
};

Hypergraph build_hypergraph(const Geometry& geometry,
                            const std::array<int, 3>& anchor,
                            const std::vector<int>& local_to_pool,
                            int arity) {
    require(local_to_pool.size() == 27, "optional-label map has wrong size");
    require(arity == 7, "THM-4179 audits only deletion arity seven");
    const u32 anchor_mask = pool_mask_of_anchor(anchor);

    std::unordered_map<u32, i64> failure_mass;
    failure_mass.reserve(geometry.atoms.size());
    for (const Atom& atom : geometry.atoms) {
        if (atom.failed_newcomer || (atom.failed_pool & anchor_mask) != 0) {
            continue;
        }
        u32 local_failure = 0;
        for (int local = 0; local < 27; ++local) {
            if ((atom.failed_pool & (u32{1} << local_to_pool[local])) != 0) {
                local_failure |= u32{1} << local;
            }
        }
        if (std::popcount(local_failure) <= arity) {
            failure_mass[local_failure] += atom.right - atom.left;
        }
    }

    Hypergraph hypergraph;
    hypergraph.arity = arity;
    Fnv1a64 edge_hash;
    edge_hash.add_u64_le(static_cast<u64>(arity));
    for (int pool_index : local_to_pool) {
        edge_hash.add_u64_le(static_cast<u64>(POOL[pool_index]));
    }

    const u32 limit = u32{1} << 27;
    u32 deletion = (u32{1} << arity) - 1;
    while (deletion != 0 && deletion < limit) {
        ++hypergraph.candidate_count;
        i64 mass = 0;
        u32 subset = deletion;
        while (true) {
            const auto found = failure_mass.find(subset);
            if (found != failure_mass.end()) mass += found->second;
            if (subset == 0) break;
            subset = (subset - 1) & deletion;
        }
        const i64 delta = static_cast<i64>(
            static_cast<i128>(63) * mass -
            static_cast<i128>(4) * geometry.denominator);
        if (delta == 0) ++hypergraph.equality_count;
        if (delta >= 0) {
            hypergraph.edges.push_back(deletion);
            hypergraph.minimum_edge_delta =
                std::min(hypergraph.minimum_edge_delta, delta);
            edge_hash.add_u64_le(
                local_to_global_mask(deletion, local_to_pool));
            edge_hash.add_i64_le(mass);
            edge_hash.add_i64_le(delta);
        } else {
            hypergraph.maximum_nonedge_delta =
                std::max(hypergraph.maximum_nonedge_delta, delta);
        }

        const u32 next = next_combination(deletion);
        if (next <= deletion) break;
        deletion = next;
    }
    require(hypergraph.candidate_count == binomial(27, arity),
            "deletion combination count mismatch");
    require(!hypergraph.edges.empty(), "repair hypergraph is empty");
    hypergraph.edge_hash = edge_hash.value();
    return hypergraph;
}

struct SearchResult {
    bool found = false;
    u32 cover = 0;
    u64 nodes = 0;
    u64 packing_prunes = 0;
    u64 memo_prunes = 0;
    std::size_t dead_states = 0;
};

class ExactCoverSearch {
  public:
    ExactCoverSearch(const std::vector<u32>& edges, int budget)
        : edges_(edges), budget_(budget) {}

    SearchResult run() {
        SearchResult result;
        result.found = search(0, budget_);
        result.cover = solution_;
        result.nodes = nodes_;
        result.packing_prunes = packing_prunes_;
        result.memo_prunes = memo_prunes_;
        result.dead_states = dead_.size();
        return result;
    }

  private:
    bool search(u32 cover, int remaining) {
        ++nodes_;
        if (dead_.find(cover) != dead_.end()) {
            ++memo_prunes_;
            return false;
        }

        std::array<int, 27> frequencies{};
        bool has_uncovered = false;
        for (u32 edge : edges_) {
            if ((edge & cover) != 0) continue;
            has_uncovered = true;
            for (u32 bits = edge; bits != 0; bits &= bits - 1) {
                ++frequencies[std::countr_zero(bits)];
            }
        }
        if (!has_uncovered) {
            solution_ = cover;
            return true;
        }
        if (remaining == 0) {
            dead_.insert(cover);
            return false;
        }

        u32 packed_vertices = 0;
        int packing_size = 0;
        for (u32 edge : edges_) {
            if ((edge & cover) != 0 || (edge & packed_vertices) != 0) continue;
            packed_vertices |= edge;
            ++packing_size;
            if (packing_size > remaining) {
                ++packing_prunes_;
                dead_.insert(cover);
                return false;
            }
        }

        u32 pivot = 0;
        u64 best_score = std::numeric_limits<u64>::max();
        for (u32 edge : edges_) {
            if ((edge & cover) != 0) continue;
            u64 score = 0;
            for (u32 bits = edge; bits != 0; bits &= bits - 1) {
                score += static_cast<u64>(
                    frequencies[std::countr_zero(bits)]);
            }
            if (score < best_score ||
                (score == best_score && (pivot == 0 || edge < pivot))) {
                best_score = score;
                pivot = edge;
            }
        }
        require(pivot != 0, "cover recursion failed to select a pivot edge");

        std::vector<int> branch_vertices;
        for (u32 bits = pivot; bits != 0; bits &= bits - 1) {
            branch_vertices.push_back(std::countr_zero(bits));
        }
        std::sort(branch_vertices.begin(), branch_vertices.end(),
                  [&](int left, int right) {
                      if (frequencies[left] != frequencies[right]) {
                          return frequencies[left] > frequencies[right];
                      }
                      return left < right;
                  });
        for (int vertex : branch_vertices) {
            if (search(cover | (u32{1} << vertex), remaining - 1)) {
                return true;
            }
        }
        dead_.insert(cover);
        return false;
    }

    const std::vector<u32>& edges_;
    int budget_;
    std::unordered_set<u32> dead_;
    u32 solution_ = 0;
    u64 nodes_ = 0;
    u64 packing_prunes_ = 0;
    u64 memo_prunes_ = 0;
};

bool is_cover(u32 cover, const std::vector<u32>& edges) {
    return std::all_of(edges.begin(), edges.end(),
                       [&](u32 edge) { return (cover & edge) != 0; });
}

struct MinimumCoverResult {
    int tau = -1;
    u32 cover = 0;
    std::vector<std::pair<int, SearchResult>> searches;
};

MinimumCoverResult exact_minimum_through_eight(
    const std::vector<u32>& edges) {
    MinimumCoverResult result;
    for (int budget = 0; budget <= 8; ++budget) {
        ExactCoverSearch search(edges, budget);
        SearchResult attempt = search.run();
        result.searches.emplace_back(budget, attempt);
        if (!attempt.found) continue;
        require(std::popcount(attempt.cover) == budget,
                "first successful budget returned a smaller cover");
        require(is_cover(attempt.cover, edges),
                "returned cover does not meet every repair edge");
        result.tau = budget;
        result.cover = attempt.cover;
        break;
    }
    require(result.tau >= 0, "transversal number exceeds eight");
    for (const auto& [budget, attempt] : result.searches) {
        if (budget < result.tau) {
            require(!attempt.found,
                    "smaller budget unexpectedly returned a cover");
        }
    }
    return result;
}

std::string search_summary(
    const std::vector<std::pair<int, SearchResult>>& searches) {
    std::ostringstream out;
    bool first = true;
    for (const auto& [budget, result] : searches) {
        if (!first) out << ';';
        first = false;
        out << budget << ':' << (result.found ? 'Y' : 'N') << ':'
            << result.nodes << ':' << result.packing_prunes << ':'
            << result.memo_prunes << ':' << result.dead_states;
    }
    return out.str();
}

std::string histogram_text(const std::array<u64, 31>& histogram) {
    std::ostringstream out;
    bool first = true;
    for (int size = 0; size <= 30; ++size) {
        if (histogram[size] == 0) continue;
        if (!first) out << ',';
        first = false;
        out << size << ':' << histogram[size];
    }
    return out.str();
}

void run_audit() {
    force_binary_stdout();
    require(POOL.size() == 30, "pool size changed");
    require(std::adjacent_find(POOL.begin(), POOL.end()) == POOL.end(),
            "pool has adjacent duplicate labels");

    const Geometry geometry = build_joint_geometry();
    const auto anchors = enumerate_anchors();
    require(anchors.size() == 6,
            "primitive divisor-complete anchor count is not six");
    require(std::find(anchors.begin(), anchors.end(), ORIGINAL_ANCHOR) !=
                anchors.end(),
            "original anchor triple is absent");

    std::cout << "AUDIT thm4179_q50_seventh_deletion_cleanroom_cpp20\n";
    std::cout << "POOL";
    for (int label : POOL) std::cout << ' ' << label;
    std::cout << "\n";
    std::cout << "NEWCOMER " << NEWCOMER << "\n";
    std::cout << "JOINT_GEOMETRY denominator=" << geometry.denominator
              << " walls=" << geometry.walls.size()
              << " atoms=" << geometry.atoms.size()
              << " labelled_atom_fnv1a64_le=" << hex64(geometry.atom_hash)
              << "\n";
    std::cout << "ATOM_HIST newcomer_safe_pool_failures="
              << histogram_text(geometry.safe_hist) << "\n";
    std::cout << "ATOM_HIST newcomer_failed_pool_failures="
              << histogram_text(geometry.failed_newcomer_hist) << "\n";
    std::cout << "PRIMITIVE_ANCHOR_COUNT " << anchors.size() << "\n";

    int audited_rows = 0;
    for (const auto& anchor : anchors) {
        const bool target =
            anchor == std::array<int, 3>{80, 143, 252} ||
            anchor == std::array<int, 3>{143, 240, 252};
        if (!target) continue;
        ++audited_rows;
        require(anchor_gcd(anchor) == 1,
                "enumerated anchor is not primitive");
        require(is_divisor_complete(anchor),
                "enumerated anchor is not divisor-complete");
        std::cout << "ANCHOR labels=" << triple_text(anchor)
                  << " original="
                  << (anchor == ORIGINAL_ANCHOR ? "yes" : "no")
                  << " gcd=" << anchor_gcd(anchor)
                  << " divisor_witnesses=";
        for (int divisor = 2; divisor <= 14; ++divisor) {
            if (divisor != 2) std::cout << ';';
            std::cout << divisor << ':';
            bool first = true;
            for (int label : anchor) {
                if (label % divisor != 0) continue;
                if (!first) std::cout << ',';
                first = false;
                std::cout << label;
            }
            require(!first, "divisor witness unexpectedly absent");
        }
        std::cout << "\n";

        const u32 anchor_mask = pool_mask_of_anchor(anchor);
        std::vector<int> optional_labels;
        std::vector<int> local_to_pool;
        for (int pool_index = 0; pool_index < 30; ++pool_index) {
            if ((anchor_mask & (u32{1} << pool_index)) != 0) continue;
            optional_labels.push_back(POOL[pool_index]);
            local_to_pool.push_back(pool_index);
        }
        require(optional_labels.size() == 27 && local_to_pool.size() == 27,
                "optional vertex set does not have size 27");

        for (int arity : {7}) {
            const Hypergraph hypergraph = build_hypergraph(
                geometry, anchor, local_to_pool, arity);
            const MinimumCoverResult cover =
                exact_minimum_through_eight(hypergraph.edges);
            require(cover.tau <= 8, "minimum cover is not at most eight");
            require(is_cover(cover.cover, hypergraph.edges),
                    "final cover verification failed");

            const std::size_t expected_edges =
                anchor == std::array<int, 3>{80, 143, 252}
                    ? 298279
                    : 286291;
            require(hypergraph.edges.size() == expected_edges &&
                        hypergraph.equality_count == 0 && cover.tau == 8,
                    "THM-4179 depth-seven ledger changed");

            std::cout << "ROW anchor=" << triple_text(anchor)
                      << " original="
                      << (anchor == ORIGINAL_ANCHOR ? "yes" : "no")
                      << " d=" << arity
                      << " candidates=" << hypergraph.candidate_count
                      << " edges=" << hypergraph.edges.size()
                      << " equalities=" << hypergraph.equality_count
                      << " edge_fnv1a64_le=" << hex64(hypergraph.edge_hash)
                      << " min_edge_delta="
                      << hypergraph.minimum_edge_delta
                      << " max_nonedge_delta="
                      << hypergraph.maximum_nonedge_delta
                      << " tau=" << cover.tau
                      << " tau_gt_7=" << (cover.tau > 7 ? "yes" : "no")
                      << " cover="
                      << labels_of_local_mask(cover.cover, optional_labels)
                      << " cover_verified=yes"
                      << " searches_b:Y:N:P:M:D="
                      << search_summary(cover.searches) << "\n";
        }
    }
    require(audited_rows == 2, "did not audit exactly two depth-seven rows");

    const u64 per_anchor = binomial(25, 7);
    const u64 pair_intersection = binomial(24, 6);
    const u64 triple_intersection = binomial(23, 5);
    const u64 union_count =
        3 * per_anchor - 3 * pair_intersection + triple_intersection;
    const u64 alternate_union = binomial(26, 8) - binomial(23, 8);
    require(per_anchor == 480700 && pair_intersection == 134596 &&
                triple_intersection == 33649 && union_count == 1071961 &&
                alternate_union == union_count &&
                union_count - per_anchor == 591261,
            "exact body-union count changed");
    std::cout << "BODY_UNION per_anchor=" << per_anchor
              << " pair_intersection=" << pair_intersection
              << " triple_intersection=" << triple_intersection
              << " distinct_exactly_one_original_anchor=" << union_count
              << " anchor_presentations=" << 3 * per_anchor
              << " increment_beyond_40_anchor=" << union_count - per_anchor
              << " alternate_formula=" << alternate_union << "\n";
    std::cout << "VERDICT PASS\n";
}

}  // namespace

int main() {
    try {
        run_audit();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
