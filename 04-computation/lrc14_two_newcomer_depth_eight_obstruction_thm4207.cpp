// Primary exact joint-wall audit for THM-4207.
//
// This path builds the literal wall arrangement of P union {50,51}, records
// every labelled pool-failure atom, enumerates the complete deletion decks,
// and exhausts every seven-set when proving tau(E_8(50,51))=8.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using i64 = std::int64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

i64 checked_lcm(i64 a, i64 b) {
    const i64 g = std::gcd(a, b);
    const i128 value = static_cast<i128>(a / g) * b;
    require(value <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(value);
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    if (negative) value = -value;
    std::string out;
    while (value != 0) {
        out.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}

bool safe_midpoint(int speed, i64 left, i64 right, i64 denominator) {
    const i128 period = static_cast<i128>(2) * denominator;
    const i128 raw = static_cast<i128>(speed) * (left + right);
    const i128 residue = raw % period;
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * denominator;
}

struct ExactMass {
    i64 denominator;
    i64 numerator;
    std::size_t wall_count;
};

struct PairAtom {
    i64 length;
    std::uint32_t failed_pool;
};

struct PairGeometry {
    i64 denominator;
    std::vector<PairAtom> atoms;
};

PairGeometry pair_geometry(int q, int r) {
    i64 denominator = 1;
    for (int speed : POOL) denominator = checked_lcm(denominator, 14LL * speed);
    denominator = checked_lcm(denominator, 14LL * q);
    denominator = checked_lcm(denominator, 14LL * r);

    std::vector<int> all(POOL.begin(), POOL.end());
    all.push_back(q);
    all.push_back(r);
    std::vector<i64> walls = {0, denominator};
    for (int speed : all) {
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    std::vector<PairAtom> atoms;
    for (std::size_t i = 0; i + 1 < walls.size(); ++i) {
        if (!safe_midpoint(q, walls[i], walls[i + 1], denominator) ||
            !safe_midpoint(r, walls[i], walls[i + 1], denominator)) continue;
        std::uint32_t failed = 0;
        for (std::size_t vertex = 0; vertex < POOL.size(); ++vertex) {
            if (!safe_midpoint(POOL[vertex], walls[i], walls[i + 1], denominator)) {
                failed |= std::uint32_t{1} << vertex;
            }
        }
        atoms.push_back({walls[i + 1] - walls[i], failed});
    }
    return {denominator, std::move(atoms)};
}

i64 body_mass(std::uint32_t body, const PairGeometry& geometry) {
    i64 mass = 0;
    for (const PairAtom& atom : geometry.atoms) {
        if ((atom.failed_pool & body) == 0) mass += atom.length;
    }
    return mass;
}

std::string body_labels(std::uint32_t body) {
    std::ostringstream out;
    bool first = true;
    for (std::size_t vertex = 0; vertex < POOL.size(); ++vertex) {
        if ((body & (std::uint32_t{1} << vertex)) == 0) continue;
        if (!first) out << ',';
        first = false;
        out << POOL[vertex];
    }
    return out.str();
}

std::uint32_t mask_of_values(std::initializer_list<int> values) {
    std::uint32_t result = 0;
    for (int value : values) {
        const auto found = std::find(POOL.begin(), POOL.end(), value);
        require(found != POOL.end(), "mask label absent from pool");
        result |= std::uint32_t{1} << std::distance(POOL.begin(), found);
    }
    return result;
}

void optimize_body(int q, int r) {
    const PairGeometry geometry = pair_geometry(q, r);
    std::uint32_t body = 0;
    for (int chosen = 0; chosen < 10; ++chosen) {
        i64 best_mass = std::numeric_limits<i64>::max();
        int best_vertex = -1;
        for (int vertex = 0; vertex < 30; ++vertex) {
            const std::uint32_t bit = std::uint32_t{1} << vertex;
            if ((body & bit) != 0) continue;
            const i64 mass = body_mass(body | bit, geometry);
            if (mass < best_mass) {
                best_mass = mass;
                best_vertex = vertex;
            }
        }
        require(best_vertex >= 0, "greedy selection failed");
        body |= std::uint32_t{1} << best_vertex;
    }

    bool improved = true;
    while (improved) {
        improved = false;
        i64 best_mass = body_mass(body, geometry);
        std::uint32_t best_body = body;
        for (int remove = 0; remove < 30; ++remove) {
            const std::uint32_t remove_bit = std::uint32_t{1} << remove;
            if ((body & remove_bit) == 0) continue;
            for (int add = 0; add < 30; ++add) {
                const std::uint32_t add_bit = std::uint32_t{1} << add;
                if ((body & add_bit) != 0) continue;
                const std::uint32_t candidate = (body ^ remove_bit) | add_bit;
                const i64 mass = body_mass(candidate, geometry);
                if (mass < best_mass) {
                    best_mass = mass;
                    best_body = candidate;
                }
            }
        }
        if (best_body != body) {
            body = best_body;
            improved = true;
        }
    }
    const i64 numerator = body_mass(body, geometry);
    const i128 delta = static_cast<i128>(63) * numerator -
                       static_cast<i128>(4) * geometry.denominator;
    std::cout << "OPT_PAIR q=" << q << " r=" << r
              << " atoms=" << geometry.atoms.size()
              << " denominator=" << geometry.denominator
              << " numerator=" << numerator
              << " delta63=" << decimal(delta)
              << " body={" << body_labels(body) << '}'
              << " verdict=" << (delta >= 0 ? "SAFE" : "FAIL") << '\n';
}

std::uint32_t next_combination(std::uint32_t mask) {
    const std::uint32_t low = mask & (~mask + 1u);
    const std::uint32_t ripple = mask + low;
    if (ripple == 0) return 0;
    return ripple | (((mask ^ ripple) >> 2) / low);
}

void pair_layer_counts(int q, int r, int max_arity) {
    const PairGeometry geometry = pair_geometry(q, r);
    require(q == 50 && r == 51 && max_arity == 8,
            "unexpected primary audit universe");
    require(geometry.denominator == INT64_C(91205797082400),
            "pair denominator changed");
    constexpr std::array<std::uint64_t, 9> expected_edges = {
        0, 0, 0, 0, 0, 0, 39, 10114, 311544};
    std::unordered_map<std::uint32_t, i64> atom_mass;
    for (const PairAtom& atom : geometry.atoms) atom_mass[atom.failed_pool] += atom.length;
    for (int arity = 0; arity <= max_arity; ++arity) {
        std::uint64_t edges = 0;
        std::uint64_t equalities = 0;
        i128 best_delta = -static_cast<i128>(geometry.denominator) * 4;
        std::uint32_t best_edge = 0;
        std::uint32_t deletion = arity == 0 ? 0 :
            (std::uint32_t{1} << arity) - 1;
        const std::uint32_t limit = std::uint32_t{1} << 30;
        while (true) {
            i64 mass = 0;
            std::uint32_t subset = deletion;
            while (true) {
                const auto found = atom_mass.find(subset);
                if (found != atom_mass.end()) mass += found->second;
                if (subset == 0) break;
                subset = (subset - 1) & deletion;
            }
            const i128 delta = static_cast<i128>(63) * mass -
                               static_cast<i128>(4) * geometry.denominator;
            if (delta > best_delta) {
                best_delta = delta;
                best_edge = deletion;
            }
            if (delta == 0) ++equalities;
            if (delta >= 0) ++edges;
            if (arity == 0) break;
            deletion = next_combination(deletion);
            if (deletion == 0 || deletion >= limit) break;
        }
        require(edges == expected_edges[arity], "pair edge count changed");
        require(equalities == 0, "pair threshold equality appeared");
        std::cout << "PAIR_LAYER q=" << q << " r=" << r
                  << " d=" << arity << " edges=" << edges
                  << " equalities=" << equalities
                  << " best_delta63=" << decimal(best_delta)
                  << " best={" << body_labels(best_edge) << "}" << '\n';
    }
}

std::vector<std::uint32_t> pair_edges(int q, int r, int arity) {
    const PairGeometry geometry = pair_geometry(q, r);
    std::unordered_map<std::uint32_t, i64> atom_mass;
    for (const PairAtom& atom : geometry.atoms) atom_mass[atom.failed_pool] += atom.length;
    std::vector<std::uint32_t> edges;
    std::uint32_t deletion = (std::uint32_t{1} << arity) - 1;
    const std::uint32_t limit = std::uint32_t{1} << 30;
    while (deletion != 0 && deletion < limit) {
        i64 mass = 0;
        std::uint32_t subset = deletion;
        while (true) {
            const auto found = atom_mass.find(subset);
            if (found != atom_mass.end()) mass += found->second;
            if (subset == 0) break;
            subset = (subset - 1) & deletion;
        }
        if (static_cast<i128>(63) * mass >= static_cast<i128>(4) * geometry.denominator) {
            edges.push_back(deletion);
        }
        deletion = next_combination(deletion);
    }
    return edges;
}

std::uint32_t greedy_cover(const std::vector<std::uint32_t>& edges) {
    std::uint32_t chosen = 0;
    while (true) {
        std::array<std::uint64_t, 30> scores{};
        std::uint64_t uncovered = 0;
        for (std::uint32_t edge : edges) {
            if ((edge & chosen) != 0) continue;
            ++uncovered;
            for (int vertex = 0; vertex < 30; ++vertex) {
                if ((edge & (std::uint32_t{1} << vertex)) != 0) ++scores[vertex];
            }
        }
        if (uncovered == 0) return chosen;
        int best = -1;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if ((chosen & (std::uint32_t{1} << vertex)) != 0) continue;
            if (best < 0 || scores[vertex] > scores[best]) best = vertex;
        }
        require(best >= 0 && scores[best] > 0, "greedy cover stalled");
        chosen |= std::uint32_t{1} << best;
    }
}

bool is_cover(std::uint32_t candidate, const std::vector<std::uint32_t>& edges) {
    for (std::uint32_t edge : edges) {
        if ((candidate & edge) == 0) return false;
    }
    return true;
}

void cover_probe(int q, int r, int arity) {
    const std::vector<std::uint32_t> edges = pair_edges(q, r, arity);
    const std::uint32_t greedy = greedy_cover(edges);
    std::cout << "PAIR_GREEDY q=" << q << " r=" << r << " d=" << arity
              << " edges=" << edges.size()
              << " cover_size=" << __builtin_popcount(greedy)
              << " cover={" << body_labels(greedy) << "}"
              << " verified=" << (is_cover(greedy, edges) ? "YES" : "NO") << '\n';

    const int upper = __builtin_popcount(greedy);
    std::uint64_t candidates = 0;
    std::uint64_t edge_checks = 0;
    bool exact_found = false;
    for (int size = 0; size <= upper; ++size) {
        std::uint32_t candidate = size == 0 ? 0 :
            (std::uint32_t{1} << size) - 1;
        const std::uint32_t limit = std::uint32_t{1} << 30;
        bool found = false;
        std::uint32_t witness = 0;
        while (true) {
            ++candidates;
            bool cover = true;
            for (std::uint32_t edge : edges) {
                ++edge_checks;
                if ((candidate & edge) == 0) {
                    cover = false;
                    break;
                }
            }
            if (cover) {
                found = true;
                witness = candidate;
                break;
            }
            if (size == 0) break;
            candidate = next_combination(candidate);
            if (candidate == 0 || candidate >= limit) break;
        }
        if (found) {
            const int expected = arity == 6 ? 2 : 5;
            require(size == expected, "shallow transversal number changed");
            std::cout << "PAIR_EXACT_TAU q=" << q << " r=" << r
                      << " d=" << arity << " tau=" << size
                      << " witness={" << body_labels(witness) << "}"
                      << " candidates=" << candidates
                      << " edge_checks=" << edge_checks << '\n';
            exact_found = true;
            break;
        }
    }
    require(exact_found, "shallow transversal witness not found");
}

std::uint64_t edge_key(std::uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

void target_seven_probe(int q, int r, int arity) {
    std::vector<std::uint32_t> edges = pair_edges(q, r, arity);
    require(q == 50 && r == 51 && arity == 8,
            "unexpected depth-eight target universe");
    require(edges.size() == 311544, "depth-eight edge count changed");
    const std::uint32_t exact_cover =
        mask_of_values({16,88,95,143,168,193,240,290});
    require(is_cover(exact_cover, edges), "exact eight-cover changed");
    const std::uint32_t greedy = greedy_cover(edges);
    std::sort(edges.begin(), edges.end(), [](std::uint32_t a, std::uint32_t b) {
        const std::uint64_t ka = edge_key(a);
        const std::uint64_t kb = edge_key(b);
        if (ka != kb) return ka < kb;
        return a < b;
    });

    std::uint64_t candidates = 0;
    std::uint64_t checks = 0;
    std::uint64_t max_checks = 0;
    std::uint32_t closest = 0;
    std::uint32_t closest_miss = 0;
    std::uint32_t cover = 0;
    std::uint32_t body = (std::uint32_t{1} << 7) - 1;
    const std::uint32_t limit = std::uint32_t{1} << 30;
    while (body < limit) {
        std::uint64_t row_checks = 0;
        std::uint32_t missed = 0;
        for (std::uint32_t edge : edges) {
            ++row_checks;
            if ((body & edge) == 0) {
                missed = edge;
                break;
            }
        }
        ++candidates;
        checks += row_checks;
        if (row_checks > max_checks) {
            max_checks = row_checks;
            closest = body;
            closest_miss = missed;
        }
        if (missed == 0) {
            cover = body;
            break;
        }
        const std::uint32_t next_body = next_combination(body);
        if (next_body <= body) break;
        body = next_body;
    }
    require(candidates == UINT64_C(2035800),
            "seven-candidate universe changed");
    require(cover == 0, "unexpected depth-eight seven-cover");
    std::cout << "TARGET7 q=" << q << " r=" << r << " d=" << arity
              << " edges=" << edges.size()
              << " greedy_size=" << __builtin_popcount(greedy)
              << " greedy={" << body_labels(greedy) << "}"
              << " greedy_verified=" << (is_cover(greedy, edges) ? "YES" : "NO")
              << " candidates=" << candidates
              << " expected=" << UINT64_C(2035800)
              << " checks=" << checks
              << " max_checks=" << max_checks
              << " closest={" << body_labels(closest) << "}"
              << " missed={" << body_labels(closest_miss) << "}"
              << " cover=" << (cover == 0 ? "NONE" : body_labels(cover))
              << " verdict=" << (cover == 0 ? "TAU_GT_7" : "TAU_LE_7")
              << '\n';
}

ExactMass exact_mass(std::vector<int> speeds) {
    std::sort(speeds.begin(), speeds.end());
    speeds.erase(std::unique(speeds.begin(), speeds.end()), speeds.end());
    i64 denominator = 1;
    for (int speed : speeds) denominator = checked_lcm(denominator, 14LL * speed);

    std::vector<i64> walls = {0, denominator};
    for (int speed : speeds) {
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    i64 numerator = 0;
    for (std::size_t i = 0; i + 1 < walls.size(); ++i) {
        bool safe = true;
        for (int speed : speeds) {
            if (!safe_midpoint(speed, walls[i], walls[i + 1], denominator)) {
                safe = false;
                break;
            }
        }
        if (safe) numerator += walls[i + 1] - walls[i];
    }
    return {denominator, numerator, walls.size()};
}

std::vector<int> pool_minus(const std::vector<int>& deletion) {
    std::vector<int> out;
    for (int value : POOL) {
        if (std::find(deletion.begin(), deletion.end(), value) == deletion.end()) {
            out.push_back(value);
        }
    }
    return out;
}

i128 report(const std::string& name, const std::vector<int>& speeds) {
    const ExactMass mass = exact_mass(speeds);
    const i128 delta = static_cast<i128>(63) * mass.numerator -
                       static_cast<i128>(4) * mass.denominator;
    std::cout << name << " size=" << speeds.size()
              << " denominator=" << mass.denominator
              << " numerator=" << mass.numerator
              << " walls=" << mass.wall_count
              << " delta63=" << decimal(delta)
              << " verdict=" << (delta >= 0 ? "SAFE" : "FAIL") << '\n';
    return delta;
}

}  // namespace

int main() {
    try {
        const std::vector<int> body = {80,85,88,95,143,145,168,193,240,252};
        const std::vector<int> cover_body =
            {10,80,95,120,132,143,145,168,193,240};
        const std::vector<int> d8_cover_body9 =
            {8,16,88,95,143,168,193,240,290};
        const std::vector<int> repair = {8,16,42,63,132,170,290};
        std::vector<int> ambient = pool_minus(repair);
        const std::vector<int> repair5 = {8,10,95,193,240};
        std::vector<int> ambient5 = pool_minus(repair5);
        const std::vector<int> repair4 = {88,95,176,193};
        std::vector<int> ambient4 = pool_minus(repair4);
        const std::vector<int> d8_missed_repair =
            {8,42,88,132,145,170,176,264};
        std::vector<int> ambient8 = pool_minus(d8_missed_repair);

        auto with = [](std::vector<int> base, std::initializer_list<int> extra) {
            base.insert(base.end(), extra.begin(), extra.end());
            return base;
        };

        report("BODY", body);
        report("BODY_Q50", with(body, {50}));
        report("BODY_Q51", with(body, {51}));
        report("BODY_Q50_R51", with(body, {50,51}));
        report("COVER_BODY_Q50", with(cover_body, {50}));
        report("COVER_BODY_Q51", with(cover_body, {51}));
        report("COVER_BODY_Q50_R51", with(cover_body, {50,51}));
        const i128 cover_body_delta =
            report("D8_COVER_BODY9_Q50_R51", with(d8_cover_body9, {50,51}));
        require(cover_body_delta == static_cast<i128>(204816647179284LL),
                "depth-eight cover-body control changed");
        report("AMBIENT", ambient);
        report("AMBIENT_Q50", with(ambient, {50}));
        report("AMBIENT_Q51", with(ambient, {51}));
        report("AMBIENT_Q50_R51", with(ambient, {50,51}));
        report("AMBIENT5_Q50", with(ambient5, {50}));
        report("AMBIENT5_Q51", with(ambient5, {51}));
        report("AMBIENT5_Q50_R51", with(ambient5, {50,51}));
        report("AMBIENT4_Q50", with(ambient4, {50}));
        report("AMBIENT4", ambient4);
        report("AMBIENT4_Q51", with(ambient4, {51}));
        report("AMBIENT4_Q50_R51", with(ambient4, {50,51}));
        const i128 missed_repair_delta =
            report("D8_MISSED_REPAIR_Q50_R51", with(ambient8, {50,51}));
        require(missed_repair_delta == static_cast<i128>(1334721427452LL),
                "depth-eight missed-repair control changed");
        optimize_body(50, 51);
        pair_layer_counts(50, 51, 8);
        cover_probe(50, 51, 6);
        cover_probe(50, 51, 7);
        target_seven_probe(50, 51, 8);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
