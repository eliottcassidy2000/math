// Independent endpoint-595 audit for the THM-4302 exchanged carrier.
//
// This is deliberately self-contained: it reconstructs the carrier from the
// frozen component ledgers and implements the literal-wall geometry directly.
// Its body replay enumerates combinations recursively, sorts them into numeric
// order, and stops on the first covering coordinate.  It neither includes nor
// calls the maintained audit_pair617 routine.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr u64 kFnvOffset = UINT64_C(0xcbf29ce484222325);
constexpr u64 kFnvPrime = UINT64_C(0x100000001b3);
constexpr u64 kBodyCount = UINT64_C(14307150);
constexpr i64 kFixedGrid = INT64_C(18241159416480);

constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};

constexpr std::array<u32, 14> kPriorDelete = {
    0x00003e1a, 0x000132a3, 0x00017464, 0x00033388,
    0x000a16c2, 0x000f8118, 0x00142a1a, 0x00154348,
    0x00184ba0, 0x001aa260, 0x00202c2b, 0x002066a4,
    0x002b018a, 0x0030c2a2};
constexpr std::array<u32, 14> kPriorAdd = {
    0x18468880, 0x080e8281, 0x22081017, 0x08422a82,
    0x004cac40, 0x19c04044, 0x00c08ec0, 0x10443016,
    0x01609124, 0x10413209, 0x01611640, 0x00606449,
    0x0128d084, 0x08806449};
constexpr std::array<u32, 3> kRepairA = {
    0x2040c641, 0x00508325, 0x002a8641};
constexpr std::array<u32, 3> kRepairB = {
    0x00619324, 0x201813a4, 0x21888126};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

struct Fnv {
    u64 state = kFnvOffset;
    void add(u64 word) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (word >> (8 * byte)) & UINT64_C(0xff);
            state *= kFnvPrime;
        }
    }
};

std::string hex16(u64 value) {
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

u32 parse_mask(const std::string& token) {
    const u64 wide = std::stoull(token, nullptr, 16);
    require(wide < (UINT64_C(1) << 30), "mask escaped 30 labels");
    return static_cast<u32>(wide);
}

u64 fingerprint(const std::vector<u32>& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> read_masks(const std::filesystem::path& path,
                            std::size_t expected_count, u64 expected_fnv,
                            unsigned minimum_rank, unsigned maximum_rank) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask ledger " + path.string());
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask(token);
        const unsigned rank = std::popcount(mask);
        require(rank >= minimum_rank && rank <= maximum_rank,
                "mask rank changed in " + path.string());
        require(distinct.insert(mask).second,
                "duplicate mask in " + path.string());
        masks.push_back(mask);
    }
    require(masks.size() == expected_count,
            "mask count changed in " + path.string());
    require(fingerprint(masks) == expected_fnv,
            "mask FNV changed in " + path.string());
    return masks;
}

void append_distinct(std::vector<u32>& target, std::set<u32>& distinct,
                     const std::vector<u32>& additions) {
    for (u32 mask : additions) {
        require(distinct.insert(mask).second, "carrier addition overlaps");
        target.push_back(mask);
    }
}

template <std::size_t N>
void append_distinct(std::vector<u32>& target, std::set<u32>& distinct,
                     const std::array<u32, N>& additions) {
    append_distinct(target, distinct,
                    std::vector<u32>(additions.begin(), additions.end()));
}

std::vector<u32> reconstruct_carrier(
    const std::filesystem::path& base_path,
    const std::filesystem::path& add45_path,
    const std::filesystem::path& suffix_path,
    const std::filesystem::path& repairs76_path,
    const std::filesystem::path& delete73_path,
    const std::filesystem::path& additions4_path) {
    std::vector<u32> carrier = read_masks(
        base_path, 8951, UINT64_C(0x188f82ab9dd1695a), 8, 8);
    std::set<u32> distinct(carrier.begin(), carrier.end());
    append_distinct(carrier, distinct, read_masks(
        add45_path, 45, UINT64_C(0xec083b65cc8c34e3), 8, 8));
    append_distinct(carrier, distinct, std::vector<u32>{0x014c9084});
    append_distinct(carrier, distinct, read_masks(
        suffix_path, 9, UINT64_C(0x02b936529030e4bc), 8, 8));

    const std::set<u32> prior_delete(kPriorDelete.begin(), kPriorDelete.end());
    std::vector<u32> post_prior_exchange;
    for (u32 mask : carrier)
        if (!prior_delete.contains(mask)) post_prior_exchange.push_back(mask);
    require(post_prior_exchange.size() + prior_delete.size() == carrier.size(),
            "prior deletions not all present");
    distinct = std::set<u32>(post_prior_exchange.begin(),
                             post_prior_exchange.end());
    append_distinct(post_prior_exchange, distinct, kPriorAdd);
    append_distinct(post_prior_exchange, distinct, kRepairA);
    append_distinct(post_prior_exchange, distinct, kRepairB);
    require(post_prior_exchange.size() == 9012,
            "pre-THM4300 carrier count changed");

    append_distinct(post_prior_exchange, distinct, read_masks(
        repairs76_path, 76, UINT64_C(0x64ce5f9d1ec8c4c2), 8, 9));
    require(post_prior_exchange.size() == 9088 &&
                fingerprint(post_prior_exchange) ==
                    UINT64_C(0x55e8588798885ae5),
            "augmented THM4300 carrier changed");

    const std::vector<u32> delete73 = read_masks(
        delete73_path, 73, UINT64_C(0x9240b264ab65aa62), 8, 8);
    const std::set<u32> deletion_set(delete73.begin(), delete73.end());
    std::vector<u32> exchanged;
    for (u32 mask : post_prior_exchange)
        if (!deletion_set.contains(mask)) exchanged.push_back(mask);
    require(exchanged.size() + deletion_set.size() == post_prior_exchange.size(),
            "THM4302 deletions not all present");
    distinct = std::set<u32>(exchanged.begin(), exchanged.end());
    append_distinct(exchanged, distinct, read_masks(
        additions4_path, 4, UINT64_C(0xdc0eebaebf688c65), 9, 9));
    require(exchanged.size() == 9019 &&
                fingerprint(exchanged) == UINT64_C(0x892fef44a9e6b37e),
            "THM4302 exchanged carrier changed");
    require(std::count_if(exchanged.begin(), exchanged.end(), [](u32 mask) {
                return std::popcount(mask) == 8;
            }) == 8961,
            "exchanged rank-eight count changed");
    return exchanged;
}

struct Pair {
    int q = 0;
    int r = 0;
    auto operator<=>(const Pair&) const = default;
};

std::vector<Pair> read_pairs(const std::filesystem::path& path,
                             std::size_t expected_count, u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pair ledger " + path.string());
    std::vector<Pair> rows;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos, "malformed pair ledger");
        const Pair pair{std::stoi(line.substr(0, comma)),
                        std::stoi(line.substr(comma + 1))};
        require(pair.q > 0 && pair.q < pair.r, "invalid pair");
        require(rows.empty() || rows.back() < pair,
                "pair ledger not strictly ordered");
        rows.push_back(pair);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(rows.size() == expected_count && ledger.state == expected_fnv,
            "pair ledger identity changed");
    return rows;
}

u64 pair_fingerprint(const std::set<Pair>& rows) {
    Fnv ledger;
    for (Pair pair : rows) {
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    return ledger.state;
}

std::vector<Pair> derive_top_rows(const std::filesystem::path& universe_path,
                                  const std::filesystem::path& old_union_path,
                                  std::set<Pair>& universe,
                                  std::set<Pair>& typed_union) {
    const auto universe_vector = read_pairs(
        universe_path, 22647, UINT64_C(0xdf5374d4aca67677));
    const auto old_union_vector = read_pairs(
        old_union_path, 1624, UINT64_C(0x11414a33ab91fef6));
    universe = {universe_vector.begin(), universe_vector.end()};
    const std::set<Pair> old_union(old_union_vector.begin(), old_union_vector.end());
    require(std::includes(universe.begin(), universe.end(),
                          old_union.begin(), old_union.end()),
            "old typed union escaped universe");
    typed_union = old_union;
    for (Pair pair : universe)
        if (pair.r >= 596) typed_union.insert(pair);
    require(typed_union.size() == 1633 &&
                pair_fingerprint(typed_union) == UINT64_C(0xb1c8ecf1dd4a71c5),
            "THM4302 typed union changed");
    std::vector<Pair> residual;
    std::set_difference(universe.begin(), universe.end(),
                        typed_union.begin(), typed_union.end(),
                        std::back_inserter(residual));
    require(residual.size() == 21014, "THM4302 residual count changed");
    const int maximum = std::max_element(
        residual.begin(), residual.end(),
        [](Pair left, Pair right) { return left.r < right.r; })->r;
    std::vector<Pair> top;
    std::copy_if(residual.begin(), residual.end(), std::back_inserter(top),
                 [=](Pair pair) { return pair.r == maximum; });
    const std::set<Pair> top_set(top.begin(), top.end());
    require(maximum == 595 && top.size() == 28 &&
                pair_fingerprint(top_set) == UINT64_C(0x47981ce64825ef2a),
            "derived residual boundary changed");
    return top;
}

i64 checked_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive LCM input");
    const i128 wide = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(wide <= std::numeric_limits<i64>::max(), "grid overflow");
    return static_cast<i64>(wide);
}

i64 fixed_grid() {
    i64 grid = 1;
    for (int speed : kPool) grid = checked_lcm(grid, 14LL * speed);
    require(grid == kFixedGrid, "fixed grid changed");
    return grid;
}

bool safe_midpoint(int speed, i64 grid, i64 left, i64 right) {
    i128 residue = static_cast<i128>(speed) *
                   (static_cast<i128>(left) + right);
    residue %= static_cast<i128>(2) * grid;
    if (residue < 0) residue += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * residue >= grid &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * grid;
}

struct Geometry {
    i64 grid = 0;
    std::vector<std::pair<u32, i64>> classes;
};

Geometry build_geometry(Pair pair) {
    i64 grid = fixed_grid();
    grid = checked_lcm(grid, 14LL * pair.q);
    grid = checked_lcm(grid, 14LL * pair.r);
    std::vector<i64> walls = {0, grid};
    auto add_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        require(grid % divisor == 0, "nonintegral wall unit");
        const i64 unit = grid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(pair.q);
    add_walls(pair.r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::map<u32, i64> by_failure;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(pair.q, grid, left, right) ||
            !safe_midpoint(pair.r, grid, left, right)) continue;
        u32 failure = 0;
        for (unsigned label = 0; label < kPool.size(); ++label)
            if (!safe_midpoint(kPool[label], grid, left, right))
                failure |= u32{1} << label;
        by_failure[failure] += right - left;
    }
    Geometry geometry;
    geometry.grid = grid;
    for (const auto& [failure, width] : by_failure)
        if (std::popcount(failure) <= 9)
            geometry.classes.emplace_back(failure, width);
    return geometry;
}

bool active(const Geometry& geometry, u32 mask) {
    i64 mass = 0;
    for (const auto& [failure, width] : geometry.classes)
        if ((failure & ~mask) == 0) mass += width;
    return static_cast<i128>(63) * mass -
               static_cast<i128>(4) * geometry.grid >= 0;
}

void generate_bodies_rec(unsigned start, unsigned needed, u32 mask,
                         std::vector<u32>& bodies) {
    if (needed == 0) {
        bodies.push_back(mask);
        return;
    }
    for (unsigned bit = start; bit + needed <= 30; ++bit)
        generate_bodies_rec(bit + 1, needed - 1,
                            mask | (u32{1} << bit), bodies);
}

std::vector<u32> generate_bodies() {
    std::vector<u32> bodies;
    bodies.reserve(kBodyCount);
    generate_bodies_rec(0, 9, 0, bodies);
    require(bodies.size() == kBodyCount, "body count changed");
    std::sort(bodies.begin(), bodies.end());
    require(std::adjacent_find(bodies.begin(), bodies.end()) == bodies.end(),
            "body enumeration duplicated");
    require(std::all_of(bodies.begin(), bodies.end(), [](u32 body) {
                return std::popcount(body) == 9 && body < (u32{1} << 30);
            }), "body enumeration escaped universe");
    return bodies;
}

struct Audit {
    Pair pair;
    u64 active_count = 0;
    u64 active_fnv = 0;
    u64 active_joint = 0;
    u64 active_nonjoint = 0;
    u64 exposed = 0;
    u64 exposed_fnv = 0;
    u64 failures = 0;
    u64 failure_fnv = 0;
};

Audit audit_pair(Pair pair, const std::vector<u32>& bodies,
                 const std::vector<u32>& joint,
                 const std::unordered_set<u32>& joint_set,
                 const std::vector<u32>& carrier) {
    const Geometry geometry = build_geometry(pair);
    std::vector<u32> active_joint;
    std::vector<u32> active_nonjoint;
    Fnv active_ledger;
    Audit audit;
    audit.pair = pair;
    for (u32 mask : joint)
        if (active(geometry, mask)) active_joint.push_back(mask);
    for (u32 mask : carrier) {
        if (!active(geometry, mask)) continue;
        ++audit.active_count;
        active_ledger.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    audit.active_fnv = active_ledger.state;
    audit.active_joint = active_joint.size();
    audit.active_nonjoint = active_nonjoint.size();
    require(audit.active_count == audit.active_joint + audit.active_nonjoint,
            "active partition changed");

    Fnv exposed_ledger;
    Fnv failure_ledger;
    for (u32 body : bodies) {
        bool covered = false;
        for (u32 mask : active_joint) {
            if ((mask & body) == 0) {
                covered = true;
                break;
            }
        }
        if (covered) continue;
        ++audit.exposed;
        exposed_ledger.add(body);
        for (u32 mask : active_nonjoint) {
            if ((mask & body) == 0) {
                covered = true;
                break;
            }
        }
        if (!covered) {
            ++audit.failures;
            failure_ledger.add(body);
        }
    }
    audit.exposed_fnv = exposed_ledger.state;
    audit.failure_fnv = failure_ledger.state;
    return audit;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 11,
                "usage: independent595 JOINT BASE8951 ADD45 SUFFIX9 "
                "REPAIRS76 DELETE73 ADDITIONS4 UNIVERSE OLD_TYPED_UNION "
                "WORKERS");
        const std::vector<u32> joint = read_masks(
            argv[1], 421, UINT64_C(0x20d63dd42fe8150e), 8, 8);
        const std::vector<u32> carrier = reconstruct_carrier(
            argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        for (u32 mask : joint)
            require(std::find(carrier.begin(), carrier.end(), mask) != carrier.end(),
                    "joint coordinate absent from exchanged carrier");

        std::set<Pair> universe;
        std::set<Pair> typed_union;
        const std::vector<Pair> top = derive_top_rows(
            argv[8], argv[9], universe, typed_union);
        const std::vector<u32> bodies = generate_bodies();

        std::vector<Audit> audits(top.size());
        std::atomic<std::size_t> next{0};
        const unsigned workers = std::stoul(argv[10]);
        require(workers >= 1 && workers <= top.size(),
                "worker count escaped row count");
        std::vector<std::thread> pool;
        for (unsigned worker = 0; worker < workers; ++worker) {
            pool.emplace_back([&] {
                while (true) {
                    const std::size_t index = next.fetch_add(1);
                    if (index >= top.size()) return;
                    audits[index] = audit_pair(top[index], bodies, joint,
                                               joint_set, carrier);
                }
            });
        }
        for (std::thread& worker : pool) worker.join();

        std::set<Pair> closed;
        std::set<Pair> surviving_top;
        const std::map<Pair, std::pair<u64, u64>> expected_failures = {
            {{96, 595}, {116, UINT64_C(0xfedacdbff3f31981)}},
            {{100, 595}, {13, UINT64_C(0x3ac9ac8b4b9ad93f)}},
            {{210, 595}, {16, UINT64_C(0xa6a226f12c168d3a)}}};
        std::cout << "THM4303_ENDPOINT595_INDEPENDENT_REPLAY_V1\n"
                  << "CARRIER " << carrier.size() << " FNV "
                  << hex16(fingerprint(carrier)) << " RANK8 8961 RANK9 58\n"
                  << "ROWS_DERIVED_FROM_TYPED_RESIDUAL " << top.size()
                  << " FNV 47981ce64825ef2a BODY_UNIVERSE_PER_ROW "
                  << bodies.size() << " WORKERS " << workers << '\n';
        for (const Audit& audit : audits) {
            std::cout << "PAIR " << audit.pair.q << ',' << audit.pair.r
                      << " ACTIVE " << audit.active_count
                      << " ACTIVE_FNV " << hex16(audit.active_fnv)
                      << " JOINT " << audit.active_joint
                      << " NONJOINT " << audit.active_nonjoint
                      << " EXPOSED " << audit.exposed
                      << " EXPOSED_FNV " << hex16(audit.exposed_fnv)
                      << " FAILURES " << audit.failures
                      << " FAILURE_FNV " << hex16(audit.failure_fnv) << '\n';
            const auto expected = expected_failures.find(audit.pair);
            if (expected == expected_failures.end()) {
                require(audit.failures == 0,
                        "unexpected independent failing row");
                closed.insert(audit.pair);
            } else {
                require(audit.failures == expected->second.first &&
                            audit.failure_fnv == expected->second.second,
                        "independent failure identity changed");
                surviving_top.insert(audit.pair);
            }
        }
        require(closed.size() == 25 && surviving_top.size() == 3,
                "unexpected endpoint-595 boundary split");
        require(pair_fingerprint(closed) == UINT64_C(0x1ad4f3c2ab6ea09d) &&
                    pair_fingerprint(surviving_top) ==
                        UINT64_C(0x9853e590efc73022),
                "independent boundary identities changed");
        std::set<Pair> enlarged_union = typed_union;
        enlarged_union.insert(closed.begin(), closed.end());
        std::set<Pair> enlarged_residual;
        std::set_difference(universe.begin(), universe.end(),
                            enlarged_union.begin(), enlarged_union.end(),
                            std::inserter(enlarged_residual,
                                          enlarged_residual.end()));
        require(enlarged_union.size() == 1658 &&
                    pair_fingerprint(enlarged_union) ==
                        UINT64_C(0x43317f1aee06e8bd) &&
                    enlarged_residual.size() == 20989 &&
                    pair_fingerprint(enlarged_residual) ==
                        UINT64_C(0xb0fbaa28440a118f),
                "enlarged typed partition count changed");
        std::cout << "BOUNDARY_CLOSED " << closed.size()
                  << " FNV " << hex16(pair_fingerprint(closed))
                  << " BOUNDARY_SURVIVORS " << surviving_top.size()
                  << " FNV " << hex16(pair_fingerprint(surviving_top)) << '\n'
                  << "SURVIVING_TOP_ROWS";
        for (Pair pair : surviving_top)
            std::cout << ' ' << '(' << pair.q << ',' << pair.r << ')';
        std::cout << '\n'
                  << "ENLARGED_TYPED_UNION " << enlarged_union.size()
                  << " FNV " << hex16(pair_fingerprint(enlarged_union))
                  << " ENLARGED_RESIDUAL " << enlarged_residual.size()
                  << " FNV " << hex16(pair_fingerprint(enlarged_residual)) << '\n'
                  << "SCOPE FINITE_EXACT_FIXED_CARRIER_ROW_CONSEQUENCES_ONLY_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "INDEPENDENT595_ERROR " << error.what() << '\n';
        return 1;
    }
}
