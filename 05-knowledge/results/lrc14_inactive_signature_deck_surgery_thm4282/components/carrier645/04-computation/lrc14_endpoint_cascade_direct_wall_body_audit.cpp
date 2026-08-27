// Clean-room exact auditor for the post-THM-4252 endpoint cascade.
//
// Input is only the directory of one-pair primary probe transcripts.  This
// auditor rebuilds a literal joint wall arrangement for the thirty fixed pool
// speeds and each named endpoint pair, integrates every emitted prefix repair
// by summing joint-cell widths, and recursively generates every labelled
// nine-body.  It does not evaluate the primitive endpoint cocycle, group pool
// cells into atoms, rank repairs in colex order, or run a superset transform.
//
// All consequence-bearing gates use require(), including under NDEBUG.  No
// floating point, sampling, or C/C++ assert participates.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};

constexpr u64 EXPECTED_BODIES = UINT64_C(14307150);
constexpr u64 EXPECTED_BAND_FNV = UINT64_C(0xb3d54b78babbcaec);

constexpr std::array<std::pair<int, int>, 59> EXPECTED_PAIRS = {{
    {616,755}, {616,756}, {616,757}, {616,758}, {616,759},
    {616,760}, {616,761}, {616,762}, {616,763}, {616,764},
    {616,765}, {616,766}, {616,767}, {616,768}, {698,755},
    {698,757}, {704,755}, {704,757}, {704,758}, {704,759},
    {704,761}, {704,762}, {704,763}, {704,764}, {704,765},
    {721,755}, {721,757}, {721,758}, {721,759}, {721,761},
    {721,762}, {721,763}, {721,764}, {721,765}, {721,766},
    {721,767}, {721,768}, {726,755}, {726,757}, {726,758},
    {726,761}, {732,755}, {732,757}, {732,761}, {732,762},
    {732,763}, {744,762}, {744,763}, {744,765}, {744,766},
    {744,768}, {750,762}, {750,763}, {750,765}, {750,766},
    {750,768}, {765,766}, {765,768}, {766,768}
}};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    if (negative) value = -value;
    std::string answer;
    while (value != 0) {
        answer.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

i128 gcd128(i128 left, i128 right) {
    if (left < 0) left = -left;
    if (right < 0) right = -right;
    while (right != 0) {
        const i128 remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

std::string fraction(i128 numerator, i128 denominator) {
    require(denominator > 0, "nonpositive fraction denominator");
    const i128 divisor = gcd128(numerator, denominator);
    return decimal(numerator / divisor) + "/" + decimal(denominator / divisor);
}

i64 checked_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive LCM input");
    const i128 value = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(value <= std::numeric_limits<i64>::max(), "joint grid overflow");
    return static_cast<i64>(value);
}

bool safe_midpoint(int speed, i64 grid, i64 left, i64 right) {
    // Cast before addition: for the (698,757) band edge, 2*grid exceeds
    // INT64_MAX even though every individual wall coordinate fits in i64.
    i128 residue = static_cast<i128>(speed) *
                   (static_cast<i128>(left) + right);
    residue %= static_cast<i128>(2) * grid;
    if (residue < 0) residue += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * residue >= grid &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * grid;
}

struct Cell {
    i64 width = 0;
    u32 failed_pool = 0;
    bool pair_safe = false;
};

struct Geometry {
    i64 grid = 0;
    i64 pair_ticks = 0;
    std::vector<Cell> cells;
};

Geometry build_joint_geometry(int q, int r) {
    i64 grid = 1;
    for (int speed : POOL) grid = checked_lcm(grid, 14LL * speed);
    require(grid == INT64_C(18241159416480),
            "fixed-pool denominator changed");
    grid = checked_lcm(grid, 14LL * q);
    grid = checked_lcm(grid, 14LL * r);

    std::vector<i64> walls = {0, grid};
    auto add_walls = [&](int speed) {
        const i64 unit = grid / (14LL * speed);
        require(unit * (14LL * speed) == grid, "nonintegral wall unit");
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : POOL) add_walls(speed);
    add_walls(q);
    add_walls(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.front() == 0 && walls.back() == grid,
            "joint boundary lost");

    Geometry answer;
    answer.grid = grid;
    answer.cells.reserve(walls.size() - 1);
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        require(left < right, "nonpositive joint cell");
        u32 failed = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if (!safe_midpoint(POOL[vertex], grid, left, right)) {
                failed |= u32{1} << vertex;
            }
        }
        const bool pair_safe = safe_midpoint(q, grid, left, right) &&
                               safe_midpoint(r, grid, left, right);
        if (pair_safe) answer.pair_ticks += right - left;
        answer.cells.push_back({right - left, failed, pair_safe});
    }
    return answer;
}

struct Fnv {
    u64 state = UINT64_C(0xcbf29ce484222325);
    void add(u64 word) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (word >> (8 * byte)) & UINT64_C(0xff);
            state *= UINT64_C(0x100000001b3);
        }
    }
};

u32 mask_from_labels(const std::string& token) {
    require(token.size() >= 2 && token.front() == '{' && token.back() == '}',
            "malformed label set");
    if (token == "{}") return 0;
    std::stringstream stream(token.substr(1, token.size() - 2));
    std::string item;
    u32 mask = 0;
    while (std::getline(stream, item, ',')) {
        const int label = std::stoi(item);
        const auto found = std::find(POOL.begin(), POOL.end(), label);
        require(found != POOL.end(), "unknown pool label in body");
        const unsigned bit = static_cast<unsigned>(found - POOL.begin());
        require((mask & (u32{1} << bit)) == 0, "duplicate pool label");
        mask |= u32{1} << bit;
    }
    return mask;
}

std::string labels(u32 mask) {
    std::string answer;
    for (unsigned bit = 0; bit < POOL.size(); ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        if (!answer.empty()) answer.push_back(',');
        answer += std::to_string(POOL[bit]);
    }
    return answer;
}

i128 parse_i128(const std::string& text) {
    require(!text.empty(), "empty exact integer");
    i128 answer = 0;
    for (char c : text) {
        require('0' <= c && c <= '9', "bad integer digit");
        const int digit = c - '0';
        require(answer <=
                    (std::numeric_limits<i128>::max() - digit) / 10,
                "exact integer overflow");
        answer = 10 * answer + digit;
    }
    return answer;
}

struct Parsed {
    int q = 0;
    int r = 0;
    i64 g = 0;
    i64 primitive_u = 0;
    i64 primitive_v = 0;
    i64 primitive_grid = 0;
    i64 primitive_safe_ticks = 0;
    u64 declared_prefix_size = 0;
    u64 declared_prefix_fnv = 0;
    std::vector<u32> deck;
    u64 declared_bodies = 0;
    u64 declared_failures = 0;
    u64 declared_checks = 0;
    u64 declared_max_checks = 0;
    u32 declared_worst_body = 0;
    i128 primary_mass_num = 0;
    i128 primary_mass_den = 0;
    bool saw_verdict = false;
};

Parsed parse_probe(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open probe transcript");
    Parsed answer;
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("PAIR ", 0) == 0) {
            std::istringstream stream(line);
            std::string label;
            std::string pair;
            std::string gcd_label;
            std::string primitive_label;
            std::string primitive;
            std::string grid_label;
            std::string safe_ticks_label;
            require(static_cast<bool>(
                        stream >> label >> pair >> gcd_label >> answer.g >>
                        primitive_label >> primitive >> grid_label >>
                        answer.primitive_grid >> safe_ticks_label >>
                        answer.primitive_safe_ticks),
                    "truncated pair row");
            require(label == "PAIR" && gcd_label == "G" &&
                        primitive_label == "PRIMITIVE" && grid_label == "GRID" &&
                        safe_ticks_label == "SAFE_TICKS",
                    "bad pair row labels");
            const std::size_t comma = pair.find(',');
            require(comma != std::string::npos, "bad pair token");
            answer.q = std::stoi(pair.substr(0, comma));
            answer.r = std::stoi(pair.substr(comma + 1));
            const std::size_t primitive_comma = primitive.find(',');
            require(primitive_comma != std::string::npos,
                    "bad primitive token");
            answer.primitive_u =
                std::stoll(primitive.substr(0, primitive_comma));
            answer.primitive_v =
                std::stoll(primitive.substr(primitive_comma + 1));
        } else if (line.rfind("PREFIX_CERT ", 0) == 0) {
            std::istringstream stream(line);
            std::string prefix_label;
            std::string size_label;
            std::string fnv_label;
            require(static_cast<bool>(stream >> prefix_label >> size_label >>
                                      answer.declared_prefix_size >> fnv_label >>
                                      std::hex >> answer.declared_prefix_fnv),
                    "truncated prefix declaration");
            require(prefix_label == "PREFIX_CERT" && size_label == "SIZE" &&
                        fnv_label == "FNV",
                    "bad prefix declaration labels");
        } else if (line.rfind("PREFIX_MASKS_HEX ", 0) == 0) {
            const std::string rest =
                line.substr(std::string("PREFIX_MASKS_HEX ").size());
            std::stringstream stream(rest);
            std::string token;
            while (std::getline(stream, token, ',')) {
                std::size_t consumed = 0;
                const unsigned long value = std::stoul(token, &consumed, 16);
                require(consumed == token.size() && value < (1UL << 30),
                        "bad repair mask token");
                answer.deck.push_back(static_cast<u32>(value));
            }
        } else if (line.rfind("BODY_SCAN ", 0) == 0) {
            std::istringstream stream(line);
            std::string body_scan;
            std::string bodies_label;
            std::string failures_label;
            std::string checks_label;
            std::string max_label;
            std::string worst_label;
            std::string worst_token;
            std::string first_label;
            std::string first_token;
            require(static_cast<bool>(
                stream >> body_scan >> bodies_label >> answer.declared_bodies >>
                failures_label >> answer.declared_failures >> checks_label >>
                answer.declared_checks >> max_label >> answer.declared_max_checks >>
                worst_label >> worst_token >> first_label >> first_token),
                "truncated body row");
            require(body_scan == "BODY_SCAN" && bodies_label == "BODIES" &&
                        failures_label == "FAILURES" && checks_label == "CHECKS" &&
                        max_label == "MAX_CHECKS" && worst_label == "WORST_BODY" &&
                        first_label == "FIRST_FAILURE" && first_token == "{}",
                    "bad body row labels or claimed failure");
            answer.declared_worst_body = mask_from_labels(worst_token);
        } else if (line.rfind("DECISIVE_MASS_NUM ", 0) == 0) {
            std::istringstream stream(line);
            std::string numerator_label;
            std::string numerator;
            std::string denominator_label;
            std::string denominator;
            require(static_cast<bool>(stream >> numerator_label >> numerator >>
                                      denominator_label >> denominator),
                    "truncated decisive mass row");
            require(numerator_label == "DECISIVE_MASS_NUM" &&
                        denominator_label == "DEN",
                    "bad decisive mass labels");
            answer.primary_mass_num = parse_i128(numerator);
            answer.primary_mass_den = parse_i128(denominator);
        } else if (line.rfind("VERDICT ", 0) == 0) {
            require(line.find("VERDICT EVERY_BODY_CLOSED") == 0,
                    "primary transcript is not a closure");
            answer.saw_verdict = true;
        }
    }
    require(answer.q > 0 && answer.q < answer.r, "missing or unordered pair");
    require(answer.g == std::gcd(answer.q, answer.r) &&
                answer.primitive_u == answer.q / answer.g &&
                answer.primitive_v == answer.r / answer.g &&
                answer.primitive_grid ==
                    14 * answer.primitive_u * answer.primitive_v &&
                answer.primitive_safe_ticks > 0 &&
                answer.primitive_safe_ticks <= answer.primitive_grid,
            "declared primitive scaling is inconsistent");
    require(answer.declared_prefix_size == answer.deck.size() &&
                !answer.deck.empty(),
            "prefix list/declaration mismatch");
    require(answer.declared_bodies == EXPECTED_BODIES &&
                answer.declared_failures == 0 &&
                answer.declared_max_checks == answer.deck.size(),
            "primary body declaration changed");
    require(answer.primary_mass_num > 0 && answer.primary_mass_den > 0 &&
                answer.saw_verdict,
            "incomplete primary transcript");
    require(answer.primary_mass_den ==
                static_cast<i128>(answer.primitive_grid) * answer.g *
                    INT64_C(18241159416480),
            "primary cancelled denominator is not g*D*N");
    return answer;
}

std::map<std::pair<int, int>, Parsed>
parse_probe_directory(const std::filesystem::path& directory) {
    require(std::filesystem::is_directory(directory),
            "probe path is not a directory");
    std::vector<std::filesystem::path> paths;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".out") {
            paths.push_back(entry.path());
        }
    }
    std::sort(paths.begin(), paths.end());
    require(paths.size() == EXPECTED_PAIRS.size(),
            "probe output count is not the frozen 59-edge band");
    std::map<std::pair<int, int>, Parsed> answer;
    for (const auto& path : paths) {
        Parsed parsed = parse_probe(path);
        const auto key = std::make_pair(parsed.q, parsed.r);
        require(answer.emplace(key, std::move(parsed)).second,
                "duplicate pair transcript");
    }
    for (const auto& pair : EXPECTED_PAIRS) {
        require(answer.count(pair) == 1, "frozen band pair missing");
    }
    return answer;
}

void enumerate_bodies(int next, int need, u32 mask, std::vector<u32>& bodies) {
    if (need == 0) {
        bodies.push_back(mask);
        return;
    }
    for (int bit = next; bit <= 30 - need; ++bit) {
        enumerate_bodies(bit + 1, need - 1, mask | (u32{1} << bit), bodies);
    }
}

struct BodyAudit {
    u64 bodies = 0;
    u64 failures = 0;
    u64 checks = 0;
    u64 max_checks = 0;
    u32 worst_body = 0;
};

BodyAudit audit_bodies(const std::vector<u32>& bodies,
                       const std::vector<u32>& deck) {
    const unsigned hardware = std::thread::hardware_concurrency();
    const unsigned thread_count =
        std::max(1u, std::min(8u, hardware == 0 ? 1u : hardware));
    std::vector<BodyAudit> locals(thread_count);
    std::vector<std::thread> threads;
    for (unsigned thread = 0; thread < thread_count; ++thread) {
        threads.emplace_back([&, thread]() {
            BodyAudit local;
            const std::size_t begin = bodies.size() * thread / thread_count;
            const std::size_t end = bodies.size() * (thread + 1) / thread_count;
            for (std::size_t index = begin; index < end; ++index) {
                const u32 body = bodies[index];
                u64 used = 0;
                bool covered = false;
                for (u32 repair : deck) {
                    ++used;
                    if ((body & repair) == 0) {
                        covered = true;
                        break;
                    }
                }
                ++local.bodies;
                local.checks += used;
                if (!covered) ++local.failures;
                if (used > local.max_checks ||
                    (used == local.max_checks && body < local.worst_body)) {
                    local.max_checks = used;
                    local.worst_body = body;
                }
            }
            locals[thread] = local;
        });
    }
    for (auto& thread : threads) thread.join();
    BodyAudit answer;
    for (const BodyAudit& local : locals) {
        answer.bodies += local.bodies;
        answer.failures += local.failures;
        answer.checks += local.checks;
        if (local.max_checks > answer.max_checks ||
            (local.max_checks == answer.max_checks &&
             local.worst_body < answer.worst_body)) {
            answer.max_checks = local.max_checks;
            answer.worst_body = local.worst_body;
        }
    }
    return answer;
}

struct PairAudit {
    u64 joint_cells = 0;
    u64 prefix_size = 0;
    u64 body_checks = 0;
};

PairAudit audit_pair(const Parsed& parsed, const std::vector<u32>& bodies,
                     Fnv& pair_labelled_prefix_incidence_ledger) {
    Fnv prefix_hash;
    std::vector<u32> unique = parsed.deck;
    for (u32 repair : parsed.deck) {
        require(std::popcount(repair) == 8, "repair has wrong weight");
        prefix_hash.add(repair);
        // Retain repeated masks across distinct pair prefixes: this is an
        // incidence ledger, not the set-theoretic union of repair masks.
        pair_labelled_prefix_incidence_ledger.add(repair);
    }
    std::sort(unique.begin(), unique.end());
    require(std::adjacent_find(unique.begin(), unique.end()) == unique.end(),
            "duplicate repair in prefix");
    require(prefix_hash.state == parsed.declared_prefix_fnv,
            "prefix FNV mismatch");

    const Geometry geometry = build_joint_geometry(parsed.q, parsed.r);
    const i128 direct_beta_divisor =
        gcd128(geometry.pair_ticks, geometry.grid);
    const i128 primitive_beta_divisor =
        gcd128(parsed.primitive_safe_ticks, parsed.primitive_grid);
    require(static_cast<i128>(geometry.pair_ticks) / direct_beta_divisor ==
                    static_cast<i128>(parsed.primitive_safe_ticks) /
                        primitive_beta_divisor &&
                static_cast<i128>(geometry.grid) / direct_beta_divisor ==
                    static_cast<i128>(parsed.primitive_grid) /
                        primitive_beta_divisor,
            "direct pair mass disagrees with primitive beta");
    i128 minimum_margin = 0;
    u32 minimum_repair = 0;
    i64 final_mass = 0;
    for (std::size_t index = 0; index < parsed.deck.size(); ++index) {
        const u32 repair = parsed.deck[index];
        i64 mass = 0;
        for (const Cell& cell : geometry.cells) {
            if (cell.pair_safe && (cell.failed_pool & ~repair) == 0) {
                mass += cell.width;
            }
        }
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * geometry.grid;
        require(margin >= 0, "inactive repair in emitted prefix");
        if (index == 0 || margin < minimum_margin) {
            minimum_margin = margin;
            minimum_repair = repair;
        }
        if (index + 1 == parsed.deck.size()) final_mass = mass;
    }
    const i128 direct_divisor = gcd128(final_mass, geometry.grid);
    const i128 primary_divisor =
        gcd128(parsed.primary_mass_num, parsed.primary_mass_den);
    require(static_cast<i128>(final_mass) / direct_divisor ==
                    parsed.primary_mass_num / primary_divisor &&
                static_cast<i128>(geometry.grid) / direct_divisor ==
                    parsed.primary_mass_den / primary_divisor,
            "direct final mass disagrees with primary cocycle mass");

    const BodyAudit body = audit_bodies(bodies, parsed.deck);
    require(body.bodies == EXPECTED_BODIES && body.failures == 0,
            "recursive body quantifier failed");
    require(body.checks == parsed.declared_checks &&
                body.max_checks == parsed.declared_max_checks &&
                body.worst_body == parsed.declared_worst_body,
            "independent body ledger disagrees with primary");
    require(body.max_checks == parsed.deck.size(),
            "prefix is not order-minimal");
    require((body.worst_body & parsed.deck.back()) == 0,
            "final repair does not miss minimality witness");
    for (std::size_t index = 0; index + 1 < parsed.deck.size(); ++index) {
        require((body.worst_body & parsed.deck[index]) != 0,
                "earlier repair misses minimality witness");
    }

    std::cout << "PAIR " << parsed.q << ',' << parsed.r
              << " G " << parsed.g
              << " PRIMITIVE " << parsed.primitive_u << ','
              << parsed.primitive_v
              << " N " << parsed.primitive_grid
              << " S " << parsed.primitive_safe_ticks
              << " PRIMARY_DEN " << decimal(parsed.primary_mass_den)
              << " GRID " << geometry.grid
              << " JOINT_CELLS " << geometry.cells.size()
              << " PAIR_TICKS " << geometry.pair_ticks << '/' << geometry.grid
              << " PREFIX " << parsed.deck.size()
              << " PREFIX_FNV " << std::hex << prefix_hash.state << std::dec
              << " MIN_MARGIN_NUM " << decimal(minimum_margin)
              << " MIN_REPAIR {" << labels(minimum_repair) << "}\n";
    std::cout << "BODY_RECURSION COUNT " << body.bodies
              << " FAILURES " << body.failures
              << " CHECKS " << body.checks
              << " MAX_CHECKS " << body.max_checks
              << " MINIMALITY_WITNESS {" << labels(body.worst_body) << "}\n";
    std::cout << "FINAL_REPAIR_DIRECT_MASS " << final_mass << '/'
              << geometry.grid << " SCALED_MARGIN_63MU_MINUS4 "
              << fraction(static_cast<i128>(63) * final_mass -
                              static_cast<i128>(4) * geometry.grid,
                          geometry.grid)
              << " CROSS_PRIMARY PASS\n";

    return {static_cast<u64>(geometry.cells.size()),
            static_cast<u64>(parsed.deck.size()), body.checks};
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: cleanroom-audit replay_band_directory");
        const auto parsed = parse_probe_directory(argv[1]);

        Fnv band_ledger;
        for (const auto& [q, r] : EXPECTED_PAIRS) {
            band_ledger.add(static_cast<u64>(q));
            band_ledger.add(static_cast<u64>(r));
        }
        require(band_ledger.state == EXPECTED_BAND_FNV,
                "hard-coded band fingerprint changed");

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "recursive body universe changed");

        std::cout << "ENDPOINT_CASCADE_CLEANROOM_DIRECT_WALL_BODY_V1\n";
        std::cout << "BAND PAIRS " << EXPECTED_PAIRS.size()
                  << " ENDPOINTS 755..768 FNV " << std::hex
                  << band_ledger.state << std::dec
                  << " BODIES_PER_PAIR " << EXPECTED_BODIES << '\n';

        Fnv pair_labelled_prefix_incidence_ledger;
        u64 total_cells = 0;
        u64 total_prefix_repairs = 0;
        u64 total_body_checks = 0;
        for (const auto& pair : EXPECTED_PAIRS) {
            const auto found = parsed.find(pair);
            require(found != parsed.end(), "expected pair vanished after parse");
            pair_labelled_prefix_incidence_ledger.add(
                static_cast<u64>(pair.first));
            pair_labelled_prefix_incidence_ledger.add(
                static_cast<u64>(pair.second));
            pair_labelled_prefix_incidence_ledger.add(
                found->second.deck.size());
            const PairAudit audit =
                audit_pair(found->second, bodies,
                           pair_labelled_prefix_incidence_ledger);
            total_cells += audit.joint_cells;
            total_prefix_repairs += audit.prefix_size;
            total_body_checks += audit.body_checks;
        }

        std::cout << "GLOBAL_LEDGER JOINT_CELLS " << total_cells
                  << " PREFIX_OCCURRENCES " << total_prefix_repairs
                  << " BODY_CASES "
                  << EXPECTED_BODIES * EXPECTED_PAIRS.size()
                  << " BODY_CHECKS " << total_body_checks
                  << " PAIR_LABELLED_PREFIX_INCIDENCE_FNV " << std::hex
                  << pair_labelled_prefix_incidence_ledger.state << std::dec
                  << " ALL_PREFIX_MARGINS_STRICTLY_POSITIVE 1\n";
        std::cout << "CHECKS PASS fresh_joint_walls,direct_midpoint_mass,"
                     "recursive_bodies,prefix_order_minimality,cross_primary\n";
        std::cout << "GLOBAL_VERDICT PASS FRESH_JOINT_WALLS 59_PREFIXES "
                     "FULL_30_POOL_9_BODY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT_CASCADE_CLEANROOM_ERROR " << error.what() << '\n';
        return 1;
    }
}
