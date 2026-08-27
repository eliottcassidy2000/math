// Primary exact certificate for THM-4252.
//
// This path keeps the fixed-pool wall cells, evaluates the primitive centered
// observable at every signed cell endpoint, groups exact cell contributions by
// their labelled failure mask, applies an ordinary-colex superset transform to
// all eight-deletion repairs, and exhausts every labelled nine-body.
//
// All consequence-bearing gates use require(), which remains active under
// NDEBUG.  No floating point, sampling, or C/C++ assert participates.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 COMMON = INT64_C(18241159416480);

constexpr u64 EXPECTED_REPAIRS = UINT64_C(5852925);
constexpr u64 EXPECTED_BODIES = UINT64_C(14307150);
constexpr i128 THRESHOLD_DEN = static_cast<i128>(8281) * 3467;
constexpr i128 THRESHOLD_NUM = 1650;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

i64 exact_lcm(i64 left, i64 right) {
    const i64 divisor = std::gcd(left, right);
    const i128 value = static_cast<i128>(left / divisor) * right;
    require(value <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(value);
}

u32 next_combination(u32 mask) {
    const u32 low = mask & (~mask + 1u);
    const u32 ripple = mask + low;
    if (ripple == 0) return 0;
    return ripple | (((mask ^ ripple) >> 2) / low);
}

bool pool_midpoint_safe(int speed, i64 left, i64 right) {
    const i64 period = 2 * COMMON;
    i128 residue = static_cast<i128>(speed) * (left + right);
    residue %= period;
    if (residue < 0) residue += period;
    return static_cast<i128>(7) * residue >= COMMON &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * COMMON;
}

struct Cell {
    i64 left = 0;
    i64 right = 0;
    u32 failed_pool = 0;
};

std::vector<Cell> build_pool_cells() {
    i64 check = 1;
    for (int speed : POOL) check = exact_lcm(check, 14LL * speed);
    require(check == COMMON, "fixed-pool common denominator changed");
    std::vector<i64> walls = {0, COMMON};
    for (int speed : POOL) {
        const i64 unit = COMMON / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.size() == 7134, "fixed-pool wall count changed");
    std::vector<Cell> cells;
    cells.reserve(walls.size() - 1);
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        u32 failed = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if (!pool_midpoint_safe(POOL[vertex], left, right)) {
                failed |= u32{1} << vertex;
            }
        }
        cells.push_back({left, right, failed});
    }
    return cells;
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

std::string labels(u32 mask) {
    std::string answer;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((mask & (u32{1} << vertex)) == 0) continue;
        if (!answer.empty()) answer += ',';
        answer += std::to_string(POOL[vertex]);
    }
    return answer;
}

std::array<std::array<u64, 9>, 31> choose8_local{};

void init_choose8_local() {
    for (int n = 0; n <= 30; ++n) {
        choose8_local[n][0] = 1;
        for (int k = 1; k <= 8; ++k) {
            choose8_local[n][k] = n == 0 ? 0 :
                choose8_local[n - 1][k] + choose8_local[n - 1][k - 1];
        }
    }
    require(choose8_local[30][8] == EXPECTED_REPAIRS,
            "rank-eight universe changed");
}

u64 colex_rank8_local(u32 mask) {
    u64 rank = 0;
    int ordinal = 1;
    for (int bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        rank += choose8_local[bit][ordinal];
        ++ordinal;
    }
    require(ordinal == 9 && rank < EXPECTED_REPAIRS, "bad rank-eight mask");
    return rank;
}

struct Reduced {
    i128 num = 0;
    i128 den = 1;
};

i128 gcd_i128(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        const i128 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

Reduced reduced(i128 num, i128 den) {
    require(den != 0, "zero denominator");
    if (den < 0) { num = -num; den = -den; }
    const i128 d = gcd_i128(num, den);
    return {num / d, den / d};
}

std::string show(Reduced x) {
    return decimal(x.num) + "/" + decimal(x.den);
}

struct PrimitiveCell {
    i64 left = 0;
    i64 right = 0;
    bool safe = false;
};

struct PrimitivePair {
    i64 u = 0;
    i64 v = 0;
    i64 grid = 0;
    i64 safe_ticks = 0;
    std::vector<i64> points;
    std::vector<PrimitiveCell> cells;
    std::vector<i64> safe_prefix_ticks;
    i128 min_raw = 0;
    i128 max_raw = 0;
    i64 min_at = 0;
    i64 max_at = 0;
};

i64 mod_tick(i128 x, i64 n) {
    x %= n;
    if (x < 0) x += n;
    return static_cast<i64>(x);
}

bool primitive_midpoint_safe(i64 speed, i64 grid, i64 left, i64 right) {
    i128 residue = static_cast<i128>(speed) * (left + right);
    residue %= static_cast<i128>(2) * grid;
    if (residue < 0) residue += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * residue >= grid &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * grid;
}

PrimitivePair build_primitive(i64 input_u, i64 input_v) {
    i64 u = input_u;
    i64 v = input_v;
    if (u > v) std::swap(u, v);
    require(u > 0 && u < v && std::gcd(u, v) == 1,
            "primitive pair invalid");
    const i128 n128 = static_cast<i128>(14) * u * v;
    require(n128 <= std::numeric_limits<i64>::max(), "primitive grid overflow");
    const i64 n = static_cast<i64>(n128);
    PrimitivePair p;
    p.u = u;
    p.v = v;
    p.grid = n;
    p.points.reserve(static_cast<std::size_t>(2 * (u + v) + 2));
    p.points.push_back(0);
    p.points.push_back(n);
    for (i64 i = 0; i < u; ++i) {
        p.points.push_back(mod_tick(static_cast<i128>(v) * (14 * i - 1), n));
        p.points.push_back(mod_tick(static_cast<i128>(v) * (14 * i + 1), n));
    }
    for (i64 i = 0; i < v; ++i) {
        p.points.push_back(mod_tick(static_cast<i128>(u) * (14 * i - 1), n));
        p.points.push_back(mod_tick(static_cast<i128>(u) * (14 * i + 1), n));
    }
    std::sort(p.points.begin(), p.points.end());
    p.points.erase(std::unique(p.points.begin(), p.points.end()), p.points.end());
    require(p.points.front() == 0 && p.points.back() == n,
            "primitive boundary lost");
    p.cells.reserve(p.points.size() - 1);
    p.safe_prefix_ticks.assign(p.points.size(), 0);
    for (std::size_t i = 1; i < p.points.size(); ++i) {
        const i64 left = p.points[i - 1];
        const i64 right = p.points[i];
        const bool safe = primitive_midpoint_safe(u, n, left, right) &&
                          primitive_midpoint_safe(v, n, left, right);
        p.cells.push_back({left, right, safe});
        p.safe_prefix_ticks[i] = p.safe_prefix_ticks[i - 1] +
                                 (safe ? right - left : 0);
    }
    p.safe_ticks = p.safe_prefix_ticks.back();
    for (std::size_t i = 0; i < p.points.size(); ++i) {
        const i128 raw = static_cast<i128>(p.safe_prefix_ticks[i]) * n -
                         static_cast<i128>(p.safe_ticks) * p.points[i];
        if (raw < p.min_raw) { p.min_raw = raw; p.min_at = p.points[i]; }
        if (raw > p.max_raw) { p.max_raw = raw; p.max_at = p.points[i]; }
    }
    return p;
}

// F_A(p/COMMON) numerator on denominator grid*COMMON.
i128 primitive_safe_prefix_num(const PrimitivePair& p, i64 phase_tick) {
    require(0 <= phase_tick && phase_tick < COMMON, "phase out of range");
    const i128 y = static_cast<i128>(p.grid) * phase_tick;
    std::size_t lo = 0;
    std::size_t hi = p.points.size();
    while (lo + 1 < hi) {
        const std::size_t mid = lo + (hi - lo) / 2;
        if (static_cast<i128>(p.points[mid]) * COMMON <= y) lo = mid;
        else hi = mid;
    }
    if (lo + 1 == p.points.size()) {
        return static_cast<i128>(p.safe_ticks) * COMMON;
    }
    i128 answer = static_cast<i128>(p.safe_prefix_ticks[lo]) * COMMON;
    if (p.cells[lo].safe) {
        answer += y - static_cast<i128>(p.points[lo]) * COMMON;
    }
    return answer;
}

// H_A(g*t/COMMON) numerator on denominator grid*COMMON.
i128 endpoint_h_num(const PrimitivePair& p, i64 g, i64 tick) {
    const i128 product = static_cast<i128>(g) * tick;
    const i64 phase = static_cast<i64>(product % COMMON);
    return primitive_safe_prefix_num(p, phase) -
           static_cast<i128>(p.safe_ticks) * phase;
}

// Integral_0^(tick/COMMON) 1_A(g y)dy numerator on grid*g*COMMON.
i128 actual_safe_prefix_num(const PrimitivePair& p, i64 g, i64 tick) {
    const i128 product = static_cast<i128>(g) * tick;
    const i128 whole = product / COMMON;
    const i64 phase = static_cast<i64>(product % COMMON);
    return whole * p.safe_ticks * COMMON + primitive_safe_prefix_num(p, phase);
}

struct AtomData {
    std::map<u32, i128> mass;
    u64 cells_checked = 0;
    i128 total_cocycle = 0;
};

AtomData build_cocycle_atoms(const std::vector<Cell>& pool_cells,
                             const PrimitivePair& p, i64 g) {
    AtomData out;
    i128 prior_direct = actual_safe_prefix_num(p, g, 0);
    for (const Cell& cell : pool_cells) {
        const i128 h_left = endpoint_h_num(p, g, cell.left);
        const i128 h_right = endpoint_h_num(p, g, cell.right);
        const i128 cocycle = h_right - h_left;
        const i128 contribution = static_cast<i128>(g) * p.safe_ticks *
                                      (cell.right - cell.left) + cocycle;
        const i128 current_direct = actual_safe_prefix_num(p, g, cell.right);
        require(contribution == current_direct - prior_direct,
                "cocycle/direct cell disagreement");
        require(contribution >= 0, "negative literal cell mass");
        if (std::popcount(cell.failed_pool) <= 8) {
            out.mass[cell.failed_pool] += contribution;
        }
        out.total_cocycle += cocycle;
        prior_direct = current_direct;
        ++out.cells_checked;
    }
    require(out.total_cocycle == 0, "full-circle cocycle did not telescope");
    require(prior_direct == static_cast<i128>(g) * p.safe_ticks * COMMON,
            "pair total mass changed");
    return out;
}

void add_supersets_pair(u32 atom, int need, int start, u32 extra, i128 value,
                        std::vector<i128>& masses, u64& operations) {
    if (need == 0) {
        masses[colex_rank8_local(atom | extra)] += value;
        ++operations;
        return;
    }
    for (int bit = start; bit <= 30 - need; ++bit) {
        const u32 flag = u32{1} << bit;
        if ((atom & flag) != 0) continue;
        add_supersets_pair(atom, need - 1, bit + 1, extra | flag,
                           value, masses, operations);
    }
}

u64 mix64(u64 x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

struct FnvLocal {
    u64 state = UINT64_C(0xcbf29ce484222325);
    void add(u64 x) {
        for (unsigned j = 0; j < 8; ++j) {
            state ^= (x >> (8 * j)) & 0xffu;
            state *= UINT64_C(0x100000001b3);
        }
    }
};

struct ScanResult {
    u64 bodies = 0;
    u64 failures = 0;
    u64 checks = 0;
    u64 max_checks = 0;
    u32 worst_body = 0;
    u32 worst_repair = 0;
    u32 first_failure = std::numeric_limits<u32>::max();
};

ScanResult scan_bodies(const std::vector<u32>& active) {
    std::vector<u32> bodies;
    bodies.reserve(EXPECTED_BODIES);
    u32 body = (u32{1} << 9) - 1;
    const u32 limit = u32{1} << 30;
    while (body != 0 && body < limit) {
        bodies.push_back(body);
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(bodies.size() == EXPECTED_BODIES, "body universe changed");

    const unsigned hardware = std::thread::hardware_concurrency();
    const unsigned thread_count = std::max(1u, std::min(8u, hardware ? hardware : 1u));
    std::vector<ScanResult> locals(thread_count);
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (unsigned thread = 0; thread < thread_count; ++thread) {
        workers.emplace_back([&, thread]() {
            ScanResult local;
            const std::size_t begin = bodies.size() * thread / thread_count;
            const std::size_t end = bodies.size() * (thread + 1) / thread_count;
            for (std::size_t index = begin; index < end; ++index) {
                const u32 b = bodies[index];
                u64 used = 0;
                u32 first = 0;
                for (u32 repair : active) {
                    ++used;
                    if ((b & repair) == 0) { first = repair; break; }
                }
                ++local.bodies;
                local.checks += used;
                if (first == 0) {
                    ++local.failures;
                    local.first_failure = std::min(local.first_failure, b);
                }
                if (used > local.max_checks ||
                    (used == local.max_checks && b < local.worst_body)) {
                    local.max_checks = used;
                    local.worst_body = b;
                    local.worst_repair = first;
                }
            }
            locals[thread] = local;
        });
    }
    for (auto& worker : workers) worker.join();
    ScanResult all;
    all.first_failure = std::numeric_limits<u32>::max();
    for (const ScanResult& local : locals) {
        all.bodies += local.bodies;
        all.failures += local.failures;
        all.checks += local.checks;
        all.first_failure = std::min(all.first_failure, local.first_failure);
        if (local.max_checks > all.max_checks ||
            (local.max_checks == all.max_checks && local.worst_body < all.worst_body)) {
            all.max_checks = local.max_checks;
            all.worst_body = local.worst_body;
            all.worst_repair = local.worst_repair;
        }
    }
    return all;
}

struct ComponentReport {
    i64 components = 0;
    i64 positive = 0;
    i64 negative = 0;
    i64 zero = 0;
    i128 base_mass_ticks = 0;
    i128 total_hdiff = 0;
    i128 exact_mass_num = 0;
    i128 max_abs_hdiff = 0;
    i64 max_abs_left = 0;
    i64 max_abs_right = 0;
    FnvLocal endpoint_ledger;
};

ComponentReport component_report(const std::vector<Cell>& cells,
                                 const PrimitivePair& p, i64 g, u32 repair) {
    std::vector<unsigned char> safe(cells.size(), 0);
    for (std::size_t i = 0; i < cells.size(); ++i) {
        safe[i] = static_cast<unsigned char>((cells[i].failed_pool & ~repair) == 0);
    }
    require(!safe.front() && !safe.back(),
            "repair unexpectedly has a component crossing the chosen origin");
    ComponentReport out;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const std::size_t prev = (i + cells.size() - 1) % cells.size();
        if (!safe[i] || safe[prev]) continue;
        ++out.components;
        const i64 left = cells[i].left;
        std::size_t j = i;
        i128 width = 0;
        while (j < cells.size() && safe[j]) {
            width += cells[j].right - cells[j].left;
            ++j;
        }
        require(j > i, "empty component");
        const i64 right = cells[j - 1].right;
        const i128 hdiff = endpoint_h_num(p, g, right) -
                           endpoint_h_num(p, g, left);
        if (hdiff > 0) ++out.positive;
        else if (hdiff < 0) ++out.negative;
        else ++out.zero;
        out.base_mass_ticks += width;
        out.total_hdiff += hdiff;
        out.exact_mass_num += static_cast<i128>(g) * p.safe_ticks * width + hdiff;
        const i128 abs_hdiff = hdiff < 0 ? -hdiff : hdiff;
        if (abs_hdiff > out.max_abs_hdiff) {
            out.max_abs_hdiff = abs_hdiff;
            out.max_abs_left = left;
            out.max_abs_right = right;
        }
        out.endpoint_ledger.add(static_cast<u64>(left));
        out.endpoint_ledger.add(static_cast<u64>(right));
        out.endpoint_ledger.add(static_cast<u64>(hdiff));
    }
    require(out.components > 0, "repair has no components");
    return out;
}

u32 mask_for_rank8(u64 wanted) {
    u32 mask = (u32{1} << 8) - 1;
    const u32 limit = u32{1} << 30;
    u64 rank = 0;
    while (mask != 0 && mask < limit) {
        if (rank == wanted) return mask;
        ++rank;
        const u32 next = next_combination(mask);
        if (next <= mask) break;
        mask = next;
    }
    require(false, "rank not found");
    return 0;
}

struct Expected {
    i64 q;
    i64 r;
    i64 g;
    i64 u;
    i64 v;
    i64 grid;
    i64 safe_ticks;
    i128 min_raw;
    i64 min_at;
    i128 max_raw;
    i64 max_at;
    u64 active;
    u64 active_ledger;
    u64 body_checks;
    u64 prefix_size;
    u64 prefix_ledger;
    u32 worst_body;
    u32 decisive_repair;
    i128 decisive_mass;
    i64 components;
    i64 positive;
    i64 negative;
    i128 total_hdiff;
    u64 endpoint_ledger;
};

constexpr std::array<Expected, 3> EXPECTED = {{
    {466, 699, 233, 2, 3, 84, 64,
     -268, 58, 268, 26,
     UINT64_C(5114702), UINT64_C(0x33a4fecc55e46462),
     UINT64_C(413551394), UINT64_C(481),
     UINT64_C(0xa5461a33cd4e1e8c), UINT32_C(0x0f722000),
     UINT32_C(0x00890586), static_cast<i128>(UINT64_C(22736853453460272)),
     196, 94, 102, -static_cast<i128>(UINT64_C(1102595351093584)),
     UINT64_C(0x848edf9116a8f904)},
    {616, 769, 1, 616, 769, INT64_C(6631856), INT64_C(4872384),
     -static_cast<i128>(UINT64_C(22034815264)), INT64_C(4588584),
     static_cast<i128>(UINT64_C(22034815264)), INT64_C(2043272),
     UINT64_C(4632061), UINT64_C(0xe7fb7af64ae2edcc),
     UINT64_C(418577253), UINT64_C(709),
     UINT64_C(0xf07e55a1f5c2b1b6), UINT32_C(0x0ce07400),
     UINT32_C(0x200e0289), static_cast<i128>(UINT64_C(8868855487806352192)),
     212, 100, 112, -static_cast<i128>(UINT64_C(545177584978334528)),
     UINT64_C(0x1a1e4c9470b48089)},
    {721, 769, 1, 721, 769, INT64_C(7762286), INT64_C(5702904),
     -static_cast<i128>(UINT64_C(19807136076)), INT64_C(6612291),
     static_cast<i128>(UINT64_C(19807136076)), INT64_C(1149995),
     UINT64_C(4762470), UINT64_C(0x059c1899134692d7),
     UINT64_C(414845236), UINT64_C(672),
     UINT64_C(0x88a16a55ce12cb5a), UINT32_C(0x06cc3001),
     UINT32_C(0x01320318), static_cast<i128>(UINT64_C(9128326347957259072)),
     186, 98, 88, -static_cast<i128>(UINT64_C(316010002749365312)),
     UINT64_C(0x36d48e458b4a3f2c)}
}};

void run_pair(const Expected& expected, const std::vector<Cell>& pool_cells) {
        const i64 q = expected.q;
        const i64 r = expected.r;
        const i64 g = std::gcd(q, r);
        const PrimitivePair primitive = build_primitive(q / g, r / g);
        require(g == expected.g && primitive.u == expected.u &&
                    primitive.v == expected.v && primitive.grid == expected.grid &&
                    primitive.safe_ticks == expected.safe_ticks &&
                    primitive.min_raw == expected.min_raw &&
                    primitive.min_at == expected.min_at &&
                    primitive.max_raw == expected.max_raw &&
                    primitive.max_at == expected.max_at,
                "primitive expected control changed");
        const AtomData atoms = build_cocycle_atoms(pool_cells, primitive, g);
        std::vector<i128> masses(EXPECTED_REPAIRS, 0);
        u64 zeta_operations = 0;
        for (const auto& [mask, value] : atoms.mass) {
            add_supersets_pair(mask, 8 - std::popcount(mask), 0, 0,
                               value, masses, zeta_operations);
        }
        require(zeta_operations == UINT64_C(152170690),
                "zeta operation count changed");
        const i128 full_den = static_cast<i128>(primitive.grid) * g * COMMON;
        std::vector<u32> active;
        active.reserve(EXPECTED_REPAIRS);
        i128 minimum_active_margin = 0;
        i128 maximum_inactive_margin = 0;
        bool saw_inactive = false;
        u64 minimum_active_rank = 0;
        u64 maximum_inactive_rank = 0;
        u64 equalities = 0;
        u32 repair = (u32{1} << 8) - 1;
        const u32 limit = u32{1} << 30;
        u64 rank = 0;
        while (repair != 0 && repair < limit) {
            const i128 margin = static_cast<i128>(63) * masses[rank] -
                                static_cast<i128>(4) * full_den;
            if (margin >= 0) {
                if (active.empty() || margin < minimum_active_margin) {
                    minimum_active_margin = margin;
                    minimum_active_rank = rank;
                }
                if (margin == 0) ++equalities;
                active.push_back(repair);
            } else if (!saw_inactive || margin > maximum_inactive_margin) {
                saw_inactive = true;
                maximum_inactive_margin = margin;
                maximum_inactive_rank = rank;
            }
            ++rank;
            const u32 next = next_combination(repair);
            if (next <= repair) break;
            repair = next;
        }
        require(rank == EXPECTED_REPAIRS, "repair enumeration incomplete");
        require(saw_inactive, "inactive repair control disappeared");

        constexpr u64 ORDER_SEED = UINT64_C(0x4245422842334245);
        std::sort(active.begin(), active.end(), [](u32 a, u32 b) {
            constexpr u64 seed = UINT64_C(0x4245422842334245);
            const u64 ka = mix64(static_cast<u64>(a) ^ seed);
            const u64 kb = mix64(static_cast<u64>(b) ^ seed);
            return ka != kb ? ka < kb : a < b;
        });
        FnvLocal active_ledger;
        for (u32 mask : active) active_ledger.add(mask);
        const ScanResult scan = scan_bodies(active);

        FnvLocal prefix_ledger;
        for (u64 i = 0; i < scan.max_checks && i < active.size(); ++i) {
            prefix_ledger.add(active[i]);
        }

        require(active.size() == expected.active &&
                    active_ledger.state == expected.active_ledger &&
                    equalities == 0,
                "active repair hypergraph changed");
        require(scan.bodies == EXPECTED_BODIES && scan.failures == 0 &&
                    scan.checks == expected.body_checks &&
                    scan.max_checks == expected.prefix_size &&
                    scan.worst_body == expected.worst_body &&
                    scan.worst_repair == expected.decisive_repair,
                "full nine-body scan changed");
        require(prefix_ledger.state == expected.prefix_ledger,
                "prefix certificate ledger changed");

        const u32 closest_active = mask_for_rank8(minimum_active_rank);
        u32 decisive_repair = scan.worst_repair;
        const i128 decisive_mass = masses[colex_rank8_local(decisive_repair)];
        const ComponentReport comp = component_report(pool_cells, primitive, g,
                                                      decisive_repair);
        require(comp.exact_mass_num == decisive_mass,
                "component sum disagrees with zeta mass");
        require(decisive_repair == expected.decisive_repair &&
                    decisive_mass == expected.decisive_mass &&
                    comp.components == expected.components &&
                    comp.positive == expected.positive &&
                    comp.negative == expected.negative && comp.zero == 0 &&
                    comp.total_hdiff == expected.total_hdiff &&
                    comp.endpoint_ledger.state == expected.endpoint_ledger,
                "decisive component certificate changed");

        const Reduced beta = reduced(primitive.safe_ticks, primitive.grid);
        const Reduced omega = reduced(primitive.max_raw - primitive.min_raw,
                                      static_cast<i128>(primitive.grid) * primitive.grid);
        const Reduced omega_over_g = reduced(omega.num, omega.den * g);
        const bool beta_gate = static_cast<i128>(91) * beta.num >=
                               static_cast<i128>(66) * beta.den;
        const bool omega_gate = omega_over_g.num * THRESHOLD_DEN <=
                                THRESHOLD_NUM * omega_over_g.den;
        const Reduced exact_margin = reduced(static_cast<i128>(63) * decisive_mass -
                                                 static_cast<i128>(4) * full_den,
                                             full_den);
        const Reduced total_error = reduced(comp.total_hdiff, full_den);
        const Reduced scalar_budget = reduced(static_cast<i128>(comp.components) *
                                                  (primitive.max_raw - primitive.min_raw),
                                              static_cast<i128>(primitive.grid) *
                                                  primitive.grid * g);
        const Reduced base_term = reduced(static_cast<i128>(primitive.safe_ticks) *
                                              comp.base_mass_ticks,
                                          static_cast<i128>(primitive.grid) * COMMON);
        require(beta_gate && !omega_gate,
                "named edge no longer has the declared scalar-hostile status");

        std::cout << "BEGIN_PAIR " << q << ',' << r << '\n';
        std::cout << "PAIR " << q << ',' << r << " G " << g
                  << " PRIMITIVE " << primitive.u << ',' << primitive.v
                  << " GRID " << primitive.grid << " SAFE_TICKS "
                  << primitive.safe_ticks << " BETA " << show(beta)
                  << " OMEGA " << show(omega) << " OMEGA_OVER_G "
                  << show(omega_over_g) << " BETA_GATE " << beta_gate
                  << " SCALAR_GATE " << (beta_gate && omega_gate) << '\n';
        std::cout << "PRIMITIVE_EXTREMA MIN_RAW " << decimal(primitive.min_raw)
                  << " AT " << primitive.min_at << " MAX_RAW "
                  << decimal(primitive.max_raw) << " AT " << primitive.max_at
                  << " THRESHOLD 1650/28710227\n";
        std::cout << "COCYCLE_CELLS " << atoms.cells_checked
                  << " ATOMS " << atoms.mass.size() << " TELESCOPING "
                  << decimal(atoms.total_cocycle) << " DIRECT_CELL_AUDIT PASS\n";
        std::cout << "REPAIRS " << EXPECTED_REPAIRS << " ACTIVE " << active.size()
                  << " INACTIVE " << EXPECTED_REPAIRS - active.size()
                  << " EQUALITIES " << equalities << " ZETA_OPS "
                  << zeta_operations << " ACTIVE_LEDGER " << std::hex
                  << active_ledger.state << std::dec << '\n';
        std::cout << "BOUNDARY_REPAIRS MIN_ACTIVE {" << labels(closest_active)
                  << "} MARGIN_NUM " << decimal(minimum_active_margin)
                  << " MAX_INACTIVE {" << labels(mask_for_rank8(maximum_inactive_rank))
                  << "} MARGIN_NUM " << decimal(maximum_inactive_margin) << '\n';
        std::cout << "BODY_SCAN BODIES " << scan.bodies << " FAILURES "
                  << scan.failures << " CHECKS " << scan.checks << " MAX_CHECKS "
                  << scan.max_checks << " WORST_BODY {" << labels(scan.worst_body)
                  << "} FIRST_FAILURE {" <<
                     (scan.failures ? labels(scan.first_failure) : std::string())
                  << "}\n";
        std::cout << "PREFIX_CERT SIZE " << scan.max_checks << " FNV " << std::hex
                  << prefix_ledger.state << std::dec << " ORDER_SEED " << std::hex
                  << ORDER_SEED << std::dec << " DECISIVE_REPAIR {"
                  << labels(decisive_repair) << "}\n";
        std::cout << "PREFIX_MASKS_HEX";
        for (u64 i = 0; i < scan.max_checks && i < active.size(); ++i) {
            std::cout << (i == 0 ? " " : ",") << std::hex << active[i];
        }
        std::cout << std::dec << '\n';
        std::cout << "COMPONENT_REPORT COMPONENTS " << comp.components
                  << " POS " << comp.positive << " NEG " << comp.negative
                  << " ZERO " << comp.zero << " BASE_MASS_TICKS "
                  << decimal(comp.base_mass_ticks) << " BASE_TERM " << show(base_term)
                  << " TOTAL_HDIFF " << decimal(comp.total_hdiff)
                  << " ERROR " << show(total_error) << " SCALAR_BUDGET "
                  << show(scalar_budget) << " MAX_ABS_HDIFF "
                  << decimal(comp.max_abs_hdiff) << " MAX_ADDR "
                  << comp.max_abs_left << ',' << comp.max_abs_right
                  << " ENDPOINT_LEDGER " << std::hex << comp.endpoint_ledger.state
                  << std::dec << '\n';
        std::cout << "DECISIVE_MASS_NUM " << decimal(decisive_mass)
                  << " DEN " << decimal(full_den) << " SCALED_MARGIN_63MU_MINUS4 "
                  << show(exact_margin) << " COMPONENT_IDENTITY PASS\n";
        std::cout << "VERDICT EVERY_BODY_CLOSED\n";
        std::cout << "END_PAIR " << q << ',' << r << '\n';
}

} // namespace

int main() {
    try {
        init_choose8_local();
        const std::vector<Cell> pool_cells = build_pool_cells();
        require(pool_cells.size() == 7133, "pool cell count changed");
        std::cout << "THM4252_PRIMARY_EXACT_ENDPOINT_COCYCLE_V1\n";
        std::cout << "POOL_CELLS " << pool_cells.size()
                  << " REPAIRS " << EXPECTED_REPAIRS
                  << " BODIES " << EXPECTED_BODIES
                  << " PAIRS " << EXPECTED.size() << '\n';
        for (const Expected& expected : EXPECTED) run_pair(expected, pool_cells);
        std::cout << "GLOBAL_VERDICT PASS THREE_EDGES FULL_30_POOL_9_BODY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4252_PRIMARY_ERROR " << error.what() << '\n';
        return 1;
    }
}
