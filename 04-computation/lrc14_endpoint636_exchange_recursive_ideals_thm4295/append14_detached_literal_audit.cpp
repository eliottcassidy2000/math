// Direct literal-wall and all-body audit of the exact-minimum 14-mask repair
// for the THM-4287 endpoint-636 carrier failures. Activity is row-specific.
#define SINGLETON_LITERAL_LIBRARY_ONLY
#include "04-computation/lrc14_signature_response_congruence_thm4286/singleton_fibre_literal_verify.cpp"
#undef SINGLETON_LITERAL_LIBRARY_ONLY

#include <array>
#include <filesystem>
#include <fstream>
#include <set>
#include <unordered_map>

namespace {

constexpr std::size_t kInheritedCount = 8951;
constexpr std::size_t kAdditionCount = 45;
constexpr std::size_t kWitnessCount = 9;
constexpr std::size_t kAppendCount = 14;
constexpr u64 kInheritedFnv = UINT64_C(0x188f82ab9dd1695a);
constexpr u64 kAdditionsFnv = UINT64_C(0xec083b65cc8c34e3);
constexpr u64 kWitnessFnv = UINT64_C(0x02b936529030e4bc);
constexpr u64 kBase9006Fnv = UINT64_C(0xfdc1c57ae4dc1bb6);
constexpr u64 kAppendFnv = UINT64_C(0xcb6e6d8963cf54d0);
constexpr u64 kCarrier9020Fnv = UINT64_C(0x920651ae987fcdef);
constexpr u32 kPriorRepair = UINT32_C(0x014c9084);
constexpr std::size_t kFailureCount = 101;
constexpr u64 kBodyCount = UINT64_C(14307150);

struct Failure {
    int q = 0;
    int r = 0;
    u32 body = 0;
};

std::vector<u32> read_exact_masks(const std::filesystem::path& path,
                                  std::size_t count, u64 fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u64 wide = std::stoull(token, nullptr, 16);
        require(wide < (UINT64_C(1) << 30), "mask outside pool");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "mask rank/distinctness changed");
        masks.push_back(mask);
    }
    require(input.eof() && masks.size() == count && mask_fnv(masks) == fnv,
            "mask ledger identity changed");
    return masks;
}

std::vector<Failure> read_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "failure header changed");
    std::vector<Failure> rows;
    while (std::getline(input, line)) {
        const std::size_t a = line.find(',');
        const std::size_t b = line.find(',', a + 1);
        require(a != std::string::npos && b != std::string::npos &&
                    line.find(',', b + 1) == std::string::npos,
                "bad failure row");
        Failure row{std::stoi(line.substr(0, a)),
                    std::stoi(line.substr(a + 1, b - a - 1)),
                    static_cast<u32>(
                        std::stoul(line.substr(b + 1), nullptr, 16))};
        require(std::popcount(row.body) == 9 &&
                    ((row.q == 100 && row.r == 636) ||
                     (row.q == 256 && row.r == 636)),
                "failure row outside hostile pairs");
        rows.push_back(row);
    }
    require(rows.size() == kFailureCount, "failure count changed");
    std::vector<u32> a;
    std::vector<u32> b;
    for (const Failure& row : rows)
        (row.q == 100 ? a : b).push_back(row.body);
    require(a.size() == 64 && mask_fnv(a) == UINT64_C(0xd611500ea833ff83) &&
                b.size() == 37 &&
                mask_fnv(b) == UINT64_C(0xee7792a8a2fd51c9),
            "failure pair ledgers changed");
    return rows;
}

struct RowLiteral {
    int q = 0;
    int r = 0;
    LiteralPair pair;
    std::vector<i128> prefix;
};

RowLiteral build_row(int q, int r, const std::vector<Cell>& cells) {
    RowLiteral out;
    out.q = q;
    out.r = r;
    out.pair = build_literal_pair(q, r);
    std::vector<i64> walls(cells.size() + 1);
    walls.front() = cells.front().left;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        require(cells[i].left == walls[i], "pool wall discontinuity");
        walls[i + 1] = cells[i].right;
    }
    require(walls.front() == 0 && walls.back() == COMMON,
            "pool wall boundary changed");
    out.prefix.resize(walls.size());
    for (std::size_t i = 0; i < walls.size(); ++i)
        out.prefix[i] = prefix_at(
            out.pair,
            static_cast<i128>(walls[i]) * out.pair.pool_scale);
    return out;
}

struct Active {
    std::vector<u32> masks;
    u64 fnv = 0;
    u64 equalities = 0;
};

Active active_masks(const std::vector<u32>& masks,
                    const std::vector<Cell>& cells,
                    const RowLiteral& row) {
    Active out;
    FnvLocal ledger;
    for (u32 mask : masks) {
        const IndexedRepair repair = index_repair(mask, cells);
        const i128 margin = repair_margin(repair, row.prefix, row.pair.grid);
        out.equalities += margin == 0;
        if (margin >= 0) {
            out.masks.push_back(mask);
            ledger.add(mask);
        }
    }
    out.fnv = ledger.state;
    return out;
}

struct BodyAudit {
    u64 bodies = 0;
    u64 checks = 0;
    u64 max_checks = 0;
    std::vector<u32> failures;
};

BodyAudit audit_bodies_direct(const std::vector<u32>& active) {
    BodyAudit out;
    u32 body = (u32{1} << 9) - 1;
    while (body < (u32{1} << 30)) {
        ++out.bodies;
        u64 local = 0;
        bool covered = false;
        for (u32 mask : active) {
            ++local;
            if ((body & mask) == 0) {
                covered = true;
                break;
            }
        }
        out.checks += local;
        out.max_checks = std::max(out.max_checks, local);
        if (!covered) out.failures.push_back(body);
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(out.bodies == kBodyCount, "body universe changed");
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 6,
                "usage: append14-literal BASE8951 ADDITIONS45 WITNESS9 "
                "APPEND14 FAILURES101");
        auto carrier = read_exact_masks(argv[1], kInheritedCount,
                                        kInheritedFnv);
        const auto additions = read_exact_masks(argv[2], kAdditionCount,
                                                kAdditionsFnv);
        const auto witness = read_exact_masks(argv[3], kWitnessCount,
                                              kWitnessFnv);
        const auto append = read_exact_masks(argv[4], kAppendCount,
                                             kAppendFnv);
        const auto failures = read_failures(argv[5]);
        std::set<u32> distinct(carrier.begin(), carrier.end());
        for (u32 mask : additions) {
            require(distinct.insert(mask).second, "addition overlap");
            carrier.push_back(mask);
        }
        require(distinct.insert(kPriorRepair).second, "prior repair overlap");
        carrier.push_back(kPriorRepair);
        for (u32 mask : witness) {
            require(distinct.insert(mask).second, "witness overlap");
            carrier.push_back(mask);
        }
        require(carrier.size() == 9006 && mask_fnv(carrier) == kBase9006Fnv,
                "base9006 identity changed");
        for (u32 mask : append)
            require(distinct.insert(mask).second,
                    "append mask overlaps base carrier");

        const auto cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        const std::array<RowLiteral, 2> rows = {
            build_row(100, 636, cells), build_row(256, 636, cells)};
        std::array<Active, 2> active_base;
        std::array<Active, 2> active_append;
        for (std::size_t i = 0; i < rows.size(); ++i) {
            active_base[i] = active_masks(carrier, cells, rows[i]);
            active_append[i] = active_masks(append, cells, rows[i]);
        }
        require(active_base[0].masks.size() == 3556 &&
                    active_base[0].fnv == UINT64_C(0x33a5e283e7ebad6a) &&
                    active_base[1].masks.size() == 3239 &&
                    active_base[1].fnv == UINT64_C(0x218f5f1b44cc4b5) &&
                    active_append[0].masks.size() == 12 &&
                    active_append[1].masks.size() == 9,
                "literal active carrier identities changed");

        std::array<u64, 2> union_response{};
        FnvLocal response_ledger;
        std::cout << "THM4295_APPEND14_DETACHED_LITERAL_V1\n"
                  << "BASE9006_FNV " << std::hex << kBase9006Fnv
                  << " APPEND14_FNV " << kAppendFnv << " CARRIER9020_FNV "
                  << kCarrier9020Fnv << std::dec << '\n';
        for (u32 mask : append) {
            std::array<u64, 2> response{};
            std::array<i128, 2> margins{};
            std::array<bool, 2> active{};
            const IndexedRepair repair = index_repair(mask, cells);
            for (std::size_t row = 0; row < rows.size(); ++row) {
                margins[row] = repair_margin(
                    repair, rows[row].prefix, rows[row].pair.grid);
                active[row] = margins[row] >= 0;
            }
            for (std::size_t i = 0; i < failures.size(); ++i) {
                const std::size_t row = failures[i].q == 100 ? 0 : 1;
                if (active[row] && (mask & failures[i].body) == 0)
                    response[i >> 6] |= UINT64_C(1) << (i & 63);
            }
            union_response[0] |= response[0];
            union_response[1] |= response[1];
            response_ledger.add(mask);
            for (std::size_t row = 0; row < rows.size(); ++row) {
                add_i128(response_ledger, margins[row]);
                add_i128(response_ledger, rows[row].pair.grid);
                response_ledger.add(active[row]);
            }
            response_ledger.add(response[0]);
            response_ledger.add(response[1]);
            std::cout << "MASK " << std::hex << std::setw(8)
                      << std::setfill('0') << mask << std::dec
                      << std::setfill(' ') << " A100 " << active[0]
                      << " M100 " << decimal(margins[0]) << '/'
                      << decimal(rows[0].pair.grid) << " A256 " << active[1]
                      << " M256 " << decimal(margins[1]) << '/'
                      << decimal(rows[1].pair.grid) << " COVER100 "
                      << std::popcount(response[0]) << " COVER256 "
                      << std::popcount(response[1]) << " WORDS " << std::hex
                      << std::setw(16) << std::setfill('0') << response[1]
                      << ',' << std::setw(16) << response[0] << std::dec
                      << std::setfill(' ') << '\n';
        }
        const std::array<u64, 2> full = {
            UINT64_MAX, (UINT64_C(1) << 37) - 1};
        require(union_response == full, "append14 response union incomplete");

        std::vector<u32> carrier9020 = carrier;
        carrier9020.insert(carrier9020.end(), append.begin(), append.end());
        require(carrier9020.size() == 9020 &&
                    mask_fnv(carrier9020) == kCarrier9020Fnv,
                "carrier9020 identity changed");
        for (std::size_t row = 0; row < rows.size(); ++row) {
            const BodyAudit before = audit_bodies_direct(active_base[row].masks);
            std::vector<u32> expected;
            for (const Failure& failure : failures)
                if ((row == 0 && failure.q == 100) ||
                    (row == 1 && failure.q == 256))
                    expected.push_back(failure.body);
            require(before.failures == expected,
                    "literal pre-append failure ledger changed");
            std::vector<u32> active9020 = active_base[row].masks;
            active9020.insert(active9020.end(), active_append[row].masks.begin(),
                              active_append[row].masks.end());
            FnvLocal active_ledger;
            for (u32 mask : active9020) active_ledger.add(mask);
            const BodyAudit after = audit_bodies_direct(active9020);
            require(after.failures.empty(),
                    "literal post-append carrier still fails");
            std::cout << "PAIR " << rows[row].q << ',' << rows[row].r
                      << " GRID " << decimal(rows[row].pair.grid)
                      << " BASE_ACTIVE " << active_base[row].masks.size()
                      << " BASE_FNV " << std::hex << active_base[row].fnv
                      << std::dec << " APPEND_ACTIVE "
                      << active_append[row].masks.size() << " ACTIVE9020 "
                      << active9020.size() << " ACTIVE9020_FNV " << std::hex
                      << active_ledger.state << std::dec
                      << " EQUALITIES "
                      << active_base[row].equalities +
                             active_append[row].equalities
                      << " BEFORE_FAILURES " << before.failures.size()
                      << " BEFORE_CHECKS " << before.checks
                      << " BEFORE_MAX " << before.max_checks
                      << " AFTER_FAILURES " << after.failures.size()
                      << " AFTER_CHECKS " << after.checks
                      << " AFTER_MAX " << after.max_checks << '\n';
        }
        std::cout << "RESPONSE_UNION " << std::hex << union_response[1] << ','
                  << union_response[0] << " RESPONSE_LEDGER_FNV "
                  << response_ledger.state << std::dec << '\n'
                  << "VERDICT PASS ROW_SPECIFIC_ACTIVITY_AND_DIRECT_BODY_COVER\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "APPEND14_LITERAL_ERROR " << error.what() << '\n';
        return 1;
    }
}
