// Independent exact audit of the capacity packet's fixed-peel residual
// response bounds.  No signature-atlas or greedy source is imported.

#include "q50_response_common.hpp"

using namespace q50_common;

namespace {
struct CoverRow {
    u32 mask = 0;
    unsigned rank = 0;
    unsigned gain = 0;
    unsigned total = 0;
    unsigned residual = 0;
};

std::vector<u32> read_q50_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint589 failures");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    std::vector<u32> q50;
    std::set<u32> distinct;
    Fnv ledger;
    u64 q96 = 0;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0, r = 0;
        std::string token, extra;
        require(static_cast<bool>(fields >> q >> r >> token) &&
                    !(fields >> extra) && r == 589 && (q == 50 || q == 96),
                "failure row changed");
        const u32 body = parse_mask_agent(token);
        require(std::popcount(body) == 9, "failure rank changed");
        if (q == 96) {
            ++q96;
            continue;
        }
        require(distinct.insert(body).second, "duplicate q50 failure");
        q50.push_back(body); ledger.add(body);
    }
    require(input.eof() && q50.size() == 20025 && q96 == 11 &&
                ledger.state == UINT64_C(0xff421454f02d9099),
            "q50 failure identity changed");
    return q50;
}

std::vector<u32> read_packing(const std::filesystem::path& path,
                              const std::set<u32>& residual) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open capacity packing");
    std::string line;
    require(std::getline(input, line) &&
                line == "solution_index,core_index,body_hex",
            "capacity packing header changed");
    std::vector<u32> result;
    std::set<u32> distinct;
    Fnv ledger;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        std::size_t index = 0, core = 0;
        std::string token, extra;
        require(static_cast<bool>(fields >> index >> core >> token) &&
                    !(fields >> extra) && index == result.size(),
                "capacity packing row changed");
        const u32 body = parse_mask_agent(token);
        require(residual.contains(body) && distinct.insert(body).second,
                "capacity packing escaped residual");
        result.push_back(body);
        ledger.add(index); ledger.add(core); ledger.add(body);
    }
    require(input.eof() && result.size() == 19 &&
                ledger.state == UINT64_C(0x195708723b22fe7d),
            "capacity packing identity changed");
    return result;
}

std::vector<CoverRow> read_cover(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open capacity cover");
    std::string line;
    require(std::getline(input, line) &&
                line == "step,phase,mask_hex,rank,newly_covered,total_covered,residual",
            "capacity cover header changed");
    std::vector<CoverRow> result;
    std::set<u32> distinct;
    Fnv ledger;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        std::size_t step = 0;
        std::string phase, token, extra;
        CoverRow row;
        require(static_cast<bool>(fields >> step >> phase >> token >> row.rank >>
                                  row.gain >> row.total >> row.residual) &&
                    !(fields >> extra) && step == result.size() + 1 &&
                    (phase == "packing_seed" || phase == "greedy"),
                "capacity cover row changed");
        row.mask = parse_mask_agent(token);
        require((row.rank == 8 || row.rank == 9) &&
                    row.rank == std::popcount(row.mask) &&
                    distinct.insert(row.mask).second,
                "capacity cover mask changed");
        result.push_back(row);
        ledger.add(row.mask);
    }
    require(input.eof() && result.size() == 95 &&
                ledger.state == UINT64_C(0xfa44f9bfad76cfe7),
            "capacity cover identity changed");
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 4,
                "usage: audit_capacity_residual_19_95 FAILURES PACKING COVER");
        const std::vector<u32> failures = read_q50_failures(argv[1]);
        const Geometry geometry = build_geometry(50, 589);
        constexpr u32 inactive_peel = UINT32_C(0x0000b3a5);
        constexpr u32 active_peel = UINT32_C(0x0220932c);
        require(margin(geometry, inactive_peel).ticks < 0 &&
                    margin(geometry, active_peel).ticks >= 0,
                "fixed peel activity changed");

        std::vector<u32> residual;
        std::set<u32> residual_set;
        Fnv peeled_ledger, residual_ledger;
        u64 peeled = 0;
        for (std::size_t index = 0; index < failures.size(); ++index) {
            const u32 body = failures[index];
            if ((active_peel & body) == 0) {
                ++peeled;
                peeled_ledger.add(index); peeled_ledger.add(body);
            } else {
                residual.push_back(body); residual_set.insert(body);
                residual_ledger.add(index); residual_ledger.add(body);
            }
        }
        require(peeled == 891 && residual.size() == 19134 &&
                    peeled_ledger.state == UINT64_C(0xf766352d228791ed) &&
                    residual_ledger.state == UINT64_C(0x0e67c635fbc71d3b),
                "fixed peel residual changed");

        const std::vector<u32> packing = read_packing(argv[2], residual_set);
        const std::vector<CoverRow> cover = read_cover(argv[3]);

        std::vector<bool> covered(residual.size());
        unsigned total = 0;
        i128 minimum_ticks = 0;
        u32 minimum_mask = 0;
        for (std::size_t step = 0; step < cover.size(); ++step) {
            const CoverRow row = cover[step];
            const i128 ticks = margin(geometry, row.mask).ticks;
            require(ticks >= 0, "capacity cover contains inactive response");
            if (step == 0 || ticks < minimum_ticks) {
                minimum_ticks = ticks;
                minimum_mask = row.mask;
            }
            unsigned gain = 0;
            for (std::size_t index = 0; index < residual.size(); ++index) {
                if (!covered[index] && (row.mask & residual[index]) == 0) {
                    covered[index] = true;
                    ++gain;
                }
            }
            total += gain;
            require(gain == row.gain && total == row.total &&
                        row.residual == residual.size() - total,
                    "capacity greedy gain disagrees with direct replay");
        }
        require(total == residual.size() &&
                    std::ranges::all_of(covered, [](bool bit) { return bit; }) &&
                    minimum_ticks ==
                        static_cast<i128>(INT64_C(1550209054968)) &&
                    minimum_mask == UINT32_C(0x0a624049),
                "capacity 95-cover misses a residual body");

        std::array<u64, 2> universe{};
        std::array<u64, 2> active{};
        std::array<u64, 2> packing_responders{};
        std::array<unsigned, 2> maximum_load{};
        Fnv active_ledger, hit_ledger;
        for (unsigned rank : {8U, 9U}) {
            const std::size_t slot = rank - 8;
            const u32 limit = UINT32_C(1) << 30;
            for (u32 mask = (UINT32_C(1) << rank) - 1; mask < limit;
                 mask = next_combination(mask)) {
                ++universe[slot];
                if (margin(geometry, mask).ticks < 0) continue;
                ++active[slot];
                active_ledger.add(rank); active_ledger.add(mask);
                unsigned hits = 0;
                for (u32 body : packing) hits += (mask & body) == 0;
                require(hits <= 1,
                        "active response hits two capacity packing bodies");
                maximum_load[slot] = std::max(maximum_load[slot], hits);
                if (hits == 0) continue;
                ++packing_responders[slot];
                hit_ledger.add(rank); hit_ledger.add(mask); hit_ledger.add(hits);
            }
        }
        require(universe == std::array<u64, 2>{5852925, 14307150} &&
                    active == std::array<u64, 2>{480409, 4112383} &&
                    packing_responders == std::array<u64, 2>{172, 17049} &&
                    maximum_load == std::array<unsigned, 2>{1, 1} &&
                    active_ledger.state ==
                        UINT64_C(0xb79e52255c5b3522) &&
                    hit_ledger.state == UINT64_C(0x430fdb51e2ee1fa1),
                "capacity packing full-stream audit changed");

        std::cout << "ENDPOINT589_Q50_CAPACITY_19_95_INDEPENDENT_AUDIT_V1\n"
                  << "FAILURES " << failures.size() << " PEEL_COVERED "
                  << peeled << " RESIDUAL " << residual.size()
                  << " RESIDUAL_FNV " << std::hex << residual_ledger.state
                  << std::dec << '\n'
                  << "UPPER COVER_SIZE " << cover.size() << " COVER_FNV "
                  << std::hex << UINT64_C(0xfa44f9bfad76cfe7) << std::dec
                  << " COVERED " << total << " MIN_ACTIVITY_TICKS "
                  << static_cast<long long>(minimum_ticks) << " MIN_MASK "
                  << hex8(minimum_mask) << '\n'
                  << "LOWER PACKING_SIZE " << packing.size()
                  << " RESPONDERS_RANK8 " << packing_responders[0]
                  << " RESPONDERS_RANK9 " << packing_responders[1]
                  << " MAX_LOAD "
                  << std::max(maximum_load[0], maximum_load[1])
                  << " HIT_FNV " << std::hex << hit_ledger.state << std::dec
                  << '\n'
                  << "STREAM ACTIVE_RANK8 " << active[0] << " ACTIVE_RANK9 "
                  << active[1] << " ACTIVE_FNV " << std::hex
                  << active_ledger.state << std::dec << '\n'
                  << "CONSEQUENCE 19_LE_RESIDUAL_RESPONSE_COVER_NUMBER_LE_95\n"
                  << "SCOPE FIXED_Q96_PEEL_AND_Q50_RESPONSE_REPRESENTATION_ONLY\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CAPACITY_19_95_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}

