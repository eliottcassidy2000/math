#include <filesystem>
// Detached direct-wall audit for THM-4283's exact 63-row carrier descent.
#define JOINT421_LITERAL_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/verify_joint421_literal_r670.cpp"
#undef JOINT421_LITERAL_LIBRARY_ONLY

namespace {
std::vector<u32> read_audit_masks(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask file");
    std::vector<u32> masks;
    std::set<u32> seen;
    std::string word;
    while (input >> word) {
        const u64 wide = std::stoull(word, nullptr, 16);
        require(wide < (UINT64_C(1) << 30), "mask outside pool");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && seen.insert(mask).second,
                "mask rank/duplicate changed");
        masks.push_back(mask);
    }
    return masks;
}
}

int main(int argc, char** argv) {
    try {
        require(argc == 3, "usage: literal BASE ADDITIONS");
        std::vector<u32> carrier = read_audit_masks(argv[1]);
        require(carrier.size() == 8951 &&
                    fnv_literal421(carrier) ==
                        UINT64_C(0x188f82ab9dd1695a),
                "base carrier changed");
        std::set<u32> seen(carrier.begin(), carrier.end());
        const std::vector<u32> additions = read_audit_masks(argv[2]);
        require(additions.size() == 45 &&
                    fnv_literal421(additions) ==
                        UINT64_C(0xec083b65cc8c34e3),
                "addition ledger changed");
        for (u32 mask : additions) {
            require(seen.insert(mask).second, "addition overlap");
            carrier.push_back(mask);
        }
        constexpr u32 repair = UINT32_C(0x014c9084);
        require(std::popcount(repair) == 8 && seen.insert(repair).second,
                "endpoint repair invalid/duplicate");
        carrier.push_back(repair);
        require(carrier.size() == 8997,
                "augmented carrier count changed");
        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "body universe changed");
        const std::array<std::pair<int,int>, 63> pairs{{
            {220,644}, {256,644}, {258,644}, {294,644},
            {366,644}, {416,644}, {512,644},
            {220,643}, {256,643}, {282,643}, {294,643}, {332,643},
            {338,643}, {366,643}, {372,643}, {384,643}, {412,643},
            {416,643}, {420,643}, {520,643}, {530,643},
            {220,642}, {256,642}, {258,642}, {294,642}, {384,642},
            {416,642}, {512,642}, {516,642}, {520,642},
            {220,641}, {256,641}, {282,641}, {294,641}, {366,641},
            {416,641}, {520,641},
            {220,640}, {248,640}, {282,640}, {294,640}, {335,640},
            {338,640}, {350,640}, {366,640}, {372,640}, {383,640},
            {416,640}, {420,640}, {440,640}, {520,640},
            {256,639}, {294,639}, {384,639}, {416,639},
            {100,638}, {282,638}, {294,638}, {366,638}, {416,638},
            {420,638}, {512,638}, {520,638}}};
        Fnv ledger;
        u64 total_failures = 0;
        for (const auto& [q,r] : pairs) {
            const Geometry geometry = build_joint_geometry(q, r);
            i64 repair_mass = 0;
            for (const Cell& cell : geometry.cells)
                if (cell.pair_safe && (cell.failed_pool & ~repair) == 0)
                    repair_mass += cell.width;
            const i128 repair_margin =
                static_cast<i128>(63) * repair_mass -
                static_cast<i128>(4) * geometry.grid;
            const std::vector<u32> active =
                active_literal421(q, r, carrier);
            const std::vector<u32> failures =
                failures_literal421(bodies, active);
            Fnv active_ledger;
            for (u32 mask : active) active_ledger.add(mask);
            Fnv failure_ledger;
            for (u32 body : failures) failure_ledger.add(body);
            total_failures += failures.size();
            ledger.add(q); ledger.add(r); ledger.add(active.size());
            ledger.add(active_ledger.state); ledger.add(failures.size());
            ledger.add(failure_ledger.state);
            ledger.add(static_cast<u64>(repair_margin));
            std::cout << "PAIR " << q << ',' << r
                      << " GRID " << decimal(geometry.grid)
                      << " ACTIVE " << active.size()
                      << " ACTIVE_FNV " << std::hex << active_ledger.state
                      << std::dec << " FAILURES " << failures.size()
                      << " FAILURE_FNV " << std::hex << failure_ledger.state
                      << std::dec << " REPAIR_MARGIN "
                      << decimal(repair_margin) << '\n';
        }
        std::cout << "LEDGER " << std::hex << ledger.state << std::dec
                  << " TOTAL_FAILURES " << total_failures << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR " << e.what() << '\n';
        return 1;
    }
}
