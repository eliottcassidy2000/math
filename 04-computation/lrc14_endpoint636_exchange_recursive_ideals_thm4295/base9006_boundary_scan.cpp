// Exact all-layer boundary scan for THM-4287's inherited 9,006-mask carrier.

#include "04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/carrier_scan_support.cpp"

#include <set>
#include <tuple>

namespace {

constexpr u32 kRepair644 = UINT32_C(0x014c9084);
constexpr u64 kCarrier9006Fnv = UINT64_C(0xfdc1c57ae4dc1bb6);
constexpr std::size_t kOld9Count = 9;
constexpr u64 kOld9Fnv = UINT64_C(0x02b936529030e4bc);

std::vector<Pair> read_live_band(const std::filesystem::path& path,
                                 int lower_endpoint) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open live residual");
    std::vector<Pair> all, band;
    FnvLocal ledger;
    std::string line;
    while (std::getline(in, line)) {
        require(!line.empty(), "blank live residual row");
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed live residual row");
        std::size_t uq = 0, ur = 0;
        Pair p{std::stoi(line.substr(0, comma), &uq),
               std::stoi(line.substr(comma + 1), &ur)};
        require(uq == comma && ur == line.size() - comma - 1 && p.q > 0 &&
                    p.q < p.r &&
                    (all.empty() || std::tie(all.back().q, all.back().r) <
                                        std::tie(p.q, p.r)),
                "invalid live residual row");
        all.push_back(p);
        ledger.add(p.q);
        ledger.add(p.r);
        if (p.r >= lower_endpoint) band.push_back(p);
    }
    require(all.size() == 22647 &&
                ledger.state == UINT64_C(0xdf5374d4aca67677),
            "THM-4287 live residual identity changed");
    std::sort(band.begin(), band.end(), [](Pair a, Pair b) {
        return a.r != b.r ? a.r > b.r : a.q < b.q;
    });
    return band;
}

PairAudit audit_fast(const std::vector<Cell>& cells,
                     const std::vector<u32>& joint,
                     const std::vector<u32>& carrier,
                     const std::unordered_set<u32>& joint_set,
                     Pair pair) {
    const ActiveUniverse active = build_active_universe(cells, pair.q, pair.r);
    std::vector<u32> active_joint, active_nonjoint;
    PairAudit audit;
    FnvLocal active_carrier_ledger, inactive_joint_ledger;
    for (std::size_t i = 0; i < joint.size(); ++i) {
        const u32 mask = joint[i];
        if (active.active[colex_rank8_local(mask)]) {
            active_joint.push_back(mask);
        } else {
            ++audit.inactive_joint;
            inactive_joint_ledger.add(i);
            inactive_joint_ledger.add(mask);
        }
    }
    for (u32 mask : carrier) {
        if (!active.active[colex_rank8_local(mask)]) continue;
        ++audit.active_carrier;
        active_carrier_ledger.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    audit.active_joint = active_joint.size();
    audit.active_nonjoint = active_nonjoint.size();
    FnvLocal exposed_ledger, response_ledger;
    u32 body = (u32{1} << 9) - 1;
    for (std::size_t ordinal = 0; ordinal < kBodyCount; ++ordinal) {
        bool joint_hit = false;
        for (u32 mask : active_joint)
            if ((body & mask) == 0) {
                joint_hit = true;
                break;
            }
        if (!joint_hit) {
            ++audit.exposed;
            exposed_ledger.add(body);
            u64 hits = 0;
            u32 least = 0;
            for (u32 mask : active_nonjoint) {
                if ((body & mask) != 0) continue;
                if (hits == 0 || mask < least) least = mask;
                ++hits;
            }
            response_ledger.add(body);
            response_ledger.add(hits);
            response_ledger.add(least);
            audit.hit_incidences += hits;
            if (hits == 0) {
                ++audit.failures;
                audit.failure_bodies.push_back(body);
            } else {
                if (hits < audit.minimum_hits) {
                    audit.minimum_hits = hits;
                    audit.minimum_hit_body = body;
                }
                if (hits > audit.maximum_hits) {
                    audit.maximum_hits = hits;
                    audit.maximum_hit_body = body;
                }
            }
        }
        if (ordinal + 1 < kBodyCount) body = next_same_popcount(body);
    }
    require(body == UINT32_C(0x3fe00000), "body enumeration changed");
    audit.active_carrier_fnv = active_carrier_ledger.state;
    audit.inactive_joint_fnv = inactive_joint_ledger.state;
    audit.exposed_fnv = exposed_ledger.state;
    audit.response_fnv = response_ledger.state;
    if (audit.exposed == 0) audit.minimum_hits = 0;
    FnvLocal failure_ledger;
    for (u32 body_mask : audit.failure_bodies) failure_ledger.add(body_mask);
    audit.failure_fnv = failure_ledger.state;
    return audit;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 8,
                "usage: base9006-scan JOINT BASE8951 ADD45 OLD9 LIVE22647 "
                "LOWER FAILURES");
        init_choose8_local();
        const std::vector<u32> joint =
            read_masks(argv[1], kJointCount, kJointFnv);
        std::vector<u32> carrier =
            read_masks(argv[2], kCarrierCount, kCarrierFnv);
        const std::vector<u32> additions = read_additions(argv[3]);
        const std::vector<u32> old9 =
            read_masks(argv[4], kOld9Count, kOld9Fnv);
        std::set<u32> distinct(carrier.begin(), carrier.end());
        for (u32 m : additions) {
            require(distinct.insert(m).second, "duplicate add45");
            carrier.push_back(m);
        }
        require(distinct.insert(kRepair644).second, "duplicate prior repair");
        carrier.push_back(kRepair644);
        for (u32 m : old9) {
            require(distinct.insert(m).second, "duplicate old9");
            carrier.push_back(m);
        }
        require(carrier.size() == 9006 &&
                    mask_fnv(carrier) == kCarrier9006Fnv,
                "base repaired carrier identity changed");
        const auto band = read_live_band(argv[5], std::stoi(argv[6]));
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<Cell> cells = build_pool_cells();
        std::ofstream failures(argv[7]);
        require(static_cast<bool>(failures), "cannot create failure ledger");
        failures << "q,r,body_hex\n";
        std::cout << "THM4295_BASE9006_BOUNDARY_SCAN_V1\n"
                  << "CARRIER " << carrier.size() << " FNV " << std::hex
                  << mask_fnv(carrier) << std::dec << " SELECTED "
                  << band.size() << '\n';
        std::size_t index = 0;
        u64 total_failures = 0;
        while (index < band.size()) {
            const int endpoint = band[index].r;
            u64 layer_failures = 0;
            std::size_t layer_rows = 0;
            FnvLocal layer_ledger;
            while (index < band.size() && band[index].r == endpoint) {
                const Pair p = band[index++];
                const PairAudit a = audit_fast(
                    cells, joint, carrier, joint_set, p);
                for (u32 body : a.failure_bodies)
                    failures << p.q << ',' << p.r << ',' << std::hex
                             << std::setw(8) << std::setfill('0') << body
                             << std::dec << std::setfill(' ') << '\n';
                layer_failures += a.failures;
                total_failures += a.failures;
                ++layer_rows;
                layer_ledger.add(p.q);
                layer_ledger.add(p.r);
                std::cout << "PAIR " << p.q << ',' << p.r
                          << " ACTIVE " << a.active_carrier
                          << " EXPOSED " << a.exposed
                          << " FAILURES " << a.failures
                          << " FAILURE_FNV " << std::hex << a.failure_fnv
                          << std::dec << '\n';
            }
            std::cout << "LAYER " << endpoint << " ROWS " << layer_rows
                      << " FNV " << std::hex << layer_ledger.state << std::dec
                      << " FAILURES " << layer_failures << '\n';
        }
        require(failures.good(), "failure ledger write failed");
        std::cout << "TOTAL_FAILURES " << total_failures << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "THM4295_BASE9006_BOUNDARY_ERROR " << e.what() << '\n';
        return 1;
    }
}
