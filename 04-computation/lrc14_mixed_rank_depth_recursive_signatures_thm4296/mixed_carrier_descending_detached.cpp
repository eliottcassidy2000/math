// Detached literal-wall descent of the carrier after the r=636 exchange and
// the six-mask mixed-rank r=632 repair.

#define main r632_detached_hostile_main
#include "r632_detached_hostile_survivor.cpp"
#undef main

#include <unordered_set>

namespace {

constexpr std::size_t kJointCountAgent = 421;
constexpr u64 kJointFnvAgent = UINT64_C(0x20d63dd42fe8150e);
constexpr std::size_t kBaseCountAgent = 8951;
constexpr u64 kBaseFnvAgent = UINT64_C(0x188f82ab9dd1695a);
constexpr std::size_t kAdd45CountAgent = 45;
constexpr u64 kAdd45FnvAgent = UINT64_C(0xec083b65cc8c34e3);
constexpr std::size_t kSuffixCountAgent = 9;
constexpr u64 kSuffixFnvAgent = UINT64_C(0x02b936529030e4bc);
constexpr std::size_t kResidualCountAgent = 22647;
constexpr u64 kResidualFnvAgent = UINT64_C(0xdf5374d4aca67677);
constexpr u64 kBodyCountAgent = UINT64_C(14307150);

constexpr std::array<u32, 14> kDeleteAgent = {
    0x00003e1a, 0x000132a3, 0x00017464, 0x00033388,
    0x000a16c2, 0x000f8118, 0x00142a1a, 0x00154348,
    0x00184ba0, 0x001aa260, 0x00202c2b, 0x002066a4,
    0x002b018a, 0x0030c2a2};
constexpr std::array<u32, 14> kExchangeAddAgent = {
    0x18468880, 0x080e8281, 0x22081017, 0x08422a82,
    0x004cac40, 0x19c04044, 0x00c08ec0, 0x10443016,
    0x01609124, 0x10413209, 0x01611640, 0x00606449,
    0x0128d084, 0x08806449};
constexpr std::array<u32, 3> kRepairAAgent = {
    0x2040c641, 0x00508325, 0x002a8641};
constexpr std::array<u32, 3> kRepairBAgent = {
    0x00619324, 0x201813a4, 0x21888126};

struct PairAgent {
    int q = 0;
    int r = 0;
};

struct PairAuditAgent {
    u64 active = 0;
    u64 active_fnv = 0;
    u64 active_joint = 0;
    u64 active_nonjoint = 0;
    u64 exposed = 0;
    u64 exposed_fnv = 0;
    u64 hit_incidences = 0;
    u64 minimum_hits = std::numeric_limits<u64>::max();
    u64 maximum_hits = 0;
    u64 failures = 0;
    u64 failure_fnv = 0;
    std::vector<u32> failure_bodies;
};

u32 parse_mask_agent(const std::string& token) {
    const u64 wide = std::stoull(token, nullptr, 16);
    require(wide < (UINT64_C(1) << 30), "mask escaped labels");
    return static_cast<u32>(wide);
}

std::vector<u32> read_masks_agent(const std::filesystem::path& path,
                                  std::size_t count, u64 expected_fnv,
                                  unsigned rank) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require(std::popcount(mask) == rank && distinct.insert(mask).second,
                "mask rank/distinctness changed");
        masks.push_back(mask);
        ledger.add(mask);
    }
    require(masks.size() == count && ledger.state == expected_fnv,
            "mask ledger identity changed");
    return masks;
}

std::vector<u32> build_mixed_carrier(const std::filesystem::path& base_path,
                                     const std::filesystem::path& add45_path,
                                     const std::filesystem::path& suffix_path) {
    std::vector<u32> carrier = read_masks_agent(
        base_path, kBaseCountAgent, kBaseFnvAgent, 8);
    const std::vector<u32> add45 = read_masks_agent(
        add45_path, kAdd45CountAgent, kAdd45FnvAgent, 8);
    const std::vector<u32> suffix = read_masks_agent(
        suffix_path, kSuffixCountAgent, kSuffixFnvAgent, 8);
    std::set<u32> distinct(carrier.begin(), carrier.end());
    for (u32 mask : add45) {
        require(distinct.insert(mask).second, "add45 overlap");
        carrier.push_back(mask);
    }
    require(distinct.insert(UINT32_C(0x014c9084)).second,
            "prior repair overlap");
    carrier.push_back(UINT32_C(0x014c9084));
    for (u32 mask : suffix) {
        require(distinct.insert(mask).second, "suffix overlap");
        carrier.push_back(mask);
    }
    const std::set<u32> deleted(kDeleteAgent.begin(), kDeleteAgent.end());
    std::vector<u32> exchanged;
    for (u32 mask : carrier)
        if (!deleted.contains(mask)) exchanged.push_back(mask);
    for (u32 mask : kExchangeAddAgent) {
        require(!distinct.contains(mask), "exchange addition overlap");
        exchanged.push_back(mask);
    }
    require(exchanged.size() == 9006, "exchange size changed");
    distinct = std::set<u32>(exchanged.begin(), exchanged.end());
    for (u32 mask : kRepairAAgent) {
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "A repair invalid");
        exchanged.push_back(mask);
    }
    for (u32 mask : kRepairBAgent) {
        require(std::popcount(mask) == 9 && distinct.insert(mask).second,
                "B repair invalid");
        exchanged.push_back(mask);
    }
    require(exchanged.size() == 9012, "mixed carrier size changed");
    return exchanged;
}

u64 masks_fnv_agent(const std::vector<u32>& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<PairAgent> read_band_agent(const std::filesystem::path& path,
                                       int lower) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open residual");
    std::vector<PairAgent> all;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos, "malformed residual row");
        PairAgent pair{std::stoi(line.substr(0, comma)),
                       std::stoi(line.substr(comma + 1))};
        require(pair.q > 0 && pair.q < pair.r, "invalid residual pair");
        if (!all.empty())
            require(std::tie(all.back().q, all.back().r) <
                        std::tie(pair.q, pair.r),
                    "residual order changed");
        all.push_back(pair);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(all.size() == kResidualCountAgent &&
                ledger.state == kResidualFnvAgent,
            "THM4287 residual changed");
    std::vector<PairAgent> band;
    for (PairAgent pair : all)
        if (pair.r >= lower) band.push_back(pair);
    std::sort(band.begin(), band.end(), [](PairAgent left, PairAgent right) {
        if (left.r != right.r) return left.r > right.r;
        return left.q < right.q;
    });
    return band;
}

PairAuditAgent audit_pair_agent(
    PairAgent pair, const std::vector<u32>& joint,
    const std::vector<u32>& carrier,
    const std::unordered_set<u32>& joint_set) {
    const Geometry geometry = build_geometry(pair.q, pair.r);
    std::vector<u32> active_joint;
    std::vector<u32> active_nonjoint;
    PairAuditAgent audit;
    Fnv active_ledger;
    for (u32 mask : joint)
        if (margin(geometry, mask).ticks >= 0) active_joint.push_back(mask);
    for (u32 mask : carrier) {
        if (margin(geometry, mask).ticks < 0) continue;
        ++audit.active;
        active_ledger.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    audit.active_fnv = active_ledger.state;
    audit.active_joint = active_joint.size();
    audit.active_nonjoint = active_nonjoint.size();
    require(audit.active == audit.active_joint + audit.active_nonjoint,
            "active partition changed");
    Fnv exposed_ledger;
    Fnv failure_ledger;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        bool joint_hit = false;
        for (u32 mask : active_joint)
            if ((mask & body) == 0) {
                joint_hit = true;
                break;
            }
        if (joint_hit) continue;
        ++audit.exposed;
        exposed_ledger.add(body);
        u64 hits = 0;
        for (u32 mask : active_nonjoint)
            if ((mask & body) == 0) ++hits;
        audit.hit_incidences += hits;
        if (hits == 0) {
            ++audit.failures;
            audit.failure_bodies.push_back(body);
            failure_ledger.add(body);
        } else {
            audit.minimum_hits = std::min(audit.minimum_hits, hits);
            audit.maximum_hits = std::max(audit.maximum_hits, hits);
        }
    }
    audit.exposed_fnv = exposed_ledger.state;
    audit.failure_fnv = failure_ledger.state;
    if (audit.exposed == 0) audit.minimum_hits = 0;
    return audit;
}

}  // namespace


int main(int argc, char** argv) {
    try {
        require(argc >= 8,
                "usage: mixed JOINT BASE8951 ADD45 SUFFIX9 RESIDUAL LOWER FAILURES [REPAIR_MASK ...]");
        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        std::vector<u32> carrier =
            build_mixed_carrier(argv[2], argv[3], argv[4]);
        for (int argument = 8; argument < argc; ++argument) {
            const u32 repair = parse_mask_agent(argv[argument]);
            require((std::popcount(repair) == 8 ||
                     std::popcount(repair) == 9) &&
                        std::find(carrier.begin(), carrier.end(), repair) ==
                            carrier.end(),
                    "repair invalid or duplicate");
            carrier.push_back(repair);
        }
        const int lower = std::stoi(argv[6]);
        const std::vector<PairAgent> band = read_band_agent(argv[5], lower);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        for (u32 mask : joint)
            require(std::find(carrier.begin(), carrier.end(), mask) !=
                        carrier.end(),
                    "joint mask absent from mixed carrier");
        std::ofstream failures_out(argv[7]);
        require(static_cast<bool>(failures_out), "cannot create failures");
        failures_out << "q,r,body_hex\n";
        const u64 rank8_count = std::count_if(
            carrier.begin(), carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        const u64 rank9_count = carrier.size() - rank8_count;
        std::cout << "MIXED_CARRIER_DESCENDING_DETACHED_V2\n"
                  << "CARRIER_SIZE " << carrier.size() << " FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " RANK8 "
                  << rank8_count
                  << " RANK9 " << rank9_count << " LOWER " << lower << " SELECTED "
                  << band.size() << '\n';
        std::size_t index = 0;
        u64 total_failures = 0;
        Fnv ledger;
        bool stopped = false;
        while (index < band.size() && !stopped) {
            const int endpoint = band[index].r;
            u64 layer_rows = 0;
            u64 layer_failures = 0;
            Fnv layer_ledger;
            while (index < band.size() && band[index].r == endpoint) {
                const PairAgent pair = band[index++];
                const PairAuditAgent audit =
                    audit_pair_agent(pair, joint, carrier, joint_set);
                for (u32 body : audit.failure_bodies)
                    failures_out << pair.q << ',' << pair.r << ','
                                 << hex8(body) << '\n';
                ++layer_rows;
                layer_failures += audit.failures;
                total_failures += audit.failures;
                layer_ledger.add(pair.q);
                layer_ledger.add(pair.r);
                ledger.add(pair.q);
                ledger.add(pair.r);
                ledger.add(audit.active);
                ledger.add(audit.active_fnv);
                ledger.add(audit.exposed);
                ledger.add(audit.exposed_fnv);
                ledger.add(audit.failures);
                ledger.add(audit.failure_fnv);
                std::cout << "PAIR " << pair.q << ',' << pair.r
                          << " ACTIVE " << audit.active << " ACTIVE_FNV "
                          << std::hex << audit.active_fnv << std::dec
                          << " ACTIVE_JOINT " << audit.active_joint
                          << " ACTIVE_NONJOINT " << audit.active_nonjoint
                          << " EXPOSED " << audit.exposed << " HIT_RANGE "
                          << audit.minimum_hits << ".." << audit.maximum_hits
                          << " HIT_INCIDENCES " << audit.hit_incidences
                          << " FAILURES " << audit.failures
                          << " FAILURE_FNV " << std::hex << audit.failure_fnv
                          << std::dec << '\n';
            }
            std::cout << "LAYER " << endpoint << " ROWS " << layer_rows
                      << " FNV " << std::hex << layer_ledger.state << std::dec
                      << " FAILURES " << layer_failures << '\n';
            stopped = layer_failures != 0;
        }
        require(failures_out.good(), "failure ledger write failed");
        std::cout << "COMPLETED_ROWS " << index << " TOTAL_FAILURES "
                  << total_failures << " STOPPED " << stopped
                  << " LEDGER_FNV " << std::hex << ledger.state << std::dec
                  << '\n'
                  << "SCOPE DETACHED_LITERAL_FIXED_POOL_THM4287_RESIDUAL_"
                     "MIXED_RANK8_9_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_MIXED_CARRIER_DESCENT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "MIXED_DESCENT_ERROR " << error.what() << '\n';
        return 1;
    }
}
