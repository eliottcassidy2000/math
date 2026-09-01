// Scratch activity-only scan for masks of the THM-4296 9,019-mask carrier
// that are strictly inactive on every one of the complete 72 current rows
// with endpoint at least 626.  Raw margins across pairs are deliberately not
// compared (MISTAKE-532); only exact signs are consumed.

#define main r632_hostile_base_main
#include "04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_detached_hostile_survivor.cpp"
#undef main

#include <array>
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
    for (u32 mask : kExchangeAddAgent) exchanged.push_back(mask);
    require(exchanged.size() == 9006, "exchange size changed");
    distinct = std::set<u32>(exchanged.begin(), exchanged.end());
    for (u32 mask : kRepairAAgent) {
        require(distinct.insert(mask).second, "A repair overlap");
        exchanged.push_back(mask);
    }
    for (u32 mask : kRepairBAgent) {
        require(distinct.insert(mask).second, "B repair overlap");
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
            "residual identity changed");
    std::vector<PairAgent> band;
    for (PairAgent pair : all)
        if (pair.r >= lower) band.push_back(pair);
    std::sort(band.begin(), band.end(), [](PairAgent left, PairAgent right) {
        if (left.r != right.r) return left.r > right.r;
        return left.q < right.q;
    });
    return band;
}

constexpr std::array<u32, 7> kPost632Repairs = {
    0x0010e125, 0x002ac4c0, 0x3882a082, 0x0041c325,
    0x08c28e40, 0x02008327, 0x0006e281};
constexpr std::array<u32, 3> kExchangeDelete626 = {
    0x0041212e, 0x0041303c, 0x0045106c};
constexpr std::array<u32, 3> kExchangeAdd626 = {
    0x210c1096, 0x02458324, 0x088a0ac1};

}  // namespace

#ifndef ENDPOINT626_EXCHANGE_MAIN
#define ENDPOINT626_EXCHANGE_MAIN main
#endif

int ENDPOINT626_EXCHANGE_MAIN(int argc, char** argv) {
    try {
        require(argc == 6,
                "usage: exchange_scan JOINT BASE8951 ADD45 SUFFIX9 RESIDUAL");
        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        std::vector<u32> carrier = build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(carrier.begin(), carrier.end());
        for (u32 mask : kPost632Repairs) {
            require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                        distinct.insert(mask).second,
                    "post-632 repair changed");
            carrier.push_back(mask);
        }
        require(carrier.size() == 9019 &&
                    masks_fnv_agent(carrier) == UINT64_C(0xd7f0e06e154e78c2),
                "final THM-4296 carrier identity changed");
        const std::vector<PairAgent> band = read_band_agent(argv[5], 626);
        require(band.size() == 72, "endpoint-626 prefix changed");

        std::vector<unsigned char> ever_active(carrier.size(), 0);
        std::vector<u64> active_rows(carrier.size(), 0);
        u64 equality_cells = 0;
        Fnv sign_ledger;
        for (PairAgent pair : band) {
            const Geometry geometry = build_geometry(pair.q, pair.r);
            for (std::size_t index = 0; index < carrier.size(); ++index) {
                const Margin value = margin(geometry, carrier[index]);
                const bool active = value.ticks >= 0;
                equality_cells += value.ticks == 0;
                ever_active[index] |= active;
                active_rows[index] += active;
                sign_ledger.add(carrier[index]);
                sign_ledger.add(active);
            }
        }

        std::vector<u32> inactive_all;
        std::vector<u32> inactive_nonjoint;
        std::array<u64, 73> row_histogram{};
        for (std::size_t index = 0; index < carrier.size(); ++index) {
            require(active_rows[index] <= band.size(), "active-row count escaped");
            ++row_histogram[active_rows[index]];
            if (ever_active[index]) continue;
            inactive_all.push_back(carrier[index]);
            if (!joint_set.contains(carrier[index]))
                inactive_nonjoint.push_back(carrier[index]);
        }
        std::sort(inactive_all.begin(), inactive_all.end());
        std::sort(inactive_nonjoint.begin(), inactive_nonjoint.end());
        Fnv all_ledger;
        for (u32 mask : inactive_all) all_ledger.add(mask);
        Fnv nonjoint_ledger;
        for (u32 mask : inactive_nonjoint) nonjoint_ledger.add(mask);
        const u64 inactive_rank8 = std::count_if(
            inactive_nonjoint.begin(), inactive_nonjoint.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        const u64 inactive_rank9 = inactive_nonjoint.size() - inactive_rank8;
        const std::set<u32> selected_delete(kExchangeDelete626.begin(),
                                            kExchangeDelete626.end());
        for (u32 mask : kExchangeDelete626)
            require(std::binary_search(inactive_nonjoint.begin(),
                                       inactive_nonjoint.end(), mask),
                    "selected deletion is not strictly inactive");
        std::vector<u32> exchanged;
        for (u32 mask : carrier)
            if (!selected_delete.contains(mask)) exchanged.push_back(mask);
        std::set<u32> exchanged_set(exchanged.begin(), exchanged.end());
        for (u32 mask : kExchangeAdd626) {
            require(std::popcount(mask) == 9 && exchanged_set.insert(mask).second,
                    "selected addition invalid or duplicate");
            exchanged.push_back(mask);
        }
        require(exchanged.size() == carrier.size(),
                "size-preserving exchange changed size");
        const u64 exchanged_rank8 = std::count_if(
            exchanged.begin(), exchanged.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        const u64 exchanged_rank9 = exchanged.size() - exchanged_rank8;

        std::cout << "ENDPOINT626_SIZE_PRESERVING_EXCHANGE_SCAN_V1\n"
                  << "CARRIER " << carrier.size() << " FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " PAIRS "
                  << band.size() << " SIGN_CELLS "
                  << carrier.size() * band.size() << " SIGN_FNV " << std::hex
                  << sign_ledger.state << std::dec << " EQUALITIES "
                  << equality_cells << '\n'
                  << "INACTIVE_ALL " << inactive_all.size() << " FNV "
                  << std::hex << all_ledger.state << std::dec
                  << " INACTIVE_NONJOINT " << inactive_nonjoint.size()
                  << " FNV " << std::hex << nonjoint_ledger.state << std::dec
                  << " RANK8 " << inactive_rank8 << " RANK9 "
                  << inactive_rank9 << '\n'
                  << "INACTIVE_NONJOINT_HEAD";
        for (std::size_t index = 0;
             index < std::min<std::size_t>(40, inactive_nonjoint.size()); ++index)
            std::cout << ' ' << hex8(inactive_nonjoint[index]);
        std::cout << '\n' << "ACTIVE_ROW_HISTOGRAM";
        for (std::size_t rows = 0; rows < row_histogram.size(); ++rows)
            if (row_histogram[rows] != 0)
                std::cout << ' ' << rows << ':' << row_histogram[rows];
        std::cout << '\n' << "SELECTED_EXCHANGE DELETE";
        for (u32 mask : kExchangeDelete626) std::cout << ' ' << hex8(mask);
        std::cout << " ADD";
        for (u32 mask : kExchangeAdd626) std::cout << ' ' << hex8(mask);
        std::cout << " RESULT_SIZE " << exchanged.size() << " RESULT_FNV "
                  << std::hex << masks_fnv_agent(exchanged) << std::dec
                  << " RANK8 " << exchanged_rank8 << " RANK9 "
                  << exchanged_rank9 << '\n'
                  << "SCOPE EXACT_LITERAL_SIGN_SCAN_COMPLETE_R_GE_626_PREFIX_"
                     "NO_CROSS_GRID_MARGIN_ORDER_NO_BELOW626_CLAIM_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_SIZE_PRESERVING_EXCHANGE_POOL\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT626_EXCHANGE_SCAN_ERROR " << error.what() << '\n';
        return 1;
    }
}
