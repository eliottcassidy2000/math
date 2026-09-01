// Full-universe exact audit for THM-4302.  It enumerates every rank-eight
// and rank-nine mask, reconstructs the complete response quotient of the 24
// THM-4300 failures at (210,596), verifies exact primal/dual certificates,
// and audits a strictly inactive size-preserving exchange on all r>=596 rows.

#define ENDPOINT626_EXCHANGE_MAIN thm4302_hidden_main
#include "../lrc14_size_preserving_response_staircase_thm4300/endpoint_exchange_primitives.cpp"
#undef ENDPOINT626_EXCHANGE_MAIN

#include <fstream>
#include <unordered_map>

namespace {

constexpr u64 kAugmentedFnv = UINT64_C(0x55e8588798885ae5);
constexpr u64 kRepairFnv = UINT64_C(0x64ce5f9d1ec8c4c2);
constexpr u64 kFailureFnv = UINT64_C(0x3dbd5b39673070ff);
constexpr u64 kResponder8Fnv = UINT64_C(0x2dddbe3405491cdd);
constexpr u64 kResponder9Fnv = UINT64_C(0x2202e93c739926df);
constexpr u64 kInactiveFnv = UINT64_C(0xfa143e58f59119f8);
constexpr u64 kDeleteFnv = UINT64_C(0x9240b264ab65aa62);
constexpr u64 kAdditionFnv = UINT64_C(0xdc0eebaebf688c65);
constexpr u64 kPrefixSignFnv = UINT64_C(0xaf07f32bd2453cd0);
constexpr u64 kExchangedFnv = UINT64_C(0x892fef44a9e6b37e);

constexpr std::array<u32, 4> kAdditionMasks = {
    0x185f0200, 0x28070a88, 0x02710e02, 0x0010c574};
constexpr std::array<u32, 4> kAdditionResponses = {
    0x00d3c000, 0x00002ec3, 0x000f4180, 0x00e0343c};
constexpr std::array<u32, 6> kRank8ControlMasks = {
    0x2010c125, 0x12980602, 0x183a0280,
    0x20014a89, 0x082244c8, 0x08720c08};

struct ResponseStats {
    u64 count8 = 0;
    u64 count9 = 0;
    u32 least8 = UINT32_MAX;
    u32 least9 = UINT32_MAX;
};

std::vector<u32> read_mixed(const std::filesystem::path& path,
                            std::size_t expected, u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mixed-mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "mixed rank/distinctness changed");
        masks.push_back(mask);
        ledger.add(mask);
    }
    require(input.eof() && masks.size() == expected &&
                ledger.state == expected_fnv,
            "mixed-mask identity changed");
    return masks;
}

std::vector<u32> read_failures596(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    std::vector<u32> bodies;
    Fnv ledger;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        int q = 0;
        int r = 0;
        std::string token;
        row >> q >> r >> token;
        require(row && q == 210 && r == 596, "failure pair changed");
        const u32 body = parse_mask_agent(token);
        require(std::popcount(body) == 9, "failure body rank changed");
        bodies.push_back(body);
        ledger.add(body);
    }
    require(bodies.size() == 24 &&
                std::set<u32>(bodies.begin(), bodies.end()).size() == 24 &&
                ledger.state == kFailureFnv,
            "failure universe changed");
    return bodies;
}

u32 response(u32 mask, const std::vector<u32>& failures) {
    u32 answer = 0;
    for (std::size_t i = 0; i < failures.size(); ++i)
        if ((mask & failures[i]) == 0) answer |= u32{1} << i;
    return answer;
}

u64 masks_fnv(const auto& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 10,
                "usage: audit JOINT BASE8951 ADD45 SUFFIX9 RESIDUAL "
                "REPAIRS76 FAILURES24 ADDITIONS4 DELETE73");
        const auto joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const auto repairs = read_mixed(argv[6], 76, kRepairFnv);
        const auto failures = read_failures596(argv[7]);
        const auto serialized_additions = read_mixed(argv[8], 4, kAdditionFnv);
        const auto serialized_delete = read_mixed(argv[9], 73, kDeleteFnv);
        require(std::equal(serialized_additions.begin(),
                           serialized_additions.end(),
                           kAdditionMasks.begin()),
                "addition ledger changed");

        std::vector<u32> augmented =
            build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(augmented.begin(), augmented.end());
        for (u32 repair : repairs) {
            require(distinct.insert(repair).second, "repair overlap");
            augmented.push_back(repair);
        }
        require(augmented.size() == 9088 &&
                    masks_fnv_agent(augmented) == kAugmentedFnv,
                "augmented carrier changed");
        const std::unordered_set<u32> augmented_set(
            augmented.begin(), augmented.end());

        const Geometry local = build_geometry(210, 596);
        std::unordered_map<u32, ResponseStats> atlas;
        u64 active8 = 0;
        u64 active9 = 0;
        u64 responders8 = 0;
        u64 responders9 = 0;
        Fnv responder8_ledger;
        Fnv responder9_ledger;
        auto scan_rank = [&](unsigned rank) {
            const u32 limit = u32{1} << 30;
            u64 universe = 0;
            for (u32 mask = (u32{1} << rank) - 1; mask < limit;
                 mask = next_combination(mask)) {
                ++universe;
                if (margin(local, mask).ticks < 0) continue;
                if (rank == 8) ++active8; else ++active9;
                const u32 reply = response(mask, failures);
                if (reply == 0) continue;
                require(!augmented_set.contains(mask),
                        "existing carrier responds to a recorded failure");
                ResponseStats& stats = atlas[reply];
                if (rank == 8) {
                    ++stats.count8;
                    stats.least8 = std::min(stats.least8, mask);
                    ++responders8;
                    responder8_ledger.add(mask);
                } else {
                    ++stats.count9;
                    stats.least9 = std::min(stats.least9, mask);
                    ++responders9;
                    responder9_ledger.add(mask);
                }
            }
            require(universe == (rank == 8 ? UINT64_C(5852925)
                                           : UINT64_C(14307150)),
                    "mask universe changed");
        };
        scan_rank(8);
        scan_rank(9);
        require(active8 == 1267366 && active9 == 7100734 &&
                    responders8 == 9090 && responders9 == 138019 &&
                    responder8_ledger.state == kResponder8Fnv &&
                    responder9_ledger.state == kResponder9Fnv &&
                    atlas.size() == 718,
                "response atlas changed");

        std::vector<u32> types;
        for (const auto& [reply, stats] : atlas) types.push_back(reply);
        std::sort(types.begin(), types.end(), [](u32 left, u32 right) {
            if (std::popcount(left) != std::popcount(right))
                return std::popcount(left) > std::popcount(right);
            return left < right;
        });
        std::vector<u32> maximal;
        for (u32 type : types) {
            bool dominated = false;
            for (u32 keep : maximal)
                if ((type | keep) == keep) {
                    dominated = true;
                    break;
                }
            if (!dominated) maximal.push_back(type);
        }
        require(maximal.size() == 82, "maximal response antichain changed");

        // Half-integral dual: weight 1/2 on seven obligations.  Every
        // realized response has scaled load at most two, so the fractional
        // value is 7/2 and every integral cover has size at least four.
        constexpr u32 mixed_dual_support = 0x00288146;
        std::array<u64, 3> mixed_load_hist{};
        for (u32 type : types) {
            const unsigned load = std::popcount(type & mixed_dual_support);
            require(load <= 2, "mixed dual overload");
            ++mixed_load_hist[load];
        }
        require(mixed_load_hist == std::array<u64, 3>{255, 355, 108},
                "mixed dual histogram changed");

        u32 covered = 0;
        for (std::size_t i = 0; i < kAdditionMasks.size(); ++i) {
            const u32 mask = kAdditionMasks[i];
            require(std::popcount(mask) == 9 &&
                        margin(local, mask).ticks > 0 &&
                        response(mask, failures) == kAdditionResponses[i] &&
                        !augmented_set.contains(mask),
                    "four-mask cover changed");
            covered |= kAdditionResponses[i];
        }
        require(covered == 0x00ffffff &&
                    masks_fnv(kAdditionMasks) == kAdditionFnv,
                "four-mask response cover changed");

        // Hostile control: rank eight alone has fractional dual value 11/2
        // and an explicit six-mask cover, hence exact minimum six.
        constexpr u32 rank8_half_support = 0x00d88143;
        constexpr u32 rank8_double_support = 0x00000004;
        std::array<u64, 3> rank8_load_hist{};
        u64 rank8_types = 0;
        for (u32 type : types) {
            if (atlas.at(type).count8 == 0) continue;
            ++rank8_types;
            const unsigned load =
                std::popcount(type & rank8_half_support) +
                2 * std::popcount(type & rank8_double_support);
            require(load <= 2, "rank-eight dual overload");
            ++rank8_load_hist[load];
        }
        require(rank8_types == 220 &&
                    rank8_load_hist == std::array<u64, 3>{78, 104, 38},
                "rank-eight dual changed");
        u32 rank8_cover = 0;
        for (u32 mask : kRank8ControlMasks) {
            require(std::popcount(mask) == 8 && margin(local, mask).ticks > 0,
                    "rank-eight control inactive");
            rank8_cover |= response(mask, failures);
        }
        require(rank8_cover == 0x00ffffff,
                "rank-eight six-mask control fails cover");

        const auto band = read_band_agent(argv[5], 596);
        require(band.size() == 363, "endpoint-at-least-596 prefix changed");
        std::vector<Geometry> geometries;
        geometries.reserve(band.size());
        for (PairAgent pair : band)
            geometries.push_back(build_geometry(pair.q, pair.r));
        std::vector<u32> inactive;
        u64 equalities = 0;
        Fnv sign_ledger;
        for (u32 mask : augmented) {
            bool ever_active = false;
            for (const Geometry& geometry : geometries) {
                const Margin exact = margin(geometry, mask);
                const bool active = exact.ticks >= 0;
                ever_active |= active;
                equalities += exact.ticks == 0;
                sign_ledger.add(mask);
                sign_ledger.add(active);
            }
            if (!ever_active) {
                require(!joint_set.contains(mask) && std::popcount(mask) == 8,
                        "common-inactive mask type changed");
                inactive.push_back(mask);
            }
        }
        std::sort(inactive.begin(), inactive.end());
        require(inactive.size() == 75 && masks_fnv(inactive) == kInactiveFnv &&
                    equalities == 0 && sign_ledger.state == kPrefixSignFnv,
                "prefix sign audit changed");
        require(std::equal(serialized_delete.begin(), serialized_delete.end(),
                           inactive.begin()),
                "delete ledger is not the lexicographic inactive prefix");

        const std::set<u32> delete_set(
            serialized_delete.begin(), serialized_delete.end());
        std::vector<u32> exchanged;
        for (u32 mask : augmented)
            if (!delete_set.contains(mask)) exchanged.push_back(mask);
        for (u32 mask : kAdditionMasks) exchanged.push_back(mask);
        const u64 rank8 = std::count_if(
            exchanged.begin(), exchanged.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(exchanged.size() == 9019 && rank8 == 8961 &&
                    masks_fnv_agent(exchanged) == kExchangedFnv,
                "size-preserving exchange changed");

        std::cout << "THM4302_ENDPOINT596_RESPONSE_EXCHANGE_AUDIT_V1\n"
                  << "FAILURES 24 FNV " << std::hex << kFailureFnv
                  << std::dec << " PAIR 210,596\n"
                  << "RANK8_UNIVERSE 5852925 ACTIVE " << active8
                  << " RESPONDERS " << responders8 << " FNV " << std::hex
                  << responder8_ledger.state << std::dec << '\n'
                  << "RANK9_UNIVERSE 14307150 ACTIVE " << active9
                  << " RESPONDERS " << responders9 << " FNV " << std::hex
                  << responder9_ledger.state << std::dec << '\n'
                  << "RESPONSE_TYPES " << atlas.size() << " MAXIMAL "
                  << maximal.size() << '\n'
                  << "MIXED_DUAL HALF_SUPPORT 00288146 LOAD_HIST 0:255 1:355 "
                     "2:108 VALUE 7/2 INTEGER_LOWER 4\n"
                  << "FOUR_MASK_COVER";
        for (std::size_t i = 0; i < kAdditionMasks.size(); ++i)
            std::cout << ' ' << hex8(kAdditionMasks[i]) << ':'
                      << hex8(kAdditionResponses[i]);
        std::cout << " FNV " << std::hex << kAdditionFnv << std::dec
                  << " UNION 00ffffff EXACT_MINIMUM 4\n"
                  << "RANK8_HOSTILE TYPES 220 DUAL_VALUE 11/2 INTEGER_LOWER 6 "
                     "SIX_MASK_COVER YES EXACT_MINIMUM 6\n"
                  << "PREFIX_ROWS " << band.size() << " SIGN_CELLS "
                  << augmented.size() * band.size() << " SIGN_FNV " << std::hex
                  << sign_ledger.state << std::dec << " EQUALITIES " << equalities
                  << '\n'
                  << "COMMON_INACTIVE_NONJOINT " << inactive.size() << " FNV "
                  << std::hex << kInactiveFnv << std::dec
                  << " RANK8 75 RANK9 0\n"
                  << "DELETE 73 FNV " << std::hex << kDeleteFnv
                  << " ADD 4 FNV " << kAdditionFnv << std::dec << '\n'
                  << "EXCHANGED 9019 FNV " << std::hex << kExchangedFnv
                  << std::dec << " RANK8 " << rank8 << " RANK9 "
                  << exchanged.size() - rank8 << " UNUSED_INACTIVE 2\n"
                  << "TRANSFER deleting masks inactive at every prefix row "
                     "preserves every old active subfamily; the four additions "
                     "cover the complete recorded boundary failure set\n"
                  << "SCOPE FINITE_EXACT_FIXED_POOL_CARRIER_CERTIFICATE_ONLY_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_MINIMUM_FOUR_AND_SIZE_PRESERVING_"
                     "ENDPOINT596_EXCHANGE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4302_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
