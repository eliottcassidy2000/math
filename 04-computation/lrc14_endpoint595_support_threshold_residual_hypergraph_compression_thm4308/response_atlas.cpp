// Exact THM-4307 audit for the three THM-4303 endpoint-595 carrier failures.
// It constructs the complete pair-tagged response quotient of all active
// rank-eight and rank-nine masks and audits strict common inactivity on the
// complete endpoint-at-least-595 prefix.  No claim beyond the frozen pool is
// intended.

#define ENDPOINT617_RAW_VERIFY_MAIN endpoint595_response_hidden_main
#include "04-computation/lrc14_size_preserving_response_staircase_thm4300/endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

#include <fstream>
#include <unordered_map>

namespace {

constexpr u64 kRepairFnv = UINT64_C(0x64ce5f9d1ec8c4c2);
constexpr u64 kAdditionFnv = UINT64_C(0xdc0eebaebf688c65);
constexpr u64 kDeleteFnv = UINT64_C(0x9240b264ab65aa62);
constexpr u64 kAugmentedFnv = UINT64_C(0x55e8588798885ae5);
constexpr u64 kCarrierFnv = UINT64_C(0x892fef44a9e6b37e);
constexpr u64 kFailurePairFnv = UINT64_C(0xf3d7f95fc38e7b49);
constexpr u64 kOldUnionFnv = UINT64_C(0x11414a33ab91fef6);
constexpr u64 kResponder8Fnv = UINT64_C(0xcc926b13c719225d);
constexpr u64 kResponder9Fnv = UINT64_C(0x8ff584ab870b72a1);
constexpr u64 kTypeLedgerFnv = UINT64_C(0xa5f5ebcdeb03ad34);
constexpr u64 kAugmentedPrefixSignFnv = UINT64_C(0x8608ee11775be3fc);
constexpr u64 kCarrierPrefixSignFnv = UINT64_C(0x0853edf378527886);
constexpr u64 kCarrierInactiveFnv = UINT64_C(0x3b5ca775eedae38b);

struct Response {
    // 96-row obligations occupy words 0--1; 100-row and 210-row obligations
    // occupy words 2 and 3.  Unused high bits are always zero.
    std::array<u64, 4> word{};
    bool operator==(const Response&) const = default;
};

struct ResponseHash {
    std::size_t operator()(const Response& response) const noexcept {
        u64 h = UINT64_C(0x9e3779b97f4a7c15);
        for (u64 x : response.word) {
            x ^= x >> 30;
            x *= UINT64_C(0xbf58476d1ce4e5b9);
            x ^= x >> 27;
            x *= UINT64_C(0x94d049bb133111eb);
            x ^= x >> 31;
            h ^= x + UINT64_C(0x9e3779b97f4a7c15) + (h << 6) + (h >> 2);
        }
        return static_cast<std::size_t>(h);
    }
};

struct ResponseStats {
    u64 count8 = 0;
    u64 count9 = 0;
    u32 least8 = UINT32_MAX;
    u32 least9 = UINT32_MAX;
};

struct FailureRow {
    PairAgent pair;
    std::vector<u32> bodies;
};

bool empty(const Response& response) {
    return std::all_of(response.word.begin(), response.word.end(),
                       [](u64 value) { return value == 0; });
}

unsigned weight(const Response& response) {
    unsigned answer = 0;
    for (u64 value : response.word) answer += std::popcount(value);
    return answer;
}

std::vector<FailureRow> read_failures595(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint-595 failures");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    std::vector<FailureRow> rows;
    Fnv pair_ledger;
    std::set<std::tuple<int, int, u32>> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0;
        int r = 0;
        std::string token;
        fields >> q >> r >> token;
        require(fields && q > 0 && q < r, "malformed failure row");
        const u32 body = parse_mask_agent(token);
        require(std::popcount(body) == 9 &&
                    distinct.insert({q, r, body}).second,
                "failure rank/distinctness changed");
        if (rows.empty() || rows.back().pair.q != q || rows.back().pair.r != r)
            rows.push_back(FailureRow{PairAgent{q, r}, {}});
        rows.back().bodies.push_back(body);
        pair_ledger.add(q);
        pair_ledger.add(r);
        pair_ledger.add(body);
    }
    require(input.eof() && rows.size() == 3 &&
                rows[0].pair.q == 96 && rows[0].pair.r == 595 &&
                rows[0].bodies.size() == 116 &&
                rows[1].pair.q == 100 && rows[1].pair.r == 595 &&
                rows[1].bodies.size() == 13 &&
                rows[2].pair.q == 210 && rows[2].pair.r == 595 &&
                rows[2].bodies.size() == 16 &&
                distinct.size() == 145 && pair_ledger.state == kFailurePairFnv,
            "failure universe identity changed");
    return rows;
}

Response response(u32 mask, const std::array<bool, 3>& active,
                  const std::vector<FailureRow>& rows) {
    Response answer;
    if (active[0]) {
        for (std::size_t index = 0; index < rows[0].bodies.size(); ++index)
            if ((mask & rows[0].bodies[index]) == 0)
                answer.word[index / 64] |= UINT64_C(1) << (index % 64);
    }
    if (active[1]) {
        for (std::size_t index = 0; index < rows[1].bodies.size(); ++index)
            if ((mask & rows[1].bodies[index]) == 0)
                answer.word[2] |= UINT64_C(1) << index;
    }
    if (active[2]) {
        for (std::size_t index = 0; index < rows[2].bodies.size(); ++index)
            if ((mask & rows[2].bodies[index]) == 0)
                answer.word[3] |= UINT64_C(1) << index;
    }
    return answer;
}

std::string word16(u64 value) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << value;
    return out.str();
}

u64 sorted_mask_fnv(std::vector<u32> masks) {
    std::sort(masks.begin(), masks.end());
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::set<std::pair<int, int>> read_pair_set595(
    const std::filesystem::path& path, std::size_t expected, u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open old typed union");
    std::set<std::pair<int, int>> answer;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos, "malformed typed row");
        const int q = std::stoi(line.substr(0, comma));
        const int r = std::stoi(line.substr(comma + 1));
        require(q > 0 && q < r && answer.insert({q, r}).second,
                "typed row invalid/duplicate");
        ledger.add(q);
        ledger.add(r);
    }
    require(input.eof() && answer.size() == expected &&
                ledger.state == expected_fnv,
            "old typed union identity changed");
    return answer;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 13,
                "usage: response_atlas JOINT BASE8951 ADD45 SUFFIX9 RESIDUAL "
                "OLD_UNION1624 REPAIRS76 ADDITIONS4 DELETE73 FAILURES145 "
                "ATLAS_OUT INACTIVE_OUT");

        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> repairs =
            read_mixed617(argv[7], 76, kRepairFnv);
        const std::vector<u32> additions =
            read_mixed617(argv[8], 4, kAdditionFnv);
        const std::vector<u32> deletions =
            read_mixed617(argv[9], 73, kDeleteFnv);
        const std::vector<FailureRow> failure_rows = read_failures595(argv[10]);

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
        const std::set<u32> deletion_set(deletions.begin(), deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : augmented)
            if (!deletion_set.contains(mask)) carrier.push_back(mask);
        for (u32 addition : additions) {
            require(distinct.insert(addition).second, "addition overlap");
            carrier.push_back(addition);
        }
        require(carrier.size() == 9019 &&
                    masks_fnv_agent(carrier) == kCarrierFnv,
                "C596 carrier changed");
        const std::unordered_set<u32> carrier_set(carrier.begin(), carrier.end());

        std::array<Geometry, 3> geometries = {
            build_geometry(96, 595), build_geometry(100, 595),
            build_geometry(210, 595)};
        std::unordered_map<Response, ResponseStats, ResponseHash> atlas;
        atlas.reserve(100000);
        std::array<std::array<u64, 3>, 2> active{};
        std::array<u64, 2> responders{};
        std::array<Fnv, 2> responder_ledgers{};
        for (unsigned rank : {8U, 9U}) {
            const unsigned slot = rank - 8;
            u64 universe = 0;
            const u32 limit = UINT32_C(1) << 30;
            for (u32 mask = (UINT32_C(1) << rank) - 1; mask < limit;
                 mask = next_combination(mask)) {
                ++universe;
                std::array<bool, 3> is_active{};
                for (unsigned row = 0; row < 3; ++row) {
                    is_active[row] = margin(geometries[row], mask).ticks >= 0;
                    active[slot][row] += is_active[row];
                }
                if (!is_active[0] && !is_active[1] && !is_active[2]) continue;
                const Response reply = response(mask, is_active, failure_rows);
                if (empty(reply)) continue;
                require(!carrier_set.contains(mask),
                        "existing carrier mask responds to a frozen failure");
                ++responders[slot];
                responder_ledgers[slot].add(mask);
                ResponseStats& stats = atlas[reply];
                if (rank == 8) {
                    ++stats.count8;
                    stats.least8 = std::min(stats.least8, mask);
                } else {
                    ++stats.count9;
                    stats.least9 = std::min(stats.least9, mask);
                }
            }
            require(universe == (rank == 8 ? UINT64_C(5852925)
                                           : UINT64_C(14307150)),
                    "mask universe changed");
        }
        require(active[0] == std::array<u64, 3>{928827, 1116371, 1328767} &&
                    active[1] == std::array<u64, 3>{6076461, 6735949, 7133375} &&
                    responders[0] == 51271 && responders[1] == 679004 &&
                    responder_ledgers[0].state == kResponder8Fnv &&
                    responder_ledgers[1].state == kResponder9Fnv,
                "response census changed");

        std::vector<std::pair<Response, ResponseStats>> ordered(
            atlas.begin(), atlas.end());
        std::sort(ordered.begin(), ordered.end(), [](const auto& left,
                                                     const auto& right) {
            const unsigned wl = weight(left.first);
            const unsigned wr = weight(right.first);
            if (wl != wr) return wl > wr;
            return left.first.word > right.first.word;
        });
        std::ofstream atlas_out(argv[11]);
        require(static_cast<bool>(atlas_out), "cannot create response atlas");
        atlas_out << "w0,w1,w2,w3,count8,least8,count9,least9,weight\n";
        Fnv type_ledger;
        u64 rank8_types = 0;
        u64 rank9_types = 0;
        for (const auto& [reply, stats] : ordered) {
            for (u64 word : reply.word) type_ledger.add(word);
            type_ledger.add(stats.count8);
            type_ledger.add(stats.count9);
            rank8_types += stats.count8 != 0;
            rank9_types += stats.count9 != 0;
            atlas_out << word16(reply.word[0]) << ',' << word16(reply.word[1])
                      << ',' << word16(reply.word[2]) << ','
                      << word16(reply.word[3]) << ',' << stats.count8 << ',';
            if (stats.count8) atlas_out << hex8(stats.least8);
            atlas_out << ',' << stats.count9 << ',';
            if (stats.count9) atlas_out << hex8(stats.least9);
            atlas_out << ',' << weight(reply) << '\n';
        }
        require(atlas_out.good(), "response atlas write failed");
        require(ordered.size() == 14619 && rank8_types == 2210 &&
                    rank9_types == 14619 && type_ledger.state == kTypeLedgerFnv,
                "response type quotient changed");

        const std::vector<PairAgent> raw_band = read_band_agent(argv[5], 595);
        require(raw_band.size() == 394, "raw endpoint-at-least-595 band changed");
        const std::set<std::pair<int, int>> old_union =
            read_pair_set595(argv[6], 1624, kOldUnionFnv);
        std::vector<PairAgent> band;
        for (PairAgent pair : raw_band)
            if (pair.r >= 596 || !old_union.contains({pair.q, pair.r}))
                band.push_back(pair);
        require(band.size() == 391, "carrier-preservation prefix changed");
        std::vector<Geometry> prefix_geometries;
        prefix_geometries.reserve(band.size());
        for (PairAgent pair : band)
            prefix_geometries.push_back(build_geometry(pair.q, pair.r));

        struct InactiveAudit {
            std::vector<u32> masks;
            u64 equalities = 0;
            Fnv signs;
        };
        auto inactive_audit = [&](const std::vector<u32>& family) {
            InactiveAudit audit;
            for (u32 mask : family) {
                bool ever_active = false;
                for (const Geometry& geometry : prefix_geometries) {
                    const Margin exact = margin(geometry, mask);
                    const bool here = exact.ticks >= 0;
                    ever_active |= here;
                    audit.equalities += exact.ticks == 0;
                    audit.signs.add(mask);
                    audit.signs.add(here);
                }
                if (!ever_active) audit.masks.push_back(mask);
            }
            std::sort(audit.masks.begin(), audit.masks.end());
            return audit;
        };
        const InactiveAudit augmented_inactive = inactive_audit(augmented);
        const InactiveAudit carrier_inactive = inactive_audit(carrier);
        require(augmented_inactive.equalities == 0 &&
                    augmented_inactive.masks.size() == 75 &&
                    sorted_mask_fnv(augmented_inactive.masks) ==
                        UINT64_C(0xfa143e58f59119f8) &&
                    augmented_inactive.signs.state == kAugmentedPrefixSignFnv,
                "augmented prefix inactivity changed");
        require(carrier_inactive.equalities == 0 &&
                    carrier_inactive.masks ==
                        std::vector<u32>{UINT32_C(0x380086a0),
                                         UINT32_C(0x388088c0)} &&
                    sorted_mask_fnv(carrier_inactive.masks) ==
                        kCarrierInactiveFnv &&
                    carrier_inactive.signs.state == kCarrierPrefixSignFnv,
                "C596 prefix inactivity changed");
        std::ofstream inactive_out(argv[12]);
        require(static_cast<bool>(inactive_out), "cannot create inactive list");
        inactive_out << "family,mask_hex,rank,joint\n";
        for (u32 mask : augmented_inactive.masks)
            inactive_out << "augmented9088," << hex8(mask) << ','
                         << std::popcount(mask) << ','
                         << joint_set.contains(mask) << '\n';
        for (u32 mask : carrier_inactive.masks)
            inactive_out << "carrier9019," << hex8(mask) << ','
                         << std::popcount(mask) << ','
                         << joint_set.contains(mask) << '\n';
        require(inactive_out.good(), "inactive list write failed");

        std::cout << "LRC14_ENDPOINT595_PAIR_TAGGED_RESPONSE_ATLAS_V1\n"
                  << "FAILURES 145 PAIRS 96,595:116 100,595:13 210,595:16 "
                     "PAIR_TAGGED_FNV "
                  << std::hex << kFailurePairFnv << std::dec << '\n'
                  << "CARRIER 9019 FNV " << std::hex << kCarrierFnv
                  << std::dec << " RANK8 8961 RANK9 58\n";
        for (unsigned slot = 0; slot < 2; ++slot) {
            const unsigned rank = slot + 8;
            std::cout << "RANK " << rank << " UNIVERSE "
                      << (rank == 8 ? UINT64_C(5852925)
                                    : UINT64_C(14307150))
                      << " ACTIVE_BY_PAIR " << active[slot][0] << ','
                      << active[slot][1] << ',' << active[slot][2]
                      << " RESPONDERS " << responders[slot] << " FNV "
                      << std::hex << responder_ledgers[slot].state << std::dec
                      << '\n';
        }
        std::cout << "RESPONSE_TYPES " << ordered.size() << " RANK8_TYPES "
                  << rank8_types << " RANK9_TYPES " << rank9_types
                  << " TYPE_LEDGER_FNV " << std::hex << type_ledger.state
                  << std::dec << '\n'
                  << "PREFIX_ROWS " << band.size()
                  << " ENDPOINT_LOWER 595\n"
                  << "AUGMENTED9088_SIGN_CELLS "
                  << augmented.size() * band.size() << " SIGN_FNV "
                  << std::hex << augmented_inactive.signs.state << std::dec
                  << " EQUALITIES " << augmented_inactive.equalities
                  << " COMMON_INACTIVE " << augmented_inactive.masks.size()
                  << " FNV " << std::hex
                  << sorted_mask_fnv(augmented_inactive.masks) << std::dec
                  << '\n'
                  << "CARRIER9019_SIGN_CELLS " << carrier.size() * band.size()
                  << " SIGN_FNV " << std::hex
                  << carrier_inactive.signs.state << std::dec << " EQUALITIES "
                  << carrier_inactive.equalities << " COMMON_INACTIVE "
                  << carrier_inactive.masks.size() << " FNV " << std::hex
                  << sorted_mask_fnv(carrier_inactive.masks) << std::dec << '\n'
                  << "OUTPUT_ATLAS " << argv[11] << " OUTPUT_INACTIVE "
                  << argv[12] << '\n'
                  << "SCOPE FINITE_EXACT_FROZEN_POOL_PAIR_TAGGED_RESPONSES_"
                     "AND_STRICT_INACTIVITY_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT595_RESPONSE_ATLAS_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
