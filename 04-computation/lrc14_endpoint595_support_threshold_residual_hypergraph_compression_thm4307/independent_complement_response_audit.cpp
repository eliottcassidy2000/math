// Independent complement-generated reconstruction of the complete endpoint-
// 595 pair-tagged response quotient.  Unlike response_atlas.cpp, this does not
// scan either full rank universe: it generates responder candidates as
// subsets of each failure body's 21-label complement.

#define ENDPOINT617_RAW_VERIFY_MAIN endpoint595_complement_hidden_main
#include "04-computation/lrc14_size_preserving_response_staircase_thm4300/endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

#include <fstream>
#include <unordered_map>

namespace {
constexpr u64 kFailurePairFnv = UINT64_C(0xf3d7f95fc38e7b49);
constexpr u64 kResponder8Fnv = UINT64_C(0xcc926b13c719225d);
constexpr u64 kResponder9Fnv = UINT64_C(0x8ff584ab870b72a1);
constexpr u64 kTypeLedgerFnv = UINT64_C(0xa5f5ebcdeb03ad34);

struct ResponseI {
    std::array<u64, 4> word{};
    bool operator==(const ResponseI&) const = default;
};
struct ResponseIHash {
    std::size_t operator()(const ResponseI& value) const noexcept {
        u64 h = UINT64_C(0x9e3779b97f4a7c15);
        for (u64 x : value.word) {
            x ^= x >> 30; x *= UINT64_C(0xbf58476d1ce4e5b9);
            x ^= x >> 27; x *= UINT64_C(0x94d049bb133111eb);
            x ^= x >> 31;
            h ^= x + UINT64_C(0x9e3779b97f4a7c15) + (h << 6) + (h >> 2);
        }
        return static_cast<std::size_t>(h);
    }
};
struct StatsI {
    u64 count8 = 0, count9 = 0;
    u32 least8 = UINT32_MAX, least9 = UINT32_MAX;
};
struct FailureRowI { PairAgent pair; std::vector<u32> bodies; };

std::vector<FailureRowI> read_failures_independent(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failures");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    std::vector<FailureRowI> rows;
    Fnv ledger;
    std::set<std::tuple<int, int, u32>> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q, r; std::string token;
        fields >> q >> r >> token;
        const u32 body = parse_mask_agent(token);
        require(fields && std::popcount(body) == 9 &&
                    distinct.insert({q, r, body}).second,
                "failure row invalid/duplicate");
        if (rows.empty() || rows.back().pair.q != q || rows.back().pair.r != r)
            rows.push_back(FailureRowI{PairAgent{q, r}, {}});
        rows.back().bodies.push_back(body);
        ledger.add(q); ledger.add(r); ledger.add(body);
    }
    require(input.eof() && rows.size() == 3 && rows[0].bodies.size() == 116 &&
                rows[1].bodies.size() == 13 && rows[2].bodies.size() == 16 &&
                ledger.state == kFailurePairFnv,
            "failure universe changed");
    return rows;
}

template <class Consumer>
void enumerate_complement_i(u32 body, unsigned rank, Consumer consume) {
    std::array<unsigned, 21> free{};
    unsigned count = 0;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((body & (UINT32_C(1) << bit)) == 0) free[count++] = bit;
    require(count == 21, "body complement changed");
    const u32 limit = UINT32_C(1) << 21;
    for (u32 local = (UINT32_C(1) << rank) - 1; local < limit;
         local = next_combination(local)) {
        u32 mask = 0;
        for (unsigned index = 0; index < 21; ++index)
            if (local & (UINT32_C(1) << index))
                mask |= UINT32_C(1) << free[index];
        consume(mask);
    }
}

ResponseI direct_response_i(u32 mask,
                            const std::array<Geometry, 3>& geometries,
                            const std::vector<FailureRowI>& rows) {
    ResponseI answer;
    for (unsigned row = 0; row < 3; ++row) {
        if (margin(geometries[row], mask).ticks < 0) continue;
        for (std::size_t index = 0; index < rows[row].bodies.size(); ++index) {
            if (mask & rows[row].bodies[index]) continue;
            if (row == 0) answer.word[index / 64] |= UINT64_C(1) << (index % 64);
            else answer.word[row + 1] |= UINT64_C(1) << index;
        }
    }
    return answer;
}

unsigned weight_i(const ResponseI& reply) {
    unsigned answer = 0;
    for (u64 word : reply.word) answer += std::popcount(word);
    return answer;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: independent_complement FAILURES145");
        const auto rows = read_failures_independent(argv[1]);
        const std::array<Geometry, 3> geometries = {
            build_geometry(96, 595), build_geometry(100, 595),
            build_geometry(210, 595)};
        std::array<std::unordered_map<u32, ResponseI>, 2> responders;
        responders[0].reserve(60000);
        responders[1].reserve(700000);
        std::array<u64, 2> generated_events{};
        for (unsigned slot = 0; slot < 2; ++slot) {
            const unsigned rank = slot + 8;
            for (unsigned row = 0; row < 3; ++row) {
                for (std::size_t obligation = 0;
                     obligation < rows[row].bodies.size(); ++obligation) {
                    enumerate_complement_i(rows[row].bodies[obligation], rank,
                                           [&](u32 mask) {
                        ++generated_events[slot];
                        if (margin(geometries[row], mask).ticks < 0) return;
                        ResponseI& reply = responders[slot][mask];
                        if (row == 0)
                            reply.word[obligation / 64] |=
                                UINT64_C(1) << (obligation % 64);
                        else
                            reply.word[row + 1] |= UINT64_C(1) << obligation;
                    });
                }
            }
        }
        require(generated_events[0] == UINT64_C(29506050) &&
                    generated_events[1] == UINT64_C(42619850),
                "complement event census changed");

        std::unordered_map<ResponseI, StatsI, ResponseIHash> atlas;
        atlas.reserve(16000);
        std::array<u64, 2> responder_fnv{};
        for (unsigned slot = 0; slot < 2; ++slot) {
            std::vector<u32> masks;
            masks.reserve(responders[slot].size());
            for (const auto& [mask, reply] : responders[slot]) {
                require(reply == direct_response_i(mask, geometries, rows),
                        "accumulated response differs from direct response");
                masks.push_back(mask);
            }
            std::sort(masks.begin(), masks.end());
            Fnv ledger;
            for (u32 mask : masks) {
                ledger.add(mask);
                StatsI& stats = atlas[responders[slot].at(mask)];
                if (slot == 0) {
                    ++stats.count8; stats.least8 = std::min(stats.least8, mask);
                } else {
                    ++stats.count9; stats.least9 = std::min(stats.least9, mask);
                }
            }
            responder_fnv[slot] = ledger.state;
        }
        require(responders[0].size() == 51271 && responders[1].size() == 679004 &&
                    responder_fnv[0] == kResponder8Fnv &&
                    responder_fnv[1] == kResponder9Fnv && atlas.size() == 14619,
                "independent responder quotient changed");

        std::vector<std::pair<ResponseI, StatsI>> ordered(atlas.begin(), atlas.end());
        std::sort(ordered.begin(), ordered.end(), [](const auto& left,
                                                     const auto& right) {
            if (weight_i(left.first) != weight_i(right.first))
                return weight_i(left.first) > weight_i(right.first);
            return left.first.word > right.first.word;
        });
        Fnv type_ledger;
        u64 rank8_types = 0;
        for (const auto& [reply, stats] : ordered) {
            for (u64 word : reply.word) type_ledger.add(word);
            type_ledger.add(stats.count8); type_ledger.add(stats.count9);
            rank8_types += stats.count8 != 0;
        }
        require(rank8_types == 2210 && type_ledger.state == kTypeLedgerFnv,
                "independent type ledger changed");
        std::cout << "LRC14_ENDPOINT595_INDEPENDENT_COMPLEMENT_RESPONSE_V1\n"
                  << "FAILURES 145 PAIR_TAGGED_FNV " << std::hex
                  << kFailurePairFnv << std::dec << '\n'
                  << "ENUMERATION complement_subsets_of_each_tagged_body\n"
                  << "GENERATED_EVENTS_RANK8 " << generated_events[0]
                  << " RANK9 " << generated_events[1] << '\n'
                  << "RESPONDERS_RANK8 " << responders[0].size() << " FNV "
                  << std::hex << responder_fnv[0] << std::dec
                  << " RESPONDERS_RANK9 " << responders[1].size() << " FNV "
                  << std::hex << responder_fnv[1] << std::dec << '\n'
                  << "RESPONSE_TYPES " << atlas.size() << " RANK8_TYPES "
                  << rank8_types << " TYPE_LEDGER_FNV " << std::hex
                  << type_ledger.state << std::dec << '\n'
                  << "DIRECT_RESPONSE_CROSSCHECK_MASKS "
                  << responders[0].size() + responders[1].size() << '\n'
                  << "SCOPE INDEPENDENT_COMPLETE_PAIR_TAGGED_RESPONSE_QUOTIENT_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "INDEPENDENT_COMPLEMENT_RESPONSE_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
