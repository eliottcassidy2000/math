#include "independent_common.hpp"

namespace {
using namespace ray_audit;

constexpr std::array<i64, 4> SCALES = {146, 256, 1015, 1016};

struct Choose {
    std::array<std::array<u64, 9>, 31> c{};
    Choose() {
        for (int n = 0; n <= 30; ++n) {
            c[n][0] = 1;
            for (int k = 1; k <= 8; ++k) {
                c[n][k] = n == 0 ? 0 : c[n - 1][k] + c[n - 1][k - 1];
            }
        }
        ensure(c[30][8] == REPAIR_COUNT, "binomial table mismatch");
    }
    u64 rank(u32 mask) const {
        u64 answer = 0;
        int ordinal = 1;
        for (int bit = 0; bit < 30; ++bit) {
            if ((mask & (u32{1} << bit)) == 0) continue;
            answer += c[bit][ordinal++];
        }
        ensure(ordinal == 9 && answer < REPAIR_COUNT, "rank-eight failure");
        return answer;
    }
};

template <class Callback>
void complete(u32 fixed, int missing, int next_bit, u32 extra,
              const Callback& callback) {
    if (missing == 0) {
        callback(fixed | extra);
        return;
    }
    for (int bit = next_bit; bit < 30; ++bit) {
        const u32 flag = u32{1} << bit;
        if ((fixed & flag) != 0) continue;
        complete(fixed, missing - 1, bit + 1, extra | flag, callback);
    }
}

using FourMasses = std::array<i128, SCALES.size()>;

}  // namespace

int main() {
    try {
        using namespace ray_audit;
        check_primitive_constants();
        const std::vector<Cell> arrangement = build_pool_arrangement();
        const Choose choose;
        std::map<u32, FourMasses> atom;
        for (const Cell& cell : arrangement) {
            if (std::popcount(cell.failed) > 8) continue;
            FourMasses& target = atom[cell.failed];
            for (std::size_t lane = 0; lane < SCALES.size(); ++lane) {
                const i64 g = SCALES[lane];
                target[lane] +=
                    primitive_integral_numerator(static_cast<i128>(g) * cell.right)
                  - primitive_integral_numerator(static_cast<i128>(g) * cell.left);
            }
        }

        std::array<std::vector<i128>, SCALES.size()> mass;
        for (auto& lane : mass) lane.assign(REPAIR_COUNT, 0);
        u64 transform_ops = 0;
        for (const auto& [fixed, values] : atom) {
            complete(fixed, 8 - std::popcount(fixed), 0, 0, [&](u32 repair) {
                const u64 rank = choose.rank(repair);
                for (std::size_t lane = 0; lane < SCALES.size(); ++lane) {
                    mass[lane][rank] += values[lane];
                }
                ++transform_ops;
            });
        }
        ensure(transform_ops == UINT64_C(152170690),
               "complete-deck transform operation mismatch");

        std::array<std::vector<u32>, SCALES.size()> active;
        u32 repair = (u32{1} << 8) - 1;
        const u32 limit = u32{1} << 30;
        u64 rank = 0;
        while (repair != 0 && repair < limit) {
            ensure(choose.rank(repair) == rank, "repair rank traversal mismatch");
            for (std::size_t lane = 0; lane < SCALES.size(); ++lane) {
                const i64 g = SCALES[lane];
                const i128 margin = static_cast<i128>(63) * mass[lane][rank]
                                  - static_cast<i128>(4) * N * D * g;
                if (margin >= 0) active[lane].push_back(repair);
            }
            ++rank;
            const u32 following = next_mask(repair);
            if (following <= repair) break;
            repair = following;
        }
        ensure(rank == REPAIR_COUNT, "complete repair universe truncated");
        for (auto& deck : active) {
            std::sort(deck.begin(), deck.end(), [](u32 a, u32 b) {
                const u64 ka = splitmix64(static_cast<u64>(a) ^ ORDER_SEED);
                const u64 kb = splitmix64(static_cast<u64>(b) ^ ORDER_SEED);
                return ka != kb ? ka < kb : a < b;
            });
        }
        for (auto& lane : mass) {
            lane.clear();
            lane.shrink_to_fit();
        }

        const std::vector<u32> bodies = all_bodies();
        Fnv64 transcript;
        std::cout << "LRC_23_RAY_INDEPENDENT_COMPLETE_DECK_CONTROLS_V1\n";
        std::cout << "TRANSFORM_OPS " << transform_ops << " ATOMS "
                  << atom.size() << '\n';
        for (std::size_t lane = 0; lane < SCALES.size(); ++lane) {
            Fnv64 deck_ledger;
            for (u32 mask : active[lane]) deck_ledger.word(mask);
            const BodyScan scan = scan_bodies(active[lane], bodies);
            ensure(scan.failures == 0, "complete exact deck leaves a body uncovered");
            Fnv64 prefix;
            for (u64 i = 0; i < scan.max_checks; ++i) prefix.word(active[lane][i]);
            transcript.word(static_cast<u64>(SCALES[lane]));
            transcript.word(active[lane].size());
            transcript.word(deck_ledger.value());
            transcript.word(scan.failures);
            transcript.word(scan.checks);
            transcript.word(scan.max_checks);
            transcript.word(scan.worst_body);
            transcript.word(prefix.value());
            std::cout << "G " << SCALES[lane] << " ACTIVE "
                      << active[lane].size() << " ACTIVE_FNV "
                      << hex64(deck_ledger.value()) << " BODIES " << scan.bodies
                      << " FAILURES " << scan.failures << " CHECKS "
                      << scan.checks << " MAX_CHECKS " << scan.max_checks
                      << " WORST_BODY {" << labels(scan.worst_body)
                      << "} PREFIX_FNV " << hex64(prefix.value()) << '\n';
        }
        ensure(active[1].size() == 3879889,
               "g=256 complete active census mismatch");
        std::cout << "CONTROL_TRANSCRIPT_FNV " << hex64(transcript.value()) << '\n';
        std::cout << "VERDICT ALL_FOUR_COMPLETE_EXACT_DECK_CONTROLS_CLOSE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
