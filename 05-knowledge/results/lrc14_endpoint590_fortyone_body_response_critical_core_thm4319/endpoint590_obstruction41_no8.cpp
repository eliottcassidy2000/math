// Exact no-eight search on a 41-obligation critical endpoint-590 subinstance.
// This is an exploratory sidecar over the already-frozen complete response
// atlas, not a physical lonely-runner statement.

#define main endpoint590_full_search_hidden_main
#include "endpoint590_exact_no8_search.cpp"
#undef main

int main(int argc, char** argv) {
    try {
        if (argc != 3)
            throw std::runtime_error("usage: obstruction41 SIGNATURES DUAL");
        const std::array<unsigned, 41> ordinals{
            0, 5, 9, 12, 14, 18, 20, 22, 24, 25, 29, 32, 35, 37,
            38, 40, 43, 47, 48, 53, 55, 57, 63, 65, 68, 71, 73, 76,
            77, 79, 80, 83, 84, 88, 89, 90, 94, 95, 96, 97, 99,
        };
        Bits universe{};
        for (unsigned ordinal : ordinals) {
            if (ordinal < 64) universe.lo |= u64{1} << ordinal;
            else universe.hi |= u64{1} << (ordinal - 64);
        }

        auto traces = read_signatures(argv[1]);
        for (Bits& trace : traces) trace = intersect(trace, universe);
        traces.erase(std::remove(traces.begin(), traces.end(), Bits{}),
                     traces.end());
        std::sort(traces.begin(), traces.end(), [](Bits left, Bits right) {
            if (weight(left) != weight(right)) return weight(left) > weight(right);
            if (left.hi != right.hi) return left.hi < right.hi;
            return left.lo < right.lo;
        });
        traces.erase(std::unique(traces.begin(), traces.end()), traces.end());
        if (traces.size() != 2041)
            throw std::runtime_error("restricted-signature count changed");

        Search search;
        search.dual = read_dual(argv[2]);
        for (Bits candidate : traces) {
            bool dominated = false;
            for (Bits kept : search.sets) {
                if (subset(candidate, kept)) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) search.sets.push_back(candidate);
        }
        if (search.sets.size() != 395)
            throw std::runtime_error("maximal-trace quotient changed");
        for (unsigned index = 0; index < search.sets.size(); ++index) {
            if (search.dual_weight(search.sets[index]) > 3)
                throw std::runtime_error("dual load exceeded three");
            for (unsigned bit = 0; bit < 100; ++bit) {
                if (contains(search.sets[index], bit))
                    search.coverers[bit].push_back(index);
            }
        }
        const bool satisfiable = search.dfs(universe, 8);
        u64 deletion_nodes = 0;
        u64 deletion_dead = 0;
        u64 deletion_sum_prunes = 0;
        u64 deletion_dual_prunes = 0;
        unsigned deletion_sat = 0;
        for (unsigned omitted : ordinals) {
            Search reduced;
            reduced.sets = search.sets;
            reduced.coverers = search.coverers;
            reduced.dual = search.dual;
            Bits reduced_universe = universe;
            if (omitted < 64) reduced_universe.lo &= ~(u64{1} << omitted);
            else reduced_universe.hi &= ~(u64{1} << (omitted - 64));
            if (reduced.dfs(reduced_universe, 7)) ++deletion_sat;
            deletion_nodes += reduced.nodes;
            deletion_dead += reduced.dead.size();
            deletion_sum_prunes += reduced.sum_bound_prunes;
            deletion_dual_prunes += reduced.dual_bound_prunes;
        }
        if (satisfiable || search.nodes != UINT64_C(1968373) ||
            search.dead.size() != 1698303 ||
            search.sum_bound_prunes != UINT64_C(530603) ||
            search.dual_bound_prunes != UINT64_C(1009466) ||
            deletion_sat != 0 || deletion_nodes != UINT64_C(278730) ||
            deletion_dead != UINT64_C(270531) ||
            deletion_sum_prunes != UINT64_C(52383) ||
            deletion_dual_prunes != UINT64_C(198392))
            throw std::runtime_error("critical-core search identity changed");
        std::cout << "ENDPOINT590_OBSTRUCTION41_NO8_V1\n"
                  << "UNIVERSE 41 RESTRICTED_SIGNATURES " << traces.size()
                  << " MAXIMAL_TRACES " << search.sets.size() << '\n'
                  << "DUAL_WEIGHT 22 CAPACITY_PER_RESPONSE 3 DEPTH 8\n"
                  << "NODES " << search.nodes << " DEAD " << search.dead.size()
                  << " SUM_PRUNES " << search.sum_bound_prunes
                  << " DUAL_PRUNES " << search.dual_bound_prunes << '\n'
                  << "RESULT " << (satisfiable ? "SAT" : "UNSAT") << '\n'
                  << "SINGLE_DELETION_DEPTH7 SAT " << deletion_sat
                  << " UNSAT " << (ordinals.size() - deletion_sat)
                  << " NODES " << deletion_nodes << " DEAD " << deletion_dead
                  << " SUM_PRUNES " << deletion_sum_prunes
                  << " DUAL_PRUNES " << deletion_dual_prunes << '\n'
                  << "SCOPE EXACT_COMPLETE_RESPONSE_ATLAS_RESTRICTED_TO_"
                     "FIXED_41_OBLIGATION_SUBINSTANCE_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT " << (satisfiable ? "FAIL" : "PASS") << '\n';
        return (satisfiable || deletion_sat != 0) ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "OBSTRUCTION41_ERROR " << error.what() << '\n';
        return 1;
    }
}
