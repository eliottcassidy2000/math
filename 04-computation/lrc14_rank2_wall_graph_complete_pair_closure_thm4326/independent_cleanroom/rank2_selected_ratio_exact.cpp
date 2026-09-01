// Selected exact replay used only to close the normalized-margin filter.
// Reuses this packet's clean-room i128 event graph and exact optimizer, not any
// source or result from the primary theorem packet.

#define main rank2_wide_full_audit_main_not_used_here
#include "rank2_event_sweep_wide_cleanroom.cpp"
#undef main

int main() {
    try {
        const std::vector<BaseEvent> base_events = make_base_events();
        constexpr i128 candidate_ticks = INT64_C(245428469244);
        constexpr i128 candidate_grid = INT64_C(91205797082400);
        const std::array<Pair, 4> selected = {
            Pair{50, 70}, Pair{50, 212}, Pair{50, 274}, Pair{100, 110}};

        std::cout << "LRC14_RANK2_SELECTED_NORMALIZED_EXACT_AUDIT_V1\n";
        unsigned challengers_above = 0;
        for (const Pair pair : selected) {
            const Graph graph = event_sweep_graph(pair.q, pair.r, base_events);
            const Coarse coarse = coarse_bound(graph);
            ExactOptimizer optimizer(graph);
            optimizer.solve();
            const i128 minimum_mass = graph.total - optimizer.best_reward();
            const i128 minimum_ticks = 63 * minimum_mass - 4 * graph.grid;
            const i128 direct = survivor_direct(graph, optimizer.best_body());
            demand(direct == minimum_mass, "selected direct replay mismatch");
            const i128 ratio_cross =
                minimum_ticks * candidate_grid - candidate_ticks * graph.grid;
            const bool is_candidate = pair.q == 50 && pair.r == 70;
            if (!is_candidate) challengers_above += ratio_cross > 0;
            if (is_candidate) demand(ratio_cross == 0, "candidate ratio drift");

            std::cout << "PAIR " << pair.q << ',' << pair.r
                      << " GRID " << dec(graph.grid)
                      << " COARSE_TICKS " << dec(coarse.ticks)
                      << " MIN_MASS " << dec(minimum_mass)
                      << " MIN_TICKS " << dec(minimum_ticks)
                      << " MIN_BODY " << hex8(optimizer.best_body())
                      << " DIRECT_MASS " << dec(direct)
                      << " RATIO_CROSS_VS_50_70 " << dec(ratio_cross)
                      << " NODES " << optimizer.nodes()
                      << " PRUNES " << optimizer.prunes() << "\n";
        }
        std::cout << "NORMALIZED_CHALLENGERS 3 STRICTLY_ABOVE_CANDIDATE "
                  << challengers_above << "\n";
        std::cout << "SCOPE FINITE_EXACT_SELECTED_FILTER_ROWS_NO_PHYSICAL_ENTRY_NO_LRC14\n";
        std::cout << "VERDICT " << (challengers_above == 3 ? "PASS" : "FAIL")
                  << "\n";
        return challengers_above == 3 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
