// Detached exact literal-wall audit of the endpoint-632 unresponsive body.
#define SINGLETON_LITERAL_LIBRARY_ONLY
#include "04-computation/lrc14_signature_response_congruence_thm4286/singleton_fibre_literal_verify.cpp"
#undef SINGLETON_LITERAL_LIBRARY_ONLY

int main() {
    try {
        std::cout.setf(std::ios::unitbuf);
        constexpr int q = 256;
        constexpr int r = 632;
        constexpr u32 body = UINT32_C(0x1d106401);
        require(std::popcount(body) == 9, "body rank changed");
        const auto cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        std::vector<i64> walls(cells.size() + 1);
        walls.front() = cells.front().left;
        for (std::size_t i = 0; i < cells.size(); ++i) {
            require(cells[i].left == walls[i], "pool wall discontinuity");
            walls[i + 1] = cells[i].right;
        }
        require(walls.front() == 0 && walls.back() == COMMON,
                "pool wall boundary changed");
        const LiteralPair pair = build_literal_pair(q, r);
        std::vector<i128> prefix(walls.size());
        for (std::size_t i = 0; i < walls.size(); ++i)
            prefix[i] = prefix_at(
                pair, static_cast<i128>(walls[i]) * pair.pool_scale);

        u64 universe = 0;
        u64 disjoint = 0;
        u64 positive = 0;
        u64 equal = 0;
        i128 maximum = 0;
        u32 maximum_mask = 0;
        bool have_maximum = false;
        FnvLocal ledger;
        u32 candidate = (u32{1} << 8) - 1;
        while (candidate < (u32{1} << 30)) {
            ++universe;
            if ((candidate & body) == 0) {
                ++disjoint;
                const IndexedRepair repair = index_repair(candidate, cells);
                const i128 margin = repair_margin(repair, prefix, pair.grid);
                positive += margin > 0;
                equal += margin == 0;
                if (!have_maximum || margin > maximum ||
                    (margin == maximum && candidate < maximum_mask)) {
                    have_maximum = true;
                    maximum = margin;
                    maximum_mask = candidate;
                }
                ledger.add(candidate);
                add_i128(ledger, margin);
                add_i128(ledger, pair.grid);
            }
            const u32 next = next_combination(candidate);
            if (next <= candidate) break;
            candidate = next;
        }
        require(universe == UINT64_C(5852925) &&
                    disjoint == UINT64_C(203490) && have_maximum,
                "repair universe changed");
        const i128 divisor = gcd_i128(maximum < 0 ? -maximum : maximum,
                                      pair.grid);
        std::cout << "THM4295_ENDPOINT632_UNRESPONSIVE_LITERAL_AUDIT_V1\n"
                  << "PAIR " << q << ',' << r << " BODY " << std::hex
                  << std::setw(8) << std::setfill('0') << body << std::dec
                  << std::setfill(' ') << " POOL_CELLS " << cells.size()
                  << " GRID " << decimal(pair.grid) << '\n'
                  << "RANK8_UNIVERSE " << universe << " DISJOINT " << disjoint
                  << " POSITIVE " << positive << " EQUAL " << equal
                  << " NEGATIVE " << disjoint - positive - equal << '\n'
                  << "MAXIMUM_MASK " << std::hex << std::setw(8)
                  << std::setfill('0') << maximum_mask << std::dec
                  << std::setfill(' ') << " MAXIMUM_MARGIN "
                  << decimal(maximum) << " DEN " << decimal(pair.grid)
                  << " REDUCED " << decimal(maximum / divisor) << '/'
                  << decimal(pair.grid / divisor) << '\n'
                  << "MARGIN_LEDGER_FNV " << std::hex << ledger.state
                  << std::dec << '\n'
                  << "VERDICT PASS DETACHED_ALL_DISJOINT_MASKS_INACTIVE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "UNRESPONSIVE_LITERAL_ERROR " << error.what() << '\n';
        return 1;
    }
}
