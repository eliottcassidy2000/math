#define ENDPOINT592_RESPONSE_POOL_MAIN endpoint592_pool_hidden_for_pricing
#include "endpoint592_retained_pool_export.cpp"
#undef ENDPOINT592_RESPONSE_POOL_MAIN

#include <unordered_set>

int main(int argc, char** argv) {
    try {
        require(argc == 5,
                "usage: pricing FAILURES2468 INTEGER_WEIGHTS2468 "
                "POOL_BITS_CSV VIOLATIONS_OUT");
        constexpr u64 kDualScale = UINT64_C(1000000000);
        const FailureUniverse592 failures = read_failures592(argv[1]);
        std::array<Geometry, 7> geometry;
        for (std::size_t row = 0; row < 7; ++row)
            geometry[row] = build_geometry(kQ592[row], 592);

        std::ifstream dual_in(argv[2]);
        require(static_cast<bool>(dual_in), "cannot open dual");
        std::vector<u64> dual;
        u64 value = 0;
        while (dual_in >> value) dual.push_back(value);
        require(dual_in.eof() && dual.size() == kTotalFailures592,
                "dual length changed");
        u64 dual_sum = 0;
        for (u64 item : dual) {
            dual_sum += item;
        }
        require(dual_sum == UINT64_C(35261737295),
                "integer dual numerator changed");
        std::array<std::array<u64, 64>, kTotalWords592> word_dual{};
        std::size_t ordinal = 0;
        for (std::size_t row = 0; row < 7; ++row) {
            const std::size_t base = failures.word_offset[row];
            for (std::size_t index = 0; index < kCounts592[row]; ++index)
                word_dual[base + index / 64][index % 64] = dual[ordinal++];
        }
        require(ordinal == dual.size(), "dual partition changed");

        std::ifstream pool_in(argv[3]);
        require(static_cast<bool>(pool_in), "cannot open pool csv");
        std::string line;
        require(std::getline(pool_in, line) && line.starts_with("mask_hex,"),
                "pool header changed");
        std::unordered_set<u32> pool;
        while (std::getline(pool_in, line)) {
            if (line.empty()) continue;
            const auto comma = line.find(',');
            require(comma == 8, "pool mask field changed");
            pool.insert(parse_mask_agent(line.substr(0, comma)));
        }
        require(pool_in.eof() && pool.size() == 37497, "pool size changed");

        std::ofstream violation_out(argv[4]);
        require(static_cast<bool>(violation_out),
                "cannot create violation ledger");
        violation_out << "mask_hex,rank,dual_numerator,in_retained_pool\n";
        std::array<u64, 2> responders{};
        std::array<u64, 2> omitted{};
        std::array<u64, 2> violating{};
        std::array<u64, 2> pool_violating{};
        u64 maximum_all = 0;
        u64 maximum_omitted = 0;
        u32 maximum_all_mask = 0, maximum_omitted_mask = 0;
        u32 first_violation_mask = 0;
        u64 first_violation_score = 0;
        Fnv violation_ledger;
        for (unsigned rank : {8U, 9U}) {
            const std::size_t slot = rank - 8;
            const u32 limit = UINT32_C(1) << 30;
            for (u32 mask = (UINT32_C(1) << rank) - 1; mask < limit;
                 mask = next_combination(mask)) {
                const Candidate592 candidate = response592(mask, geometry, failures);
                if (candidate.weight == 0) continue;
                ++responders[slot];
                u64 score = 0;
                for (std::size_t word_index = 0; word_index < kTotalWords592;
                     ++word_index) {
                    u64 word = candidate.response[word_index];
                    while (word) {
                        const unsigned bit = std::countr_zero(word);
                        score += word_dual[word_index][bit];
                        word &= word - 1;
                    }
                }
                if (score > maximum_all ||
                    (score == maximum_all && mask < maximum_all_mask)) {
                    maximum_all = score; maximum_all_mask = mask;
                }
                const bool in_pool = pool.contains(mask);
                if (!in_pool) {
                    ++omitted[slot];
                    if (score > maximum_omitted ||
                        (score == maximum_omitted && mask < maximum_omitted_mask)) {
                        maximum_omitted = score; maximum_omitted_mask = mask;
                    }
                }
                if (score > kDualScale) {
                    ++violating[slot];
                    if (in_pool) ++pool_violating[slot];
                    if (first_violation_mask == 0) {
                        first_violation_mask = mask;
                        first_violation_score = score;
                    }
                    violation_ledger.add(mask);
                    violation_ledger.add(score);
                    violation_out << hex8(mask) << ',' << rank << ',' << score
                                  << ',' << (in_pool ? 1 : 0) << '\n';
                }
            }
        }
        require(violation_out.good(), "violation ledger write failed");
        require(pool_violating[0] == 0 && pool_violating[1] == 0,
                "retained pool integer certificate contradicted");
        std::cout << "ENDPOINT592_FULL_UNIVERSE_INTEGER_DUAL_PRICING_V1\n"
                  << "DUAL_NUMERATOR " << dual_sum << " SCALE " << kDualScale
                  << " POSITIVE "
                  << std::count_if(dual.begin(), dual.end(), [](u64 x) {
                         return x > 0;
                     }) << "\n"
                  << "RANK8_RESPONDERS " << responders[0] << " OMITTED "
                  << omitted[0] << " VIOLATING " << violating[0] << "\n"
                  << "RANK9_RESPONDERS " << responders[1] << " OMITTED "
                  << omitted[1] << " VIOLATING " << violating[1] << "\n"
                  << "MAX_ALL_NUMERATOR " << maximum_all << " MASK "
                  << hex8(maximum_all_mask) << "\n"
                  << "MAX_OMITTED_NUMERATOR " << maximum_omitted << " MASK "
                  << hex8(maximum_omitted_mask) << "\n"
                  << "VIOLATION_FNV " << std::hex << violation_ledger.state
                  << std::dec << " FIRST_MASK " << hex8(first_violation_mask)
                  << " FIRST_NUMERATOR " << first_violation_score << "\n"
                  << "SCOPE EXACT_INTEGER_DUAL_PRICING_COMPLETE_RANK8_RANK9_"
                     "FIXED_FAILURES_NO_OPTIMIZATION\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PRICING_ERROR " << error.what() << '\n';
        return 1;
    }
}
