// Complete all-body witness multiplicity through two for the endpoint-591
// carrier. Unlike the protected-joint census, this includes joint witnesses.

#define ENDPOINT591_TWOHIT_MAIN endpoint591_twohit_hidden_for_all_witness
#include "endpoint591_twohit_census.cpp"
#undef ENDPOINT591_TWOHIT_MAIN

namespace {

struct LowWitness {
    u32 body = 0;
    unsigned hits = 0;
    u32 first = 0;
    u32 second = 0;
};

struct AllWitnessRow {
    u64 zero = 0;
    std::vector<LowWitness> low;
};

AllWitnessRow census_all_witnesses(PairAgent pair,
                                   const std::vector<u32>& carrier) {
    const Geometry geometry = build_geometry(pair.q, pair.r);
    std::vector<u32> active;
    for (u32 mask : carrier)
        if (margin(geometry, mask).ticks >= 0) active.push_back(mask);
    AllWitnessRow result;
    u64 bodies = 0;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++bodies;
        std::array<u32, 3> witnesses{};
        unsigned hits = 0;
        for (u32 mask : active) {
            if ((mask & body) != 0) continue;
            witnesses[hits++] = mask;
            if (hits == 3) break;
        }
        if (hits == 0) ++result.zero;
        if (hits <= 2) {
            if (hits == 2 && witnesses[1] < witnesses[0])
                std::swap(witnesses[0], witnesses[1]);
            result.low.push_back({body, hits, witnesses[0], witnesses[1]});
        }
    }
    require(bodies == kBodyCount617, "body universe changed");
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 19,
                "usage: all-witness JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "ADDITION593 DELETE593 TOP591 COVER43 DELETE43 LOW_CSV WORKERS");
        const unsigned worker_count = std::stoul(argv[18]);
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped range");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::vector<u32> old_carrier = reconstruct_c3925_593(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9], argv[10],
            argv[11]);
        (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        const auto add593 = read_flexible_masks593(argv[12], 1,
                                                   kAddition593Fnv);
        const auto del593 = read_flexible_masks593(argv[13], 1,
                                                   kDeletion593Fnv);
        std::vector<u32> carrier593;
        for (u32 mask : old_carrier)
            if (mask != del593.front()) carrier593.push_back(mask);
        carrier593.push_back(add593.front());
        const auto additions = read_cover591(argv[15]);
        const auto deletions = read_deletions591(argv[16]);
        const std::set<u32> deletion_set(deletions.begin(), deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : carrier593)
            if (!deletion_set.contains(mask)) carrier.push_back(mask);
        carrier.insert(carrier.end(), additions.begin(), additions.end());
        require(carrier.size() == 3925 &&
                    masks_fnv_agent(carrier) == kFinalCarrier591Fnv,
                "final carrier changed");
        for (u32 mask : joint)
            require(std::find(carrier.begin(), carrier.end(), mask) !=
                        carrier.end(), "joint mask lost");

        const auto row_set = read_pairs_repaired(argv[14], 13, kTop591Fnv);
        std::vector<PairAgent> rows;
        for (const auto& [q, r] : row_set) rows.push_back({q, r});
        std::vector<AllWitnessRow> audits(rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= rows.size()) return;
                    audits[index] = census_all_witnesses(rows[index], carrier);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
        };
        std::vector<std::thread> workers;
        for (unsigned index = 0; index < worker_count; ++index)
            workers.emplace_back(worker);
        for (auto& thread : workers) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        std::ofstream output(argv[17]);
        require(static_cast<bool>(output), "cannot create low ledger");
        output << "q,r,body_hex,hits,first_mask_hex,second_mask_hex\n";
        u64 total_zero = 0, total_one = 0, total_two = 0;
        Fnv ledger;
        std::map<std::pair<u32, u32>, u64> pair_load;
        for (std::size_t index = 0; index < rows.size(); ++index)
            for (const LowWitness& item : audits[index].low) {
                output << rows[index].q << ',' << rows[index].r << ','
                       << hex8(item.body) << ',' << item.hits << ',';
                if (item.hits >= 1) output << hex8(item.first);
                output << ',';
                if (item.hits >= 2) output << hex8(item.second);
                output << '\n';
                total_zero += item.hits == 0;
                total_one += item.hits == 1;
                total_two += item.hits == 2;
                ledger.add(rows[index].q); ledger.add(rows[index].r);
                ledger.add(item.body); ledger.add(item.hits);
                ledger.add(item.first); ledger.add(item.second);
                if (item.hits == 2) ++pair_load[{item.first, item.second}];
            }
        require(output.good(), "low ledger write failed");
        require(total_zero == 0, "carrier failure reappeared");
        u64 maximum_pair_load = 0;
        std::pair<u32, u32> least_maximum{};
        for (const auto& [pair, load] : pair_load)
            if (load > maximum_pair_load) {
                maximum_pair_load = load;
                least_maximum = pair;
            }

        std::cout << "LRC14_ENDPOINT591_ALL_WITNESS_CENSUS_V1\n"
                  << "CARRIER 3925 FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec
                  << " ROWS 13 BODY_TESTS " << 13 * kBodyCount617 << '\n';
        for (std::size_t index = 0; index < rows.size(); ++index) {
            u64 one = 0, two = 0;
            for (const auto& item : audits[index].low) {
                one += item.hits == 1;
                two += item.hits == 2;
            }
            std::cout << "ROW " << rows[index].q << ",591 ZERO "
                      << audits[index].zero << " ONE " << one << " TWO "
                      << two << '\n';
        }
        std::cout << "SUMMARY ZERO " << total_zero << " ONE " << total_one
                  << " TWO " << total_two << " LOW_FNV " << std::hex
                  << ledger.state << std::dec << " DISTINCT_TWO_PAIRS "
                  << pair_load.size() << " MAX_PAIR_LOAD " << maximum_pair_load
                  << " LEAST_MAX_PAIR " << hex8(least_maximum.first) << ','
                  << hex8(least_maximum.second) << '\n'
                  << "SCOPE COMPLETE_ALL_CARRIER_WITNESSES_ENDPOINT591_ONLY_"
                     "NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT591_ALL_WITNESS_ERROR " << error.what() << '\n';
        return 1;
    }
}
