// Complete rank-8/rank-9 response census for the 2,468 endpoint-592 failures,
// plus a reproducible high-coverage candidate-pool cover.  The cover is only
// an upper bound: the reduced pool is not claimed to contain an exact optimum.

#define ENDPOINT593_RESPONSE_CAPACITY_MAIN endpoint593_response_hidden_for_endpoint592_pool
#include "04-computation/lrc14_endpoint593_response_exchange_thm4311/endpoint593_response_capacity.cpp"
#undef ENDPOINT593_RESPONSE_CAPACITY_MAIN

#include <queue>
#include <unordered_set>

namespace {
constexpr u64 kFailureFnv592 = UINT64_C(0x2209b8d6760280cc);
constexpr u64 kExchangedCarrierFnv592 = UINT64_C(0xc9e5faef52ca5707);
constexpr u64 kAddition1Fnv592 = UINT64_C(0x60873ef7a2b4ab90);
constexpr u64 kDelete1Fnv592 = UINT64_C(0x4c14214a64ec202c);
constexpr std::array<int, 7> kQ592 = {96, 100, 105, 192, 210, 256, 416};
constexpr std::array<std::size_t, 7> kCounts592 = {13, 3, 48, 1, 288, 2101, 14};
constexpr std::array<std::size_t, 7> kWords592 = {1, 1, 1, 1, 5, 33, 1};
constexpr std::size_t kTotalWords592 = 43;
constexpr std::size_t kTotalFailures592 = 2468;
constexpr std::size_t kGlobalKeep592 = 20000;
constexpr std::size_t kRowKeep592 = 3000;

using Bits592 = std::array<u64, kTotalWords592>;

struct Candidate592 {
    u32 mask = 0;
    Bits592 response{};
    unsigned weight = 0;
};

struct FailureUniverse592 {
    std::array<std::vector<u32>, 7> bodies;
    std::array<std::size_t, 7> word_offset{};
    std::array<std::array<Bits592, 30>, 7> incidence{};
    Bits592 valid{};
};

FailureUniverse592 read_failures592(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint592 failures");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "endpoint592 failure header changed");
    FailureUniverse592 universe;
    std::size_t offset = 0;
    for (std::size_t row = 0; row < 7; ++row) {
        universe.word_offset[row] = offset;
        offset += kWords592[row];
    }
    require(offset == kTotalWords592, "word partition changed");
    Fnv ledger;
    std::set<std::tuple<int, int, u32>> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0, r = 0;
        std::string token;
        fields >> q >> r >> token;
        require(fields && r == 592, "malformed endpoint592 failure");
        const auto found = std::find(kQ592.begin(), kQ592.end(), q);
        require(found != kQ592.end(), "unexpected endpoint592 failed row");
        const std::size_t row = found - kQ592.begin();
        const u32 body = parse_mask_agent(token);
        require(std::popcount(body) == 9 &&
                    distinct.insert({q, r, body}).second,
                "endpoint592 body rank/distinctness changed");
        universe.bodies[row].push_back(body);
        ledger.add(q); ledger.add(r); ledger.add(body);
    }
    require(input.eof() && distinct.size() == kTotalFailures592 &&
                ledger.state == kFailureFnv592,
            "endpoint592 obligation identity changed");
    for (std::size_t row = 0; row < 7; ++row) {
        require(universe.bodies[row].size() == kCounts592[row],
                "endpoint592 row count changed");
        const std::size_t base = universe.word_offset[row];
        for (std::size_t index = 0; index < universe.bodies[row].size(); ++index) {
            const std::size_t word = base + index / 64;
            const u64 bit = UINT64_C(1) << (index % 64);
            universe.valid[word] |= bit;
            u32 body = universe.bodies[row][index];
            while (body) {
                const unsigned label = std::countr_zero(body);
                universe.incidence[row][label][word] |= bit;
                body &= body - 1;
            }
        }
    }
    return universe;
}

Candidate592 response592(u32 mask,
                          const std::array<Geometry, 7>& geometry,
                          const FailureUniverse592& universe) {
    Candidate592 candidate;
    candidate.mask = mask;
    std::array<unsigned, 9> labels{};
    unsigned label_count = 0;
    for (u32 work = mask; work; work &= work - 1)
        labels[label_count++] = std::countr_zero(work);
    require(label_count == 8 || label_count == 9,
            "candidate escaped ranks 8/9");
    for (std::size_t row = 0; row < 7; ++row) {
        if (margin(geometry[row], mask).ticks < 0) continue;
        const std::size_t base = universe.word_offset[row];
        for (std::size_t local = 0; local < kWords592[row]; ++local) {
            const std::size_t word = base + local;
            u64 blocked = 0;
            for (unsigned index = 0; index < label_count; ++index)
                blocked |= universe.incidence[row][labels[index]][word];
            const u64 reply = universe.valid[word] & ~blocked;
            candidate.response[word] = reply;
            candidate.weight += std::popcount(reply);
        }
    }
    return candidate;
}

unsigned row_weight592(const Candidate592& candidate, std::size_t row,
                       const FailureUniverse592& universe) {
    unsigned answer = 0;
    const std::size_t base = universe.word_offset[row];
    for (std::size_t local = 0; local < kWords592[row]; ++local)
        answer += std::popcount(candidate.response[base + local]);
    return answer;
}

unsigned intersection_weight592(const Bits592& left, const Bits592& right) {
    unsigned answer = 0;
    for (std::size_t word = 0; word < kTotalWords592; ++word)
        answer += std::popcount(left[word] & right[word]);
    return answer;
}

bool empty592(const Bits592& bits) {
    return std::all_of(bits.begin(), bits.end(),
                       [](u64 word) { return word == 0; });
}

void subtract592(Bits592& left, const Bits592& right) {
    for (std::size_t word = 0; word < kTotalWords592; ++word)
        left[word] &= ~right[word];
}

using Heap592 = std::priority_queue<std::pair<unsigned, u32>,
                                    std::vector<std::pair<unsigned, u32>>,
                                    std::greater<std::pair<unsigned, u32>>>;

void retain592(Heap592& heap, std::size_t limit, unsigned score, u32 mask) {
    if (score == 0) return;
    const std::pair<unsigned, u32> entry{score, mask};
    if (heap.size() < limit) {
        heap.push(entry);
    } else if (entry > heap.top()) {
        heap.pop();
        heap.push(entry);
    }
}
}

#ifndef ENDPOINT592_RESPONSE_POOL_MAIN
#define ENDPOINT592_RESPONSE_POOL_MAIN main
#endif

int ENDPOINT592_RESPONSE_POOL_MAIN(int argc, char** argv) {
    try {
        require(argc == 18,
                "usage: pool JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE OLD_UNION "
                "REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE ADDITION1 DELETE1 "
                "FAILURES2468 COVER_OUT RESPONSE_OUT POOL_BITS_OUT");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> base = reconstruct_c3925_593(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9], argv[10],
            argv[11]);
        const auto addition = read_flexible_masks593(
            argv[12], 1, kAddition1Fnv592);
        const auto deletion = read_flexible_masks593(
            argv[13], 1, kDelete1Fnv592);
        std::vector<u32> carrier;
        for (u32 mask : base)
            if (mask != deletion.front()) carrier.push_back(mask);
        require(carrier.size() + 1 == base.size(), "deletion absent");
        carrier.push_back(addition.front());
        require(carrier.size() == 3925 &&
                    masks_fnv_agent(carrier) == kExchangedCarrierFnv592,
                "exchanged carrier changed");
        const std::unordered_set<u32> carrier_set(carrier.begin(), carrier.end());
        (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        const FailureUniverse592 failures = read_failures592(argv[14]);
        std::array<Geometry, 7> geometry;
        for (std::size_t row = 0; row < 7; ++row)
            geometry[row] = build_geometry(kQ592[row], 592);

        Heap592 global;
        std::array<Heap592, 7> per_row;
        std::array<u32, kTotalFailures592> first_response{};
        std::size_t missing_first = kTotalFailures592;
        std::array<u64, 2> responder_count{};
        std::array<Fnv, 2> responder_ledger;
        std::array<unsigned, 2> maximum_weight{};
        std::array<std::array<unsigned, 7>, 2> maximum_row_weight{};

        std::array<std::pair<std::size_t, std::size_t>, 7> ordinal_range{};
        std::size_t ordinal = 0;
        for (std::size_t row = 0; row < 7; ++row) {
            ordinal_range[row] = {ordinal, ordinal + kCounts592[row]};
            ordinal += kCounts592[row];
        }
        require(ordinal == kTotalFailures592, "ordinal partition changed");

        for (unsigned rank : {8U, 9U}) {
            const std::size_t slot = rank - 8;
            u64 mask_count = 0;
            const u32 limit = UINT32_C(1) << 30;
            for (u32 mask = (UINT32_C(1) << rank) - 1; mask < limit;
                 mask = next_combination(mask)) {
                ++mask_count;
                const Candidate592 candidate = response592(mask, geometry, failures);
                if (candidate.weight == 0) continue;
                require(!carrier_set.contains(mask),
                        "existing carrier mask responds to a frozen failure");
                ++responder_count[slot];
                responder_ledger[slot].add(mask);
                maximum_weight[slot] =
                    std::max(maximum_weight[slot], candidate.weight);
                retain592(global, kGlobalKeep592, candidate.weight, mask);
                std::size_t response_ordinal = 0;
                for (std::size_t row = 0; row < 7; ++row) {
                    const unsigned score = row_weight592(candidate, row, failures);
                    maximum_row_weight[slot][row] =
                        std::max(maximum_row_weight[slot][row], score);
                    retain592(per_row[row], kRowKeep592, score, mask);
                    if (missing_first != 0) {
                        const std::size_t base_word = failures.word_offset[row];
                        for (std::size_t index = 0; index < kCounts592[row]; ++index) {
                            const std::size_t global_index = response_ordinal + index;
                            if (first_response[global_index] != 0) continue;
                            if ((candidate.response[base_word + index / 64] >>
                                 (index % 64)) & 1U) {
                                first_response[global_index] = mask;
                                --missing_first;
                            }
                        }
                    }
                    response_ordinal += kCounts592[row];
                }
            }
            require(mask_count == (rank == 8 ? UINT64_C(5852925)
                                              : UINT64_C(14307150)),
                    "rank universe changed");
        }
        require(missing_first == 0, "some endpoint592 obligation has no response");

        std::set<u32> pool_masks;
        while (!global.empty()) {
            pool_masks.insert(global.top().second);
            global.pop();
        }
        for (Heap592& heap : per_row)
            while (!heap.empty()) {
                pool_masks.insert(heap.top().second);
                heap.pop();
            }
        for (u32 mask : first_response) pool_masks.insert(mask);
        std::vector<Candidate592> pool;
        pool.reserve(pool_masks.size());
        for (u32 mask : pool_masks)
            pool.push_back(response592(mask, geometry, failures));

        // Export the entire finite retained pool, not merely the subsequent
        // greedy choices.  One row is an exact 2,468-obligation response
        // signature represented by the same 43-word partition used above.
        std::ofstream pool_bits_out(argv[17]);
        require(static_cast<bool>(pool_bits_out),
                "cannot create pool signature output");
        pool_bits_out << "mask_hex,rank,total_weight";
        for (std::size_t word = 0; word < kTotalWords592; ++word)
            pool_bits_out << ",w" << word;
        pool_bits_out << '\n';
        for (const Candidate592& candidate : pool) {
            pool_bits_out << hex8(candidate.mask) << ','
                          << std::popcount(candidate.mask) << ','
                          << candidate.weight;
            for (u64 word : candidate.response)
                pool_bits_out << ',' << std::hex << std::setw(16)
                              << std::setfill('0') << word << std::dec;
            pool_bits_out << '\n';
        }
        require(pool_bits_out.good(), "pool signature output write failed");

        Bits592 uncovered = failures.valid;
        std::vector<Candidate592> cover;
        while (!empty592(uncovered)) {
            const Candidate592* best = nullptr;
            unsigned best_gain = 0;
            for (const Candidate592& candidate : pool) {
                const unsigned gain =
                    intersection_weight592(candidate.response, uncovered);
                if (gain > best_gain ||
                    (gain == best_gain && gain != 0 &&
                     candidate.mask < best->mask)) {
                    best = &candidate;
                    best_gain = gain;
                }
            }
            require(best != nullptr && best_gain != 0, "pool cover stalled");
            cover.push_back(*best);
            subtract592(uncovered, best->response);
        }

        // Delete any now-redundant greedy choices, scanning in reverse order.
        for (std::size_t index = cover.size(); index-- > 0;) {
            Bits592 union_without{};
            for (std::size_t other = 0; other < cover.size(); ++other) {
                if (other == index) continue;
                for (std::size_t word = 0; word < kTotalWords592; ++word)
                    union_without[word] |= cover[other].response[word];
            }
            bool complete = true;
            for (std::size_t word = 0; word < kTotalWords592; ++word)
                complete &= (union_without[word] & failures.valid[word]) ==
                            failures.valid[word];
            if (complete) cover.erase(cover.begin() + index);
        }

        Bits592 covered{};
        for (const Candidate592& candidate : cover)
            for (std::size_t word = 0; word < kTotalWords592; ++word)
                covered[word] |= candidate.response[word];
        for (std::size_t word = 0; word < kTotalWords592; ++word)
            require((covered[word] & failures.valid[word]) == failures.valid[word],
                    "reduced cover incomplete");

        std::ofstream cover_out(argv[15]);
        std::ofstream response_out(argv[16]);
        require(cover_out && response_out, "cannot create response outputs");
        cover_out << "mask_hex,rank,total_weight";
        for (int q : kQ592) cover_out << ",q" << q;
        cover_out << '\n';
        response_out << "mask_hex,q,r,body_hex\n";
        Fnv cover_ledger;
        for (const Candidate592& candidate : cover) {
            cover_ledger.add(candidate.mask);
            cover_out << hex8(candidate.mask) << ',' << std::popcount(candidate.mask)
                      << ',' << candidate.weight;
            for (std::size_t row = 0; row < 7; ++row)
                cover_out << ',' << row_weight592(candidate, row, failures);
            cover_out << '\n';
            for (std::size_t row = 0; row < 7; ++row) {
                const std::size_t base_word = failures.word_offset[row];
                for (std::size_t index = 0; index < kCounts592[row]; ++index)
                    if ((candidate.response[base_word + index / 64] >>
                         (index % 64)) & 1U)
                        response_out << hex8(candidate.mask) << ',' << kQ592[row]
                                     << ",592," << hex8(failures.bodies[row][index])
                                     << '\n';
            }
        }
        require(cover_out.good() && response_out.good(), "output write failed");

        std::cout << "LRC14_ENDPOINT592_RESPONSE_POOL_SCOUT_V1\n"
                  << "OBLIGATIONS " << kTotalFailures592 << " FNV " << std::hex
                  << kFailureFnv592 << std::dec << " ROWS 7\n";
        for (unsigned rank : {8U, 9U}) {
            const std::size_t slot = rank - 8;
            std::cout << "RANK " << rank << " RESPONDERS "
                      << responder_count[slot] << " RESPONDER_FNV " << std::hex
                      << responder_ledger[slot].state << std::dec
                      << " MAX_TOTAL " << maximum_weight[slot] << " MAX_BY_ROW";
            for (unsigned value : maximum_row_weight[slot])
                std::cout << ' ' << value;
            std::cout << '\n';
        }
        std::cout << "CANDIDATE_POOL " << pool.size() << " GLOBAL_KEEP "
                  << kGlobalKeep592 << " ROW_KEEP " << kRowKeep592
                  << " FIRST_RESPONSE_SIDECAR " << kTotalFailures592 << '\n'
                  << "GREEDY_REDUCED_COVER " << cover.size() << " FNV "
                  << std::hex << cover_ledger.state << std::dec << '\n';
        for (const Candidate592& candidate : cover) {
            std::cout << "SELECT " << hex8(candidate.mask) << " RANK "
                      << std::popcount(candidate.mask) << " WEIGHT "
                      << candidate.weight << " ROW_WEIGHTS";
            for (std::size_t row = 0; row < 7; ++row)
                std::cout << ' ' << row_weight592(candidate, row, failures);
            std::cout << '\n';
        }
        std::cout << "SCOPE COMPLETE_RESPONDER_CENSUS_BUT_REDUCED_POOL_"
                     "UPPER_BOUND_ONLY_FIXED_FAILURES_NO_EXACT_MINIMUM_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT592_RESPONSE_POOL_ERROR " << error.what() << '\n';
        return 1;
    }
}
