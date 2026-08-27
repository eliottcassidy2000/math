// Exact response-fibre census for a rank-eight deck on the labelled
// nine-subsets of a thirty-element pool.  Input is one hexadecimal mask per
// token.  The program has no repository-relative or /tmp dependency.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace {

using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr std::size_t DECK_SIZE = 421;
constexpr std::size_t RESPONSE_WORDS = 7;
constexpr std::size_t BODY_COUNT = 14307150;
constexpr u64 EXPECTED_DECK_FNV = UINT64_C(0x20d63dd42fe8150e);

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

struct Fnv {
    u64 state = UINT64_C(0xcbf29ce484222325);
    void add(u64 value) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (value >> (8 * byte)) & UINT64_C(0xff);
            state *= UINT64_C(0x100000001b3);
        }
    }
};

u64 mix64(u64 value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

std::vector<u32> read_deck(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid deck token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "deck mask rank/distinctness changed");
        deck.push_back(mask);
    }
    Fnv ledger;
    for (u32 mask : deck) ledger.add(mask);
    require(deck.size() == DECK_SIZE && ledger.state == EXPECTED_DECK_FNV,
            "deck identity changed");
    return deck;
}

u32 next_same_popcount(u32 value) {
    const u32 low = value & (0u - value);
    const u32 ripple = value + low;
    return ripple | (((value ^ ripple) >> 2) / low);
}

using Signature = std::array<u64, RESPONSE_WORDS>;

Signature response_signature(u32 body, const std::vector<u32>& deck) {
    Signature words{};
    for (std::size_t index = 0; index < deck.size(); ++index) {
        if ((body & deck[index]) == 0)
            words[index / 64] |= UINT64_C(1) << (index % 64);
    }
    return words;
}

std::pair<u64, u64> signature_hash(const Signature& words) {
    Fnv first;
    u64 second = UINT64_C(0x243f6a8885a308d3);
    for (std::size_t index = 0; index < words.size(); ++index) {
        first.add(words[index]);
        second = mix64(second ^ mix64(words[index] +
                         UINT64_C(0x9e3779b97f4a7c15) * (index + 1)));
    }
    return {first.state, second};
}

unsigned signature_size(const Signature& words) {
    unsigned count = 0;
    for (u64 word : words) count += std::popcount(word);
    return count;
}

unsigned singleton_index(const Signature& words) {
    for (unsigned word = 0; word < words.size(); ++word)
        if (words[word] != 0)
            return 64 * word + std::countr_zero(words[word]);
    return std::numeric_limits<unsigned>::max();
}

struct HashRow {
    u64 first = 0;
    u64 second = 0;
    u32 body = 0;
    u16 multiplicity = 0;
};
static_assert(sizeof(HashRow) == 24);

struct ExactRow {
    Signature words{};
    u32 body = 0;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: response-fibre-census JOINT_DECK");
        const std::vector<u32> deck = read_deck(argv[1]);

        std::vector<HashRow> rows(BODY_COUNT);
        u32 body = (UINT32_C(1) << 9) - 1;
        for (std::size_t index = 0; index < rows.size(); ++index) {
            require(body < (UINT32_C(1) << 30) && std::popcount(body) == 9,
                    "body enumeration escaped universe");
            rows[index].body = body;
            if (index + 1 < rows.size()) body = next_same_popcount(body);
        }
        require(rows.back().body == UINT32_C(0x3fe00000),
                "body enumeration endpoint changed");

        // Fix the lane count so the frozen transcript is byte-identical on
        // hosts with different reported hardware concurrency.
        constexpr unsigned thread_count = 8;
        std::vector<std::array<u64, DECK_SIZE + 1>> local_hist(thread_count);
        std::vector<std::array<u64, DECK_SIZE>> local_private(thread_count);
        std::vector<std::thread> workers;
        for (unsigned lane = 0; lane < thread_count; ++lane) {
            workers.emplace_back([&, lane]() {
                const std::size_t begin = rows.size() * lane / thread_count;
                const std::size_t end = rows.size() * (lane + 1) / thread_count;
                for (std::size_t index = begin; index < end; ++index) {
                    const Signature words =
                        response_signature(rows[index].body, deck);
                    const unsigned count = signature_size(words);
                    require(count <= DECK_SIZE, "response size overflow");
                    const auto [first, second] = signature_hash(words);
                    rows[index].first = first;
                    rows[index].second = second;
                    rows[index].multiplicity = static_cast<u16>(count);
                    ++local_hist[lane][count];
                    if (count == 1) {
                        const unsigned only = singleton_index(words);
                        require(only < DECK_SIZE, "bad singleton response");
                        ++local_private[lane][only];
                    }
                }
            });
        }
        for (auto& worker : workers) worker.join();

        std::array<u64, DECK_SIZE + 1> body_hist{};
        std::array<u64, DECK_SIZE> private_count{};
        for (unsigned lane = 0; lane < thread_count; ++lane) {
            for (std::size_t count = 0; count < body_hist.size(); ++count)
                body_hist[count] += local_hist[lane][count];
            for (std::size_t index = 0; index < private_count.size(); ++index)
                private_count[index] += local_private[lane][index];
        }
        u64 body_sum = 0;
        u64 incidence_sum = 0;
        unsigned minimum = DECK_SIZE + 1;
        unsigned maximum = 0;
        for (unsigned count = 0; count <= DECK_SIZE; ++count) {
            body_sum += body_hist[count];
            incidence_sum += static_cast<u64>(count) * body_hist[count];
            if (body_hist[count] != 0) {
                minimum = std::min(minimum, count);
                maximum = std::max(maximum, count);
            }
        }
        require(body_sum == BODY_COUNT && minimum == 1,
                "body coverage census changed");

        Fnv row_ledger;
        for (const HashRow& row : rows) {
            row_ledger.add(row.body);
            row_ledger.add(row.first);
            row_ledger.add(row.second);
            row_ledger.add(row.multiplicity);
        }

        std::sort(rows.begin(), rows.end(), [](const HashRow& left,
                                               const HashRow& right) {
            return std::tie(left.first, left.second, left.body) <
                   std::tie(right.first, right.second, right.body);
        });

        u64 hash_groups = 0;
        u64 exact_fibres = 0;
        u64 collision_splits = 0;
        u64 maximum_fibre_size = 0;
        u64 maximum_hash_group = 0;
        std::vector<u64> fibre_body_hist(BODY_COUNT + 1);
        Fnv fibre_ledger;
        std::size_t begin = 0;
        while (begin < rows.size()) {
            std::size_t end = begin + 1;
            while (end < rows.size() && rows[end].first == rows[begin].first &&
                   rows[end].second == rows[begin].second)
                ++end;
            ++hash_groups;
            maximum_hash_group = std::max<u64>(maximum_hash_group, end - begin);
            if (end - begin == 1) {
                ++exact_fibres;
                ++fibre_body_hist[1];
                maximum_fibre_size = std::max<u64>(maximum_fibre_size, 1);
                fibre_ledger.add(rows[begin].first);
                fibre_ledger.add(rows[begin].second);
                fibre_ledger.add(0);
                fibre_ledger.add(1);
            } else {
                std::vector<ExactRow> exact(end - begin);
                for (std::size_t offset = 0; offset < exact.size(); ++offset) {
                    exact[offset].body = rows[begin + offset].body;
                    exact[offset].words =
                        response_signature(exact[offset].body, deck);
                    require(signature_hash(exact[offset].words) ==
                                std::pair{rows[begin].first, rows[begin].second},
                            "response hash replay changed");
                }
                std::sort(exact.begin(), exact.end(),
                          [](const ExactRow& left, const ExactRow& right) {
                    return std::tie(left.words, left.body) <
                           std::tie(right.words, right.body);
                });
                std::size_t exact_begin = 0;
                u64 split_index = 0;
                while (exact_begin < exact.size()) {
                    std::size_t exact_end = exact_begin + 1;
                    while (exact_end < exact.size() &&
                           exact[exact_end].words == exact[exact_begin].words)
                        ++exact_end;
                    const u64 size = exact_end - exact_begin;
                    ++exact_fibres;
                    ++fibre_body_hist[size];
                    maximum_fibre_size = std::max(maximum_fibre_size, size);
                    fibre_ledger.add(rows[begin].first);
                    fibre_ledger.add(rows[begin].second);
                    fibre_ledger.add(split_index++);
                    fibre_ledger.add(size);
                    exact_begin = exact_end;
                }
                collision_splits += split_index - 1;
            }
            begin = end;
        }
        require(exact_fibres == hash_groups + collision_splits,
                "exact fibre accounting changed");

        u64 private_masks = 0;
        u64 private_bodies = 0;
        std::vector<std::size_t> redundant;
        Fnv private_ledger;
        for (std::size_t index = 0; index < deck.size(); ++index) {
            private_bodies += private_count[index];
            private_masks += private_count[index] != 0;
            if (private_count[index] == 0) redundant.push_back(index);
            private_ledger.add(index);
            private_ledger.add(deck[index]);
            private_ledger.add(private_count[index]);
        }
        require(private_masks == 420 && private_bodies == body_hist[1] &&
                    private_bodies == 3512 && redundant.size() == 1 &&
                    redundant.front() == 318 &&
                    deck[redundant.front()] == UINT32_C(0x003c900c) &&
                    private_ledger.state == UINT64_C(0x6afd7072b1c0127f),
                "private-fibre accounting changed");
        const u64 minimal_fibres = private_masks;

        std::cout << "JOINT421_RESPONSE_FIBRE_CENSUS_V1\n"
                  << "DECK " << deck.size() << " FNV " << std::hex
                  << EXPECTED_DECK_FNV << std::dec << " THREADS "
                  << thread_count << '\n'
                  << "BODIES " << body_sum << " INCIDENCES " << incidence_sum
                  << " MIN " << minimum << " MAX " << maximum
                  << " SINGLETON_BODIES " << body_hist[1]
                  << " RESPONSE_ROW_FNV " << std::hex << row_ledger.state
                  << std::dec << '\n';
        for (unsigned count = minimum; count <= maximum; ++count)
            if (body_hist[count] != 0)
                std::cout << "BODY_MULTIPLICITY " << count << " BODIES "
                          << body_hist[count] << '\n';
        std::cout << "HASH_PAIR_GROUPS " << hash_groups
                  << " EXACT_DISTINCT_FIBRES " << exact_fibres
                  << " HASH_COLLISION_SPLITS " << collision_splits
                  << " DUPLICATE_BODY_EXCESS " << body_sum - exact_fibres
                  << " MAX_BODIES_PER_FIBRE " << maximum_fibre_size
                  << " MAX_HASH_GROUP " << maximum_hash_group
                  << " FIBRE_FNV " << std::hex << fibre_ledger.state
                  << std::dec << '\n';
        for (std::size_t size = 1; size < fibre_body_hist.size(); ++size)
            if (fibre_body_hist[size] != 0)
                std::cout << "FIBRE_REALIZATIONS " << size
                          << " DISTINCT_FIBRES " << fibre_body_hist[size]
                          << '\n';
        std::cout << "PRIVATE_MASKS " << private_masks
                  << " PRIVATE_BODIES " << private_bodies
                  << " PRIVATE_FNV " << std::hex << private_ledger.state
                  << std::dec << '\n'
                  << "BODY_REDUNDANT_MASK INDEX " << redundant.front()
                  << " MASK " << std::hex << deck[redundant.front()]
                  << std::dec << '\n'
                  << "INCLUSION_MINIMAL_FIBRES " << minimal_fibres
                  << " ALL_SINGLETON YES\n"
                  << "CRITERION ACTIVE_ALL_" << private_masks
                  << "_BODY_ESSENTIAL_MASKS IFF_DECK_CERTIFIES_ALL_BODIES\n"
                  << "VERDICT PASS EXACT_RESPONSE_FIBRE_CENSUS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "JOINT421_RESPONSE_FIBRE_CENSUS_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
