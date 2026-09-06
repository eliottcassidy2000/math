// Independent complete K8 Hamilton base: Gray-code bitsets, no Walsh transform.
// Every root-gauged signing is visited. Cycle signs are toggled as packed bits.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

static uint64_t gates = 0;
static void need(bool ok, const char* what) {
    ++gates;
    if (!ok) throw std::runtime_error(what);
}

int main() {
    constexpr int n = 8, dim = 21, total = 2520, words = (total + 63) / 64;
    constexpr uint32_t end = uint32_t(1) << dim;
    int index[n][n] = {}, canonical[n][n] = {};
    int d = 0;
    // Different edge ordering from the primary: colexicographic in vertices.
    for (int v = 2; v < n; ++v)
        for (int u = 1; u < v; ++u) index[u][v] = index[v][u] = d++;
    d = 0;
    for (int u = 1; u < n; ++u)
        for (int v = u + 1; v < n; ++v) canonical[u][v] = canonical[v][u] = d++;
    std::array<uint32_t, dim> canonical_bit{};
    for (int v = 2; v < n; ++v)
        for (int u = 1; u < v; ++u)
            canonical_bit[index[u][v]] = uint32_t(1) << canonical[u][v];

    std::array<std::array<uint64_t, words>, dim> contains{};
    std::vector<int> order = {1,2,3,4,5,6,7};
    std::set<uint32_t> unique;
    int count = 0;
    do {
        if (order.front() > order.back()) continue;
        uint32_t mask = 0;
        for (int j = 0; j < n-2; ++j)
            mask |= uint32_t(1) << index[order[j]][order[j+1]];
        need(unique.insert(mask).second, "anchored cycle reversal quotient is injective");
        for (int bit = 0; bit < dim; ++bit)
            if (mask & (uint32_t(1) << bit))
                contains[bit][count / 64] |= uint64_t(1) << (count % 64);
        ++count;
    } while (std::next_permutation(order.begin(), order.end()));
    need(count == total, "complete unoriented Hamilton cycle universe");
    for (const auto& row : contains) {
        int incidence = 0;
        for (uint64_t block : row) incidence += __builtin_popcountll(block);
        need(incidence == 720, "every fixed nongauge edge lies in 6! cycles");
    }

    std::set<uint32_t> expected;
    for (int bit = 0; bit < dim; ++bit) expected.insert(uint32_t(1) << bit);
    for (int v = 1; v < n; ++v) {
        uint32_t star = 0;
        for (int u = 1; u < n; ++u)
            if (u != v) star |= uint32_t(1) << index[u][v];
        expected.insert(star);
    }
    auto single = expected;
    for (uint32_t mask : single) expected.insert((end-1) ^ mask);
    need(single.size() == 28 && expected.size() == 56, "labelled equality masks are distinct");

    std::array<uint64_t, words> negative{};
    std::vector<int16_t> spectrum(end);
    std::set<uint32_t> ties;
    int minimum = total, zero_count = 0;
    uint32_t canonical_mask = 0;
    for (uint32_t step = 0; step < end; ++step) {
        int weight = 0;
        if (step) {
            const int flip = __builtin_ctz(step);
            canonical_mask ^= canonical_bit[flip];
            for (int word = 0; word < words; ++word) {
                negative[word] ^= contains[flip][word];
                weight += __builtin_popcountll(negative[word]);
            }
        }
        const uint32_t mask = step ^ (step >> 1);
        spectrum[canonical_mask] = int16_t(total - 2*weight);
        if (!weight) {
            ++zero_count;
            need(mask == 0 || mask == end-1, "complete zero set is B,A");
            continue;
        }
        need(mask != 0 && mask != end-1, "B,A must have zero weight");
        if (weight < minimum) { minimum = weight; ties.clear(); }
        if (weight == minimum) ties.insert(mask);
    }
    need(zero_count == 2, "exactly two zero switching classes");
    need(minimum == 720, "complete nonzero K8 Hamilton minimum");
    need(ties == expected, "complete equality classification");
    uint64_t fingerprint = 1469598103934665603ULL;
    for (int16_t value : spectrum) {
        fingerprint ^= uint16_t(value);
        fingerprint *= 1099511628211ULL;
    }
    need(fingerprint == 0xe2dfba14125e7983ULL, "entire canonical spectrum matches frozen producer digest");
    std::cout << "method=packed_cycle_signs_gray_code; no_Walsh_or_repository_imports\n";
    std::cout << "n=8 switching_classes=" << end << " Hamilton_cycles=" << count
              << " minimum=" << minimum << " zero_classes=" << zero_count
              << " equality_classes=" << ties.size() << "\n";
    std::cout << "semantic_fnv64=" << std::hex << fingerprint << std::dec << "\n";
    std::cout << "spectrum_memory_bytes=" << spectrum.size()*sizeof(int16_t)
              << " cycle_incidence_memory_bytes=" << sizeof(contains) << "\n";
    std::cout << "checks=" << gates << "\nRESULT=PASS\n";
}
