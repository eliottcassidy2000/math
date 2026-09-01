#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

using u64 = std::uint64_t;

constexpr u64 kFnvOffset = UINT64_C(0xcbf29ce484222325);
constexpr u64 kFnvPrime = UINT64_C(0x100000001b3);

void fnv_add(u64& state, u64 value) {
    for (unsigned shift = 0; shift < 64; shift += 8) {
        state ^= (value >> shift) & 0xff;
        state *= kFnvPrime;
    }
}

struct Bits {
    u64 lo = 0;
    u64 hi = 0;
    bool operator==(const Bits&) const = default;
};

struct State {
    Bits uncovered;
    unsigned depth = 0;
    bool operator==(const State&) const = default;
};

struct StateHash {
    std::size_t operator()(const State& state) const noexcept {
        u64 x = state.uncovered.lo ^ std::rotl(state.uncovered.hi, 23)
                ^ (u64{state.depth} << 57);
        x ^= x >> 30;
        x *= UINT64_C(0xbf58476d1ce4e5b9);
        x ^= x >> 27;
        x *= UINT64_C(0x94d049bb133111eb);
        x ^= x >> 31;
        return static_cast<std::size_t>(x);
    }
};

unsigned weight(Bits bits) {
    return std::popcount(bits.lo) + std::popcount(bits.hi);
}

Bits intersect(Bits left, Bits right) {
    return {left.lo & right.lo, left.hi & right.hi};
}

Bits subtract(Bits left, Bits right) {
    return {left.lo & ~right.lo, left.hi & ~right.hi};
}

bool subset(Bits left, Bits right) {
    return (left.lo & ~right.lo) == 0 && (left.hi & ~right.hi) == 0;
}

bool contains(Bits bits, unsigned index) {
    return index < 64 ? ((bits.lo >> index) & 1U)
                      : ((bits.hi >> (index - 64)) & 1U);
}

u64 parse_hex(const std::string& text) {
    std::size_t used = 0;
    const u64 value = std::stoull(text, &used, 16);
    if (used != text.size()) throw std::runtime_error("bad hex word");
    return value;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t comma = line.find(',', begin);
        fields.push_back(line.substr(begin, comma - begin));
        if (comma == std::string::npos) return fields;
        begin = comma + 1;
    }
}

std::vector<Bits> read_signatures(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open signatures");
    std::string line;
    if (!std::getline(input, line) ||
        line != "w0,w1,weight,count8,least8,count9,least9")
        throw std::runtime_error("signature header changed");
    std::vector<Bits> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split(line);
        if (fields.size() != 7) throw std::runtime_error("signature row width");
        Bits bits{parse_hex(fields[0]), parse_hex(fields[1])};
        if (bits.hi >> 36 || weight(bits) != std::stoul(fields[2]))
            throw std::runtime_error("invalid signature bits/weight");
        result.push_back(bits);
    }
    if (!input.eof() || result.size() != 14368)
        throw std::runtime_error("signature count changed");
    std::sort(result.begin(), result.end(), [](Bits left, Bits right) {
        if (weight(left) != weight(right)) return weight(left) > weight(right);
        if (left.hi != right.hi) return left.hi < right.hi;
        return left.lo < right.lo;
    });
    for (std::size_t index = 1; index < result.size(); ++index)
        if (result[index] == result[index - 1])
            throw std::runtime_error("duplicate signature");
    return result;
}

std::array<unsigned, 100> read_dual(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open dual weights");
    std::string line;
    if (!std::getline(input, line) || line != "failure_ordinal,weight")
        throw std::runtime_error("dual header changed");
    std::array<unsigned, 100> result{};
    unsigned total = 0;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split(line);
        if (fields.size() != 2) throw std::runtime_error("dual row width");
        const unsigned index = std::stoul(fields[0]);
        const unsigned value = std::stoul(fields[1]);
        if (index >= 100 || result[index] != 0 || value == 0)
            throw std::runtime_error("invalid dual row");
        result[index] = value;
        total += value;
    }
    if (total != 22) throw std::runtime_error("dual total changed");
    return result;
}

struct Search {
    std::vector<Bits> sets;
    std::array<std::vector<unsigned>, 100> coverers;
    std::array<unsigned, 100> dual{};
    std::unordered_set<State, StateHash> dead;
    std::vector<unsigned> witness;
    u64 nodes = 0;
    u64 sum_bound_prunes = 0;
    u64 dual_bound_prunes = 0;
    u64 trace_fnv = kFnvOffset;

    unsigned dual_weight(Bits bits) const {
        unsigned total = 0;
        for (u64 work = bits.lo; work; work &= work - 1)
            total += dual[std::countr_zero(work)];
        for (u64 work = bits.hi; work; work &= work - 1)
            total += dual[64 + std::countr_zero(work)];
        return total;
    }

    bool possible_by_sum(Bits uncovered, unsigned depth) {
        std::array<unsigned, 8> top{};
        for (Bits response : sets) {
            const unsigned gain = weight(intersect(response, uncovered));
            if (gain <= top.back()) continue;
            auto position = std::upper_bound(top.begin(), top.end(), gain,
                                             std::greater<unsigned>());
            std::move_backward(position, top.end() - 1, top.end());
            *position = gain;
        }
        unsigned sum = 0;
        for (unsigned index = 0; index < depth; ++index) sum += top[index];
        return sum >= weight(uncovered);
    }

    bool dfs(Bits uncovered, unsigned depth) {
        ++nodes;
        fnv_add(trace_fnv, uncovered.lo);
        fnv_add(trace_fnv, uncovered.hi);
        fnv_add(trace_fnv, depth);
        if (uncovered.lo == 0 && uncovered.hi == 0) return true;
        if (depth == 0) return false;
        const State state{uncovered, depth};
        if (dead.contains(state)) return false;
        if (dual_weight(uncovered) > 3 * depth) {
            ++dual_bound_prunes;
            dead.insert(state);
            return false;
        }
        if (!possible_by_sum(uncovered, depth)) {
            ++sum_bound_prunes;
            dead.insert(state);
            return false;
        }

        unsigned pivot = 100;
        std::size_t least = sets.size() + 1;
        for (unsigned bit = 0; bit < 100; ++bit) {
            if (!contains(uncovered, bit)) continue;
            if (coverers[bit].size() < least) {
                least = coverers[bit].size();
                pivot = bit;
            }
        }
        if (pivot == 100) throw std::runtime_error("missing pivot");

        std::vector<std::pair<Bits, unsigned>> branches;
        for (unsigned index : coverers[pivot]) {
            const Bits gain = intersect(sets[index], uncovered);
            bool dominated = false;
            for (const auto& [kept, ignored] : branches)
                if (subset(gain, kept)) {
                    dominated = true;
                    break;
                }
            if (dominated) continue;
            branches.erase(std::remove_if(branches.begin(), branches.end(),
                [&](const auto& entry) { return subset(entry.first, gain); }),
                branches.end());
            branches.push_back({gain, index});
        }
        std::sort(branches.begin(), branches.end(), [&](const auto& left,
                                                        const auto& right) {
            const unsigned left_dual = dual_weight(left.first);
            const unsigned right_dual = dual_weight(right.first);
            if (left_dual != right_dual) return left_dual > right_dual;
            if (weight(left.first) != weight(right.first))
                return weight(left.first) > weight(right.first);
            return left.second < right.second;
        });
        for (const auto& [gain, index] : branches) {
            witness.push_back(index);
            if (dfs(subtract(uncovered, gain), depth - 1)) return true;
            witness.pop_back();
        }
        dead.insert(state);
        return false;
    }
};

int main(int argc, char** argv) {
    try {
        if (argc != 4)
            throw std::runtime_error("usage: search SIGNATURES DUAL DEPTH");
        const auto raw = read_signatures(argv[1]);
        Search search;
        search.dual = read_dual(argv[2]);
        u64 raw_fnv = kFnvOffset;
        u64 dual_fnv = kFnvOffset;
        for (unsigned value : search.dual) fnv_add(dual_fnv, value);
        for (Bits candidate : raw) {
            fnv_add(raw_fnv, candidate.lo);
            fnv_add(raw_fnv, candidate.hi);
            if (search.dual_weight(candidate) > 3)
                throw std::runtime_error("raw dual constraint violated");
        }
        if (raw_fnv != UINT64_C(0x98a18b04fa42e318))
            throw std::runtime_error("complete signature-set FNV changed");
        if (dual_fnv != UINT64_C(0x0f34be558e560d67))
            throw std::runtime_error("dual-weight FNV changed");
        for (Bits candidate : raw) {
            bool dominated = false;
            for (Bits kept : search.sets)
                if (subset(candidate, kept)) {
                    dominated = true;
                    break;
                }
            if (!dominated) search.sets.push_back(candidate);
        }
        if (search.sets.size() != 1165)
            throw std::runtime_error("maximal quotient count changed");
        u64 maximal_fnv = kFnvOffset;
        for (unsigned index = 0; index < search.sets.size(); ++index) {
            fnv_add(maximal_fnv, search.sets[index].lo);
            fnv_add(maximal_fnv, search.sets[index].hi);
            if (search.dual_weight(search.sets[index]) > 3)
                throw std::runtime_error("dual constraint violated");
            for (unsigned bit = 0; bit < 100; ++bit)
                if (contains(search.sets[index], bit))
                    search.coverers[bit].push_back(index);
        }
        if (maximal_fnv != UINT64_C(0x40cc79e6322deab0))
            throw std::runtime_error("maximal signature quotient FNV changed");
        const unsigned depth = std::stoul(argv[3]);
        if (depth != 8) throw std::runtime_error("audit depth must equal eight");
        for (const auto& choices : search.coverers)
            if (choices.empty()) throw std::runtime_error("uncovered obligation");
        const Bits full{UINT64_MAX, (UINT64_C(1) << 36) - 1};
        const bool satisfiable = search.dfs(full, depth);
        Bits covered{};
        for (unsigned index : search.witness) {
            covered.lo |= search.sets[index].lo;
            covered.hi |= search.sets[index].hi;
        }
        if (satisfiable && !(covered == full))
            throw std::runtime_error("internal witness failed");
        if (satisfiable || search.nodes != UINT64_C(7163197) ||
            search.dead.size() != 6775810 ||
            search.sum_bound_prunes != UINT64_C(943159) ||
            search.dual_bound_prunes != UINT64_C(5575606) ||
            search.trace_fnv != UINT64_C(0xc2354a48f73aef96))
            throw std::runtime_error("exact no-eight traversal identity changed");
        std::cout << "ENDPOINT590_EXACT_DEPTH_SEARCH_CPP_V1\n"
                  << "SIGNATURES " << raw.size() << " MAXIMAL "
                  << search.sets.size() << " DEPTH " << depth << '\n'
                  << "RAW_FNV " << std::hex << raw_fnv << " MAXIMAL_FNV "
                  << maximal_fnv << " DUAL_FNV " << dual_fnv << std::dec << '\n'
                  << "NODES " << search.nodes << " DEAD " << search.dead.size()
                  << " SUM_PRUNES " << search.sum_bound_prunes
                  << " DUAL_PRUNES " << search.dual_bound_prunes
                  << " TRACE_FNV " << std::hex << search.trace_fnv << std::dec
                  << '\n';
        if (satisfiable) {
            std::cout << "RESULT SAT SIZE " << search.witness.size() << '\n';
            for (unsigned index : search.witness)
                std::cout << "SELECT " << std::hex << std::setw(16)
                          << std::setfill('0') << search.sets[index].lo << ','
                          << std::setw(16) << search.sets[index].hi << std::dec
                          << " WEIGHT " << weight(search.sets[index]) << '\n';
        } else {
            std::cout << "RESULT UNSAT\n";
        }
        std::cout << "SCOPE EXACT_COMPLETE_SIGNATURE_QUOTIENT_DOMINANCE_"
                     "AND_DEPTH_FIRST_SEARCH\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SEARCH590_ERROR " << error.what() << '\n';
        return 1;
    }
}
