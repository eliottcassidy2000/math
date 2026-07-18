// r4_finite_horn_tail_referee_codex_S66.cpp
//
// Independent exact referee for the published tail replay (core indices
// 160,...,219) of r4_finite_horn_v2_kps_S128c61.py.  Note that the old progress
// line labelled core 160 was printed after processing that core, so the old
// prefix and the published tail overlap at index 160; the coverage conclusion
// is unaffected, but adding their work counts double-counts that one core.
// This implementation does
// not enumerate quadruples in the primary program's fraction-sorted order.
// Instead it regards every killer as a hyperedge on the core-safe residue
// states and runs a memoised, uncovered-state-first set-cover search.  Thus it
// proves the stronger statement that NO set of at most four available killers
// covers the safe-state universe.
//
// Exact predicate.  A state (q,a), 15 <= q <= 40 and 1 <= a < q, is safe for
// speed v when min(v*a mod q, q-(v*a mod q)) >= ceil(q/14).  For a 9-speed
// core P, killers are k = 13*max(P)+1,...,399.  A four-killer family could be
// uncertified by this finite horn only if its four unsafe-state hyperedges
// covered every state safe for P.  The search below asks that cover question
// directly, with integer arithmetic and fixed-width bit masks.
//
// Tournament/carrier audit (AGENTS.md).
//   * Candidate vertices considered: runners, moduli, residue states (q,a),
//     killer support patterns, and wall-crossing events.
//   * The exact carrier selected is the bipartite incidence hypergraph between
//     killers and core-safe residue states.  It preserves the required
//     four-union predicate.  Quotienting to pairwise runner data destroys the
//     triple/quadruple overlap information and cannot certify the implication.
//   * For a diagnostic killer tournament, the pairwise observable is
//       D(i,j) = |K_i| - |K_j|,
//     where K_i is killer i's unsafe-state edge.  Orient i -> j when D(i,j)>0;
//     the tie gauge is smaller speed -> larger speed.  This is a transitive
//     tournament: score histogram {0,1,...,m-1}, no directed cycles, singleton
//     SCCs, and exactly one Hamiltonian path.  We report its edge flips against
//     the increasing-speed gauge.  This fingerprint preserves marginal edge
//     sizes but deliberately does NOT stand in for the exact hypergraph test.
//
// Standard library only.  Prints deterministic data only.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr int kFirstCore = 160;
constexpr int kPastLastCore = 220;
constexpr int kMaxWords = 11;  // sum_{q=15}^{40}(q-1) = 689 < 11*64

struct Mask {
  std::array<std::uint64_t, kMaxWords> word{};

  bool operator==(const Mask& other) const noexcept {
    return word == other.word;
  }
};

struct MaskHash {
  std::size_t operator()(const Mask& mask) const noexcept {
    std::uint64_t h = UINT64_C(1469598103934665603);
    for (std::uint64_t x : mask.word) {
      h ^= x;
      h *= UINT64_C(1099511628211);
      h ^= x >> 32;
      h *= UINT64_C(1099511628211);
    }
    return static_cast<std::size_t>(h);
  }
};

struct State {
  int q;
  int a;
};

struct Killer {
  int speed;
  Mask unsafe;
  int unsafe_count;
};

int circular_distance(int value, int modulus) {
  const int residue = value % modulus;
  return std::min(residue, modulus - residue);
}

bool empty(const Mask& mask, int words) {
  for (int i = 0; i < words; ++i) {
    if (mask.word[i] != 0) return false;
  }
  return true;
}

int cardinality(const Mask& mask, int words) {
  int result = 0;
  for (int i = 0; i < words; ++i) {
    result += std::popcount(mask.word[i]);
  }
  return result;
}

int intersection_cardinality(const Mask& lhs, const Mask& rhs, int words) {
  int result = 0;
  for (int i = 0; i < words; ++i) {
    result += std::popcount(lhs.word[i] & rhs.word[i]);
  }
  return result;
}

Mask subtract_edge(const Mask& residual, const Mask& unsafe, int words) {
  Mask result;
  for (int i = 0; i < words; ++i) {
    result.word[i] = residual.word[i] & ~unsafe.word[i];
  }
  return result;
}

void set_bit(Mask& mask, int index) {
  mask.word[index / 64] |= UINT64_C(1) << (index % 64);
}

bool has_bit(const Mask& mask, int index) {
  return ((mask.word[index / 64] >> (index % 64)) & UINT64_C(1)) != 0;
}

void fnv_integer(std::uint64_t& digest, std::uint64_t value) {
  for (int shift = 0; shift < 64; shift += 8) {
    digest ^= (value >> shift) & UINT64_C(0xff);
    digest *= UINT64_C(1099511628211);
  }
}

std::string hex64(std::uint64_t value) {
  std::ostringstream out;
  out << std::hex << std::setfill('0') << std::setw(16) << value;
  return out.str();
}

std::vector<std::array<int, 9>> make_cores() {
  std::vector<std::array<int, 9>> cores;
  std::array<int, 9> core{};
  const auto extend = [&](auto&& self, int depth, int next) -> void {
    if (depth == 9) {
      cores.push_back(core);
      return;
    }
    const int needed_after = 8 - depth;
    for (int value = next; value <= 12 - needed_after; ++value) {
      core[depth] = value;
      self(self, depth + 1, value + 1);
    }
  };
  extend(extend, 0, 1);
  return cores;
}

std::string core_string(const std::array<int, 9>& core) {
  std::ostringstream out;
  out << '[';
  for (int i = 0; i < 9; ++i) {
    if (i != 0) out << ',';
    out << core[i];
  }
  out << ']';
  return out.str();
}

class CoverSearch {
 public:
  CoverSearch(const std::vector<Killer>& killers,
              const std::vector<std::vector<int>>& killers_by_state,
              int state_count)
      : killers_(killers),
        killers_by_state_(killers_by_state),
        state_count_(state_count),
        words_((state_count + 63) / 64) {
    for (auto& table : failed_) table.reserve(4096);
  }

  bool run(std::vector<int>& cover) {
    Mask all;
    for (int state = 0; state < state_count_; ++state) set_bit(all, state);
    path_.clear();
    const bool found = search(all, 4);
    if (found) cover = path_;
    return found;
  }

  const std::array<std::uint64_t, 5>& nodes() const { return nodes_; }
  std::uint64_t memo_hits() const { return memo_hits_; }
  std::uint64_t bound_prunes() const { return bound_prunes_; }
  std::size_t memo_states() const {
    std::size_t total = 0;
    for (const auto& table : failed_) total += table.size();
    return total;
  }

 private:
  bool search(const Mask& residual, int remaining) {
    ++nodes_[remaining];
    if (empty(residual, words_)) return true;
    if (remaining == 0) return false;

    auto& memo = failed_[remaining];
    if (memo.find(residual) != memo.end()) {
      ++memo_hits_;
      return false;
    }

    const int residual_count = cardinality(residual, words_);
    std::array<int, 4> top_counts{};
    for (const Killer& killer : killers_) {
      const int count = intersection_cardinality(residual, killer.unsafe, words_);
      for (int pos = 0; pos < remaining; ++pos) {
        if (count > top_counts[pos]) {
          for (int move = remaining - 1; move > pos; --move) {
            top_counts[move] = top_counts[move - 1];
          }
          top_counts[pos] = count;
          break;
        }
      }
    }
    const int optimistic =
        std::accumulate(top_counts.begin(), top_counts.begin() + remaining, 0);
    if (optimistic < residual_count) {
      ++bound_prunes_;
      memo.insert(residual);
      return false;
    }

    int pivot = -1;
    std::size_t fewest_candidates = std::numeric_limits<std::size_t>::max();
    for (int state = 0; state < state_count_; ++state) {
      if (!has_bit(residual, state)) continue;
      const std::size_t candidate_count = killers_by_state_[state].size();
      if (candidate_count < fewest_candidates) {
        pivot = state;
        fewest_candidates = candidate_count;
      }
    }
    if (pivot < 0) return true;
    if (fewest_candidates == 0) {
      memo.insert(residual);
      return false;
    }

    std::vector<std::pair<int, int>> branches;
    branches.reserve(fewest_candidates);
    for (int killer_index : killers_by_state_[pivot]) {
      branches.emplace_back(
          -intersection_cardinality(residual, killers_[killer_index].unsafe, words_),
          killer_index);
    }
    std::sort(branches.begin(), branches.end(), [&](const auto& lhs, const auto& rhs) {
      if (lhs.first != rhs.first) return lhs.first < rhs.first;
      return killers_[lhs.second].speed < killers_[rhs.second].speed;
    });

    for (const auto& branch : branches) {
      const int killer_index = branch.second;
      const Mask next = subtract_edge(residual, killers_[killer_index].unsafe, words_);
      path_.push_back(killers_[killer_index].speed);
      if (search(next, remaining - 1)) return true;
      path_.pop_back();
    }

    memo.insert(residual);
    return false;
  }

  const std::vector<Killer>& killers_;
  const std::vector<std::vector<int>>& killers_by_state_;
  int state_count_;
  int words_;
  std::array<std::unordered_set<Mask, MaskHash>, 5> failed_;
  std::array<std::uint64_t, 5> nodes_{};
  std::uint64_t memo_hits_ = 0;
  std::uint64_t bound_prunes_ = 0;
  std::vector<int> path_;
};

}  // namespace

int main() {
  const std::vector<std::array<int, 9>> cores = make_cores();
  if (cores.size() != 220) {
    std::cerr << "internal error: expected 220 cores\n";
    return 2;
  }

  std::cout << "### r=4 finite-horn tail independent hypergraph referee ###\n";
  std::cout << "scope: lexicographic core indices [160,220), KB=400, q=15..40\n";
  std::cout << "scope_note: old prefix and published tail overlap at core 160\n";
  std::cout << "claim tested: no <=4 killer unsafe-edges cover the core-safe residue states\n";
  std::cout << "search: uncovered-state-first exact DPLL + memo + cardinality upper bound\n";
  std::cout << "tournament: D(i,j)=|K_i|-|K_j|; ties use increasing speed\n\n";

  std::uint64_t global_digest = UINT64_C(1469598103934665603);
  std::array<std::uint64_t, 5> total_nodes{};
  std::uint64_t total_memo_hits = 0;
  std::uint64_t total_bound_prunes = 0;
  std::uint64_t total_memo_states = 0;
  std::uint64_t total_flips = 0;
  std::uint64_t total_safe_states = 0;
  int min_safe_states = std::numeric_limits<int>::max();
  int max_safe_states = 0;
  int covers_found = 0;

  for (int core_index = kFirstCore; core_index < kPastLastCore; ++core_index) {
    const auto& core = cores[core_index];
    const int lower_killer = 13 * core.back() + 1;

    std::vector<State> states;
    for (int q = 15; q <= 40; ++q) {
      const int threshold = (q + 13) / 14;
      for (int a = 1; a < q; ++a) {
        bool safe = true;
        for (int speed : core) {
          if (circular_distance(speed * a, q) < threshold) {
            safe = false;
            break;
          }
        }
        if (safe) states.push_back({q, a});
      }
    }

    const int words = (static_cast<int>(states.size()) + 63) / 64;
    if (words > kMaxWords) {
      std::cerr << "internal error: mask too small\n";
      return 2;
    }

    std::vector<Killer> killers;
    killers.reserve(400 - lower_killer);
    for (int speed = lower_killer; speed < 400; ++speed) {
      Killer killer{speed, Mask{}, 0};
      for (int state_index = 0; state_index < static_cast<int>(states.size());
           ++state_index) {
        const State& state = states[state_index];
        const int threshold = (state.q + 13) / 14;
        if (circular_distance(speed * state.a, state.q) < threshold) {
          set_bit(killer.unsafe, state_index);
          ++killer.unsafe_count;
        }
      }
      killers.push_back(killer);
    }

    std::vector<std::vector<int>> killers_by_state(states.size());
    for (int killer_index = 0; killer_index < static_cast<int>(killers.size());
         ++killer_index) {
      for (int state_index = 0; state_index < static_cast<int>(states.size());
           ++state_index) {
        if (has_bit(killers[killer_index].unsafe, state_index)) {
          killers_by_state[state_index].push_back(killer_index);
        }
      }
    }

    std::uint64_t core_digest = UINT64_C(1469598103934665603);
    fnv_integer(core_digest, static_cast<std::uint64_t>(core_index));
    for (int speed : core) fnv_integer(core_digest, static_cast<std::uint64_t>(speed));
    fnv_integer(core_digest, static_cast<std::uint64_t>(states.size()));
    for (const State& state : states) {
      fnv_integer(core_digest, static_cast<std::uint64_t>(state.q));
      fnv_integer(core_digest, static_cast<std::uint64_t>(state.a));
    }
    for (const Killer& killer : killers) {
      fnv_integer(core_digest, static_cast<std::uint64_t>(killer.speed));
      for (int word = 0; word < words; ++word) {
        fnv_integer(core_digest, killer.unsafe.word[word]);
      }
    }
    fnv_integer(global_digest, core_digest);

    std::uint64_t flips = 0;
    for (std::size_t i = 0; i < killers.size(); ++i) {
      for (std::size_t j = i + 1; j < killers.size(); ++j) {
        if (killers[i].unsafe_count < killers[j].unsafe_count) ++flips;
      }
    }

    CoverSearch search(killers, killers_by_state, static_cast<int>(states.size()));
    std::vector<int> cover;
    const bool has_cover = search.run(cover);
    if (has_cover) ++covers_found;

    for (int remaining = 0; remaining <= 4; ++remaining) {
      total_nodes[remaining] += search.nodes()[remaining];
    }
    total_memo_hits += search.memo_hits();
    total_bound_prunes += search.bound_prunes();
    total_memo_states += search.memo_states();
    total_flips += flips;
    total_safe_states += states.size();
    min_safe_states = std::min(min_safe_states, static_cast<int>(states.size()));
    max_safe_states = std::max(max_safe_states, static_cast<int>(states.size()));

    std::cout << "core=" << core_index << " P=" << core_string(core)
              << " safe=" << states.size() << " killers=" << killers.size()
              << " digest=" << hex64(core_digest) << " nodes=[";
    for (int remaining = 4; remaining >= 0; --remaining) {
      if (remaining != 4) std::cout << ',';
      std::cout << search.nodes()[remaining];
    }
    std::cout << "] memo=" << search.memo_states()
              << " hits=" << search.memo_hits()
              << " bound=" << search.bound_prunes()
              << " tournament_flips=" << flips
              << " cover=" << (has_cover ? "YES" : "NO");
    if (has_cover) {
      std::cout << " example=[";
      for (std::size_t i = 0; i < cover.size(); ++i) {
        if (i != 0) std::cout << ',';
        std::cout << cover[i];
      }
      std::cout << ']';
    }
    std::cout << '\n';
  }

  std::cout << "\n### aggregate ###\n";
  std::cout << "cores_checked: " << (kPastLastCore - kFirstCore) << '\n';
  std::cout << "covers_by_at_most_four_killers: " << covers_found << '\n';
  std::cout << "safe_states: total=" << total_safe_states
            << " min=" << min_safe_states << " max=" << max_safe_states << '\n';
  std::cout << "incidence_digest_fnv1a64: " << hex64(global_digest) << '\n';
  std::cout << "dpll_nodes_remaining_[4,3,2,1,0]: [";
  for (int remaining = 4; remaining >= 0; --remaining) {
    if (remaining != 4) std::cout << ',';
    std::cout << total_nodes[remaining];
  }
  std::cout << "]\n";
  std::cout << "memo_states: " << total_memo_states << '\n';
  std::cout << "memo_hits: " << total_memo_hits << '\n';
  std::cout << "cardinality_bound_prunes: " << total_bound_prunes << '\n';
  std::cout << "tournament_edge_flips_vs_increasing_speed: " << total_flips << '\n';
  std::cout << "tournament_fingerprint_each_core: transitive; cycles=0; "
               "SCCs=singletons; Hamiltonian_paths=1\n";
  std::cout << "carrier_verdict: pairwise tournament loses overlap; exact hypergraph "
               "preserves four-cover predicate\n";
  std::cout << "small_modulus_filter: moot here, since even unrestricted killers "
               "have no four-cover\n";
  std::cout << "RESULT: "
            << (covers_found == 0
                    ? "ALL 60 TAIL CORES CERTIFIED (NO FOUR-COVER)"
                    : "FAILURE: A FOUR-COVER EXISTS")
            << '\n';
  std::cout << "DONE\n";
  return covers_found == 0 ? 0 : 1;
}
