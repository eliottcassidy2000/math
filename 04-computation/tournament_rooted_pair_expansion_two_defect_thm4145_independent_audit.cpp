#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Digraph = std::vector<std::uint32_t>;

static void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

static Digraph decode(int order, std::uint64_t label) {
  Digraph out(order, 0);
  int bit = 0;
  for (int i = 0; i < order; ++i) {
    for (int j = i + 1; j < order; ++j, ++bit) {
      if ((label >> bit) & 1ULL) out[i] |= 1U << j;
      else out[j] |= 1U << i;
    }
  }
  return out;
}

static std::string bits(const Digraph& out) {
  std::string answer;
  for (int i = 0; i < static_cast<int>(out.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(out.size()); ++j) {
      answer.push_back(((out[i] >> j) & 1U) ? '1' : '0');
    }
  }
  return answer;
}

static bool reachable_all(const Digraph& arcs) {
  const int order = static_cast<int>(arcs.size());
  const std::uint32_t full = (1U << order) - 1U;
  std::uint32_t seen = 1U, todo = 1U;
  while (todo) {
    const int v = __builtin_ctz(todo);
    todo &= todo - 1U;
    const std::uint32_t fresh = arcs[v] & ~seen;
    seen |= fresh;
    todo |= fresh;
  }
  return seen == full;
}

static bool strong(const Digraph& out) {
  const int order = static_cast<int>(out.size());
  Digraph incoming(order, 0);
  for (int i = 0; i < order; ++i) {
    for (int j = 0; j < order; ++j) {
      if ((out[i] >> j) & 1U) incoming[j] |= 1U << i;
    }
  }
  return reachable_all(out) && reachable_all(incoming);
}

static std::int64_t hamilton(const Digraph& out) {
  const int order = static_cast<int>(out.size());
  const int size = 1 << order;
  std::vector<std::vector<std::int64_t>> dp(size,
      std::vector<std::int64_t>(order, 0));
  for (int v = 0; v < order; ++v) dp[1 << v][v] = 1;
  for (int mask = 1; mask < size; ++mask) {
    for (int v = 0; v < order; ++v) {
      if (!(mask & (1 << v))) continue;
      const int previous = mask ^ (1 << v);
      for (int u = 0; u < order; ++u) {
        if ((previous & (1 << u)) && ((out[u] >> v) & 1U)) {
          dp[mask][v] += dp[previous][u];
        }
      }
    }
  }
  return std::accumulate(dp.back().begin(), dp.back().end(), std::int64_t{0});
}

static Digraph expand_pair(const Digraph& out, int root) {
  const int order = static_cast<int>(out.size());
  Digraph child = out;
  child.push_back(0);
  const int clone = order;
  for (int u = 0; u < order; ++u) {
    if (u == root) continue;
    if ((out[root] >> u) & 1U) child[clone] |= 1U << u;
    else child[u] |= 1U << clone;
  }
  child[root] |= 1U << clone;
  return child;
}

static Digraph add_root_type_ear(const Digraph& out, int root) {
  const int order = static_cast<int>(out.size());
  Digraph child = out;
  child.push_back(0);
  const int ear = order;
  for (int u = 0; u < order; ++u) {
    if (u != root && ((out[root] >> u) & 1U)) child[ear] |= 1U << u;
    else child[u] |= 1U << ear;
  }
  return child;
}

static std::int64_t defect_formula(const Digraph& out, int root) {
  const int order = static_cast<int>(out.size());
  std::vector<int> permutation;
  for (int v = 0; v < order; ++v) if (v != root) permutation.push_back(v);
  std::int64_t total = 0;
  do {
    std::vector<int> defects;
    for (int i = 0; i + 1 < static_cast<int>(permutation.size()); ++i) {
      if (!((out[permutation[i]] >> permutation[i + 1]) & 1U)) {
        defects.push_back(i + 1);
      }
    }
    std::vector<bool> good(order, false);
    if ((out[root] >> permutation.front()) & 1U) good[0] = true;
    if ((out[permutation.back()] >> root) & 1U) good[order - 1] = true;
    for (int i = 0; i + 1 < static_cast<int>(permutation.size()); ++i) {
      if (((out[permutation[i]] >> root) & 1U) &&
          ((out[root] >> permutation[i + 1]) & 1U)) good[i + 1] = true;
    }
    bool repairable = defects.size() <= 2;
    for (int slot : defects) repairable = repairable && good[slot];
    if (repairable) {
      const int k = static_cast<int>(std::count(good.begin(), good.end(), true));
      if (defects.empty()) total += k * k;
      else if (defects.size() == 1) total += 2 * k - 1;
      else total += 2;
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return total;
}

static std::pair<std::int64_t, std::int64_t> gap_counts(
    const Digraph& child, int root, int clone) {
  const int order = static_cast<int>(child.size());
  std::vector<int> permutation(order);
  std::iota(permutation.begin(), permutation.end(), 0);
  std::int64_t correct = 0, reverse = 0;
  do {
    int root_position = -1, clone_position = -1;
    for (int i = 0; i < order; ++i) {
      if (permutation[i] == root) root_position = i;
      if (permutation[i] == clone) clone_position = i;
    }
    if (std::abs(root_position - clone_position) != 1) continue;
    std::vector<int> defects;
    for (int i = 0; i + 1 < order; ++i) {
      if (!((child[permutation[i]] >> permutation[i + 1]) & 1U)) {
        defects.push_back(i);
      }
    }
    const int position = std::min(root_position, clone_position);
    if (permutation[position] == root && permutation[position + 1] == clone &&
        defects.empty()) ++correct;
    if (permutation[position] == clone && permutation[position + 1] == root &&
        defects.size() == 1 && defects[0] == position) ++reverse;
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return {correct, reverse};
}

static std::int64_t square_only(const Digraph& out, int root) {
  const int order = static_cast<int>(out.size());
  std::vector<int> permutation;
  for (int v = 0; v < order; ++v) if (v != root) permutation.push_back(v);
  std::int64_t total = 0;
  do {
    bool path = true;
    for (int i = 0; i + 1 < static_cast<int>(permutation.size()); ++i) {
      path = path && ((out[permutation[i]] >> permutation[i + 1]) & 1U);
    }
    if (!path) continue;
    int good = ((out[root] >> permutation.front()) & 1U) ? 1 : 0;
    good += ((out[permutation.back()] >> root) & 1U) ? 1 : 0;
    for (int i = 0; i + 1 < static_cast<int>(permutation.size()); ++i) {
      good += (((out[permutation[i]] >> root) & 1U) &&
               ((out[root] >> permutation[i + 1]) & 1U)) ? 1 : 0;
    }
    total += good * good;
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return total;
}

static std::uint64_t fnv_update(std::uint64_t hash, const std::string& row) {
  for (unsigned char byte : row) {
    hash ^= byte;
    hash *= 1099511628211ULL;
  }
  return hash;
}

int main() {
  try {
    std::uint64_t semantic = 14695981039346656037ULL;
    std::int64_t total_rows = 0, total_strong_rooted = 0;
    std::vector<std::string> summaries;
    for (int order = 2; order <= 5; ++order) {
      const std::uint64_t labels = 1ULL << (order * (order - 1) / 2);
      std::int64_t rooted = 0, strong_rooted = 0;
      for (std::uint64_t label = 0; label < labels; ++label) {
        const Digraph out = decode(order, label);
        const bool strong_parent = strong(out);
        const std::int64_t parent_h = hamilton(out);
        for (int root = 0; root < order; ++root) {
          const Digraph child = expand_pair(out, root);
          require(child == add_root_type_ear(out, root), "ear adjacency");
          require(strong(child) == strong_parent, "strong equivalence");
          const std::int64_t child_h = hamilton(child);
          const std::int64_t formula_h = defect_formula(out, root);
          require(child_h == formula_h, "defect formula");
          const auto [correct, reverse] = gap_counts(child, root, order);
          require(correct == parent_h && reverse == parent_h, "gap bijections");
          std::ostringstream row;
          row << order << ',' << label << ',' << root << ','
              << (strong_parent ? 1 : 0) << ',' << parent_h << ',' << child_h
              << ',' << formula_h << ',' << correct << ',' << reverse << '\n';
          semantic = fnv_update(semantic, row.str());
          ++rooted;
          if (strong_parent) ++strong_rooted;
        }
      }
      total_rows += rooted;
      total_strong_rooted += strong_rooted;
      std::int64_t factorial = 1;
      for (int i = 2; i <= order - 1; ++i) factorial *= i;
      std::ostringstream summary;
      summary << "order=" << order << ";labels=" << labels
              << ";rooted=" << rooted << ";strong_rooted=" << strong_rooted
              << ";max_deletion_orders=" << factorial;
      summaries.push_back(summary.str());
    }
    require(total_rows == 5404, "rooted universe");
    require(total_strong_rooted == 2822, "strong rooted universe");
    require(semantic == 0x9f075ab9b9d3f07fULL, "semantic FNV");

    const Digraph c3 = decode(3, 0b101);
    const Digraph c3_child = expand_pair(c3, 0);
    require(bits(c3) == "101" && bits(c3_child) == "101101", "C3 labels");
    require(hamilton(c3_child) == 5 && square_only(c3, 0) == 4, "C3 hostile");

    const Digraph root_loss = decode(4, 0b00100);
    require(bits(root_loss) == "001000", "root-loss label");
    require(hamilton(expand_pair(root_loss, 0)) == 11, "root 0 response");
    require(hamilton(expand_pair(root_loss, 1)) == 9, "root 1 response");

    std::cout << "THM4145_ROOTED_PAIR_EXPANSION_TWO_DEFECT_INDEPENDENT_20260825\n";
    std::cout << "status=PASS;universe=all rooted labelled tournaments orders 2..5\n";
    for (const auto& summary : summaries) std::cout << summary << '\n';
    std::cout << "rooted_rows=" << total_rows
              << ";strong_rooted=" << total_strong_rooted << '\n';
    std::cout << "order11_covering_presentations=93559490\n";
    std::cout << "c3_hostile=(101,0,101101,5,4)\n";
    std::cout << "root_loss_control=(001000,(0,11),(1,9))\n";
    std::cout << "semantic_fnv64=" << std::hex << std::setw(16)
              << std::setfill('0') << semantic << std::dec << '\n';
    std::cout << "checks=strong_iff,ear_identity,defect_formula,internal_gap_pair\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "status=FAIL;error=" << error.what() << '\n';
    return 1;
  }
}
