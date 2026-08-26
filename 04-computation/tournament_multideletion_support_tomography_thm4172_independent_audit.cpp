#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using i64 = long long;
using u64 = std::uint64_t;

[[noreturn]] static void fail(const std::string& message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

static void need(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

static i64 choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  i64 answer = 1;
  for (int j = 1; j <= k; ++j) answer = answer * (n - j + 1) / j;
  return answer;
}

struct Fraction {
  i64 num = 0;
  i64 den = 1;

  Fraction(i64 a = 0, i64 b = 1) {
    need(b != 0, "zero denominator");
    if (b < 0) { a = -a; b = -b; }
    i64 g = std::gcd(a < 0 ? -a : a, b);
    num = a / g;
    den = b / g;
  }

  friend bool operator==(const Fraction& a, const Fraction& b) {
    return a.num == b.num && a.den == b.den;
  }

  std::string str() const {
    return std::to_string(num) + "/" + std::to_string(den);
  }
};

using Adj = std::vector<u64>;
using Tensor = std::vector<std::vector<i64>>;

static Adj decode(u64 code, int n) {
  Adj adj(n, 0);
  int bit = 0;
  for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j, ++bit) {
    if ((code >> bit) & 1ULL) adj[i] |= 1ULL << j;
    else adj[j] |= 1ULL << i;
  }
  return adj;
}

static u64 encode(const Adj& adj) {
  u64 code = 0;
  int bit = 0;
  for (int i = 0; i < (int)adj.size(); ++i) {
    for (int j = i + 1; j < (int)adj.size(); ++j, ++bit) {
      code |= ((adj[i] >> j) & 1ULL) << bit;
    }
  }
  return code;
}

struct Restriction {
  Adj adj;
  Tensor tensor;
  std::vector<int> keep;
  std::vector<int> where;
};

static Restriction restrict_set(const Adj& adj, const Tensor& z, u64 removed) {
  int n = (int)adj.size();
  Restriction answer;
  answer.where.assign(n, -1);
  for (int v = 0; v < n; ++v) if (!((removed >> v) & 1ULL)) {
    answer.where[v] = (int)answer.keep.size();
    answer.keep.push_back(v);
  }
  int m = (int)answer.keep.size();
  answer.adj.assign(m, 0);
  answer.tensor.assign(m, std::vector<i64>(m, 0));
  for (int a = 0; a < m; ++a) for (int b = a + 1; b < m; ++b) {
    int i = answer.keep[a], j = answer.keep[b];
    if ((adj[i] >> j) & 1ULL) answer.adj[a] |= 1ULL << b;
    else answer.adj[b] |= 1ULL << a;
    if (!z.empty()) answer.tensor[a][b] = answer.tensor[b][a] = z[i][j];
  }
  return answer;
}

static Tensor exposed_capacity(const Adj& adj) {
  int n = (int)adj.size();
  Tensor c(n, std::vector<i64>(n, 0));
  std::vector<int> word(n);
  std::iota(word.begin(), word.end(), 0);
  do {
    int reversed = 0, reversed_at = -1;
    for (int k = 0; k + 1 < n; ++k) {
      if (!((adj[word[k]] >> word[k + 1]) & 1ULL)) {
        ++reversed;
        reversed_at = k;
      }
    }
    if (reversed == 0) {
      for (int k = 0; k + 1 < n; ++k) {
        int i = word[k], j = word[k + 1];
        ++c[i][j]; ++c[j][i];
      }
    } else if (reversed == 1) {
      int i = word[reversed_at], j = word[reversed_at + 1];
      ++c[i][j]; ++c[j][i];
    }
  } while (std::next_permutation(word.begin(), word.end()));
  return c;
}

static i64 hamilton_paths(const Adj& adj) {
  int n = (int)adj.size();
  int size = 1 << n;
  std::vector<i64> dp((std::size_t)size * n, 0);
  for (int v = 0; v < n; ++v) dp[(std::size_t)(1 << v) * n + v] = 1;
  for (int mask = 1; mask < size; ++mask) {
    for (int last = 0; last < n; ++last) {
      i64 ways = dp[(std::size_t)mask * n + last];
      if (!ways) continue;
      u64 todo = adj[last] & (u64)(size - 1) & ~(u64)mask;
      while (todo) {
        u64 bit = todo & (~todo + 1);
        todo ^= bit;
        int next = __builtin_ctzll(bit);
        dp[(std::size_t)(mask | (int)bit) * n + next] += ways;
      }
    }
  }
  i64 answer = 0;
  for (int last = 0; last < n; ++last) answer += dp[(std::size_t)(size - 1) * n + last];
  return answer;
}

static Adj add_ear(const Adj& adj, u64 outgoing) {
  int n = (int)adj.size();
  Adj child = adj;
  child.push_back(0);
  for (int u = 0; u < n; ++u) {
    if ((outgoing >> u) & 1ULL) child[n] |= 1ULL << u;
    else child[u] |= 1ULL << n;
  }
  return child;
}

static Tensor response_capacity(const Adj& adj) {
  int n = (int)adj.size();
  Tensor c(n, std::vector<i64>(n, 0));
  i64 empty = hamilton_paths(add_ear(adj, 0));
  std::vector<i64> singleton(n);
  for (int i = 0; i < n; ++i) singleton[i] = hamilton_paths(add_ear(adj, 1ULL << i));
  for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
    i64 pair = hamilton_paths(add_ear(adj, (1ULL << i) | (1ULL << j)));
    c[i][j] = c[j][i] = singleton[i] + singleton[j] - empty - pair;
    need(c[i][j] >= 0, "capacity positivity");
  }
  return c;
}

static std::pair<i64, i64> packet(const Adj& adj, const Tensor& z) {
  int n = (int)adj.size();
  std::vector<i64> degree(n, 0), field(n, 0);
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (i != j) {
    degree[i] += z[i][j];
    field[i] += ((adj[i] >> j) & 1ULL) ? z[i][j] : -z[i][j];
  }
  i64 c_value = 0, d_value = 0;
  for (int i = 0; i < n; ++i) c_value += degree[i] * field[i];
  for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
    for (int k = i + 1; k < n; ++k) for (int l = k + 1; l < n; ++l) {
      if (i != k && i != l && j != k && j != l) d_value += z[i][j] * z[k][l];
    }
  }
  return {c_value, d_value};
}

static Fraction tilt(int n, i64 c_value, i64 d_value) {
  need(d_value > 0, "tilt denominator");
  return Fraction((n - 3) * c_value, (n % 2 ? 4 : 2) * d_value);
}

static std::pair<i64, i64> capacity_tomography() {
  std::vector<std::vector<Tensor>> cache(6);
  for (int n = 2; n <= 5; ++n) {
    int count = 1 << (n * (n - 1) / 2);
    cache[n].reserve(count);
    for (int code = 0; code < count; ++code) cache[n].push_back(exposed_capacity(decode(code, n)));
  }

  i64 rows = 0, gates = 0;
  for (int n = 2; n <= 5; ++n) {
    int count = 1 << (n * (n - 1) / 2);
    for (int code = 0; code < count; ++code) {
      Adj adj = decode(code, n);
      Tensor none;
      for (int x = 0; x < n; ++x) for (int y = x + 1; y < n; ++y) {
        std::vector<i64> layers(n - 1, 0);
        u64 outside = ((1ULL << n) - 1) ^ (1ULL << x) ^ (1ULL << y);
        for (u64 removed = outside;; removed = (removed - 1) & outside) {
          int r = __builtin_popcountll(removed);
          Restriction card = restrict_set(adj, none, removed);
          const Tensor& cc = cache[card.adj.size()][encode(card.adj)];
          layers[r] += cc[card.where[x]][card.where[y]];
          if (removed == 0) break;
        }
        need(layers[0] == cache[n][code][x][y], "layer zero");
        int N = n - 2;
        for (int s = 0; s <= N; ++s) {
          int degree = N - s;
          i64 weight = 0;
          for (int r = degree; r <= N; ++r) {
            i64 sign = ((r - degree) & 1) ? -1 : 1;
            weight += sign * choose(r, degree) * layers[r];
          }
          need(weight >= 0, "inverted atom weight");
        }
        ++rows;
        gates += layers.size();
      }
    }
  }
  need(rows == 10650 && gates == 42162, "tomography universe size");
  return {rows, gates};
}

static std::pair<i64, i64> multideletion_algebra() {
  i64 rows = 0, restrictions = 0;
  for (int n = 5; n <= 9; ++n) {
    Adj adj = decode(((1ULL << (n * (n - 1) / 2)) - 1) / 3, n);
    Tensor z(n, std::vector<i64>(n, 0));
    int index = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j, ++index) {
      z[i][j] = z[j][i] = 1 + (29 * index + 11 * n) % 31;
    }
    auto [parent_c, parent_d] = packet(adj, z);
    for (int r = 1; r <= n - 4; ++r) {
      i64 sum_c = 0, sum_d = 0, count = 0;
      for (u64 removed = 0; removed < (1ULL << n); ++removed) {
        if (__builtin_popcountll(removed) != r) continue;
        Restriction child = restrict_set(adj, z, removed);
        auto [child_c, child_d] = packet(child.adj, child.tensor);
        sum_c += child_c; sum_d += child_d; ++count;
      }
      need(sum_c == choose(n - 3, r) * parent_c, "C support degree");
      need(sum_d == choose(n - 4, r) * parent_d, "D support degree");
      int m = n - r;
      Fraction parent = tilt(n, parent_c, parent_d);
      Fraction mean((m - 3) * sum_c, (m % 2 ? 4 : 2) * sum_d);
      int kappa_parent = n % 2 ? 4 : 2;
      int kappa_child = m % 2 ? 4 : 2;
      Fraction transported(kappa_child * mean.num, kappa_parent * mean.den);
      need(parent == transported, "multideletion holonomy");
      if ((r & 1) == 0) need(parent == mean, "same-parity barycenter");
      ++rows; restrictions += count;
    }
  }
  need(rows == 15 && restrictions == 632, "algebra universe size");
  return {rows, restrictions};
}

static void named_controls() {
  Adj prime11 = decode(3169369058263173ULL, 11);
  Tensor z11 = response_capacity(prime11);
  auto [c11, d11] = packet(prime11, z11);
  Fraction tau11 = tilt(11, c11, d11);
  need(hamilton_paths(prime11) == 23685, "prime11 H");
  need(tau11 == Fraction(1055017002LL, 11090656697LL), "prime11 tilt");
  int central11 = 0;
  Fraction worst11;
  std::pair<int, int> worst_pair11{-1, -1};
  for (int i = 0; i < 11; ++i) for (int j = i + 1; j < 11; ++j) {
    Restriction child = restrict_set(prime11, z11, (1ULL << i) | (1ULL << j));
    auto [c, d] = packet(child.adj, child.tensor);
    Fraction t = tilt(9, c, d);
    if ((__int128)(t.num < 0 ? -t.num : t.num) < t.den) ++central11;
    Fraction absolute(t.num < 0 ? -t.num : t.num, t.den);
    if ((__int128)absolute.num * worst11.den > (__int128)worst11.num * absolute.den) {
      worst11 = absolute; worst_pair11 = {i, j};
    }
  }
  need(central11 == 55, "prime11 corrected pair centrality");
  need(worst_pair11 == std::pair<int, int>(7, 10), "prime11 worst pair label");
  need(worst11 == Fraction(1629373665LL, 9484374388LL), "prime11 worst pair tilt");

  Adj hostile12 = {3070, 3644, 3704, 3824, 4064, 4032,
                   3970, 3846, 3598, 1024, 2049, 512};
  Tensor z12 = response_capacity(hostile12);
  auto [c12, d12] = packet(hostile12, z12);
  Fraction tau12 = tilt(12, c12, d12);
  need(hamilton_paths(hostile12) == 27759, "hostile12 H");
  need(tau12 == Fraction(-53092739331LL, 40435524866LL), "hostile12 tilt");
  int central_one = 0;
  for (int i = 0; i < 12; ++i) {
    Restriction child = restrict_set(hostile12, z12, 1ULL << i);
    auto [c, d] = packet(child.adj, child.tensor);
    Fraction t = tilt(11, c, d);
    central_one += (__int128)(t.num < 0 ? -t.num : t.num) < t.den;
  }
  need(central_one == 12, "hostile12 one-deletion false negative");
  int central_two = 0;
  Fraction unique_central, worst12;
  std::pair<int, int> central_pair{-1, -1}, worst_pair12{-1, -1};
  for (int i = 0; i < 12; ++i) for (int j = i + 1; j < 12; ++j) {
    Restriction child = restrict_set(hostile12, z12, (1ULL << i) | (1ULL << j));
    auto [c, d] = packet(child.adj, child.tensor);
    Fraction t = tilt(10, c, d);
    Fraction absolute(t.num < 0 ? -t.num : t.num, t.den);
    if ((__int128)absolute.num < absolute.den) {
      ++central_two; unique_central = t; central_pair = {i, j};
    }
    if ((__int128)absolute.num * worst12.den > (__int128)worst12.num * absolute.den) {
      worst12 = absolute; worst_pair12 = {i, j};
    }
  }
  need(central_two == 1 && central_pair == std::pair<int, int>(9, 11),
       "hostile12 unique central pair");
  need(unique_central == Fraction(-4461429LL, 4143186457LL), "central pair tilt");
  need(worst_pair12 == std::pair<int, int>(9, 10), "hostile12 worst pair label");
  need(worst12 == Fraction(770665315LL, 458646486LL), "hostile12 worst pair tilt");

  std::cout << "prime11 tau=" << tau11.str() << " central_pairs=" << central11
            << " worst_pair=7,10 worst_abs=" << worst11.str() << '\n';
  std::cout << "hostile12 tau=" << tau12.str() << " central_one=" << central_one
            << " noncentral_pairs=" << 66 - central_two
            << " central_pair=9,11 central_tau=" << unique_central.str()
            << " worst_pair=9,10 worst_abs=" << worst12.str() << '\n';
}

int main() {
  auto [tomography_rows, tomography_gates] = capacity_tomography();
  std::cout << "capacity_tomography " << tomography_rows << " " << tomography_gates << '\n';
  auto [algebra_rows, restrictions] = multideletion_algebra();
  std::cout << "multideletion_algebra " << algebra_rows << " " << restrictions << '\n';
  named_controls();
  std::cout << "THM4172_INDEPENDENT_AUDIT_ACCEPT\n";
  return 0;
}
