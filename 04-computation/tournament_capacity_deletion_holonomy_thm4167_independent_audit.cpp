#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
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
  friend bool operator<(const Fraction& a, const Fraction& b) {
    return (__int128)a.num * b.den < (__int128)b.num * a.den;
  }
  friend Fraction operator+(const Fraction& a, const Fraction& b) {
    return Fraction((i64)((__int128)a.num * b.den + (__int128)b.num * a.den),
                    a.den * b.den);
  }
  friend Fraction operator*(const Fraction& a, const Fraction& b) {
    return Fraction(a.num * b.num, a.den * b.den);
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

static std::tuple<Adj, std::vector<int>, std::vector<int>>
delete_vertex(const Adj& adj, int v) {
  int n = (int)adj.size();
  std::vector<int> keep, where(n, -1);
  for (int u = 0; u < n; ++u) if (u != v) {
    where[u] = (int)keep.size();
    keep.push_back(u);
  }
  Adj card(n - 1, 0);
  for (int a = 0; a < n - 1; ++a) for (int b = a + 1; b < n - 1; ++b) {
    int i = keep[a], j = keep[b];
    if ((adj[i] >> j) & 1ULL) card[a] |= 1ULL << b;
    else card[b] |= 1ULL << a;
  }
  return {card, keep, where};
}

static i64 hamilton_paths(const Adj& adj) {
  int n = (int)adj.size();
  int size = 1 << n;
  std::vector<i64> dp((std::size_t)size * n, 0);
  for (int v = 0; v < n; ++v) dp[((1 << v) * n) + v] = 1;
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
  std::vector<i64> one(n);
  for (int i = 0; i < n; ++i) one[i] = hamilton_paths(add_ear(adj, 1ULL << i));
  for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
    i64 pair = hamilton_paths(add_ear(adj, (1ULL << i) | (1ULL << j)));
    c[i][j] = c[j][i] = one[i] + one[j] - empty - pair;
    need(c[i][j] >= 0, "response capacity sign");
  }
  return c;
}

static Tensor exposed_capacity(const Adj& adj) {
  int n = (int)adj.size();
  Tensor c(n, std::vector<i64>(n, 0));
  std::vector<int> word(n);
  std::iota(word.begin(), word.end(), 0);
  do {
    int bad = 0, bad_at = -1;
    for (int k = 0; k + 1 < n; ++k) {
      if (!((adj[word[k]] >> word[k + 1]) & 1ULL)) {
        ++bad;
        bad_at = k;
      }
    }
    if (bad == 0) {
      for (int k = 0; k + 1 < n; ++k) {
        int a = word[k], b = word[k + 1];
        ++c[a][b]; ++c[b][a];
      }
    } else if (bad == 1) {
      int a = word[bad_at], b = word[bad_at + 1];
      ++c[a][b]; ++c[b][a];
    }
  } while (std::next_permutation(word.begin(), word.end()));
  return c;
}

static std::pair<i64, i64> packet(const Adj& adj, const Tensor& z) {
  int n = (int)adj.size();
  std::vector<i64> degree(n, 0), field(n, 0);
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (i != j) {
    degree[i] += z[i][j];
    field[i] += ((adj[i] >> j) & 1ULL) ? z[i][j] : -z[i][j];
  }
  i64 cval = 0, dval = 0;
  for (int i = 0; i < n; ++i) cval += field[i] * degree[i];
  std::vector<std::pair<int, int>> edges;
  for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) edges.push_back({i, j});
  for (std::size_t a = 0; a < edges.size(); ++a) {
    auto [i, j] = edges[a];
    for (std::size_t b = a + 1; b < edges.size(); ++b) {
      auto [k, l] = edges[b];
      if (i != k && i != l && j != k && j != l) dval += z[i][j] * z[k][l];
    }
  }
  return {cval, dval};
}

static std::pair<Adj, Tensor> restrict_tensor(const Adj& adj, const Tensor& z, int v) {
  auto [card, keep, where] = delete_vertex(adj, v);
  Tensor out(keep.size(), std::vector<i64>(keep.size(), 0));
  for (std::size_t a = 0; a < keep.size(); ++a) for (std::size_t b = a + 1; b < keep.size(); ++b) {
    out[a][b] = out[b][a] = z[keep[a]][keep[b]];
  }
  return {card, out};
}

static std::vector<std::vector<int>> exposed_words(const Adj& adj, int x, int y) {
  int n = (int)adj.size();
  std::vector<int> word(n);
  std::iota(word.begin(), word.end(), 0);
  std::vector<std::vector<int>> answer;
  do {
    int marked = 0;
    bool good = true;
    for (int k = 0; k + 1 < n; ++k) {
      bool is_marked = (word[k] == x && word[k + 1] == y)
                    || (word[k] == y && word[k + 1] == x);
      marked += is_marked;
      if (!is_marked && !((adj[word[k]] >> word[k + 1]) & 1ULL)) good = false;
    }
    if (good && marked == 1) answer.push_back(word);
  } while (std::next_permutation(word.begin(), word.end()));
  return answer;
}

static int marked_audit_through_five() {
  int rows = 0;
  for (int n = 3; n <= 5; ++n) {
    int bits = n * (n - 1) / 2;
    for (u64 code = 0; code < (1ULL << bits); ++code) {
      Adj adj = decode(code, n);
      for (int v = 0; v < n; ++v) {
        auto [card, keep, where] = delete_vertex(adj, v);
        for (int x = 0; x < n; ++x) for (int y = x + 1; y < n; ++y) if (x != v && y != v) {
          auto local = exposed_words(card, where[x], where[y]);
          std::set<std::vector<int>> card_words;
          i64 redundancy = 0, theft = 0, extensions = 0;
          for (const auto& p : local) {
            std::vector<int> word;
            for (int u : p) word.push_back(keep[u]);
            card_words.insert(word);
            std::vector<int> signature;
            for (int u : word) signature.push_back((adj[v] >> u) & 1ULL);
            int drops = 0, mark = -1;
            for (int k = 0; k + 1 < n - 1; ++k) {
              drops += signature[k] == 1 && signature[k + 1] == 0;
              if ((word[k] == x && word[k + 1] == y) || (word[k] == y && word[k + 1] == x)) mark = k;
            }
            need(mark >= 0, "independent marked gap");
            int stolen = signature[mark] == 0 && signature[mark + 1] == 1;
            int literal = ((adj[v] >> word.front()) & 1ULL)
                        + ((adj[word.back()] >> v) & 1ULL);
            for (int k = 0; k + 1 < n - 1; ++k) if (k != mark) {
              literal += ((adj[word[k]] >> v) & 1ULL) && ((adj[v] >> word[k + 1]) & 1ULL);
            }
            need(literal == 1 + drops - stolen, "independent insertion formula");
            redundancy += drops; theft += stolen; extensions += literal;
          }
          auto full = exposed_words(adj, x, y);
          i64 orphans = 0;
          for (const auto& word : full) {
            std::vector<int> deleted;
            for (int u : word) if (u != v) deleted.push_back(u);
            orphans += !card_words.count(deleted);
          }
          need((i64)full.size() == extensions + orphans, "independent deletion fibres");
          need((i64)full.size() - (i64)local.size() == redundancy + orphans - theft,
               "independent marked identity");
          ++rows;
        }
      }
    }
  }
  return rows;
}

static std::tuple<i64, i64, i64> monotonicity_n6() {
  std::vector<Tensor> card_cache(1 << 10);
  for (u64 code = 0; code < (1ULL << 10); ++code) card_cache[code] = exposed_capacity(decode(code, 5));
  i64 gates = 0, zeros = 0, minimum = -1;
  for (u64 code = 0; code < (1ULL << 15); ++code) {
    Adj adj = decode(code, 6);
    Tensor c = exposed_capacity(adj);
    for (int v = 0; v < 6; ++v) {
      auto [card, keep, where] = delete_vertex(adj, v);
      u64 card_code = 0;
      int bit = 0;
      for (int i = 0; i < 5; ++i) for (int j = i + 1; j < 5; ++j, ++bit) {
        card_code |= ((card[i] >> j) & 1ULL) << bit;
      }
      const Tensor& cc = card_cache[card_code];
      for (int a = 0; a < 5; ++a) for (int b = a + 1; b < 5; ++b) {
        i64 delta = c[keep[a]][keep[b]] - cc[a][b];
        if (minimum < 0 || delta < minimum) minimum = delta;
        zeros += delta == 0;
        ++gates;
      }
    }
  }
  need(minimum == 0, "independent n6 monotonicity");
  return {gates, zeros, minimum};
}

static int restriction_algebra() {
  int rows = 0;
  for (int n = 5; n <= 8; ++n) {
    Adj adj = decode((1ULL << (n * (n - 1) / 2)) / 3, n);
    Tensor z(n, std::vector<i64>(n, 0));
    int index = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j, ++index) {
      z[i][j] = z[j][i] = 1 + (17 * index + 5 * n) % 23;
    }
    auto [cval, dval] = packet(adj, z);
    i64 sum_c = 0, sum_d = 0;
    std::vector<std::pair<i64, i64>> deck;
    for (int v = 0; v < n; ++v) {
      auto [card, rz] = restrict_tensor(adj, z, v);
      deck.push_back(packet(card, rz));
      sum_c += deck.back().first; sum_d += deck.back().second;
    }
    need(sum_c == (n - 3) * cval, "independent C deck");
    need(sum_d == (n - 4) * dval, "independent D deck");
    Fraction parent = n % 2 ? Fraction((n - 3) * cval, 4 * dval)
                            : Fraction((n - 3) * cval, 2 * dval);
    Fraction mean = (n - 1) % 2 ? Fraction((n - 4) * sum_c, 4 * sum_d)
                                : Fraction((n - 4) * sum_c, 2 * sum_d);
    Fraction transported = Fraction(n % 2 ? 1 : 2, n % 2 ? 2 : 1) * mean;
    need(parent == transported, "independent parity transport");
    ++rows;
  }
  return rows;
}

static void named_hostiles() {
  Adj prime11 = decode(3169369058263173ULL, 11);
  Tensor c11 = response_capacity(prime11);
  auto p11 = packet(prime11, c11);
  need(hamilton_paths(prime11) == 23685, "independent prime11 H");
  need(p11 == std::pair<i64, i64>(4220068008LL, 88725253576LL), "independent prime11 packet");
  auto [bad10, keep, where] = delete_vertex(prime11, 10);
  Tensor c10 = response_capacity(bad10);
  auto p10 = packet(bad10, c10);
  need(hamilton_paths(bad10) == 2037, "independent card H");
  need(p10 == std::pair<i64, i64>(755197384LL, 1016215996LL), "independent card packet");

  Adj order12 = {3070, 3644, 3704, 3824, 4064, 4032,
                 3970, 3846, 3598, 1024, 2049, 512};
  Tensor c12 = response_capacity(order12);
  auto p12 = packet(order12, c12);
  need(hamilton_paths(order12) == 27759, "independent order12 H");
  need(p12 == std::pair<i64, i64>(-94387092144LL, 323484198928LL), "independent order12 packet");
  i64 sum_d = 0;
  Fraction maximum;
  std::vector<std::tuple<i64, i64, Fraction>> deck;
  for (int v = 0; v < 12; ++v) {
    auto [card, rz] = restrict_tensor(order12, c12, v);
    auto [cv, dv] = packet(card, rz);
    Fraction tau(8 * cv, 4 * dv);
    deck.push_back({cv, dv, tau});
    sum_d += dv;
    Fraction absolute(tau.num < 0 ? -tau.num : tau.num, tau.den);
    if (maximum < absolute) maximum = absolute;
  }
  i64 sum_c = 0;
  for (auto [cv, dv, tau] : deck) sum_c += cv;
  Fraction mean(8 * sum_c, 4 * sum_d);
  need(mean == Fraction(-53092739331LL, 80871049732LL), "independent order12 mean");
  need(maximum == Fraction(9073595176LL, 12026131621LL), "independent order12 maximum");
  need(Fraction(2, 1) * mean == Fraction(-53092739331LL, 40435524866LL),
       "independent order12 amplification");
  std::cout << "named prime11=" << p11.first << "," << p11.second
            << " bad10=" << p10.first << "," << p10.second
            << " order12_mean=" << mean.str()
            << " order12_max=" << maximum.str() << '\n';
}

int main() {
  std::cout << "marked_rows " << marked_audit_through_five() << '\n';
  auto [gates, zeros, minimum] = monotonicity_n6();
  std::cout << "n6_monotonicity " << gates << " " << zeros << " " << minimum << '\n';
  std::cout << "restriction_rows " << restriction_algebra() << '\n';
  named_hostiles();
  std::cout << "THM4167_INDEPENDENT_AUDIT_PASS\n";
  return 0;
}
