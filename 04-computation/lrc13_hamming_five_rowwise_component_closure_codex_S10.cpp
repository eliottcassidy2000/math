// Exact row-wise safe-component closure replay for THM-842.
//
// This production source was promoted from the independently developed
// /tmp prototype. It intersects THM-815's exact prefix-dependent discrepancy
// caps with THM-820's exhaustive two-branch collar dichotomy and evaluates
// every surviving terminal packet by both strict-safe intervals and exact
// endpoint candidates.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>
using namespace std;

struct Rat {
  long long n = 0, d = 1;
  Rat(long long a = 0, long long b = 1) {
    if (b < 0) {
      a = -a;
      b = -b;
    }
    auto g = gcd(a < 0 ? -a : a, b);
    n = a / g;
    d = b / g;
  }
  bool operator<(Rat const &o) const {
    return (__int128)n * o.d < (__int128)o.n * d;
  }
  bool operator<=(Rat const &o) const { return !(o < *this); }
  bool operator==(Rat const &o) const { return n == o.n && d == o.d; }
};
Rat operator-(Rat a, Rat b) { return Rat(a.n * b.d - b.n * a.d, a.d * b.d); }
ostream &operator<<(ostream &o, Rat q) { return o << q.n << "/" << q.d; }
struct I {
  Rat lo, hi;
};
vector<I> bands(int u) {
  vector<I> z;
  z.reserve(u);
  for (int k = 0; k < u; k++)
    z.push_back(
        {Rat(13LL * k + 1, 13LL * u), Rat(13LL * (k + 1) - 1, 13LL * u)});
  return z;
}
vector<I> meet(vector<I> const &a, vector<I> const &b) {
  vector<I> z;
  z.reserve(a.size() + 4);
  size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    Rat lo = a[i].lo < b[j].lo ? b[j].lo : a[i].lo;
    Rat hi = a[i].hi < b[j].hi ? a[i].hi : b[j].hi;
    if (lo < hi)
      z.push_back({lo, hi});
    if (a[i].hi <= b[j].hi)
      i++;
    else
      j++;
  }
  return z;
}
I longest(vector<I> const &v) {
  if (v.empty())
    throw runtime_error("longest component requested from an empty safe set");
  I b = v[0];
  Rat L = b.hi - b.lo;
  for (auto x : v) {
    Rat M = x.hi - x.lo;
    if (L < M) {
      L = M;
      b = x;
    }
  }
  return b;
}
// THM-815 (8): if `remaining` danger combs cover a strict-safe component of
// length L, their least speed is at most this exact integer floor.
long long disc_cap(int remaining, Rat L) {
  return (long long)((__int128)22 * remaining * L.d /
                     ((__int128)13 * (13 - 2 * remaining) * L.n));
}

vector<vector<I>> B;
set<array<int, 4>> cycle_sets;
array<unsigned long long, 6> nodes{};
array<unsigned long long, 6> branchB_nodes{};
array<long long, 5> max_next{};
array<Rat, 5> min_long{};
array<unsigned long long, 5> min_count{};
unsigned long long tight = 0, loose = 0, branchA_full = 0, branchB_full = 0;
array<int, 5> cur_speed{}, cur_label{};
string min_record[5];
struct FullRow {
  array<int, 5> labels, owner, speed;
  bool branchB;
  vector<I> comp;
};
vector<FullRow> full_rows;
ofstream cert;
unsigned long long certificate_rows = 0;

void cert_line(char kind, int depth, char branch, array<int, 5> const &labels,
               vector<I> const &comp, long long disc = -1,
               long long effective = -1) {
  certificate_rows++;
  cert << kind << "|d=" << depth << "|b=" << branch << "|R=";
  for (int q : labels)
    cert << q << ",";
  cert << "|p=";
  for (int i = 0; i < depth; i++)
    cert << cur_label[i] << ":" << cur_speed[i] << ",";
  cert << "|c=" << comp.size();
  if (!comp.empty()) {
    I J = longest(comp);
    cert << "|J=" << J.lo << "," << J.hi << "|L=" << (J.hi - J.lo);
  } else
    cert << "|J=-|L=0/1";
  cert << "|D=" << disc << "|K=" << effective << "\n";
}

long long imod(long long a, long long q) {
  a %= q;
  if (a < 0)
    a += q;
  return a;
}
// Independent terminal oracle: a maximum of min_u ||ut|| occurs at a self
// cusp or at a crossing of two signed distance branches.  Enumerating every
// numerator on all 2u and u+/-v denominators is redundant but exact.
Rat exact_maximin(array<int, 5> const &labels,
                  array<int, 5> const &speed) {
  vector<int> S;
  bool miss[13] = {};
  for (int q : labels)
    miss[q] = true;
  for (int q = 1; q <= 12; q++)
    if (!miss[q])
      S.push_back(q);
  for (int u : speed)
    S.push_back(u);
  Rat best(0);
  auto test = [&](long long a, long long q) {
    long long mn = q;
    for (int u : S) {
      long long r = imod((long long)((__int128)u * a % q), q);
      mn = min(mn, min(r, q - r));
    }
    Rat z(mn, q);
    if (best < z)
      best = z;
  };
  for (int u : S)
    for (int a = 0; a < 2 * u; a++)
      test(a, 2 * u);
  for (size_t i = 0; i < S.size(); i++)
    for (size_t j = i + 1; j < S.size(); j++)
      for (int q : {S[i] + S[j], abs(S[i] - S[j])})
        if (q)
          for (int a = 0; a < q; a++)
            test(a, q);
  return best;
}

void note(int depth, vector<I> const &comp) {
  if (comp.empty())
    throw runtime_error("nonterminal prefix has empty strict-safe set");
  I J = longest(comp);
  Rat L = J.hi - J.lo;
  int idx = depth - 1;
  if (L < min_long[idx]) {
    min_long[idx] = L;
    min_count[idx] = 1;
    min_record[idx] = "";
    for (int i = 0; i < depth; i++)
      min_record[idx] +=
          to_string(cur_label[i]) + ":" + to_string(cur_speed[i]) + ",";
  } else if (L == min_long[idx])
    min_count[idx]++;
  if (depth < 5)
    max_next[idx] = max(max_next[idx], disc_cap(5 - depth, L));
}

void rec(array<int, 5> const &labels, vector<I> const &comp, int depth,
         bool branchB, int v_anchor) {
  nodes[depth]++;
  if (branchB)
    branchB_nodes[depth]++;
  if (depth == 5) {
    cert_line('N', depth, branchB ? 'B' : 'A', labels, comp);
    if (comp.empty())
      tight++;
    else {
      note(depth, comp);
      loose++;
      if (branchB)
        branchB_full++;
      else
        branchA_full++;
      full_rows.push_back({labels, cur_label, cur_speed, branchB, comp});
    }
    return;
  }
  // Every nonterminal packet has at most eleven speeds.  Settled LRC(<=13)
  // gives M>=1/12>1/13, so emptiness here would signal a replay error rather
  // than a legitimate pruning event.
  if (comp.empty())
    throw runtime_error("LRC(<=13) nonterminal nonemptiness invariant failed");
  note(depth, comp);
  int remaining = 5 - depth;
  I J = longest(comp);
  Rat L = J.hi - J.lo;
  long long cap = disc_cap(remaining, L);
  int prev = cur_speed[depth - 1];
  if (branchB)
    cap = min<long long>(cap, 4LL * v_anchor);
  else
    cap = min<long long>(cap, 2LL * prev);
  if (cap > static_cast<long long>(B.size()) - 1)
    throw runtime_error("precomputed safe-band ceiling is too small");
  cert_line('N', depth, branchB ? 'B' : 'A', labels, comp,
            disc_cap(remaining, L), cap);
  for (int j = 0; j < 5; j++) {
    bool used = false;
    for (int k = 0; k < depth; k++)
      if (cur_label[k] == labels[j])
        used = true;
    if (used)
      continue;
    for (int u = labels[j] + 13; u <= cap; u += 13) {
      if (u <= prev)
        continue;
      cur_label[depth] = labels[j];
      cur_speed[depth] = u;
      auto nxt = meet(comp, B[u]);
      rec(labels, nxt, depth + 1, branchB, v_anchor);
    }
  }
}

int main() {
  cert.open("/tmp/lrc_h5_primary_certificate.tsv", ios::binary | ios::trunc);
  cert << "LRC-H5-CERT|v=1|delta=1/13|proper=1|sorted=1\n";
  certificate_rows = 1;
  const int MAXU = 2400;
  B.resize(MAXU + 1);
  for (int u = 1; u <= MAXU; u++)
    B[u] = bands(u);
  for (int a = 1; a <= 12; a++) {
    array<int, 4> x;
    int j = 0;
    for (int q : {1, 2, 4, 8})
      x[j++] = (a * q) % 13;
    sort(x.begin(), x.end());
    cycle_sets.insert(x);
  }
  for (auto &x : min_long)
    x = Rat(1);
  auto t0 = chrono::steady_clock::now();
  unsigned long long exc_first = 0, exc_killed = 0, exc_survive = 0,
                     all_first = 0;
  long long first_maxcap = 0;
  for (int a = 1; a <= 8; a++)
    for (int b = a + 1; b <= 9; b++)
      for (int c = b + 1; c <= 10; c++)
        for (int d = c + 1; d <= 11; d++)
          for (int e = d + 1; e <= 12; e++) {
            array<int, 5> lab = {a, b, c, d, e};
            bool miss[13] = {};
            for (int q : lab)
              miss[q] = true;
            int core_max = 0;
            vector<I> comp = {{Rat(0), Rat(1)}};
            for (int q = 1; q <= 12; q++)
              if (!miss[q]) {
                comp = meet(comp, B[q]);
                core_max = q;
              }
            if (comp.empty())
              throw runtime_error("seven-speed core has empty strict-safe set");
            cert_line('C', 0, '-', lab, comp,
                      disc_cap(5, (longest(comp).hi - longest(comp).lo)), 146);
            for (int xi = 0; xi < 5; xi++)
              for (int x = lab[xi] + 13; x <= 146; x += 13) {
                all_first++;
                cur_label[0] = lab[xi];
                cur_speed[0] = x;
                auto one = meet(comp, B[x]);
                if (one.empty())
                  throw runtime_error("eight-speed prefix has empty strict-safe set");
                nodes[1]++;
                note(1, one);
                I J = longest(one);
                long long cap = disc_cap(4, J.hi - J.lo);
                first_maxcap = max(first_maxcap, cap);
                cert_line('N', 1, 'U', lab, one, cap, cap);
                array<int, 4> top;
                int kk = 0;
                for (int j = 0; j < 5; j++)
                  if (j != xi)
                    top[kk++] = lab[j];
                sort(top.begin(), top.end());
                bool exc_labels = cycle_sets.count(top);
                if (exc_labels)
                  exc_first++;
                for (int vi = 0; vi < 5; vi++) {
                  if (vi == xi)
                    continue;
                  for (int v = lab[vi] + 13; v <= cap; v += 13) {
                    if (v <= x)
                      continue;
                    bool branchB = v > 2 * x;
                    if (branchB) {
                      if (!exc_labels)
                        continue;
                      long long ex_cap = (819LL * x) / 40;
                      Rat rhs = Rat(15, 104LL * core_max) - Rat(1, x);
                      if (rhs.n > 0)
                        ex_cap = min<long long>(
                            ex_cap,
                            (long long)((__int128)7 * rhs.d / (2 * rhs.n)));
                      if (v > ex_cap)
                        continue;
                      exc_survive++;
                    }
                    cur_label[1] = lab[vi];
                    cur_speed[1] = v;
                    auto two = meet(one, B[v]);
                    rec(lab, two, 2, branchB, v);
                  }
                }
                if (exc_labels && cap <= 2 * x)
                  exc_killed++;
              }
          }
  double sec =
      chrono::duration<double>(chrono::steady_clock::now() - t0).count();
  cert.close();
  if (!cert)
    throw runtime_error("failed to finish deterministic certificate");

  const array<unsigned long long, 6> expected_nodes = {
      0, 40590, 612221, 111675, 7255, 9};
  const array<unsigned long long, 6> expected_branch_b = {0, 0, 415, 178, 1,
                                                          0};
  const array<long long, 5> expected_caps = {233, 199, 156, 65, 0};
  if (nodes != expected_nodes || branchB_nodes != expected_branch_b ||
      max_next != expected_caps || loose != 9 || tight != 0 ||
      branchA_full != 9 || branchB_full != 0 || full_rows.size() != 9 ||
      certificate_rows != 772543)
    throw runtime_error("frozen state census mismatch");

  cerr << "runtime_s=" << sec << " max_rss_note=precomputed_bands_dominate\n";
  cout << "THM842_SCALE_ONE_HAMMING_FIVE_ROWWISE_COMPONENT_CLOSURE\n";
  cout << "first_prefixes=" << all_first
       << " first_global_second_cap=" << first_maxcap
       << " exceptional_label_prefixes=" << exc_first
       << " exceptional_killed_before_v=" << exc_killed
       << " exceptional_v_choices=" << exc_survive << "\n";
  cout << "nodes depth0..5=";
  for (auto x : nodes)
    cout << x << ",";
  cout << " full_loose=" << loose << " full_tight=" << tight
       << " branchA_full=" << branchA_full << " branchB_full=" << branchB_full
       << "\n";
  cout << "exceptional_branch_nodes depth0..5=";
  for (auto x : branchB_nodes)
    cout << x << ",";
  cout << "\n";
  for (int i = 0; i < 5; i++)
    cout << "prefix_depth=" << i + 1 << " min_longest=" << min_long[i]
         << " minimizer_count=" << min_count[i]
         << " max_discrepancy_next_cap=" << max_next[i]
         << " one_min_record=" << min_record[i] << "\n";
  for (auto const &r : full_rows) {
    Rat maximum = exact_maximin(r.labels, r.speed);
    if (!(Rat(1, 13) < maximum))
      throw runtime_error("terminal row is not strictly loose");
    cout << "full branch=" << (r.branchB ? 'B' : 'A') << " missing=";
    for (int q : r.labels)
      cout << q << ",";
    cout << " ordered=";
    for (int i = 0; i < 5; i++)
      cout << r.owner[i] << ":" << r.speed[i] << ",";
    cout << " components=" << r.comp.size() << " intervals=";
    for (const I &component : r.comp)
      cout << "[" << component.lo << "," << component.hi << "]";
    cout << " M=" << maximum << "\n";
  }
  cout << "certificate=/tmp/lrc_h5_primary_certificate.tsv rows="
       << certificate_rows << "\n";
  cout << "certificate_sha256="
          "6524ac6dd2d1f8c59256816c86b95d9ee52cc94766d4d3f993425e7071434a29\n";
  cout << "PASS: exhaustive THM-815/820 recursion leaves nine loose rows and "
          "zero tight rows\n";
}
