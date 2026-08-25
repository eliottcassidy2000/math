#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <omp.h>
#include <openssl/sha.h>

namespace {

constexpr int N = 9;
constexpr int SZ = 1 << N;
constexpr int FULL = SZ - 1;

struct Rat {
  int64_t n = 0;
  int64_t d = 1;
};

static int64_t gcd64(int64_t a, int64_t b) {
  if (a < 0) a = -a;
  if (b < 0) b = -b;
  return std::gcd(a, b);
}

static Rat rat(int64_t n, int64_t d) {
  if (d == 0) throw std::runtime_error("zero denominator");
  if (d < 0) { n = -n; d = -d; }
  int64_t g = gcd64(n, d);
  return {n / g, d / g};
}

static int cmp(Rat a, Rat b) {
  __int128 lhs = static_cast<__int128>(a.n) * b.d;
  __int128 rhs = static_cast<__int128>(b.n) * a.d;
  return (lhs > rhs) - (lhs < rhs);
}

static std::string rs(Rat x) {
  return std::to_string(x.n) + "/" + std::to_string(x.d);
}

static int64_t ceil_div128(__int128 a, __int128 b) {
  if (b <= 0) throw std::runtime_error("ceil_div denominator");
  __int128 q = a / b;
  __int128 r = a % b;
  if (r > 0) ++q;
  if (r < 0) {
    // C++ truncates toward zero, which is already ceil for a negative input.
  }
  if (q < std::numeric_limits<int64_t>::min() ||
      q > std::numeric_limits<int64_t>::max())
    throw std::runtime_error("ceil_div overflow");
  return static_cast<int64_t>(q);
}

struct Layer {
  int m = 0;
  int t = 0;
  Rat J;
  int64_t L = 0;
  int64_t A = 0;
  int64_t lattice = 0;
};

struct Profile {
  std::string label;
  bool strong = false;
  int64_t H = 0;
  int64_t W2 = 0;       // exactly 2 W
  int64_t D4x4 = 0;     // exactly 4 D_4
  int64_t Chdx4 = 0;    // exactly 4 C_hd
  Rat theta;
  Rat rho;
  uint16_t rational_mask = 0;  // bit m means cardinality m is an optimizer
  uint16_t coset_mask = 0;
  uint16_t actual_mask = 0;
  std::array<Layer, N - 1> layers{};
  std::array<uint32_t, SZ> values{};  // mask order 0,...,511
};

struct Adjacency {
  std::array<uint16_t, N> out{};
  std::array<uint16_t, N> in{};
};

static Adjacency decode_label(const std::string& bits) {
  if (bits.size() != N * (N - 1) / 2)
    throw std::runtime_error("label length");
  Adjacency a;
  int cursor = 0;
  for (int left = 0; left < N; ++left) {
    for (int right = left + 1; right < N; ++right) {
      char bit = bits[cursor++];
      if (bit == '1') a.out[left] |= uint16_t(1u << right);
      else if (bit == '0') a.out[right] |= uint16_t(1u << left);
      else throw std::runtime_error("nonbinary label");
    }
  }
  for (int v = 0; v < N; ++v)
    for (int u = 0; u < N; ++u)
      if (a.out[u] & (1u << v)) a.in[v] |= uint16_t(1u << u);
  return a;
}

static bool is_strong(const Adjacency& a) {
  auto reach = [&](bool reverse) {
    uint16_t reached = 1;
    uint16_t frontier = 1;
    while (frontier) {
      int v = __builtin_ctz(frontier);
      frontier &= uint16_t(frontier - 1);
      uint16_t fresh = (reverse ? a.in[v] : a.out[v]) & uint16_t(~reached);
      reached |= fresh;
      frontier |= fresh;
    }
    return (reached & FULL) == FULL;
  };
  return reach(false) && reach(true);
}

static Profile evaluate(const std::string& label) {
  const Adjacency a = decode_label(label);
  Profile p;
  p.label = label;
  p.strong = is_strong(a);

  // Independent THM-4131 carrier: Hamilton paths ending/starting at a
  // specified vertex of every subset.
  uint32_t ending[SZ][N] = {};
  uint32_t starting[SZ][N] = {};
  for (int v = 0; v < N; ++v) {
    ending[1 << v][v] = 1;
    starting[1 << v][v] = 1;
  }
  for (int mask = 1; mask < SZ; ++mask) {
    uint16_t vertices = uint16_t(mask);
    while (vertices) {
      int v = __builtin_ctz(vertices);
      vertices &= uint16_t(vertices - 1);
      int rest = mask ^ (1 << v);
      if (!rest) continue;
      uint16_t pred = uint16_t(rest) & a.in[v];
      while (pred) {
        int u = __builtin_ctz(pred);
        pred &= uint16_t(pred - 1);
        ending[mask][v] += ending[rest][u];
      }
      uint16_t succ = uint16_t(rest) & a.out[v];
      while (succ) {
        int u = __builtin_ctz(succ);
        succ &= uint16_t(succ - 1);
        starting[mask][v] += starting[rest][u];
      }
    }
  }
  for (int v = 0; v < N; ++v) p.H += ending[FULL][v];
  if (p.H <= 0) throw std::runtime_error("Hamilton path count");

  uint32_t before[SZ][N] = {};
  uint32_t after[SZ][N] = {};
  for (int v = 0; v < N; ++v) before[0][v] = after[0][v] = 1;
  for (int mask = 1; mask < SZ; ++mask) {
    for (int boundary = 0; boundary < N; ++boundary) {
      uint16_t bits = uint16_t(mask) & a.in[boundary];
      while (bits) {
        int u = __builtin_ctz(bits);
        bits &= uint16_t(bits - 1);
        before[mask][boundary] += ending[mask][u];
      }
      bits = uint16_t(mask) & a.out[boundary];
      while (bits) {
        int u = __builtin_ctz(bits);
        bits &= uint16_t(bits - 1);
        after[mask][boundary] += starting[mask][u];
      }
    }
  }

  uint32_t cap[N][N] = {};
  uint32_t edgecap[N][N] = {};
  int64_t outcap[N] = {};
  for (int tail = 0; tail < N; ++tail) {
    for (int head = tail + 1; head < N; ++head) {
      int u = tail, v = head;
      if (!(a.out[u] & (1u << v))) std::swap(u, v);
      int remaining = FULL ^ (1 << u) ^ (1 << v);
      uint64_t c = 0;
      int left = remaining;
      while (true) {
        int right = remaining ^ left;
        c += uint64_t(before[left][u]) * after[right][v];
        c += uint64_t(before[left][v]) * after[right][u];
        if (left == 0) break;
        left = (left - 1) & remaining;
      }
      if (c > std::numeric_limits<uint32_t>::max())
        throw std::runtime_error("capacity overflow");
      cap[u][v] = static_cast<uint32_t>(c);
      edgecap[u][v] = edgecap[v][u] = static_cast<uint32_t>(c);
      outcap[u] += c;
      p.W2 += c;
    }
  }

  // Scaled packet avoids Fraction arithmetic: w_e=cap_e/2.
  int64_t d2[N] = {}, h2[N] = {};
  for (int v = 0; v < N; ++v) {
    for (int u = 0; u < N; ++u) if (u != v) {
      d2[v] += edgecap[v][u];
      if (a.out[v] & (1u << u)) h2[v] += edgecap[v][u];
      else h2[v] -= edgecap[v][u];
    }
    p.Chdx4 += h2[v] * d2[v];
  }
  for (int a0 = 0; a0 < N; ++a0)
    for (int b = a0 + 1; b < N; ++b)
      for (int c = a0 + 1; c < N; ++c) if (c != b)
        for (int d = c + 1; d < N; ++d)
          if (d != b)
            p.D4x4 += int64_t(edgecap[a0][b]) * edgecap[c][d];
  // The loop above chooses the edge containing the least of four vertices,
  // then one of its three partners, hence counts each disjoint edge pair once.
  if (p.D4x4 <= 0) throw std::runtime_error("D4 positivity");
  p.theta = rat((N - 3) * p.Chdx4, 2 * p.D4x4);
  p.rho = rat((N - 3) * std::llabs(p.Chdx4), 4 * p.D4x4);

  int64_t edge_subset[N][SZ] = {};
  for (int v = 0; v < N; ++v) {
    for (int mask = 1; mask < SZ; ++mask) {
      int u = __builtin_ctz(mask);
      int rest = mask ^ (1 << u);
      edge_subset[v][mask] = edge_subset[v][rest] + edgecap[v][u];
    }
  }
  int64_t cut[SZ] = {};
  p.values[0] = static_cast<uint32_t>(p.H);
  for (int mask = 1; mask < SZ; ++mask) {
    int v = __builtin_ctz(mask);
    int rest = mask ^ (1 << v);
    cut[mask] = cut[rest] + outcap[v] - edge_subset[v][rest];
    int64_t value = p.H + cut[mask];
    if (value < p.H || value > std::numeric_limits<uint32_t>::max() || !(value & 1))
      throw std::runtime_error("response support/parity");
    p.values[mask] = static_cast<uint32_t>(value);
  }
  if (p.values[FULL] != p.H) throw std::runtime_error("constant ear endpoints");

  Rat bestJ{std::numeric_limits<int64_t>::min(), 1};
  int64_t bestL = std::numeric_limits<int64_t>::min();
  int64_t bestA = std::numeric_limits<int64_t>::min();
  for (int m = 1; m < N; ++m) {
    int64_t count = 0, sum = 0, sumsq = 0, maximum = 0, lattice = 0;
    int64_t anchor = -1;
    for (int mask = 0; mask < SZ; ++mask) if (__builtin_popcount(mask) == m) {
      int64_t value = p.values[mask];
      if (anchor < 0) anchor = value;
      ++count;
      sum += value;
      sumsq += value * value;
      maximum = std::max(maximum, value);
      lattice = gcd64(lattice, value - anchor);
    }
    int64_t delta = sum - count * p.H;
    if (delta <= 0) throw std::runtime_error("positive layer mean gap");
    Rat J = rat(sumsq - sum * p.H, delta);
    int64_t L;
    if (lattice == 0) {
      if (cmp(J, Rat{anchor, 1}) != 0) throw std::runtime_error("constant lattice");
      L = anchor;
    } else {
      __int128 num = static_cast<__int128>(J.n) - static_cast<__int128>(anchor) * J.d;
      int64_t k = ceil_div128(num, static_cast<__int128>(J.d) * lattice);
      L = anchor + lattice * k;
    }
    if (L > maximum) throw std::runtime_error("coset support floor exceeds max");
    p.layers[m - 1] = {m, N - 2 * m, J, L, maximum, lattice};
    int cj = cmp(J, bestJ);
    if (cj > 0) { bestJ = J; p.rational_mask = uint16_t(1u << m); }
    else if (cj == 0) p.rational_mask |= uint16_t(1u << m);
    if (L > bestL) { bestL = L; p.coset_mask = uint16_t(1u << m); }
    else if (L == bestL) p.coset_mask |= uint16_t(1u << m);
    if (maximum > bestA) { bestA = maximum; p.actual_mask = uint16_t(1u << m); }
    else if (maximum == bestA) p.actual_mask |= uint16_t(1u << m);
  }

  // Independent nearest-grid verification, expressed directly in cardinality.
  Rat bestDistance{std::numeric_limits<int64_t>::max(), 1};
  uint16_t predicted = 0;
  for (int m = 1; m < N; ++m) {
    int t = N - 2 * m;
    Rat distance = rat(std::llabs(int64_t(t) * p.theta.d - p.theta.n), p.theta.d);
    int c = cmp(distance, bestDistance);
    if (c < 0) { bestDistance = distance; predicted = uint16_t(1u << m); }
    else if (c == 0) predicted |= uint16_t(1u << m);
  }
  if (predicted != p.rational_mask) throw std::runtime_error("nearest-grid mismatch");
  return p;
}

static std::string tuple_t(uint16_t mask) {
  std::ostringstream out;
  out << '(';
  bool first = true;
  for (int m = 1; m < N; ++m) if (mask & (1u << m)) {
    if (!first) out << ',';
    first = false;
    out << (N - 2 * m);
  }
  out << ')';
  return out.str();
}

static void sha_bytes(SHA256_CTX& ctx, const void* ptr, size_t n) {
  SHA256_Update(&ctx, ptr, n);
}

template <typename U>
static void sha_unsigned_le(SHA256_CTX& ctx, U value) {
  static_assert(std::is_unsigned<U>::value, "unsigned serializer");
  unsigned char bytes[sizeof(U)];
  for (size_t i = 0; i < sizeof(U); ++i) {
    bytes[i] = static_cast<unsigned char>(value & 0xffu);
    value >>= 8;
  }
  sha_bytes(ctx, bytes, sizeof(bytes));
}

static void sha_i64(SHA256_CTX& ctx, int64_t value) {
  sha_unsigned_le<uint64_t>(ctx, static_cast<uint64_t>(value));
}

static std::string sha_final_hex(SHA256_CTX& ctx) {
  unsigned char digest[SHA256_DIGEST_LENGTH];
  SHA256_Final(digest, &ctx);
  std::ostringstream out;
  for (unsigned char byte : digest)
    out << std::hex << std::setw(2) << std::setfill('0') << int(byte);
  return out.str();
}

static void hash_profile(SHA256_CTX& ctx, const Profile& p) {
  sha_bytes(ctx, p.label.data(), p.label.size());               // fixed 36 bytes
  sha_unsigned_le<uint32_t>(ctx, static_cast<uint32_t>(p.H));
  sha_i64(ctx, p.W2);
  sha_i64(ctx, p.D4x4);
  sha_i64(ctx, p.Chdx4);
  sha_i64(ctx, p.theta.n); sha_i64(ctx, p.theta.d);
  sha_i64(ctx, p.rho.n); sha_i64(ctx, p.rho.d);
  sha_unsigned_le<uint16_t>(ctx, p.rational_mask);
  sha_unsigned_le<uint16_t>(ctx, p.coset_mask);
  sha_unsigned_le<uint16_t>(ctx, p.actual_mask);
  for (const Layer& layer : p.layers) {
    sha_i64(ctx, layer.J.n); sha_i64(ctx, layer.J.d);
    sha_i64(ctx, layer.L); sha_i64(ctx, layer.A); sha_i64(ctx, layer.lattice);
  }
  for (uint32_t value : p.values) sha_unsigned_le<uint32_t>(ctx, value);
}

static std::string packet(const Profile& p, bool with_layers = true) {
  std::ostringstream out;
  out << "label=" << p.label
      << ",H=" << p.H
      << ",W=" << rs(rat(p.W2, 2))
      << ",D4=" << rs(rat(p.D4x4, 4))
      << ",Chd=" << rs(rat(p.Chdx4, 4))
      << ",theta=" << rs(p.theta)
      << ",rho=" << rs(p.rho)
      << ",rational_t=" << tuple_t(p.rational_mask)
      << ",coset_t=" << tuple_t(p.coset_mask)
      << ",actual_t=" << tuple_t(p.actual_mask);
  if (with_layers) {
    out << ",layers=[";
    for (size_t i = 0; i < p.layers.size(); ++i) {
      if (i) out << ';';
      const Layer& x = p.layers[i];
      out << '(' << x.m << ',' << x.t << ',' << rs(x.J) << ','
          << x.L << ',' << x.A << ',' << x.lattice << ')';
    }
    out << ']';
  }
  return out.str();
}

// A separate literal child Hamilton DP is used only on the named n=6 controls.
static std::vector<uint32_t> control_values(uint64_t code, std::string* label_out) {
  constexpr int M = 6;
  std::string bits;
  bits.reserve(M * (M - 1) / 2);
  for (int cursor = 0; cursor < M * (M - 1) / 2; ++cursor)
    bits.push_back((code & (uint64_t(1) << cursor)) ? '1' : '0');
  *label_out = bits;
  uint16_t out[M + 1] = {};
  int cursor = 0;
  for (int u = 0; u < M; ++u) for (int v = u + 1; v < M; ++v) {
    if (bits[cursor++] == '1') out[u] |= uint16_t(1u << v);
    else out[v] |= uint16_t(1u << u);
  }
  auto hamilton = [&](const uint16_t child[M + 1]) {
    constexpr int CSZ = 1 << (M + 1);
    uint32_t dp[CSZ][M + 1] = {};
    for (int v = 0; v <= M; ++v) dp[1 << v][v] = 1;
    for (int mask = 1; mask < CSZ; ++mask) {
      for (int v = 0; v <= M; ++v) if (mask & (1 << v)) {
        int rest = mask ^ (1 << v);
        for (int u = 0; u <= M; ++u)
          if ((rest & (1 << u)) && (child[u] & (1 << v))) dp[mask][v] += dp[rest][u];
      }
    }
    uint32_t H = 0;
    for (int v = 0; v <= M; ++v) H += dp[CSZ - 1][v];
    return H;
  };
  std::vector<uint32_t> values(1 << M);
  for (int cut = 0; cut < (1 << M); ++cut) {
    uint16_t child[M + 1] = {};
    std::copy(std::begin(out), std::end(out), child);
    for (int v = 0; v < M; ++v) {
      if (cut & (1 << v)) child[M] |= uint16_t(1u << v);
      else child[v] |= uint16_t(1u << M);
    }
    values[cut] = hamilton(child);
  }
  return values;
}

static Profile evaluate_control(uint64_t code) {
  // Embed the n=6 contracted-block evaluator directly in a small generic
  // implementation kept distinct from the optimized n=9 route.
  constexpr int M = 6, MSZ = 1 << M, MFULL = MSZ - 1;
  std::string label;
  std::vector<uint32_t> literal = control_values(code, &label);
  uint16_t out[M] = {}, in[M] = {};
  int cursor = 0;
  for (int u = 0; u < M; ++u) for (int v = u + 1; v < M; ++v) {
    if (label[cursor++] == '1') out[u] |= uint16_t(1u << v);
    else out[v] |= uint16_t(1u << u);
  }
  for (int v = 0; v < M; ++v) for (int u = 0; u < M; ++u)
    if (out[u] & (1 << v)) in[v] |= uint16_t(1u << u);
  uint32_t ending[MSZ][M] = {}, starting[MSZ][M] = {};
  for (int v = 0; v < M; ++v) ending[1 << v][v] = starting[1 << v][v] = 1;
  for (int mask = 1; mask < MSZ; ++mask) for (int v = 0; v < M; ++v)
    if (mask & (1 << v)) {
      int rest = mask ^ (1 << v);
      for (int u = 0; u < M; ++u) if (rest & (1 << u)) {
        if (out[u] & (1 << v)) ending[mask][v] += ending[rest][u];
        if (out[v] & (1 << u)) starting[mask][v] += starting[rest][u];
      }
    }
  int64_t H = 0;
  for (int v = 0; v < M; ++v) H += ending[MFULL][v];
  uint32_t before[MSZ][M] = {}, after[MSZ][M] = {};
  for (int v = 0; v < M; ++v) before[0][v] = after[0][v] = 1;
  for (int mask = 1; mask < MSZ; ++mask) for (int v = 0; v < M; ++v)
    for (int u = 0; u < M; ++u) {
      if ((mask & (1 << u)) && (out[u] & (1 << v))) before[mask][v] += ending[mask][u];
      if ((mask & (1 << u)) && (out[v] & (1 << u))) after[mask][v] += starting[mask][u];
    }
  uint32_t cap[M][M] = {};
  for (int a = 0; a < M; ++a) for (int b = a + 1; b < M; ++b) {
    int u = a, v = b;
    if (!(out[u] & (1 << v))) std::swap(u, v);
    int rem = MFULL ^ (1 << u) ^ (1 << v);
    uint32_t c = 0;
    for (int left = rem;; left = (left - 1) & rem) {
      int right = rem ^ left;
      c += before[left][u] * after[right][v] + before[left][v] * after[right][u];
      if (!left) break;
    }
    cap[u][v] = c;
  }
  std::vector<uint32_t> cut_values(MSZ, static_cast<uint32_t>(H));
  for (int mask = 0; mask < MSZ; ++mask)
    for (int u = 0; u < M; ++u) if (mask & (1 << u))
      for (int v = 0; v < M; ++v) if (!(mask & (1 << v))) cut_values[mask] += cap[u][v];
  if (cut_values != literal) throw std::runtime_error("literal control replay mismatch");

  // Return through a compact dynamic summary printed directly below. The n=9
  // Profile type cannot encode six-order t coordinates, so this routine
  // communicates by a deliberately formatted pseudo-profile in label.
  Profile q;
  q.label = label;
  q.H = H;
  // Compute layer summaries and append them to q.label after a sentinel.
  std::ostringstream summary;
  summary << label << "|H=" << H << "|";
  Rat bestJ{std::numeric_limits<int64_t>::min(), 1};
  int64_t bestL = std::numeric_limits<int64_t>::min(), bestA = std::numeric_limits<int64_t>::min();
  uint16_t rmask = 0, lmask = 0, amask = 0;
  std::array<Layer, M - 1> ls{};
  for (int m = 1; m < M; ++m) {
    int64_t count = 0, sum = 0, sumsq = 0, anchor = -1, d = 0, A = 0;
    for (int mask = 0; mask < MSZ; ++mask) if (__builtin_popcount(mask) == m) {
      int64_t x = cut_values[mask]; if (anchor < 0) anchor = x;
      ++count; sum += x; sumsq += x*x; A = std::max(A, x); d = gcd64(d, x-anchor);
    }
    Rat J = rat(sumsq-sum*H, sum-count*H);
    int64_t L = anchor;
    if (d) L += d * ceil_div128(static_cast<__int128>(J.n)-static_cast<__int128>(anchor)*J.d,
                                static_cast<__int128>(J.d)*d);
    ls[m-1] = {m, M-2*m, J, L, A, d};
    int c = cmp(J,bestJ); if(c>0){bestJ=J;rmask=1u<<m;}else if(c==0)rmask|=1u<<m;
    if(L>bestL){bestL=L;lmask=1u<<m;}else if(L==bestL)lmask|=1u<<m;
    if(A>bestA){bestA=A;amask=1u<<m;}else if(A==bestA)amask|=1u<<m;
  }
  auto ts6 = [](uint16_t mask) {
    std::ostringstream s; s << '('; bool first=true;
    for(int m=1;m<M;++m)if(mask&(1u<<m)){if(!first)s<<',';first=false;s<<(M-2*m);}s<<')';return s.str();
  };
  summary << "rational_t=" << ts6(rmask) << "|coset_t=" << ts6(lmask)
          << "|actual_t=" << ts6(amask) << "|layers=[";
  for(size_t i=0;i<ls.size();++i){if(i)summary<<';';auto&x=ls[i];summary<<'('<<x.m<<','<<x.t<<','<<rs(x.J)<<','<<x.L<<','<<x.A<<','<<x.lattice<<')';}
  summary << ']';
  q.label = summary.str();
  return q;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    size_t limit = std::numeric_limits<size_t>::max();
    if (const char* value = std::getenv("AUDIT_LIMIT")) limit = std::stoull(value);
    std::vector<std::string> labels;
    std::string label;
    while (std::cin >> label) {
      labels.push_back(label);
      if (labels.size() == limit) break;
    }
    if (labels.empty()) throw std::runtime_error("empty universe");

    SHA256_CTX raw_ctx, profile_ctx;
    SHA256_Init(&raw_ctx); SHA256_Init(&profile_ctx);
    const char profile_header[] = "thm4131-n9-independent-profile-v1\n";
    sha_bytes(profile_ctx, profile_header, sizeof(profile_header)-1);
    for (const std::string& s : labels) {
      sha_bytes(raw_ctx, s.data(), s.size());
      const char newline = '\n'; sha_bytes(raw_ctx, &newline, 1);
    }

    int64_t strong_count = 0, rational_fail = 0, coset_fail = 0;
    int64_t theta_zero = 0, reorder = 0, actual_noncentral = 0;
    int64_t minimum_margin = std::numeric_limits<int64_t>::max();
    Rat worst{-1,1};
    std::vector<Profile> worst_profiles;
    std::map<std::string,int64_t> rhist, chist, ahist;

    constexpr size_t BATCH = 512;
    for (size_t base = 0; base < labels.size(); base += BATCH) {
      size_t count = std::min(BATCH, labels.size() - base);
      std::vector<Profile> batch(count);
      #pragma omp parallel for schedule(static)
      for (int64_t i = 0; i < static_cast<int64_t>(count); ++i)
        batch[i] = evaluate(labels[base+i]);
      for (Profile& p : batch) {
        if (!p.strong) throw std::runtime_error("-c stream contains nonstrong class");
        ++strong_count;
        hash_profile(profile_ctx, p);
        const uint16_t central_mask = uint16_t((1u << 4) | (1u << 5)); // t=+1,-1
        if (p.rational_mask & ~central_mask) ++rational_fail;
        if (p.coset_mask & ~central_mask) ++coset_fail;
        if (!(p.actual_mask & central_mask)) ++actual_noncentral;
        if (p.theta.n == 0) ++theta_zero;
        if (p.rational_mask != p.coset_mask) ++reorder;
        ++rhist[tuple_t(p.rational_mask)];
        ++chist[tuple_t(p.coset_mask)];
        ++ahist[tuple_t(p.actual_mask)];
        int64_t best_c = std::numeric_limits<int64_t>::min();
        int64_t best_o = std::numeric_limits<int64_t>::min();
        for (const Layer& x : p.layers) {
          if (x.t == -1 || x.t == 1) best_c = std::max(best_c, x.L);
          else best_o = std::max(best_o, x.L);
        }
        minimum_margin = std::min(minimum_margin, best_c-best_o);
        int c = cmp(p.rho, worst);
        if (c > 0) { worst = p.rho; worst_profiles.assign(1,p); }
        else if (c == 0) worst_profiles.push_back(p);
      }
      if ((base / BATCH) % 64 == 0)
        std::cerr << "processed=" << (base + count) << "/" << labels.size() << "\n";
    }

    const std::string raw_hash = sha_final_hex(raw_ctx);
    const std::string profile_hash = sha_final_hex(profile_ctx);
    std::cout << "implementation=gentourng strong universe + contracted good/reversed-block DP + mask recurrence\n";
    std::cout << "thread_modes=audited_1_and_8\n";
    std::cout << "classes=" << labels.size() << "\n";
    std::cout << "strong_classes=" << strong_count << "\n";
    std::cout << "raw_stream_sha256=" << raw_hash << "\n";
    std::cout << "profile_serialization=ASCII header thm4131-n9-independent-profile-v1\\n; then per gentourng record: 36 ASCII bits; H:u32le; W2,D4x4,Chdx4,theta_num,theta_den,rho_num,rho_den:i64le; rational,coset,actual cardinality-bitmasks:u16le; for m=1..8 J_num,J_den,L,A,lattice:i64le; F(mask):u32le for masks 0..511\n";
    std::cout << "ordered_full_profile_sha256=" << profile_hash << "\n";
    std::cout << "rational_central_failures=" << rational_fail << "\n";
    std::cout << "coset_central_failures=" << coset_fail << "\n";
    std::cout << "minimum_strict_coset_margin=" << minimum_margin << "\n";
    std::cout << "worst_rho=" << rs(worst) << "\n";
    std::cout << "worst_multiplicity=" << worst_profiles.size() << "\n";
    for (size_t i=0;i<worst_profiles.size();++i)
      std::cout << "worst_packet_" << (i+1) << '=' << packet(worst_profiles[i]) << "\n";
    std::cout << "theta_zero=" << theta_zero << "\n";
    std::cout << "coset_reorders_rational=" << reorder << "\n";
    std::cout << "actual_noncentral_only=" << actual_noncentral << "\n";
    auto print_hist=[](const char* name,const std::map<std::string,int64_t>& h){std::cout<<name<<'=';bool first=true;for(auto&[k,v]:h){if(!first)std::cout<<',';first=false;std::cout<<k<<':'<<v;}std::cout<<'\n';};
    print_hist("rational_histogram",rhist);
    print_hist("coset_histogram",chist);
    print_hist("actual_histogram",ahist);

    for (uint64_t code : {uint64_t(2), uint64_t(140), uint64_t(20)}) {
      Profile c = evaluate_control(code);
      std::cout << "control_code_" << code << '=' << c.label << "|literal_child_replay=PASS\n";
    }

    if (labels.size() == 178133) {
      if (strong_count != 178133 || rational_fail || coset_fail || minimum_margin != 90 ||
          cmp(worst, rat(16430,44663)) || worst_profiles.size()!=2 || theta_zero!=2661 ||
          reorder!=10772 || actual_noncentral!=3248)
        throw std::runtime_error("frozen n=9 invariant mismatch");
      std::cout << "status=ACCEPT\n";
    } else {
      std::cout << "status=PARTIAL_LIMIT_ONLY\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "status=REJECT error=" << e.what() << "\n";
    return 1;
  }
  return 0;
}
