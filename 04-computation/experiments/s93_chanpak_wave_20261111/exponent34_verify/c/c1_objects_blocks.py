#!/usr/bin/env python3
"""Independent check (verifier C) of the Objects section and of the block law / three-block formula of THM-5141 (C).

Exact integer arithmetic only.  Nothing is imported from the coordinator's scripts.

Checks
  O1  g(v) = v + 1 - 2 q(v) - (v mod 2), x(h) = g(2h)                         all v < 2^21
  O2  0 <= q(h) <= (h-1)/2 for h >= 1;  h + 2 <= x(h) <= 2h + 1 for h >= 1; x(0) = 1
  O3  g(word) = number of Hamiltonian paths of the one-apex tournament T(word) (built from the definition:
      chain 1 -> 2 -> ... transitive, apex x -> i iff bit b_i = 1, b_1 = most significant bit), by permutations
      (chain length <= 6) and by subset DP (chain length <= 11)
  O4  all tournaments on N <= 6 labelled vertices: the H-values of those with a transitive one-vertex deletion are
      exactly {g(v) : v < 2^(N-1)}
  O5  {g(v) : v < 2^m} = {x(h) : h < 2^(m-1)} (m <= 20) ; g(odd) reduces to g(even)
  B1  unrolled formula q(h) = sum_i h_i ceil(floor(h/2^(i+1))/2)
  B2  block law 2 q(A 2^m + B) = 2 q(A) + 2 q(B) + A rev_m(B) + (A mod 2) B_(m-1), exhaustive small + random large
  B3  three-block formula with explicit K(P, Z), all (t', w, s) small, s = 0 and w = 1 included; Z_(s-1) = c mod 2;
      Z -> c bijection onto (2^s, 2^(s+1)]
  G1  max |gamma| over ALL odd pairs and parities, w <= 10 (bound 2^w)
"""
import itertools
import random
import sys

sys.setrecursionlimit(10000)
LIM = 1 << 21

# ---------------------------------------------------------------- recursions, as quoted in the theorem file
q = [0] * LIM
for h in range(1, LIM):
    q[h] = q[h >> 1] + ((((h >> 1) + 1) >> 1) if (h & 1) else 0)   # q(2h+1) = q(h) + ceil(h/2), q(2h) = q(h)

g = [0] * LIM
g[0] = 1
for v in range(1, LIM):
    u = v >> 1
    g[v] = g[u] if (v & 1) else g[u] + u + (u & 1)                 # g(2u+1) = g(u), g(2u) = g(u) + u + (u mod 2)


def xval(h):
    return 2 * h + 1 - 2 * q[h]


def rev(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r


def qbig(h):
    """q for arbitrary size integers (unrolled recursion, exact)."""
    t = 0
    while h:
        if h & 1:
            t += ((h >> 1) + 1) >> 1
        h >>= 1
    return t


def out(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- O1, O2
bad = 0
for v in range(LIM):
    if g[v] != v + 1 - 2 * q[v] - (v & 1):
        bad += 1
out("O1a g(v) = v + 1 - 2q(v) - (v mod 2) for all v < 2^21: violations", bad)
bad = 0
for h in range(LIM // 2):
    if xval(h) != g[2 * h]:
        bad += 1
out("O1b x(h) = g(2h) for all h < 2^20: violations", bad)
bad = 0
tight_low = tight_high = 0
for h in range(1, LIM):
    if not (0 <= 2 * q[h] <= h - 1):
        bad += 1
    xv = xval(h)
    if not (h + 2 <= xv <= 2 * h + 1):
        bad += 1
    tight_low += (xv == h + 2)
    tight_high += (xv == 2 * h + 1)
out("O2 0 <= q(h) <= (h-1)/2 and h+2 <= x(h) <= 2h+1 for 1 <= h < 2^21: violations", bad,
    "| x(0) =", xval(0), "(h = 0 is excluded: x(0) = 1 < 0 + 2)",
    "| equality x = h+2 at", tight_low, "values of h, x = 2h+1 at", tight_high)
out("   equality cases x(h) = h + 2 (h < 64):", [h for h in range(1, 64) if xval(h) == h + 2])
out("   equality cases x(h) = 2h + 1 (h < 64):", [h for h in range(1, 64) if xval(h) == 2 * h + 1])

# ---------------------------------------------------------------- O3 one-apex tournaments from the definition


def one_apex_rows(bits):
    """vertex 0 = apex, vertices 1..n = chain; rows[i] = bitmask of out-neighbours.
    chain: i -> j for i < j; apex -> i iff bits[i-1] == 1, else i -> apex."""
    n = len(bits)
    rows = [0] * (n + 1)
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            rows[i] |= 1 << j
        if bits[i - 1]:
            rows[0] |= 1 << i
        else:
            rows[i] |= 1
    return rows


def ham_perm(rows):
    N = len(rows)
    c = 0
    for p in itertools.permutations(range(N)):
        ok = True
        for i in range(N - 1):
            if not (rows[p[i]] >> p[i + 1]) & 1:
                ok = False
                break
        c += ok
    return c


def ham_dp(rows):
    N = len(rows)
    full = (1 << N) - 1
    dp = [[0] * N for _ in range(full + 1)]
    for i in range(N):
        dp[1 << i][i] = 1
    for mask in range(1, full + 1):
        row = dp[mask]
        for i in range(N):
            ways = row[i]
            if ways:
                avail = rows[i] & ~mask
                while avail:
                    b = avail & -avail
                    avail ^= b
                    dp[mask | b][b.bit_length() - 1] += ways
    return sum(dp[full])


cnt_perm = cnt_dp = bad = 0
for n in range(0, 7):
    for bits in itertools.product((0, 1), repeat=n):
        v = int("".join(map(str, bits)), 2) if n else 0
        Hp = ham_perm(one_apex_rows(bits))
        cnt_perm += 1
        if Hp != g[v]:
            bad += 1
out("O3a g(word) = H(T(word)) by literal permutations, all words of chain length <= 6:", cnt_perm, "words, violations", bad)
bad = 0
for n in range(0, 10):
    for bits in itertools.product((0, 1), repeat=n):
        v = int("".join(map(str, bits)), 2) if n else 0
        Hd = ham_dp(one_apex_rows(bits))
        cnt_dp += 1
        if Hd != g[v]:
            bad += 1
out("O3b g(word) = H(T(word)) by subset DP (pure Python), all words of chain length <= 9:", cnt_dp, "words, violations", bad)

try:
    import numpy as np
    from numba import njit

    @njit(cache=False)
    def ham_dp_nb(rows, N):
        full = (1 << N) - 1
        dp = np.zeros((full + 1, N), dtype=np.int64)
        for i in range(N):
            dp[1 << i, i] = 1
        for mask in range(1, full + 1):
            for i in range(N):
                ways = dp[mask, i]
                if ways != 0:
                    avail = rows[i] & (full ^ mask)
                    for j in range(N):
                        if (avail >> j) & 1:
                            dp[mask | (1 << j), j] += ways
        tot = 0
        for i in range(N):
            tot += dp[full, i]
        return tot

    @njit(cache=False)
    def check_len(n, garr):
        """all words of chain length n; returns number of violations."""
        N = n + 1
        rows = np.zeros(N, dtype=np.int64)
        badc = 0
        for v in range(1 << n):
            for i in range(N):
                rows[i] = 0
            for i in range(1, n + 1):
                for j in range(i + 1, n + 1):
                    rows[i] |= 1 << j
                b = (v >> (n - i)) & 1              # b_1 = most significant bit
                if b == 1:
                    rows[0] |= 1 << i
                else:
                    rows[i] |= 1
            if ham_dp_nb(rows, N) != garr[v]:
                badc += 1
        return badc

    garr = np.array(g[:1 << 15], dtype=np.int64)
    tot = 0
    badn = 0
    for n in range(0, 15):
        badn += check_len(n, garr)
        tot += 1 << n
    out("O3b' g(word) = H(T(word)) by subset DP (numba), all words of chain length <= 14:", tot, "words, violations", badn)
except Exception as e:  # pragma: no cover
    out("O3b' numba check skipped:", repr(e))

# the opposite conventions (to see that the convention matters for the function but not for the value set)
bad_comp = bad_revc = 0
for n in range(1, 9):
    for bits in itertools.product((0, 1), repeat=n):
        v = int("".join(map(str, bits)), 2)
        comp = tuple(1 - b for b in bits)
        vc = int("".join(map(str, comp)), 2)
        revc = tuple(1 - b for b in reversed(bits))
        vr = int("".join(map(str, revc)), 2)
        if g[vc] != g[v]:
            bad_comp += 1
        if g[vr] != g[v]:
            bad_revc += 1
out("O3c (information) words of length <= 8 with g(complement) != g: ", bad_comp,
    "; with g(reverse-complement) != g:", bad_revc, "(converse invariance => 0)")

# ---------------------------------------------------------------- O4 all tournaments of order N <= 6


def all_tournament_values(N):
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    vals_one_apex = set()
    n_one_apex = 0
    for code in range(1 << len(pairs)):
        rows = [0] * N
        for k, (i, j) in enumerate(pairs):
            if (code >> k) & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        # transitive one-vertex deletion?  T - a is transitive iff its score sequence is 0,1,...,N-2
        ok = False
        for a in range(N):
            keep = [i for i in range(N) if i != a]
            km = sum(1 << i for i in keep)
            sc = sorted(bin(rows[i] & km).count("1") for i in keep)
            if sc == list(range(N - 1)):
                ok = True
                break
        if ok:
            n_one_apex += 1
            vals_one_apex.add(ham_dp(rows))
    return n_one_apex, vals_one_apex


for N in range(1, 7):
    n1, vals = all_tournament_values(N)
    gv = set(g[v] for v in range(1 << (N - 1)))
    out("O4 order", N, ": labelled tournaments with a transitive one-vertex deletion:", n1,
        "| H-values", sorted(vals), "| equal to {g(v): v < 2^(N-1)}:", vals == gv)

# ---------------------------------------------------------------- O5 value sets
for m in (4, 8, 12, 16, 20):
    Sg = set(g[v] for v in range(1 << m))
    Sx = set(xval(h) for h in range(1 << (m - 1)))
    out("O5 m =", m, ": {g(v): v < 2^m} == {x(h): h < 2^(m-1)}:", Sg == Sx, "| size", len(Sg))

# ---------------------------------------------------------------- B1 unrolled formula
bad = 0
for h in range(1 << 16):
    t = 0
    i = 0
    while (h >> i):
        if (h >> i) & 1:
            t += ((h >> (i + 1)) + 1) // 2          # ceil(floor(h/2^(i+1))/2)
        i += 1
    if t != q[h]:
        bad += 1
out("B1 q(h) = sum_i h_i ceil(floor(h/2^(i+1))/2), h < 2^16: violations", bad)

# ---------------------------------------------------------------- B2 block law
bad = cnt = 0
for A in range(0, 1 << 8):
    for m in range(1, 10):
        for B in range(1 << m):
            lhs = 2 * q[(A << m) + B]
            rhs = 2 * q[A] + 2 * q[B] + A * rev(B, m) + (A & 1) * ((B >> (m - 1)) & 1)
            cnt += 1
            if lhs != rhs:
                bad += 1
out("B2a block law, A < 2^8 (A = 0 and odd A included), 1 <= m <= 9, all B < 2^m:", cnt, "cases, violations", bad)
random.seed(20261111)
bad = cnt = 0
for _ in range(20000):
    a = random.randrange(0, 200)
    m = random.randrange(1, 200)
    A = random.getrandbits(a) if a else 0
    B = random.getrandbits(m)
    lhs = 2 * qbig((A << m) + B)
    rhs = 2 * qbig(A) + 2 * qbig(B) + A * rev(B, m) + (A & 1) * ((B >> (m - 1)) & 1)
    cnt += 1
    if lhs != rhs:
        bad += 1
out("B2b block law, random A (< 2^200), m < 200, B < 2^m:", cnt, "cases, violations", bad)
# m = 0 is excluded in the statement; what would happen at m = 0?  (B = 0, rev_0 = 0, 'B_(-1)' undefined)
out("B2c (information) the law is stated for m >= 1; at m = 0 it would read 2q(A) = 2q(A) + (A mod 2) B_(-1): the last term"
    " must be dropped, which is why s = 0 is treated separately (c = 2, c mod 2 = 0).")

# ---------------------------------------------------------------- B3 three blocks
bad = cnt = 0
badK = 0
badpar = 0
badrange = 0
for tp in range(0, 6):                       # tp = 0 means P = 0 (not needed by the theorem, formula still holds)
    Ps = [0] if tp == 0 else range(1 << (tp - 1), 1 << tp)
    for w in range(1, 7):
        for s in range(0, 6):
            cs = set()
            for Z in range(1 << s):
                c = (1 << (s + 1)) - rev(Z, s)
                cs.add(c)
                if not ((1 << s) < c <= (1 << (s + 1))):
                    badrange += 1
                if s >= 1 and ((Z >> (s - 1)) & 1) != (c & 1):
                    badpar += 1
                if s == 0 and c != 2:
                    badpar += 1
                for P in Ps:
                    Kexp = 2 * P * (1 << (w + s)) + 2 * Z + 1 - 2 * q[P] - 2 * q[Z] - P * (1 << w) * rev(Z, s)
                    K0 = None
                    for x in range(1 << w):
                        h = (P << (w + s)) + (x << s) + Z
                        val = xval(h)
                        rest = c * x - P * rev(x, w) - 2 * q[x] - (P & 1) * ((x >> (w - 1)) & 1) - (c & 1) * (x & 1)
                        K = val - rest
                        cnt += 1
                        if K0 is None:
                            K0 = K
                        if K != K0:
                            bad += 1
                        if K != Kexp:
                            badK += 1
            if cs != set(range((1 << s) + 1, (1 << (s + 1)) + 1)):
                badrange += 1
out("B3 three-block formula x(h) = K(P,Z) + c x - P rev_w(x) - 2q(x) - (P mod 2) x_(w-1) - (c mod 2) x_0,",
    "t' <= 5 (and P = 0), 1 <= w <= 6, 0 <= s <= 5:", cnt, "cases; K depends on x:", bad,
    "; K != 2P 2^(w+s) + 2Z + 1 - 2q(P) - 2q(Z) - P 2^w rev_s(Z):", badK,
    "; Z_(s-1) != c mod 2 (s >= 1) or c != 2 (s = 0):", badpar, "; Z -> c not a bijection onto (2^s, 2^(s+1)]:", badrange)

# random large three-block instances (exact big integers)
bad = cnt = 0
for _ in range(20000):
    tp = random.randrange(1, 80)
    w = random.randrange(1, 80)
    s = random.randrange(0, 80)
    P = (1 << (tp - 1)) | random.getrandbits(tp - 1) if tp > 1 else 1
    Z = random.getrandbits(s) if s else 0
    c = (1 << (s + 1)) - rev(Z, s)
    x1 = random.getrandbits(w)
    x2 = random.getrandbits(w)
    h1 = (P << (w + s)) + (x1 << s) + Z
    h2 = (P << (w + s)) + (x2 << s) + Z
    lhs = (2 * h1 + 1 - 2 * qbig(h1)) - (2 * h2 + 1 - 2 * qbig(h2))
    X = x1 - x2
    Xi = rev(x1, w) - rev(x2, w)
    gam = 2 * (qbig(x1) - qbig(x2)) + (P & 1) * (((x1 >> (w - 1)) & 1) - ((x2 >> (w - 1)) & 1)) \
        + (c & 1) * ((x1 & 1) - (x2 & 1))
    cnt += 1
    if lhs != c * X - P * Xi - gam:
        bad += 1
out("B3b difference form x(h) - x(h') = c X - P Xi - gamma for random blocks up to 80 bits each (any x, x'):",
    cnt, "cases, violations", bad)

# ---------------------------------------------------------------- G1 the bound |gamma| <= 2^w over all odd pairs
for w in range(1, 11):
    mx = 0
    arg = None
    n_odd = 0
    top = 1 << (w - 1)
    for x1 in range(1 << w):
        for x2 in range(1 << w):
            if ((x1 ^ x2) & 1) == 0 or ((x1 ^ x2) & top) == 0:
                continue
            n_odd += 1
            base = 2 * (q[x1] - q[x2])
            dt = ((x1 >> (w - 1)) & 1) - ((x2 >> (w - 1)) & 1)
            db = (x1 & 1) - (x2 & 1)
            for pp in (0, 1):
                for cp in (0, 1):
                    gam = base + pp * dt + cp * db
                    if abs(gam) > mx:
                        mx = abs(gam)
                        arg = (x1, x2, pp, cp)
    out("G1 w =", w, ": odd pairs", n_odd, "(4^(w-1) =", 4 ** (w - 1), "; 2 for w = 1)", "max |gamma| =", mx,
        "2^w =", 1 << w, "attained at (x, x', P mod 2, c mod 2) =", arg)
out("done")
