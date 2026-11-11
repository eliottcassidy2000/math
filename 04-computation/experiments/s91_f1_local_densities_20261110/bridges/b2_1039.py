#!/usr/bin/env python3
"""S91 bridges (2): the claims about 1039 (understand_sweep_numbers.txt items (f)7-(f)8), verified with own code.

Conventions (CONTEXT_S91.md section 1).
  psi(0) = 0, psi(2h) = psi(h) + h, psi(2h+1) = psi(h) + floor(h/2) + 1;  x(h) = 2 psi(h) + 1,  y(h) = 2h + 1.
  One-apex tournament R(h): u = 2h written in binary with m = bitlength(u) digits b_1 ... b_m (b_1 the leading digit);
  a transitive chain c_1 -> c_2 -> ... -> c_m (c_i -> c_j for i < j) and an apex a with a -> c_i iff b_i = 1.
  F1 join J(h1, h2, s): R(h1), a transitive chain Z of s vertices, and the converse R(h2)^op, with every arc directed
  R(h1) -> Z -> R(h2)^op and R(h1) -> R(h2)^op, except the single arc b -> a between the two apexes.
  Claimed count: N = x1 x2 + 2^s y1 y2 - [s = 0];  r(N) = number of ORDERED triples (h1, h2, s) with value N.
  Tournament masks: masks[i] has bit j set iff i -> j (out-neighbour masks).

Checks (all exact integer arithmetic):
  1. x(h) = H(R(h)) and y(h) = number of ordered covers (A1, A2) of R(h), A1 a possibly empty path, A2 a path that
     starts at the apex, for all h < 256 (brute-force DP on the explicit tournaments);
  2. H(J(h1, h2, s)) = x1 x2 + 2^s y1 y2 - [s = 0] for all h1, h2 < 12 and s <= 3 (576 explicit joins);
  3. psi(h) >= h/2 and psi(h) <= h for h < 2^20 (the induction is in the docstring of b3), hence psi(h) = t forces
     t <= h <= 2t; psi never takes the values 3, 6, 519 (so 7, 13, 1039 are not one-apex counts);
  4. r(N) by a complete search (y1 y2 2^s <= N + 1 and x >= (y + 1)/2 bound the search): r = 0 at 1039 and the
     quoted values r(3), r(5), r(7), r(11), r(13), r(19), r(21), r(27); the HYP-9657 list of zeros among N = 4q - 1;
     the number of misses among odd N in [1025, 2049);
  5. arithmetic of 1039: prime, 1039 + 1 = 2^4 * 5 * 13, residues, ord_1039(2) = 519 (so <2> = squares), 3 is a
     primitive root;
  6. least order: the maximum of H over ALL tournaments on n <= 8 vertices by exhaustive enumeration (every
     tournament has a Hamiltonian path, so it is isomorphic to one containing 0 -> 1 -> ... -> n-1: 2^C(n-1,2)
     completions); max H(8) = 661 < 1039, so order >= 9;
  7. the order-9 realiser [496, 317, 473, 449, 392, 28, 434, 290, 32]: it is a tournament and has exactly 1039
     Hamiltonian paths (bitmask DP, and an independent count over all 9! permutations).
"""
import itertools
import sys

import numpy as np
from numba import njit, prange

RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append(bool(cond))
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""), flush=True)


def psi(h):
    if h == 0:
        return 0
    if h % 2 == 0:
        return psi(h // 2) + h // 2
    k = h // 2
    return psi(k) + k // 2 + 1


def xy(h):
    return 2 * psi(h) + 1, 2 * h + 1


@njit(cache=False)
def ham_paths(adj, n):
    """number of Hamiltonian paths of the digraph with 0/1 adjacency matrix adj (n <= 20)."""
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for i in range(n):
        dp[1 << i, i] = 1
    for mask in range(1, full + 1):
        for last in range(n):
            c = dp[mask, last]
            if c == 0:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1 == 0 and adj[last, nxt]:
                    dp[mask | (1 << nxt), nxt] += c
    tot = 0
    for last in range(n):
        tot += dp[full, last]
    return tot


@njit(cache=False)
def covers_second_from(adj, n, root):
    """ordered covers (A1, A2): A1 a possibly empty path, A2 a path starting at root, A1 u A2 = everything."""
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)          # dp[mask][last]: paths with vertex set mask ending at last
    for i in range(n):
        dp[1 << i, i] = 1
    paths = np.zeros(1 << n, np.int64)            # number of paths with vertex set mask
    from_root = np.zeros(1 << n, np.int64)        # ... starting at root
    dpr = np.zeros((1 << n, n), np.int64)
    dpr[1 << root, root] = 1
    for mask in range(1, full + 1):
        for last in range(n):
            c = dp[mask, last]
            cr = dpr[mask, last]
            if c == 0 and cr == 0:
                continue
            paths[mask] += c
            from_root[mask] += cr
            for nxt in range(n):
                if (mask >> nxt) & 1 == 0 and adj[last, nxt]:
                    dp[mask | (1 << nxt), nxt] += c
                    dpr[mask | (1 << nxt), nxt] += cr
    tot = 0
    for mask in range(1, full + 1):
        if (mask >> root) & 1:
            rest = full ^ mask
            tot += from_root[mask] * (1 if rest == 0 else paths[rest])
    return tot


def one_apex(h):
    """adjacency matrix of R(h); vertices 0..m-1 = chain c_1..c_m, vertex m = apex."""
    u = 2 * h
    m = u.bit_length()
    adj = np.zeros((m + 1, m + 1), np.uint8)
    for i in range(m):
        for j in range(i + 1, m):
            adj[i, j] = 1
        bit = (u >> (m - 1 - i)) & 1
        if bit:
            adj[m, i] = 1
        else:
            adj[i, m] = 1
    return adj, m


def f1_join(h1, h2, s):
    a1, m1 = one_apex(h1)
    a2, m2 = one_apex(h2)
    a2 = a2.T.copy()                                # converse
    n1, n2 = m1 + 1, m2 + 1
    n = n1 + s + n2
    adj = np.zeros((n, n), np.uint8)
    adj[:n1, :n1] = a1
    adj[n1 + s:, n1 + s:] = a2
    for i in range(n1):
        for j in range(n1, n):
            adj[i, j] = 1
    for i in range(n1, n1 + s):
        for j in range(i + 1, n):
            adj[i, j] = 1
    apex1, apex2 = m1, n1 + s + m2
    adj[apex1, apex2] = 0
    adj[apex2, apex1] = 1
    return adj, n


# ---------------------------------------------------------------- 1, 2: conventions
ok_x = ok_y = True
for h in range(256):
    adj, m = one_apex(h)
    x, y = xy(h)
    ok_x &= (ham_paths(adj, m + 1) == x)
    ok_y &= (covers_second_from(adj, m + 1, m) == y)
check("1a x(h) = 2 psi(h) + 1 = H(R(h)) for all h < 256 (explicit tournaments, DP)", ok_x)
check("1b y(h) = 2h + 1 = number of ordered covers (A1 possibly empty, A2 starting at the apex), h < 256", ok_y)
ok_f1 = True
njoin = 0
for h1 in range(12):
    for h2 in range(12):
        for s in range(4):
            adj, n = f1_join(h1, h2, s)
            x1, y1 = xy(h1)
            x2, y2 = xy(h2)
            njoin += 1
            ok_f1 &= (ham_paths(adj, n) == x1 * x2 + (y1 * y2 << s) - (1 if s == 0 else 0))
check("2  H(J(h1, h2, s)) = x1 x2 + 2^s y1 y2 - [s = 0] on %d explicit joins (h1, h2 < 12, s <= 3)" % njoin, ok_f1)

# ---------------------------------------------------------------- 3: psi bounds and missing values
LIM = 1 << 20
ps = np.zeros(LIM, np.int64)
for k in range(20):
    lo, hi = 1 << k, 1 << (k + 1)
    idx = np.arange(lo, hi, dtype=np.int64)
    half = idx >> 1
    ps[lo:hi] = ps[half] + np.where(idx & 1 == 1, (half >> 1) + 1, half)
hh = np.arange(LIM, dtype=np.int64)
check("3a 2 psi(h) >= h and psi(h) <= h for all h < 2^20", bool(np.all(2 * ps >= hh) and np.all(ps <= hh)))
check("3b the vectorised table agrees with the recursion for h < 4096", all(int(ps[h]) == psi(h) for h in range(4096)))
missing = {}
for t in (3, 6, 519):
    missing[t] = [h for h in range(t, 2 * t + 1) if psi(h) == t]
check("3c psi never equals 3, 6 or 519 (search h in [t, 2t], sufficient by 3a): 7, 13, 1039 are not one-apex counts",
      all(not v for v in missing.values()), str(missing))
words = {N: [h for h in range((N - 1) // 2, N) if 2 * psi(h) + 1 == N] for N in (1, 3, 5, 7, 13, 19, 21, 27, 1039)}
check("3d one-apex words with x(h) = N: 3: {1}, 5: {2,3}, 19: {10,13}, 27: {17,23}; none for 7, 13, 21, 1039",
      words[3] == [1] and words[5] == [2, 3] and words[19] == [10, 13] and words[27] == [17, 23]
      and not words[7] and not words[13] and not words[21] and not words[1039], str(words))


# ---------------------------------------------------------------- 4: r(N)
def r_table(X):
    """r(N) for all N <= X, with witnesses; complete because x >= (y + 1)/2 >= 1 and y1 y2 2^s <= N + 1."""
    H = [h for h in range(0, X // 2 + 2) if 2 * h + 1 <= X + 1]
    st = [(2 * int(ps[h]) + 1, 2 * h + 1, h) for h in H]
    r = [0] * (X + 1)
    wit = {}
    for x1, y1, h1 in st:
        for x2, y2, h2 in st:
            yy = y1 * y2
            if yy > X + 1:
                break                                 # y increases with h
            base = x1 * x2
            s = 0
            while True:
                N = base + (yy << s) - (1 if s == 0 else 0)
                if N > X:
                    break
                r[N] += 1
                wit.setdefault(N, []).append((h1, h2, s))
                s += 1
    return r, wit


X = 4096
r, wit = r_table(X)
quoted = {3: 1, 5: 3, 7: 0, 11: 2, 13: 0, 19: 4, 21: 0, 27: 5, 1039: 0}
got = {N: r[N] for N in quoted}
check("4a r(N) at 3, 5, 7, 11, 13, 19, 21, 27, 1039 equals 1, 3, 0, 2, 0, 4, 0, 5, 0", got == quoted, str(got))
for N in (3, 5, 19, 27):
    print("      witnesses (h1, h2, s) of N = %d: %s" % (N, wit[N]))


def r_single(N):
    """independent per-N search: for every (h1, s) solve for h2 by scanning."""
    cnt = 0
    for s in range(0, N.bit_length() + 1):
        target = N + (1 if s == 0 else 0)
        h1 = 0
        while (2 * h1 + 1) << s <= target:
            x1, y1 = 2 * int(ps[h1]) + 1, 2 * h1 + 1
            h2 = 0
            while x1 * (h2 + 1) + ((y1 * (2 * h2 + 1)) << s) <= target + x1 * (h2 + 1):
                x2, y2 = 2 * int(ps[h2]) + 1, 2 * h2 + 1
                if ((y1 * y2) << s) > target:
                    break
                if x1 * x2 + ((y1 * y2) << s) == target:
                    cnt += 1
                h2 += 1
            h1 += 1
    return cnt


check("4b independent per-N search gives the same r at the nine values and r(1039) = 0",
      all(r_single(N) == v for N, v in quoted.items()))


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


rab = [N for N in range(7, 3000, 4) if is_prime(N) and is_prime((N + 1) // 4)]
rab0 = [N for N in rab if r[N] == 0]
check("4c HYP-9657 (B): 30 integers N = 4q - 1 < 3000 with N, q prime; r = 0 exactly at 7, 43, 67, 163, 211, 547, 1171, 1723",
      len(rab) == 30 and rab0 == [7, 43, 67, 163, 211, 547, 1171, 1723], "%d, %s" % (len(rab), rab0))
window = [N for N in range(1025, 2049, 2)]
miss = [N for N in window if r[N] == 0]
mean_r = sum(r[N] for N in window) / len(window)
check("4d misses among the 512 odd N in [1025, 2049): 65 (so r(1039) = 0 is not exceptional)", len(miss) == 65,
      "%d misses, mean r = %.3f, 1039 in misses: %s" % (len(miss), mean_r, 1039 in miss))
near = [(N, r[N]) for N in range(1029, 1051, 2)]
print("      r near 1039:", near)
print("      is (1039 + 1)/4 = 260 prime? no: 1039 is not of the form 4q - 1 with q prime (260 = 2^2 * 5 * 13)")


# ---------------------------------------------------------------- 5: arithmetic of 1039
def order(a, p):
    k, x = 1, a % p
    while x != 1:
        x = x * a % p
        k += 1
    return k


check("5a 1039 is prime, 1040 = 2^4 * 5 * 13 = 2^4 (2^6 + 1)", is_prime(1039) and 1040 == 16 * 5 * 13 == 16 * (2 ** 6 + 1))
check("5b residues: 1039 = 1 mod 3, -1 mod 5, -1 mod 13, 3 mod 7, 13 mod 19, 13 mod 27, 7 mod 8, 7 mod 12",
      (1039 % 3, 1039 % 5, 1039 % 13, 1039 % 7, 1039 % 19, 1039 % 27, 1039 % 8, 1039 % 12) == (1, 4, 12, 3, 13, 13, 7, 7))
check("5c ord_1039(2) = 519 = (1039 - 1)/2, so <2> is the group of squares (as at p = 7)", order(2, 1039) == 519)
check("5d 3 is a primitive root modulo 1039; 519 = 3 * 173", order(3, 1039) == 1038 and 519 == 3 * 173 and is_prime(173))
p7 = [q for q in range(3, 1100) if is_prime(q) and q % 8 == 7 and order(2, q) == (q - 1) // 2]
print("      primes q < 1100 with q = 7 mod 8 and <2> = squares:", p7)
check("5e 1039 is one of %d primes below 1100 with <2> = squares and q = 7 mod 8 (not a rare property)" % len(p7), 1039 in p7 and 7 in p7)


# ---------------------------------------------------------------- 6: maximum of H at n <= 8
@njit(parallel=True, cache=False)
def max_h_fixed_path(n, nthreads):
    """all tournaments on n vertices containing the path 0 -> 1 -> ... -> n-1; returns (max H, count, occurring values)."""
    pairs_i = np.zeros(64, np.int64)
    pairs_j = np.zeros(64, np.int64)
    k = 0
    for i in range(n):
        for j in range(i + 2, n):
            pairs_i[k] = i
            pairs_j[k] = j
            k += 1
    total = 1 << k
    full = (1 << n) - 1
    best = np.zeros(nthreads, np.int64)
    bestcnt = np.zeros(nthreads, np.int64)
    seen = np.zeros((nthreads, 1024), np.uint8)
    chunk = total // nthreads + 1
    for t in prange(nthreads):
        out = np.zeros(n, np.int64)
        dp = np.zeros((1 << n, n), np.int64)
        for code in range(t * chunk, min(total, (t + 1) * chunk)):
            for i in range(n):
                out[i] = 0
            for i in range(n - 1):
                out[i] |= 1 << (i + 1)
            for b in range(k):
                if (code >> b) & 1:
                    out[pairs_i[b]] |= 1 << pairs_j[b]
                else:
                    out[pairs_j[b]] |= 1 << pairs_i[b]
            for mask in range(1, full + 1):
                for last in range(n):
                    dp[mask, last] = 0
            for i in range(n):
                dp[1 << i, i] = 1
            for mask in range(1, full + 1):
                for last in range(n):
                    c = dp[mask, last]
                    if c == 0:
                        continue
                    avail = out[last] & (full ^ mask)
                    while avail:
                        low = avail & (-avail)
                        nxt = 0
                        while (low >> nxt) != 1:
                            nxt += 1
                        dp[mask | low, nxt] += c
                        avail ^= low
            h = 0
            for last in range(n):
                h += dp[full, last]
            if h < 1024:
                seen[t, h] = 1
            if h > best[t]:
                best[t] = h
                bestcnt[t] = 1
            elif h == best[t]:
                bestcnt[t] += 1
    m = best.max()
    cnt = 0
    for t in range(nthreads):
        if best[t] == m:
            cnt += bestcnt[t]
    allseen = np.zeros(1024, np.uint8)
    for t in range(nthreads):
        for v in range(1024):
            if seen[t, v]:
                allseen[v] = 1
    return m, cnt, allseen


maxima = {}
spectrum = {}
for n in range(2, 9):
    m, cnt, seen = max_h_fixed_path(n, 3)
    maxima[n] = int(m)
    spectrum[n] = [v for v in range(1024) if seen[v]]
check("6a maximum of H over all tournaments on n = 2..8 vertices: 1, 3, 5, 15, 45, 189, 661 (exhaustive)",
      [maxima[n] for n in range(2, 9)] == [1, 3, 5, 15, 45, 189, 661], str(maxima))
check("6b hence no tournament on <= 8 vertices has H = 1039 (661 < 1039): the least order is >= 9", max(maxima.values()) < 1039)
check("6c sanity: every H value at n <= 8 is odd, and 7, 21 never occur (THM-1370 at these orders)",
      all(v % 2 == 1 for n in spectrum for v in spectrum[n]) and all(7 not in spectrum[n] and 21 not in spectrum[n] for n in spectrum))
print("      number of distinct H values at n = 8: %d; odd values <= 661 missing at n = 8: %d" %
      (len(spectrum[8]), len([v for v in range(1, 662, 2) if v not in spectrum[8]])))

# ---------------------------------------------------------------- 7: the realiser of order 9
MASKS = [496, 317, 473, 449, 392, 28, 434, 290, 32]
n = 9
is_tournament = all(((MASKS[i] >> j) & 1) + ((MASKS[j] >> i) & 1) == 1 for i in range(n) for j in range(i + 1, n)) \
    and all((MASKS[i] >> i) & 1 == 0 for i in range(n)) and all(m < (1 << n) for m in MASKS)
check("7a the masks define a tournament on 9 vertices (bit j of masks[i] set iff i -> j)", is_tournament)
adj = np.zeros((n, n), np.uint8)
for i in range(n):
    for j in range(n):
        adj[i, j] = (MASKS[i] >> j) & 1
h_dp = int(ham_paths(adj, n))
h_perm = 0
for perm in itertools.permutations(range(n)):
    okp = True
    for a, b in zip(perm, perm[1:]):
        if not (MASKS[a] >> b) & 1:
            okp = False
            break
    if okp:
        h_perm += 1
check("7b the realiser has exactly 1039 Hamiltonian paths (DP and all 9! permutations)", h_dp == 1039 and h_perm == 1039,
      "DP %d, permutations %d" % (h_dp, h_perm))
scores = sorted(bin(m).count("1") for m in MASKS)
print("      score sequence of the realiser:", scores)
check("7c least order of a tournament with H = 1039 is exactly 9", h_dp == 1039 and max(maxima.values()) < 1039)

if all(RESULTS):
    print("ALL CHECKS PASSED (%d checks)" % len(RESULTS))
    sys.exit(0)
print("CHECK FAILED (%d of %d)" % (RESULTS.count(False), len(RESULTS)))
sys.exit(1)
