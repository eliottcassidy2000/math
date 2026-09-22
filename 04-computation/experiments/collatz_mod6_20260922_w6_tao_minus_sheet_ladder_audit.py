#!/usr/bin/env python3
"""
Audit of lane tao_minus_sheet_ladder (wave 6, 2026-09-22).  Independent
recomputation of every number the note relies on, with different code paths
from the explorer's script where a different path exists, plus the boundary
and definition-mismatch probes listed in the note's audit record.

A1  Delta = 4 ladder, mod-3 obstruction, pasted Lean sets.
A2  Offset identity by the recursion F_n = (3 F_{n-1} + 1)/2^{a_n} (exact),
    both signs; Terras form of the word law: every word with |a| = s is realised
    by exactly one residue class of odd N mod 2^(s+1), on both sheets.
A3  Syrac negation mod 3^n with exact Geom(2) weights, n = 3 and n = 4; the
    'never a multiple of 3' fact is proved by the last summand 2^{-a_n}.
A4  Prime targets under two readings of 'target' (Syracuse image; raw 3n+b).
A5  PPT (77,36,85), 2^11 - 3^7, the seven-cycle, carry 2363, gate 17.
A6  196 = 14^2, (G12), A110979 recomputation (primes < 2*10^7).
A7  Basin census on odd n <= 10^6, both sheets, memoised C-form.
A8  Square-sum graph Q_n: components by union-find, low-degree vertices,
    the {1,4,6} membership test, Hamiltonian path COUNTS by subset DP for
    n <= 23 (cross-checked with A071983), path existence by an independent
    backtracker for n <= 34, cycle existence for n = 30..34 (A071984).
A9  Tao text anchors from the cached pdftotext output (Remark 5.1 line, the
    quotation line, the 'artificial' convention, the 'non-negative' line).
A10 Lean loopless witnesses.
Explicit raise everywhere.
"""
import os
import re
import sys
import time
from fractions import Fraction
from math import gcd, isqrt

from sympy import isprime, primerange

TXT = ("/private/tmp/claude-501/-Users-e-Documents-GitHub-math/"
       "e197ec98-d8f9-4475-947b-5af87889cf35/scratchpad/tao/tao.txt")


def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def hr(t):
    print()
    print("-" * 78)
    print(t)
    print("-" * 78)


def v2(m):
    return (m & -m).bit_length() - 1


def syr(N, b):
    m = 3 * N + b
    a = v2(m)
    return m >> a, a


T0 = time.time()

# ---------------------------------------------------------------- A1
hr("A1 Delta = 4 ladder")
lad = [3, 7, 11, 17]
print("differences:", [lad[i + 1] - lad[i] for i in range(3)], "AP:", len({lad[i + 1] - lad[i] for i in range(3)}) == 1)
runs = {}
for p in primerange(2, 10 ** 6):
    L = 1
    while isprime(p + 4 * L):
        L += 1
    runs[L] = runs.get(L, 0) + 1
    if L >= 3:
        print("run of length", L, "starting at", p, ":", [p + 4 * i for i in range(L)])
print("run-length histogram of maximal difference-4 prime runs (start < 10^6):", dict(sorted(runs.items())))
check(max(runs) == 3, "max run length 3")
print("shifted set:", sorted({p + 4 for p in lad}), "primality:", {p + 4: bool(isprime(p + 4)) for p in lad})

# ---------------------------------------------------------------- A2
hr("A2 offset identity by recursion; Terras form of the word law")
NMAX, nmax = 4096, 8
cnt = 0
for b in (1, -1):
    for N in range(1, NMAX, 2):
        x, F, s = N, Fraction(0), 0
        for n in range(1, nmax + 1):
            x, a = syr(x, b)
            s += a
            F = (3 * F + 1) / Fraction(2 ** a)      # F_n = (3 F_{n-1} + 1) / 2^{a_n}
            check(x == Fraction(3 ** n * N, 2 ** s) + b * F, "recursion identity b=%d N=%d n=%d" % (b, N, n))
            cnt += 1
print("recursion identity T_b^n(N) = 3^n 2^-|a| N + b F_n(a): instances checked", cnt, "(odd N <", NMAX, ", n <=", nmax, ")")
J = 10
K = J + 3          # odd N in [1, 2^K)
for b in (1, -1):
    real = {}
    for N in range(1, 1 << K, 2):
        x, s, w = N, 0, []
        while True:
            x, a = syr(x, b)
            s += a
            if s > J:
                break
            w.append(a)
            real.setdefault(tuple(w), []).append(N)
    check(len(real) == 2 ** J - 1, "number of words with |a| <= J is 2^J - 1")
    for w, Ns in real.items():
        s = sum(w)
        m = 2 ** (s + 1)
        check(len({N % m for N in Ns}) == 1, "word %s is one class mod 2^(s+1)" % (w,))
        check(len(Ns) == 2 ** (K - 1 - s), "count of word %s" % (w,))
    print("b=%+d: all %d words with |a| <= %d: each is exactly one residue class of odd N mod 2^(|a|+1), realised by 2^(%d-|a|) odd N < 2^%d"
          % (b, len(real), J, K - 1, K))
    if b == 1:
        plus = {w: sorted(N % 2 ** (sum(w) + 1) for N in Ns)[0] for w, Ns in real.items()}
    else:
        minus = {w: sorted(N % 2 ** (sum(w) + 1) for N in Ns)[0] for w, Ns in real.items()}
same_class = sum(1 for w in plus if plus[w] == minus[w])
neg_class = sum(1 for w in plus if (plus[w] + minus[w]) % 2 ** (sum(w) + 1) == 0)
print("residue classes: identical on both sheets for %d words; negatives of each other (r_- = -r_+ mod 2^(|a|+1)) for %d words of %d"
      % (same_class, neg_class, len(plus)))
check(neg_class == len(plus), "minus-sheet class is the negative of the plus-sheet class")

# ---------------------------------------------------------------- A3
hr("A3 Syrac negation mod 3^n with exact Geom(2) weights")
for n in (3, 4):
    mod = 3 ** n
    SMAX = 4 * n + 4
    hist = {1: {}, -1: {}}
    tot = Fraction(0)

    def words(n, smax):
        if n == 0:
            yield ()
            return
        for a in range(1, smax - (n - 1) + 1):
            for rest in words(n - 1, smax - a):
                yield (a,) + rest

    for w in words(n, SMAX):
        s = sum(w)
        wt = Fraction(1, 2 ** s)
        tot += wt
        Fn = Fraction(0)
        for a in w:
            Fn = (3 * Fn + 1) / Fraction(2 ** a)
        r = (Fn.numerator * pow(Fn.denominator, -1, mod)) % mod
        for b in (1, -1):
            hist[b][(b * r) % mod] = hist[b].get((b * r) % mod, 0) + wt
    neg = {(-r) % mod: c for r, c in hist[1].items()}
    check(neg == hist[-1], "negation n=%d" % n)
    m3 = [r for r in hist[1] if r % 3 == 0]
    print("n=%d, words |a| <= %d (mass %s of 1): minus histogram mod %d = negation of plus histogram: True; residues hit %d of %d; multiples of 3 hit: %d"
          % (n, SMAX, tot, mod, len(hist[1]), mod, len(m3)))
print("proof that 3 never divides F_n(a): F_n = 3*(...) + 2^{-a_n}, and 2^{-a_n} is a unit mod 3")

# ---------------------------------------------------------------- A4
hr("A4 prime targets, two readings")
for b in (1, -1):
    first_syr = next((n, syr(n, b)[0]) for n in range(1, 10 ** 4, 2) if syr(n, b)[0] not in (1, n) and isprime(syr(n, b)[0]))
    pre3 = [n for n in range(1, 1 << 16, 2) if syr(n, b)[0] == 3]
    raw_first = [(n, 3 * n + b) for n in range(1, 20) if isprime(3 * n + b)][:3]
    raw3 = [n for n in range(1, 10 ** 5) if 3 * n + b == 3]
    print("b=%+d: Syracuse reading: first odd n with prime target other than 1, n: %s; preimages of 3 (odd n < 2^16): %s"
          % (b, first_syr, pre3))
    print("b=%+d: raw reading 3n%+d, n >= 1: first prime values %s; n with 3n%+d = 3: %s (3n%+d = %+d mod 3)"
          % (b, b, raw_first, b, raw3, b, b))

# ---------------------------------------------------------------- A5
hr("A5 PPT, 2^11 - 3^7, seven-cycle")
print("77^2+36^2 =", 77 ** 2 + 36 ** 2, "85^2 =", 85 ** 2, "gcd(77,36) =", gcd(77, 36))
m, n_ = (7 + 11) // 2, (11 - 7) // 2
print("(m,n) = (%d,%d): m^2-n^2 = %d, 2mn = %d, m^2+n^2 = %d" % (m, n_, m * m - n_ * n_, 2 * m * n_, m * m + n_ * n_))
print("T_+(7) =", syr(7, 1), "(target 11, valuation 1)")
D = 2 ** 11 - 3 ** 7
print("2^11 =", 2 ** 11, "3^7 =", 3 ** 7, "Delta =", D, "; (-17)*Delta =", -17 * D, "; 17*139 =", 17 * 139, "; factor 2363 =", [p for p in range(2, 2363) if 2363 % p == 0 and isprime(p)])
x, cyc, word = 17, [17], []
for _ in range(7):
    x, a = syr(x, -1)
    cyc.append(x)
    word.append(a)
print("3n-1 odd cycle from 17:", cyc, "valuations", word, "|a| =", sum(word))
check(cyc[-1] == 17 and sum(word) == 11, "seven-cycle")
B = -D * 17
print("carry B = -(Delta)*17 =", B, "; gate n_0 = b*B/Delta with b = -1:", (-B) // D, "; q = |Delta|/gcd(B,|Delta|) =", abs(D) // gcd(B, abs(D)))
# cycle identity: 2^|a| n_0 = 3^L n_0 + b B  <=>  (2^|a| - 3^L) n_0 = b B
check((2 ** 11 - 3 ** 7) * 17 == -B, "cycle identity")
print("cycle identity (2^11 - 3^7)*17 = %d = b*B with b = -1: True" % ((2 ** 11 - 3 ** 7) * 17))

# ---------------------------------------------------------------- A6
hr("A6 196, (G12), A110979")
op = list(primerange(3, 40))[:11]
S = 1 + sum(op)
print("1 + sum of first 11 odd primes", op, "=", S, "=", isqrt(S), "^2")
check(S == 196, "196")
c = [(p - (2 * i + 1)) // 2 for i, p in enumerate(op, 1)]
comp_count = [sum(1 for q in range(3, p + 1, 2) if not isprime(q)) for p in op]
print("c_i =", c, "sum", sum(c), "; direct count of odd composites in [3, p_i]:", comp_count, "equal:", c == comp_count)
check(c == comp_count and 12 ** 2 + 2 * sum(c) == 196, "G12")
print("(G12): 12^2 + 2*%d = %d" % (sum(c), 144 + 2 * sum(c)))
tot, hits = 0, []
for mm, p in enumerate(primerange(2, 20_000_000), 1):
    tot += p
    r = isqrt(tot - 1)
    if r * r == tot - 1:
        hits.append((mm, p, tot - 1, r))
print("sum of first m primes minus 1 is a square, primes < 2*10^7:", hits)
oeis = [1, 4, 9, 16, 196, 839056, 7796654478001]
print("A110979 name (lookup 2026-09-22): 'Squares equal to the sum of the first k primes minus 1.'; data", oeis,
      "; comment: no more terms < 7472966967498 (sum of first 10^6 primes minus 1, Ray Chandler 2005)")
check([h[2] for h in hits] == oeis, "A110979 match")
print("index of 196 in A110979:", oeis.index(196) + 1, "(1-based)")

# ---------------------------------------------------------------- A7
hr("A7 basin census, memoised C-form, odd n <= 10^6")
for b in (1, -1):
    cyc_set = {1} if b == 1 else {1, 5, 17}
    root = {}
    counts = {}
    for n in range(1, 10 ** 6 + 1, 2):
        x = n
        while x >= n and x not in cyc_set:
            x = 3 * x + b
            while x % 2 == 0:
                x //= 2
        r = root[x] if x < n and x in root else x
        root[n] = r
        counts[r] = counts.get(r, 0) + 1
    print("b=%+d: basins" % b, dict(sorted(counts.items())), "fractions",
          {k: round(v / 500000, 6) for k, v in sorted(counts.items())})
    if b == -1:
        check(counts == {1: 163486, 5: 162122, 17: 174392}, "minus basins")
    else:
        check(counts == {1: 500000}, "plus basins")
print("agrees with minus_sheet_positive_control S1 row 10^6 (T-form odd counts 163486 / 162122 / 174392)")
print("[A7 %.1fs]" % (time.time() - T0))

# ---------------------------------------------------------------- A8
hr("A8 square-sum graph Q_n")


def edges(n):
    return [(x, y) for x in range(1, n + 1) for y in range(x + 1, n + 1) if isqrt(x + y) ** 2 == x + y]


def uf_components(n):
    par = list(range(n + 1))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for x, y in edges(n):
        par[f(x)] = f(y)
    comps = {}
    for v in range(1, n + 1):
        comps.setdefault(f(v), []).append(v)
    return sorted(comps.values())


cc = {n: len(uf_components(n)) for n in range(1, 41)}
print("component counts n=1..40:", cc)
check(all(cc[n] == 3 for n in range(4, 13)) and cc[13] == 2 and all(cc[n] == 1 for n in range(14, 41)), "components")
c12 = uf_components(12)
print("Q_12 components:", c12, "founders (minima):", [c[0] for c in c12])
same16 = [n for n in range(6, 41) if any(1 in c and 6 in c for c in uf_components(n))]
sep14 = [n for n in range(4, 41) if not any(1 in c and 4 in c for c in uf_components(n))]
print("n in [6,40] with 1 and 6 in the same component of Q_n:", "all of 6..40" if same16 == list(range(6, 41)) else same16)
print("n in [4,40] with 1 and 4 in different components:", sep14, "(1 and 4 join at n = 13)")
check(same16 == list(range(6, 41)) and sep14 == list(range(4, 13)), "{1,4,6} membership")
print("so {1,4,6} is NOT a transversal of the three Q_n chains: 1 ~ 3 ~ 6 (1+3 = 4, 3+6 = 9); the chain {2,7,9} contains none of 1,4,6")
adjs = {}
for n in range(1, 41):
    adj = {v: set() for v in range(1, n + 1)}
    for x, y in edges(n):
        adj[x].add(y)
        adj[y].add(x)
    adjs[n] = adj
low = {n: [v for v in adjs[n] if len(adjs[n][v]) <= 1] for n in range(14, 33)}
print("degree <= 1 vertices n=14..32:", low)
check(low[18] == [16, 17, 18] and low[19] == [16, 18] and all(low[n] == [18] for n in range(20, 31)) and low[31] == [] and low[32] == [], "low degree")
print("N(18) in Q_30:", sorted(adjs[30][18]), "in Q_31:", sorted(adjs[31][18]), "(18+18 = 36 is excluded as x != y; 18+31 = 49)")
print("degree-0 vertices for n <= 13:", {n: [v for v in adjs[n] if not adjs[n][v]] for n in range(1, 14) if any(not adjs[n][v] for v in adjs[n])})


def count_paths(n):
    """Number of Hamiltonian paths (reversals identified) by subset DP on reachable states."""
    adj = adjs[n]
    nb = {v: sum(1 << (u - 1) for u in adj[v]) for v in adj}
    dp = {}
    for v in range(1, n + 1):
        dp[(1 << (v - 1), v)] = 1
    layer = dict(dp)
    for _ in range(n - 1):
        nxt = {}
        for (mask, v), c in layer.items():
            free = nb[v] & ~mask
            while free:
                lb = free & -free
                u = lb.bit_length()
                free ^= lb
                key = (mask | lb, u)
                nxt[key] = nxt.get(key, 0) + c
        layer = nxt
        if len(layer) > 6_000_000:
            raise RuntimeError("DP too large at n=%d" % n)
    full = (1 << n) - 1
    tot = sum(c for (mask, v), c in layer.items() if mask == full)
    check(tot % 2 == 0 or n == 1, "path count parity")
    return tot // 2 if n > 1 else tot


t1 = time.time()
pc = {n: count_paths(n) for n in range(1, 26)}
print("Hamiltonian path counts (reversal identified) n=1..25 by subset DP: %s [%.1fs]" % (pc, time.time() - t1))
a071983 = {15: 1, 16: 1, 17: 1, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 3, 24: 0, 25: 10}
print("A071983 (offset 15, lookup 2026-09-22) a(15..25):", a071983, "match:", all(pc[n] == a071983[n] for n in a071983))
check(all(pc[n] == a071983[n] for n in a071983), "A071983")
check(all(pc[n] == 0 for n in range(2, 15)), "no path for 2..14")


def ham_exists(n, cycle=False, budget=60.0):
    """Independent backtracker: DFS from each start (or from 1 for cycles), pruning when the
    unvisited vertices plus the current end are disconnected or when two unvisited vertices
    other than the required endpoint have remaining degree <= 1."""
    adj = adjs[n]
    nb = {v: sum(1 << (u - 1) for u in adj[v]) for v in adj}
    full = (1 << n) - 1
    t0 = time.time()

    def connected(mask, v):
        seen = 1 << (v - 1)
        stack = [v]
        while stack:
            u = stack.pop()
            f = nb[u] & mask & ~seen
            while f:
                lb = f & -f
                f ^= lb
                seen |= lb
                stack.append(lb.bit_length())
        return seen & mask == mask

    def dfs(v, mask, start):
        if mask == full:
            return (not cycle) or (nb[v] >> (start - 1)) & 1
        if time.time() - t0 > budget:
            raise TimeoutError
        rem = full & ~mask
        if not connected(rem, v):
            return False
        ones = 0
        r = rem
        while r:
            lb = r & -r
            r ^= lb
            u = lb.bit_length()
            d = bin(nb[u] & rem).count("1") + ((nb[u] >> (v - 1)) & 1) + (1 if cycle and (nb[u] >> (start - 1)) & 1 else 0)
            if d == 0:
                return False
            if d == 1:
                ones += 1
                if ones > (0 if cycle else 1):
                    return False
        f = nb[v] & rem
        while f:
            lb = f & -f
            f ^= lb
            if dfs(lb.bit_length(), mask | lb, start):
                return True
        return False

    try:
        starts = [1] if cycle else sorted(adj, key=lambda v: len(adj[v]))
        for s in starts:
            if dfs(s, 1 << (s - 1), s):
                return "yes"
        return "no"
    except TimeoutError:
        return "TIMEOUT"


sys.setrecursionlimit(10000)
t1 = time.time()
ex = {n: ham_exists(n) for n in range(1, 35)}
print("Hamiltonian path existence n=1..34 (independent backtracker): %s [%.1fs]" % (ex, time.time() - t1))
yes = [n for n in ex if ex[n] == "yes"]
print("   yes:", yes)
print("   no :", [n for n in ex if ex[n] == "no"], "unresolved:", [n for n in ex if ex[n] == "TIMEOUT"])
check(yes == [1, 15, 16, 17, 23] + list(range(25, 35)), "path existence")
t1 = time.time()
cy = {n: ham_exists(n, cycle=True, budget=90.0) for n in range(30, 35)}
print("Hamiltonian cycle existence n=30..34: %s [%.1fs]" % (cy, time.time() - t1))
print("A071984 (offset 32, lookup 2026-09-22): a(32..34) = 1, 1, 11; comment (Dobbelaere 2018): no cycle for n <= 30")
check(cy[30] == "no" and cy[31] == "no" and all(cy[n] == "yes" for n in (32, 33, 34)), "cycles")
p15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
p23 = [18, 7, 9, 16, 20, 5, 11, 14, 2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22]
for n, p in ((15, p15), (23, p23)):
    ok = sorted(p) == list(range(1, n + 1)) and all(isqrt(p[i] + p[i + 1]) ** 2 == p[i] + p[i + 1] for i in range(n - 1))
    print("pasted path n=%d valid: %s; ends %s; sums %s" % (n, ok, (p[0], p[-1]), sorted({p[i] + p[i + 1] for i in range(n - 1)})))
    check(ok, "pasted path")
print("n=15: the unique path (count 1) has ends 8, 9 = the two leaves of Q_15; vertex 4 has neighbours", sorted(adjs[15][4]), "in Q_15")

# ---------------------------------------------------------------- A9
hr("A9 Tao text anchors (cached pdftotext of arXiv:1909.03562v7)")
if os.path.exists(TXT):
    with open(TXT, encoding="utf-8", errors="replace") as fh:
        L = fh.read().split("\n")
    print("lines:", len(L))
    find = lambda pat: [i + 1 for i, l in enumerate(L) if re.search(pat, l)]
    print("Remark 5.1 (Syr_min(N) <= N^theta, Korec recovery) at line:", find(r"^Remark 5\.1"))
    print("'eventually recover the results of Korec' at line:", find(r"eventually recover the results of Korec"))
    print("quotation 'beyond (1.19), (1.20)' at line:", find(r"beyond \(1\.19\), \(1\.20\)"))
    print("'artificial' convention Syr^infinity(N) := 1 at line:", find(r"artificial"))
    print("'negative' lines:", find(r"\bnegative\b"), "; of these 'non-negative':", find(r"non-negative"))
    print("Proposition 1.17 header line:", find(r"^Proposition 1\.17"), "; 'not divisible by 3' within 3 lines:",
          any("not divisible by 3" in L[i] for i in range(find(r"^Proposition 1\.17")[0] - 1, find(r"^Proposition 1\.17")[0] + 3)))
    print("Syrac definition (1.22) 'Fn (Geom(2)n ) mod 3n' at line:", find(r"Syrac\(Z/3n Z\) ≡ Fn \(Geom\(2\)n \) mod 3n"))
    print("footnote on negative time indexing / reversal at lines:", find(r"ancient|necessitates some reversal"))
    print("Korec exponent line (theta > log 3 / log 4 ~ 0.7924):", find(r"0\.7924"))
    print("Section 1.3 spans roughly lines", find(r"^1\.3\. ")[:1], "to", find(r"^1\.4\. ")[:1])
    print("martingale / dyadic / entropy decrement lines:", find(r"martingale"), find(r"\bdyadic\b"), find(r"entropy decrement"))
    print("black-triangle renewal lines (first three):", find(r"black")[:3], "; Pascal pairing b_j = a_(2j-1) + a_(2j) at line:", find(r"bj := a2j−1 \+ a2j"))
else:
    print("cached text absent; A9 skipped")

# ---------------------------------------------------------------- A10
hr("A10 Lean loopless witnesses")
print("v <= 40 with 2v a square:", [v for v in range(1, 41) if isqrt(2 * v) ** 2 == 2 * v], "(v = 2k^2)")
print("smallest: v = 2, Fin index 1; note also Fin n for n = 1: value 1, 1+1 = 2 not a square, so the n = 1 graph is loopless")
print()
print("next-question parameter recorded for the note: computable alpha such as 1.1 in place of 1.001")
print("audit total %.1fs" % (time.time() - T0))
