#!/usr/bin/env python3
"""collatz_crossings_20260926_audit.py -- adversarial audit of
05-knowledge/results/collatz_crossings_20260926_potential_and_seeds.md (session collatz-crossings-20260926)
and its scripts collatz_crossings_20260926_phase.py, collatz_poset_20260926_balance.py,
collatz_wythoff_20260926_interlock_probe.py.  Independent recomputation, exact where the claim is exact.

Sections (numbered as in the audit brief):
 1. Wythoff colouring: AA/B/AB = Zeckendorf lowest-index rule (exact floor(k phi) via isqrt), densities,
    golden-coordinate characterisation, doubling matrix (empirical + exact derivation), mutual information.
 2. {2,3,11}: primes p with 2p below the next odd square above p, to 10^6; the algebra 4j^2-12j+1<0.
 3. No-descent posets P(o,e): f(j), linear-extension DP vs word enumeration, balance table k<=20 vs .out,
    the four claims of section 2.2, W_k of THM-4495, width 2.
 4. Proposition 1: phase identity, discrepancies (base 2 and base 10), pure rotation comparison, shift
    convergence, and a synthetic sequence showing base-2 equidistribution does not give base 10.
 5. Proposition 2: exact Q_l identities, valuations, increments, the reciprocal-sum constant.
 6. Proposition 3: 3-adic residue formula, m_l<3^l iff 2^d_l>Q_l, locality mod 3^j.
 7. Exits lemma: X_0(b), M_k(y), stays on the T-orbits, crossing counts.
 8. Section 3.6 binary tree on odd numbers, enumeration to 10^5.
 9. Text checks (labels, Proposition 4, Benford wording, speculation markers).
Usage: python3 collatz_crossings_20260926_audit.py   (writes 05-knowledge/results/collatz_crossings_20260926_audit.out)
"""
import math, os, re, sys
from fractions import Fraction
from decimal import Decimal, getcontext
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
NOTE = os.path.join(ROOT, '05-knowledge', 'results', 'collatz_crossings_20260926_potential_and_seeds.md')
BAL_OUT = os.path.join(HERE, 'collatz_poset_20260926_balance.out')
OUT_PATH = os.path.join(ROOT, '05-knowledge', 'results', 'collatz_crossings_20260926_audit.out')

LINES = []
def out(s=""):
    print(s); LINES.append(s)

getcontext().prec = 60
SQRT5 = Decimal(5).sqrt()
PHI = (1 + SQRT5) / 2
ALPHA = math.log2(3)


def floor_kphi(k):
    """floor(k*phi) exactly: phi = (1+sqrt5)/2, so k*phi = (k + sqrt(5k^2))/2."""
    return (k + math.isqrt(5 * k * k)) // 2


def v2(x):
    c = 0
    while x % 2 == 0:
        x //= 2; c += 1
    return c


def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def syracuse(n, steps, q=3):
    """odd iterates m_0=n, m_(l+1)=(q m_l+1)/2^v_l; d_l = v_0+...+v_(l-1); stops after reaching 1 (q=3)."""
    ms = [n]; ds = [0]; vs = []; d = 0
    for _ in range(steps):
        x = q * n + 1; v = 0
        while x % 2 == 0:
            x //= 2; v += 1
        d += v; n = x
        ms.append(n); ds.append(d); vs.append(v)
        if n == 1: break
    return ms, ds, vs


def star_discrepancy(xs):
    xs = sorted(x % 1 for x in xs); N = len(xs)
    return max(max(abs((i + 1) / N - x), abs(i / N - x)) for i, x in enumerate(xs))


def entropy(counter, N):
    return -sum((c / N) * math.log2(c / N) for c in counter.values() if c)


def mutual_info(pairs):
    N = len(pairs); cx, cy, cxy = Counter(), Counter(), Counter()
    for x, y in pairs:
        cx[x] += 1; cy[y] += 1; cxy[(x, y)] += 1
    I = sum((c / N) * math.log2(c * N / (cx[x] * cy[y])) for (x, y), c in cxy.items())
    return I, entropy(cy, N), len(cx), len(cy)


# ----------------------------------------------------------------------------------------------------
def section1():
    out("=" * 100)
    out("1. Wythoff colouring (note section 1.1)")
    LIM = 2_000_010
    # exact Wythoff sets
    col = bytearray(LIM + 1)          # 1 red = AA, 2 black = B, 3 blue = AB
    hits = bytearray(LIM + 1)
    A = []
    k = 1
    while True:
        a = floor_kphi(k)
        if a > LIM: break
        A.append(a); k += 1
    B = []
    k = 1
    while True:
        b = floor_kphi(k) + k         # floor(k phi^2) = floor(k phi) + k
        if b > LIM: break
        B.append(b); k += 1
    for b in B:
        col[b] = 2; hits[b] += 1
    for a in A:
        aa = floor_kphi(a)
        if aa <= LIM: col[aa] = 1; hits[aa] += 1
    for b in B:
        ab = floor_kphi(b)
        if ab <= LIM: col[ab] = 3; hits[ab] += 1
    partition_ok = all(hits[n] == 1 for n in range(1, LIM + 1))
    out("   AA, B, AB computed exactly (floor(k phi) = (k + isqrt(5k^2))//2); they partition [1, %d]: %s" % (LIM, partition_ok))

    # Zeckendorf lowest index rule
    F = [1, 2]
    while F[-1] <= LIM: F.append(F[-1] + F[-2])
    def zeck_lowest(n):
        i = len(F) - 1; low = None
        while n > 0:
            while F[i] > n: i -= 1
            low = i + 2; n -= F[i]; i -= 2
        return low
    def zeck_colour(n):
        k = zeck_lowest(n)
        return 1 if k == 2 else (2 if k % 2 == 1 else 3)
    ZLIM = 400_000
    mism = [n for n in range(1, ZLIM + 1) if zeck_colour(n) != col[n]]
    out("   Zeckendorf rule (lowest index 2 -> red, odd -> black, even>=4 -> blue) == (AA, B, AB) for n <= %d: %s (mismatches: %d)" % (ZLIM, not mism, len(mism)))

    names = {1: 'red', 2: 'black', 3: 'blue'}
    first8 = [names[col[n]] for n in range(1, 9)]
    out("   colours of n = 1..8: %s  (note quotes: red, black, blue, red, black, red, black, blue)" % ", ".join(first8))
    out("   colour(35) = %s (Zeckendorf 35 = 34 + 1, lowest index 2); colour(37) = %s (37 = 34 + 3, lowest index 4)" % (names[col[35]], names[col[37]]))
    AB_list = [n for n in range(1, 60) if col[n] == 3]
    out("   AB (blue) up to 59: %s" % AB_list)
    gaps = [AB_list[i + 1] - AB_list[i] for i in range(len(AB_list) - 1)]
    # Fibonacci word over {5,3}: a=5, b=3 with a -> ab, b -> a
    w = "a"
    while len(w) < 40:
        w = "".join("ab" if c == "a" else "a" for c in w)
    fibword = [5 if c == "a" else 3 for c in w[:len(gaps)]]
    out("   blue gaps: %s; Fibonacci word (5=a,3=b; a->ab, b->a): %s; equal: %s" % (gaps, fibword, gaps == fibword))
    out("   NOTE: the owner's full 35-term list is not in the repository (only these 8 terms are quoted in the note), so")
    out("         'agrees for n <= 34' cannot be audited independently; the checkable part (n <= 8, 35 red, 37 blue) holds.")

    # densities
    N = 1_000_000
    cnt = Counter(col[n] for n in range(1, N + 1))
    out("   densities at N = 10^6: red %.5f black %.5f blue %.5f; phi^-2 = %.5f, phi^-3 = %.5f" % (
        cnt[1] / N, cnt[2] / N, cnt[3] / N, float(PHI ** -2), float(PHI ** -3)))

    # golden coordinate characterisation
    t1 = 1 / PHI ** 2; t2 = 1 - 1 / PHI ** 3
    def golden_colour(n):
        x = (Decimal(n) * PHI) % 1
        return 2 if x < t1 else (1 if x < t2 else 3)
    GLIM = 100_000
    gm = sum(1 for n in range(1, GLIM + 1) if golden_colour(n) != col[n])
    out("   golden-coordinate rule: x = {n phi}: black iff x < phi^-2 = %.6f, red iff phi^-2 < x < 1 - phi^-3 = %.6f, blue iff x > 1 - phi^-3;" % (float(t1), float(t2)))
    out("      agrees with (AA, B, AB) for n <= %d: %s (mismatches %d)" % (GLIM, gm == 0, gm))

    # doubling matrix, empirical and exact
    M = [[0] * 4 for _ in range(4)]
    for n in range(1, N + 1):
        M[col[n]][col[2 * n]] += 1
    out("   doubling matrix colour(n) -> colour(2n), N = 10^6, row-normalised (red, black, blue):")
    for i in (1, 2, 3):
        s = sum(M[i])
        out("      %-6s %s" % (names[i], "  ".join("%.5f" % (M[i][j] / s) for j in (1, 2, 3))))
    out("   exact derivation ({2n phi} = {2x}, x = {n phi} uniform on the colour interval):")
    out("      black x in (0, p2) -> 2x in (0, 2 p2 = 1 - p3): black iff 2x < p2 (half), red otherwise (half), never blue")
    out("      blue  x in (1 - p3, 1) -> 2x mod 1 in (1 - 2 p3, 1): red on (1 - 2 p3, 1 - p3) (half), blue on (1 - p3, 1) (half), never black")
    out("      red   x in (p2, 1 - p3) -> 2x in (1 - p3, 1) [blue] u [1, 1 + p2) [black] u (1 + p2, 2 - 2 p3) [red]:")
    p1, p2, p3 = float(PHI ** -1), float(PHI ** -2), float(PHI ** -3)
    rr = (p1 - 2 * p3) / 2 / p2; rb = (p2 / 2) / p2; rbl = (p3 / 2) / p2
    out("            red->red = (phi^-1 - 2 phi^-3)/(2 phi^-2) = phi^-2/2 = %.5f, red->black = 1/2 = %.5f, red->blue = phi^-3/(2 phi^-2) = phi^-1/2 = %.5f" % (rr, rb, rbl))
    # exact MI of the doubling channel
    P = {1: p2, 2: p2, 3: p3}
    Trans = {1: {1: p2 / 2, 2: 0.5, 3: p1 / 2}, 2: {1: 0.5, 2: 0.5, 3: 0.0}, 3: {1: 0.5, 2: 0.0, 3: 0.5}}
    I_exact = sum(P[i] * Trans[i][j] * math.log2(Trans[i][j] / P[j]) for i in P for j in P if Trans[i][j] > 0)
    out("      exact I(colour(n); colour(2n)) from this matrix = %.5f bits (note: 0.369); H(colour) = %.4f bits" % (
        I_exact, -sum(p * math.log2(p) for p in P.values())))

    # mutual information with Collatz targets
    for NMI in (100_000, 300_000):
        stop = {1: 0}
        def stopping(n):
            path = []; m = n
            while m not in stop:
                path.append(m); m = T(m)
            s = stop[m]
            for x in reversed(path):
                s += 1; stop[x] = s
            return stop[n]
        rows = []
        for n in range(2, NMI + 1):
            c = col[n]
            m = n; desc = 0
            for _ in range(20):
                m = T(m)
                if m < n: desc = 1; break
            st = stopping(n) % 3
            if n % 2 == 1:
                x = 3 * n + 1; v = v2(x); U = x >> v
                rows.append((c, n % 2, v % 3, desc, st, col[U] if U <= LIM else -1, col[2 * n]))
            else:
                rows.append((c, n % 2, -1, desc, st, -1, col[2 * n]))
        tg = [("parity", 1), ("v2(3n+1) mod 3 (odd n)", 2), ("descent<=20 T-steps", 3), ("T-steps-to-1 mod 3", 4), ("colour(U(n)) (odd n)", 5), ("colour(2n)", 6)]
        out("   mutual information I(colour(n); target) at N = %d [target entropy; plug-in noise floor (r-1)(c-1)/(2 N ln 2)]:" % NMI)
        for name, ti in tg:
            pairs = [(r[0], r[ti]) for r in rows if r[ti] >= 0]
            I, Hy, rr_, cc_ = mutual_info(pairs)
            floor_ = (rr_ - 1) * (cc_ - 1) / (2 * len(pairs) * math.log(2))
            out("      %-26s I = %.6f bits  [H = %.3f; noise ~ %.1e; n = %d]" % (name, I, Hy, floor_, len(pairs)))
        out("      (note claims: < 1e-5 for parity / v2 / descent / stopping mod 3 at N = 300000; 0.005 for colour(U(n)); 0.369 for colour(2n))")


# ----------------------------------------------------------------------------------------------------
def section2():
    out("=" * 100)
    out("2. {2, 3, 11}: primes p with 2p below the next odd square above p (note section 1.2)")
    LIM = 1_000_000
    sieve = bytearray([1]) * (LIM + 1); sieve[0] = sieve[1] = 0
    for i in range(2, math.isqrt(LIM) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    found = []
    for p in range(2, LIM + 1):
        if not sieve[p]: continue
        s = math.isqrt(p) + 1
        if s % 2 == 0: s += 1
        assert s * s > p and (s - 2) ** 2 <= p
        if 2 * p < s * s:
            found.append((p, s * s, [t for t in range(2 * p, s * s, p)]))
    out("   primes p <= 10^6 with 2p < next odd square: %s" % [(p, sq) for p, sq, _ in found])
    out("   multiples of p in [2p, square): %s" % [(p, mult) for p, _, mult in found])
    out("   algebra: 2(2j-1)^2 < (2j+1)^2  <=>  4j^2 - 12j + 1 < 0; values for j = 1..6: %s; roots %.4f, %.4f; so j in {1, 2}" % (
        [4 * j * j - 12 * j + 1 for j in range(1, 7)], (12 - math.sqrt(128)) / 8, (12 + math.sqrt(128)) / 8))
    out("   case j = 1: primes in (1, 9] with 2p < 9: %s; case j = 2: primes in (9, 25] with 2p < 25: %s" % (
        [p for p in (2, 3, 5, 7) if 2 * p < 9], [p for p in (11, 13, 17, 19, 23) if 2 * p < 25]))
    out("   ratio of consecutive odd squares (2j+1)^2/(2j-1)^2 for j = 2, 3, 4: %s (< 2 from j = 3 on)" % ["%.4f" % ((2 * j + 1) ** 2 / (2 * j - 1) ** 2) for j in (2, 3, 4)])
    out("   VERDICT: proposition and case analysis correct; the necessary condition 2(2j-1)^2 < (2j+1)^2 is then closed by the finite check in j = 1, 2.")


# ----------------------------------------------------------------------------------------------------
def f_of(j):
    o = 0
    while not (3 ** o > 2 ** (o + j)): o += 1
    return o


def nodescent_words(k):
    """words (1 = odd, 0 = halving) of length k with 3^(o') > 2^(k') for every nonempty prefix (exact integers)."""
    P3 = [3 ** i for i in range(k + 2)]; P2 = [2 ** i for i in range(k + 2)]
    res = []
    def rec(prefix, o, t):
        if t == k:
            res.append(tuple(prefix)); return
        if P3[o + 1] > P2[t + 1]:
            prefix.append(1); rec(prefix, o + 1, t + 1); prefix.pop()
        if P3[o] > P2[t + 1]:
            prefix.append(0); rec(prefix, o, t + 1); prefix.pop()
    rec([], 0, 0)
    return res


def ext_count_dp(o, e, f):
    """number of linear extensions of P(o,e): shuffles of O_1..O_o, H_1..H_e with H_j after O_(f(j))."""
    dp = [[0] * (e + 1) for _ in range(o + 1)]
    dp[0][0] = 1
    for i in range(o + 1):
        for j in range(e + 1):
            if i == 0 and j == 0: continue
            v = 0
            if i >= 1: v += dp[i - 1][j]
            if j >= 1 and i >= f[j]: v += dp[i][j - 1]
            dp[i][j] = v
    return dp[o][e]


def section3():
    out("=" * 100)
    out("3. No-descent posets, balance constants (note sections 2.1-2.2)")
    L32 = 1 / math.log2(1.5)
    f = {j: f_of(j) for j in range(1, 61)}
    formula_ok = all(f[j] == math.floor(j * L32) + 1 for j in f)
    out("   f(j) = least o' with 3^o' > 2^(o'+j): %s; equals floor(j log_(3/2) 2) + 1 for j <= 60: %s (THM-4503's m_j is the same quantity)" % (
        [f[j] for j in range(1, 13)], formula_ok))
    W_THM4495 = [1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990, 27328]
    # parse the session's .out
    reported = {}
    rx = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+O_(\d+) before H_(\d+): p = ([\d.]+)")
    with open(BAL_OUT) as fh:
        for line in fh:
            m = rx.match(line)
            if m:
                k, o, n, d, i, j, p = m.groups()
                reported[(int(k), int(o))] = (int(n), float(d), int(i), int(j), float(p))
    KMAX = 20
    table = {}      # (k,o) -> (count, delta, best pair, p, npairs)
    mism = []
    Wk = []
    for k in range(1, KMAX + 1):
        words = nodescent_words(k)
        Wk.append(len(words))
        by_o = defaultdict(list)
        for w in words: by_o[sum(w)].append(w)
        for o in sorted(by_o):
            ws = by_o[o]; e = k - o
            dpc = ext_count_dp(o, e, f)
            if dpc != len(ws): mism.append(("count", k, o, dpc, len(ws)))
            pairs = [(i, j) for i in range(1, o + 1) for j in range(1, e + 1) if i > f[j]]
            cnt = Counter()
            for w in ws:
                posO = [t for t, c in enumerate(w) if c == 1]; posH = [t for t, c in enumerate(w) if c == 0]
                for (i, j) in pairs:
                    if posO[i - 1] < posH[j - 1]: cnt[(i, j)] += 1
            best = (0.0, None, None)
            for (i, j) in pairs:
                p = cnt[(i, j)] / len(ws); bal = min(p, 1 - p)
                if bal > best[0] + 1e-15: best = (bal, (i, j), p)
            # width: an antichain of size 3 would contain two elements of one chain -> impossible; width 2 iff some pair incomparable
            width = 2 if pairs else 1
            table[(k, o)] = (len(ws), best[0], best[1], best[2], len(pairs), width, cnt)
            if k >= 4 and len(ws) >= 2:
                rep = reported.get((k, o))
                if rep is None:
                    mism.append(("missing in .out", k, o))
                else:
                    n_r, d_r, i_r, j_r, p_r = rep
                    p_mine_for_reported_pair = cnt[(i_r, j_r)] / len(ws)
                    if n_r != len(ws) or abs(d_r - best[0]) > 6e-5 or abs(p_r - p_mine_for_reported_pair) > 6e-5:
                        mism.append(("value", k, o, rep, (len(ws), best[0], best[1], best[2])))
    out("   linear-extension DP on P(o,e) == exact prefix enumeration for every cell, k <= %d: %s" % (KMAX, not any(m[0] == 'count' for m in mism)))
    out("   sum over o of e(P(o,e)) = W_k for k = 1..20: %s; equals THM-4495's list: %s" % (Wk, Wk == W_THM4495))
    out("   comparison with collatz_poset_20260926_balance.out (cells with >= 2 extensions, 4 <= k <= 20): %d cells compared, mismatches: %s" % (
        sum(1 for (k, o), t in table.items() if k >= 4 and t[0] >= 2), mism if mism else "none"))
    cells2 = [(k, o, t) for (k, o), t in table.items() if t[0] >= 2]
    out("   every cell with >= 2 extensions has width 2: %s; every such cell has delta >= 1/3: %s" % (
        all(t[5] == 2 for _, _, t in cells2), all(t[1] >= 1 / 3 - 1e-12 for _, _, t in cells2)))
    third = [(k, o, t[0]) for k, o, t in cells2 if abs(t[1] - 1 / 3) < 1e-12]
    out("   cells with delta = 1/3 exactly: %s  (note: exactly (5,4) and (6,4), three-word cells)" % third)
    for thr in (30, 100):
        sub = [(t[1], k, o, t[0]) for k, o, t in cells2 if t[0] >= thr]
        mn = min(sub)
        out("   cells with >= %d extensions: min delta = %.4f at (k,o) = (%d,%d) with %d extensions  (note: >= 0.441)" % (thr, mn[0], mn[1], mn[2], mn[3]))
    out("   best cell per k, distance from 1/2 (note claims 'for every k within 0.006 of 1/2'):")
    bad_k = []
    for k in range(4, KMAX + 1):
        cand = [(t[1], o) for (kk, o), t in table.items() if kk == k and t[0] >= 2]
        b = max(cand)
        flag = "" if abs(0.5 - b[0]) <= 0.006 else "   <-- VIOLATES the 0.006 claim"
        if flag: bad_k.append(k)
        out("      k=%2d best delta = %.4f (o=%d), |1/2 - delta| = %.4f%s" % (k, b[0], b[1], abs(0.5 - b[0]), flag))
    out("   k with best cell farther than 0.006 from 1/2: %s" % bad_k)
    # the two three-word cells by hand
    for (k, o) in ((5, 4), (6, 4)):
        ws = [w for w in nodescent_words(k) if sum(w) == o]
        out("   cell (%d,%d) words: %s" % (k, o, ["".join('O' if c else 'H' for c in w) for w in ws]))
    h = 0.9499555272
    out("   section 2.2(c) 'log_2 e(P) = h k - (3/2) log_2 k + O(1)': for the largest cell vs the whole level W_k (k, W_k, max cell, max cell/W_k, log_2 W_k - (hk - 1.5 log_2 k), same for max cell):")
    for k in (8, 12, 16, 20):
        mx = max(t[0] for (kk, o), t in table.items() if kk == k)
        base = h * k - 1.5 * math.log2(k)
        out("      k=%2d  W_k=%6d  max cell=%6d  ratio=%.3f  %.3f  %.3f  (a single cell is W_k/O(sqrt k) heuristically; the O(1) is loose for one poset)" % (
            k, Wk[k - 1], mx, mx / Wk[k - 1], math.log2(Wk[k - 1]) - base, math.log2(mx) - base))
    out("   Linial (1984): every finite width-2 poset that is not a chain has an incomparable pair with 1/3 <= P(x<y) <= 2/3; sharp on the 3-element poset a<b, c free (1/3).")
    out("   THM-4503 bridge note section 4 already states the bijection and checked the 1/3 bound in all fixed-count universes through ell = 14;")
    out("   new in the audited note: the exact delta values to k = 20, sharpness at (5,4),(6,4), the 0.441 floor for large cells.")


# ----------------------------------------------------------------------------------------------------
def section4():
    out("=" * 100)
    out("4. Proposition 1 (phases) -- recomputation and the base question")
    cases = [("3n+1 from 2^60-1", 2 ** 60 - 1, 3, 4000), ("3n+1 from 2^200-1", 2 ** 200 - 1, 3, 8000),
             ("3n+1 record 63728127", 63728127, 3, 4000), ("5n+1 from 7", 7, 5, 3000)]
    for name, n, q, steps in cases:
        ms, ds, vs = syracuse(n, steps, q)
        N = len(ms); aq = math.log2(q)
        ph2 = [math.log2(m) for m in ms]
        rot = [math.log2(n) + l * aq for l in range(N)]
        ph10 = [math.log10(m) for m in ms]
        shift = [ph2[l] - ph2[0] + ds[l] - l * aq for l in range(N)]       # = log2 C_l exactly (d_l integer)
        # independent: log2 C_l = sum_{j<l} log2(1 + 1/(q m_j))
        shift2 = [0.0]
        for j in range(N - 1): shift2.append(shift2[-1] + math.log2(1 + 1 / (q * ms[j])))
        agree = max(abs(a - b) for a, b in zip(shift, shift2))
        half = shift[N // 2:]
        acc_last20 = shift[-1] - shift[max(0, N - 21)]
        out("   %-22s N = %4d odd iterates: D_N{log_2 m} = %.4f, D_N{pure rotation log_2 n + l log_2 q} = %.4f, 1/sqrt(N) = %.4f, D_N{log_10 m} = %.4f" % (
            name, N, star_discrepancy(ph2), star_discrepancy(rot), 1 / math.sqrt(N), star_discrepancy(ph10)))
        out("        shift log2 C_l: at l=N/2: %.5f, at l=0.9N: %.5f, last: %.5f; max-min over last half: %.2e; accumulated in the last 20 odd steps: %.5f; identity vs product: %.1e" % (
            shift[N // 2], shift[int(0.9 * N)], shift[-1], max(half) - min(half), acc_last20, agree))
    out("   reading: for the finite 3n+1 orbits the 'convergence' of the shift is finiteness: it is ~0 while the iterates are huge and accumulates")
    out("   its whole value in the last few dozen odd steps (small m_j); only the (presumed divergent) 5n+1 orbit shows genuine convergence (1e-12).")
    out("   one-line reason for 'Weyl + convergent = equidistributed': for h != 0, |e(h(x_l+y_l)) - e(h(x_l+y))| <= 2 pi |h| |y_l - y| -> 0,")
    out("   so the Weyl sums of x_l + y_l and of the shifted rotation x_l + y have the same Cesaro limit 0.  CONFIRMED (base 2).")
    # base question
    out("   Base question: log_2 m_l = log_2 Q_l + l alpha - d_l with d_l an INTEGER, so {log_2 m_l} = {log_2 Q_l + l alpha} (pure rotation + convergent).")
    out("   In base 10: log_10 m_l = log_10 Q_l + l log_10 3 - d_l log_10 2, and d_l log_10 2 is NOT an integer: the phase is the linear form")
    out("   l log_10 3 - d_l log_10 2 along the orbit's own path (l, d_l); no rotation argument applies. Synthetic demonstration (a word with v_l in {1,2}")
    out("   chosen greedily, m_l := 3^l 2^(-d_l), exactly the identity's shape with Q constant):")
    l10_3 = math.log10(3); l10_2 = math.log10(2)
    Lsyn = 5000; y = 0.0; d = 0; ph2s = []; ph10s = []
    for l in range(Lsyn):
        ph2s.append((l * ALPHA - d) % 1); ph10s.append(y % 1)
        v = 1 if y + l10_3 - l10_2 <= 0.35 else 2
        y += l10_3 - v * l10_2; d += v
    out("      N = %d: D_N{log_2 m} = %.4f (equidistributed), D_N{log_10 m} = %.4f, base-10 phases confined to [%.3f, %.3f]" % (
        Lsyn, star_discrepancy(ph2s), star_discrepancy(ph10s), min(ph10s), max(ph10s)))
    out("   So base-2 significand equidistribution (what Prop. 1 proves) does not imply base-10 Benford; whether an actual divergent orbit can")
    out("   choose its valuations like this is unknown -- the note's proof does not address it. 'Benford's law' must be read as 'in base 2'. GAP.")
    out("   Attribution: Kontorovich-Miller (2005) and Lagarias-Soundararajan (2006) prove base-B Benford for the first N iterates of MOST seeds")
    out("   (unconditional, finite segments, the valuations acting as a random walk on the circle); Prop. 1 is a conditional statement about ONE")
    out("   non-eventually-periodic orbit (vacuous for 3n+1 if Collatz holds) in base 2 only, where the valuations drop out. Related but not the same content.")


# ----------------------------------------------------------------------------------------------------
def section5():
    out("=" * 100)
    out("5. Proposition 2 (the two-place clock Q_l)")
    for n in (27, 703, 26623):
        ms, ds, vs = syracuse(n, 500)
        L = len(ms) - 1
        S = 0; C = Fraction(1); prevQ = Fraction(n); ok_id = ok_C = ok_inc = ok_incr = ok_val = True
        recips = Fraction(0)
        for l in range(1, L + 1):
            S = 3 * S + 2 ** ds[l - 1]
            if 2 ** ds[l] * ms[l] != 3 ** l * n + S: ok_id = False
            C *= (1 + Fraction(1, 3 * ms[l - 1]))
            Q = Fraction(2 ** ds[l] * ms[l], 3 ** l)
            if Q != n * C: ok_C = False
            if Q - prevQ != prevQ / (3 * ms[l - 1]): ok_inc = False
            if not Q > prevQ: ok_incr = False
            num, den = Q.numerator, Q.denominator
            if v2(num) != ds[l] or den != 3 ** l or ms[l] % 3 == 0: ok_val = False
            prevQ = Q
            recips += Fraction(1, ms[l - 1])
        logC = math.log(float(C))
        ratio = float(recips) / logC
        mmin = min(ms[:L])
        out("   n=%6d (%3d odd steps): identity 2^d_l m_l = 3^l n + S_(l-1): %s; Q_l = n C_l: %s; Q_(l+1)-Q_l = Q_l/(3 m_l): %s; Q_l increasing: %s;" % (
            n, L, ok_id, ok_C, ok_inc, ok_incr))
        out("        |Q_l|_2 = 2^-d_l and |Q_l|_3 = 3^l (i.e. 3 does not divide m_l) for 1 <= l <= %d: %s; at l = 0, |Q_0|_3 = |n|_3 = 3^-%d" % (L, ok_val, 0 if n % 3 else 3))
        out("        sum_(j<L) 1/m_j = %.6f; log(Q_L/n) = %.6f; ratio = %.4f (note claims ratio <= 3; correct: 3 < ratio <= 3 + 1/min m_j = %.4f)" % (
            float(recips), logC, ratio, 3 + 1 / mmin))
    # the reciprocal constant over many orbits
    worst_lo = 10; worst_hi = 0
    for n in range(3, 4001, 2):
        ms, ds, vs = syracuse(n, 10000)
        L = len(ms) - 1
        s = sum(1 / m for m in ms[:L]); logC = sum(math.log(1 + 1 / (3 * m)) for m in ms[:L])
        r = s / logC
        worst_lo = min(worst_lo, r); worst_hi = max(worst_hi, r)
    out("   over all odd n <= 4000 (full orbits to 1): sum 1/m_j / log(Q_L/n) in [%.4f, %.4f]  -- always > 3 (log(1+x) < x), never <= 3" % (worst_lo, worst_hi))
    out("   correct inequality: x/(1+x) <= log(1+x) <= x with x = 1/(3 m) <= 1/3 gives  3 log(Q_inf/n) <= sum 1/m_l <= (3 + 1/min_l m_l) log(Q_inf/n) <= 4 log(Q_inf/n).")
    out("   The note's 'sum 1/m_l <= 3 log(Q_inf/n)' has the inequality reversed; the universal constant in the stated direction is 4. ERROR (constant), harmless for the use made of it.")
    # 3 does not divide m_l for l >= 1: enumeration
    bad = 0
    for n in range(1, 100001, 2):
        ms, ds, vs = syracuse(n, 10000)
        if any(m % 3 == 0 for m in ms[1:]): bad += 1
    out("   3 | m_l for some l >= 1 on the orbit of odd n <= 10^5: %d orbits (algebra: m_l = 2^-v (3 m_(l-1) + 1) = 2^-v mod 3 != 0)" % bad)


# ----------------------------------------------------------------------------------------------------
def section6():
    out("=" * 100)
    out("6. Proposition 3 (the past is the 3-adic address)")
    for n in (27, 703, 26623):
        ms, ds, vs = syracuse(n, 500)
        L = len(ms) - 1
        S = 0; equiv_ok = True; cong_all = True; resid_ok = True; first = None; tested = 0
        for l in range(1, L + 1):
            S = 3 * S + 2 ** ds[l - 1]
            Q = Fraction(2 ** ds[l] * ms[l], 3 ** l)
            c1 = (2 ** ds[l] > Q); c2 = (ms[l] < 3 ** l)
            if c1 != c2: equiv_ok = False
            mod = 3 ** l
            if (2 ** ds[l] * ms[l] - S) % mod != 0: cong_all = False
            if c1:
                tested += 1
                if first is None: first = l
                r = (pow(2, -ds[l], mod) * S) % mod
                if r != ms[l]: resid_ok = False
        out("   n=%6d: (2^d_l > Q_l) <=> (m_l < 3^l) for all l: %s; 2^d_l m_l = S_(l-1) mod 3^l for all l: %s; m_l = least positive residue of 2^-d_l S_(l-1) mod 3^l for all %d indices l >= %s with 2^d_l > Q_l: %s" % (
            n, equiv_ok, cong_all, tested, first, resid_ok))
    # locality: m_l mod 3^j depends only on the last j valuations
    out("   locality test: over the orbits of all odd n <= 3000, key = (v_(l-j),...,v_(l-1)) -> set of residues m_l mod 3^j:")
    orbits = {n: syracuse(n, 10000) for n in range(3, 3001, 2)}
    for j in range(1, 7):
        seen = defaultdict(set); formula_ok = True
        for n, (ms, ds, vs) in orbits.items():
            for l in range(j, len(ms)):
                key = tuple(vs[l - j:l]); mod = 3 ** j
                seen[key].add(ms[l] % mod)
                # closed form: m_l = sum_{t=1}^{j} 3^(t-1) 2^(-(v_(l-t)+...+v_(l-1))) mod 3^j
                acc = 0; s = 0
                for t in range(1, j + 1):
                    s += vs[l - t]
                    acc += 3 ** (t - 1) * pow(2, -s, mod)
                if acc % mod != ms[l] % mod: formula_ok = False
        multi = sum(1 for k_, v_ in seen.items() if len(v_) > 1)
        out("      j=%d: %5d distinct valuation windows, windows with more than one residue: %d; closed form m_l = sum_t 3^(t-1) 2^-(v_(l-t)+..+v_(l-1)) mod 3^j: %s" % (
            j, len(seen), multi, formula_ok))
    # explicit examples
    ex = []
    ms1, ds1, vs1 = orbits[27]; ms2, ds2, vs2 = orbits[703]
    for l1 in range(4, len(ms1)):
        for l2 in range(4, len(ms2)):
            if vs1[l1 - 3:l1] == vs2[l2 - 3:l2] and vs1[l1 - 4] != vs2[l2 - 4] and ms1[l1] % 81 != ms2[l2] % 81:
                ex.append((l1, l2)); break
        if ex: break
    if ex:
        l1, l2 = ex[0]
        out("   example: orbit 27 at l=%d (m=%d, last 4 valuations %s) and orbit 703 at l=%d (m=%d, last 4 valuations %s): m mod 27 = %d, %d (equal); mod 81 = %d, %d (differ, 4th-last valuation differs)" % (
            l1, ms1[l1], vs1[l1 - 4:l1], l2, ms2[l2], vs2[l2 - 4:l2], ms1[l1] % 27, ms2[l2] % 27, ms1[l1] % 81, ms2[l2] % 81))
    out("   VERDICT: Prop. 3 CONFIRMED (reduce the identity mod 3^l; m_l < 3^l iff Q_l < 2^d_l is immediate from m_l = 3^l Q_l / 2^d_l).")


# ----------------------------------------------------------------------------------------------------
def M_k(k, pred):
    """#{w in {0,1}^k : pred(o_i, i) for all 1 <= i <= k}; pred is the barrier condition."""
    dp = {0: 1}
    for i in range(1, k + 1):
        nd = defaultdict(int)
        for o, c in dp.items():
            for step in (0, 1):
                o2 = o + step
                if pred(o2, i): nd[o2] += c
        dp = nd
    return sum(dp.values())


def section7():
    out("=" * 100)
    out("7. Exits lemma and crossings (note section 3.4)")
    # X_0(b): |b| (3/2)^k <= X/4 with k = floor(log_2 X)
    for b in (1, 3, 5):
        fails = [X for X in range(2, 200000) if b * 1.5 ** math.floor(math.log2(X)) > X / 4]
        fails_sharp = [X for X in range(2, 200000) if b * (1.5 ** math.floor(math.log2(X)) - 1) > X / 4]
        out("   |b|=%d: |b|(3/2)^floor(log_2 X) <= X/4 fails exactly for X in [2, %d] (%d values), holds for all X >= %d; with the sharper carry |b|((3/2)^k - 1): holds for X >= %d" % (
            b, max(fails), len(fails), max(fails) + 1, max(fails_sharp) + 1))
    out("   (at X = 2^t the ratio is (3/4)^t |b|, decreasing in t, so the threshold is final: X_0(1) = 21; the 'X^0.585 <= X/4' reading gives X >= 4^(1/0.415) = 28.3.)")
    out("   log_2(3/5) = %.4f > -1.6: the derived reversed partial sums are > log_2(3/5), so M_k(0.737) already suffices; M_k(1.6) is a valid but looser count." % math.log2(0.6))
    pred_35 = lambda o, i: 5 * 3 ** o > 3 * 2 ** i          # 3^o/2^i > 3/5 exactly
    pred_16 = lambda o, i: o * ALPHA - i > -1.6
    out("   M_k(log_2(5/3)) and M_k(1.6) for k = 10, 14, 18, 22: %s" % [(k, M_k(k, pred_35), M_k(k, pred_16)) for k in (10, 14, 18, 22)])
    for name, n in (("63728127", 63728127), ("2^60-1", 2 ** 60 - 1)):
        xs = [n]
        while xs[-1] != 1: xs.append(T(xs[-1]))
        out("   T-orbit of %s: %d iterates to 1" % (name, len(xs) - 1))
        for t in (10, 14, 18, 22):
            X = 2 ** t; k = t
            stays = []; i = 0
            while i < len(xs):
                if xs[i] <= X:
                    j = i
                    while j + 1 < len(xs) and xs[j + 1] <= X: j += 1
                    if j + 1 < len(xs): stays.append((i, j))      # has an exit at j+1
                    i = j + 1
                else: i += 1
            long = [(i, j) for (i, j) in stays if j - i + 1 >= k]
            cond_ok = True; zs = set()
            for (i, j) in long:
                e = j + 1; zi = e - k; z = xs[zi]; zs.add(z)
                word = [xs[zi + s] % 2 for s in range(k)]
                o = [0]
                for s in range(k): o.append(o[-1] + word[s])
                for s in range(k):
                    if not (5 * 3 ** (o[k] - o[s]) > 3 * 2 ** (k - s)): cond_ok = False
                if any(xs[zi + s] > X for s in range(k)): cond_ok = False
            out("      X = 2^%d: stays with an exit: %d, of length >= k=%d: %d (bound 2 M_k(log2 5/3) = %d, 2 M_k(1.6) = %d); every long stay's z satisfies M_k/M_s > 3/5 for all s < k: %s; distinct z: %s" % (
                t, len(stays), k, len(long), 2 * M_k(k, pred_35), 2 * M_k(k, pred_16), cond_ok, len(zs) == len(long)))
    out("   Crossing counts on the odd iterates (level 2^t: #odd iterates <= 2^t / #downcrossings m_(l-1) > 2^t >= m_l):")
    for name, n in (("63728127", 63728127), ("2^60-1", 2 ** 60 - 1)):
        ms, ds, vs = syracuse(n, 5000)
        row = []
        for t in (10, 14, 18, 22, 26, 30):
            Y = 2 ** t
            below = sum(1 for m in ms if m <= Y)
            down = sum(1 for l in range(1, len(ms)) if ms[l - 1] > Y >= ms[l])
            row.append("2^%d: %d/%d" % (t, below, down))
        out("      %-10s (%d odd iterates incl. start) " % (name, len(ms)) + "  ".join(row))
    # the script's 5n+1 band pairs use ALPHA = log2 3 -- check with the right alpha
    ms, ds, vs = syracuse(7, 2000, 5)
    Y = 2 ** 18; vis = [l for l, m in enumerate(ms) if Y // 2 < m <= Y]
    pairs = [(vis[i + 1] - vis[i], ds[vis[i + 1]] - ds[vis[i]]) for i in range(len(vis) - 1)]
    out("   5n+1 band (2^17, 2^18] pairs (M, D) = %s: |M log_2 5 - D| = %s (the session script printed |M log_2 3 - D| = 2.25, 6.74 for these: wrong alpha for 5n+1; harmless, not quoted in the note)" % (
        pairs, ["%.2f" % abs(M * math.log2(5) - D) for M, D in pairs]))


# ----------------------------------------------------------------------------------------------------
def section8():
    out("=" * 100)
    out("8. The binary tree on odd numbers (note section 3.6)")
    def parent(m):
        if m % 8 == 5: return (m - 1) // 4
        x = 3 * m + 1; v = 0
        while x % 2 == 0: x //= 2; v += 1
        return x, v
    def children(m):
        ch = [4 * m + 1]
        if m % 3 == 1: c = (4 * m - 1) // 3
        elif m % 3 == 2: c = (2 * m - 1) // 3
        else: c = None
        if c is not None and c != m: ch.append(c)
        return ch
    LIM = 100_000
    alg_ok = all((v2(3 * m + 1) >= 3) == (m % 8 == 5) for m in range(1, 4 * LIM + 2, 2))
    out("   v_2(3m+1) >= 3 iff m = 5 mod 8, all odd m <= %d: %s (3m+1 = 0 mod 8 iff m = 7 * 3^-1 = 21 = 5 mod 8)" % (4 * LIM + 1, alg_ok))
    # parent well defined
    par = {}
    pv_ok = True
    for m in range(3, 4 * LIM + 2, 2):
        p = parent(m)
        if isinstance(p, tuple):
            x, v = p
            if v not in (1, 2): pv_ok = False
            par[m] = x
        else:
            par[m] = p
        if par[m] % 2 == 0 or par[m] < 1: pv_ok = False
    out("   parent(m) odd and >= 1 with v in {1,2} in the U-case for all odd 3 <= m <= %d: %s" % (4 * LIM + 1, pv_ok))
    # children consistency
    inv = defaultdict(set)
    for m, p in par.items(): inv[p].add(m)
    ch_ok = True; count_ok = True
    for m in range(1, LIM + 1, 2):
        pred = set(children(m))
        if pred != inv[m]: ch_ok = False
        expect = 1 if (m == 1 or m % 3 == 0) else 2
        if len(pred) != expect: count_ok = False
    out("   predicted children {4m+1} u {(2^v0 m - 1)/3 : 3 !| m, v0 least admissible} == {m' <= %d : parent(m') = m} for all odd m <= %d: %s" % (4 * LIM + 1, LIM, ch_ok))
    out("   child counts: 2 for 3 !| m, m > 1; 1 for 3 | m; the root 1 has one child (5; its 'root preimage' (4-1)/3 = 1 is itself): %s" % count_ok)
    # parent chains reach 1
    reach = {1: 0}
    def chain_len(m):
        path = []
        while m not in reach:
            path.append(m)
            m = par.get(m)
            if m is None:
                p = parent(path[-1]); m = p[0] if isinstance(p, tuple) else p
            if m in path: return None
        s = reach[m]
        for x in reversed(path):
            s += 1; reach[x] = s
        return reach[path[0]] if path else 0
    cyc = 0; mx = 0
    for m in range(1, LIM + 1, 2):
        r = chain_len(m)
        if r is None: cyc += 1
        else: mx = max(mx, r)
    out("   parent chains from every odd m <= %d reach 1 (no cycle): %s (cycles found: %d; longest chain %d edges); with unique parents this is a tree rooted at 1" % (LIM, cyc == 0, cyc, mx))
    out("   remark: U(4m+1) = U(m) (3(4m+1)+1 = 4(3m+1)), so the parent chain follows the Syracuse orbit with comb descents inserted;")
    out("   'the graph is a tree rooted at 1' is exactly the Collatz conjecture (no divergence, no other cycle), as the note says.  CONFIRMED (marked SPECULATION for the height idea).")


# ----------------------------------------------------------------------------------------------------
def section9():
    out("=" * 100)
    out("9. Text checks on the note")
    with open(NOTE, encoding='utf-8') as fh:
        txt = fh.read()
    props = sorted(set(re.findall(r"Proposition (\d)", txt)))
    out("   'Proposition N' labels present in the note: %s; 'Proposition 4' present: %s; status line says 'Propositions 1-4' (sections 3.1-3.3 are Props 1-3; '6' is the cited foundry Prop. 6)." % (props, '4' in props))
    out("   'PROVED' labels: %d occurrences; 'SPECULATION' markers: %d; section 2.4 header: %s; section 3.6 header: %s" % (
        txt.count("PROVED"), txt.count("SPECULATION"), "marked" if "2.4 Speculation on adjacent constructions (marked as such)" in txt else "NOT marked",
        "marked" if "3.6 A proof-idea attempt through the Pythagorean tree (SPECULATION)" in txt else "NOT marked"))
    out("   'Benford' occurrences: %d; the note never names a base; the proof gives base 2 only (section 4 above)." % txt.count("Benford"))
    out("   'for every k the best cell has delta within 0.006 of 1/2': contradicted by the session's own .out at k = 5, 7, 11, 13 (section 3 above).")


def main():
    out("collatz_crossings_20260926_audit.py -- independent audit of collatz_crossings_20260926_potential_and_seeds.md (2026-09-26)")
    section1(); section2(); section3(); section4(); section5(); section6(); section7(); section8(); section9()
    out("=" * 100)
    out("SUMMARY: 1.1 CONFIRMED (colouring, densities, doubling matrix explained by {2x}, MI values reproduced; owner's 35-term list not in repo);")
    out("  1.2 CONFIRMED; 2.1 CONFIRMED (f = THM-4503's m_j; bijection; width 2; W_k list); 2.2 CONFIRMED except the 'every k within 0.006' claim (ERROR: k=5,7,11,13);")
    out("  Prop 1 CONFIRMED in base 2, 'Benford' GAP (base 10 does not follow); Prop 2 CONFIRMED except the reciprocal-sum constant (ERROR: reversed, correct 4 / 3+1/min m);")
    out("  Prop 3 CONFIRMED; exits lemma CONFIRMED (X_0(1) = 21, 1.6 loose, needs distinct terms); crossings CONFIRMED; 3.6 tree CONFIRMED; 'Propositions 1-4' mislabel.")
    with open(OUT_PATH, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == '__main__':
    main()
