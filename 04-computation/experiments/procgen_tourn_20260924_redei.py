#!/usr/bin/env python3
"""procgen_tourn_20260924, part C: Redei as a fixed-point parity, Sperner, and the Collatz odd-count candidates.

Sections (every check raises on failure):
  C0  Redei through the OCF (THM-002, Grinberg-Stanley form): H(T) = sum over permutations sigma whose
      nontrivial cycles are directed odd cycles of T of 2^{#nontrivial cycles}.  On X = {(sigma, bits)}
      (|X| = H) the involution "flip the bit of the cycle through the least moved point" has exactly one
      fixed point, (id, ()): H = 1 (mod 2).  All tournaments n <= 5, random n = 6, 7.
  C1  Sperner's lemma in dimension 2 (the door-in/door-out parity): random Sperner labellings of a
      triangulated triangle; rainbow count odd; door graph has max degree 2.
  C2  Collatz odd-count candidate (a): depth-d preimage counts N_d(a) on Z/3^d (Haar model, inverse-tree
      note section 3.1): parity is not constant; the recursion increment is not even (contrast: the
      tournament vertex-addition increment H(T) - H(T-v) = 2 sum mu(C) is always even, THM-070).
      Second moments: E[N_d^2]/E[N_d]^2 (d <= 15) versus the tournament ratio W(n)/n! (THM-589) and the
      Galton-Watson value 1.5.
  C3  candidates (b), (c): 2-adic periodic points x_w = c_w/(2^p - 3^a) of T, p <= 16: exactly 2^p, paired
      by the FREE involution iota = Phi^-1 o complement o Phi (x_w <-> x_{w-bar}), so #Fix(T^p) is even;
      positive and integer ones; R_p(n) = #{w : x_w in [1, n]} has no parity law.
  C4  Hamiltonian cycles of the forward type graph of T_b mod 2^k: it is the de Bruijn graph B(2,k) for
      3n+1, 3n-1 and 5n+1 alike; #Ham cycles = 2^(2^(k-1) - k) (de Bruijn / Flye Sainte-Marie), checked
      by DP (k <= 4) and by the BEST / matrix-tree theorem (k <= 8).
Run: python3 04-computation/experiments/procgen_tourn_20260924_redei.py   (about 1 minute, < 400 MB)
"""
import itertools
import math
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


# ---------------------------------------------------------------- C0
def ham_brute(n, A):
    return sum(1 for p in itertools.permutations(range(n)) if all(A[p[i]][p[i + 1]] for i in range(n - 1)))


def cycles_of(perm):
    n = len(perm)
    seen = [False] * n
    cyc = []
    for s in range(n):
        if not seen[s]:
            c = [s]
            seen[s] = True
            t = perm[s]
            while t != s:
                c.append(t)
                seen[t] = True
                t = perm[t]
            cyc.append(c)
    return cyc


def ocf_configurations(n, A):
    """Permutations whose nontrivial cycles are directed odd cycles of A (v -> perm[v] is an arc)."""
    confs = []
    for perm in itertools.permutations(range(n)):
        cyc = [c for c in cycles_of(perm) if len(c) > 1]
        if all(len(c) % 2 == 1 and all(A[c[i]][c[(i + 1) % len(c)]] for i in range(len(c))) for c in cyc):
            confs.append((perm, cyc))
    return confs


def section_C0():
    print("=" * 100)
    print("C0  Redei's theorem as the fixed-point count of an involution (OCF / Grinberg-Stanley form)")
    rng = random.Random(20260924)
    tested = 0
    for n in range(1, 8):
        pairs = list(itertools.combinations(range(n), 2))
        if n <= 5:
            codes = range(1 << len(pairs))
        else:
            codes = [rng.getrandbits(len(pairs)) for _ in range(40 if n == 6 else 8)]
        for code in codes:
            A = [[0] * n for _ in range(n)]
            for k, (u, v) in enumerate(pairs):
                if (code >> k) & 1:
                    A[u][v] = 1
                else:
                    A[v][u] = 1
            H = ham_brute(n, A)
            confs = ocf_configurations(n, A)
            X = [(perm, bits) for perm, cyc in confs for bits in itertools.product((0, 1), repeat=len(cyc))]
            check(len(X) == H, "OCF (Grinberg-Stanley form): H = sum 2^psi (n=%d)" % n)
            cycmap = {perm: cyc for perm, cyc in confs}
            fixed = 0
            for perm, bits in X:
                cyc = cycmap[perm]
                if not cyc:
                    fixed += 1
                    continue
                j = min(range(len(cyc)), key=lambda i: min(cyc[i]))
                img = (perm, tuple(b ^ (1 if i == j else 0) for i, b in enumerate(bits)))
                check(img != (perm, bits), "fixed-point-free off the identity")
            check(fixed == 1, "exactly one fixed point, (id, ())")
            check(H % 2 == 1, "Redei")
            tested += 1
    print("  %d tournaments (all with n <= 5, random n = 6, 7): H = #X, X = {(sigma, bits)}; the involution" % tested)
    print("  flipping the bit of the cycle through the least moved point is fixed-point-free except at (id, ()).")
    print("  So H = 1 (mod 2): Redei = 'the identity is the only permutation without a nontrivial cycle'.")
    print("  (Higher digits: H = sum_{k<m} alpha_k 2^k mod 2^m, THM-466; checked on windows in part B.)")


# ---------------------------------------------------------------- C1
def section_C1():
    print("=" * 100)
    print("C1  Sperner's lemma (dimension 2): odd number of rainbow cells via the door graph (handshake)")
    rng = random.Random(7)
    N = 14
    pts = [(i, j) for i in range(N + 1) for j in range(N + 1 - i)]
    tris = []
    for i in range(N):
        for j in range(N - i):
            tris.append(((i, j), (i + 1, j), (i, j + 1)))
            if i + j + 1 < N:
                tris.append(((i + 1, j), (i + 1, j + 1), (i, j + 1)))
    counts = {}
    for trial in range(400):
        lab = {}
        for (i, j) in pts:
            k = N - i - j  # barycentric (i, j, k); corners (N,0,0)->0, (0,N,0)->1, (0,0,N)->2
            allowed = [c for c, coord in enumerate((i, j, k)) if coord > 0]
            lab[(i, j)] = rng.choice(allowed)
        rainbow = sum(1 for t in tris if {lab[p] for p in t} == {0, 1, 2})
        # door graph: doors = edges labelled {0,1}; degree of each cell = number of its {0,1}-edges
        deg_ok = True
        for t in tris:
            doors = sum(1 for a, b in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2])) if {lab[a], lab[b]} == {0, 1})
            deg_ok &= doors <= 2
            if {lab[p] for p in t} == {0, 1, 2}:
                deg_ok &= doors == 1
        side01 = [(N - s, s) for s in range(N + 1)]  # the side k = 0, from corner 0 to corner 1
        bdoors = sum(1 for u, v in zip(side01, side01[1:]) if {lab[u], lab[v]} == {0, 1})
        check(deg_ok, "door degrees")
        check(rainbow % 2 == 1 and bdoors % 2 == 1, "Sperner parity (rainbow and boundary doors odd)")
        counts[rainbow] = counts.get(rainbow, 0) + 1
    print("  400 random Sperner labellings of a %d-subdivided triangle: rainbow counts %s" % (N, dict(sorted(counts.items()))))
    print("  every cell has <= 2 doors (edges labelled {0,1}), rainbow cells exactly 1: paths pair the doors,")
    print("  the odd number of boundary doors on the 0-1 side forces an odd number of rainbow cells.")


# ---------------------------------------------------------------- C2
def W_odd_compositions(nmax):
    """THM-589: W(n) = sum_k k! [x^n] g(x)^k, g = x(1+x^2)/(1-x^2) = x + 2x^3 + 2x^5 + ..."""
    g = [0] * (nmax + 1)
    for e in range(1, nmax + 1, 2):
        g[e] = 1 if e == 1 else 2
    W = [0] * (nmax + 1)
    power = [1] + [0] * nmax  # g^0
    for k in range(1, nmax + 1):
        new = [0] * (nmax + 1)
        for i, c in enumerate(power):
            if c:
                for j in range(1, nmax + 1 - i):
                    if g[j]:
                        new[i + j] += c * g[j]
        power = new
        f = math.factorial(k)
        for n in range(nmax + 1):
            W[n] += f * power[n]
    return W


def section_C2():
    print("=" * 100)
    print("C2  Collatz candidate (a): preimage counts N_d(a) = |T^-d(a)| as functions on Z/3^d")
    dmax = 15
    t0 = time.time()
    rows = []
    # independent 3-type branching model: type 0 -> one type-0 child; type 1 -> one type-2 child;
    # type 2 -> one type-1 child + one child of uniform type (subtrees independent)
    m = {0: Fr(1), 1: Fr(1), 2: Fr(1)}
    s2 = {0: Fr(1), 1: Fr(1), 2: Fr(1)}
    model = {}
    for d in range(1, dmax + 1):
        mu = (m[0] + m[1] + m[2]) / 3
        su = (s2[0] + s2[1] + s2[2]) / 3
        nm = {0: Fr(1), 1: m[2], 2: m[1] + mu}
        ns = {0: Fr(1), 1: s2[2], 2: s2[1] + su + 2 * m[1] * mu}
        m, s2 = nm, ns
        allm, alls = sum(m.values()) / 3, sum(s2.values()) / 3
        um, us = (m[1] + m[2]) / 2, (s2[1] + s2[2]) / 2
        model[d] = (float(alls / allm ** 2), float(us / um ** 2))
    prev = np.ones(1, dtype=np.int32)
    for d in range(1, dmax + 1):
        M, Mp = 3 ** d, 3 ** (d - 1)
        a = np.arange(M, dtype=np.int32)
        base = prev[(2 * a) % Mp]
        leg = (a % 3) == 2
        incr = np.zeros(M, dtype=np.int32)
        incr[leg] = prev[((2 * a[leg] - 1) // 3) % Mp]
        Nd = base + incr
        check(int(Nd.max()) < 46000, "N_d^2 fits in int32")
        mean = Fr(int(Nd.sum(dtype=np.int64)), M)
        check(mean == Fr(4, 3) ** d, "E[N_d] = (4/3)^d")
        m2 = Fr(int((Nd * Nd).sum(dtype=np.int64)), M)
        unit = (a % 3) != 0
        Nu = Nd[unit]
        um = Fr(int(Nu.sum(dtype=np.int64)), int(Nu.size))
        u2 = Fr(int((Nu * Nu).sum(dtype=np.int64)), int(Nu.size))
        odd = float(np.mean(Nd % 2))
        incr_odd = float(np.mean(incr % 2))
        rows.append((d, float(mean), float(m2 / mean ** 2), float(um), float(u2 / um ** 2), odd, incr_odd,
                     int(Nu.min()), int(Nd.max())))
        prev = Nd
        del a, leg, base, incr, unit, Nu
    print("   d   E[N_d]  ratio(all)  model(all)   unit mean  ratio(units)  model(units)  P(N_d odd)  P(incr odd)  min_units  max")
    for d, mean, ratio, umean, uratio, odd, iodd, mn, mx in rows:
        print("  %2d %8.3f %10.5f %10.5f %11.3f %12.5f %12.5f %11.4f %11.4f %9d %6d"
              % (d, mean, ratio, model[d][0], umean, uratio, model[d][1], odd, iodd, mn, mx))
    check(abs(rows[3][4] - 1.070) < 0.001 and abs(rows[13][4] - 1.089) < 0.001,
          "units ratio reproduces the inverse-tree note table (1.070 at d=4, 1.089 at d=14)")
    check(all(0.2 < r[5] < 0.8 for r in rows[3:]), "N_d parity is not constant")
    check(all(r[6] > 0.05 for r in rows[1:]), "the increment is not even")
    print("  ratio(units) reproduces the inverse-tree note's table (1.070, 1.078, 1.088, 1.089 at d = 4, 8, 12, 14)")
    print("  and reaches %.4f at d = 15; the model with independent subtrees gives %.4f: the 3-adic digits shared by"
          % (rows[-1][4], model[dmax][1]))
    print("  the D- and E-subtrees make the true counts MORE concentrated, but the ratio does not decrease toward 1.")
    print("  Parity: N_d(a) is odd on a fraction ~%.2f of classes and the increment [a = 2 mod 3] N_{d-1}(E(a))" % rows[-1][5])
    print("  is odd on a fraction ~%.2f: NO Redei-type parity law, and no even-increment law (contrast: the" % rows[-1][6])
    print("  tournament vertex-addition increment H(T) - H(T-v) = 2 sum_(C ni v) mu(C) is always even, THM-070).")
    W = W_odd_compositions(40)
    check(W[1:9] == [1, 2, 8, 32, 158, 928, 6350, 49752], "THM-589 values W(1..8)")
    print("  Tournaments (THM-589): E[H^2]/E[H]^2 = W(n)/n! = %s" % ", ".join(
        "%d:%.4f" % (n, W[n] / math.factorial(n)) for n in (4, 6, 8, 12, 16, 24, 32, 40)))
    check(abs(W[40] / math.factorial(40) - 1) < 0.06, "W(n)/n! -> 1")
    print("  => tournament H is self-averaging (ratio -> 1, like 1 + 2/n); the Collatz preimage count is not")
    print("     (units ratio oscillates near 1.09 for d = 9..15, finite-exact).  E[H] = n!/2^(n-1) grows")
    print("     factorially, E[N_d] = (4/3)^d exponentially.  'N_d behaves like H' is false in growth and in")
    print("     fluctuation.  (Existence and value of the limit of the Haar ratio: OPEN, not claimed.)")
    print("  (time %.1f s)" % (time.time() - t0), file=sys.stderr)


# ---------------------------------------------------------------- C3
def x_word(w):
    """Rational periodic point of T with parity word w (length p): x = c_w/(2^p - 3^a)."""
    # compose the affine branches x -> x/2 (e=0), x -> (3x+1)/2 (e=1): T^p(x) = A x + B on the cylinder
    A, B = Fr(1), Fr(0)
    for e in w:
        if e:
            A, B = 3 * A / 2, (3 * B + 1) / 2
        else:
            A, B = A / 2, B / 2
    return B / (1 - A)


def T_rat(x):
    check(x.denominator % 2 == 1, "odd denominator")
    return x / 2 if x.numerator % 2 == 0 else (3 * x + 1) / 2


def section_C3():
    print("=" * 100)
    print("C3  Collatz candidates (b), (c): the 2-adic periodic points of T (the fixed points of T^p)")
    print("   p   #words  #distinct  T^p-fixed  iota free  #x_w>0  #x_w in Z_{>0}  R_p(1) R_p(2) R_p(10) R_p(100)")
    parities = {k: set() for k in ("pos", "R1", "R2", "R10", "R100")}
    int_points = {}
    for p in range(1, 17):
        pts = {}
        for w in itertools.product((0, 1), repeat=p):
            x = x_word(w)
            pts[w] = x
        vals = list(pts.values())
        check(len(set(vals)) == 2 ** p, "2^p distinct periodic points")
        # verify T^p(x_w) = x_w and the parity word (sample all for p <= 12, else every 7th)
        ws = list(pts)
        for w in (ws if p <= 12 else ws[::7]):
            x = pts[w]
            y = x
            for e in w:
                check(y.numerator % 2 == e, "parity word")
                y = T_rat(y)
            check(y == x, "T^p(x_w) = x_w")
        # the free involution: complement word <-> iota
        for w in ws:
            wb = tuple(1 - e for e in w)
            check(pts[wb] != pts[w], "iota has no fixed point")
        pos = sum(1 for x in vals if x > 0)
        ints = sorted({int(x) for x in vals if x.denominator == 1})
        int_points[p] = ints
        R = {n: sum(1 for x in vals if 1 <= x <= n) for n in (1, 2, 10, 100)}
        parities["pos"].add(pos % 2)
        for n, key in ((1, "R1"), (2, "R2"), (10, "R10"), (100, "R100")):
            parities[key].add(R[n] % 2)
        print("  %2d %7d %9d %10s %10s %7d %14d %7d %6d %7d %8d" % (
            p, 2 ** p, len(set(vals)), "yes", "yes", pos, sum(1 for x in ints if x > 0), R[1], R[2], R[10], R[100]))
    print("  integer periodic points (all periods p <= 16): %s" % sorted({x for v in int_points.values() for x in v}))
    check(sorted({x for v in int_points.values() for x in v}) ==
          sorted([0, -1, 1, 2, -5, -7, -10, -17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34]),
          "integer periodic points = the five known cycles")
    print("  parities seen: #positive %s, R_p(1) %s, R_p(2) %s, R_p(10) %s, R_p(100) %s" % tuple(
        sorted(parities[k]) for k in ("pos", "R1", "R2", "R10", "R100")))
    check(parities["R10"] == {0, 1} and parities["R100"] == {0, 1}, "no parity law for R_p(n)")
    print("  => #Fix(T^p) on Z_2 is 2^p: every periodic point is rational, and iota = Phi^-1 o (bit complement) o Phi")
    print("     (Phi the parity-vector isometry, which conjugates T to the shift) is a FREE involution commuting")
    print("     with T.  Even count with no fixed point: the opposite of Redei's single fixed configuration.")
    print("     Counts of cycle candidates in [1, n] take both parities: no Redei-type law.")


# ---------------------------------------------------------------- C4
def type_graph(k, q=3, b=1):
    """Forward type graph mod 2^k: r -> T(r') mod 2^k for the two lifts r' of r mod 2^(k+1)."""
    M = 2 ** k
    edges = []
    for r in range(M):
        for lift in (r, r + M):
            y = lift // 2 if lift % 2 == 0 else (q * lift + b) // 2
            edges.append((r, y % M))
    return edges


def word_of(r, k, q=3, b=1):
    w = []
    x = r
    for _ in range(k):
        w.append(x % 2)
        x = x // 2 if x % 2 == 0 else (q * x + b) // 2
    return tuple(w)


def bareiss_det(Mx):
    M = [row[:] for row in Mx]
    n = len(M)
    sign, prev = 1, 1
    for i in range(n - 1):
        if M[i][i] == 0:
            sw = next((r for r in range(i + 1, n) if M[r][i] != 0), None)
            if sw is None:
                return 0
            M[i], M[sw] = M[sw], M[i]
            sign = -sign
        for r in range(i + 1, n):
            for c in range(i + 1, n):
                M[r][c] = (M[r][c] * M[i][i] - M[r][i] * M[i][c]) // prev
            M[r][i] = 0
        prev = M[i][i]
    return sign * M[n - 1][n - 1]


def ham_cycles_dp(nv, edges):
    out = [0] * nv
    for u, v in edges:
        if u != v:
            out[u] |= 1 << v
    # count Hamiltonian cycles through vertex 0 (each directed cycle once)
    dp = {(1, 0): 1}
    for _ in range(nv - 1):
        nxt = {}
        for (mask, v), c in dp.items():
            o = out[v]
            while o:
                u = (o & -o).bit_length() - 1
                o &= o - 1
                if not (mask >> u) & 1:
                    key = (mask | (1 << u), u)
                    nxt[key] = nxt.get(key, 0) + c
        dp = nxt
    return sum(c for (mask, v), c in dp.items() if (out[v] & 1))


def section_C4():
    print("=" * 100)
    print("C4  Hamiltonian cycles ('chains') of the forward type graph of the map mod 2^k")
    maps = {"3n+1": (3, 1), "3n-1": (3, -1), "5n+1": (5, 1)}
    print("   k   de Bruijn iso (3n+1, 3n-1, 5n+1)   #Ham cycles (DP, k<=4)   #Eulerian circuits of level k-1 (BEST)   2^(2^(k-1)-k)")
    for k in range(1, 9):
        iso = []
        for name, (q, b) in maps.items():
            E = type_graph(k, q, b)
            phi = {r: word_of(r, k, q, b) for r in range(2 ** k)}
            check(len(set(phi.values())) == 2 ** k, "Terras bijection mod 2^k")
            ok = sorted((phi[u], phi[v]) for u, v in E) == sorted(
                (w, w[1:] + (e,)) for w in itertools.product((0, 1), repeat=k) for e in (0, 1))
            iso.append(ok)
        check(all(iso), "type graph = de Bruijn graph B(2,k)")
        target = 2 ** (2 ** (k - 1) - k)
        dp = None
        if k <= 4:
            counts = {ham_cycles_dp(2 ** k, type_graph(k, q, b)) for (q, b) in maps.values()}
            check(len(counts) == 1, "same Ham-cycle count for every map")
            dp = counts.pop()
            check(dp == target, "Ham cycles = 2^(2^(k-1)-k)")
        best = None
        if k >= 2:
            vals = set()
            for (q, b) in maps.values():
                E = type_graph(k - 1, q, b)
                n = 2 ** (k - 1)
                L = [[0] * n for _ in range(n)]
                for u, v in E:
                    if u != v:
                        L[u][u] += 1
                        L[u][v] -= 1
                red = [row[1:] for row in L[1:]]
                vals.add(bareiss_det(red))  # arborescences; BEST with all out-degrees 2: circuits = t_w
            check(len(vals) == 1, "same BEST count for every map")
            best = vals.pop()
            check(best == target, "BEST count = 2^(2^(k-1)-k)")
        print("  %2d   %-34s %-24s %-40s %d" % (k, str(iso), str(dp), str(best), target))
    print("  => the Collatz type graph mod 2^k is the de Bruijn graph for 3n+1, 3n-1 and 5n+1 alike; its")
    print("     Hamiltonian cycles are the binary de Bruijn sequences of order k, 2^(2^(k-1)-k) of them (even")
    print("     for k >= 3).  A doubly exponential 'chain growth' that is sheet-blind and drift-blind.")


def main():
    t0 = time.time()
    section_C0()
    section_C1()
    section_C2()
    section_C3()
    section_C4()
    print("=" * 100)
    print("PART C: ALL CHECKS PASSED")
    print("(time %.0f s)" % (time.time() - t0), file=sys.stderr)


if __name__ == "__main__":
    main()
