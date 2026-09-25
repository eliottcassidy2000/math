#!/usr/bin/env python3
"""procgen_continuous_20260925_mahler.py

Lane: natural boundaries, Mahler/Cobham and harmonic trees
(session collatz-procgen-20260922, 2026-09-25).

Part B: the Mahler / Cobham structure and its obstruction.

  MB1  Exact search for pure k-Mahler equations  sum_{i<=d} a_i(z) F(z^{k^i}) = 0,
       deg a_i <= delta, k in {2, 3}, over F_p (p = 2^31 - 1) on L = 3000 coefficients.
       Full rank => no such equation exists over Q (hence over C) with those (d, delta):
       a FINITE-EXACT non-existence statement.  Controls that must be FOUND:
       Thue-Morse (2-automatic), the base-3 Cantor set (3-automatic), an eventually periodic
       set and the Collatz basin truncation (all ones, rational).  Targets: the three 3n-1
       basins, and the 3-smooth numbers.
  MB2  k-kernel complexity (automaticity evidence): number of distinct length-W prefixes of
       n -> s(k^e n + r), 0 <= r < k^e, for e <= E.  Bounded for automatic sequences.
  MB3  The rotation obstruction (PROVED): the 3-smooth numbers S = {2^a 3^b} satisfy BOTH
       rotation-Mahler equations  F(z^2) = (F(z)+F(-z))/2  and  F(z^3) = (F(z)+F(wz)+F(w^2 z))/3,
       yet S is infinite of density zero, so F is not rational (natural boundary by Fabry).
       Hence Adamczewski-Bell / Schafke-Singer rigidity fails once roots of unity are allowed;
       the Collatz equation is of this rotation type.
  MB4  SHEET at the level of generating functions (checked): the T_+ basins on the negative
       integers, read in 1/z, are exactly the 3n-1 basins; so the two-sided 3n+1 equation
       splits into E_{3,+1} (for z) and E_{3,-1} (for 1/z), and the second has
       natural-boundary 0/1 solutions.

Every check raises on failure.  Peak memory about 250 MB; runtime about 1 min.
"""
import sys
import time
import math
import cmath
import random

import numpy as np

T0 = time.time()
P = 2**31 - 1


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)


def mem(tag):
    import resource
    print(f"[mem] {tag}: max RSS so far {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB",
          file=sys.stderr)


def hdr(s):
    mem("before " + s[:4])
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


def Tqb(n, q, b):
    return n // 2 if n % 2 == 0 else (q * n + b) // 2


def minus_basins(N):
    """Labels 1,2,3 for the 3n-1 basins of 1, 5, 17 on [1, N] (drop-below census, chunked)."""
    target = np.zeros(N + 1, dtype=np.int32)
    cyc = []
    for lo in range(1, N + 1, 1_000_000):
        idx = np.arange(lo, min(N, lo + 999_999) + 1, dtype=np.int64)
        val = idx.copy()
        steps = 0
        while idx.size:
            steps += 1
            check(steps < 20000, "minus census settles")
            odd = (val & 1).astype(bool)
            val[odd] = (3 * val[odd] - 1) // 2
            val[~odd] //= 2
            below = val < idx
            ret = val == idx
            target[idx[below]] = val[below].astype(np.int32)
            for n in idx[ret].tolist():
                cyc.append(n)
                target[n] = n
            keep = ~(below | ret)
            idx, val = idx[keep], val[keep]
    check(sorted(cyc) == [1, 5, 17], "3n-1 cycle minima")
    ptr = target
    for _ in range(40):
        nxt = ptr[ptr]
        if np.array_equal(nxt, ptr):
            break
        ptr = nxt
    lab = np.zeros(N + 1, dtype=np.int8)
    for i, m in enumerate((1, 5, 17)):
        lab[ptr == m] = i + 1
    return lab


def rank_mod_p(A):
    """Rank of an integer matrix modulo P (entries reduced mod P first)."""
    A = np.array(A, dtype=np.int64) % P
    m, n = A.shape
    r = 0
    for c in range(n):
        if r == m:
            break
        piv = np.nonzero(A[r:, c])[0]
        if piv.size == 0:
            continue
        i = r + int(piv[0])
        if i != r:
            A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), P - 2, P)
        A[r] = (A[r] * inv) % P
        rows = np.nonzero(A[:, c])[0]
        rows = rows[rows != r]
        if rows.size:
            f = A[rows, c].reshape(-1, 1)
            A[rows] = (A[rows] - (f * A[r]) % P) % P
        r += 1
    return r


def mahler_matrix(f, k, d, delta, L):
    """Rows m = 0..L-1 (coefficient of z^m); columns (i, t): coefficient of z^t in a_i.
    Entry = [z^(m-t)] F(z^(k^i)) = f[(m-t)/k^i] if k^i | (m-t) >= 0."""
    cols = []
    for i in range(d + 1):
        ki = k**i
        for t in range(delta + 1):
            col = np.zeros(L, dtype=np.int64)
            ms = np.arange(t, L)
            ok = ((ms - t) % ki) == 0
            src = (ms[ok] - t) // ki
            col[ms[ok]] = f[src]
            cols.append(col)
    return np.stack(cols, axis=1)


def has_mahler_eq(f, k, d, delta, L):
    A = mahler_matrix(f, k, d, delta, L)
    rk = rank_mod_p(A)
    return rk < A.shape[1], rk, A.shape[1]


def kernel_counts(s, k, E, W):
    out = []
    seen = set()
    for e in range(E + 1):
        ke = k**e
        check(ke * (W - 1) + ke - 1 < len(s), "kernel window inside the data")
        for r in range(ke):
            seen.add(s[r: r + ke * W: ke][:W].tobytes())
        out.append(len(seen))
    return out


def main():
    random.seed(20260925)
    N = 10**7
    t1 = time.time()
    lab = minus_basins(N)
    print(f"[3n-1 basins to 10^7 computed in {time.time()-t1:.1f}s]")

    L = 3000
    seqs = {}
    tm = np.array([bin(n).count("1") % 2 for n in range(L + 1)], dtype=np.int64)
    seqs["Thue-Morse (control, 2-automatic)"] = tm

    def cantor(n):
        while n:
            if n % 3 == 1:
                return 0
            n //= 3
        return 1
    seqs["base-3 Cantor set (control, 3-automatic)"] = np.array([cantor(n) if n else 1 for n in range(L + 1)],
                                                              dtype=np.int64)
    ep = np.array([1 if (n % 7 in (1, 2, 4) or n < 10 and n % 2 == 0) else 0 for n in range(L + 1)], dtype=np.int64)
    ep[0] = 0
    seqs["eventually periodic set (control, rational)"] = ep
    ones = np.ones(L + 1, dtype=np.int64)
    ones[0] = 0
    seqs["3n+1 basin(1) truncation = all ones (control, rational)"] = ones

    t1 = np.ones(L + 1, dtype=np.int64)
    t1[0] = 0
    m = 6
    while m <= L:
        t1[m] = 0
        m *= 2
    seqs["planted T1: basin = N minus {3*2^m} (DEFECT control)"] = t1

    def smooth3(n):
        if n == 0:
            return 0
        while n % 2 == 0:
            n //= 2
        while n % 3 == 0:
            n //= 3
        return 1 if n == 1 else 0
    sm = np.array([smooth3(n) for n in range(L + 1)], dtype=np.int64)
    seqs["3-smooth numbers {2^a 3^b}"] = sm
    for i, m in enumerate((1, 5, 17)):
        v = (lab[:L + 1] == i + 1).astype(np.int64)
        v[0] = 0
        seqs[f"3n-1 basin of {m}"] = v

    # -----------------------------------------------------------------------------------------
    hdr("MB1  Pure k-Mahler equations of small order d and degree delta (exact, mod p = 2^31-1)")
    grid = [(1, 40), (2, 30), (3, 24), (4, 19)]
    print(f"  L = {L} coefficients; (d, delta) grid {grid}; 'FOUND' = kernel nonzero mod p")
    results = {}
    for name, f in seqs.items():
        line = []
        for k in (2, 3):
            found = None
            for (d, delta) in grid:
                h, rk, nu = has_mahler_eq(f, k, d, delta, L)
                if h:
                    found = (d, delta, rk, nu)
                    break
            results[(name, k)] = found
            line.append(f"k={k}: " + (f"FOUND at (d,delta)=({found[0]},{found[1]}) rank {found[2]}/{found[3]}"
                                      if found else "none (full rank at every grid point)"))
        print(f"  {name:<52s} " + "; ".join(line))
    check(results[("Thue-Morse (control, 2-automatic)", 2)] is not None, "TM is 2-Mahler")
    check(results[("base-3 Cantor set (control, 3-automatic)", 3)] is not None, "Cantor is 3-Mahler")
    check(results[("eventually periodic set (control, rational)", 2)] is not None, "rational is 2-Mahler")
    check(results[("eventually periodic set (control, rational)", 3)] is not None, "rational is 3-Mahler")
    check(results[("3n+1 basin(1) truncation = all ones (control, rational)", 2)] is not None, "z/(1-z) 2-Mahler")
    for m in (1, 5, 17):
        for k in (2, 3):
            check(results[(f"3n-1 basin of {m}", k)] is None, f"no small {k}-Mahler equation for basin {m}")
    check(results[("planted T1: basin = N minus {3*2^m} (DEFECT control)", 2)] is not None, "B_T1 is 2-Mahler")
    check(results[("planted T1: basin = N minus {3*2^m} (DEFECT control)", 3)] is None, "no small 3-Mahler eq for B_T1")
    print("  => FINITE-EXACT: no 3n-1 basin series satisfies a 2- or 3-Mahler equation with")
    print("     (order, degree) in the grid; the solver finds the equations of every control.")
    print("  DEFECT: the planted map T1 has basin(1) = N minus {3*2^m : m>=1} (Part A, NB5); its series")
    print("  z/(1-z) - G(z), G(z) = z^6 + G(z^2), is 2-Mahler (PROVED) and irrational, while T1 has a")
    print("  divergent orbit.  So a single pure Mahler equation cannot imply Collatz by any argument that")
    print("  survives density-zero planting; the SECOND (3-Mahler) equation of Adamczewski-Bell is essential.")

    # -----------------------------------------------------------------------------------------
    hdr("MB2  k-kernel complexity (distinct length-W prefixes of n -> s(k^e n + r), cumulative over e)")
    W = 48
    ctl = {}
    Nk = 3_200_000
    tm_big = np.array([0], dtype=np.int8)
    # Thue-Morse to Nk via bit tricks
    x = np.arange(Nk, dtype=np.int64)
    pc = np.zeros(Nk, dtype=np.int8)
    y = x.copy()
    while y.any():
        pc ^= (y & 1).astype(np.int8)
        y >>= 1
    ctl["Thue-Morse"] = pc
    sm_big = np.zeros(Nk, dtype=np.int8)
    a = 1
    while a < Nk:
        b3 = a
        while b3 < Nk:
            sm_big[b3] = 1
            b3 *= 3
        a *= 2
    ctl["3-smooth"] = sm_big
    for i, m in enumerate((1, 5, 17)):
        ctl[f"3n-1 basin of {m}"] = (lab[:Nk] == i + 1).astype(np.int8)
    for k, E in ((2, 16), (3, 10)):
        print(f"  k = {k}, e = 0..{E}, W = {W}:")
        for name, s in ctl.items():
            cnt = kernel_counts(s, k, E, W)
            print(f"    {name:<22s} " + " ".join(str(c) for c in cnt))
            if name == "Thue-Morse" and k == 2:
                check(cnt[-1] == 2, "TM 2-kernel has 2 elements")
    # the Terras identity explains the k = 2 rows of the 3n-1 basins exactly
    seen = set()
    pairs = []
    for e in range(17):
        for r in range(2**e):
            x, a = r, 0
            for _ in range(e):
                if x % 2:
                    a += 1
                x = x // 2 if x % 2 == 0 else (3 * x - 1) // 2
            seen.add((a, x))
        pairs.append(len(seen))
    rows2 = {name: kernel_counts(s, 2, 16, W) for name, s in ctl.items() if name.startswith("3n-1")}
    for name, cnt in rows2.items():
        check(cnt == pairs, f"2-kernel count of {name} = number of distinct pairs (a_e(r), T^e(r))")
    print("  Terras identity (PROVED): for ANY T-invariant s, s(2^e n + r) = s(3^a n + T^e(r)), a = number of")
    print("  odd steps among the first e steps of r.  So the 2-kernel of an invariant set is a set of")
    print("  subsequences n -> s(3^a n + t); the k=2 rows of the three basins coincide EXACTLY with the")
    print(f"  cumulative number of distinct pairs (a, T^e(r)): {pairs}")
    print("  (every pair gives a different subsequence, for each basin, up to e = 16).  The statistic is")
    print("  type-level: it cannot tell the basins apart, and it is 1 for the Collatz basin (all ones).")
    print("  Reading: bounded rows are automatic; rows growing like k^e are consistent with")
    print("  non-automaticity (evidence only; automaticity of the 3n-1 basins is not decided here).")
    del x, pc, y

    # -----------------------------------------------------------------------------------------
    hdr("MB3  The rotation obstruction: 3-smooth numbers satisfy both rotation-Mahler equations")
    Ns = 10**6
    s = np.zeros(Ns + 1, dtype=np.int8)
    a = 1
    while a <= Ns:
        b3 = a
        while b3 <= Ns:
            s[b3] = 1
            b3 *= 3
        a *= 2
    n = np.arange(1, Ns // 2 + 1)
    check(np.array_equal(s[2 * n], s[n]), "s(2n) = s(n)")
    n = np.arange(1, Ns // 3 + 1)
    check(np.array_equal(s[3 * n], s[n]), "s(3n) = s(n)")
    cnt = int(s.sum())
    pred = (math.log(Ns) ** 2) / (2 * math.log(2) * math.log(3))
    print(f"  s(2n) = s(n), s(3n) = s(n) on [1, {Ns}] (exact); #S cap [1,{Ns}] = {cnt} "
          f"(~ (ln N)^2/(2 ln2 ln3) = {pred:.0f})")
    coef = s[:4001].astype(float)

    def F(z):
        acc = 0j
        for c in coef[::-1]:
            acc = acc * z + c
        return acc
    om = cmath.exp(2j * math.pi / 3)
    worst = 0.0
    for _ in range(8):
        z = random.uniform(0.2, 0.85) * cmath.exp(1j * random.uniform(0, 2 * math.pi))
        e2 = abs(F(z * z) - (F(z) + F(-z)) / 2)
        e3 = abs(F(z**3) - (F(z) + F(om * z) + F(om * om * z)) / 3)
        worst = max(worst, e2, e3)
    check(worst < 1e-10, "rotation-Mahler identities numerically")
    print(f"  numerically: max |F(z^2)-(F(z)+F(-z))/2|, |F(z^3)-mean_j F(w^j z)| over 8 points = {worst:.1e}")
    print("  S is infinite with density 0, hence not eventually periodic, so F is irrational; its exponents")
    print("  n_k satisfy n_k/k -> infinity (counting function O(log^2 N)), so |z| = 1 is a natural boundary")
    print("  (Fabry-Faber gap theorem).  Uncountably many more: S_A = {2^a 3^b m : m in A} for any A of")
    print("  integers prime to 6.  => The Adamczewski-Bell / Schafke-Singer theorem has no analogue for")
    print("  Mahler equations with roots of unity; the Collatz (Berg-Meinardus) equation is of that type.")

    # -----------------------------------------------------------------------------------------
    hdr("MB4  SHEET at the level of series: the negative half of 3n+1 is the 3n-1 problem")
    Nn = 10**5
    cyc_neg = {-1: 1, -5: 2, -7: 2, -10: 2}
    c17 = [-17]
    x = Tqb(-17, 3, 1)
    while x != -17:
        c17.append(x)
        x = Tqb(x, 3, 1)
    for v in c17:
        cyc_neg[v] = 3
    bad = 0
    for m in range(1, Nn + 1):
        x = -m
        for _ in range(100000):
            if x in cyc_neg:
                break
            x = Tqb(x, 3, 1)
        else:
            raise AssertionError("negative orbit did not settle")
        if cyc_neg[x] != int(lab[m]):
            bad += 1
    check(bad == 0, "T_+ on -m lands in the cycle of -c exactly when T_- on m lands in the cycle of c")
    print(f"  for every m <= {Nn}: the T_+-orbit of -m ends in the cycle of -1, -5 or -17 exactly when the")
    print("  T_- orbit of m ends in the cycle of 1, 5 or 17 (direct iteration of T_+ on negatives).")
    print(f"  T_+ cycle through -17: {c17}")
    print("  So for a T_+-invariant S in Z\\{0}: sum_{n in S, n>0} z^n solves E_{3,+1} and")
    print("  sum_{n in S, n<0} z^{-n} solves E_{3,-1}; the second equation has three natural-boundary")
    print("  0/1 solutions (the 3n-1 basins).  A proof that E_{3,+1} has only rational 0/1 solutions")
    print("  must therefore use what distinguishes b = +1 from b = -1 in")
    print("      h(z^3) - h(z^6) = z^{-b} (1/3) sum_j w^{bj} h(w^j z^2):")
    print("  the monomial z^{-b} and the sieved class n = -b mod 3, i.e. the sign.")

    print()
    print(f"ALL CHECKS PASSED  ({time.time() - T0:.1f}s)")
    mem("end")


if __name__ == "__main__":
    main()
