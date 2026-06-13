#!/usr/bin/env python3
"""
c3_spectrum_dp_kpc1.py -- THREAD C (HYP-2368, kind-pasteur-2026-06-10-S1 sub-agent kpc1)

EXACT spectrum of c3 (number of 3-cycles) over all n-vertex tournaments, for
n = 3..30 (and beyond while a time budget lasts), WITHOUT enumerating the
~A000571(n) Landau sequences individually (A000571(30) ~ 10^14 -- hopeless).

Method: c3(T) = C(n,3) - f(s),  f(s) = sum C(s_i,2)  (Kendall-Babington Smith /
Goodman identity, sanity-verified in c3_spectrum_kpc1.py), and Landau's theorem
(a nondecreasing integer sequence is a score sequence iff prefix sums >= C(k,2)
with equality at k = n).  So the achievable c3-set is determined by the set of
achievable f-values over Landau sequences.

DP over score VALUES v = 0..n-1 with multiplicities.  State = (c, w) where c =
number of scores chosen so far (all <= v) and w = their sum.  Value of a state =
bitset (Python int) of achievable partial f.  Transition: choose multiplicity t
of score v: (c,w) -> (c+t, w+t*v), f-bitset shifted by t*C(v,2).

Landau-condition correctness inside a block of t equal scores v: the prefix
condition at position k = c+j (0 <= j <= t) reads  w + j*v >= C(c+j, 2).
g(j) = w + j*v - C(c+j,2) is CONCAVE in j (linear minus convex), so its minimum
over 0 <= j <= t is at an endpoint; both endpoints are checked by the DP.
Hence checking the boundary states (c,w) and (c+t, w+t*v) suffices -- EXACT.
This is additionally cross-checked against the full backtracking enumeration of
c3_spectrum_kpc1.py for all n <= 13.

EXACT integer arithmetic throughout (Python bigints as bitsets). No floats.
"""

import time

def comb2(x):
    return x * (x - 1) // 2

def comb3(x):
    return x * (x - 1) * (x - 2) // 6

def max_c3_formula(n):
    if n % 2 == 1:
        return (n**3 - n) // 24
    return (n**3 - 4 * n) // 24

def near_regular_f(n):
    """f of the (near-)regular sequence = minimum of f over Landau sequences."""
    if n % 2 == 1:
        h = (n - 1) // 2
        return n * comb2(h)
    h = n // 2
    return h * (comb2(h) + comb2(h - 1))

def dp_f_bitset(n):
    """Bitset (int) U with bit q set iff some Landau sequence on n vertices has
    f(s) = sum C(s_i,2) = q."""
    N2 = comb2(n)
    layer = {(0, 0): 1}  # before any value processed: 0 scores, sum 0, f=0
    for v in range(n):
        cv = comb2(v)
        nxt = {}
        for (c, w), B in layer.items():
            maxt = n - c
            for t in range(maxt + 1):
                c2 = c + t
                w2 = w + t * v
                if w2 > N2:
                    break                       # w2 increases with t
                if w2 < comb2(c2):
                    continue                    # Landau boundary condition
                r = n - c2
                if r == 0:
                    if w2 != N2:
                        continue
                else:
                    # remaining r scores must lie in [v+1, n-1]
                    if w2 + r * (v + 1) > N2:
                        continue                # decreases with t: keep trying
                    if w2 + r * (n - 1) < N2:
                        break                   # non-increasing with t
                key = (c2, w2)
                nb = B << (t * cv)
                if key in nxt:
                    nxt[key] |= nb
                else:
                    nxt[key] = nb
        layer = nxt
    return layer.get((n, comb2(n)), 0)

def landau_fset_brute(n):
    """Backtracking enumeration (same as c3_spectrum_kpc1.py) for cross-check."""
    N2 = comb2(n)
    fset = set()

    def rec(k, last, psum, f):
        if k == n:
            if psum == N2:
                fset.add(f)
            return
        r = n - k
        for s in range(last, n):
            ps = psum + s
            if ps + (r - 1) * s > N2:
                break
            if ps < comb2(k + 1):
                continue
            if ps + (r - 1) * (n - 1) < N2:
                continue
            rec(k + 1, s, ps, f + comb2(s))

    rec(0, 0, 0, 0)
    return fset

def analyze(n):
    t0 = time.time()
    U = dp_f_bitset(n)
    elapsed = time.time() - t0
    C3 = comb3(n)
    assert U != 0, f"n={n}: empty DP result -- bug"
    # f-range
    fmin = (U & -U).bit_length() - 1
    fmax = U.bit_length() - 1
    assert fmax == C3, f"n={n}: max f {fmax} != C(n,3) {C3} (transitive missing)"
    assert fmin == near_regular_f(n), \
        f"n={n}: min f {fmin} != near-regular value {near_regular_f(n)}"
    mx = C3 - fmin
    assert mx == max_c3_formula(n), \
        f"n={n}: max c3 {mx} != formula {max_c3_formula(n)}"
    # gaps in c3 in [0, mx]  <=>  zero bits of U in [fmin, C3]
    bits = bin(U)[2:][::-1]  # bits[q] = '1' iff f=q achievable
    gaps_c3 = [C3 - q for q in range(fmin, C3 + 1) if bits[q] == '0']
    gaps_c3.sort()
    size = bits[fmin:C3 + 1].count('1')
    return mx, size, gaps_c3, elapsed

def main():
    print("=" * 78)
    print("CROSS-CHECK: DP bitset vs full backtracking enumeration, n = 3..13")
    print("=" * 78)
    for n in range(3, 14):
        U = dp_f_bitset(n)
        dp_set = {q for q in range(U.bit_length()) if (U >> q) & 1}
        brute = landau_fset_brute(n)
        status = "MATCH" if dp_set == brute else "*** MISMATCH ***"
        print(f"n={n:2d}: DP f-set size={len(dp_set):4d}  brute={len(brute):4d}  {status}")
        assert dp_set == brute
    print()
    print("=" * 78)
    print("EXACT c3 SPECTRUM, n = 3..30 (+ extension while budget lasts)")
    print("=" * 78)
    maxes, sizes = [], []
    total_t0 = time.time()
    budget_seconds = 480
    n = 3
    last_n = 0
    while True:
        if n > 30 and (time.time() - total_t0) > budget_seconds:
            break
        if n > 60:
            break
        mx, size, gaps, el = analyze(n)
        maxes.append(mx)
        sizes.append(size)
        gtxt = "NONE (gap-free)" if not gaps else f"{gaps}"
        print(f"n={n:2d}: max c3={mx:5d}  spectrum size={size:5d}  "
              f"gaps={gtxt}  [{el:.2f}s]")
        if gaps:
            print(f"   *** GAPS FOUND at n={n}: missing c3 values {gaps} ***")
        last_n = n
        n += 1
    print()
    print(f"Completed n = 3..{last_n} in {time.time()-total_t0:.1f}s total.")
    print()
    print("OEIS lookup material (search by VALUES):")
    print(f"  max c3   (n=3..{last_n}): {maxes}")
    print(f"  size     (n=3..{last_n}): {sizes}")
    print(f"  size == max+1 everywhere: {all(s == m + 1 for s, m in zip(sizes, maxes))}")

if __name__ == "__main__":
    main()
