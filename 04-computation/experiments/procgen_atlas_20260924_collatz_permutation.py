#!/usr/bin/env python3
"""
procgen_atlas_20260924_collatz_permutation.py

Implication-atlas lane (session collatz-procgen-20260922, 2026-09-24): Collatz's original
permutation (notebook 1932; open question publicised by Klamkin 1963; Lagarias's overview, arXiv
2111.02635, p. 3, writes it as U(2n)=3n, U(4n+1)=3n+1, U(4n+3)=3n+2 and asks whether the
iterates of 8 form an infinite set):

    g(3m) = 2m,  g(3m+1) = 4m+1,  g(3m+2) = 4m+3      (g = U^{-1})
    U(2n) = 3n,  U(4n+1) = 3n+1,  U(4n+3) = 3n+2

Sections:
  P1  bijection and exact identities (all checked on |n| <= 10^6 unless stated):
        g odd; g(n) = round(4n/3) for 3 !| n; U(m) = 3m/2 (m even), round(3m/4) (m odd);
        U = Mahler's W(m) = ceil(3m/2) on evens and floor(W(m)/2) on odds;
        U(k) = T^2(k) (k = 1 mod 4) and U(k) = T_-^2(k) (k = 3 mod 4)   [T_-: 3x-1 shortcut map];
        U(k) = T_+((k-1)/2) (k = 3 mod 4), U(k) = T_-((k+1)/2) (k = 1 mod 4);
        g(T(n)) = 2n+1 for every odd n, i.e. T(n) = U(2n+1) (the 3x+1 odd step IS one U-step);
        g(m) = 2 T_+^{-1}(m) + 1 (m = 2 mod 3), g(m) = 2 T_-^{-1}(m) - 1 (m = 1 mod 3),
        where T_+^{-1}(m) = (2m-1)/3 and T_-^{-1}(m) = (2m+1)/3 are the odd-branch predecessors;
        g(3m+1) = 2 g(3m) + 1 and g(3m+2) = 2 g(3m) + 3 = g(3m+1) + 2.
  P2  p-adic extensions: g is a full-branch expanding map of Z_3 (three branches, factor 3), U a
      full-branch expanding map of Z_2 (branch measures 1/2, 1/4, 1/4); both preserve Haar measure;
      the first k itinerary letters of n are determined by n mod 3^k (resp. a 2-adic prefix) and are
      uniformly distributed (exact counts).  Haar probability that the log-walk of g ever descends.
  P3  cycle census: every orbit meeting [1, N0] is one of the known cycles or exceeds B in both directions.
  P4  the orbit of 8: forward (g) and backward (U) to 10^5 steps; growth rates vs the drifts
      (1/3)log(2/3)+(2/3)log(4/3) = +0.056633 and (1/2)log(3/2)+(1/2)log(3/4) = +0.058892.
  P5  typing against the foundry controls (printed summary).
  P6  the shared gate: g-cycles live on the clocks 2^K vs 3^k of log_2 3 (signed carries), with an
      Eliahou-type cycle-length bound transferred to g.
"""
import math, sys
from fractions import Fraction as Fr

def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)

def g(n):
    m, r = divmod(n, 3)
    return 2 * m if r == 0 else (4 * m + 1 if r == 1 else 4 * m + 3)

def U(k):
    if k % 2 == 0:
        return 3 * (k // 2)
    n, r = divmod(k, 4)
    return 3 * n + 1 if r == 1 else 3 * n + 2

def T(x):   # 3x+1 shortcut
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2

def Tm(x):  # 3x-1 shortcut
    return x // 2 if x % 2 == 0 else (3 * x - 1) // 2

def W(m):   # Mahler's map ceil(3m/2)
    return -((-3 * m) // 2)

def rnd(fr):  # nearest integer to a Fraction with non-half fractional part
    f = math.floor(fr)
    return f if fr - f < Fr(1, 2) else f + 1

def P1():
    hdr("P1  bijection and exact identities")
    N = 10 ** 6
    for n in range(-N, N + 1):
        assert U(g(n)) == n and g(U(n)) == n
        assert g(-n) == -g(n)
    print(f"  U(g(n)) = g(U(n)) = n and g(-n) = -g(n) for |n| <= {N}: PASS (g is an odd permutation of Z)")
    for n in range(-3000, 3001):
        if n % 3:
            assert g(n) == rnd(Fr(4 * n, 3))
        m = n
        if m % 2 == 0:
            assert U(m) == 3 * m // 2 == W(m)
        else:
            assert U(m) == rnd(Fr(3 * m, 4)) == W(m) // 2
    print("  g(n) = round(4n/3) (3 !| n);  U(m) = 3m/2 = W(m) (m even), round(3m/4) = floor(W(m)/2) (m odd),")
    print("  with W(m) = ceil(3m/2) Mahler's map: PASS on |n| <= 3000.")
    for k in range(-N, N + 1, 2):  # odd k: start at -N (N even) -> -N is even; use k+1
        kk = k + 1
        r = kk % 4
        if r == 1:
            assert U(kk) == T(T(kk)) == Tm((kk + 1) // 2)
            assert (3 * kk + 1) % 4 == 0 and (3 * kk - 1) % 4 != 0
        else:
            assert U(kk) == Tm(Tm(kk)) == T((kk - 1) // 2)
            assert (3 * kk - 1) % 4 == 0 and (3 * kk + 1) % 4 != 0
    print("  odd k: U(k) = (3k+s)/4 with the unique sign s making 4 | 3k+s (Althoefer's one-player descending")
    print("  3n+-1 move); U(k) = T^2(k) = T_-((k+1)/2) for k = 1 mod 4 and U(k) = T_-^2(k) = T((k-1)/2) for k = 3 mod 4: PASS.")
    for n in range(-N + 1, N, 2):
        assert g(T(n)) == 2 * n + 1 and T(n) == U(2 * n + 1)
    print("  g(T(n)) = 2n+1, equivalently T(n) = U(2n+1), for every odd |n| < 10^6: PASS")
    for m in range(-N, N + 1):
        r = m % 3
        if r == 2:
            assert (2 * m - 1) % 3 == 0 and g(m) == 2 * ((2 * m - 1) // 3) + 1 and T((2 * m - 1) // 3) == m
        elif r == 1:
            assert (2 * m + 1) % 3 == 0 and g(m) == 2 * ((2 * m + 1) // 3) - 1 and Tm((2 * m + 1) // 3) == m
        else:
            assert g(m) == 2 * (m // 3)
    print("  g(m) = 2 T_+^{-1}(m) + 1 (m = 2 mod 3), g(m) = 2 T_-^{-1}(m) - 1 (m = 1 mod 3), g(3m) = 2m: PASS")
    print("  -> Collatz's permutation = (odd-branch predecessor on the sheet chosen by m mod 3) followed by x -> 2x+-1.")
    for m in range(-N, N + 1):
        assert g(3 * m + 1) == 2 * g(3 * m) + 1 and g(3 * m + 2) == 2 * g(3 * m) + 3
    print("  g(3m+1) = 2 g(3m) + 1, g(3m+2) = 2 g(3m) + 3: PASS  (the triple 3m,3m+1,3m+2 -> x, 2x+1, 2x+3, x=2m)")
    # the 4x+1 merge: T^3(4m+1) = T(m) for odd m
    for m in range(1, 200001, 2):
        assert T(T(T(4 * m + 1))) == T(m)
    print("  T^3(4m+1) = T(m) for odd m (the '4x+1 fractal recursion' of the T-tree): PASS; hence for odd m,")
    print("  T^3(g(3m+1)) = T(g(3m)/2).")

def P2():
    hdr("P2  p-adic extensions and itinerary statistics (exact)")
    # g on Z/3^(k+1) -> Z/3^k : each residue class r mod 3 is mapped onto Z/3^k bijectively
    for k in range(1, 9):
        M1, M0 = 3 ** (k + 1), 3 ** k
        for r in range(3):
            imgs = set(g(n) % M0 for n in range(r, M1, 3))
            assert len(imgs) == M0
    print("  g: for each r mod 3, {n = r mod 3} mod 3^(k+1) -> Z/3^k is a bijection (k <= 8): g extends to a")
    print("     full-branch map of Z_3 with three branches of expansion 3; Haar measure is invariant.")
    for k in range(1, 12):
        for (res, mod, e) in ((0, 2, 1), (1, 4, 2), (3, 4, 2)):
            M1, M0 = 2 ** (k + e), 2 ** k
            imgs = set(U(n) % M0 for n in range(res, M1, mod))
            assert len(imgs) == M0
    print("  U: the branches 2Z_2, 1+4Z_2, 3+4Z_2 map onto Z_2 with expansions 2, 4, 4 (k <= 11): Haar-invariant,")
    print("     itinerary letters Bernoulli(1/2, 1/4, 1/4).")
    # itinerary uniformity: first k ternary letters of g-orbit of n are uniform over n mod 3^k
    for k in range(1, 9):
        seen = {}
        for n in range(3 ** k):
            it = []
            x = n
            for _ in range(k):
                it.append(x % 3); x = g(x)
            it = tuple(it)
            seen[it] = seen.get(it, 0) + 1
        assert len(seen) == 3 ** k and all(v == 1 for v in seen.values())
    print("  the first k letters of the ternary g-itinerary are a bijection of Z/3^k (k <= 8): Terras-type theorem.")
    # Haar probability that the g log-walk (factors 2/3, 4/3, 4/3) ever goes below 0 (below its start)
    import math as m_
    l23, l43 = m_.log(2 / 3), m_.log(4 / 3)
    # states: (i2, i4) counts -> position i2*l23 + i4*l43 ; DP over steps with killing at < 0
    alive = {(0, 0): 1.0}
    dropped = 0.0
    for step in range(1, 3001):
        new = {}
        for (i2, i4), p in alive.items():
            for (a, b, w) in ((1, 0, 1 / 3), (0, 1, 2 / 3)):
                s = (i2 + a) * l23 + (i4 + b) * l43
                if s < -1e-12:
                    dropped += p * w
                elif s > 40:   # far above: cannot matter at 1e-12 level (drift positive)
                    pass
                else:
                    key = (i2 + a, i4 + b)
                    new[key] = new.get(key, 0.0) + p * w
        alive = new
    print(f"  Haar probability that the g-orbit's log-walk ever drops below its start: {dropped:.6f}")
    print("  (the analogue of Terras's theorem FAILS: a positive-density set of integers never descends in the")
    print("   model -- the drift is positive, as for 5x+1.)")

def P3(N0=100000, B=10 ** 24):
    hdr(f"P3  cycle census: all orbits meeting [1, {N0}], escape bound B = 10^{len(str(B)) - 1}")
    owner = bytearray(N0 + 1)
    cycles = []
    open_orbits = 0
    maxsteps = 0
    for n in range(1, N0 + 1):
        if owner[n]:
            continue
        owner[n] = 1
        x = g(n)
        steps = 0
        cyc = False
        seq = [n]
        while True:
            steps += 1
            if x == n:
                cyc = True
                break
            if x > B:
                break
            if x <= N0:
                owner[x] = 1
            seq.append(x)
            x = g(x)
        if cyc:
            cycles.append(sorted(seq))
            continue
        y = U(n)
        while y <= B:
            steps += 1
            if y <= N0:
                owner[y] = 1
            y = U(y)
        maxsteps = max(maxsteps, steps)
        open_orbits += 1
    for c in cycles:
        print(f"  cycle of length {len(c)}: min {c[0]}, max {c[-1]}: {c if len(c) <= 12 else c[:12]}")
    print(f"  distinct non-periodic-looking orbits meeting [1,{N0}] (exceed B both ways): {open_orbits}; max steps {maxsteps}")
    print(f"  FINITE-EXACT: every positive g-cycle with minimum <= {N0} and maximum <= B is one of the {len(cycles)} above.")
    return cycles

def P4(steps=100000):
    hdr(f"P4  the orbit of 8: {steps} steps forward (g) and backward (U)")
    for name, f in (("forward g", g), ("backward U", U)):
        x = 8
        mn = 8
        cnt = {}
        logs = []
        for s in range(1, steps + 1):
            key = x % 3 if f is g else (x % 4 if x % 2 else 0)
            cnt[key] = cnt.get(key, 0) + 1
            x = f(x)
            if x == 8:
                print(f"  {name}: RETURNED to 8 after {s} steps")
                break
            if x < mn:
                mn = x
            if s in (1000, 10000, 100000):
                logs.append((s, x.bit_length() * math.log(2)))
        rate = [(s, L / s) for s, L in logs]
        print(f"  {name}: min value seen {mn}; ln|x_s|/s at s=1e3,1e4,1e5: " +
              ", ".join(f"{r:.6f}" for _, r in rate) + f";  letter counts {dict(sorted(cnt.items()))}")
    print("  predicted drifts: g +0.056633 (letters mod 3 uniform), U +0.058892 (even : 1 mod 4 : 3 mod 4 = 2:1:1).")
    print("  8 is not periodic within 10^5 steps in either direction (known; the orbit of 8 is conjectured infinite).")

def P5():
    hdr("P5  typing of 'the orbit of 8 is infinite' against the foundry controls")
    rows = [
        ("SHEET", "g is odd (g(-n) = -g(n)): no sign asymmetry at all; but each step mixes the 3x+1 and 3x-1 odd",
         "branches (P1), so every sheet-specific Collatz tool is inapplicable"),
        ("DRIFT", "both directions expand (P2: +0.0566, +0.0589 per step); infinite orbits are the EXPECTED behaviour,",
         "as for 5x+1: the statement is an existence-of-divergence problem (reverse transfer), not a no-divergence one"),
        ("DEFECT", "a single orbit; density/Haar statements (P2) say nothing about 8",
         ""),
        ("INTEGRAL", "rational cycles of g exist for every periodic ternary itinerary (full-branch map)",
         "so an argument must isolate the integer 8"),
        ("UNIFORM", "g is an RCWA permutation; for the class of RCWA permutations cycle questions are undecidable-type",
         "(Kurtz-Simon type; see the note)"),
        ("LOGIC", "'orbit of 8 infinite' <=> 'g^k(8) != 8 for all k >= 1' is Pi^0_1 (refutable by finding a cycle)",
         "like Goldbach/RH; 'finite' is Sigma^0_1"),
    ]
    for r in rows:
        print(f"  {r[0]:8s}: {r[1]}")
        if r[2]:
            print(f"  {'':8s}  {r[2]}")

def P6(cycles, N0):
    hdr("P6  the gate shared with 3x+1: g-cycles sit on the clocks K/k of log_2 3 (exact)")
    print("  A g-step is x -> (c x + d)/3 with (c,d) = (2,0), (4,-1), (4,1) for x = 0,1,2 mod 3.  Over a word w")
    print("  of length k with a zeros and b = k-a nonzero letters, g^k(x) = (2^K x + D_w)/3^k, K = a + 2b,")
    print("  D_() = 0, D_{w r} = c_r D_w + d_r 3^{|w|}.  A cycle point satisfies  n (3^k - 2^K) = D_w,")
    print("  the same gate 2^K - 3^L as a 3x+1 cycle (Boehm-Sontacchi: n (2^K - 3^L) = B(w) > 0), but with a")
    print("  SIGNED carry D_w.  Product identity: 3^k/2^K = prod over nonzero letters of (1 + d_r/(4 n_i)).")
    from mpmath import mp, mpf, log as mlog, floor as mfloor
    mp.dps = 60
    alpha = mlog(3) / mlog(2)
    # convergents of alpha
    a_ = []; x = alpha
    for _ in range(40):
        ai = int(mfloor(x)); a_.append(ai); x = 1 / (x - ai)
    P = [a_[0], a_[0] * a_[1] + 1]; Q = [1, a_[1]]
    for n in range(2, 30):
        P.append(a_[n] * P[-1] + P[-2]); Q.append(a_[n] * Q[-1] + Q[-2])
    conv = set(zip(P, Q))
    for c in cycles:
        n0 = c[0]
        x = n0; w = []
        while True:
            w.append(x % 3); x = g(x)
            if x == n0:
                break
        k = len(w); a = w.count(0); b = k - a; K = a + 2 * b
        D = 0
        for i, r in enumerate(w):
            cr, dr = ((2, 0), (4, -1), (4, 1))[r]
            D = cr * D + dr * 3 ** i
        G = 3 ** k - 2 ** K
        assert n0 * G == D
        from math import gcd
        gg = gcd(K, k)
        red = (K // gg, k // gg)
        print(f"  cycle min {n0:3d}: k={k:2d} a={a} b={b} K={K:2d}  clock K/k={K}/{k}  gate 3^k-2^K={G:6d}  D_w={D:8d}"
              f"  n0 = D/gate: {D // G == n0}   reduced clock {red[0]}/{red[1]} convergent: {red in conv}")
    print("  (Shanks 1965 [R] already noted that the cycle lengths 1, 2, 5, 12 are convergent denominators of log_2 3;")
    print("   Atkin 1966 [R]: no further cycle of period < 200; Simons 2022 [P, arXiv 2205.10582]: finitely many m-cycles")
    print("   by Baker/Rhin, none for m <= 2 besides the known ones.  New here: the signed-carry comparison below.)")
    # unit-gap clocks (Catalan / Gersonides: |3^k - 2^K| = 1 only for (K,k) = (1,1), (2,1), (3,2))
    from itertools import product
    found = set()
    for k in (1, 2):
        for w in product((0, 1, 2), repeat=k):
            b = sum(1 for r in w if r); K = k + b
            G = 3 ** k - 2 ** K
            if abs(G) != 1:
                continue
            D = 0
            for i, r in enumerate(w):
                cr, dr = ((2, 0), (4, -1), (4, 1))[r]
                D = cr * D + dr * 3 ** i
            n0 = D // G
            assert n0 * G == D and n0 % 3 == w[0] % 3
            x = n0; cyc = []
            for _ in range(k):
                cyc.append(x); x = g(x)
            assert x == n0
            found.add(tuple(sorted(cyc)))
    print(f"  unit-gap clocks (1/1, 2/1, 3/2; |3^k - 2^K| = 1 only there by Catalan/Gersonides): every word gives an")
    print(f"  integer cycle, and the cycles are exactly {sorted(found)}: the unit gaps carry only the known small cycles.")
    print("  3x+-1 comparison (inherited gate census): the signed 3x+-1 cycles sit on 1/1, 2/1, 3/2, 11/7")
    print("  (gates -1, 1, -1, -139); Collatz's permutation uses 2/1, 3/2, 8/5, 19/12 (gates -1, 1, -13, 7153).")
    print("  The unit (Catalan) gaps 2/1 and 3/2 carry cycles of BOTH; 13 and 7153, which no 3x+-1 carry cancels")
    print("  (all 35 words on 8/5 have q=13; min q=23 on 19/12), are cancelled by g's signed carries.")
    # Eliahou-type bound for g-cycles all of whose elements exceed N0
    delta = 1 / ((4 * mpf(N0) - 1) * mlog(2))
    best = None
    for qq in range(1, 200000):
        pp = int(mfloor(alpha * qq + mpf(1) / 2))
        for p2 in (pp - 1, pp, pp + 1):
            if abs(mpf(p2) / qq - alpha) < delta:
                best = (p2, qq); break
        if best:
            break
    print(f"  Eliahou-type transfer: a g-cycle with every element > N0 = {N0} has |K/k - log_2 3| < 1/((4N0-1) ln 2)")
    print(f"  = {float(delta):.3e}, so k >= {best[1]} (smallest denominator in that window, fraction {best[0]}/{best[1]}).")
    print(f"  With P3: any g-cycle besides the four has min > {N0} (hence length >= {best[1]}) or max > 10^24.")

if __name__ == "__main__":
    P1(); P2()
    N0 = 10 ** 6
    cyc = P3(N0=N0)
    P4(); P5(); P6(cyc, N0)
    print("\nDONE collatz_permutation")
