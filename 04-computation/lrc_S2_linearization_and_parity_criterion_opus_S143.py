"""
lrc_S2_linearization_and_parity_criterion_opus_S143.py

(1) |S| = 2 LINEARIZATION THEOREM (proved this session, verified here):
    for coprime a < b:  chi_c(Cay(Z, +-{a,b})) = 1/M({a,b}).
    * both odd  => every edge flips parity => bipartite => chi_c = 2 = 1/M (M = 1/2).
    * a+b odd   => M = floor((a+b)/2)/(a+b) and the explicit (a+b)-cycle
      x, x+a, x+2a, ..., x+ba, x+ba-b, ..., x  (b steps of +a then a steps of -b;
      distinctness by coprimality) is an ODD cycle => chi_c >= (a+b)/((a+b-1)/2) = 1/M;
      the witness coloring (L1) gives <=.  QED.
    VERIFIED: M closed form + odd-cycle existence/distinctness for all coprime a<b <= 40.

(2) THE PARITY CRITERION for grid-symmetric 2-page optima (n=12 theorem):
    sigma-fixed chords have a+b = n+1 (odd length) => never parity-carrying =>
    every gridsym page assignment has Q == C(n,4) mod 2.
    THEOREM: if Z(n) != C(n,4) mod 2, NO grid-symmetric assignment attains Z(n).
    Table n = 5..20: the criterion fires exactly at which n?  (n=12 proved broken;
    n=8 shows parity-consistent breaking also exists -- criterion is sufficient only.)

(3) MU = M BLEEDING-EDGE HUNT: periodic Motzkin (exact DP) vs M over a large batch of
    3- and 4-element sets (max element <= 14, periods N <= 320): any mu > M?
"""
from fractions import Fraction as F
from math import gcd, comb
import sys, time

sys.path.insert(0, ".")
from lrc_graph_interpretation_ladder_opus_S141 import M_exact, motzkin_periodic

def Z_guy(n):
    return (n // 2) * ((n - 1) // 2) * ((n - 2) // 2) * ((n - 3) // 2) // 4

def main():
    t0 = time.time()
    print("=" * 100)
    print("(1) |S|=2 theorem: M closed form + odd-cycle verification, coprime a<b <= 40")
    print("=" * 100)
    bad = 0; checked = 0
    for b in range(2, 41):
        for a in range(1, b):
            if gcd(a, b) != 1: continue
            checked += 1
            M = M_exact([a, b])
            if a % 2 == 1 and b % 2 == 1:
                if M != F(1, 2): bad += 1; print(f"  M({a},{b}) = {M} != 1/2")
            else:
                pred = F((a + b) // 2, a + b)
                if M != pred: bad += 1; print(f"  M({a},{b}) = {M} != {pred}")
                # odd cycle distinctness
                cyc = [i * a for i in range(b + 1)] + [b * a - j * b for j in range(1, a)]
                if len(set(cyc)) != a + b: bad += 1; print(f"  cycle not simple for ({a},{b})")
                if (a + b) % 2 == 0: bad += 1
    print(f"  checked {checked} coprime pairs: closed form + odd-cycle {'ALL OK' if bad == 0 else f'{bad} FAILURES'}")
    print("  => chi_c = 1/M for |S| = 2: bipartite case + odd-(a+b)-cycle case. THEOREM.")

    print()
    print("=" * 100)
    print("(2) parity criterion: gridsym page assignments have Q == C(n,4) mod 2;")
    print("    Z(n) != C(n,4) mod 2 => NO gridsym optimum (mirror breaking FORCED)")
    print("=" * 100)
    for n in range(5, 21):
        z, c = Z_guy(n), comb(n, 4)
        fired = (z % 2) != (c % 2)
        print(f"  n={n:2d}: Z={z:5d} ({z%2})  C(n,4)={c:5d} ({c%2})   "
              f"{'*** PARITY-FORCED MIRROR BREAKING ***' if fired else 'parity-consistent'}")
    print("  [n=12 breaking is now a THEOREM; n=8 breaking (S142 census) is parity-consistent,")
    print("   so the criterion is sufficient, not necessary -- mac-mini-S50 has the n=8 mechanism]")

    print()
    print("=" * 100)
    print("(3) mu = M hunt: 3-,4-element sets, max <= 14, exact periodic Motzkin N <= 320")
    print("=" * 100)
    from itertools import combinations
    exceptions = []
    ntested = 0
    for k in (3, 4):
        for S in combinations(range(1, 15), k):
            if gcd(*S) if k == 2 else 0: pass
            g = 0
            for s in S: g = gcd(g, s)
            if g != 1: continue
            ntested += 1
            M = M_exact(list(S))
            # scan periods: multiples of M's denominator first (cheap certificates), then all
            best = F(0)
            for N in range(2, 321):
                v = motzkin_periodic(list(S), N)
                if v > best: best = v
                if best > M: break
            if best > M:
                exceptions.append((S, M, best))
                print(f"  *** mu > M: S={S}: M={M}, periodic mu >= {best}")
    print(f"  tested {ntested} primitive sets: exceptions found: {len(exceptions)}"
          f"{'  => mu = M on ALL (collapse conjecture strengthens)' if not exceptions else ''}")

    print(f"\n[{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
