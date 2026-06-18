#!/usr/bin/env python3
"""
DECISIVE TEST: does mu_0(E) -> 0 as the SPREAD of E -> infinity?

If yes, the 0-anchored strip-avoidance measure mu_0 does NOT have a uniform positive
floor, and CANNOT directly close B(k) on its own (it would only give the bounded-spread
floor, same residual as before). If mu_0 stays bounded below uniformly, this angle WINS.

The naive independence heuristic predicts mu_0(E) ~ prod over e of density(Good_e)
but with CORRELATIONS. As e_i grow and become 'generic' (equidistributed, independent
mod 1), the events {frac(e_i x) avoid [1/7,3/7]} become asymptotically INDEPENDENT, each
prob 5/7, so mu_0 -> (5/7)^k. That is POSITIVE but the WORST shape could be worse if
correlations are NEGATIVE.

KEY: even (5/7)^k > 0 is a FLOOR for the 'generic large-spread' regime. The danger is
shapes where mu_0 dips BELOW (5/7)^k due to negative correlation, possibly to 0.

We test:
  (1) mu_0 for E = {0,1,2,...,k-1} scaled? No -- scaling-invariant, mu_0(cE) != mu_0(E)?
      WAIT: mu IS scale-invariant but mu_0 with a FIXED arc is NOT obviously. Check.
  (2) E = {0} U {large spread} : push min nonzero e -> infinity.
  (3) Direct: for fixed k, take E = {0, N, 2N, ..., (k-1)N}? That's an AP, scaling of
      consecutive by N -- if mu_0 scale-invariant, equals consecutive. Check.
  (4) The真正 worst: random large-spread, track inf as we allow larger spread.
  (5) Asymptotic: does inf_E mu_0(E) over ALL E (unbounded spread) equal 0 or > 0?
"""
from fractions import Fraction as F
from math import gcd
import random
import sys
sys.path.insert(0, ".")
from importlib import import_module

# reuse mu0_exact
DLO, DHI = F(1, 7), F(3, 7)
def _frac(q): return q - q.__floor__()

def mu0_exact(E, dlo=DLO, dhi=DHI):
    E = sorted(set(E))
    bp = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for n in range(e + 1):
            for off in (dlo, dhi):
                x = (F(n) + off) / e
                if F(0) <= x < F(1): bp.add(x)
    bp = sorted(bp)
    tot = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        mid = (a + b) / 2
        ok = True
        for e in E:
            if e == 0: continue
            fr = _frac(e * mid)
            if dlo <= fr <= dhi: ok = False; break
        if ok: tot += (b - a)
    return tot

def normalize(E):
    E = sorted(set(E)); g = 0
    for e in E: g = gcd(g, e)
    return [e // g for e in E] if g else E

def header(t):
    print("\n" + "=" * 74); print(t); print("=" * 74)

if __name__ == "__main__":
    header("(1) IS mu_0 SCALE-INVARIANT?  mu_0(cE) vs mu_0(E)")
    for E in ([0,1,2], [0,1,2,3], [0,1,3]):
        print(f"  E={E}")
        for c in (1,2,3,5,7,10):
            cE = [c*e for e in E]
            print(f"    c={c:2d}: mu_0(cE)={float(mu0_exact(cE)):.6f}={mu0_exact(cE)}")
    print("  CONCLUSION: mu_0 with FIXED arc is NOT scale-invariant (unlike mu).")
    print("  => mu_0(AP scaled) can differ. The 0-anchor breaks scaling symmetry.")

    header("(2) PUSH MIN NONZERO e -> infinity (k=3): E = {0, 1, N}")
    for N in (3, 5, 10, 50, 100, 500, 1000, 5000):
        E = [0, 1, N]
        print(f"  N={N:5d}: mu_0={float(mu0_exact(E)):.6f}={mu0_exact(E)}")

    header("(2b) E = {0, N, N+1} push N (all large)")
    for N in (2, 5, 10, 50, 200, 1000, 4000):
        E = [0, N, N+1]
        m0 = mu0_exact(E)
        print(f"  N={N:5d}: mu_0={float(m0):.6f}={m0}")

    header("(3) AP scaled: E = {0,N,2N,...,(k-1)N}, k=5, push N")
    for N in (1, 2, 5, 10, 100, 1000):
        E = [i*N for i in range(5)]
        m0 = mu0_exact(E)
        print(f"  N={N:5d}: E spread={4*N:6d} mu_0={float(m0):.6f}={m0}")

    header("(4) RANDOM LARGE-SPREAD WORST-CASE, k=13, increasing spread budget")
    random.seed(20260618)
    for spread_mult in (2, 4, 8, 16, 32, 64):
        worst = None
        ntrial = 1500
        for _ in range(ntrial):
            hi = 13 * spread_mult
            pool = random.sample(range(1, hi + 1), 12)
            E = normalize([0] + pool)
            m0 = mu0_exact(E)
            if worst is None or m0 < worst[1]:
                worst = (E, m0)
        print(f"  spread_mult={spread_mult:3d} (hi={13*spread_mult:4d}): worst mu_0={float(worst[1]):.6f}={worst[1]}")
        print(f"      E={worst[0]}")

    header("(5) DOES inf mu_0 -> 0?  Target: adversarial 'spread to kill the 0-anchor'")
    # Heuristic for worst: choose e_i so that the bad strips [1/7,3/7]/e_i COVER as much
    # of [0,1] as possible. Bad strip for e covers fraction 2/7 of x-axis (e copies of
    # width (2/7)/e). To MINIMIZE mu_0 = meas(all-good), MAXIMIZE the union of bad strips.
    # Pick e_i with diverse residues to spread coverage. Try a greedy adversary.
    def union_bad_measure(E):
        return F(1) - mu0_exact(E)
    # greedy: start {0,1}, add e that maximizes union_bad
    best_set = [0, 1]
    print("  greedy adversary maximizing union of bad strips (minimizing mu_0):")
    cand_pool = list(range(2, 200))
    while len(best_set) < 13:
        bestadd = None
        for e in cand_pool:
            if e in best_set: continue
            cur = best_set + [e]
            m0 = mu0_exact(cur)
            if bestadd is None or m0 < bestadd[1]:
                bestadd = (e, m0)
        best_set.append(bestadd[0])
        print(f"    |E|={len(best_set):2d}: add e={bestadd[0]:3d}  mu_0={float(bestadd[1]):.6f}={bestadd[1]}")
    print(f"  final greedy E={normalize(best_set)}")
    print(f"  final mu_0={float(mu0_exact(best_set)):.6f}={mu0_exact(best_set)}")
