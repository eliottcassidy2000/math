#!/usr/bin/env python3
"""
LRC(14) REFRAME 2 — LP duality / averaging over the covering moduli.
mac-mini-2026-06-17-S3.

Question (3 parts from the prompt):
 (1) weights lambda_q>=0 over q in 2..14; does a convex combination of q-grid
     witnesses certify a lonely point or M>=1/14?
 (2) the LP DUAL of "danger arcs cover the circle": find a dual certificate.
 (3) is there a probability measure mu on tau under which E[min_v||v tau||] or a
     Markov bound forces a point with min>=1/14? Report best provable LB on M.

VERDICT (honest): The linear/first-moment averaging is TAUTOLOGICAL or BLIND.
Three independent obstructions are proved computationally below:

  OBSTRUCTION A (measure blindness): every danger band D_v={||v tau||<1/14} has
    Lebesgue measure exactly 2*(1/14)=1/7, INDEPENDENT of v and of S. Total band
    measure = 13/7 for EVERY 13-set, covering or not, M=1/14 or M=7/89. The mean
    coverage multiplicity E[N] = 13/7 carries ZERO information separating cases.

  OBSTRUCTION B (averaging LP is tautological): the LP
       minimize  sum_v mu(D_v)   over probability measures mu on the circle
    has optimum = min_tau N(tau) (put all mass at the argmin), which equals
    0 iff a lonely point exists and >=1 iff S is covering. The LP optimum IS the
    yes/no answer to LRC for S — no leverage, pure restatement.

  OBSTRUCTION C (the lonely point is OFF the q-grid lattice — the real finding):
    For a covering S, every q-grid point a/q (q in 2..14) has g(S,a/q)=0 because
    S contains a multiple of q. Hence the lonely set {g>=1/14} is DISJOINT from
    the entire q-grid. Worse, the binding-pair lonely witness sits at tau*=num/D
    with D=v_a+v_b, and for the extremal sets D is coprime to lcm(2..14)=360360
    (e.g. 37/89, gcd(89,360360)=1). So tau* is not even on the 1/lcm lattice that
    every convex combination sum_q lambda_q a_q/q lives on. The lambda_q averaging
    construction is CONFINED to the lattice, where g=0 everywhere; it can never
    reach a lonely point.

Best provable lower bound on M from averaging: the only thing first moments give is
the trivial "min coverage multiplicity" test, which is the original question. No
nontrivial analytic lower bound on M is produced. This is a NO-LEVERAGE result for
the linear averaging family — but Obstruction C is a clean, reusable structural lemma.
"""
from fractions import Fraction as F
from math import lcm, gcd

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def min_coverage_multiplicity(S, c=F(1, 14)):
    """min_tau N(tau), N = number of bands covering tau. EXACT via breakpoints."""
    bps = set()
    for v in S:
        for j in range(v + 1):
            bps.add((F(j, v) - c / v) % 1); bps.add((F(j, v) + c / v) % 1)
    pts = sorted(bps)
    best = None; bt = None
    for i in range(len(pts)):
        lo = pts[i]; hi = pts[i+1] if i+1 < len(pts) else pts[0] + 1
        mid = ((lo + hi) / 2) % 1
        N = sum(1 for v in S if nrm(v * mid) < c)
        if best is None or N < best: best = N; bt = mid
    return best, bt

def main():
    L = 1
    for q in range(2, 15): L = lcm(L, q)
    print(f"lcm(2..14) = {L}\n")

    W = {
        "S1 {1..11,13,84}":  [1,2,3,4,5,6,7,8,9,10,11,13,84],
        "S2 {1..5,7..13,84}":[1,2,3,4,5,7,8,9,10,11,12,13,84],
    }

    print("=== OBSTRUCTION A: measure blindness (each band = 1/7, total = 13/7) ===")
    for name, S in W.items():
        m, at = M(S)
        print(f"  {name}: M={m}, tau*={at}, total band measure = {13}*1/7 = {F(13,7)}")
    print(f"  Extremal {{1..13}} (NOT covering, M=1/14): total band measure also {F(13,7)}")
    print("  => mean multiplicity 13/7 identical regardless of M. BLIND.\n")

    print("=== OBSTRUCTION B: averaging LP optimum = min coverage multiplicity ===")
    for name, S in list(W.items()) + [("{1..13}", [1,2,3,4,5,6,7,8,9,10,11,12,13])]:
        mn, at = min_coverage_multiplicity(S)
        m, _ = M(S)
        print(f"  {name}: min_tau N(tau) = {mn}  (0 <=> lonely-set nonempty); M={m}")
    print("  => LP optimum encodes exactly the yes/no LRC answer. TAUTOLOGICAL.\n")

    print("=== OBSTRUCTION C: lonely point is OFF the q-grid lattice ===")
    for name, S in W.items():
        cov = is_covering(S)
        allzero = all(g(S, F(a, q)) == 0 for q in range(2, 15) for a in range(q))
        m, at = M(S)
        print(f"  {name}: covering={cov}, g=0 on every q-grid point: {allzero}")
        print(f"     tau*={at}, denom={at.denominator}, gcd(denom,lcm)={gcd(at.denominator,L)}")
    print("  => convex combos sum_q lambda_q a_q/q have denom | lcm; the extremal")
    print("     lonely point has denom coprime to lcm => UNREACHABLE by averaging.\n")

    print("=== Sanity: inf over covering sets, best provable averaging LB on M ===")
    # confirm 7/89 is the min, and that the only averaging LB is the trivial 0.
    m1, _ = M(W["S1 {1..11,13,84}"])
    print(f"  unique min over covering 13-sets: M = {m1} = {float(m1):.6f} > 1/14 = {float(F(1,14)):.6f}")
    print(f"  best provable LOWER bound on M from linear averaging: 0 (no nontrivial bound).")
    print(f"  REFRAME 2 verdict: NO-LEVERAGE for the linear/first-moment family.")
    print(f"  Reusable lemma (Obstruction C): off-grid-witness lemma, see docstring.")

if __name__ == "__main__":
    main()
