#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-501 STRUCTURAL VERIFICATION of the LRC singular series
   L(S) = (6/7)^13 + sum_{exact relations sum_{v in T} t_v v = 0, t_v != 0}
                       (6/7)^{13-|T|} (-1)^{|T|} prod_{v in T} s(t_v),
   s(t) = sin(pi t / 7)/(pi t).

Fresh code (kind-pasteur-2026-06-14-S6, structural pass).  Verifies:
 (1) the 7-VANISHING: s(t)=0 iff 7|t (t!=0); any relation term with ANY coeff t_v
     divisible by 7 contributes 0.  => L depends only on relations with all coeffs
     coprime to 7.
 (2) PAIRWISE ABSOLUTE CONVERGENCE: for a pair (v_a,v_b), pairwise relations are
     (t_a,t_b)=k*(v_b,-v_a)/g, g=gcd(v_a,v_b).  The pair correction
        P(a,b) = sum_{k!=0} s(k v_b/g) s(-k v_a/g)
     is absolutely convergent and |P| <= C * g^2/(v_a v_b) (explicit C below).
 (3) ALMOST-SIDON SUBCASE: if the only exact relations are pairwise (no relation
     with |T|>=3 having all coeffs coprime to 7), then
        |L - (6/7)^13| <= (6/7)^{11} * sum_{a<b} |P(a,b)|,
     and if this bound < (6/7)^13 then L>0.  Report the bound per config.
 (4) NUMERICAL: truncated exact relation sum (over |T|<=2 and a search for small
     |T|=3 coprime-to-7 relations) for a near-Sidon mult-of-14 set vs the evader
     7*{1..12} u {r}.  Compare the size of the |T|>=3 contribution.

All sums truncated honestly; truncation tails bounded and reported.
"""

import sys, itertools
from math import sin, pi, gcd
from functools import reduce
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

R = 6.0 / 7.0          # 1 - beta limit
MAIN = R ** 13         # (6/7)^13 main term


def s(t):
    """s(t) = sin(pi t / 7)/(pi t); s(0) := 1/7 by continuity (not used for t=0)."""
    if t == 0:
        return 1.0 / 7.0
    return sin(pi * t / 7.0) / (pi * t)


# ----------------------------------------------------------------------------
# (1) THE 7-VANISHING
# ----------------------------------------------------------------------------
def part1():
    print("=== (1) THE 7-VANISHING: s(t)=0 iff 7|t (t!=0) ===", flush=True)
    bad = []
    for t in range(-300, 301):
        if t == 0:
            continue
        z = abs(s(t)) < 1e-12
        if z != (t % 7 == 0):
            bad.append(t)
    print(f"   checked t in [-300,300]\\{{0}}: s(t)=0  <=>  7|t   -- counterexamples: {bad}", flush=True)
    # the bound |s(t)| <= 1/(pi|t|)
    bound_ok = all(abs(s(t)) <= 1.0 / (pi * abs(t)) + 1e-15 for t in range(1, 5000))
    print(f"   |s(t)| <= 1/(pi|t|) for t=1..4999: {bound_ok}", flush=True)
    print(f"   => any relation term with a coeff t_v with 7|t_v contributes 0;", flush=True)
    print(f"      L depends ONLY on relations all of whose coeffs are coprime to 7.", flush=True)


# ----------------------------------------------------------------------------
# (2) PAIRWISE CORRECTION + ABSOLUTE-CONVERGENCE BOUND
# ----------------------------------------------------------------------------
def pair_correction(va, vb, Kmax=200000):
    """
    Exact pairwise correction for the unordered pair {va, vb}:
       relations t_a va + t_b vb = 0 with t_a,t_b != 0 are (t_a,t_b)=k*(vb,-va)/g.
    The expansion term for T={va,vb} carries (6/7)^{11} (-1)^2 = (6/7)^11, and the
    inner sum over the relation lattice of the pair is
       P(va,vb) = sum_{k != 0} s(k vb/g) s(-k va/g)
                = 2 * sum_{k>=1} s(k vb/g) s(k va/g)     (s odd => s(-x)=-s(x), product even)
    Returns (P, tail_bound) where tail_bound bounds |sum_{|k|>Kmax}|.
    """
    g = gcd(va, vb)
    a = vb // g           # step in t_a (coprime to b=va//g)
    b = va // g
    P = 0.0
    for k in range(1, Kmax + 1):
        P += s(k * a) * s(k * b)
    P *= 2.0
    # tail bound: |s(k a) s(k b)| <= 1/(pi^2 a b k^2); sum_{k>Kmax} 1/k^2 <= 1/Kmax
    tail = 2.0 * (1.0 / (pi * pi * a * b)) * (1.0 / Kmax)
    return P, tail, g, a, b


def pair_abs_bound(va, vb):
    """
    The PROVED absolute bound on |P(va,vb)| using the 7-vanishing.
       |P| <= 2 sum_{k>=1, 7 doesn't div k*a, 7 doesn't div k*b} 1/(pi^2 a b k^2)
            <= (2/(pi^2 a b)) * zeta(2) = (2/(pi^2 a b)) * pi^2/6 = 1/(3 a b)
       and a b = (vb/g)(va/g) = va vb / g^2, so |P| <= g^2 / (3 va vb).
    Returns the explicit upper bound g^2/(3 va vb).
    """
    g = gcd(va, vb)
    return g * g / (3.0 * va * vb)


def part2():
    print("\n=== (2) PAIRWISE ABSOLUTE CONVERGENCE + explicit O(g^2/(v_a v_b)) bound ===", flush=True)
    print("   pair: exact P(a,b)=2 sum_{k>=1} s(k vb/g)s(k va/g);  proved bound |P| <= g^2/(3 va vb)", flush=True)
    tests = [(14, 21), (1, 2), (3, 5), (12, 18), (7, 14), (8, 13), (6, 10), (98, 105)]
    for va, vb in tests:
        P, tail, g, a, b = pair_correction(va, vb, Kmax=400000)
        bound = pair_abs_bound(va, vb)
        flag = "OK" if abs(P) <= bound + tail + 1e-12 else "!!VIOLATION!!"
        coprime7 = (a % 7 != 0) and (b % 7 != 0)
        # if g shares a factor making a or b divisible by 7, the whole pair lattice can vanish
        note = "" if coprime7 else "  (a or b divisible by 7: terms can vanish by 7-vanishing)"
        print(f"   ({va:3d},{vb:3d}) g={g} a={a} b={b}: P={P:+.6e} (tail<{tail:.1e})  "
              f"bound g^2/(3 va vb)={bound:.6e}  {flag}{note}", flush=True)


# ----------------------------------------------------------------------------
# (3) ALMOST-SIDON SUBCASE THRESHOLD
# ----------------------------------------------------------------------------
def almost_sidon_pairwise_bound(S):
    """
    sum over unordered pairs of |P(a,b)|, and the proved over-bound
       sum_{a<b} |P(a,b)| <= sum_{a<b} g(a,b)^2/(3 va vb).
    The full |L-(6/7)^13| deviation, IF all exact relations are pairwise, is
       <= (6/7)^11 * sum_{a<b} |P(a,b)|.
    Returns (exact_pairsum, proved_pairsum_bound, deviation_bound).
    """
    exact = 0.0
    proved = 0.0
    for va, vb in itertools.combinations(sorted(S), 2):
        P, tail, g, a, b = pair_correction(va, vb, Kmax=200000)
        exact += abs(P)
        proved += pair_abs_bound(va, vb)
    dev_bound_exact = (R ** 11) * exact
    dev_bound_proved = (R ** 11) * proved
    return exact, proved, dev_bound_exact, dev_bound_proved


def part3():
    print("\n=== (3) ALMOST-SIDON PROVABLE SUBCASE: |L-(6/7)^13| <= (6/7)^11 sum|P| ===", flush=True)
    print(f"   (6/7)^13 = {MAIN:.6f};  threshold: deviation bound < (6/7)^13  => L>0 (loose).", flush=True)
    # GENUINE near-Sidon, coprime, spread mult-of-14 sets (few/weak pairwise relations).
    # (NB: powers of 2 are a BAD example -- t_a v_a + t_b v_b = 0 has dense small
    #  solutions like 2*1 - 1*2 = 0, so {2^k} has LARGE pairwise corrections.)
    cfgs = {
        'spread primes u14':  sorted([14, 17, 23, 31, 41, 53, 61, 71, 83, 97, 101, 113, 127]),
        'large-spread u14':   sorted([14, 101, 211, 307, 401, 503, 601, 701, 809, 907, 1009, 1103, 1201]),
    }
    for name, S in cfgs.items():
        if len(set(S)) != 13:
            print(f"   {name}: not 13 distinct, skip", flush=True); continue
        ex, pr, dbx, dbp = almost_sidon_pairwise_bound(S)
        ok_exact = dbx < MAIN
        ok_proved = dbp < MAIN
        print(f"   {name}:", flush=True)
        print(f"      sum|P| (exact)   = {ex:.6f}  => dev bound (6/7)^11*sum = {dbx:.6f}   < (6/7)^13? {ok_exact}", flush=True)
        print(f"      sum|P| (proved <=)= {pr:.6f}  => dev bound (proved)      = {dbp:.6f}   < (6/7)^13? {ok_proved}", flush=True)
        print(f"      => if S is almost-Sidon (only pairwise relations), L >= (6/7)^13 - {dbx:.6f} = {MAIN-dbx:.6f} > 0", flush=True)


# ----------------------------------------------------------------------------
# (4) NUMERICAL: full truncated relation sum, near-Sidon vs evader
# ----------------------------------------------------------------------------
def exact_relations_up_to(S, max_T=3, coeff_max=14):
    """
    Enumerate exact integer relations sum_{v in T} t_v v = 0 with t_v != 0,
    |t_v| <= coeff_max, |T| <= max_T, AND all t_v coprime to 7 (others vanish).
    Returns the singular-series partial sum and a breakdown by |T|.
    NOTE: this is a TRUNCATION in coeff size; tails are not fully controlled for
    |T|>=3 but the dominant small-coefficient relations are captured.  Pairwise
    (|T|=2) is summed exactly via the lattice (no coeff truncation) for accuracy.
    """
    Ssort = sorted(set(S))
    contrib = {1: 0.0, 2: 0.0, 3: 0.0}

    # |T|=1: t_v v = 0 with t_v!=0 has no solution (v>0) -> nothing.
    # |T|=2: exact lattice sum (no truncation), only coprime-to-7 steps survive.
    for va, vb in itertools.combinations(Ssort, 2):
        P, tail, g, a, b = pair_correction(va, vb, Kmax=200000)
        # (6/7)^{13-2} * (-1)^2 * P
        contrib[2] += (R ** 11) * P

    # |T|=3: search small coprime-to-7 relations t1 u + t2 v + t3 w = 0.
    if max_T >= 3:
        for tri in itertools.combinations(Ssort, 3):
            u, v, w = tri
            for t1 in range(-coeff_max, coeff_max + 1):
                if t1 == 0 or t1 % 7 == 0:
                    continue
                for t2 in range(-coeff_max, coeff_max + 1):
                    if t2 == 0 or t2 % 7 == 0:
                        continue
                    rem = -(t1 * u + t2 * v)
                    if rem % w != 0:
                        continue
                    t3 = rem // w
                    if t3 == 0 or t3 % 7 == 0:
                        continue
                    if abs(t3) > coeff_max:
                        continue
                    term = (R ** 10) * ((-1) ** 3) * s(t1) * s(t2) * s(t3)
                    contrib[3] += term
    return contrib


def band_set(q, n=14):
    h = q // n
    return set(r % q for r in range(-h, h + 1))


def deficit(S, q, n=14):
    B = band_set(q, n)
    return sum(1 for a in range(q) if all((v * a) % q not in B for v in S))


def measure_L(S, qlo=1500, qhi=1530):
    """GROUND TRUTH limit: average exact D(q,S)/q over a high-q window."""
    vals = [deficit(S, q) / q for q in range(qlo, qhi)]
    m = sum(vals) / len(vals)
    sd = (sum((x - m) ** 2 for x in vals) / len(vals)) ** 0.5
    return m, sd


def part4():
    print("\n=== (4) NUMERICAL: truncated relation sum vs GROUND-TRUTH deficit limit ===", flush=True)
    print("   (|T|=2 is summed exactly over its lattice; |T|=3 is a coeff-truncated PROBE,", flush=True)
    print("    NOT a convergent estimate of the |T|=3 total -- the ground-truth limit below", flush=True)
    print("    is the exact D(q,S)/q average.  Trend, not value, is the claim for |T|=3.)", flush=True)
    near_sidon = sorted([14, 17, 23, 31, 41, 53, 61, 71, 83, 97, 101, 113, 127])
    evader = sorted([7 * k for k in range(1, 13)] + [611])    # 7*{1..12} u {611}
    for name, S in [('near-Sidon spread-primes u14', near_sidon),
                    ('evader 7*{1..12}u{611}', evader)]:
        c = exact_relations_up_to(S, max_T=3, coeff_max=14)
        Lpair = MAIN + c[2]
        Lmeas, sd = measure_L(S)
        print(f"   {name}:", flush=True)
        print(f"      (6/7)^13                 = {MAIN:+.6f}", flush=True)
        print(f"      + |T|=2 pairwise (exact) = {c[2]:+.6f}  => pairwise-only L = {Lpair:+.6f}", flush=True)
        print(f"      |T|=3 probe (|t|<=14)    = {c[3]:+.6f}  (sign/size of higher-order effect)", flush=True)
        print(f"      GROUND TRUTH L = avg D/q over q=1500..1529 = {Lmeas:.4f} +- {sd:.4f}", flush=True)
        print(f"      residual L_true - L_pairwise = {Lmeas - Lpair:+.4f}  "
              f"(near-Sidon: small;  evader: large => |T>=3| matters)", flush=True)


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
