#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_threadA_closed_form_kpswf5.py  (kind-pasteur 2026-06-21, THREAD A part 2)

Find the UNIFORM closed form of the reduced f_P(a,b) = S_P/P table on
a,b in {1..(P-1)/2}, valid for all odd primes P.

From kpswf5 part 1, the tables (rowdef = S/P) are:

  P=5:  a\b  1   2          P=7:  a\b 1   2   3
        1    8  12                1   12  20  24
        2   12  16                2   20  32  36
                                  3   24  36  44

  P=11: a\b 1   2   3   4   5     P=13: a\b 1  2   3   4   5   6
        1   20  36  48  56  60          1  24 44  60  72  80  84
        2   36  64  84  96 100          2  44 80 108 128 140 144
        3   48  84 108 124 136          3  60 108 144 168 188 200
        4   56  96 124 148 156          4  72 128 168 204 228 236
        5   60 100 136 156 168          5  80 140 188 228 248 264
                                        6  84 144 200 236 264 276

The DIAGONAL (a=a): max corner f(h,h)=maxS/P with h=(P-1)/2.
  P=5: 16, P=7: 44, P=11: 168, P=13: 276.  These are 4*(1,11,42,69)=4*{4,11,42,69}? no.
  Actually maxS/4 = 4,11,42,69 for P=5,7,11,13.  And maxS = (P^2-1)*something?
    P=5: (25-1)=24 -> 16 = 24*2/3 ; P=7: 48 -> 44 = 48*11/12 ; hmm.

Strategy: this is a discrete L1 discrepancy of a 1D sawtooth.  Re-derive f_P(a,b)
analytically from the cov closed form, then fit.  Key object:
  rowdef(p,q) = sum_j | P*c_0j - p*q |, c_0j = (p//P)*q + #{z in [qj,q(j+1)): z%P < p%P}.
  Subtract the integer part: let r=p%P (=a up to sign), then c_0j - (p//P)*q =
  #{z in [qj,q(j+1)): z%P < r}.  And p*q/P = (p//P)*q + r*q/P.  So
  P*c_0j - p*q = P*#{z in window: z%P<r} - r*q.
  i.e. rowdef = sum_j | P * N_j(r) - r*q |, N_j(r)=#{z in [qj,q(j+1)): z%P<r}.
This depends on (r, q%P) only.  We compute it as a clean function of (r mod P, q mod P)
and reduce by ||.||.  Then fit a closed form.
"""
from math import gcd
from fractions import Fraction as Fr


def rowdef_residue(r, qm, P):
    """rowdef as a function of (r=p%P, qm=q%P), via the residue-only structure.
    We use the exact analytic form: choose any (p,q) with p%P=r, q%P=qm, coprime,
    in-window, large enough that the period-P sawtooth is fully developed.
    Then rowdef = sum_j |P*N_j - r*q| with N_j = #{z in [qj,q(j+1)): z%P<r}."""
    # pick concrete representatives: q = qm + P*Q, p = r + P*Q2 ensuring coprime, p>q, p<Pq.
    # Use moderately large multipliers so the count is "developed".
    for Qm in range(1, 40):
        q = qm + P * Qm if qm != 0 else P * Qm
        if q < 1:
            continue
        for Pm in range(1, 40):
            p = r + P * Pm if r != 0 else P * Pm
            if p <= q or p >= P * q:
                continue
            if gcd(p, q) != 1:
                continue
            base = p // P
            rr = p % P
            rd = 0
            for j in range(P):
                Nj = sum(1 for z in range(q * j, q * (j + 1)) if z % P < rr)
                rd += abs(P * Nj - rr * q)
            return rd
    return None


def norm_P(x, P):
    m = x % P
    return min(m, P - m)


def candidate_forms(P):
    """Build the reduced f(a,b) table and TEST several uniform closed-form candidates."""
    h = (P - 1) // 2
    f = {}
    for a in range(1, h + 1):
        for b in range(1, h + 1):
            rd = rowdef_residue(a, b, P)
            f[(a, b)] = rd
    return f, h


def main():
    print("Searching for the UNIFORM closed form f_P(a,b), a,b in 1..(P-1)/2")
    print("=" * 74)
    primes = [3, 5, 7, 11, 13, 17, 19]
    F = {}
    for P in primes:
        f, h = candidate_forms(P)
        F[P] = (f, h)
        print(f"\nP={P} (h={h}): f-table")
        for a in range(1, h + 1):
            print("   " + " ".join(f"{f[(a,b)]:4d}" for b in range(1, h + 1)))

    # ---- candidate 1: f = a*b*(P - a*b/?)/?  -- product-sawtooth (triangle) form ----
    # The 1D L1 discrepancy of {k*alpha} for k=0..M-1 vs uniform is ~ a*(P-a)/P type.
    # Sawtooth sum_j |P N_j - r q| with r=a is a known "Dedekind-sum-like" object.
    # Hypothesis (Dedekind/sawtooth): f(a,b) = 2*a*b*(P-a)*(P-b)/P^2 ? test.
    print("\n" + "=" * 74)
    print("Testing closed-form hypotheses (exact rational, must match ALL entries):")

    def test(name, fn):
        allok = True
        for P in primes:
            f, h = F[P]
            ok = True
            for a in range(1, h + 1):
                for b in range(1, h + 1):
                    val = fn(a, b, P)
                    if val != f[(a, b)]:
                        ok = False
            print(f"   {name:52s} P={P:2d}: {'MATCH' if ok else 'no'}")
            allok = allok and ok
        print(f"   --> {name}: {'UNIFORM MATCH on all P' if allok else 'fails'}")
        return allok

    # H1: 2 a b (P-a)(P-b)/P^2  (separable Dedekind product)
    test("H1: 2ab(P-a)(P-b)/P^2",
         lambda a, b, P: Fr(2 * a * b * (P - a) * (P - b), P * P))
    # H2: 4 a b (P-a)(P-b)/P^2
    test("H2: 4ab(P-a)(P-b)/P^2",
         lambda a, b, P: Fr(4 * a * b * (P - a) * (P - b), P * P))
    # H3: a b (P - a) (P - b) * 4 / (P*(P-1)) ... try several normalizations later
    # Let's just compute the ratio f/(ab(P-a)(P-b)) to find the right constant.
    print("\n  ratio f(a,b) / [a b (P-a)(P-b)] (look for constant):")
    for P in primes:
        f, h = F[P]
        ratios = set()
        for a in range(1, h + 1):
            for b in range(1, h + 1):
                denom = a * b * (P - a) * (P - b)
                ratios.add(Fr(f[(a, b)], denom))
        print(f"   P={P:2d}: distinct ratios = {sorted(ratios)}")

    # H4: maybe f separates as g(a)*g(b)+something? check rank-1 structure of f-table.
    print("\n  is f(a,b) rank-1 (separable) as a matrix?  check f(a,b)*f(1,1) vs f(a,1)*f(1,b):")
    for P in primes:
        f, h = F[P]
        sep = True
        for a in range(1, h + 1):
            for b in range(1, h + 1):
                if f[(a, b)] * f[(1, 1)] != f[(a, 1)] * f[(1, b)]:
                    sep = False
        print(f"   P={P:2d}: separable(rank-1)={sep}")

    # H5: the marginal f(a,1) = the 'edge' values.  Tabulate and fit.
    print("\n  edge marginal f(a,1) vs a (the b=1 column):")
    for P in primes:
        f, h = F[P]
        col = [f[(a, 1)] for a in range(1, h + 1)]
        # compare to a*(2P-?)...
        print(f"   P={P:2d}: f(a,1)= {col}")
        # try f(a,1) = a*(P-a)*c ; print f(a,1)/(a(P-a))
        rr = [Fr(f[(a, 1)], a * (P - a)) for a in range(1, h + 1)]
        print(f"        f(a,1)/[a(P-a)] = {rr}")


if __name__ == "__main__":
    main()
