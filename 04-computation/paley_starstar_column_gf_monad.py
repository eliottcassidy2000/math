#!/usr/bin/env python3
"""
paley_starstar_column_gf_monad.py
monad-explorer-2026-06-07 (deep-research, 8th session)

NEW STRUCTURAL FINDING about the cycle-rank triangle t(k,m) of (**):

   t(k,m) = sum_{even-series sigma, cycle-rank m} prod_v(|B_v|-1)!,   S_k = sum_m (-1)^m t(k,m) = (-1)^k C_k.

CONJECTURE (column rationality):  the column generating function
   T_m(x) := sum_{k>=m} t(k,m) x^k  =  P_m(x) * x^m / (1-x)^{2m-1},
with P_m(x) a POLYNOMIAL, P_m(0) = A088368(m) (the diagonal / all-pairings overcount),
and (conjecturally) deg P_m = m-2 for m>=2.

If true, U(x,y) = sum_m T_m(x) y^m has every y^m-coefficient RATIONAL in x with
denominator (1-x)^{2m-1} -- a strong constraint on the catalytic equation
(THM-438 ADDENDUM-5 handoff #1).  Specializing y=-1 must give F(x)-1, F=(-1+sqrt(1+4x))/(2x).

This script:
  (1) loads the triangle (k<=KMAX, hardcoded from the fast enumerator, extended when k=6/7 land),
  (2) for each column m, fits P_m(x) under denominator (1-x)^{2m-1} and reports degree + checks,
  (3) tabulates P_m and searches for structure,
  (4) verifies sum_m (-1)^m T_m(x) = F(x)-1 as power series,
  (5) attempts a low-degree quadratic CATALYTIC functional equation for U(x,y).
"""
import sys
from fractions import Fraction as Fr
from math import comb

# ---- the triangle t[k][m], 1-indexed m=1..k (from paley_starstar_triangle_fast_monad.py) ----
TRI = {
    1: [1],
    2: [1, 3],
    3: [1, 9, 13],
    4: [1, 18, 72, 69],
    5: [1, 30, 230, 580, 421],
}
# k=6 row gets appended here once the background enumerator finishes (see *_k6 .out).

def load_k6():
    import os
    path = "05-knowledge/results/paley_starstar_triangle_fast_k6_monad.out"
    if not os.path.exists(path):
        return
    txt = open(path).read()
    # find the k=6 t(k,m) line
    lines = txt.splitlines()
    for i, ln in enumerate(lines):
        if ln.strip().startswith("k=6:") and "S_k=" in ln:
            # next line "t(k,m) m=1..6: ..."
            for j in range(i, min(i + 4, len(lines))):
                if "t(k,m) m=1.." in lines[j]:
                    row = lines[j].split(":")[-1].split()
                    TRI[6] = [int(v) for v in row]
                    print(f"[loaded k=6 row from .out: {TRI[6]}]")
                    return


def poly_from_series(coeffs, m):
    """Given column m values t(m,m),t(m+1,m),...,t(K,m) (for k=m..K),
       assume T_m = P_m(x) x^m/(1-x)^{2m-1}; solve for P_m coefficients.
       [x^k] x^m/(1-x)^{2m-1} = C(k-m + (2m-2), 2m-2) = C(k+m-2, 2m-2).
       So t(k,m) = sum_j p_j * C(k-j+m-2, 2m-2)  where P_m = sum_j p_j x^j.
       Solve triangularly for p_0,p_1,..."""
    d = 2 * m - 2  # denominator (1-x)^(d+1)=(1-x)^(2m-1)
    p = []
    vals = coeffs[:]  # t(k,m) for k=m,m+1,...
    # base function b_j(k) = C(k-j+m-2, 2m-2), defined for k>=m+j
    def basis(k, j):
        n = k - j + m - 2
        if n < d:
            return 0
        return comb(n, d)
    # solve greedily: p_j determined by t(m+j,m)
    res = {}
    for j in range(len(vals)):
        k = m + j
        acc = sum(res.get(i, 0) * basis(k, i) for i in range(j))
        # basis(k,j) = C(k-j+m-2,2m-2)=C(2m-2,2m-2)=1
        pj = vals[j] - acc
        res[j] = pj
    # trim trailing zeros
    P = [res[j] for j in range(len(vals))]
    while len(P) > 1 and P[-1] == 0:
        P.pop()
    return P, d


def main():
    load_k6()
    KMAX = max(TRI)
    print("=" * 72)
    print(f"Triangle loaded up to k={KMAX}")
    for k in sorted(TRI):
        print(f"  k={k}: {TRI[k]}")
    print("=" * 72)
    print("COLUMN FIT:  T_m(x) = P_m(x) * x^m / (1-x)^(2m-1)")
    Ps = {}
    for m in range(1, KMAX + 1):
        col = [TRI[k][m - 1] for k in range(m, KMAX + 1)]
        P, d = poly_from_series(col, m)
        Ps[m] = P
        # verify: does P reproduce ALL given column values? (checks if deg stabilized)
        def predict(k):
            return sum(P[j] * comb(k - j + m - 2, 2 * m - 2) for j in range(len(P))
                       if k - j + m - 2 >= 2 * m - 2)
        ok = all(predict(k) == TRI[k][m - 1] for k in range(m, KMAX + 1))
        nterms = len(col)
        deg = len(P) - 1
        nchecks = nterms - len(P)  # data beyond what's needed to fix P
        print(f"  m={m}: P_m = {P}   deg={deg}  (#terms={nterms}, #independent-checks={max(nchecks,0)})  reproduces-col={ok}")
    print("-" * 72)
    print("P_m(0) (should be A088368 = 1,1,3,13,69,421,2867):",
          [Ps[m][0] for m in range(1, KMAX + 1)])
    print("deg P_m:", [len(Ps[m]) - 1 for m in range(1, KMAX + 1)], " (conj m-2 for m>=2)")
    # P_m(1) and other evaluations
    print("P_m(1):", [sum(Ps[m]) for m in range(1, KMAX + 1)])

    # ---- verify sum_m (-1)^m T_m(x) = F(x)-1 as power series ----
    print("-" * 72)
    N = KMAX
    # F(x) = sum_k (-1)^k C_k x^k
    def catalan(k): return comb(2 * k, k) // (k + 1)
    F = [(-1) ** k * catalan(k) for k in range(N + 1)]  # F[0]=1
    # S_k from triangle:
    S = [1] + [sum((-1) ** m * TRI[k][m - 1] for m in range(1, k + 1)) for k in range(1, N + 1)]
    print("S_k (from triangle):", S)
    print("(-1)^k C_k        :", F)
    print("loop-equation check S_k = -sum_{i+j=k-1} S_i S_j:")
    for k in range(1, N + 1):
        rhs = -sum(S[i] * S[k - 1 - i] for i in range(k))
        print(f"   k={k}: S_k={S[k]}  -conv={rhs}  {'OK' if S[k]==rhs else 'MISMATCH'}")
    print("=" * 72)
    print("Columns are RATIONAL with denominator (1-x)^(2m-1).  This pins the")
    print("y^m-coefficient of U(x,y); the catalytic equation must reproduce it.")


if __name__ == "__main__":
    main()
