#!/usr/bin/env python3
"""
paley_starstar_core_decomp_monad.py
monad-explorer-2026-06-07 (deep-research, 9th session)

CORE/SUBDIVISION DECOMPOSITION of the cycle-rank column GF (proves column rationality).

Claim (THM-438 ADDENDUM-7):
   T_m(x) = sum_k t(k,m) x^k = sum_{e=m}^{2m-1} R(m,e) * (x/(1-x))^e
where R(m,e) = sum over REDUCED even-series patterns (every line length exactly 2,
hence k=e) of cycle-rank m with exactly e lines, weighted by W=prod(|B_v|-1)!.

Consequences (all checked here):
   * pole order = max e = 2m-1  (Eulerian: <=2 odd-deg vertices; reduced core: min deg>=3
     => V_core<=m => e=V_core+m-1<=2m-1).
   * R(m,m) = A088368(m)  (single-vertex m-loop "rosette"; minimal e=m, V_core=1).
   * P_m(x) = sum_e R(m,e) x^{e-m} (1-x)^{2m-1-e}, a polynomial of degree<=m-1.
   * deg P_m = m-2  <=>  p_{m-1} = -sum_e (-1)^e R(m,e) = 0  <=>  Q_m(-1)=0,
     where Q_m(t)=sum_e R(m,e) t^e = t^m (1+t)^{m-1} P_m(t/(1+t)).
   * lead P_m = p_{m-2} = sum_e (-1)^e (2m-1-e) R(m,e) = 2^m - 1.

This script BRUTE-extracts R(m,e) (from k=e enumeration) to VALIDATE the decomposition,
then reconstructs P_m and the full column, cross-checking the known triangle.
"""
import sys
from math import comb, factorial
from collections import defaultdict

# import the validated analyzer from the fast enumerator
sys.path.insert(0, "04-computation")
from paley_starstar_triangle_fast_monad import rgs_iter, analyze, catalan


def brute_R(EMAX):
    """R[(m,e)] = sum of W over even-series patterns at k=e with cycle-rank m and e lines.
       (At k=e, #lines=e forces every line length exactly 2 = 'reduced'.)
       We bucket ALL even-series patterns at each k by (m,e,W); R(m,e)=bucket at k=e."""
    R = defaultdict(int)
    full = {}  # full[k] = dict (m,e)->W-sum  (for cross-checks / #lines distribution)
    for k in range(1, EMAX + 1):
        L = 2 * k
        n = L + 1
        tab = defaultdict(int)
        for a in rgs_iter(n):
            res = analyze(a, L)
            if res is None:
                continue
            m, W, e = res
            tab[(m, e)] += W
        full[k] = tab
        # R(m,e) is exactly the e=k slice
        for (m, e), w in tab.items():
            if e == k:
                R[(m, e)] = w
    return R, full


def reconstruct_column(R, m):
    """P_m(x) = sum_{e=m}^{2m-1} R(m,e) x^{e-m} (1-x)^{2m-1-e}, return coeff list."""
    deg_bound = m - 1
    P = [0] * (deg_bound + 1)
    for e in range(m, 2 * m):
        Rme = R.get((m, e), 0)
        if Rme == 0:
            continue
        # x^{e-m} * (1-x)^{2m-1-e}
        p = 2 * m - 1 - e  # power of (1-x)
        for i in range(p + 1):
            c = Rme * comb(p, i) * ((-1) ** i)
            P[(e - m) + i] += c
    while len(P) > 1 and P[-1] == 0:
        P.pop()
    return P


def column_value(P, m, k):
    """t(k,m) from P_m via x^k coeff of P_m(x) x^m/(1-x)^{2m-1}."""
    s = 0
    for j, pj in enumerate(P):
        nn = k - j + m - 2
        if nn >= 2 * m - 2:
            s += pj * comb(nn, 2 * m - 2)
    return s


def main():
    EMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    A088368 = {1: 1, 2: 3, 3: 13, 4: 69, 5: 421, 6: 2867, 7: 22417}
    TRI = {1:[1],2:[1,3],3:[1,9,13],4:[1,18,72,69],5:[1,30,230,580,421],
           6:[1,45,560,2626,4845,2867]}
    print("=" * 74)
    print("CORE/SUBDIVISION DECOMPOSITION  T_m = sum_e R(m,e) (x/(1-x))^e")
    print("=" * 74)
    R, full = brute_R(EMAX)

    print("\nR(m,e) table (brute, k=e slice):  e = m .. 2m-1")
    mmax = EMAX
    for m in range(1, mmax + 1):
        row = []
        for e in range(m, 2 * m):
            if e <= EMAX:
                row.append((e, R.get((m, e), "?")))
        print(f"  m={m}: " + "  ".join(f"R(m,{e})={v}" for e, v in row))

    print("\nStructural checks:")
    for m in range(1, mmax + 1):
        # only do full checks if all e=m..2m-1 are within reach (2m-1<=EMAX)
        have_all = (2 * m - 1) <= EMAX
        diag = R.get((m, m), None)
        ok_diag = (diag == A088368.get(m))
        line = f"  m={m}: R(m,m)={diag}  A088368={A088368.get(m)}  diag-match={ok_diag}"
        if have_all:
            alt = sum((-1) ** e * R.get((m, e), 0) for e in range(m, 2 * m))
            P = reconstruct_column(R, m)
            lead = P[-1] if P else 0
            line += (f"  | sum_e(-1)^e R = {alt} (want 0)"
                     f"  | P_m={P}  deg={len(P)-1}(want {m-2})"
                     f"  lead={lead}(want {2**m-1})")
        print(line)

    print("\nReconstruct columns from R and cross-check known triangle:")
    for m in range(1, mmax + 1):
        if (2 * m - 1) > EMAX:
            print(f"  m={m}: need e up to {2*m-1}>EMAX={EMAX} -- partial")
            continue
        P = reconstruct_column(R, m)
        # predict t(k,m) for k=m..6 and compare to TRI
        preds = []
        ok = True
        for k in range(m, 7):
            pv = column_value(P, m, k)
            preds.append(pv)
            if k in TRI and m <= len(TRI[k]):
                if pv != TRI[k][m - 1]:
                    ok = False
        print(f"  m={m}: P_m={P}  predicts t(k,m) k={m}..6: {preds}  match-triangle={ok}")


if __name__ == "__main__":
    main()
