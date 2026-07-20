#!/usr/bin/env python3
"""
klein-2026-07-20-S359 -- TNC BY EXACT ELIMINATION, span by span: the induction whose base
cases are Groebner computations and whose payoff (via the Gamma bridge) is GMC(2) on bounded
CHARGE SPAN with UNBOUNDED coefficient degree.

WHERE THIS SITS.  The assembled chain is  TNC => NC2 => GMC(2):
  * NC2 => GMC(2)                      PROVED (THM-1510 SS C, weight counting).
  * TNC => NC2                         the Gamma bridge (klein-S351): psi_m is a polynomial
                                       in r whose degree GROWS with m and whose TOP
                                       coefficient is the toral quantity, so the k! = Gamma
                                       weights make the top term dominate.
  * TNC at extreme weight +-1          PROVED (THM-1530, Lagrange-Buermann, exact).
  * TNC at M, N >= 2                   THE REMAINING LINK.  Attacked here.

WHY ELIMINATION IS THE RIGHT INDUCTION HERE.  TNC is a statement about a Laurent polynomial
with CONSTANT coefficients -- Lam = sum_{q=-M}^{N} c_q u^q -- so for each span (M,N) it is a
FINITE system in M+N+1 unknowns:  CT(Lam^m) = 0 for m = 1..K.  Unlike THM-1590, no coefficient
DEGREE is bounded: the span alone bounds the problem.  So proving TNC span-by-span and feeding
the Gamma bridge gives GMC(2) for that charge span at EVERY coefficient degree.

THE CORRECT TEST, and it is NOT the one from THM-1590.  Here the degenerate locus is not the
origin: it is "the span collapses", i.e. c_{-M} = 0 or c_N = 0.  So TNC for span (M,N) says

        V(I)  is contained in  {c_{-M} = 0} u {c_N = 0},

which by the Nullstellensatz is exactly:  (c_{-M} * c_N)^k lies in I for some k.
That is the ideal-membership question computed below.
"""
import sympy as sp
from sympy import groebner
import itertools, math

def CT_pow(cs, M, N, m):
    """CT(Lam^m) where Lam = sum_{q=-M}^{N} cs[q+M] u^q -- exact symbolic constant term"""
    cur = {0: sp.Integer(1)}
    for _ in range(m):
        nxt = {}
        for e1, v1 in cur.items():
            for j, cq in enumerate(cs):
                q = j - M
                if cq == 0: continue
                nxt[e1 + q] = nxt.get(e1 + q, 0) + v1 * cq
        cur = nxt
    return sp.expand(cur.get(0, sp.Integer(0)))

print("=" * 86)
print("TNC BY EXACT ELIMINATION, span by span")
print("=" * 86)
print("  test: is (c_{-M} * c_N)^k in the ideal <CT(Lam^m) : m=1..K> ?")
print("  YES  =>  every solution has a collapsed span  =>  TNC HOLDS for that (M,N).\n")
print(f"{'(M,N)':>8} {'unknowns':>9} {'m range':>10} {'k with (c-M cN)^k in I':>24} {'TNC?':>8}")
import sys
for (M, N) in [(1,1), (1,2), (2,1), (2,2), (3,1), (1,3), (2,3), (3,2), (3,3)]:
    d = M + N
    cs = list(sp.symbols(f'c0:{d+1}'))
    # torus normalisation: CT is invariant under u -> lam u (c_q -> lam^q c_q) and Lam -> mu Lam.
    # Use it to set the two extreme coefficients to 1 -- legitimate on the open set where both
    # are nonzero, which is exactly the set TNC claims is empty.
    csn = cs[:]; csn[0] = sp.Integer(1); csn[d] = sp.Integer(1)
    unk = [c for c in cs[1:d]]
    m0 = (M + N) // math.gcd(M, N)            # first m admitting a weight-0 monomial
    K = m0 + len(unk) + 1                     # start lean; escalate below only if needed
    eqs = []
    for m in range(1, K + 1):
        e = CT_pow(csn, M, N, m)
        if e != 0: eqs.append(e)
    if not unk:
        # no free unknowns: TNC holds iff some equation is a nonzero constant
        nz = [e for e in eqs if e.is_number and e != 0]
        print(f"{str((M,N)):>8} {0:>9} {f'{1}..{K}':>10} {'(no free unknowns)':>24} "
              f"{('YES' if nz else 'CHECK'):>8}")
        continue
    inconsistent = None; kk = 'not reached'; Kused = K
    for Ktry in (K, K + 2, K + 4):
        eqs = []
        for m in range(1, Ktry + 1):
            e = CT_pow(csn, M, N, m)
            if e != 0: eqs.append(e)
        try:
            G = groebner(eqs, *unk, order='grevlex')
            gens = list(G.exprs)
            inconsistent = (len(gens) == 1 and gens[0] == 1)
            kk = '<1>  INCONSISTENT' if inconsistent else f'basis size {len(gens)}'
            Kused = Ktry
            if inconsistent: break
        except Exception as ex:
            kk = f'groebner failed'; break
    print(f"{str((M,N)):>8} {len(unk):>9} {f'{1}..{Kused}':>10} {kk:>24} "
          f"{('YES' if inconsistent else ('OPEN' if inconsistent is not None else '?')):>8}")
    sys.stdout.flush()
print("""
  NOTE ON THE NORMALISATION.  CT(Lam^m) is invariant under Lam -> mu*Lam only up to mu^m, and
  under u -> lam*u exactly; together these let me set c_{-M} = c_N = 1 on the open set where
  BOTH are nonzero.  TNC for span (M,N) is precisely the claim that this open set contains no
  solution -- i.e. that the normalised system is INCONSISTENT, Groebner basis <1>.  So the
  test above is the honest one: <1> means TNC holds for that span, anything else means the
  span is not yet settled by these K moments.
""")
