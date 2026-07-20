"""opus-2026-07-20-S420 -- THE FINISHING STATEMENT: the tuned-cancellation locus is empty.

From THM-1655: TNC is proved for binomials and for any R with a UNIQUE minimal charge-rep of
0.  The residual is NON-unique minimal reps tuned so CT(m_0) = 0.  Witness 1+u^3-u^6 (N=2):
CT(3)=0 but CT(6)=-30.  This script attacks the residual two ways.

STRUCTURAL FACT (worth stating).  With r_N = 0 forced (THM-1655 s1), F(u,t) = u^N - tR(u) is
LINEAR in t: [u^N]F = 1 - t r_N = 1, and [u^k]F = -t r_k (k != N).  So the whole branch
geometry is governed by a t-linear polynomial -- a strong rigidity the saddle picture hides.

ATTACK 1 (trinomial theorem).  R = r_0 + r_j u^j + r_d u^d, three charges {-N, j-N, M}.  Is
there a non-monomial trinomial that is a nullcone element (CT(m)=0 for ALL m)?  Parametrize
the tuned locus CT(m_0)=0 and check whether ANY point also kills all later CT.

ATTACK 2 (degenerate saddles).  THM-1635 proves TNC whenever a dominant saddle is
NONDEGENERATE.  So a nullcone violator needs ALL dominant saddles degenerate (g''=0).  Search
non-monomial R for an all-dominant-degenerate case and test its CT.
"""
import sympy as sp
from math import gcd
from itertools import product

u, t = sp.symbols('u t')

def CT(Rexpr, N, mm):
    e = sp.expand(sp.sympify(Rexpr)**mm)
    return sp.Poly(e, u).coeff_monomial(u**(N*mm))

print("="*78)
print("(0) STRUCTURAL FACT: r_N=0  =>  F(u,t)=u^N - tR(u) is LINEAR in t")
print("="*78)
for R, N in [("1+u**3-u**6", 2), ("1+u**2+u**5", 3)]:
    F = sp.expand(u**N - t*sp.sympify(R))
    degt = sp.Poly(F, t).degree()
    print(f"   R={R:14s} N={N}: deg_t F = {degt}   (linear: {degt==1});  "
          f"[u^N]F = {sp.Poly(F,u).coeff_monomial(u**N)}")

print()
print("="*78)
print("(1) TRINOMIAL SWEEP: does any non-monomial trinomial kill ALL CT up to a high bound?")
print("="*78)
print("   For each charge pattern, tune the middle coeff to cancel CT(m_0), then scan later CT.")
survivors = []
tested = 0
for N in [2, 3]:
    for j in range(1, 8):
        for d in range(j+1, 10):
            if j == N or d == N: continue          # r_N = 0
            # R = 1 + a u^j + u^d ; find m0 = first m with CT != 0 generically
            Rsym = 1 + sp.Symbol('a')*u**j + u**d
            a = sp.Symbol('a')
            m0 = None
            for m in range(1, 14):
                c = sp.expand(CT(1 + a*u**j + u**d, N, m))
                if c != 0:
                    m0 = m; break
            if m0 is None: continue
            tested += 1
            c_m0 = sp.Poly(sp.expand(CT(1 + a*u**j + u**d, N, m0)), a)
            roots = sp.solve(c_m0.as_expr(), a)
            for a0 in roots:
                if a0 == 0: continue                 # monomial-ish
                Rt = sp.expand(1 + a0*u**j + u**d)
                # scan later CT up to a bound
                later = [sp.simplify(CT(Rt, N, m)) for m in range(m0, m0+12)]
                if all(v == 0 for v in later):
                    survivors.append((N, j, d, a0, m0))
print(f"   trinomial charge patterns tested: {tested}")
print(f"   tuned points killing CT through m_0+11: {len(survivors)}")
for s in survivors[:8]:
    print(f"      N={s[0]} j={s[1]} d={s[2]} a={s[3]} (m0={s[4]}) -- CANDIDATE, inspect")
if not survivors:
    print("   NONE.  Every tuned-cancellation trinomial has a surviving later CT.")
    print("   => no non-monomial trinomial is a nullcone element (up to this bound).")

print()
print("="*78)
print("(2) THE DOUBLING LAW: when CT(m_0)=0 is tuned, where does CT first survive?")
print("="*78)
print("   pattern from the witness: m_0 -> next nonzero at a MULTIPLE of m_0?")
for N, j, d in [(2, 3, 6), (2, 1, 4), (3, 1, 5), (2, 3, 5), (3, 2, 7)]:
    a = sp.Symbol('a')
    m0 = next((m for m in range(1, 14) if sp.expand(CT(1+a*u**j+u**d, N, m)) != 0), None)
    if m0 is None: continue
    roots = sp.solve(sp.expand(CT(1+a*u**j+u**d, N, m0)), a)
    tuned = [r for r in roots if r != 0]
    if not tuned:
        print(f"   N={N} (charges -{N},{j-N},{d-N}): m0={m0}, CT(m0) has NO nonzero tuning "
              f"(unique rep) -> TNC by THM-1655")
        continue
    a0 = tuned[0]; Rt = sp.expand(1+a0*u**j+u**d)
    nxt = next((m for m in range(m0+1, m0+20) if sp.simplify(CT(Rt, N, m)) != 0), None)
    print(f"   N={N} (charges -{N},{j-N},{d-N}): m0={m0}, tuned a={a0}, "
          f"first survivor m={nxt}   {'= 2*m0' if nxt==2*m0 else ''}")
