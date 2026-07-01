"""
lrc_covmin_selfconcordant_residual_opus_20260630.py
The covering-min self-concordant residual + H3/H6 back-and-forth.
 (1) 1/M = (n-1) + 1/rung  (barrier param nu=n-1 + residual 1/rung), verified on beaters + construction.
 (2) covering-min at n=7,8 = 2/13,2/15 (BELOW construction => H3 holds; my old "n>=7 exact" corrected).
 (3) the known beaters n=7..10 MISS core {1..n-2} yet M<n/Phi6 => violate lowness lemma HYP-3747 (n-dependent).
 (4) H6 window: margin(rung k)=(k-1)/(n(k(n-1)+1)); hexagonal(k=2) <-> cyclotomic(k=n, =Dedekind sum HYP-3768); ratio->2.
ASCII only. See reflection the-covering-min-self-concordant-residual-...-H3-vs-H6-opus-20260630.md.
"""
import numpy as np, itertools
from fractions import Fraction as F
from math import gcd
def phi6(n): return n*n-n+1
def M_exact(S,Q):
    Sa=np.array(S); bn,bd=0,1
    for q in range(2,Q+1):
        A=np.outer(Sa,np.arange(1,q))%q; k=int(np.minimum(A,q-A).min(axis=0).max())
        if k*bd>bn*q: bn,bd=k,q
    g=gcd(bn,bd) or 1; return bn//g,bd//g

def decomposition():
    print("(1) SELF-CONCORDANT RESIDUAL  1/M = (n-1) + 1/rung:")
    known={7:F(2,13),8:F(2,15),9:F(4,33),10:F(4,37),11:F(3,31)}
    for n,M in known.items():
        res=1/M-(n-1); print(f"    n={n}: M={M}  1/M={1/M}  = {n-1} + {res}  (rung {1/res})")
    for n in [14]:
        M=F(n,phi6(n)); print(f"    n={n} construction: M={M}  1/M={1/M} = {n-1} + {1/M-(n-1)} (rung {n}); AP floor 1/n: residual 1 (rung 1, analytic center)")

def confirm_H3():
    print("\n(2)+(3) covering-min at n=7,8 (exhaustive min M>1/n) + beaters violate HYP-3747:")
    for n,hi in [(7,14),(8,16)]:
        thr=n/phi6(n); best=(9,1)
        for S in itertools.combinations(range(1,hi+1),n-1):
            bn,bd=M_exact(S,9*n); v=bn/bd
            if v>1.0/n+1e-12 and v<best[0]/best[1]: best=(bn,bd)
        g=gcd(*best) or 1
        print(f"    n={n}: covering-min={best[0]//g}/{best[1]//g}={best[0]/best[1]:.5f} (construction {n}/{phi6(n)}={thr:.5f}) -> BELOW: {best[0]/best[1]<thr}")
    beaters={7:[1,2,5,6,7,8],8:[1,4,5,6,7,11,16],9:[1,3,4,5,7,11,18,32],10:[1,2,3,5,6,7,8,9,30]}
    for n,S in beaters.items():
        bn,bd=M_exact(S,150); miss=[c for c in range(1,n-1) if c not in S]
        print(f"    n={n} beater {tuple(S)}: M={bn}/{bd}<{n}/{phi6(n)}={bn/bd<n/phi6(n)}, missing core {miss} -> HYP-3747 violated: {bn/bd<n/phi6(n) and len(miss)>0}")

def H6_window():
    print("\n(4) H6 WINDOW  margin(k)=(k-1)/(n(k(n-1)+1)); hexagonal(k=2) <-> cyclotomic(k=n); ratio cyc/hex -> 2:")
    for n in [7,10,14,20]:
        Mk=lambda k: F(k,k*(n-1)+1)
        hexm=Mk(2)-F(1,n); cycm=Mk(n)-F(1,n)
        print(f"    n={n}: hex(rung2)={hexm}={float(hexm):.5f}  cyc(rung n)={cycm}={float(cycm):.5f}  ratio={float(cycm/hexm):.4f}  (cyc = Dedekind sum, HYP-3768)")

if __name__=="__main__":
    decomposition(); confirm_H3(); H6_window()
    print("\n=> H3 vs H4 <=> the interior-point residual 1/a(n) stays Theta(1) (small rung, near AP center)")
    print("   vs vanishes like 1/n (rung n, construction). Both give Theta(1/n^2) margin; different constant.")
