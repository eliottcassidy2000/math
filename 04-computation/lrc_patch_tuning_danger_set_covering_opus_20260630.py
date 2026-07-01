"""
PATCH TUNING via the danger-set covering. M(S)<=1/n  <=>  union_{v in S} D_v = [0,1),
D_v = {t: ||vt||<=1/n} (v intervals of half-width 1/(vn) around j/v).
AP {1..n-1} TILES [0,1). Sporadic = retile: remove D_k, patch hole H_k = D_k minus union_{v!=k} D_v by D_g.
TEST: (1) AP covers; (2) the hole H_k; (3) which g patch it -- is g BOUNDED? (=> bounded-speed => finiteness).
"""
from fractions import Fraction
def frac(x): x=x%1.0; return min(x,1-x)
def indanger(v,t,n): return frac(v*t) <= 1.0/n + 1e-12
def covered(S,t,n): return any(indanger(v,t,n) for v in S)
def AP(n): return list(range(1,n))
# grid fine enough to resolve widths ~1/(n*maxv)
def hole_and_patches(n,k,gmax):
    base=AP(n); others=[v for v in base if v!=k]
    G=40*n*n  # grid
    # hole H_k = points in D_k not covered by others
    holepts=[i for i in range(G) if indanger(k,i/G,n) and not covered(others,i/G,n)]
    holefrac=len(holepts)/G
    # which g patches: D_g covers all hole points
    patches=[]
    for g in range(1,gmax+1):
        if g in base: continue
        if all(indanger(g,i/G,n) for i in holepts): patches.append(g)
    return holefrac, patches
print("Danger-set covering: AP tiles [0,1); remove D_k, find hole H_k and its patches g (bounded?):")
for (n,k) in [(6,2),(8,6),(14,12),(14,2),(14,6),(10,5)]:
    hf,patches=hole_and_patches(n,k,20*n)
    print(f"  n={n}, remove k={k}: |H_k|={hf:.4f}; patches g (<=20n) = {patches}  {'[g=2k='+str(2*k)+']' if 2*k in patches else ''}")
print()
print("=> if patches are always bounded (small, ~2k), single-swap sporadics have BOUNDED speed -> the")
print("   single-swap tight locus is a FINITE search. Extending: multi-patch (remove/add several).")
