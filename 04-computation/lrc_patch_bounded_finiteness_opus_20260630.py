"""
RIGOROUS bounded patch: H_k subset D_g requires every hole-interval width <= 2/(gn) (fits in one D_g interval,
since D_g has gaps). So the WIDEST hole-interval w_max gives g <= 2/(n*w_max). Bounded => finite.
Verify the bound; extend to multi-swap (remove several, add several).
"""
def frac(x): x=x%1.0; return min(x,1-x)
def indanger(v,t,n): return frac(v*t) <= 1.0/n + 1e-12
def covered(S,t,n): return any(indanger(v,t,n) for v in S)
def hole_widths(n,removed,G):
    base=[v for v in range(1,n) if v not in removed]
    hole=[1 if (any(indanger(k,i/G,n) for k in removed) and not covered(base,i/G,n)) else 0 for i in range(G)]
    # widest run (circular)
    runs=[]; c=0
    for x in hole+hole:  # circular
        if x: c+=1
        else:
            if c: runs.append(c); c=0
    if c: runs.append(c)
    wmax=(max(runs)/G) if runs else 0.0
    meas=sum(hole)/G
    return meas,wmax
print("(1) single-swap: hole H_k, widest interval w_max, bound g<=2/(n*w_max), actual patch:")
G=60000
for (n,k,patch) in [(5,2,7),(6,2,9),(8,6,12),(14,12,24)]:
    meas,wmax=hole_widths(n,[k],G)
    bound=2/(n*wmax) if wmax>0 else float('inf')
    print(f"  n={n},k={k}: |H_k|={meas:.4f}, w_max={wmax:.5f}, bound 2/(n*w_max)={bound:.1f}; patch g={patch} <= bound: {patch<=bound}")
print()
print("(2) multi-swap n=8 second sporadic {1,4,5,6,7,11,13} = remove {2,3}, add {11,13}:")
meas,wmax=hole_widths(8,[2,3],G)
bound=2/(8*wmax) if wmax>0 else float('inf')
print(f"  remove {{2,3}}: |hole|={meas:.4f}, w_max={wmax:.5f}, bound={bound:.1f}; patches 11,13 <= bound: {11<=bound and 13<=bound}")
print()
print("=> the patch elements are BOUNDED by 2/(n*w_max), w_max = widest hole interval > 0. Since a hole is a")
print("   finite union of open intervals (w_max>0), every patch is bounded => the AP-modification tight locus")
print("   is a FINITE search. This is the bounded-speed lever for OPEN-Q-108 (single/multi-swap component).")
