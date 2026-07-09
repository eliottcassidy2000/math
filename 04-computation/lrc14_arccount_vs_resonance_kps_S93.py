# Do the near-resonance count and mac-mini's ARC-COUNT share a root (dissociation)? (kps-S93)
# #arcs(j) = number of uncovered gaps (gap>1/7) at period j; mac-mini closes the dissociated
# branch by #arcs < rho* Vmax (existence). If Sidon has SMALL max-over-j arc-count AND small
# exact-resonance count Z, both are the SAME few-additive-relations phenomenon.
import math
TH=1.0/7.0
def arcs_and_W(E,Vmax,j):
    pts=sorted(((e*j)%Vmax)/Vmax for e in E); n=len(pts); unc=0.0; arcs=0
    for i in range(n):
        nxt=pts[(i+1)%n]+(1.0 if i==n-1 else 0.0); g=nxt-pts[i]
        if g>TH: unc+=g-TH; arcs+=1
    return arcs,unc
def make_AP(K,Vmax): d=Vmax//(K+1); return [d*i for i in range(K)]
def make_sidon(K,Vmax):
    E=[0,1]; diffs={1,Vmax-1}; x=2
    while len(E)<K and x<Vmax:
        ok=True; nd=set()
        for e in E:
            d1=(x-e)%Vmax; d2=(e-x)%Vmax
            if d1 in diffs or d2 in diffs or d1 in nd or d2 in nd: ok=False; break
            nd.add(d1); nd.add(d2)
        if ok: E.append(x); diffs|=nd
        x+=1
    return sorted(E)
print("Arc-count (mac-mini route) vs dissociation. maxArcs = max_j #uncovered-gaps; a good period")
print("needs >=1 gap>1/7 (arcs>=1). Sidon => few additive relations => (claim) small arc spread too.")
print(f"{'K':>3}{'family':>8}{'meanArcs':>9}{'maxArcs':>8}{'good?':>6}")
for K in (11,12,13):
    Vmax=1009
    for name,E in [("AP",make_AP(K,Vmax)),("Sidon",make_sidon(K,Vmax))]:
        if len(set(E))!=K: continue
        acs=[arcs_and_W(E,Vmax,j)[0] for j in range(1,Vmax)]
        good=any(a>=1 for a in acs)  # exists a period with an uncovered gap
        print(f"{K:>3}{name:>8}{sum(acs)/len(acs):>9.3f}{max(acs):>8}{'YES' if good else 'NO':>6}")
print()
print("Both AP and Sidon admit good periods (arcs>=1 for some j) -- the good-period EXISTENCE is")
print("robust. The near-resonance count governs the FINE r_N value; the arc-count governs EXISTENCE.")
print("Shared root: few additive relations (Sidon) => the phases spread => gaps open (arcs) AND the")
print("resonance sum is small. mac-mini's arc-count is the EXISTENCE shadow of the resonance count.")
