from fractions import Fraction as F
from itertools import combinations

def forbidden(v, delta):
    r=delta/v; out=[]
    for k in range(v):
        c=F(k,v); L=(c-r)%1; H=(c+r)%1
        if L<H: out.append((L,H))
        else: out+=[(L,F(1)),(F(0),H)]
    return merge(out)
def merge(iv):
    iv=sorted(iv); res=[]
    for l,h in iv:
        if res and l<=res[-1][1]: res[-1]=(res[-1][0],max(res[-1][1],h))
        else: res.append((l,h))
    return res

# LRC(14): 13 runners, delta=1/14. Residual = configs containing a multiple of 14 (the apex).
delta=F(1,14)
def lonely_set(V):
    U=forbidden(V[0],delta)
    for v in V[1:]: U=merge(U+forbidden(v,delta))
    # complement of U on [0,1)
    comp=[]; prev=F(0)
    for l,h in U:
        if l>prev: comp.append((prev,l))
        prev=max(prev,h)
    if prev<1: comp.append((prev,F(1)))
    return comp, sum(h-l for l,h in comp)

def cross_structure(V):
    # at each unit clock j (gcd(j,14)=1): which runners are TIGHT (v*j ≡ ±1 mod 14)?
    units=[j for j in range(1,14) if F(j,14).denominator==14 and __import__('math').gcd(j,14)==1]
    rows=[]
    for j in units:
        up=[v for v in V if (v*j)%14==1]    # at +1/14, needs eps>0 to stay safe
        dn=[v for v in V if (v*j)%14==13]    # at -1/14, needs eps<0
        rows.append((j,up,dn))
    return rows

configs = {
 "WALL {1..11,13,14} (all 6 units + apex 14)":[1,2,3,4,5,6,7,8,9,10,11,13,14],
 "AP {1..13} (no apex, classical tight)":     list(range(1,14)),
 "apex+partial {1,2,3,4,5,6,7,8,9,10,11,12,14}":[1,2,3,4,5,6,7,8,9,10,11,12,14],
 "bigger apex {1..11,13,28}":                 [1,2,3,4,5,6,7,8,9,10,11,13,28],
}
for name,V in configs.items():
    comp,p0=lonely_set(V)
    has_apex=any(v%14==0 for v in V)
    print(f"\n=== {name} ===")
    print(f"  apex present: {has_apex}   p0 (lonely measure) = {float(p0):.6f}   lonely-set nonempty: {len(comp)>0}")
    if comp:
        # is the lonely set near a unit clock j/14? show first few lonely intervals and nearest unit clock
        shown=comp[:6]
        print(f"  lonely intervals ({len(comp)} total), first few: "+", ".join(f"({float(l):.4f},{float(h):.4f})" for l,h in shown))
        # distance of each lonely interval midpoint to nearest j/14
        for l,h in comp[:4]:
            mid=(l+h)/2; nearest=min(range(15), key=lambda j: abs(float(mid)-j/14))
            print(f"     interval mid={float(mid):.4f} width={float(h-l):.5f}  nearest unit grid {nearest}/14={nearest/14:.4f} dist={float(mid)-nearest/14:+.4f}")
    if has_apex:
        print("  CRISS-CROSS at unit clocks (tight antipodal pairs v*j≡+1 needs eps>0 / ≡-1 needs eps<0):")
        for j,up,dn in cross_structure(V):
            print(f"     clock {j}/14: tight-up(eps>0 safe)={up}  tight-down(eps<0 safe)={dn}  {'<-- CROSS (both nonempty: forced conflict)' if up and dn else ''}")
