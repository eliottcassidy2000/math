from fractions import Fraction as F
from math import gcd
def nd_res(v,a,q): return (v*a)%q
def val_q_a(V):
    # exact M via THM-999: M=val/q, q | v_i+v_j for some pair, q<=2max; maximizer t=a/q.
    # We scan candidate q from all pairwise sums' divisors up to 2*max, and a in 1..q-1,
    # M(t=a/q)=min_i |v_i a mod q|_q / q ; take the global max over (q,a). Return best.
    mx=max(V); best=None
    Qs=set()
    for i in range(len(V)):
        for j in range(len(V)):
            s=V[i]+V[j]
            for d in range(1,2*mx+1):
                if s% d==0 and d<=2*mx: Qs.add(d)
    Qs=[q for q in Qs if q>=2]
    for q in Qs:
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q, q-((v*a)%q)) for v in V)  # val = min circular residue
            if best is None or F(m,q)>F(best[0],best[1]):
                best=(m,q,a)
    val,q,a=best
    return val,q,a
def gaps_struct(V,val,q,a):
    R=sorted(((v*a)%q) for v in V)
    # residues should lie in [val,q-val]; reflect if needed so min touches val side
    gaps=[R[i+1]-R[i] for i in range(len(R)-1)]
    nsmall=sum(1 for g in gaps if g<val)
    return R,gaps,nsmall
# deep wells
for m in [1,2,3]:
    V=list(range(1,13))+[182*m]
    val,q,a=val_q_a(V)
    R,gaps,nsmall=gaps_struct(V,val,q,a)
    isAP = sorted(V[:-1])==[i for i in range(1,13)] # trivially AP core here
    print("deepwell m=%d V_max=%d: M=%s(=%.5f<1/13?%s) val=%d q=%d a=%d | #gaps<val=%d | span=%d=11val+%d? g0=q-13val=%d"
          %(m,V[-1],F(val,q),val/q, val/q<1/13, val,q,a, nsmall, R[-1]-R[0], (R[-1]-R[0])-11*val, q-13*val))
    print("   sorted gaps:",gaps)
