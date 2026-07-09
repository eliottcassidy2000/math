"""kps-S116: extend the pair-sum native_decide leg to Vmax<=20 + the ratio tiling census.
Enumerates the 6084 primitive covering [1,20] 13-sets, finds a grid-free pair-sum band witness
(p,q) for each (q=v_i+v_j, all residues (v_l*p) mod q in [q/14,13q/14]); reports the ratio split
(ratio<=13 = covered non-enumeratively by mreach_ge_of_pairsum_ratioBand; ratio>13 = enumeration
residual = mac-mini's C1/C2/C3 territory). Emits LRCCoveringVmax20.lean (chunked, native_decide)."""
from itertools import combinations
from math import gcd
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def prim(S):
    g=0
    for s in S:
        g=gcd(g,s)
        if g==1: break
    return g==1
def psmod(v): return sorted({v[i]+v[j] for i in range(len(v)) for j in range(i,len(v))})
def band_witness(v):
    for q in psmod(v):
        if q<2: continue
        for p in range(1,q):
            if all(q<=14*((x*p)%q)<=13*q for x in v): return (p,q)
    return None
if __name__=="__main__":
    for VB in [18,20]:
        cov=[S for S in combinations(range(1,VB+1),13) if is_cov(S) and prim(S)]
        wit=[(S,band_witness(list(S))) for S in cov]
        miss=sum(1 for _,w in wit if w is None)
        r13=sum(1 for S,_ in wit if max(S)<=13*min(S))
        mq=max(w[1] for _,w in wit if w); mp=max(w[0] for _,w in wit if w)
        print(f"[1,{VB}]: covering={len(cov)} witnessed={len(cov)-miss} missing={miss} | "
              f"ratio<=13(ratio-band cert)={r13} ratio>13(residual)={len(cov)-r13} | max q={mq} p={mp}")
