from fractions import Fraction as Fr
def unsafe_arcs(S, thr=Fr(1,14)):
    arcs=[]
    for s in S:
        hw=thr/s  # half-width of U_s arc
        for k in range(s):
            c=Fr(k,s); lo=c-hw; hi=c+hw
            lo%=1
            if lo+2*hw<=1: arcs.append((lo,lo+2*hw))
            else: arcs.append((lo,Fr(1))); arcs.append((Fr(0),lo+2*hw-1))
    return arcs
def union(arcs):
    arcs=sorted(arcs); m=[]; clo=chi=None
    for lo,hi in arcs:
        if chi is None: clo,chi=lo,hi
        elif lo<=chi: chi=max(chi,hi)
        else: m.append((clo,chi)); clo,chi=lo,hi
    if chi is not None: m.append((clo,chi))
    return m
def complement(m):
    safe=[]; prev=Fr(0)
    for lo,hi in m:
        if lo>prev: safe.append((prev,lo))
        prev=max(prev,hi)
    if prev<1: safe.append((prev,Fr(1)))
    return safe
def safe_set(S): return complement(union(unsafe_arcs(S)))
def meas(iv): return sum(b-a for a,b in iv)
def intersect(A,B):
    out=[]; i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: out.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return out
seed=list(range(1,12))+[13]
Ss=safe_set(seed)
print(f"seed {{1..11,13}}: meas(Safe@1/14)={float(meas(Ss)):.5f} (S43 said 0.0122 -- should match now), #int={len(Ss)}")
if Ss: print(f"  smallest safe interval = {float(min(b-a for a,b in Ss)):.2e}")
for v in [30030, 60060, 510510]:
    both=intersect(Ss, safe_set([v]))
    print(f"  v={v}: meas(Safe(seed) ∩ Safe(v))={float(meas(both)):.6f}; M(seed∪{{v}})>=1/14: {meas(both)>0}; ~equidist pred {float(meas(Ss)*Fr(6,7)):.6f}")
