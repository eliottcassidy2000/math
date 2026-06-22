"""
HONEST scope of Bonferroni-3 (kps-S31u): does it hold on the BINDING configs (large base + TIGHT far,
p0 near cap = the actual dangerous case, HYP-2797 doublet) even though it fails on spread-far slack?
"""
from itertools import combinations
def p0(E):
    E=sorted(set(e for e in E if e!=0))
    if not E: return 0.0
    bset={0.0,1.0}
    for e in E:
        for j in range(8):
            b=j/7.0; m=0
            while True:
                xv=(b+m)/e
                if xv>=1: break
                if xv>=0: bset.add(xv)
                m+=1
    Bs=sorted(bset); tot=0.0
    for lo,hi in zip(Bs,Bs[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2
        if len(set(int((e*mid)%1*7) for e in E)&set(range(1,7)))==6: tot+=hi-lo
    return tot
def s3(base,far):
    pf={S:p0(base+[far[i] for i in S]) for r in range(min(4,len(far)+1)) for S in combinations(range(len(far)),r)}
    tot=0.0
    for r in range(1,min(4,len(far)+1)):
        for S in combinations(range(len(far)),r):
            for r2 in range(r+1):
                for T in combinations(S,r2):
                    if T in pf: tot+=((-1)**(r-r2))*pf[T]
    return tot
cap=0.3815
print("BINDING configs (large consec base + TIGHT far cluster -- the near-cap/doublet case):")
binding=[([0,1,2,3,4,5],[f,f+1]) for f in [12,15,20,30]] + [([0,1,2,3,4,5],[f,f+1,f+2]) for f in [12,20,30]]
for base,far in binding:
    full=p0(base+far); b3=s3(base,far)
    print(f"  base={base} far={far}: p0={full:.4f} (cap={cap}) Bonf3 bound={b3:.4f}  p0<=Bonf3? {full<=b3+1e-9}")
print("\nSLACK configs (small base + SPREAD far -- p0<<cap):")
slack=[([0,3,9],[27,28,42,43,47]),([0,3,7],[50,53,54,57,64]),([0,2,5],[33,41,52,58,61])]
for base,far in slack:
    full=p0(base+far); b3=s3(base,far)
    print(f"  base={base} far={far}: p0={full:.4f} (cap={cap}) Bonf3 bound={b3:.4f}  p0<=Bonf3? {full<=b3+1e-9}")
print("\n=> Bonferroni-3 holds on BINDING (tight-far, near-cap) configs -- the dangerous ones; fails only")
print("   on SLACK (spread-far) where p0<<cap. So it IS a valid handle on the binding doublet leg, but")
print("   is NOT a standalone universal closure (the slack leg needs the separate p0<<cap argument).")
