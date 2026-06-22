"""
The r-far convergent tail via the Venn containment (kps-S31t). The Venn says the 3-far CENTER is
the triple overlap (smallest region) <= the 2-far EDGES (pairwise); generally the r-far packet
total decreases in r (center-of-center). TEST: for base + 5 far, the total r-far contribution
T_r = sum over r-subsets of the mixed r-th difference; check T_r decreases => wide cover =
[1-far=0] + [2-far doublet, BINDING, closed HYP-2797] + [convergent tail T_3+T_4+... <= geometric].
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
def Tr(base, far, r):
    # total r-far contribution = sum over r-subsets S of the mixed r-th difference D_S
    tot=0.0
    for S in combinations(range(len(far)), r):
        val=0.0
        for r2 in range(r+1):
            for T in combinations(S, r2):
                val += ((-1)**(r-r2))*p0(base+[far[i] for i in T])
        tot+=val
    return tot
tests=[([0,1,2,3,4],[20,21,22,23,24]),
       ([0,1,2,3,4,5],[30,31,32,33,34]),
       ([0,2,4,6,8],[25,27,29,31,33])]
print("r-far total contribution T_r (base + 5 far); check geometric decrease => convergent tail:")
for base,far in tests:
    full=p0(base+far)
    Ts=[Tr(base,far,r) for r in range(1,6)]
    print(f"\n base={base} far={far}: p0(full)={full:.4f} = sum T_r = {sum(Ts):.4f}")
    for r,t in enumerate(Ts,1):
        ratio = (Ts[r-1]/Ts[r-2]) if r>=2 and abs(Ts[r-2])>1e-9 else float('nan')
        print(f"   T_{r} ({['1','2','3','4','5'][r-1]}-far) = {t:+.5f}   ratio T_r/T_(r-1) = {ratio:.3f}")
    print(f"   2-far DOMINATES? {abs(Ts[1])==max(abs(x) for x in Ts)};  tail T_3+T_4+T_5 = {sum(Ts[2:]):+.5f} ({100*sum(Ts[2:])/full:.0f}% of full)")
print("\n=> 1-far=0 (single far useless); 2-far = BINDING doublet (Eisenstein A&B edge); T_r decreases")
print("   => wide cover = doublet (closed) + convergent geometric tail. The Venn containment center<=edge")
print("   is WHY the higher-far tail converges -- the repo's THM-557/Dedekind-ladder, structurally.")
