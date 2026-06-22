"""
Applying the 3-set Venn understanding to the WIDE BOUND (kps-S31t).
The LRC multi-far cover is the 3-set Venn (THM-548): p0(B u {u,v,w}) = p0(B) + (1-far corners)
- (2-far edges) + (3-far center), MOBIUS signs. The Newton packets:
  D_u = p0(B+u)-p0(B)                                   (corner, 1-far)
  D_uv = p0(B+u+v)-p0(B+u)-p0(B+v)+p0(B)                (edge, 2-far, mixed 2nd diff)
  D_uvw = p0(B+uvw)-sum(pairs)+sum(singles)-p0(B)       (center, 3-far)
TEST: (1) inclusion-exclusion holds; (2) are the 2-far EDGES negative (=> BONFERRONI:
p0(full) <= p0(B)+sum(1-far), reducing wide bound to single-far THM-563); (3) is 3-far center
SUB-DOMINANT (the even/degenerate prediction)?
"""
from itertools import combinations
def p0(E):
    E=sorted(set(e for e in E if e!=0)); 
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
    B=sorted(bset); tot=0.0
    for lo,hi in zip(B,B[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2
        if len(set(int((e*mid)%1*7) for e in E)&set(range(1,7)))==6: tot+=hi-lo
    return tot
def packets(base, far):
    # far = [u,v,w]; return dict of Newton mixed differences over subsets of far
    pf={}
    for r in range(len(far)+1):
        for S in combinations(range(len(far)), r):
            E=base+[far[i] for i in S]
            pf[S]=p0(E)
    D={}
    for r in range(1,len(far)+1):
        for S in combinations(range(len(far)), r):
            # mixed r-th difference
            val=0.0
            for r2 in range(r+1):
                for T in combinations(S, r2):
                    val += ((-1)**(r-r2))*pf[T]
            D[S]=val
    return pf, D
print("3-FAR VENN packets for wide configs (base + 3 far runners):  cap_8~0.382")
tests = [([0,1,2,3,4], [20,21,22]),    # consec base + tight far triplet
         ([0,1,2,3,4], [20,40,60]),    # spread far
         ([0,2,4,6,8], [25,26,27]),    # even base + far
         ([0,1,2,3,4,5], [30,31,32])]
for base,far in tests:
    pf,D=packets(base,far)
    one=sum(D[S] for S in D if len(S)==1)
    two=sum(D[S] for S in D if len(S)==2)
    three=sum(D[S] for S in D if len(S)==3)
    full=pf[(0,1,2)]; b0=pf[()]
    print(f"\n base={base} far={far}:")
    print(f"   p0(B)={b0:.4f}  p0(full)={full:.4f}  inc-excl: B+1far-2far+3far={b0+one+two+three:.4f} {'OK' if abs(b0+one+two+three-full)<1e-6 else 'X'}")
    print(f"   1-far(corners) sum={one:+.4f}   2-far(edges) sum={two:+.4f}   3-far(center)={three:+.4f}")
    print(f"   BONFERRONI-1 upper: p0(B)+1far = {b0+one:.4f}  >= p0(full)={full:.4f}? {b0+one>=full-1e-9}")
    print(f"   2-far edges negative (sub-additive)? {two<0};  |3-far| << |2-far| (center sub-dom)? {abs(three)<abs(two) if two!=0 else 'na'}")
