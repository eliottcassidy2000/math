# opus-2026-07-17-S377 -- WHY 12, AND ONLY 12?
#
# {1,...,13} with speed i replaced by r is tight iff the union still covers, i.e.
#     E_i  SUBSET OF  D_r        where   E_i := [0,1] \ union_{j != i} D_j
# is the ESSENTIAL REGION of speed i -- the part of the circle that ONLY speed i
# covers.  A small essential region is easy to re-cover, so the swap i -> r
# should be possible exactly when E_i is small and sits inside D_r.
from fractions import Fraction as F
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def complement(u):
    out=[]; prev=F(0)
    for a,b in u:
        if a>prev: out.append((prev,a))
        prev=b
    if prev<1: out.append((prev,F(1)))
    return out
def essential(V,i):
    allv=[]
    for x in V:
        if x!=i: allv.extend(teeth01(x))
    return complement(union(allv))
def contained(E,D):
    Du=union(D)
    for a,b in E:
        ok=False
        for c,d in Du:
            if c<=a and b<=d: ok=True; break
        if not ok: return False
    return True
def meas(iv): return sum(b-a for a,b in iv)

base=list(range(1,14))
print("(10) THE ESSENTIAL REGION E_i of each speed in {1,...,13}")
print("     i   |E_i|        #components   E_i inside D_2i ?")
rows=[]
for i in base:
    E=essential(base,i)
    ok=contained(E,teeth01(2*i))
    rows.append((i,meas(E),len(E),ok))
    print(f"     {i:2d}  {float(meas(E)):.6f}       {len(E):3d}          {'YES  <-- swap works' if ok else 'no'}")

print()
print("(11) CROSS-CHECK: the criterion should predict the swap results exactly")
print("     (S377 found i->2i tight ONLY for i=12)")
pred=[i for i,m,c,ok in rows if ok]
print(f"     criterion predicts i->2i works for: {pred}")
print(f"     search found:                       [12]")
print(f"     MATCH: {pred==[12]}")

print()
print("(12) WHY 12 -- the essential regions ranked by size")
for i,m,c,ok in sorted(rows,key=lambda r:r[1]):
    print(f"     i={i:2d}: |E_i| = {float(m):.6f}  ({c} components){'   <-- the swappable one' if ok else ''}")
print()
print("     A speed is swappable when its essential region is small enough AND")
print("     positioned inside the finer comb D_2i.  Both conditions are needed:")
print("     the ranking alone does not determine swappability.")
