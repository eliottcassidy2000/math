# opus-2026-07-17-S369 -- THE PAYOFF: does the 13-term AP have a limiting
# uncovered measure?  The lattice argument predicts YES -- for an AP the
# resonance lattice's s=0 part is independent of (a,d), and the s != 0 part
# needs |sum n_i| >= d, so it dies as d grows.  If the exact uncovered measure
# converges to a POSITIVE constant, the additively-richest (extremal) corner
# is settled for all large d at once.
from fractions import Fraction as F
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    lo=0
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0], max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    u=union(allv)
    return 1-sum(b-a for a,b in u)

print("13-TERM APs {a, a+d, ..., a+12d}: exact uncovered measure vs d")
print("    a    d    max speed   uncovered (exact)     float")
prev=None
for (a,d) in [(1,1),(1,2),(1,3),(1,5),(1,8),(1,13),(1,21),(1,34),(1,55),(1,89),
              (2,3),(3,5),(5,8),(2,21),(3,34),(7,55)]:
    V=[a+i*d for i in range(13)]
    u=uncovered(V)
    print(f"    {a:3d}  {d:3d}   {max(V):8d}   {str(u):20s} {float(u):.8f}")
print()
print("If these converge to a common positive value as d grows, every")
print("sufficiently-spread AP is lonely, and the extremal corner closes.")
