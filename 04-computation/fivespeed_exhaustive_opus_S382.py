# opus-2026-07-17-S382 -- HYP-7670: THE FIVE-SPEED CASE, EXHAUSTIVE.
#
# k=5 needs a second prune.  THM-1165's per-level bound uses only ell_max:
#     each component must satisfy  ell <= j/(r*(7-j)).
# THERE IS A GENUINELY ADDITIONAL GLOBAL CONSTRAINT.  If E is a union of c
# intervals, then for one comb
#     mu(D_w n E) = sum_i mu(D_w n I_i) <= sum_i (2*lam*|I_i| + 2*lam/w)
#                 = 2*lam*mu(E) + 2*lam*c/w,
# so covering E by j combs all >= r forces
#     mu(E)*(1 - 2*j*lam) <= 2*lam*c*sum(1/w_i) <= 2*lam*c*j/r,
# i.e.  mu(E) <= c(E) * j / (r*(7-j))   at lam = 1/14.
# The two prunes are complementary: the ell_max one binds when a single long
# component dominates, the measure one when many components share the budget.
# Both are sound (each is implied by full covering), so applying both only
# removes nodes that provably cannot lead to a tight family.
from fractions import Fraction as F
from itertools import combinations
import sys, time
LAM=F(1,14)
TEETH={}; UT={}
def T(x):
    if x not in TEETH:
        w=LAM/x; out=[]
        for j in range(0,x+1):
            a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
            if a<b: out.append((a,b))
        TEETH[x]=out
    return TEETH[x]
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def UT_(x):
    if x not in UT: UT[x]=union(T(x))
    return UT[x]
def complement(u):
    out=[]; prev=F(0)
    for a,b in u:
        if a>prev: out.append((prev,a))
        prev=b
    if prev<1: out.append((prev,F(1)))
    return out
def subtract(E,D):
    out=[]; j=0; n=len(D)
    for a,b in E:
        cur=a
        while j<n and D[j][1]<=cur: j+=1
        jj=j
        while jj<n and D[jj][0]<b:
            c,d=D[jj]
            if c>cur: out.append((cur,min(c,b)))
            if d>cur: cur=d
            if cur>=b: break
            jj+=1
        if cur<b: out.append((cur,b))
    return out
def circ_stats(c):
    """(longest component treating [0,1] circularly, total measure, count)."""
    if not c: return F(0),F(0),0
    tot=sum(b-a for a,b in c)
    if len(c)>1 and c[0][0]==0 and c[-1][1]==1:
        wrap=(c[0][1]-c[0][0])+(c[-1][1]-c[-1][0])
        inner=max((b-a) for a,b in c[1:-1]) if len(c)>2 else F(0)
        return max(wrap,inner), tot, len(c)-1
    return max(b-a for a,b in c), tot, len(c)
def ess(V,drop):
    allv=[]
    for x in V:
        if x not in drop: allv.extend(T(x))
    return complement(union(allv))
def coef(j, r):
    """per-component budget with j combs all >= r; sharper single-arc at j=1."""
    if j==1: return 2*LAM/r
    return F(j, r*(7-j))

base=list(range(1,14))
st={'nodes':0,'p_max':0,'p_meas':0,'empty':0,'leaf':0}
tight=set(); t0=time.time()

def rec(E, j, rmin, chosen, dropped):
    st['nodes']+=1
    L, mu, c = circ_stats(E)
    B = coef(j, rmin)
    if L > B: st['p_max']+=1; return            # longest component uncoverable
    if j>1 and mu > c*F(j, rmin*(7-j)):         # global measure budget exceeded
        st['p_meas']+=1; return
    wb = int((2*LAM)/L) if j==1 else int(F(j,(7-j))/L)
    for w in range(rmin, wb+1):
        E2=subtract(E, UT_(w))
        if j==1:
            st['leaf']+=1
            if not E2:
                V=tuple(sorted(set([x for x in base if x not in dropped]+chosen+[w])))
                if len(V)==13: tight.add(V)
        else:
            if not E2: st['empty']+=1; continue
            rec(E2, j-1, w, chosen+[w], dropped)

quints=list(combinations(base,5))
for idx,q in enumerate(quints):
    E=ess(base,set(q))
    if E: rec(E, 5, 1, [], set(q))
    if idx%150==0:
        print(f"  [{idx}/{len(quints)}] nodes={st['nodes']} pruned={st['p_max']}+{st['p_meas']} "
              f"leaf={st['leaf']} t={time.time()-t0:.0f}s", flush=True)
print(f"nodes visited: {st['nodes']}")
print(f"pruned by ell_max: {st['p_max']}   pruned by measure: {st['p_meas']}")
print(f"empty residues (would contradict LRC(10)..LRC(13)): {st['empty']}")
print(f"leaf checks: {st['leaf']}   elapsed {time.time()-t0:.0f}s")
print(f"DISTINCT TIGHT FAMILIES: {len(tight)}")
B0=tuple(base); T2=(1,2,3,4,5,6,7,8,9,10,11,13,24)
for V in sorted(tight):
    tag = "= {1..13}" if V==B0 else ("= {1..11,13,24}" if V==T2 else "*** NEW ***")
    print(f"    {list(V)}   [{tag}]")
