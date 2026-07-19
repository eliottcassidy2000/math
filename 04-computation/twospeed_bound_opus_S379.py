# opus-2026-07-17-S379 -- HYP-7630: MAKING THE TWO-SPEED CASE EXHAUSTIVE.
#
# THM-1125 made SINGLE substitutions exhaustive via  r <= 2*lam/ell_max(E_i).
# Two-speed is harder: a component of the joint essential region
#     E_{i,j} = I \ union_{k != i,j} D_k
# may be covered by arcs of BOTH new speeds, so "fits inside one arc" fails.
#
# THE ROUTE, in three steps.
# (1) DENSITY BOUND.  The window lemma (FragmentationLemma) gives
#         mu(D_w n I) <= 2*lam*|I| + 2*lam/w
#     for any interval I.  Full covering of a component of length ell needs
#         ell <= (2*lam*ell + 2*lam/r) + (2*lam*ell + 2*lam/s),
#     i.e.  ell*(1 - 4*lam) <= 2*lam*(1/r + 1/s).  At lam = 1/14 that is
#         ell <= (1/5)*(1/r + 1/s).
#     With r <= s this gives  r <= 2/(5*ell_max(E_{i,j}))  -- the SMALLER speed
#     is bounded.  (The larger is NOT bounded by this alone.)
# (2) LRC(13) MAKES THE RECURSION TERMINATE.  After fixing r, the residue
#     E' = E_{i,j} \ D_r must satisfy E' subset D_s.  E' is NONEMPTY: otherwise
#     the 12 speeds ({1..13} minus {i,j}) plus r would cover everything at
#     radius 1/14, contradicting LRC(13) -- 12 speeds always leave a point of
#     gap >= 1/13 > 1/14.  (LRC(<=13) is a citation hypothesis for this project.)
# (3) SINGLE-SPEED BOUND ON THE SECOND.  E' nonempty => s <= 2*lam/ell_max(E').
# So both speeds are bounded and the search is FINITE.
from fractions import Fraction as F
from itertools import combinations
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
def circ_lmax(comp):
    """longest component, treating [0,1] as a circle (merge wrap-around piece)."""
    if not comp: return F(0)
    c=list(comp)
    if len(c)>1 and c[0][0]==0 and c[-1][1]==1:
        wrap=(c[0][1]-c[0][0])+(c[-1][1]-c[-1][0])
        inner=max((b-a) for a,b in c[1:-1]) if len(c)>2 else F(0)
        return max(wrap,inner)
    return max(b-a for a,b in c)
def ess(V, drop):
    allv=[]
    for x in V:
        if x not in drop: allv.extend(teeth01(x))
    return complement(union(allv))
def diff_arcs(E, r):
    """E minus badArcs r, as intervals."""
    D=union(teeth01(r)); out=[]
    for a,b in E:
        cur=[(a,b)]
        for c,d in D:
            nxt=[]
            for x,y in cur:
                if d<=x or c>=y: nxt.append((x,y)); continue
                if x<c: nxt.append((x,min(c,y)))
                if y>d: nxt.append((max(d,x),y))
            cur=nxt
        out.extend([(x,y) for x,y in cur if y>x])
    return out
def contained(E,D):
    Du=union(D)
    return all(any(c<=a and b<=d for c,d in Du) for a,b in E)
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))

base=list(range(1,14))
print("(1) STEP 1 -- bound on the SMALLER new speed:  r <= 2/(5*ell_max(E_ij))")
print("    over all 78 pairs (i,j) of {1,...,13}")
rows=[]
for i,j in combinations(base,2):
    E=ess(base,{i,j}); L=circ_lmax(E)
    rb=int(F(2,5)/L) if L>0 else None
    rows.append((i,j,L,rb))
rbs=[rb for *_,rb in rows if rb]
print(f"    ell_max ranges {float(min(r[2] for r in rows)):.6f} .. {float(max(r[2] for r in rows)):.6f}")
print(f"    bound on smaller speed ranges {min(rbs)} .. {max(rbs)}   (MAX = {max(rbs)})")
worst=max(rows,key=lambda t:t[3] or 0)
print(f"    worst pair: (i,j)=({worst[0]},{worst[1]})  ell_max={float(worst[2]):.6f}  r <= {worst[3]}")

print()
print("(2) STEP 2 -- verify E' = E_ij \ D_r is NEVER empty (LRC(13) predicts this)")
empt=0; checked=0
for i,j,L,rb in rows:
    E=ess(base,{i,j})
    for r in range(1,(rb or 0)+1):
        Ep=diff_arcs(E,r); checked+=1
        if not Ep: empt+=1; print(f"    EMPTY at (i,j,r)=({i},{j},{r}) -- would contradict LRC(13)!")
print(f"    checked {checked} (pair, r) combinations; empty residues: {empt}")
print("    (0 empties confirms the LRC(13) argument and the recursion terminates)")
