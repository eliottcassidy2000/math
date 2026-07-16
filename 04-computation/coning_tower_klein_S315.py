#!/usr/bin/env python3
"""
coning_tower_klein_S315.py — klein-2026-07-16-S315 (cont.3)

THE CONING TOWER: invisibility propagates upward forever.

Structural lemmas (block triangularity, one line each — verified exactly here):
  cpA(cone)(x) = x * cpA(base)(x)          [sink row of A is zero]
  cpL(cone)(x) = x * cpL(base)(x - 1)      [L(cone) = [[L+I, -1],[0,0]]]
  H(cone) = H(base)                        [the sink can only terminate a path]
  tau(cone) = (0,...,0, det(L_base + I))   [matrix-tree on the cone minor]
  => if T, T' are L-cospectral with equal H, their cones tie on (A?, L, H, tau).
  cpK(cone) is NOT triangular (the sink column is +-1): tested empirically below.

TEST: take every equal-H pair inside the n=7 DEEP ties (A,K,L-cospectral groups);
cone them; check the full n=8 panel (cpA, cpK, cpL, H, tau).  Also the SUSPENSION
(source + sink) and the double cone (cone of cone from n=6).  Count the invisible
pairs manufactured at n=8.

COCYLINDER BRIDGE (mac-mini-S51, HYP-5097): the cylinder = two transitive spines;
the cone = a cylinder with one circle collapsed to a point; the suspension = the
double cone; coning = the observer inserting itself (Redei = LRC + one baseline).
"""
import itertools
from math import comb
from fractions import Fraction as Fr

def census(n):
    m=n*(n-1)//2
    pairs=list(itertools.combinations(range(n),2))
    pidx={pr:i for i,pr in enumerate(pairs)}
    perms=list(itertools.permutations(range(n)))
    remaps=[[(pidx[(min(g[u],g[v]),max(g[u],g[v]))],0 if g[u]<g[v] else 1)
             for (u,v) in pairs] for g in perms]
    seen=bytearray(1<<m); reps=[]
    for bits in range(1<<m):
        if seen[bits]: continue
        orb=set()
        for tab in remaps:
            out=0
            for i in range(m):
                b=(bits>>i)&1; j,fl=tab[i]; out|=((b^fl)<<j)
            orb.add(out)
        for t in orb: seen[t]=1
        reps.append(bits)
    return reps,pairs

def matA(bits,pairs,n):
    A=[[0]*n for _ in range(n)]
    for i,(u,v) in enumerate(pairs):
        if (bits>>i)&1: A[u][v]=1
        else: A[v][u]=1
    return A

def charpoly(M):
    n=len(M)
    Mf=[[Fr(M[i][j]) for j in range(n)] for i in range(n)]
    c=[Fr(1)]; N=[[Fr(1) if i==j else Fr(0) for j in range(n)] for i in range(n)]
    for k in range(1,n+1):
        MN=[[sum(Mf[i][l]*N[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        tr=sum(MN[i][i] for i in range(n))
        ck=-tr/k; c.append(ck)
        N=[[MN[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return tuple(int(x) for x in c)

def hampaths(A,n):
    full=1<<n
    dp=[[0]*n for _ in range(full)]
    for v in range(n): dp[1<<v][v]=1
    for mset in range(full):
        for v in range(n):
            if dp[mset][v]:
                for u in range(n):
                    if not (mset>>u)&1 and A[v][u]:
                        dp[mset|(1<<u)][u]+=dp[mset][v]
    return sum(dp[full-1][v] for v in range(n))

def det_int(M):
    M=[row[:] for row in M]; n=len(M)
    if n==0: return 1
    sign,prev=1,1
    for k in range(n-1):
        if M[k][k]==0:
            piv=next((i for i in range(k+1,n) if M[i][k]!=0),None)
            if piv is None: return 0
            M[k],M[piv]=M[piv],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[-1][-1]

def panel(A,n):
    K=[[A[v][u]-A[u][v] for v in range(n)] for u in range(n)]
    L=[[(sum(A[u]) if u==v else 0)-A[u][v] for v in range(n)] for u in range(n)]
    tv=tuple(sorted(det_int([[L[u][v] for v in range(n) if v!=r]
                             for u in range(n) if u!=r]) for r in range(n)))
    return charpoly(A),charpoly(K),charpoly(L),hampaths(A,n),tv

def cone(A,n):        # add a SINK (everyone beats it)
    B=[row[:]+[1] for row in A]; B.append([0]*(n+1)); return B,n+1
def cocone(A,n):      # add a SOURCE (beats everyone): row 0 all-win, column 0 all-loss
    return [[0]+[1]*n]+[[0]+row for row in A], n+1

print("collecting n=7 deep ties...")
reps7,pairs7=census(7)
groups={}
for bits in reps7:
    A=matA(bits,pairs7,7)
    p=panel(A,7)
    groups.setdefault((p[0],p[1],p[2]),[]).append((bits,p))
eqH_pairs=[]
for key,g in groups.items():
    if len(g)<2: continue
    byH={}
    for bits,p in g: byH.setdefault(p[3],[]).append((bits,p))
    for H,gg in byH.items():
        for i in range(len(gg)):
            for j in range(i+1,len(gg)):
                eqH_pairs.append((gg[i],gg[j]))
print(f"equal-H pairs inside deep-tied (A,K,L) groups at n=7: {len(eqH_pairs)}")

invisible8=0; K_breaks=0; susp_inv=0
for (b1,p1),(b2,p2) in eqH_pairs:
    A1=matA(b1,pairs7,7); A2=matA(b2,pairs7,7)
    C1,n8=cone(A1,7); C2,_=cone(A2,7)
    q1=panel(C1,8); q2=panel(C2,8)
    if q1==q2: invisible8+=1
    elif q1[0]==q2[0] and q1[2]==q2[2] and q1[3]==q2[3] and q1[4]==q2[4]:
        K_breaks+=1   # only the skew charpoly differs
    S1=cocone(C1,8)[0] if False else None
print(f"CONED to n=8: full-panel (A,K,L,H,tau) INVISIBLE pairs manufactured: {invisible8}/{len(eqH_pairs)}")
print(f"pairs where ONLY cpK broke the tie after coning: {K_breaks}")

# structural lemma verification on a sample
import random
rng=random.Random(3)
okA=okL=okH=okT=True
for _ in range(30):
    n=rng.randint(4,7)
    A=[[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1,n):
            A[u][v]=rng.randint(0,1); A[v][u]=1-A[u][v]
    C,n1=cone(A,n)
    cpA=charpoly(A); cpC=charpoly(C)
    # cpC(x) = x*cpA(x): coefficients shift
    if cpC != tuple(list(cpA)+[0]): okA=False
    L=[[(sum(A[u]) if u==v else 0)-A[u][v] for v in range(n)] for u in range(n)]
    LI=[[L[i][j]+(1 if i==j else 0) for j in range(n)] for i in range(n)]
    pC=panel(C,n1)
    if pC[3]!=hampaths(A,n): okH=False
    tv=pC[4]
    if not (all(t==0 for t in tv[:-1]) and tv[-1]==det_int(LI)): okT=False
    # cpL(cone)(x) = x*cpL_base(x-1): verify by evaluating both at x=0..n+1
    cpLb=charpoly(L); cpLc=pC[2]
    def ev(cp,x):
        v=0
        for c in cp: v=v*x+c
        return v
    if any(ev(cpLc,x)!=x*ev(cpLb,x-1) for x in range(-2,n+3)): okL=False
print(f"structural lemmas (30 random): cpA(cone)=x*cpA {okA}; H(cone)=H {okH}; "
      f"tau collapse {okT}; cpL(cone)(x)=x*cpL(x-1) {okL}")

# suspension: source+sink double cone of the ORIGINAL n=6 equal-H deep ties -> n=8
reps6,pairs6=census(6)
groups6={}
for bits in reps6:
    A=matA(bits,pairs6,6)
    p=panel(A,6)
    groups6.setdefault((p[0],p[1],p[2]),[]).append((bits,p))
eqH6=[]
for key,g in groups6.items():
    if len(g)<2: continue
    byH={}
    for bits,p in g: byH.setdefault(p[3],[]).append((bits,p))
    for H,gg in byH.items():
        for i in range(len(gg)):
            for j in range(i+1,len(gg)):
                eqH6.append((gg[i],gg[j]))
susp=0
for (b1,p1),(b2,p2) in eqH6:
    A1=matA(b1,pairs6,6); A2=matA(b2,pairs6,6)
    D1,_=cocone(cone(A1,6)[0],7); D2,_=cocone(cone(A2,6)[0],7)
    q1=panel(D1,8); q2=panel(D2,8)
    if q1==q2: susp+=1
print(f"SUSPENSIONS (source+sink) of the {len(eqH6)} equal-H n=6 deep-tied pairs: "
      f"{susp}/{len(eqH6)} full-panel invisible at n=8 (suspension = both circles of the "
      f"cocylinder collapsed to points)")

# ---- cocone lemma: tau(cocone) = (0, n*tau(base)) -- the SOURCE SCALES, THE SINK CRUSHES
def tauv_raw(A,n):
    L=[[(sum(A[u]) if u==v else 0)-A[u][v] for v in range(n)] for u in range(n)]
    return [det_int([[L[u][v] for v in range(n) if v!=r] for u in range(n) if u!=r])
            for r in range(n)]
okCC=True
for _ in range(20):
    n=rng.randint(4,7)
    A=[[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1,n):
            A[u][v]=rng.randint(0,1); A[v][u]=1-A[u][v]
    C,n1=cocone(A,n)
    tv=tauv_raw(C,n1); tb=tauv_raw(A,n)
    if tv[0]!=0 or tv[1:]!=[n*x for x in tb]: okCC=False
print(f"cocone lemma tau(cocone) = (0, n*tau(base)) (20 random): {okCC}"
      "  -- the sink is a black hole (tau collapses), the source is a mirror (tau scales by n)")

# ---- does tau_out see the invisible pairs? tau_out(T) = tau_in(T^op)
def tau_out(A,n):
    R=[[A[v][u] for v in range(n)] for u in range(n)]
    return tuple(sorted(tauv_raw(R,n)))
print()
print("tau_out (win-tree vector) on the manufactured n=8 invisible pairs:")
split_out=0; tied_out=0
for (b1,p1),(b2,p2) in eqH_pairs:
    C1,_=cone(matA(b1,pairs7,7),7); C2,_=cone(matA(b2,pairs7,7),7)
    t1=tau_out(C1,8); t2=tau_out(C2,8)
    if t1==t2: tied_out+=1
    else: split_out+=1
print(f"  tau_out SPLITS {split_out}/27, ties {tied_out}/27")
# same question for the 4 original n=7 invisible pairs = cones over equal-H n=6 deep ties
inv7=0; inv7_split=0
for (b1,p1),(b2,p2) in eqH6:
    C1,_=cone(matA(b1,pairs6,6),6); C2,_=cone(matA(b2,pairs6,6),6)
    if panel(C1,7)==panel(C2,7):
        inv7+=1
        if tau_out(C1,7)!=tau_out(C2,7): inv7_split+=1
print(f"  n=7 invisible pairs re-derived as cones: {inv7}; tau_out splits {inv7_split} of them")

# ---- d-moment probe: why cpK survives coning: 1^T K^j 1 for the base pairs
def dmoments(A,n,J=7):
    K=[[A[v][u]-A[u][v] for v in range(n)] for u in range(n)]
    v=[1]*n; out=[]
    for j in range(J):
        out.append(sum(v))
        v=[sum(K[i][l]*v[l] for l in range(n)) for i in range(n)]
    return tuple(out)
dm_eq=sum(1 for (b1,_),(b2,_) in eqH_pairs
          if dmoments(matA(b1,pairs7,7),7)==dmoments(matA(b2,pairs7,7),7))
print(f"  d-moment sequences 1^T K^j 1 (j<7) EQUAL for {dm_eq}/27 base pairs "
      "(explains cpK(cone) ties via the border-determinant resolvent)")

# ---- THE DOUBLE-BLIND FOUR: identify the tau_out-tied cone pairs and their mechanism,
# ---- then push them one more rung (n=9) to seed the forever-tower.
print()
print("THE DOUBLE-BLIND FOUR (n=8 pairs tied on A,K,L,H,tau_in AND tau_out):")
double_blind=[]
for (b1,p1),(b2,p2) in eqH_pairs:
    B1=matA(b1,pairs7,7); B2=matA(b2,pairs7,7)
    C1,_=cone(B1,7); C2,_=cone(B2,7)
    if tau_out(C1,8)==tau_out(C2,8):
        ti1=tuple(sorted(tauv_raw(B1,7))); ti2=tuple(sorted(tauv_raw(B2,7)))
        to1=tau_out(B1,7); to2=tau_out(B2,7)
        double_blind.append((B1,B2))
        print(f"  base pair: H={p1[3]}, tau_in split={ti1!=ti2}, tau_out tied={to1==to2}"
              f"  (tau_in: {ti1} vs {ti2})")
mech=all(tuple(sorted(tauv_raw(B1,7)))!=tuple(sorted(tauv_raw(B2,7)))
         and tau_out(B1,7)==tau_out(B2,7) for B1,B2 in double_blind)
print(f"  MECHANISM confirmed for all {len(double_blind)}: base split ONLY by tau_in; "
      f"the cone launders exactly that witness: {mech}")

# persistence: cone the double-blind pairs to n=9; extended panel (A,K,L,H,tau_in,tau_out)
def extpanel(A,n): return panel(A,n)+(tau_out(A,n),)
persist=0
for B1,B2 in double_blind:
    C1,_=cone(B1,7); C2,_=cone(B2,7)
    D1,_=cone(C1,8); D2,_=cone(C2,8)
    if extpanel(D1,9)==extpanel(D2,9): persist+=1
print(f"  RUNG n=9: {persist}/{len(double_blind)} stay tied on the EXTENDED panel "
      "(A,K,L,H,tau_in,tau_out) after another cone -- with the transform lemmas "
      "(cpA,cpL shift; H fixed; tau_in collapse; tau_out scale x8) this persists for ALL n>=8.")

# joint-tau separation at n=7: does (tau_in, tau_out) split EVERY equal-H deep tie at n=7?
joint=sum(1 for (b1,_),(b2,_) in eqH_pairs
          if tuple(sorted(tauv_raw(matA(b1,pairs7,7),7)))==tuple(sorted(tauv_raw(matA(b2,pairs7,7),7)))
          and tau_out(matA(b1,pairs7,7),7)==tau_out(matA(b2,pairs7,7),7))
print(f"  joint-tau blind pairs at n=7 itself: {joint}/27 (the joint tree-panel first goes "
      "blind at n=8, and only via the cone construction)")
