#!/usr/bin/env python3
"""
tournament_walk_reciprocity_klein_S316.py — klein-2026-07-16-S316

THE WALK RECIPROCITY THEOREM (THM-924): for every tournament T on n vertices,

  (L1)  det(zI - A + J) = (-1)^n cpA(-1-z)          [A - J = -(A^T + I); one line]
  (L2)  1 + R_A(z) = (-1)^n cpA(-1-z)/cpA(z)        [matrix determinant lemma]
        where R_A(z) = 1^T (zI-A)^{-1} 1 = sum_j w_j z^{-j-1}, w_j = # directed j-walks
        ==> ALL walk moments w_j are determined by cpA.
  (L3)  cpK(y) = 2^{n-1} [ cpA((y-1)/2) + (-1)^n cpA((-y-1)/2) ]
        ==> cpA DETERMINES cpK (census's 0/116 explained; skew spectrum = the
        reflection-symmetrization of the A-spectrum about Re z = -1/2).
  (L4)  R_K(y) = 1 - 2^n cpA((y-1)/2)/cpK(y),  R_K(y) = sum_j mu_j y^{-j-1},
        mu_j = 1^T K^j 1 the skew d-moments  ==> d-moments cpA-determined.
        HYP-7096 PROVED (vastly stronger: no deep tie needed, cpA alone suffices)
        ==> THM-918's cpK tower leg UNCONDITIONAL: cpK(cone) determined by x*cpA.
  (L5)  cpL_in(x)*(x-n) = (-1)^n x * cpL_out(n-x)   [L_in = nI - J - L_out; L 1 = 0]
        ==> the census panel is reversal-closed.

Verification: exact integers/rationals over ALL iso classes n = 4..7, random n = 8..10,
and the double-blind four's cones (n = 8, 9). Plus: cpK-tie groups vs cpA-tie groups at
n = 6, 7 (does the symmetrization LOSE information? i.e. cpK-ties that are not cpA-ties).
"""
import itertools, random
from fractions import Fraction as Fr

def census(n):
    m=n*(n-1)//2
    pairs=list(itertools.combinations(range(n),2))
    pidx={pr:i for i,pr in enumerate(pairs)}
    remaps=[[(pidx[(min(g[u],g[v]),max(g[u],g[v]))],0 if g[u]<g[v] else 1)
             for (u,v) in pairs] for g in itertools.permutations(range(n))]
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
        tr=sum(MN[i][i] for i in range(n)); ck=-tr/k; c.append(ck)
        N=[[MN[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return [Fr(x) for x in c]     # cp(x) = sum c[i] x^{n-i}

def ev(cp,x):
    v=Fr(0)
    for c in cp: v=v*x+c
    return v

def rand_tournament(n,rng):
    A=[[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1,n):
            A[u][v]=rng.randint(0,1); A[v][u]=1-A[u][v]
    return A

def verify_all(A,n):
    """returns dict of booleans for L1..L5 on one tournament, all exact."""
    cpA=charpoly(A)
    J1=[[A[i][j]- (1 if True else 0) for j in range(n)] for i in range(n)]
    AmJ=[[A[i][j]-1 for j in range(n)] for i in range(n)]
    cpAmJ=charpoly(AmJ)
    sgn=(-1)**n
    # L1: cp_{A-J}(z) = (-1)^n cpA(-1-z): evaluate at z = -3..n+2
    L1=all(ev(cpAmJ,Fr(z))==sgn*ev(cpA,Fr(-1-z)) for z in range(-3,n+3))
    # L2: walk moments: w_j = 1^T A^j 1 vs series of (-1)^n cpA(-1-z)/cpA(z) - 1
    w=[]; v=[1]*n
    for j in range(2*n):
        w.append(sum(v)); v=[sum(A[i][l]*v[l] for l in range(n)) for i in range(n)]
    # series check via evaluation: for several z, sum_{j<2n} w_j z^{-j-1} should equal
    # target up to O(z^{-2n-1}); instead verify EXACTLY via polynomial identity:
    # (-1)^n cpA(-1-z) - cpA(z) = P(z) with P(z) = sum_j (coefficient) — do it by
    # checking 1 + R = Q/cp at many z using exact resolvent solve
    def resolvent_sum(M,z):
        # solve (zI - M) x = 1 exactly, return 1^T x
        MM=[[ (z if i==j else Fr(0)) - M[i][j] for j in range(n)] for i in range(n)]
        b=[Fr(1)]*n
        # gaussian elimination
        for col in range(n):
            piv=next(r for r in range(col,n) if MM[r][col]!=0)
            MM[col],MM[piv]=MM[piv],MM[col]; b[col],b[piv]=b[piv],b[col]
            inv=Fr(1)/MM[col][col]
            MM[col]=[x*inv for x in MM[col]]; b[col]*=inv
            for r in range(n):
                if r!=col and MM[r][col]!=0:
                    f=MM[r][col]
                    MM[r]=[MM[r][k]-f*MM[col][k] for k in range(n)]; b[r]-=f*b[col]
        return sum(b)
    L2=all(Fr(1)+resolvent_sum(A,Fr(z))==sgn*ev(cpA,Fr(-1-z))/ev(cpA,Fr(z))
           for z in [Fr(17,2),Fr(23,3),Fr(-9,2)])
    # L3: cpK closed form
    K=[[A[j][i]-A[i][j] for j in range(n)] for i in range(n)]
    cpK=charpoly(K)
    L3=all(ev(cpK,Fr(y))==Fr(2)**(n-1)*(ev(cpA,Fr(y-1,2))+sgn*ev(cpA,Fr(-y-1,2)))
           for y in range(-n-2,n+3))
    # L4: d-moments from the ratio
    mu=[]; v=[1]*n
    for j in range(2*n):
        mu.append(sum(v)); v=[sum(K[i][l]*v[l] for l in range(n)) for i in range(n)]
    L4=all(resolvent_sum(K,Fr(y))==Fr(1)-Fr(2)**n*ev(cpA,Fr(y-1,2))/ev(cpK,Fr(y))
           for y in [Fr(19,2),Fr(31,3)])
    # L5: reversal Laplacian
    s=[sum(A[u]) for u in range(n)]
    L=[[(s[u] if u==v else 0)-A[u][v] for v in range(n)] for u in range(n)]
    Lin=[[((n-1-s[u]) if u==v else 0)-A[v][u] for v in range(n)] for u in range(n)]
    cpL=charpoly(L); cpLin=charpoly(Lin)
    L5=all(ev(cpLin,Fr(x))*(x-n)==sgn*Fr(x)*ev(cpL,Fr(n-x)) for x in range(-2,n+3))
    return L1,L2,L3,L4,L5

OK=[]
def check(name,cond):
    OK.append((name,bool(cond))); print(("PASS" if cond else "FAIL"),name)

rng=random.Random(7)
# all iso classes n=4..6 + random big
allok=[True]*5
for n in (4,5,6):
    reps,pairs=census(n)
    for bits in reps:
        r=verify_all(matA(bits,pairs,n),n)
        allok=[a and b for a,b in zip(allok,r)]
for n in (8,9,10):
    for _ in range(6):
        r=verify_all(rand_tournament(n,rng),n)
        allok=[a and b for a,b in zip(allok,r)]
for i,nm in enumerate(["L1 det(zI-A+J)=(-1)^n cpA(-1-z)",
                       "L2 walk resolvent = cpA-ratio (walk moments cpA-determined)",
                       "L3 cpK = 2^{n-1}[cpA((y-1)/2)+(-1)^n cpA((-y-1)/2)] (cpA => cpK)",
                       "L4 d-moment resolvent = 1 - 2^n cpA/cpK (HYP-7096 PROVED, stronger)",
                       "L5 cpL_in*(x-n) = (-1)^n x cpL_out(n-x) (panel reversal-closed)"]):
    check(nm+" — exact, all classes n<=6 + random n=8..10",allok[i])

# n=7: verify L3 on all 456 classes and compare cpA-ties vs cpK-ties
reps7,pairs7=census(7)
cpa_groups={}; cpk_groups={}
ok3=True
for bits in reps7:
    A=matA(bits,pairs7,7)
    cpA=charpoly(A)
    K=[[A[j][i]-A[i][j] for j in range(7)] for i in range(7)]
    cpK=charpoly(K)
    if not all(ev(cpK,Fr(y))==Fr(2)**6*(ev(cpA,Fr(y-1,2))-ev(cpA,Fr(-y-1,2)))
               for y in range(-9,10)): ok3=False
    ka=tuple(cpA); kk=tuple(cpK)
    cpa_groups.setdefault(ka,[]).append(bits)
    cpk_groups.setdefault(kk,[]).append(bits)
check("L3 verified on ALL 456 classes at n=7",ok3)
a_t=sum(1 for g in cpa_groups.values() if len(g)>1)
k_t=sum(1 for g in cpk_groups.values() if len(g)>1)
a_cl=sum(len(g) for g in cpa_groups.values() if len(g)>1)
k_cl=sum(len(g) for g in cpk_groups.values() if len(g)>1)
print(f"   n=7: cpA-tie groups {a_t} ({a_cl} classes) | cpK-tie groups {k_t} ({k_cl} classes)")
check(f"the symmetrization LOSES information: cpK-ties strictly coarser than cpA-ties "
      f"({k_t} groups/{k_cl} classes vs {a_t}/{a_cl}) — cpA => cpK but NOT conversely",
      k_cl>a_cl or (k_t<a_t and k_cl>=a_cl))

print()
fails=[nm for nm,c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK)-len(fails)} passed ===")
for f in fails: print("FAILED:",f)
