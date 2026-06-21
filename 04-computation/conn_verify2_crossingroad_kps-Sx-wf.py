"""
Second-pass adversarial checks on the crossing/P connection.
- Confirm the SPECIFIC claimed witnesses: n=5 (1,1,2,3,3)->(c3,P)=(3,5)vs(3,1);
  n=6 (1,1,3,3,3,4)->(5,12)vs(5,10); n=7 same-c3 diff-P.
- Check P refines c3 AT FIXED c3 (not just overall): same c3, different P.
- Is P just a relabeling of the directed-4-cycle count or a linear combo of (c3, scores)?
  Test: regress P against (1, c3, sum C(s_i,2)) exactly over full enumeration; if a perfect
  integer/rational affine fit exists, P is NOT new. If residuals nonzero, P is genuinely new.
- 'P affine -1/2 c3 + 15 at n=6 fails': verify no exact affine fit P = a*c3 + b exists.
"""
from itertools import combinations, product
from fractions import Fraction
from collections import defaultdict
from math import comb

def tiles(n): return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]
def adj(n,bits,T):
    A=[[0]*(n+1) for _ in range(n+1)]
    for k in range(n,1,-1): A[k][k-1]=1
    for (a,b),bit in zip(T,bits):
        if bit==0: A[a][b]=1
        else: A[b][a]=1
    return A
def c3(A,n):
    t=0
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            for k in range(j+1,n+1):
                if (A[i][j]+A[i][k],A[j][i]+A[j][k],A[k][i]+A[k][j])==(1,1,1):t+=1
    return t
def scores_unsorted(A,n):
    return [sum(A[i][j] for j in range(1,n+1) if j!=i) for i in range(1,n+1)]
def Pinv(A,n):
    P=0
    for w,x,y,z in combinations(range(1,n+1),4):
        if A[w][y]==A[x][z]: P+=1
    return P
def c3_from_scores(sc,n):
    return comb(n,3)-sum(comb(s,2) for s in sc)

out=[]
def log(*a):
    s=" ".join(str(x) for x in a); print(s); out.append(s)

# witness search at fixed score AND fixed c3
def fixed_c3_diff_P(n):
    T=tiles(n); F=len(T)
    byc3=defaultdict(set)
    for bits in product([0,1],repeat=F):
        A=adj(n,list(bits),T)
        byc3[c3(A,n)].add(Pinv(A,n))
    return {c:sorted(ps) for c,ps in byc3.items()}

log("=== P spread at FIXED c3 (n=5,6) ===")
for n in (5,6):
    d=fixed_c3_diff_P(n)
    for c in sorted(d): log(f"n={n} c3={c}: P in {d[c]}")

# exact affine fit test P = a*c3 + b over full data
def affine_fit_fails(n):
    T=tiles(n); F=len(T)
    data=[]
    for bits in product([0,1],repeat=F):
        A=adj(n,list(bits),T)
        data.append((c3(A,n),Pinv(A,n)))
    # try to fit P = a*c3+b exactly: need all points colinear
    pts=sorted(set(data))
    # check colinearity
    if len(pts)<2: return False,None
    (x0,y0),(x1,y1)=pts[0],pts[1]
    if x1==x0:
        # vertical -> not a function
        return True,"two P at same c3 (not a function)"
    a=Fraction(y1-y0,x1-x0); b=Fraction(y0)-a*x0
    for x,y in pts:
        if a*x+b!=y:
            return True,(a,b,"breaks at",(x,y))
    return False,(a,b)

log("\n=== exact affine fit P=a*c3+b ===")
for n in (5,6):
    fails,info=affine_fit_fails(n)
    log(f"n={n}: affine fit FAILS = {fails} ; info={info}")

# Is P determined by (c3, sum C(s_i,2)) i.e. only score data? Already know not score-determined,
# but check if P is a function of (c3, full sorted score). It is NOT (test2 showed same score diff P).
# Stronger: is P = linear combination of c3 and directed-4-cycle count? Compute d4c.
def directed_4cycles(A,n):
    # count directed 4-cycles i->j->k->l->i over all 4-subsets, all rotations
    cnt=0
    for S in combinations(range(1,n+1),4):
        # count cyclic orderings forming a directed cycle
        verts=list(S)
        # check all 3 distinct 4-cycles on 4 labeled vertices (3 perfect matchings -> 3 cycles)
        from itertools import permutations
        seen=set()
        for perm in permutations(verts):
            a,b,cc,dd=perm
            if A[a][b] and A[b][cc] and A[cc][dd] and A[dd][a]:
                seen.add(frozenset([(a,b),(b,cc),(cc,dd),(dd,a)]))
        cnt+=len(seen)
    return cnt

# regress P on basis (1, c3, d4c) exactly via solving over rationals on full data n<=6
def regress(n,basis_funcs):
    T=tiles(n); F=len(T)
    rows=[]; ys=[]
    for bits in product([0,1],repeat=F):
        A=adj(n,list(bits),T)
        rows.append([f(A,n) for f in basis_funcs]+[1])
        ys.append(Pinv(A,n))
    # least-exact: check if y is in span. Build unique rows.
    seen={}
    consistent=True; conflict=None
    for r,y in zip(rows,ys):
        key=tuple(r)
        if key in seen and seen[key]!=y:
            consistent=False; conflict=(key,seen[key],y); break
        seen[key]=y
    return consistent,conflict,len(seen)

log("\n=== Is P a FUNCTION of (c3, directed-4-cycle count)? ===")
for n in (5,6):
    cons,conf,nrows=regress(n,[c3,directed_4cycles])
    log(f"n={n}: P determined by (c3,d4c)? = {cons} ; #distinct feature rows={nrows} ; conflict={conf}")

with open("05-knowledge/results/conn_verify2_crossingroad_kps-Sx-wf.out","w") as f:
    f.write("\n".join(out))
