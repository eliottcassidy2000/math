"""
lrc14_tourmap_diophantine_m7robust_kps-S2-wf.py

ROBUSTNESS CHECK for the one interesting signal found:
M7 at n=5 argmax-image = 11/12 (misses class #3cyc=1, scores=(0,1,3,3,3))
while allcand (binding crossings) and unconstrained both = 12/12.

Question: is the argmax-miss ROBUST as we widen the speed range, or does it
just fill in?  If it persists across maxspeed 12,14,16,18, that is a genuine
(weak) "the lonely OPTIMUM tau cannot express this class" signal -- a soft
loneliness->tournament obstruction. If it fills in, M7 is just slow to saturate
and there is NO forbidding.

We scan argmax-image of M7 at n=5 for maxspeed in {12,14,16,18,20}.
We also print, when the class is realized, the SMALLEST speed set + argmax tau
that first realizes it, to show how marginal it is.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def signed_frac(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else r-1
def cont_frac(p,q):
    a=[]
    while q: a.append(p//q); p,q=q,p%q
    return a
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def g(S,t): return min(nrm(v*t) for v in S)
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def gcd_list(xs):
    g0=0
    for x in xs: g0=gcd(g0,x)
    return g0

n=5
PAIRS=list(combinations(range(n),2))
def arcs_to_mask(arcs):
    idx={pr:k for k,pr in enumerate(PAIRS)}; mask=0
    for (i,j) in arcs:
        a,b=(i,j) if i<j else (j,i); k=idx[(a,b)]
        if (i,j)==(a,b): mask|=(1<<k)
    return mask
def mask_to_adj(mask):
    A=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(PAIRS):
        if (mask>>k)&1: A[i][j]=1
        else: A[j][i]=1
    return A
def adj_num_3cyc(A):
    c=0
    for a,b,cc in combinations(range(n),3):
        outs={a:0,b:0,cc:0}
        for (x,y) in combinations([a,b,cc],2):
            if A[x][y]==1: outs[x]+=1
            else: outs[y]+=1
        if sorted(outs.values())==[1,1,1]: c+=1
    return c
def adj_score(A): return tuple(sorted(sum(A[i]) for i in range(n)))
def relabel_mask(mask,perm):
    A=mask_to_adj(mask); B=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n): B[perm[i]][perm[j]]=A[i][j]
    m2=0
    for k,(i,j) in enumerate(PAIRS):
        if B[i][j]==1: m2|=(1<<k)
    return m2
PERMS=list(permutations(range(n)))
def class_key(mask):
    return min(relabel_mask(mask,p) for p in PERMS)

def method7(S,tau):
    dist=[nrm(v*tau) for v in S]
    side=[1 if signed_frac(v*tau)>0 else (-1 if signed_frac(v*tau)<0 else 0) for v in S]
    arcs=[]
    for i in range(n):
        for j in range(i+1,n):
            if side[i]!=side[j]:
                if dist[i]==dist[j]:
                    a,b=max(S[i],S[j]),min(S[i],S[j]); P=sum(cont_frac(a,b))
                    big=i if S[i]>S[j] else j; sml=j if big==i else i
                    arcs.append((sml,big) if P%2==0 else (big,sml)); continue
                arcs.append((i,j) if dist[i]<dist[j] else (j,i))
            else:
                if dist[i]==dist[j]:
                    a,b=max(S[i],S[j]),min(S[i],S[j]); P=sum(cont_frac(a,b))
                    big=i if S[i]>S[j] else j; sml=j if big==i else i
                    arcs.append((sml,big) if P%2==0 else (big,sml)); continue
                arcs.append((i,j) if dist[i]>dist[j] else (j,i))
    return arcs

# build identity of the target class (0,1,3,3,3) with 1 three-cycle
def find_target_key():
    # construct one such tournament: score (0,1,3,3,3) means a sink (score0),
    # then a near-transitive. Build transitive then flip one arc to make c3=1.
    # transitive scores would be (0,1,2,3,4). We need (0,1,3,3,3) with 1 cycle.
    # brute search:
    from itertools import product as prod
    for bits in prod([0,1],repeat=len(PAIRS)):
        arcs=[]
        for (bt,(i,j)) in zip(bits,PAIRS):
            arcs.append((i,j) if bt==0 else (j,i))
        A=mask_to_adj(arcs_to_mask(arcs))
        if adj_score(A)==(0,1,3,3,3) and adj_num_3cyc(A)==1:
            return class_key(arcs_to_mask(arcs))
    return None

def main():
    target=find_target_key()
    print(f"Target class (0,1,3,3,3),#3cyc=1 key found: {target is not None}\n")
    print("M7 argmax-image at n=5 vs maxspeed; does it ever realize the target?\n")
    for ms in (12,14,16,18,20):
        realized=set()
        first_target=None
        cnt=0
        for S in combinations(range(1,ms+1),n):
            if gcd_list(S)!=1: continue
            cnt+=1
            tau=M(S)[1]
            ck=class_key(arcs_to_mask(method7(S,tau)))
            realized.add(ck)
            if ck==target and first_target is None:
                first_target=(S,tau)
        hits_target = target in realized
        print(f" maxspeed={ms}: speedsets={cnt}  argmax-image={len(realized)}/12  "
              f"target-realized={hits_target}"
              + (f"  FIRST at S={first_target[0]} tau={first_target[1]}" if first_target else ""))

if __name__=="__main__":
    main()
