#!/usr/bin/env python3
"""Finite, exact controls for slot modules and pointed Young-lattice transfer."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial, comb
import json


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def tournament(n, bits):
    out=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if bits>>k&1: out[i]|=1<<j
        else: out[j]|=1<<i
    return tuple(out)


def induced(a, vs):
    return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)


@lru_cache(None)
def hstarts(a):
    n=len(a)
    if n==0: return ()
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1 and mask!=(1<<v):
                rest=mask^(1<<v)
                dp[mask][v]=sum(dp[rest][w] for w in range(n) if a[v]>>w&1)
    return tuple(dp[-1])


@lru_cache(None)
def covers(a):
    n=len(a)
    if not n: return (1,)
    out=[0]*(n+1); full=(1<<n)-1
    for sub in range(1,1<<n,2):
        h=sum(hstarts(induced(a,[i for i in range(n) if sub>>i&1])))
        q=covers(induced(a,[i for i in range(n) if not(sub>>i&1)]))
        for k,v in enumerate(q): out[k+1]+=h*v
    return tuple(out)


def startcovers(a,q):
    n=len(a);out=[0]*(n+1)
    for sub in range(1,1<<n):
        if not(sub>>q&1):continue
        vs=[i for i in range(n) if sub>>i&1]
        h=hstarts(induced(a,vs))[vs.index(q)]
        p=covers(induced(a,[i for i in range(n) if not(sub>>i&1)]))
        for k,v in enumerate(p):out[k+1]+=h*v
    return tuple(out)


def pointed(a,q):
    pred=[v for v in range(len(a)) if a[v]>>q&1]
    require(len(pred)==1,'invalid guard')
    p=covers(a);s=startcovers(a,q)
    vs=[v for v in range(len(a)) if v!=q]
    d=startcovers(induced(a,vs),vs.index(pred[0]))
    return s,tuple(x-y for x,y in zip(p,s)),d


def one_attach(a,q):
    n=len(a);out=list(a)+[0];out[q]|=1<<n
    out[n]=((1<<n)-1)^(1<<q)
    return tuple(out)


def triangle(a):
    n=len(a)
    return tuple([v|(1<<n) for v in a]+[1<<(n+1),(1<<n)-1])


def ordinal(a,b):
    n=len(a);m=len(b)
    return tuple([v|(((1<<m)-1)<<n) for v in a]+[v<<n for v in b])


def slots(p,r):
    falling=1;out=0
    for k,v in enumerate(p):
        if k:falling*=r-k+1
        out+=v*falling
    return out


@lru_cache(None)
def partitions(n, cap=None):
    if not n:return ((),)
    cap=n if cap is None else min(cap,n)
    return tuple((k,)+tail for k in range(cap,0,-1) for tail in partitions(n-k,k))


@lru_cache(None)
def predecessors(lam):
    out=[]
    for i,v in enumerate(lam):
        if i+1==len(lam) or v>lam[i+1]:
            q=list(lam);q[i]-=1
            if not q[i]:q.pop(i)
            out.append(tuple(q))
    return tuple(out)


@lru_cache(None)
def skew_tableaux(lam, row):
    """f^(lambda/(row)), counted by corner removal; row=0 means empty."""
    if sum(lam)<row:return 0
    if sum(lam)==row:return int(lam==((row,) if row else ()))
    return sum(skew_tableaux(mu,row) for mu in predecessors(lam))


def multiplicities(p,r):
    return tuple(sum(p[k]*skew_tableaux(lam,r-k) for k in range(min(r,len(p)-1)+1)) for lam in partitions(r))


def matrix(r):
    parts=partitions(r)
    return tuple(tuple(len(set(predecessors(a))&set(predecessors(b))) for b in parts) for a in parts)


def matvec(K,v):return tuple(sum(a*b for a,b in zip(row,v)) for row in K)
def add(*vs):return tuple(map(sum,zip(*vs)))
def neg(v):return tuple(-x for x in v)


def transfer(v,r,move):
    s,u,d=v;K=matrix(r);ks=matvec(K,s);ku=matvec(K,u);kd=matvec(K,d)
    if move=='C':return add(ks,neg(s),ku),add(s,u,d),s
    return add(ks,d),add(ku,u),add(kd,neg(d))


def compact(p):
    return (p[1] if len(p)>1 else 0, sum(p[2:4]))


def compact_transfer(v,move):
    (hs,bs),(hu,bu),(hd,bd)=v
    if move=='C':return ((hu,hs+2*bs+hu+3*bu),(hs+hu+hd,bs+bu+bd),(hs,bs))
    return ((hs+hd,hs+3*bs+bd),(2*hu,hu+4*bu),(0,hd+2*bd))


def odd_cycles(a):
    n=len(a);out=[]
    for k in range(3,n+1,2):
        for vs in combinations(range(n),k):
            for tail in permutations(vs[1:]):
                cyc=(vs[0],)+tail
                if all(a[cyc[i]]>>cyc[(i+1)%k]&1 for i in range(k)):
                    out.append((sum(1<<v for v in cyc),k))
    return out


def odd_cycle_monomials(a):
    cyc=odd_cycles(a);out=Counter();n=len(a)
    def visit(i,used,excess,weight):
        if i==len(cyc):out[n-excess]+=weight;return
        visit(i+1,used,excess,weight)
        mask,k=cyc[i]
        if not(mask&used):visit(i+1,used|mask,excess+k-1,2*weight)
    visit(0,0,0,1)
    return out


@lru_cache(None)
def stirling(n,k):
    if n==0:return int(k==0)
    if k==0:return 0
    return stirling(n-1,k-1)+k*stirling(n-1,k)


def fixed_slots(a,r,sigma):
    """Direct path-and-cut enumeration, independent of subset DP."""
    n=len(a);fixed=0
    for order in permutations(range(n)):
        for cutbits in range(1<<(n-1)):
            blocks=[];start=0;valid=True
            for i in range(n-1):
                if cutbits>>i&1:blocks.append(order[start:i+1]);start=i+1
                elif not(a[order[i]]>>order[i+1]&1):valid=False;break
            if not valid:continue
            blocks.append(order[start:])
            # Canonical order chooses one occurrence of each unordered cover.
            if tuple(blocks)!=tuple(sorted(blocks)):continue
            k=len(blocks)
            if k>r:continue
            for assignment in permutations(range(r),k):
                if all(sigma[j]==j for j in assignment):fixed+=1
    return fixed


def main():
    report={}
    # Character statement and S_r Specht branching, small independent controls.
    character_cases=0
    for n in range(1,4):
        for bits in range(1<<(n*(n-1)//2)):
            a=tournament(n,bits);p=covers(a)
            for r in range(1,5):
                ms=multiplicities(p,r)
                dims=[skew_tableaux(lam,0) for lam in partitions(r)]
                require(sum(x*y for x,y in zip(ms,dims))==slots(p,r),'dimension decomposition')
                for sigma in permutations(range(r)):
                    f=sum(i==v for i,v in enumerate(sigma))
                    require(fixed_slots(a,r,sigma)==slots(p,f),'fixed-cover character')
                    character_cases+=1
    report['direct_character_cases']=character_cases
    # Young-lattice transfer against independently counted literal graphs.
    guard_cases=0
    for n in range(2,5):
        for bits in range(1<<(n*(n-1)//2)):
            a=tournament(n,bits)
            for q in range(n):
                if sum(v>>q&1 for v in a)!=1:continue
                state=pointed(a,q);b=one_attach(a,q)
                for r in range(2,7):
                    v=tuple(multiplicities(p,r) for p in state)
                    for move,newq in [('C',n),('D',q)]:
                        actual=tuple(multiplicities(p,r) for p in pointed(b,newq))
                        require(transfer(v,r,move)==actual,'pointed Specht transfer')
                        if r==3:
                            compact_actual=tuple(compact(p) for p in pointed(b,newq))
                            require(compact_transfer(tuple(compact(p) for p in state),move)==compact_actual,'compact six-coordinate transfer')
                        guard_cases+=1
    triangle_cases=0
    for n in range(1,5):
        for bits in range(1<<(n*(n-1)//2)):
            a=tournament(n,bits);p=covers(a)
            ps,pu,pd=pointed(triangle(a),n+1)
            for r in range(1,7):
                require(slots(pd,r)==r*slots(p,r-1),'triangle induction')
                require(slots(pu,r)==r*slots(p,r+1),'triangle restriction')
                require(slots(ps,r)+slots(pd,r)==r*r*slots(p,r),'triangle complement')
                triangle_cases+=1
    report['triangle_functor_cases']=triangle_cases
    report['pointed_Specht_cases']=guard_cases
    report['compact_six_coordinate_cases']=guard_cases//5
    for r in range(2,11):
        K=matrix(r)
        require(all(K[i][j]-(i==j)>=0 for i in range(len(K)) for j in range(len(K))),'K-I positivity')
    # Odd-cycle tensor powers against independent path-cover DP, all tournaments <=5.
    ocf=0
    for n in range(1,6):
        for bits in range(1<<(n*(n-1)//2)):
            a=tournament(n,bits);coeff=odd_cycle_monomials(a);p=covers(a)
            actual=tuple(sum(weight*stirling(b,k) for b,weight in coeff.items()) for k in range(n+1))
            require(actual==p,'Stirling-positive cycle lift')
            ocf+=1
    report['cycle_tensor_cases']=ocf
    # Tensor/ordinal action including non-isomorphic inputs.
    small=[tournament(n,b) for n in range(1,4) for b in range(1<<(n*(n-1)//2))]
    for a in small:
        for b in small:
            for r in range(6):
                require(slots(covers(ordinal(a,b)),r)==slots(covers(a),r)*slots(covers(b),r),'ordinal tensor character')
    report['ordinal_character_cases']=len(small)**2*6
    # Literal C-chain and exact S3 decomposition; recurrence through n=9.
    a=(2,0);q=1;records=[];profiles=[]
    for j in range(8):
        p=covers(a);profiles.append(p)
        rec={'order':len(a),'H':p[1],'F3':slots(p,3),'S3_multiplicities':multiplicities(p,3)}
        records.append(rec)
        if j>=3:
            for r in range(1,9):
                require(slots(profiles[j],r)==r*slots(profiles[j-1],r)+slots(profiles[j-2],r)+r*slots(profiles[j-3],r),'full character recurrence')
        a=one_attach(a,q);q=len(a)-1
    report['C_chain']=records
    # R5 obstruction pertains to vertex-degree symmetric function, not slot module.
    a=tuple(sum(1<<j for j in range(5) if (j-i)%5 in (1,2)) for i in range(5))
    cyc=Counter(k for _,k in odd_cycles(a))
    require(cyc=={3:5,5:2},'R5 directed cycle counts')
    schur32=5+2*cyc[3]*(-1)+2*cyc[5]*0
    require(schur32==-5,'R5 Schur obstruction')
    report['R5']={'cycles':dict(cyc),'cover_profile':covers(a),'s32_coefficient':schur32}
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
