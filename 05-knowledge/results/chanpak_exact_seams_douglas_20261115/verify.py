#!/usr/bin/env python3
"""Independent controls for the rigidity and unique-cycle exclusion."""
from math import comb
from itertools import product, combinations, chain
from fractions import Fraction as Q

def merge(D):
    p,q=len(D),len(D[0]); A=[[0]*(q+1) for _ in range(p+1)]; B=[[0]*(q+1) for _ in range(p+1)]
    for i in range(1,p+1): A[i][0]=1
    for j in range(1,q+1): B[0][j]=1
    for i,j in product(range(1,p+1),range(1,q+1)):
        A[i][j]=A[i-1][j]+(not D[i-1][j-1])*B[i-1][j]
        B[i][j]=B[i][j-1]+D[i-1][j-1]*A[i][j-1]
    F=[[x+y for x,y in zip(a,b)] for a,b in zip(A,B)]
    F[0][0]=1
    return F

def ferrers_control():
    checks=0
    for p,q in [(2,2),(2,3),(3,2),(3,3),(2,5),(5,2),(3,4),(4,3),(4,4)]:
        smallest=float('inf')
        for bits in product((0,1),repeat=p*q):
            D=[bits[i*q:(i+1)*q] for i in range(p)]
            R=[[D[p-i-1][q-j-1] for i in range(p)] for j in range(q)]
            F,G=merge(D),merge(R)
            L=sum(F[i][j]*G[q-j][p-i]-1 for i,j in product(range(p+1),range(q+1)))
            good=all(D[i][j]<=D[i][j+1] for i,j in product(range(p),range(q-1))) and all(D[i][j]>=D[i+1][j] for i,j in product(range(p-1),range(q)))
            if good: assert L==0
            else: smallest=min(smallest,L)
            checks+=1
        assert smallest==p+q-1, (p,q,smallest)
    return checks

def hp(T):
    n=len(T); q=1<<n
    dp=[[0]*n for _ in range(q)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,q):
        for i in range(n):
            if dp[s][i]:
                for j in range(n):
                    if not(s>>j&1) and T[i][j]: dp[s|1<<j][j]+=dp[s][i]
    return sum(dp[-1])

def backward_spine(g):
    n=1+sum(1+x for x in g); S=[0]
    for x in g:S.append(S[-1]+1+x)
    T=[[int(i<j) for j in range(n)] for i in range(n)]
    for i,j in zip(S,S[1:]): T[i][j],T[j][i]=0,1
    return T

def douglas(g,r):
    # g lists natural gap indices of transitive D in its given order.
    n=r+1+len(g); A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if j<=r: bit=int(j!=i+1)
            elif i>r: bit=1
            else: bit=int(i<=g[j-r-1])
            A[i][j],A[j][i]=bit,1-bit
    return A

def gap_profile(g,r):
    return tuple(g.count(i) for i in range(r))

def compositions(n):
    if n==0: yield ()
    for i in range(n):
        for tail in compositions(n-i-1): yield (i,)+tail

def transfer(g):
    H,S,Z=1,1,1
    for a in g:
        t=(Z if a==0 else 0)
        H,S,Z=H+(1<<a)*S-t,H+(1<<a)*S,H-t
    return H

def floor_d(d):
    return (3,5,9)[(d-2)%3]*5**((d-2)//3)

def cone():
    def V(k,t,z):
        H=S=1; Z=z
        for _ in range(k): H,S,Z=H+S-Z,H+S,H-Z
        return H+t*S
    for k in (0,1,2):
        for z in (Q(0),Q(1,3)):
            for t in (0,1):
                assert V(k+3,t,z)-5*V(k,t,z)>=0
            assert V(k,2,0)>=V(k,2,Q(1,3))
            assert V(k,3,0)-V(k,3,Q(1,3))>=V(k,2,0)-V(k,2,Q(1,3))
    R=min(V(k,1<<b,Q(1,3))/floor_d(k+b+1) for k,b in product((1,2,3),repeat=2))
    R2=min((1+(1<<a))*V(k,1<<b,Q(1,1+(1<<a)))/floor_d(k+a+b+2) for k,a,b in product((1,2,3),repeat=3))
    assert (R,R2)==(Q(29,27),Q(49,45))

def spine_control():
    exact=0; order=0
    for n in range(3,11):
        for g in compositions(n-1):
            assert hp(backward_spine(g))==transfer(g)
            h=transfer(g); f=floor_d(n-1)
            assert h==f or 25*h>=27*f
            exact+=1
    for r in range(2,7):
        for m in range(min(9-r,5)):
            for g in product(range(r),repeat=m):
                if any(g[i]>g[j]+1 for i,j in combinations(range(m),2)):continue
                T=douglas(g,r)
                h=hp(T); s=transfer(gap_profile(g,r))
                assert h>=s
                for i,j in combinations(range(m),2):
                    if g[i]>g[j]:
                        M={g[i],r+1+i,r+1+j}
                        U=[[x for k,x in enumerate(row) if k not in M] for l,row in enumerate(T) if l not in M]
                        assert h>=s+2*hp(U),(r,g,i,j,h,s,hp(U))
                order+=1
    return exact,order

def n8_boundary():
    hist={}
    for r in range(2,8):
        for g in product(range(r),repeat=7-r):
            if any(g[i]>g[j]+1 for i,j in combinations(range(len(g)),2)):continue
            h=hp(douglas(g,r))
            hist[h]=hist.get(h,0)+1
    assert sum(hist.values())==144 and sorted(hist)[:2]==[45,49], hist
    assert hp(backward_spine([1,3]))==27
    return hist

if __name__=='__main__':
    print(ferrers_control(),'cross matrices; abstract rigidity sharp at p+q-1 in each shape',flush=True)
    cone()
    print(spine_control(),'backward-spine and Douglas sorting controls',flush=True)
    print('order-eight Douglas histogram:',sorted(n8_boundary().items()),flush=True)
