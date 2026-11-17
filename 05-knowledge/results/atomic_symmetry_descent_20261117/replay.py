"""Literal controls for symmetry, ordered covers, and exact universal divisors."""
from itertools import product, permutations
from math import gcd, factorial, comb, isqrt
from functools import cache
from collections import Counter

def need(ok):
    if not ok: raise ArithmeticError("certificate failed")

def counts(T):
    n=len(T);f=(1<<n)-1
    need(all(0<=a<=f and not a>>i&1 for i,a in enumerate(T)))
    need(all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(n) for j in range(i)))
    d=[[0]*n for _ in range(f+1)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,f+1):
        for i,v in enumerate(d[s]):
            if not v:continue
            B=T[i]&~s
            while B:
                j=B&-B; B^=j
                d[s|j][j.bit_length()-1]+=v
    return [1]+[sum(L) for L in d[1:]]

def profile(T):
    H=counts(T); f=len(H)-1
    @cache
    def C(s,r):
        if r<=1:return H[s] if r==1 and s else int(r==0 and s==0)
        A=s; out=0
        while A:
            out+=H[A]*C(s^A,r-1); A=(A-1)&s
        return out
    c=tuple(C(f,r) for r in range(1,len(T)+1))
    p=tuple(sum(comb(q,i)*val for i,val in enumerate(c,1)) for q in (1,2,3))
    D=gcd(*(factorial(i)*x for i,x in enumerate(c)))
    return p,c,D

def cyclic(n,S):
    return tuple(sum(1<<((i+a)%n) for a in S) for i in range(n))

def strong(T):
    f=(1<<len(T))-1
    for i in range(len(T)):
        s=1<<i
        for _ in T:
            for j,x in enumerate(T):
                if s>>j&1:s|=x
        if s!=f:return False
    return True

def substitute(T,i,B):
    n=len(T); k=len(B); groups=[(j,) if j<i else (j+k-1,) for j in range(n)]
    groups[i]=tuple(range(i,i+k)); out=[]
    for p in range(n):
        ext=sum(1<<z for q in range(n) if T[p]>>q&1 for z in groups[q])
        if p==i:out.extend(ext+(b<<i) for b in B)
        else:out.append(ext)
    return tuple(out)

def matching(r):
    n=2*r-1
    def arc(i,j):
        if i==0:return j<r
        if j==0:return i>=r
        if (i<r)==(j<r):return i<j
        return j==i+r-1 if i<r else i!=j+r-1
    return tuple(sum(1<<j for j in range(n) if i!=j and arc(i,j)) for i in range(n))

def slots(T,w):
    ans=Counter()
    for P in permutations(i for i in range(len(T)) if i!=w):
        req=0; opt=(T[w]>>P[0]&1)+(T[P[-1]]>>w&1)
        for a,b in zip(P,P[1:]):
            can=(T[a]>>w&1)&(T[w]>>b&1)
            if T[a]>>b&1:opt+=can
            elif can:req+=1
            else:break
        else:
            for j in range(opt+1):
                if req+j:ans[req+j]+=comb(opt,j)
    return ans

def run():
    I=cyclic(7,(1,2,3));P=cyclic(7,(1,2,4))
    bucket=Counter()
    for es in product((-1,1),repeat=3):
        a=cyclic(7,tuple(e*(i+1) for i,e in enumerate(es)))
        p,c,D=profile(a);bucket[p[0],D]+=1
        need(c==(p[0],1638,6678,15120,20160,15120,5040))
    need(bucket=={(175,7):6,(189,63):2})
    for n,S in ((5,(1,2)),(5,(1,-2))):need(profile(cyclic(n,S))[0][0]==15)
    Q=(2,4,9,3)
    need([slots(Q,i) for i in range(4)]==[{1:5,2:2},{1:5,2:3},{1:5,2:3},{1:5,2:2}])
    T=substitute(Q,0,I)
    need(strong(T) and counts(T)[-1]==4151 and all(593%d for d in range(2,isqrt(593)+1)))
    for r in range(2,5):
        X=matching(r); a=slots(X,0)
        need(strong(X) and max(a)==r and a[r]==factorial(r-1))
    examples=(
        ((0,),0),
        ((30,60,120,112,96,65,3),1),
        ((126,92,120,112,96,66,0),2),
        ((8,25,51,116,225,195,135,15),3),
        ((62,124,104,112,100,64,1),4),
        ((110,92,56,112,97,66,4),5),
        ((126,60,120,112,96,64,2),6),
        ((118,108,120,49,98,64,8),7),
    )
    for T,z in examples:
        q,_,D=profile(T)
        need(sum((a%7==0)<<i for i,a in enumerate(q))==z and ((D%7==0)==(z==7)))
        print(z,q,D)
    A=(118,108,120,49,98,64,8)
    B=(2,20,17,7,105,15,47)
    need(profile(A)[0]==(49,812,6363) and profile(A)[2]==7 and strong(A))
    need(profile(B)[0]==(49,872,6843) and profile(B)[2]==1 and strong(B))
    bad=examples[3][0]
    cut=matching(3); total=counts(substitute(cut,0,bad))[-1]
    need(total%7==1 and profile(bad)[0]==(427,6496,49851))
    X=(4,49,90,99,169,197,147,15)
    val,c,d=profile(X)
    need(c==(459,5778,31860,99000,187200,216000,141120,40320) and gcd(*c)==9 and d==27)
    print('593 symmetry anchor:',substitute(Q,0,I))
    print('three-visit separating read:',total)
    print('all passport and orbit fixtures: PASS')

if __name__=='__main__':run()
