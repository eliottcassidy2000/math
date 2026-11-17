from fractions import Fraction as Q
from itertools import combinations
from math import factorial,comb
import json
from pathlib import Path

def ck(x,y):
    if x!=y:raise AssertionError((x,y))
def pars(n):return list(combinations(range(n),2))
def TT(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def rows(n,bits):
    a=[0]*n
    for k,(i,j) in enumerate(pars(n)):
        if bits>>k&1:a[j]|=1<<i
        else:a[i]|=1<<j
    return tuple(a)
def R(n):
    m=n if n%2 else n+1;inds=list(range(m)) if n%2 else list(range(1,m))
    return tuple(sum((0<(j-i)%m<m/2)<<v for v,j in enumerate(inds)) for i in inds)
def count(a):
    n=len(a);dp=[[0]*n for _ in range(1<<n)];h=[1]+[0]*(len(dp)-1)
    for i in range(n):dp[1<<i][i]=1
    for V in range(1,len(dp)):
        h[V]=sum(dp[V])
        for i,z in enumerate(dp[V]):
            if not z:continue
            E=a[i]&~V
            while E:
                b=E&-E;E-=b;dp[V|b][b.bit_length()-1]+=z
    return h
def P2(a):
    h=count(a);f=len(h)-1
    return sum(h[s]*h[f^s] for s in range(f+1))
def sub(a,S):
    return tuple(sum((a[i]>>j&1)<<k for k,j in enumerate(S)) for i in S)
def second(n):
    Z=1
    for k in range(1,(n-1)//2+1):
        for c in range(1,min(k,n-2*k)+1):
            Z+=Q(2**c*comb(k-1,c-1)*comb(n-2*k,c)*factorial(n-2*k),factorial(n))
    return Z
def walk(n):
    a=list(TT(n));target=R(n)
    yield tuple(a),None
    for i,j in pars(n):
        if (a[i]^target[i])>>j&1:
            a[i]^=1<<j;a[j]^=1<<i
            yield tuple(a),(i,j)
    ck(tuple(a),target)
def main():
    for n in range(1,6):
        nlab=1<<len(pars(n));X=[count(rows(n,k))[-1] for k in range(nlab)]
        mu=Q(factorial(n),1<<(n-1))
        ck(Q(sum(X),nlab),mu)
        ck(Q(sum(x*x for x in X),nlab)/(mu*mu),second(n))
        if n>=2:
            z,s2=0,0
            for k,x in enumerate(X):
                deg=sum((2*a.bit_count()-(n-1))**2 for a in rows(n,k))
                p=2-Q(deg,n*(n-1)); R=x/mu-p
                z+=R*R;s2+=R*(p-1)
            ck(z/nlab,second(n)-1-Q(2*(n-2),n*(n-1)))
            ck(s2,0)
    print('exact second moments n=1..5:',*(str(second(n)) for n in range(1,6)))
    path=Path(__file__).resolve().parents[2]/'birth-atlas-now'/'atlas.json'
    known={1:1,2:1,3:3,4:5,5:15,6:45,7:189,8:661}
    for n in range(2,11):
        history=[];old_h=1
        for a,arc in walk(n):
            h=count(a)[-1];history.append(h)
            if arc:
                q=P2(sub(a,[i for i in range(n) if i not in arc]))
                ck(abs(old_h-h)<=q,True)
            old_h=h
        if path.exists() and n<=9:
            atlas=json.loads(path.read_text())
            bound=max(int(h) for h,j in atlas['births'].items() if j['order']<n)
            ck(h>bound if n>=3 else h==bound,True)
            print('n, circular endpoint, earlier maximum from supplied atlas, witnessed new values on path:',n,h,bound,len({x for x in history if x>bound}))
    print('all 9 deterministic paths independently verify the corrected two-cover edge bound')
if __name__=='__main__':main()
