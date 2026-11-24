"""Compute exact ordered path-cover polynomials from odd-cycle packings."""
from itertools import combinations
from functools import cache
from math import gcd
import subprocess

def rows(n,forward):
    T=[0]*n
    for i,j in combinations(range(n),2):
        a,b=(i,j) if forward(i,j) else (j,i)
        T[a]|=1<<b
    return tuple(T)

def paired(n,A,order):
    return rows(2*n,lambda i,j:j<n or (i<n and (order[i]==order[j-n] or A(order[i],order[j-n]))))

def mask(N,m):
    C={(i,j) for k,(i,j) in enumerate((i,j) for i,j in combinations(range(N),2) if j>i+1) if m>>k&1}
    return rows(N,lambda i,j:(i,j) not in C)

def profile(T):
    n=len(T);size=1<<n; full=size-1
    C=[0]*size
    for v in range(n):
        U=(size-1)^((1<<v)-1);start=1<<v
        D={start:{v:1}}
        for b in range(start,size,2*start):
            if b&U!=b or b not in D:continue
            E=D[b]
            if b.bit_count()>=3 and b.bit_count()%2:
                C[b]=sum(c for u,c in E.items() if T[u]&start)
            for u,c in E.items():
                w=T[u]&U&~b
                while w:
                    q=w&-w;w^=q
                    j=q.bit_length()-1;X=D.setdefault(b|q,{})
                    X[j]=X.get(j,0)+c
    @cache
    def res(S):
        if not S:return (1,)
        a=S&-S;R=S^a
        z=[0]+list(res(R))
        if R:
            q=R
            while q:
                if C[q|a]:
                    for j,x in enumerate(res(R^q)):
                        z[j+1]+=2*C[q|a]*x
                q=(q-1)&R
        return tuple(z)
    return res(full)

def evaluate(f,k):return sum(a*k**i for i,a in enumerate(f))

def factors(N):
    if N==1:return "1"
    return subprocess.check_output(['factor',str(N)],text=True).strip().split(':',1)[1].strip()

def fixtures():
    A=lambda i,j:(i-j)%5 in (1,4)
    return {
        '9 cell A':(16,1,3,7,14),
        '9 cell B':(20,1,2,7,14),
        'C5 cyclic':rows(5,lambda i,j:(j-i)%5 in (1,2)),
        'C5 paired':paired(5,A,(0,1,4,3,2)),
        'Petersen sparse':mask(10,4304406596),
        'Petersen compact':rows(10,lambda i,j:((j-i)%5 in (1,2)) if i//5==j//5 else i!=j-5),
        'Petersen hidden':(43198,28988,57464,21680,36128,7232,923,19254,37421,9279,63943,62286,58069,49656,33649,746),
    }

def main():
    for name,T in fixtures().items():
        f=profile(T);e=next(i for i,c in enumerate(f) if c);F=f[e:]
        D=gcd(*(evaluate(F,k) for k in range(len(F))))
        raw=gcd(*(evaluate(f,k) for k in range(len(f))))
        assert raw=={'9 cell A':3,'9 cell B':3,'C5 cyclic':15,'C5 paired':1,'Petersen sparse':1,'Petersen compact':15,'Petersen hidden':3}[name]
        print(name,'profile',f,'minpower',e,'fixed divisors raw/residual',raw,D,flush=True)
        for k in range(1,7):
            z=evaluate(f,k);r=evaluate(F,k)
            print(k,z,' residual',r,factors(r),flush=True)
if __name__=='__main__':main()
