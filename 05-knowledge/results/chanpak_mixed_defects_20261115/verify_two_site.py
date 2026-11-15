"""All-order two-site rectangle check, after the exact finite reduction in the note."""
from itertools import product, combinations, permutations
from math import gcd

def c(d):return 1 if d<=2 else 1<<(d-2)
def inner(w):return sum(c(j-i)*w[i]*(1-w[j]) for i in range(len(w)) for j in range(i+1,len(w)))
def kernel(X,Y,s,t):
    return sum(c(abs(i-j))*a*b for i,a in zip(X,(1,s)) for j,b in zip(Y,(1,t)))

def exgcd(a,b):
    if not b:return abs(a), (1 if a>=0 else -1), 0
    g,x,y=exgcd(b,a%b)
    return g,y,x-(a//b)*y

def integer_solution(A,B,F,C,D,G):
    # A x + B y=F; C x + D y=G, x,y nonnegative ints.
    det=A*D-B*C
    if det:
        x=F*D-B*G;y=A*G-F*C
        return not (x%det or y%det) and x//det>=0 and y//det>=0
    if A*G!=C*F:return False
    g,x,y=exgcd(A,B)
    if F%g:return False
    x*=F//g;y*=F//g
    u=B//g;v=-A//g
    lo=None;hi=None
    for p,d in ((x,u),(y,v)):
        if d>0:lo=-(p//d) if lo is None else max(lo,-(p//d))
        if d<0:hi=p//(-d) if hi is None else min(hi,p//(-d))
        if d==0 and p<0:return False
    return lo is None or hi is None or lo<=hi

def require(x):
    if not x:raise AssertionError

for r,h,t in product(range(1,15),repeat=3):
    a,b,d,e=0,r,r+h,r+h+t
    for k,(X,Y) in enumerate((((a,b),(d,e)),((a,d),(b,e)),((a,e),(b,d)))):
        for s1,s2 in product((1,-1),repeat=2):
            mix=kernel(X,Y,s1,s2)
            expected=k==2 and ((s1*s2==-1 and r==t) or (s1==s2==-1 and r==h==t==1))
            require((mix==0)==expected)
            if b+h+t>11 and mix:require(abs(mix)>4)

cases={}
for p,q,r,s in permutations(range(4)):
    cases.setdefault(p+s-q-r,set()).add((q-p,r-p))
checked=0
for r,h,t in product(range(1,10),repeat=3):
    L=r+h+t
    if L>11:continue
    a,b,d,e=0,r,r+h,L
    gap=[i for i in range(L+1) if i not in (a,b,d,e)]
    for X,Y in (((a,b),(d,e)),((a,d),(b,e)),((a,e),(b,d))):
        for s1,s2 in product((1,-1),repeat=2):
            mix=-kernel(X,Y,s1,s2)
            if not mix or mix not in cases:continue
            def controls(T,s):
                return (-sum(k*(1<<i) for i,k in zip(T,(1,s))),
                        sum(k*(1<<(L-i)) for i,k in zip(T,(1,s))),
                       -sum(k*c(i+1) for i,k in zip(T,(1,s))),
                        sum(k*c(L-i+1) for i,k in zip(T,(1,s))))
            A,B,e1,f1=controls(X,s1);C,D,e2,f2=controls(Y,s2)
            for w0 in product(range(2),repeat=len(gap)):
                w=[0]*(L+1)
                for i,k in zip(gap,w0):w[i]=k
                w[X[1]]=int(s1==-1);w[Y[1]]=int(s2==-1)
                z=inner(w)
                for i in X:w[i]^=1
                delta1=inner(w)-z
                for i in X:w[i]^=1
                for i in Y:w[i]^=1
                delta2=inner(w)-z
                for u0,v0 in product(range(2),repeat=2):
                    for x,y in cases[mix]:
                        checked+=1
                        F=x-delta1-e1*u0-f1*v0
                        G=y-delta2-e2*u0-f2*v0
                        if integer_solution(A,B,F,C,D,G):
                            raise AssertionError((X,Y,s1,s2,w0,u0,v0,x,y))
print(f"{checked} residual two-site systems; no four consecutive odds")
