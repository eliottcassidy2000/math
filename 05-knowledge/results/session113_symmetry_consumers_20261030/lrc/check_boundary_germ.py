from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce

CHECKS=0
def check(ok, label):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(label)

def survivor(V,r):
    S=[(F(0),F(1))]
    for v in sorted(V):
        A=[((k+r)/v,(k+1-r)/v) for k in range(v)]
        T=[]; i=j=0
        while i<len(S) and j<len(A):
            a,b=S[i]; c,d=A[j]
            lo=max(a,c); hi=min(b,d)
            if lo<=hi:
                if T and lo<=T[-1][1]: T[-1]=(T[-1][0],max(hi,T[-1][1]))
                else: T.append((lo,hi))
            if b<d:i+=1
            elif d<b:j+=1
            else:i+=1;j+=1
        S=T
    return S

def event_height(V):
    return max([2]+[(v+w)//gcd(v,w) for v,w in combinations(V,2)])

def mass(S):return sum((b-a for a,b in S),F(0))
def owners(V,t,r,sign):
    target=r if sign==1 else 1-r
    return tuple(v for v in V if (v*t)%1==target)
def germ(V,r):
    S=survivor(V,r); B=F(0); marked=[]
    for a,b in S:
        P=owners(V,a,r,1); N=owners(V,b,r,-1)
        if not P or not N:raise RuntimeError('lost endpoint owner')
        B+=F(1,max(P))+F(1,max(N))
        marked.append((a,b,P,N))
    return mass(S),B,marked

def audit(V,r):
    A,B,marked=germ(V,r)
    eta=min(r,F(1,r.denominator*event_height(V)))
    for eps in (eta/2,eta/3,eta*F(9,10)):
        S=survivor(V,r-eps)
        check(mass(S)==A+eps*B,('affinity',V,r,eps,A,B))
        check(bool(S)==bool(marked),('finite bridge',V,r))
        check((mass(S)>=2*eps*reduce(gcd,V)/max(V))==bool(marked),('finite quantitative gate',V,r))
    # Genuine source symmetry: reflection swaps endpoint owners.
    for a,b,P,N in marked:
        check((1-b,1-a,N,P) in marked,('reflection',V,r))
    return A,B,marked

def main():
    rows=0
    for n in range(1,5):
        for V in combinations(range(1,10),n):
            for r in (F(1,14),F(1,7),F(1,3),F(2,5),F(17,50)):
                audit(V,r); rows+=1
    for n in range(3,21):
        V=tuple(range(1,n)); r=F(1,n)
        A,B,S=audit(V,r)
        expected=2*sum((F(1,a) for a in range(1,n) if gcd(a,n)==1),F(0))
        check(A==0 and B==expected,('tight AP',n))
        check(len(S)==sum(gcd(a,n)==1 for a in range(1,n)),('AP count',n))
    A0,B0,S0=audit((1,2),F(1,3))
    A1,B1,S1=audit((1,2,4),F(1,3))
    check(A0==A1==0 and len(S0)==len(S1)==2 and B0==3 and B1==F(3,2),'ties matter')
    V=tuple(range(1,14)); A,B,S=audit(V,F(1,14))
    check(B!=6*(F(1)+F(1,13)),'unit orbit weight is not constant')
    print('tight AP14: mass=',A,'boundary coefficient=',B,'witnesses=',len(S))
    for d in (2,3,7):
        Ad,Bd,Sd=audit(tuple(d*v for v in V),F(1,14))
        check((Ad,Bd,len(Sd))==(A,B,d*len(S)),('dilation',d))
        check(event_height(tuple(d*v for v in V))==event_height(V),('event height dilation',d))
    for V in (tuple(26*k for k in range(1,13))+(339,),tuple(range(1,13))+(5460,)):
        A,B,S=audit(V,F(1,14))
        print('hostile',V,'mass=',A,'coefficient=',B,'components=',len(S))
    r=F(17,50); V=(1,2)
    check(not survivor(V,r),'near miss original')
    check(bool(survivor(V,r-F(1,100))),'unguarded thickening false positive')
    check(not survivor(V,r-F(1,400)),'guarded thickening')
    print('finite rows=',rows,'checks=',CHECKS)
if __name__=='__main__':main()
