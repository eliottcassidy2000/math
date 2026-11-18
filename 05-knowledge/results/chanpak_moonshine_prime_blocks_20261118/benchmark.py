"""Exact prequential source-bank benchmark; freeze train before holdout."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ans_inverse'))
from fixtures import dp,rev,state,forward,H
from math import isqrt
from fractions import Fraction as F
from hashlib import sha256

BASE=((29,11,31,13),(23,2,50,9),(25,2,50,9))
SHELLS=((4096,16384),(16384,65536),(65536,262144))
TRAIN=(1024,4096)
K=8
FROZEN=(
    (3,1,3,1),(5,2,6,3),(5,1,9,3),(9,3,11,5),
    (11,4,12,5),(9,4,10,5),(13,2,18,5),(11,3,17,5))
WITNESSES=(
    ((126,124,120,112,32,64,16),4),
    ((126,124,120,48,96,64,8),5),
    ((126,124,120,48,96,64,8),4),
    ((126,124,40,112,100,64,4),5),
    ((126,124,24,112,96,68,4),6),
    ((126,124,56,112,96,64,4),5),
    ((126,124,88,48,96,68,8),4),
    ((126,124,24,112,96,68,4),5))
UPPER_BASE=BASE+(tuple(7*x for x in BASE[0]),)
EXPECTED=((962,5451,238,1184,1336),
          (3923,21931,714,4098,4642),
          (15463,87585,2546,14458,16458))
CAT_SHA='3989273b1634117975f2a88e3f11c36e5ae263812e5f071f26d549fb25b7edc1'

def profiles(catalog):
    dct={}
    for line in catalog.open():
        T=tuple(map(int,line.split()))
        E=dp(T);F=dp(rev(T));M=len(E)-1;H=[1]+list(map(sum,E[1:]))
        for r in range(len(T)):
            h=H[M];e=E[M][r];t=sum(F[X][r]*H[M^X] for X in range(1,M+1));d=H[M^(1<<r)]
            if h>1 and t>=2*d:dct.setdefault((h,e,t,d),(T,r))
    return dct

def image(S,lo,hi,trace=None,mass=None):
    stack=[(S,'')]
    hits=set()
    while stack:
        (h,e,t,d),word=stack.pop()
        if lo<=h<hi and (not word or word[-1]=='1'):
            hits.add(h)
            if mass is not None:mass[0]+=bool(word)
            if trace is not None:trace.setdefault(h,word)
        if h>=hi or t>=hi:continue
        stack.append(((h,0,2*t,d),word+'0'))
        stack.append(((h+t-e,d,2*t-d,d),word+'1'))
    return hits

def union(B,lo,hi):
    return set().union(*(image(S,lo,hi) for S in B))

def pick_pool(P):
    lo,hi=TRAIN
    covered=union(BASE,lo,hi)
    training={S:image(S,lo,hi) for S in P if S not in BASE}
    bank=[]
    for j in range(K):
        S=min(training,key=lambda S:(-len(training[S]-covered),S))
        fresh=training.pop(S)-covered
        print('training selection',j+1,S,len(fresh),P[S])
        bank.append(S);covered|=fresh
    return bank

def primes_to(N):
    a=bytearray([1])*(N+1);a[:2]=b'\0\0'
    for i in range(2,isqrt(N)+1):
        if a[i]:a[i*i::i]=b'\0'*((N-i*i)//i+1)
    return a

def score(B,S):
    A=primes_to(S[-1][-1]);caps=[]
    for j,(lo,hi) in enumerate(S):
        U=union(BASE,lo,hi);V=union(B,lo,hi)
        added=V-U
        C=U|V;P=[q for q in range(lo|1,hi,2) if A[q]]
        newp=sum(A[q] for q in added)
        G=[sum(1 for q in P if q in X) for X in (U,C)]
        print(lo,hi,len(U),len(C),len(added),*G,len(P),newp,'first-new-primes',
            sorted(q for q in added if A[q])[:5])
        assert (len(U),len(C),*G,len(P))==EXPECTED[j]
        upper0=union(UPPER_BASE,7*lo,7*hi);upper1=upper0|union(B,7*lo,7*hi)
        pair0=sum(q in U and 7*q in upper0 for q in P)
        pair1=sum(q in C and 7*q in upper1 for q in P)
        print('both axes',pair0,pair1,'/',len(P))
        mass=[0]
        for src in BASE+B:image(src,lo,hi,mass=mass)
        print('nonempty terminal-word mass',mass[0], 'distinct',len(C))
        caps.append((lo,hi,added))
    return caps

def cap_upper(t,d,n=9):
    # ln(t/(t-d)) = 2 atanh(d/(2t-d)); geometric upper bound on tail.
    x=F(d,2*t-d)
    z=2*sum(x**(2*j+1)/F(2*j+1) for j in range(n))
    z+=2*x**(2*n+1)/(F(2*n+1)*(1-x*x))
    return 2*z/d

def main():
    if len(sys.argv)>1:
        cat=Path(sys.argv[1]);assert sha256(cat.read_bytes()).hexdigest()==CAT_SHA
        P=profiles(cat)
        assert len(P)==2806
        print('pool',len(P),'train',TRAIN,'holdouts',SHELLS)
        assert tuple(pick_pool(P))==FROZEN
    for (A,r),S in zip(WITNESSES,FROZEN):
        assert state(A,r)==S
    assert sum(cap_upper(S[2],S[3]) for S in BASE+FROZEN[:1]) < F(983,1000)
    T=(6,28,24,17,1);r=4
    assert state(T,r)==(11,4,12,5)
    for N,w in [(4111,'001111101'),(28777,'100100100011')]:
        assert H(forward(T,r,tuple(map(int,w))))==N
    assert primes_to(4111)[4111]
    print('literal 4111 and 7 * 4111')
    score(FROZEN,SHELLS)

if __name__=='__main__':main()
