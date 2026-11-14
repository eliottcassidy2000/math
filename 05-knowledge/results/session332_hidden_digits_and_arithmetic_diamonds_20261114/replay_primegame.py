"""Exact PRIMEGAME entrance, checkpoints, timings, and literal block inventory."""
from itertools import product
from math import gcd, isqrt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from profile_tools import profile, strong

P=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),
   (77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
Alt=P[:11]+((15,14),(15,2),(55,1))

def fact(x,p):
    k=0
    while not x%p:
        x//=p;k+=1
    return k

def next_step(x,P=P):
    for a,b in P:
        if x%b==0:
            return x//b*a
    raise RuntimeError('original end instruction lost')

def get_trial(x):
    if x%13 or x%169==0:
        return None
    n=fact(x,5)
    d=fact(x,7)
    return (n,d) if x==13*5**n*7**d and n else None

def trial(n,d,P=P):
    x=5**n*7**d*13
    t=0;p=[]
    while True:
        x=next_step(x,P)
        t+=1
        if x&(x-1)==0:
            p.append(x.bit_length()-1)
        end=get_trial(x)
        if end is not None:
            return end,t,p

def enter(seed,P=P):
    for _ in range(10000):
        n=get_trial(seed)
        if n is not None:
            return n
        seed=next_step(seed,P)
    raise RuntimeError('initialisation not recovered')

def spf(n):
    return next((d for d in range(2,isqrt(n)+1) if not n%d),n)

for u,v,w,t in product(range(5),repeat=4):
    s=2**u*3**v*5**w*7**t
    assert enter(s)==(u+w+1,u+v)
    assert enter(s,Alt)==(u+w+1,v+max(u,t))
for n,d in product(range(1,23),range(29)):
    end,t,p=trial(n,d)
    alt,ta,pa=trial(n,d,Alt)
    if d==0:
        assert end==(n,0) and t==2 and p==[]
        continue
    q,r=divmod(n,d)
    if r:
        assert end==(n,d-1) and t==6*n+2*q+2 and p==[]
        assert ta==t
    else:
        assert end==(n+1,n) and t==7*n+2*q+d+3
        assert p==([n] if d==1 else [])
        assert ta==t-(d-1)
    assert alt==end and p==pa

def epoch(n,d):
    work=0;p=[]
    while True:
        end,t,hit=trial(n,d)
        work+=t;p+=hit
        if end[0]>n or d==0:
            return work,p,end
        _,d=end

for n in range(2,36):
    g=n//spf(n)
    prescribed=sum(6*n+2*(n//d)+2 for d in range(g+1,n))+7*n+2*n//g+g+3
    time,p,end=epoch(n,n-1)
    assert time==prescribed and p==([n] if g==1 else []) and end==(n+1,n)
assert epoch(10,9)==(338,[],(11,10))
for s in [7,8,9,10,12,14,15,16,18,25,50]:
    n,d=enter(s)
    time,p,end=epoch(n,d)
    print('seed',s,'entry',(n,d),'first pure registers',p,'next',end)
for s in range(1,51):
    u,v,w,t=(fact(s,j) for j in (2,3,5,7))
    if s!=2**u*3**v*5**w*7**t:
        continue
    n,d=enter(s)
    _,power,_=epoch(n,d)
    assert (bool(power) and n>1 and spf(n)<n)==(s==50)
print('canonical candidate 10 steps',epoch(10,9)[0])
B={
    2:((0,),1,2),
    3:((2,4,1),3,12),5:((2,4,9,3),5,32),
    7:((2,4,25,35,75,23,47),49,824),
    11:((2,20,9,3,13),11,100),13:((2,12,9,17,7),13,104),
    17:((2,4,9,19,39,15),17,224),19:((2,4,41,19,7,27),19,256),
    23:((2,4,25,19,35,15),23,272),29:((2,12,17,21,35,15),29,320),
}
for p,(U,h,q) in B.items():
    assert profile(U)==(h,q) and strong(U)
def join(A,B):
    return tuple(x|(((1<<len(B))-1)<<len(A)) for x in A)+tuple(x<<len(A) for x in B)
for s in (2,6,15,7,10,50,4):
    graph=(); h=q=1
    for p,(U,a,b) in B.items():
        for _ in range(fact(s,p)):
            graph=join(graph,U)
            h*=a;q*=b
    assert profile(graph)==(h,q)
print('PRIMEGAME return, timings, and tournament register inventory: PASS')
