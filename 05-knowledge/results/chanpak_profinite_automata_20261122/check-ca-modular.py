"""Replay exact DFA state counts and the suffix difference certificate."""
from collections import deque
from math import gcd, prod

def f(y):
    if y == 0:
        return 0
    x = y//2
    return f(x) + (not y%2)*(x+(x%2))

def L(y):
    return 149+162*y+25*f(y)

def check(M):
    d=M//gcd(M,50)
    def step(s,b):
        x,e,z=s
        return ((2*x+b)%M, b if d>1 else 0,
                (z+(1-b)*(x+e)*pow(2,-1,d))%d)
    def out(s):
        x,e,z=s
        return (149+162*x+50*z)%M
    def o(s,w):
        for b in w:s=step(s,int(b))
        return out(s)
    Q={(0,0,0)}
    q=deque(Q)
    while q:
        s=q.popleft()
        for b in (0,1):
            t=step(s,b)
            if t not in Q:
                Q.add(t)
                q.append(t)
    assert len(Q)==(M if d==1 else 2*M*d)
    assert len({tuple(o(s,w) for w in ('','1','0','00')) for s in Q})==len(Q)
    for y in range(230):
        assert o((0,0,0),bin(y)[2:])==L(y)%M

for m in (1,3,5,7,9,15,25,27,35,49,75):
    check(m)

for D in range(1,80):
    b=(2*D).bit_length()
    B=1<<b
    u,v=B//2-D,B//2
    rev=lambda n:int(f'{n:0{b}b}'[::-1],2)
    d=rev(u)-rev(v)
    vals=[L(B*z+v)-L(B*z+u) for z in range(10)]
    assert vals==[vals[0]+25*(d*z-z%2) for z in range(10)]
    assert [vals[i+1]-vals[i] for i in range(2)]==[25*(d-1),25*(d+1)]

m=prod((3,5,7,11,13,17,19,23,29,31,41,47,59,71))
print("modular reader replay PASS; odd Ogg product:",m,"DFA states:",2*m*m//5)
