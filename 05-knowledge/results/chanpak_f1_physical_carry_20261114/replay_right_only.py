from math import isqrt
from replay_f1 import Core, audit, path_arrays, make
from compile_absolute import binary_h,joined

low3=Core((4,1,2,7),7,3,24,9,3,3,3,3)
hi3=Core((10,0,3,6),3,3,24,9,6,3,2,5)
low9=Core((20,1,2,7,14),22,9,84,41,22,7,4,23)
hi9=Core((16,1,3,7,14),22,9,84,41,25,7,5,25)
for C in [low3,hi3,low9,hi9]:audit(C)
assert tuple(a^b for a,b in zip(low9.masks,hi9.masks))==(4,0,1,0,0)

def strong(R):
    n=len(R)
    for s in range(n):
        v=1<<s
        while True:
            u=v
            for i in range(n):
                if v>>i&1:u|=R[i]
            if u==v:break
            v=u
        if v!=(1<<n)-1:return False
    return True

for C,N in [(low9,40241),(hi9,40243)]:
    R=make([C],0,925)
    assert len(R)==16 and strong(R)
    assert path_arrays(R)[0][-1]==N
    assert C.K+C.h*(binary_h(925)-1)+C.P*925+C.S==N
    assert all(binary_h(y)!=N for y in range(1,4*(N-1)))
assert 40243==7*5749
for p in (40241,5749):
    assert all(p%d for d in range(2,isqrt(p)+1))

samples=(Core((0,),0,1,2,0,1,0,1,1),
         Core((0,),1,1,2,1,0,1,0,1),low3,hi3,low9,hi9)
for X in samples:
    for Y in samples:
        for A,B,s in [(low3,hi3,-1),(low9,hi9,1)]:
            diff=joined(joined(X,B),Y).K-joined(joined(X,A),Y).K
            assert diff==Y.h*(2*X.h+3*X.P+s*X.S)
            assert (diff==2)==(X.h==Y.h==1 and X.P==0)
print("PASS: four profiles, literal 16-vertex counts, all-order transitive-spine exclusions, two-sided limits")

for C in samples:
    f=lambda x,y:C.K+C.h*(binary_h(x)+binary_h(y)-2)+C.Q*x+C.P*y+C.f*x*y+C.T*(x%2)+C.S*(y%2)
    K,P1,P3,P2,Q1,Q3,B=f(0,0),f(0,1),f(0,3),f(0,2),f(1,0),f(3,0),f(1,1)
    P=(P3-P1)//2;Q=(Q3-Q1)//2
    S=P1-K-P;T=Q1-K-Q
    h=(P2-K-2*P)//2
    fc=B-K-P-Q-S-T
    assert (K,h,fc,P,Q,S,T)==(C.K,C.h,C.f,C.P,C.Q,C.S,C.T)
print("PASS: four and seven probe recovery")
