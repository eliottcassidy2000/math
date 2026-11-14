from itertools import product, combinations
from math import prod, gcd
from fractions import Fraction

def params(h):
    J=prod(h); W=tuple(J//m for m in h)
    u=tuple(pow(w,-1,m) for m,w in zip(h,W))
    return J,W,u

def digits(h,t):
    J,W,u=params(h)
    return tuple((v*t)%m for v,m in zip(u,h))

def representatives(h):
    J,W,u=params(h)
    F=[sum(w*d for w,d in zip(W,digits(h,t))) for t in range(J)]
    B=(sum(w*v for w,v in zip(W,u))-1)//J
    good=0
    for t in range(J):
        C=sum(d+v>=m for d,v,m in zip(digits(h,t),u,h))
        assert F[(t+1)%J]-F[t]==1+J*(B-C)
        good+=B==C
    p=[Fraction(v,m) for v,m in zip(u,h)]
    prob=sum(prod((a if b else 1-a) for a,b in zip(p,S))
             for S in product((0,1),repeat=len(h)) if sum(S)==B)
    assert good==J*prob
    return F,good

def forward(h,D):
    W=params(h)[1]
    return {sum(w*d for w,d in zip(W,ds)) for ds in product(*D)}

def inverse(h,D,t):
    J,W,u=params(h); r=digits(h,t)
    c=(t-sum(w*a for w,a in zip(W,r)))//J
    Z=[{(d-a)//m for d in pool if d%m==a} for pool,a,m in zip(D,r,h)]
    return any(sum(z)==c for z in product(*Z))

assert representatives((3,5))[1]==8
assert sorted(representatives((3,5))[0])==[0,3,5,6,8,9,10,11,12,13,14,16,17,19,22]
assert set(range(15,30))<=forward((3,5),(range(6),range(5)))
for k in range(2,4):
    for h in combinations(range(3,12,2),k):
        if any(gcd(p,q)>1 for p,q in combinations(h,2)):continue
        J,W,u=params(h)
        F,_=representatives(h)
        assert max(sum(a<=f<a+J for f in F) for a in F)*36<=35*J
        D=tuple(range(2*m) if i<k-1 else range(m) for i,m in enumerate(h))
        S=forward(h,D)
        assert set(range((k-1)*J,k*J))<=S
        for t in range(-5,max(S)+5):
            assert inverse(h,D,t)==(t in S)

def path_arrays(R):
    n=len(R);F=(1<<n)-1
    P=[[0]*n for _ in range(F+1)]
    for j in range(n):P[1<<j][j]=1
    for s in range(1,F+1):
        for j,a in enumerate(P[s]):
            if a:
                out=R[j]&(F^s)
                while out:
                    b=out&-out;out^=b
                    P[s|b][b.bit_length()-1]+=a
    return P
def profile(R,U):
    n=len(R);F=(1<<n)-1
    E=path_arrays(R);S=path_arrays(tuple(F^r^(1<<j) for j,r in enumerate(R)))
    H=[1]+[sum(q) for q in E[1:]]
    ps=[sum(v for j,v in enumerate(q) if U>>j&1) for q in S]
    es=[sum(v for j,v in enumerate(q) if U>>j&1) for q in E]
    qs=[sum(v for j,v in enumerate(q) if not U>>j&1) for q in E]
    h=H[-1]
    f=sum(H[s]*H[F^s] for s in range(F+1))
    P=sum(ps[s]*H[F^s] for s in range(F+1))
    Q=sum(qs[s]*H[F^s] for s in range(F+1))
    T=sum(v for j,v in enumerate(S[-1]) if not U>>j&1)
    e=sum(es[s]*ps[F^s] for s in range(F+1))
    aug=tuple(r|((not U>>j&1)<<n) for j,r in enumerate(R))+(U,)
    K=sum(path_arrays(aug)[-1])
    assert K==h+P-es[-1]-e
    return (h,f,P,Q,es[-1],T),e

cases=[
 (((2,20,9,35,77,23,47),17),((2,20,33,7,45,75,31),33)),
 (((130,172,169,209,135,217,151,0),130),((0,9,147,101,171,199,23,75),9)),
 (((2,12,17,37,139,151,63,79),130),((2,12,17,37,139,151,63,79),36))
]
h=[]
for a,b in cases:
    A,e=profile(*a)
    B,f=profile(*b)
    assert A==B and f-e==2
    h.append(A[0])
assert h==[65,109,161]
S=forward(h,((0,1),)*3)
assert sorted(S)==[0,7085,10465,17549,17550,24634,28014,35099]
assert [v for v in S if v+1 in S]==[17549]
print("PASS: exact CRT wraps, height strips, redundant quotients, and six literal source profiles")
