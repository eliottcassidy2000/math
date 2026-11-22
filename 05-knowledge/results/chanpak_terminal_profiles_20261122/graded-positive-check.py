"""Literal directed subset recurrences for the graded crossing bridge."""
from itertools import product

def paths(A, tags=None):
    m=1<<len(A);n=len(A)
    dp=[[0]*n for _ in range(m)]
    for u in range(n):
        if tags is None or tags>>u&1:dp[1<<u][u]=1
    for s in range(1,m):
        for u,v in enumerate(dp[s]):
            bits=A[u]&~s
            while bits and v:
                q=bits&-bits;dp[s|q][q.bit_length()-1]+=v;bits^=q
    return dp

def H(A):
    return sum(paths(A)[-1]) if A else 1

def grade(A,tags=0):
    seq=list(map(sum,paths(A)));seq[0]=1
    full=(1<<len(A))-1
    val=seq[full]
    covers=sum(seq[i]*seq[full^i] for i in range(full+1))
    assert sum(seq[i]*seq[full^i] for i in range(full+1) if i.bit_count()%2==0)==covers//2
    tagged=paths(A,tags)
    mass=sum(sum(v)*seq[full^i] for i,v in enumerate(tagged))
    E=sum(sum(q for j,q in enumerate(v) if tags>>j&1)*seq[full^i]
          for i,v in enumerate(tagged) if i.bit_count()%2==0)
    return val,covers,mass,E

def jo(*B):
    n=sum(map(len,B));p=0;ans=[]
    for A in B:
        e=p+len(A);tail=(1<<n)-(1<<e)
        ans.extend(tail|r<<p for r in A);p=e
    return tuple(ans)

def turn(A,u,v):
    A=list(A);A[u]^=1<<v;A[v]^=1<<u;return tuple(A)

def bridge(A,X,S):
    T=jo(A,S,(0,))
    for u in range(len(A)):
        if X>>u&1:T=turn(T,u,len(T)-1)
    return T

def strong(A):
    full=(1<<len(A))-1
    for i in range(len(A)):
        s=1<<i
        while True:
            t=s
            for j,r in enumerate(A):
                if s>>j&1:t|=r
            if s==t:break
            s=t
        if s!=full:return False
    return True

def trans(n):
    return tuple((1<<n)-(2<<i) for i in range(n))

def collar(D):
    return turn(jo((0,),D,(0,)),0,len(D)+1)

B=((0,), (2,0), (2,4,1), (2,4,9,3))
SEEDS=(((0,),1,1),((2,4,1),3,6),((2,4,9,3),3,11),
       ((2,4,9,3),6,12),((18,4,9,19,6),18,29))
for A,X,C in SEEDS:
    h,Q,c,e=grade(A,X)
    assert (c,h-2*e)==(C,1)
    for D,S in product(B,repeat=2):
        hD,qD,*_=grade(D);hS,qS,*_=grade(S)
        T=bridge(jo(A,D),X,S)
        assert strong(T) and H(T)==hD*hS+C*qD*qS
for D,S in product(B,repeat=2):
    hD,qD,*_=grade(D);hS,qS,*_=grade(S)
    A=collar(D);X=1|(1<<(len(A)-1))
    assert grade(A,X)[::2]==(hD+qD,3*qD)
    T=bridge(A,X,S)
    assert strong(T) and H(T)==hD*hS+3*qD*qS
for a,b in product(range(1,5),repeat=2):
    S=collar(trans(b))
    h,q,*_=grade(S)
    assert (h,q)==(1+2**b,2**(b+2)+2*3**b-2)
    T=bridge(collar(trans(a)),1|(1<<(a+1)),S)
    assert H(T)==1+2**b+3*2**a*q
print("graded crossing bridge, five fixed grades, and symmetric 3-bridge: literal PASS")
P=(2,4,9,3)
R=(38,4,8,19,39,14)
assert grade(P)[:2]==(5,32)
assert grade(R)[:2]==(19,256)
assert strong(P) and strong(R)
for i,j,pad in ((1,0,0),(0,1,0),(1,1,0),(1,0,2),(0,0,4)):
    k=5*i+8*j+pad
    S=jo(*([P]*i+[R]*j+[(0,)]*pad))
    assert grade(S)[:2]==(5**i*19**j,2**k)
    for A,X,L in SEEDS[:3]:
        T=bridge(A,X,S)
        assert H(T)==5**i*19**j+L*2**k and strong(T)
print("slot-neutral 5 and 19 modules: literal PASS")
A=(244,81,210,39,8,150,56,90)
B=turn(A,0,5)
assert grade(A)[:2]==(383,5960) and grade(B)[:2]==(389,5960)
for S,X,L in SEEDS[:3]:
    for k in (0,1):
        X0,Y0=tuple(bridge(S,X,jo(D,trans(k))) for D in (A,B))
        assert strong(X0) and strong(Y0)
        assert Y0==turn(X0,len(S),len(S)+5)
        assert (H(X0),H(Y0))==(383+5960*L*2**k,389+5960*L*2**k)
print("literal Q-neutral +6 edge and grafts PASS")
