from itertools import combinations
from math import gcd

def orient(n,t=0,spine=False):
    E=[(i,j) for i,j in combinations(range(n),2) if not spine or j>i+1]
    A=[0]*n
    if spine:
        for i in range(n-1): A[i]|=1<<(i+1)
    for k,(i,j) in enumerate(E):
        a,b=(j,i) if t>>k&1 else (i,j)
        A[a]|=1<<b
    return tuple(A)

def counts(A):
    n=len(A);W=(1<<n)-1
    dp=[[0]*n for _ in range(W+1)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,W):
        for i,v in enumerate(dp[s]):
            if not v:continue
            e=A[i]&~s
            while e:
                p=e&-e;e^=p;dp[s|p][p.bit_length()-1]+=v
    return [1]+[sum(row) for row in dp[1:]]

def profile(A,maxk):
    h=counts(A);z=[1]+[0]*(len(h)-1);out=[z[-1]]
    while maxk:
        q=[]
        for s in range(len(h)):
            t=s;u=0
            while True:
                u+=h[t]*z[s^t]
                if not t:break
                t=(t-1)&s
            q.append(u)
        z=q;out.append(z[-1]);maxk-=1
    return out

def join(A,B):
    a,b=len(A),len(B)
    return tuple(x|(((1<<b)-1)<<a) for x in A)+tuple(x<<a for x in B)

def triangle(A,B,C):
    a,b,c=len(A),len(B),len(C)
    x=(1<<a)-1;y=((1<<b)-1)<<a;z=((1<<c)-1)<<(a+b)
    return tuple(v|y for v in A)+tuple((v<<a)|z for v in B)+tuple((v<<(a+b))|x for v in C)

def invert(N,H,B):
    if N<=2*H:return []
    D=N-2*H;m=1;out=[]
    while (d:=(1<<m)-1)<=D:
        if D%d==0 and ((D//d)+H)%B==0:
            h=((D//d)+H)//B
            if h>0 and not (h&(h-1)):out.append((m,h.bit_length()-1))
        m+=1
    return out

small=[orient(n,j) for n in range(1,4) for j in range(1<<(n*(n-1)//2))]
V={A:profile(A,5) for A in small}
tests=0
for A in small:
    for B in small:
        a,b=V[A],V[B]
        x=profile(triangle(A,B,(0,)),4)
        for k in range(1,5):
            want=k*(a[k]*b[k]+(a[k+1]-a[k])*(b[k+1]-b[k])+(a[k]-a[k-1])*(b[k]-b[k-1]))
            assert x[k]==want
            tests+=1
print(tests,'literal profile recurrence checks')
F=[
 ('C5',orient(5,22,True),15,120),
 ('Petersen chord',orient(10,4304406596,True),5213,104408),
 ('Petersen chord first descent',orient(10,3041897544,True),1747,50656),
 ('doubled C5',(254,380,696,848,896,24,52,106,229,483),6473,134288),
 ('compressed Petersen',(966,940,888,753,483,193,386,772,552,112),8265,155040)
]
for name,A,H,B in F:
    assert profile(A,2)==[0,H,B]
    print(name,H,B,gcd(H,B))
    for m,t in [(1,0),(1,1),(1,2),(1,3),(2,0),(2,1),(3,0)]:
        N=2*H+((1<<m)-1)*((B<<t)-H)
        assert (m,t) in invert(N,H,B)
        literal=triangle(join(A,orient(t)),orient(m),(0,))
        assert counts(literal)[-1]==N
assert not invert(109623,5213,104408)
