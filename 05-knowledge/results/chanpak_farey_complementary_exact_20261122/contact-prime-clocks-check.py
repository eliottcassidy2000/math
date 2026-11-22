from itertools import permutations
from collections import defaultdict

A=1<<60
F,G=35183298379775,10504481267175
CO={
53:(36,158450767610363602280225177600,346084103871386661670010152033395225),
56:(55,41538374868278620992215173614796800,-344221346644057575),
}

def rows(s,t):
    n=s+t+3;B=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    for i,j in ((0,s+1),(s+1,n-1)):B[i]^=1<<j;B[j]^=1<<i
    return B

def paths(B,O,end):
    n=len(B);Ans=0
    V=[v for v in range(n) if v not in O]
    for o in permutations(O):
        S={(0,0,0,0):1}
        for v in V:
            T=defaultdict(int)
            for d,count in S.items():
                for i in range(4):
                    if not d[i] and i and not (B[o[i-1]]>>v&1):continue
                    leave=((B[v]>>o[i]&1) if i<3 else (end>>v&1))
                    e=d[:i]+(2-leave,)+d[i+1:]
                    T[e]+=count
            S=T
        for d,z in S.items():
            if 2 in d:continue
            if any(not d[i] and not (B[o[i-1]]>>o[i]&1) for i in (1,2)):continue
            if d[3] or end>>o[-1]&1:Ans+=z
    return Ans

def reads(s,t,a=None):
    B=rows(s,t);n=len(B);f=(1<<n)-1;O=(0,s+1,n-1)
    if a is not None:B[0]^=1<<a;B[a]^=1
    k=paths(B,O,f)
    end=(f^1)&~B[0]
    return k,paths([f^1]+[v|1 for v in B[1:]],O,end)

def phases(i,p,r):
    s,C,D=CO[i];m=p**r;d=D*pow(2,-s,m)%m
    X,P=A%m,1;S=set();R=[]
    while (X,P) not in S:
        if (C+d*X)*P%m==2%m:R.append(len(S))
        S.add((X,P));X,P=X*X%m,P*(X+1)%m
    assert (X,P)==(A%m,1)
    return len(S),R

def powers(a,m):
    p=1;R=[]
    while p not in R:
        R.append(p);p=p*a%m
    assert p==1
    return R

def mod_delta(i,p,r,N):
    s,C,D=CO[i];m=p**r;den=A-1;t=1
    while den%p==0:den//=p;t*=p
    x=pow(A,N,m*t)
    return ((x-1)//t)*pow(den,-1,m)*(C+D*pow(2,-s,m)*(x%m))%m

def lifts(i,p,d,j0,up):
    k=0;step=1;R=[]
    for r in range(2,up+1):
        M=p**r
        vals=[x for a in range(p) if mod_delta(i,p,r,pow(2,j0+d*(x:=k+a*step),M))==2]
        assert len(vals)==1
        k=vals[0];step*=p
        R.append((j0+d*k,d*step))
    return R

for i in (53,56):
    s0,C,D=CO[i]
    for w in range(1,5):
        s,t=(15,37+w) if i==53 else (56+w,15)
        k,z=reads(s,t);u,v=reads(s,t,i)
        got=(u-k,v-z)
        exp=(4503530907893760+9836602018824134393855*(1<<w),
             4503599627370496+9836602018824134393856*(1<<w)) if i==53 else (
             1180627649514430267392-32769*(1<<w),1180627649514430267392)
        assert got==exp,(i,w,got,exp)
        assert G*got[0]+(F-G)*got[1]==C+D*(1<<w)
assert phases(53,3,1)==(2,[0]) and phases(53,3,2)==(6,[2])
for i,p,pr in ((53,13,(12,[0])),(53,19,(2,[0])),(53,47,(11,[0])),(56,13,(12,[0])),
               (56,59,(28,[0,3])),(56,71,(3,[0])),(53,11,(10,[1])),(56,11,(10,[]))):
    assert phases(i,p,1)==pr
assert lifts(53,13,132,0,3)==[(528,1716),(3960,22308)]
assert lifts(56,13,84,0,3)==[(1008,1092),(4284,14196)]
assert lifts(53,11,10,1,3)==[(1,110),(881,1210)]
for i,p,m,E in ((53,19,3,[34]),(53,47,23,[868,1013]),(56,59,29,[182,929]),(56,71,7,[183,444])):
    s,C,D=CO[i];M=p*p;T=D*pow(2,-s,M)%M
    R=[x for x in range(M) if ((x-1)*(C+T*x)-2*(A-1))%M==0]
    V=powers(A,M)
    assert R and [k for k,x in enumerate(V) if x in R]==E and len(V)==m*p
    assert all(n not in powers(2,m*p) for n in E)
print('contact affine laws, individual clocks, forbidden squares, and example Hensel digits pass')
