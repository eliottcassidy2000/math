from fractions import Fraction
from math import gcd
from itertools import combinations
from sys import argv

NONTRIVIAL=((5,32),(11,100),(17,224),(19,256),(23,272),(23,296),(25,304),(29,320),(31,328),(37,352),(37,376),(41,416),(43,424))
PROFILES=((1,2),)+NONTRIVIAL
L=144
SCALE=10**9
A=(1000,1362,1374,1227)
B=(116,194,195,186)
TYPES=((0,1),(1,0),(1,1),(1,2))
def ceil(n,d):return (n+d-1)//d
def z(t):
    D,U,V=t
    return D+2*U,V,D+U
def project(x,y):
    x%=3;y%=3
    if x==0:return 0 if y else None
    if y==0:return 1
    return 2 if x==y else 3

def pressure():
    X,Y=(1,0,0),(0,0,1)
    blocks=[]
    for j in range(1,49):
        X,Y=z(X),z(Y)
        a,b,c,d=sum(X),sum(Y),X[1],Y[1]
        for h,q in PROFILES:
            for k in range(14):
                Q=q<<k
                nxt=[]
                for x,y in TYPES:
                    F=a*x+b*y;U=c*x+d*y
                    nxt.append(project(Q*(F+U),h*F))
                blocks.append((a,b,c,d,h,Q,nxt,(h,q)!=(1,2)))
    worst=Fraction(0);worst_unweighted=0
    for i in range(48):
        sums=[0]*4;unweighted=0
        for a,b,c,d,h,Q,nxt,nontriv in blocks:
            f0=a*(L-i)+b*i;f1=f0+b-a
            u0=c*(L-i)+d*i;u1=u0+d-c
            g0=h*f0+Q*(f0+u0);g1=h*f1+Q*(f1+u1)
            minimum=min(g0,g1)
            assert max(g0,g1)<L*2**80
            f,g=(f0,g0) if f0*g1>=f1*g0 else (f1,g1)
            for t,w in enumerate(nxt):
                if w is not None:
                    sums[t]+=ceil(SCALE*L*(A[w]*g+B[w]*h*f),1000*g*minimum)
            if nontriv:unweighted+=ceil(SCALE*L,minimum)
        for t,s in enumerate(sums):
            assert 10**6*L*s <966*SCALE*(A[t]*L+B[t]*i),(i,t,s)
            worst=max(worst,Fraction(1000*L*s,SCALE*(A[t]*L+B[t]*i)))
        worst_unweighted=max(worst_unweighted,unweighted)
    assert worst_unweighted==157998488
    M=tuple((i==j)*1 for i in range(3) for j in range(3))
    for k in range(3):
        t=M[3*k:3*k+3];assert z(z(z(t)))[k]>=2
    mass=sum(Fraction(1,q) for h,q in NONTRIVIAL)
    assert mass<Fraction(77,1000) and mass+Fraction(1,2)<Fraction(577,1000)
    tail=Fraction(10,65536)+Fraction(5,8192)
    assert Fraction(144,100)*(mass+Fraction(1,2))*tail<Fraction(1,1000)
    assert Fraction(worst_unweighted,SCALE)+mass*tail<Fraction(159,1000)
    beta=sum(ceil(3*SCALE,h+2*q) for h,q in NONTRIVIAL)
    assert beta==108211408 and beta<Fraction(109,1000)*SCALE
    x=Fraction(109,1000)
    total=Fraction(967,1000)+Fraction(144,100)*Fraction(159,1000)*x/(1-x)
    assert total<Fraction(996,1000)
    assert all(1000<=A[t] and 3*A[t]+B[t]<=4320 for t in range(4))
    # Explicit exponent s=65535/65536; all finite terms inflate by
    # at most f because their growth is < 2^80. Extend the tails
    # with the larger exact common ratio r for 2^{-s}.
    f=Fraction(1001,1000)
    r=f/2
    assert 1001**65536 > 1000**65536*2**80
    jt=2*r**16 + 3*r**17/(1-r)
    kt=r**14/(1-r)
    tail_s=jt/(1-r)+(2+3*r/(1-r))*kt
    single_s=f*(Fraction(966,1000)+Fraction(144,100)*Fraction(577,1000)*tail_s)
    nontriv_s=f*(Fraction(worst_unweighted,SCALE)+Fraction(77,1000)*tail_s)
    beta_s=f*Fraction(beta,SCALE)
    total_s=single_s+Fraction(144,100)*nontriv_s*beta_s/(1-beta_s)
    assert total_s<Fraction(996,1000)
    print(f'Explicit s=65535/65536 pressure upper {float(total_s):.10f}.')
    print(f'192 exact comparisons pass; actual rectangle weighted max {float(worst):.10f}; unweighted numerator {worst_unweighted}.')
    print(f'Exact total upper {total} = {float(total):.10f}; beta ceiling numerator {beta}.')

def spectrum():
    sets={}
    for n in range(3,7):
        edges=list(combinations(range(n),2)); found=set()
        N=(1<<n)-1
        for code in range(1<<len(edges)):
            rows=[0]*n
            for k,(i,j) in enumerate(edges):
                if code>>k&1:rows[i]|=1<<j
                else:rows[j]|=1<<i
            s=1
            while True:
                x=s
                for i in range(n):
                    if s>>i&1:x|=rows[i]
                if x==s:break
                s=x
            if s!=N:continue
            s=1
            while True:
                x=s
                for i in range(n):
                    if rows[i]&s:x|=1<<i
                if x==s:break
                s=x
            if s!=N:continue
            dp=[[0]*n for _ in range(N+1)]
            for i in range(n):dp[1<<i][i]=1
            for S in range(1,N):
                for i,m in enumerate(dp[S]):
                    if not m:continue
                    z=rows[i]&~S
                    while z:
                        bit=z&-z;z-=bit
                        dp[S|bit][bit.bit_length()-1]+=m
            h=sum(dp[N]);hs=[1]+[sum(d)for d in dp[1:]]
            Q=sum(x*y for x,y in zip(hs,reversed(hs)))
            if gcd(h,Q)==1:found.add((h,Q))
        sets[n]=found
    assert set().union(*sets.values())==set(NONTRIVIAL),sets
    print("Exhaustive coprime strong profiles by order:",sorted((n,sorted(v)) for n,v in sets.items()))

if __name__=='__main__':
    pressure()
    if '--profiles' in argv:spectrum()
