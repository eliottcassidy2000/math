"""Conditional simple-cycle census from the diameter-two Moore identity."""
from fractions import Fraction as F
from math import comb, prod, factorial
from collections import Counter
from itertools import combinations
import csv
from decimal import Decimal
import sys

SPECS={2:None,3:(1,-2,5,4),7:(2,-3,28,21),57:(7,-8,1729,1520)}

def trace(k,d):
    r,s,a,b=SPECS[d]; q=d-1; n=d*d+1
    def powers(l):
        if k==0:return 2
        x,y=2,l
        for _ in range(1,k): x,y=y,l*y-q*x
        return y
    return 1+q**k+a*powers(r)+b*powers(s)+n*(d-2)//2*(1+(-1)**k)

def route_channels(d,K):
    q=d-1
    # Matrices aI+bA+cJ; multiply by A via its quadratic relation.
    def times(v):
        a,b,c=v
        return q*b,a-b,b+d*c
    prev=(0,0,0); now=(1,0,0)
    diag=[0]*(K+1); off=diag.copy()
    old=old2=(0,0,0)
    for j in range(K):
        def entry(t,dist):
            return t[2]+(t[dist] if dist<2 else 0)
        if j>=1:
            for a,dist in [(diag,0),(off,2)]:
                l=j+1
                a[l]=a[l-2]+entry(old,dist)-2*entry(old2,1)-(j==1 and dist==0)
                if j>=3:
                    a[l]+=entry(older,0)
        nxt=tuple(a-q*b for a,b in zip(times(now),prev))
        older,old2,old=old2,old,now
        prev,now=now,nxt
    assert not any(diag[1:5]) and not any(off[1:5])
    return diag,off

def coll(d,K,H):
    n=d*d+1; q=d-1
    D,O=H
    L=[d*(a+q*b) for a,b in zip(D,O)]
    def cross(a,b):
        return L[a]*L[b]*(d-2)//d+d*(D[a]*D[b]+q*O[a]*O[b])
    return sum(n*cross(a,K-a) for a in range(5,K-4))*K//2

def exact(k,d):
    if d==2:return int(k==5)
    if k<5 or k>d*d+1:return 0
    tr=trace(k,d); n=d*d+1; q=d-1
    if k<=9:return tr//(2*k)
    P=n*d*q*q//10
    if k==10:
        s3=n*d*q//2*comb(q,2)
        s2=n*d//2*comb(q*q,2)-2*s3
        s1=n*comb(d*q*q//2,2)-2*s2-3*s3
        return (tr-10*P-40*s1-20*(s2+s3))//20
    if k==11:
        s4=n*d*q*q//2*(q-1)
        s3=n*d*q//2*q*q*(q-1)-2*s4
        s2=n*d//2*q**4*(q-1)-2*s3-3*s4
        iso=n*d*q//2*q*comb(q,2)-s4
        s1=n*(d*q*q//2)**2*(q-1)-2*s2-3*s3-4*s4-2*iso
        counts=(s1,s2,s3,s4,iso)
        assert min(counts)>=0, counts
        return (tr-44*(s1+iso)-22*(s2+s3+s4))//22
    return None

def branch(k,d):
    n=d*d+1;q=d-1;best=F(0)
    for m in range(3,d+1):
        if not (m<k<=3*m and (k-m)%2==0): continue
        j=(k-m)//2
        avg=F(q**j-1) if j%2==0 else F(q**j+1)-F(q,d-m+1)
        if j==1: r=n
        else:
            r=min(n,comb(k,2)//comb(j,2))
            if j*j>k:r=min(r,k*(j-1)//(j*j-k))
        order=prod(range(d-m+1,d+1))//(2*m)
        best=max(best,F(n*order*comb(m,j),r)*avg)
    return -(-best.numerator//best.denominator)

def ceildiv(a,b):return -(-a//b)
def bounds(k,d,H):
    v=exact(k,d)
    if v is not None:return v,v
    tr=trace(k,d); c=coll(d,k,H)
    # Every bad word contributes between one and M equal-position pairs.
    p=k//5; M=(k//p)*comb(p,2)+comb(k%p,2)
    return max(0,ceildiv(tr-c,2*k),branch(k,d)),(tr-ceildiv(c,M))//(2*k)

def graph(d):
    if d==3:
        V=list(combinations(range(5),2))
        return [set(j for j,v in enumerate(V) if not set(u)&set(v)) for u in V]
    V=list(range(50)); A=[set() for _ in V]
    def add(x,y):A[x].add(y);A[y].add(x)
    for i in range(5):
        for x in range(5):
            add(5*i+x,5*i+(x+1)%5)
            add(25+5*i+x,25+5*i+(x+2)%5)
            for j in range(5):add(5*i+x,25+5*j+(x+i*j)%5)
    return A

def controls():
    for d in (3,7):
        H=route_channels(d,13)
        A=graph(d); n=len(A); root=0
        for j in range(n):
            assert len(A[j])==d
            for k in range(j):
                assert len(A[j]&A[k])==int(j not in A[k])
        # Explicit oriented-edge powers, independent of the resolvent.
        x,y=sorted(A[root])[:2]
        step={(root,x):1}
        for l in range(1,14):
            assert (step.get((x,root),0),step.get((y,root),0))==(H[0][l],H[1][l])
            nxt=Counter()
            for (a,b),w in step.items():
                for c in A[b]-{a}:nxt[b,c]+=w
            step=nxt
        cnt=[0]*12
        # All length <=11 on a branching-seven graph can take a while;
        # verify via one root and use transitivity in the known controls.
        def rooted(a,p,seen,k):
            for b in A[p]:
                if b==a:
                    if k>=5:cnt[k]+=1
                elif k<min(11,n) and not (seen>>b&1):
                    rooted(a,b,seen|(1<<b),k+1)
        rooted(0,0,1,1)
        assert all(cnt[k]*n==exact(k,d)*2*k for k in range(5,12)),(d,cnt)
        print('Control',d,[exact(k,d) for k in range(5,12)],flush=True)

def factorial_envelope(d):
    Q=10**6
    U={}
    for j in range(1,d+1):
        target=factorial(j)*Q**j
        lo,hi=0,d*Q
        while hi-lo>1:
            m=(hi+lo)//2
            if m**j>=target: hi=m
            else:lo=m
        U[j]=hi
    hull=[]
    def slope(a,b):
        return F(U[b]-U[a],b-a)
    for j in U:
        while len(hull)>=2 and slope(hull[-2],hull[-1])<=slope(hull[-1],j):
            hull.pop()
        hull.append(j)
    af=[(F(U[a])-slope(a,b)*a,slope(a,b))
        for a,b in zip(hull,hull[1:])]
    assert all(s>=0 and all(i+s*j>=U[j] for j in U) for i,s in af)
    return lambda x:min((i+s*x)/Q for i,s in af)

def permanent_bound(k, n, majorant):
    avg=F(2*((7*65*k+k*k)//130),k)
    v=majorant(avg)
    return comb(n,k)*v.numerator**k//(2*v.denominator**k)

def trace_table(d, K):
    r,s,a,b=SPECS[d];q=d-1;n=d*d+1
    rr,ss=2,2;R,S=r,s
    power=1
    out=[0]
    for k in range(1,K+1):
        power*=q
        out.append(1+power+a*R+b*S+n*(d-2)//2*(1+(-1)**k))
        rr,R=R,r*R-q*rr
        ss,S=S,s*S-q*ss
    return out

def table(outfile=None):
    d=57; H=route_channels(d,171)
    n=d*d+1; tr=trace_table(d,n); majorant=factorial_envelope(d)
    schedule={5,6,7,8,9,10,11,12,21,57,63,81,90,100,140,170,171,200,500,1000,1420,1430,1431,2000,3250}
    old=0; done=[]
    for k in range(5,3251):
        lo,hi=bounds(k,d,H) if k<=171 else (0,tr[k]//(2*k))
        if k>=1000:
            hi=min(hi,permanent_bound(k,n,majorant))
        if 5<k<=1432:
            j=k-1
            # 2*floor((7+50*j/3250)*j/2) bounds the sum of internal degrees.
            induced_degree_sum=2*((7*65*j+j*j)//130)
            lifts=max(0,58*j-2*induced_degree_sum)
            lo=max(lo,ceildiv(old*lifts,k))
        assert 0<=lo<=hi and (lo>0 or k>1430),(k,lo,hi)
        old=lo
        if outfile:
            done.append([k,lo,hi])
        if k in schedule:
            print(k, lo if k<=11 else (f'{Decimal(lo):.8E}' if lo else 0), hi if k<=11 else f'{Decimal(hi):.8E}',flush=True)
    if outfile:
        sys.set_int_max_str_digits(0)
        with open(outfile,'w',newline='') as f:
            w=csv.writer(f,lineterminator='\n');w.writerow(['length','lower','upper']);w.writerows(done)
if __name__=='__main__':
    if '--controls-only' not in sys.argv:
        table(next((s for s in sys.argv[1:] if not s.startswith('--')),None))
    if '--bounds-only' not in sys.argv: controls()
