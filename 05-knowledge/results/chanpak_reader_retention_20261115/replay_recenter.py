"""Literal painted-core profiles, exact affine obstruction and dyadic recentering certificate."""
from replay_reader import count
from math import gcd
M,Q,h=16384,5760,45
B0,T0,S0,b0=2855,4034,3914,11158

def require(x):
    if not x:
        raise AssertionError("certificate mismatch")

def p(n):
    out=0
    while n:
        n,b=divmod(n,2)
        if not b:out+=n+(n%2)
    return out

def revs(q):
    out=[0]*q
    for i in range(1,q):
        out[i]=(out[i>>1]>>1)+(i&1)*(q>>1)
    return out

D=[]
for m in range(256):
    B=count(m,0)
    D.append((m,B,count(m,-1)-B,count(m,1)-B))
require(D[35]==(35,B0,T0,S0))
require(all(2131<=B<=4671 and 2880<=T<=5670 and 2880<=S<=5670 for m,B,T,S in D))
R=revs(M)
P=[p(b) for b in range(M)]
best=(10**12,)
exact=[]
for m,B,T,S in D:
    if S==S0:continue
    for b in range(M):
        a=T-T0+M*(S-S0)+Q*(b-b0)+h*(R[b0]-R[b])
        c=B-B0+S*b-S0*b0+h*(P[b]-P[b0])
        best=min(best,(abs(a),a,m,b,c))
        for e in (0,1):
            ce=c+h*e*(b<M//2)
            if not a:raise AssertionError("unbounded repaint")
            if ce%a==0 and -ce//a>=0 and (-ce//a)%2==e:
                exact.append((-ce//a,m,b))
require(best==(14,14,177,9347,-1184510))
require(exact==[(0,233,8886)])
print("Fixed-row block: sharp slope 14, unique equality on row zero; 2-step near-return on row 84608.")

motifs=[]
for k in range(15):
    d=1<<k;e=M//d
    rd,re=revs(d),revs(e)
    for m,B,T,S in D:
        if m==35:continue
        t=T*d+S*e-T0-S0*M
        if t%h:continue
        t=t//h +d+e-M-1+R[b0]
        for c in range(-4,5):
            for j in (-1,0,1):
                for a in range(d):
                    g=t-rd[a]+128*(M*c+e*a-b0)+j
                    # Need 0 <= g+128*d*b = reverse_e(b) <= e-1.
                    lo=max(0,-(g//(128*d)))
                    hi=min(e-1,(e-1-g)//(128*d))
                    for b in range(lo,hi+1):
                        if g+128*d*b==re[b]:
                            motifs.append((c,j,d,e,m,a,b))
require(motifs==[(0,0,4096,4,195,a,2) for a in (2840,2842,2845,2847)])

def F(m,x,y):
    _,B,T,S=D[m]
    return B+T*x+S*y+Q*x*y+h*(p(x)+p(y))

steps=[]
for c,j,d,e,m,a,b in motifs:
    K=F(m,a,b)-F(35,0,b0)
    steps.append(K)
    for u in [0,1,2,3,128,177,84608,2**128+1]:
        require(F(m,d*u+a,e*u+b)-F(35,u,M*u+b0)==K)
require(steps==[-11126,-12118,-13786,-14958] and gcd(*steps)==2)
print("All-order dilation certificate: four motifs",steps)
