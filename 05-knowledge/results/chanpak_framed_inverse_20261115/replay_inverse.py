from math import isqrt, factorial
from collections import defaultdict
from itertools import combinations

def need(ok,*args):
    if not ok:raise AssertionError(args)
def p(e):
    z=0
    while e:
        e,b=divmod(e,2)
        if b==0:z+=e+(e&1)
    return z
def N(e):return 1+3*e+p(e)+(e&1)
def energy(e):
    sites=[j for j in range(e.bit_length()) if e>>j&1]
    return sum(max(2,1<<(j-i-1)) for i,j in combinations(sites,2))
def decode(t,width=None,r=None):
    q=1<<max(1,(isqrt(t)).bit_length()-1)
    l=q.bit_length()-1
    out=[]
    for b in range(q):
        Y=t-1-3*b-p(b)-(b&1)
        if Y<0:continue
        L=4*q-1-int(f'{b:0{l}b}'[::-1],2)
        for a in range(max(0,-(-(Y-1)//(L+1))),Y//L+1):
            e=q*a+b
            if e and N(e)==t and (width is None or e.bit_length()==width) and (r is None or e.bit_count()==r):
                out.append(e)
    return sorted(out)
def paths(rows):
    n=len(rows);M=1<<n
    D=[[0]*n for _ in range(M)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,M-1):
        for i,v in enumerate(D[s]):
            if not v:continue
            c=rows[i]&~s
            while c:
                b=c&-c;c-=b;D[s|b][b.bit_length()-1]+=v
    return sum(D[-1])
def chart(e,plus=False):
    m=e.bit_length()+1
    n=m+1+int(plus)
    # Ascending order is spine then u, or spine then v,u.
    U=n-1
    rows=[0]*n
    for i,j in combinations(range(n),2):
        rev=(j==U and i<m-1 and bool(e>>(m-2-i)&1))
        a,b=(j,i) if rev else (i,j)
        rows[a]|=1<<b
    f=[(1<<(n-1))-1-(1<<j) for j in range(n-1)]
    for i,j in combinations(range(n),2):
        if j>=i+2 and rows[j]>>i&1:
            f[i]&=~(1<<(j-1));f[j-1]&=~(1<<i)
    return rows,f

for k in range(4,7):
    for Y,X,sgn in ((27,6,1),(29,5,-1)):
        e=(Y<<k)+X
        t=(783<<(k-3))+sgn
        need(N(e)==t and e in decode(t,k+5,6),k,e)
        for plus in (False,True):
            R,G=chart(e,plus)
            d=len(G)-1-6
            need(paths(R)==(t if plus else t-2*e))
            need(paths(G)==d*(d+1)*factorial(len(G)-2))
table=defaultdict(list)
for e in range(1,4001):
    n=N(e)
    need(n==1+4*e-energy(e) and 1+3*e<=n<=1+4*e,e)
    table[n].append(e)
    for j in range(e.bit_length()):
        if e>>j&3==1:
            need(N(e+(1<<j))-n>=(3*(1<<j)+2 if j else 4))
for n in range(1,3000,2):need(decode(n)==table[n],n)
for n,e in ((63,19),(69,18)):
    R,G=chart(e,True);need(paths(R)==n and e in decode(n),n)
need(not decode(15))
print('exact inverse controls: 1500 targets, six two-scale physical outputs, both undirected charts')
print('first same-chart adjacent pair: 1565, 1567; target 63 -> (H,E)=(25,19); 69 -> (33,18)')

for t in range(0,21):
    A=1<<(t+2)
    for k in range(4,16):
        for x,y,shift in ((6,8*A-5,5),(5,8*A-3,3)):
            e=(y<<k)+x
            z=(8*A-5)*(8*A-3)*(1<<(k-3))+shift-A
            need(N(e<<t)==z,t,k)
            need(e.bit_count()==t+6)
            if t==1 and k==4:
                R,G=chart(e<<t,True)
                need(paths(R)==z)
                d=k+t
                need(paths(G)==d*(d+1)*factorial(len(G)-2))
print('504 two-parameter exact responses; literal double-scale check at t=1,k=4')
