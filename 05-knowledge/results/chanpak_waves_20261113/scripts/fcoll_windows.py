# fcoll_windows.py (wave 22): check the window collision equation on every pair listed by fcollide.c and
# print the distributions of the lowest differing bit j, the window length n1 and log2(V/2^j).
#   A rho + B V = 917 2^j (A/9) - (Psi(M')-Psi(M)) - v0 - sigma*eps
# Usage: python3 fcoll_windows.py pairs.txt
import sys, collections
def Psi(x):
    b=[r for r in range(x.bit_length()) if (x>>r)&1]; s=0
    for i in range(len(b)):
        for k in range(i+1,len(b)):
            d=b[k]-b[i]; s+= 2 if d==1 else 1<<(d-1)
    return s
def rev(x,n): return int(format(x,'0%db'%n)[::-1],2) if n>0 else 0
pairs=[tuple(map(int,l.split())) for l in open(sys.argv[1] if len(sys.argv)>1 else 'coll24.txt')]
stats=collections.Counter(); ok=0; rows=[]
for u,up in pairs:
    x=u^up; m=x.bit_length(); j=(x&-x).bit_length()-1
    V=u>>m; W=u&((1<<m)-1); Wp=up&((1<<m)-1)
    assert (Wp>>(m-1))&1==1 and (W>>(m-1))&1==0
    n1=m-j; M=W>>j; Mp=Wp>>j; L=u&((1<<j)-1)
    A=Mp-M; B=rev(Mp,n1)-rev(M,n1); rho=rev(L,j)
    sigma=(Mp&1)-(M&1); eps=((L>>(j-1))&1 if j>0 else 1)
    lhs=A*rho+B*V; rhs=917*(1<<j)*(A//9)-(Psi(Mp)-Psi(M))-(V&1)-sigma*eps
    ok+= (lhs==rhs and A%9==0)
    rows.append((u,up,m,j,n1,V,L,M,Mp,A,B))
print("identity (***) holds on",ok,"of",len(pairs))
c=collections.Counter((r[3]) for r in rows); print("j distribution:",sorted(c.items())[:12])
c=collections.Counter((r[4]) for r in rows); print("n1 distribution:",sorted(c.items()))
import math
# V relative to 2^j
c=collections.Counter(round(math.log2(max(r[5],1)/2**r[3])) for r in rows); print("log2(V/2^j):",sorted(c.items()))
c=collections.Counter((r[7]==0) for r in rows); print("M==0 share:",c)
c=collections.Counter(bin(r[7]).count('1') for r in rows); print("pop(M):",sorted(c.items())[:15])
c=collections.Counter(bin(r[8]).count('1') for r in rows); print("pop(M'):",sorted(c.items())[:15])
c=collections.Counter(r[9]//9 for r in rows); print("A/9 most common:",c.most_common(12))
c=collections.Counter(r[10] for r in rows); print("B most common:",c.most_common(12))
