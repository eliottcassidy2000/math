#!/usr/bin/env python3
"""Finite checks of the reverse-prefix decomposition; not used for asymptotics."""
from collections import Counter, defaultdict

def reversal(z, width):
    ans = 0
    for _ in range(width):
        z, bit = divmod(z, 2)
        ans = (ans << 1) + bit
    assert z == 0
    return ans

def chart(z, z1, radix):
    u, b = divmod(z,radix)
    u1, b1 = divmod(z1,radix)
    assert abs(u-u1) <= 1
    a = radix
    if u != u1:
        first = min(u,u1)
        a <<= (first ^ (first+1)).bit_length()
    top, bot = divmod(z,a)
    top1, bot1 = divmod(z1,a)
    assert top == top1
    return a,top,bot,bot1

L=15
N=1<<L
p=[0]*N
for n in range(1,N):
    w,b=divmod(n,2)
    p[n]=p[w]+(0 if b else w+(w&1))
    assert p[n]%2 ==0 and p[n]<=n

for c,k,phi in [(18,3,2),(189,11,0),(2,1,0)]:
    fibers=defaultdict(list)
    G=[c*n+k*p[n]+phi*(n&1) for n in range(N)]
    for n in range(N):
        fibers[G[n]].append(n)
    chart_counts=Counter(); count=0; annulus_count=0; smallreverse=0; short=0
    for grp in fibers.values():
        for i,n in enumerate(grp):
            for n1 in grp[:i]:
                ell=(n ^ n1).bit_length(); V=1<<ell
                w,b=divmod(n,V); w1,b1=divmod(n1,V)
                assert w==w1 and b >= V//2 > b1
                z,z1=reversal(b,ell),reversal(b1,ell); d=z-z1
                assert d&1 and abs(d)*(1+w)*k <= (3*k+c+abs(phi))*V
                assert k*d*w+k*(w&1) == c*(b-b1)+k*(p[b]-p[b1])+phi*((b&1)-(b1&1))
                assert d!=-1
                if d==1 and ell>=2:
                    S=V//2
                    z0=reversal(b1,ell-1)
                    assert k*(w+(w&1)) == c*S+k*(S-1-z0+(z0%2==0))
                D=min(V,(3*k+c+abs(phi))*V//(k*(w+1)))
                radix=1<<(D-1).bit_length()
                a,top,bot,bot1=chart(z,z1,radix)
                aa=a.bit_length()-1; s=ell-aa; S=1<<s
                t=reversal(top,s)
                U=a*w+reversal(bot,aa); U1=a*w+reversal(bot1,aa)
                assert (n,n1)==(S*U+t,S*U1+t) and U!=U1
                if s:
                    rhs = c*S*(U-U1)+k*(p[U]-p[U1]+(S-1-reversal(t,s))*(U-U1)
                                              +(t<S//2)*((U&1)-(U1&1)))
                    assert rhs==0
                key=(ell,w,radix,a,bot,bot1,int(s and t<S//2))
                chart_counts[key]+=1
                assert chart_counts[key]==1
                count+=1
                if n-n1 <= 1<<(L//2): short+=1
                if n1>=N//2:
                    annulus_count+=1
                    if abs(d)<=1<<(L//3): smallreverse+=1
    print(c,k,phi,'unordered collisions',count,'short',short,'dyadic-outer-half',annulus_count,
          'small reverse within outer half',smallreverse,'unique half-charts',len(chart_counts))
