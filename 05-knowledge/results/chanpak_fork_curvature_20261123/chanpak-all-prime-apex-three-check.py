#!/usr/bin/env python3
"""Exact mod-3-power two-reading channels for long guarded macro tiles."""
from itertools import product

def inv(w): return [-a for a in reversed(w)]
def mul(x,y,q):
    a,b,d,g=x; A,B,D,G=y
    return ((A+a*D)%q,(b+d*B)%q,(d*D)%q,(g+G+a*B)%q)
def mpow(X,n,q):
    z=(0,0,1,0)
    for i in range(n.bit_length()):
        if n>>i&1:z=mul(z,X,q)
        X=mul(X,X,q)
    return z
def exact_g(w,L,q):
    return sum(pow(L,j-i-1,q) for i,c in enumerate(w) if c=='O'
               for j in range(i+1,len(w)) if w[j]=='I')%q

visits=0
for e in range(1,5):
    q=3**e; E=q*(q//3*2)
    C='O'*e+'I'*e
    tiles=(C,'O'+C,C+'I')
    O=[2,-1];I=[-1,3]
    Z=[1]+inv(O*e+I*e)
    u=O+inv(I)
    U=O+u+inv(O)
    K=u+U+inv(u)+inv(U)
    h=pow(2,-1,q)
    channels=(2%q,0,(3*h)%q,h)
    for lam in (0,(3*h)%q):
        for w in product(range(3),repeat=3):
            assert exact_g(''.join(tiles[j] for j in w),lam,q)==sum(exact_g(tiles[j],lam,q) for j in w)%q
    def g(word,L):
        if L%3==0:
            return sum((1 if j>0 else -1)*exact_g(tiles[abs(j)-1],L,q) for j in word)%q
        z=(0,0,1,0)
        for j in word:
            y=(0,0,1,0)
            for t in tiles[abs(j)-1]:y=mul(y,(1,0,L,0) if t=='O' else (0,1,L,0),q)
            z=mul(z,y if j>0 else mpow(y,E-1,q),q)
        return z[-1]
    assert all(g(O,L)==0 for L in channels)
    def pair(w):
        a,b,c,d=(g(w,L) for L in channels)
        return ((a+b)%q,(c+d)*h*h%q)
    a=pair(Z); b=pair(K)
    assert a==(1,1) and b==(-3*h*h%q,3*h*h%q)
    for parity in (0,1):
        scale=pow(2,2-parity,q)
        for H,Q in product(range(q),repeat=2):
            if (Q-(-1)**parity*H)%3:continue
            X=(H-1)%q;Y=(Q*pow(scale,-1,q)-1)%q
            # Solve using difference, dividing by 3 before inversion.
            k=(((Y-X)%q)//3*2)%(q//3)
            z=(X-k*b[0])%q
            assert (1+z*a[0]+k*b[0])%q==H
            assert scale*(1+z*a[1]+k*b[1])%q==Q
            visits+=1
print(visits,'all compatible two-reading residues checked through 3^4, both parities')
