"""Exact arithmetic inverse of the symmetric physically constructed factor deck."""
from math import isqrt

def invert(N,U,C,k):
    P=N-C
    u=-((P-U*U)//(2*U))
    a=P-U*U+2*U*u
    if a%4: return None
    v=a//4
    D=u*u-4*v
    if D<0: return None
    m=isqrt(D)
    if m*m!=D or (u-m)%2: return None
    i=(u-m)//2; j=(u+m)//2
    if 0<=i<=j<=k: return i,j

def profile(k):
    S=1957; v=0
    for i in range(k):
        v+=S; S=4096*S+12645
    return 898+S+1655*v+5175*k,2*S

checked=0
for k in range(16):
    U,S=profile(k); C=(2*(45+S))**2
    for i in range(k+1):
        for j in range(i,k+1):
            N=(U-2*i)*(U-2*j)+C
            if invert(N,U,C,k)!=(i,j): raise AssertionError((k,i,j,N))
            if any(invert(N+d,U,C,k) is not None for d in (-2,2)):
                raise AssertionError((k,i,j,N,'unexpected adjacency'))
            checked+=1
U,S=profile(1); C=(2*(45+S))**2; N=U*(U-2)+C
print('Exact deck inversions and adjacent controls:',checked)
print('Height one target',N,'slide positions',invert(N,U,C,1),'next odd',invert(N+2,U,C,1))
