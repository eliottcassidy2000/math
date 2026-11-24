"""Physical height overlap witnesses and exact bulk readout arithmetic."""
from math import comb, isqrt, prod
from itertools import combinations

def rows(n, mask=0, reverse_last=False):
    free=[(i,j) for i,j in combinations(range(n),2) if j>i+1]
    back={ij for k,ij in enumerate(free) if mask>>k&1}
    if reverse_last: back.add((0,n-1))
    out=[0]*n
    for i,j in combinations(range(n),2):
        if (i,j) in back:out[j]|=1<<i
        else:out[i]|=1<<j
    return tuple(out)

def height(A):
    n=len(A); M=(1<<n)-1
    F=[[0]*n for _ in range(M+1)]
    for i in range(n):F[1<<i][i]=1
    for s,E in enumerate(F):
        for i,count in enumerate(E):
            if not count:continue
            t=A[i]&~s
            while t:
                b=t&-t;t^=b
                F[s|b][b.bit_length()-1]+=count
    return sum(F[-1])

def prime(n):
    return n>=2 and (n==2 or n%2!=0 and all(n%p for p in range(3,isqrt(n)+1,2)))

if __name__=='__main__':
    chosen=None
    for n in range(3,7):
        for mask in range(1<<((n-1)*(n-2)//2)):
            T=rows(n,mask)
            if height(T)==43:
                chosen=(n,mask,T);break
        if chosen:break
    assert chosen
    print('alternate factor 43 witness',chosen)
    for d in range(1,9):assert height(rows(d+2,reverse_last=True))==2**d+1
    for d in (2,3,7,57):
        n=d*d+1; path=n*comb(d,2);one=n*d//2*(n-2*d);center=n*comb(d,3)
        print('degree',d,'powers',3**d,3**(2*d-2),3**(2*d),'extreme arc',2**d+1,
              'triple census',path,one,center,comb(n,3)-path-one-center)
    for a,deck in [
      (57,[3,3,571,174763,160465489]),
      (112,[449,2689,65537,183076097,358429848460993]),
      (114,[5,13,229,457,131101,160969,525313,275415303169])]:
        assert prod(deck)==2**a+1 and all(map(prime,deck))
        print('factorization passed',a,deck)
