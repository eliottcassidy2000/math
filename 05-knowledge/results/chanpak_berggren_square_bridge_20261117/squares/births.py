from math import gcd,isqrt
from heapq import heappop,heappush
def birth(m,n):
    a,b=m*m-n*n,2*m*n
    return max(a,b)**2
def tree(N):
    heap=[(16,2,1)];out=[]
    while heap:
        w,m,n=heappop(heap)
        if w>N:break
        out.append((m,n))
        for a,b in ((2*m-n,m),(2*m+n,m),(m+2*n,n)):
            t=birth(a,b)
            if t<=w:raise ArithmeticError((m,n,a,b))
            heappush(heap,(t,a,b))
    return out
def euclid(N):
    # max(a²,b²)≤N implies m⁴≤2N because hypotenuse=m²+n².
    return {(m,n) for m in range(2,isqrt(isqrt(2*N))+1) for n in range(1,m)
            if (m+n)%2 and gcd(m,n)==1 and birth(m,n)<=N}
if __name__=='__main__':
    N=10**8;q=tree(N)
    if len(q)!=1788 or set(q)!=euclid(N):raise ArithmeticError(len(q))
    if [birth(*v) for v in ((3,2),(5,2),(4,1))]!=[144,441,225]:raise ArithmeticError()
    if not 63**2+16**2==33**2+56**2==65**2:raise ArithmeticError()
    print('1788 identical literal first-birth edges through 10^8; child stamps and norm collision pass.')
