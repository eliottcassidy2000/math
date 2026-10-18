# S43 referee: valuations of A(n), s7(n) at Frobenius points for p = 11 (non-ordinary for eta(2t)^4 eta(4t)^4: a(11) = -44)
from math import comb
def s7(n): return sum(comb(n,k)**2*comb(n+k,k)*comb(2*k,n) for k in range((n+1)//2, n+1))
def A(n): return sum((comb(n,k)*comb(n+k,k))**2 for k in range(n+1))
def vp(x,p):
    v=0
    while x%p==0: x//=p; v+=1
    return v
for p in (5, 7, 11, 13):
    for r in (1, 2, 3):
        n=(p**r-1)//2
        print('p=%d r=%d n=%d: v_p(A(n)) = %d, v_p(s7(n)) = %d, v_p(s7 - (-1)^n A) = %d' % (p, r, n, vp(A(n),p), vp(s7(n),p), vp(s7(n)-(-1)**n*A(n),p)))
