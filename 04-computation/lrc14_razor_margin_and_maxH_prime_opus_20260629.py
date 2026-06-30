from fractions import Fraction as F
from math import gcd
def norm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def M_exact(S):
    C=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): C.add(F(m,d))
        for m in range(2*S[i]+1): C.add(F(m,2*S[i]))
    best=F(0)
    for t in C:
        if 0<t<1:
            v=min(norm(s*t) for s in S)
            if v>best: best=v
    return best
thr=F(1,14)
print("Divisor-loaded covering family {1..11,13, 84*m}: does M -> 1/14 (razor) or stay (margin)?")
for m in [1,2,3,4,5,6,8]:
    S=[1,2,3,4,5,6,7,8,9,10,11,13,84*m]
    M=M_exact(S)
    print(f"  m={m:2d} (last={84*m:4d}): M={M}={float(M):.6f}  margin above 1/14 = {float(M-thr):+.6f}")
print(f"  (1/14={float(thr):.6f}, 7/89={float(F(7,89)):.6f})")
print()

# (b) WHY max-H is a single strongly-connected prime: H(X(+)Y)=H(X)H(Y); single prime beats any product
def Hadj(n,adj):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c or not (mask>>last)&1: continue
            av=adj[last]&~mask
            while av:
                nb=av&-av;nx=nb.bit_length()-1;dp[mask|nb][nx]+=c;av^=nb
    return sum(dp[(1<<n)-1])
# max H over SC primes (single) vs max product over compositions, per size
# H-atoms (max SC-prime H) from earlier: size 3->3,4->5,5->15,6->45 ; products multiply
maxprime={1:1,2:0,3:3,4:5,5:15,6:45}   # max H of a strongly-connected tournament of that size (0=none)
def best_product(n):
    # max over compositions of n into parts>=1 (parts with a prime) of prod maxprime
    from functools import lru_cache
    @lru_cache(None)
    def f(m):
        if m==0: return 1
        best=0
        for k in range(1,m+1):
            if maxprime.get(k,0)>0:
                best=max(best, maxprime[k]*f(m-k))
        return best
    return f(n)
print("(b) max-H = single SC prime, because condensation imposes a NO-RETURN constraint (H multiplies):")
for n in [3,4,5,6]:
    sp=maxprime[n]; bp=best_product(n)
    print(f"  n={n}: max single SC prime H={sp};  max H over ANY composition (products)={bp};  "
          f"single wins: {sp>=bp} (a connected tournament beats every X(+)Y split)")
print("  => H super-multiplicative in prime size: a single irreducible (SC) tournament lets Ham paths")
print("     INTERLEAVE freely; condensation X(+)Y forces 'finish X before Y' => H(X)*H(Y) <= H(connected).")
