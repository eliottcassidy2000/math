# S43 referee post-repair spot check: does Almkvist-Zudilin agree with Cooper s18 at (p-1)/2 mod p (twist (-1/p))?
from math import comb, factorial
AZ = lambda n: sum((-1)**k * 3**(n-3*k) * comb(n, 3*k) * comb(n+k, n) * factorial(3*k)//factorial(k)**3 for k in range(n//3+1))
S18 = lambda n: 1 if n == 0 else sum((-1)**k * comb(n,k) * comb(2*k,k) * comb(2*n-2*k,n-k) * (comb(2*n-3*k-1,n) + comb(2*n-3*k,n)) for k in range(n//3+1))
ps = [p for p in range(5, 800) if all(p % d for d in range(2, int(p**.5)+1))]
res = {}
for name, eps in (('+1', lambda p: 1), ('(-1/p)', lambda p: 1 if p % 4 == 1 else -1), ('(-3/p)', lambda p: 1 if p % 3 == 1 else -1), ('(3/p)', lambda p: pow(3, (p-1)//2, p) if pow(3,(p-1)//2,p)==1 else -1)):
    bad1 = [p for p in ps if (AZ((p-1)//2) - eps(p) * S18((p-1)//2)) % p]
    bad2 = [p for p in ps if (AZ((p-1)//2) - eps(p) * S18((p-1)//2)) % (p*p)]
    res[name] = ('mod p failures: %d (first %s)' % (len(bad1), bad1[:5]), 'mod p^2 failures: %d' % len(bad2))
for k, v in res.items():
    print('AZ((p-1)/2) vs %s * s18((p-1)/2), primes 5..800 (p=5,7 incl.):' % k, v)
