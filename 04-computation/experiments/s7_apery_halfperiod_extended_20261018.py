import sys, json
from math import comb
def A(n): return sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
def s7(n): return sum(comb(n,k)**2*comb(n+k,k)*comb(2*k,n) for k in range((n+1)//2,n+1))
primes=[p for p in range(5,1500) if all(p%d for d in range(2,int(p**.5)+1))]
fails=[]
for p in primes:
    n=(p-1)//2
    if (s7(n)-(-1)**n*A(n))%(p*p): fails.append(p)
print('half-period: primes tested',len(primes),'max',primes[-1],'failures',fails,flush=True)
res={'half':dict(tested=len(primes),pmax=primes[-1],fails=fails)}
hi=[]
for p in (5,7,11,13):
    for r in (2,3):
        n=(p**r-1)//2
        if n>1300: continue
        ok=(s7(n)-(-1)**n*A(n))%(p*p)==0
        hi.append((p,r,n,ok)); print('p=%d r=%d n=%d: s7(n) == (-1)^n A(n) mod p^2: %s'%(p,r,n,ok),flush=True)
res['higher']=hi
json.dump(res,open(sys.argv[1],'w'),indent=1)
