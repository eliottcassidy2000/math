# Multiply perfect numbers sigma(n) = k n up to a bound (sieve), their prime factors, and two numerology controls.
import json
B = 2_000_000
sig=[0]*(B+1)
for d in range(1,B+1):
    for m in range(d,B+1,d): sig[m]+=d
mp={}
for n in range(1,B+1):
    if sig[n]%n==0 and sig[n]//n>=2: mp.setdefault(sig[n]//n,[]).append(n)
def factor(n):
    f=[]; p=2
    while p*p<=n:
        while n%p==0: f.append(p); n//=p
        p+=1
    if n>1: f.append(n)
    return f
ss={2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}   # supersingular primes
out={}
for k,ns in sorted(mp.items()):
    out[k]={n: factor(n) for n in ns}
print(json.dumps({'bound':B,'k_perfect':{k:list(v) for k,v in out.items()},
  'all_prime_factors_supersingular': all(set(factor(n))<=ss for ns in mp.values() for n in ns),
  'first_nonsupersingular_factor': next(((n,p) for ns in mp.values() for n in ns for p in factor(n) if p not in ss), None)}, indent=0))
