"""Exact-integer, dependency-free fixed-kernel amicable-pair certificates."""
from math import gcd, isqrt
import json

def require(condition):
    if not condition: raise AssertionError("certificate check failed")

def sigma_independent(n):
    total=0
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            total += d
            if d*d != n: total += n//d
    return total

def factor(n):
    out = {}; p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p,0)+1; n //= p
        p += 1
    if n > 1: out[n] = out.get(n,0)+1
    return out

def sigma(n):
    v = 1
    for p,e in factor(n).items(): v *= (p**(e+1)-1)//(p-1)
    return v

def prime(n):
    return n >= 2 and all(n%d for d in range(2,isqrt(n)+1))

def divisors(n):
    ds = [1]
    for p,e in factor(n).items(): ds = [d*p**i for d in ds for i in range(e+1)]
    return sorted(ds)

def certificates(a):
    S = sigma(a); D = 2*a-S; E = S-a
    if D <= 0: return []
    out = []
    for d in divisors(a*a):
        e = a*a//d
        if d >= e: break
        if (d+E)%D or (e+E)%D: continue
        p,q = (d+E)//D,(e+E)//D
        r = p*q+p+q
        if not (prime(p) and prime(q) and prime(r) and gcd(a,p*q*r)==1): continue
        n,m = a*p*q,a*r
        require(sigma(n)==sigma(m)==n+m)
        require((D*p-E)*(D*q-E)==a*a)
        if p>3 and q>3: require(r%12==11)
        out.append(dict(a=a,S=S,D=D,E=E,d=d,e=e,p=p,q=q,r=r,n=n,m=m))
    return out

def main():
    # Exhaustive over the stated kernel range, NOT over all amicable pairs.
    found = [c for a in range(1,5001) for c in certificates(a)]
    require(any((c['n'],c['m'])==(220,284) for c in found))
    for c in found:
        require(sigma_independent(c['n'])==sigma_independent(c['m'])==c['n']+c['m'])
    # Independent brute-force equivalence on small kernels and p,q <= 100.
    ps = [p for p in range(2,101) if prime(p)]
    checked = 0
    for a in range(1,101):
        produced = {(c['p'],c['q']) for c in certificates(a)}
        for p in ps:
            for q in ps:
                if p >= q: continue
                r=p*q+p+q
                if not prime(r) or gcd(a,p*q*r)!=1: continue
                n,m=a*p*q,a*r
                good=sigma(n)==sigma(m)==n+m
                require(good == ((p,q) in produced))
                checked += 1
    # Independent modular sieve exhaustive over reduced residue classes.
    residue_results=[]
    for p in (1,5,7,11):
        for q in (1,5,7,11):
            r=(p*q+p+q)%12
            require(r in (3,11))
            residue_results.append([p,q,r])
    print(json.dumps(dict(kernel_bound=5000,certificate_count=len(found),
                         odd_kernel_count=sum(c['a']%2 for c in found),
                         brute_force_cases=checked,residue_table=residue_results,
                         certificates=found),indent=2))

if __name__=='__main__': main()
