"""Part 4: local density product. d_p = lim_k f(p^k). Compare prod_p d_p to actual
lonely measure mu = lim over q of (# lonely a in units)/phi(q), and the true continuous
measure of the safe set under the orbit.

Continuous lonely measure: mu = measure{ t in [0,1] : all runners safe }. By equidistribution
of (v_i t) when gcd structure generic, mu relates to product of local densities.
We estimate mu by f(N) for highly composite large N, and compare to prod over small primes
of d_p (using deepest f(p^k)).
"""
from lrc import lonely_count_at_q, f_pk
from fractions import Fraction

def sieve_primes(limit):
    s=[True]*(limit+1); s[0]=s[1]=False
    for i in range(2,int(limit**0.5)+1):
        if s[i]:
            for j in range(i*i,limit+1,i): s[j]=False
    return [i for i in range(2,limit+1) if s[i]]

SETS = {
    "m=3 {2,3,7}":[2,3,7],
    "m=3 {1,3,4}":[1,3,4],
    "m=3 {3,5,7}":[3,5,7],
    "m=4 {2,3,5,7}":[2,3,5,7],
    "m=4 {1,3,4,7}":[1,3,4,7],
}
SMALL_PRIMES=[2,3,5,7,11,13]

def d_p_estimate(speeds,p,maxq=4000):
    # deepest level
    k=1; last=None
    while p**k<=maxq:
        frac,_,_=f_pk(speeds,p,k); last=frac; k+=1
    return last

def mu_estimate(speeds, N):
    # use a highly composite N; measure over units mod N approximates local-global product
    cnt,tot=lonely_count_at_q(speeds,N)
    return Fraction(cnt,tot)

# highly composite N coprime-rich
N_list=[2520, 5040, 27720]  # 2^3*3^2*5*7 etc
for name,speeds in SETS.items():
    print(f"\n=== {name} ===")
    dps={}
    prod=Fraction(1)
    for p in SMALL_PRIMES:
        d=d_p_estimate(speeds,p)
        dps[p]=d
        prod*=d
        print(f"  d_{p} ~ {float(d):.5f}")
    # product over primes dividing the relevant range
    print(f"  PRODUCT prod_p d_p (p<=13) = {float(prod):.6f}")
    for N in N_list:
        mu=mu_estimate(speeds,N)
        print(f"  mu via units mod N={N}: {float(mu):.6f}")
