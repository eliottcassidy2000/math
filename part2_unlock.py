"""Part 2: which prime power unlocks vs composite q (CRT necessity).
Find smallest lonely q over (a) prime powers only, (b) all q. Uses ALL a (not just units)
since loneliness is about existence of any time a/q."""
from lrc import lonely_count_all_a, is_safe
from math import gcd

def first_prime_powers(limit):
    """list of prime powers <= limit, sorted."""
    primes = []
    n = 2
    while n <= limit:
        isp = all(n % p for p in primes if p*p <= n)
        # better sieve below
        n += 1
    return None

def sieve_primes(limit):
    sieve = [True]*(limit+1)
    sieve[0]=sieve[1]=False
    for i in range(2,int(limit**0.5)+1):
        if sieve[i]:
            for j in range(i*i,limit+1,i):
                sieve[j]=False
    return [i for i in range(2,limit+1) if sieve[i]]

def prime_powers_upto(limit):
    pps = []
    for p in sieve_primes(limit):
        q = p
        while q <= limit:
            pps.append((q,p))
            q *= p
    return sorted(pps)

def smallest_lonely_q_anytype(speeds, limit):
    """smallest q (any type) with a lonely time."""
    for q in range(2, limit+1):
        if lonely_count_all_a(speeds, q) > 0:
            return q
    return None

def smallest_lonely_prime_power(speeds, limit):
    for q,p in prime_powers_upto(limit):
        if lonely_count_all_a(speeds, q) > 0:
            return q,p
    return None

SETS = {
    "m=3 {1,2,3}": [1,2,3],
    "m=3 {1,3,4}": [1,3,4],
    "m=3 {2,3,7}": [2,3,7],
    "m=3 {1,4,5}": [1,4,5],
    "m=4 {1,2,3,4}": [1,2,3,4],
    "m=4 {1,3,4,7}": [1,3,4,7],
    "m=4 {2,3,5,7}": [2,3,5,7],
    "m=4 {1,5,6,8}": [1,5,6,8],
    "m=5 {1,2,3,4,5}": [1,2,3,4,5],
    "m=5 {1,3,4,5,9}": [1,3,4,5,9],
    "m=5 {1,2,3,4,6}": [1,2,3,4,6],
}
LIMIT = 5000

print(f"{'set':22s} {'smallest_q(any)':16s} {'smallest_prime_pow':18s} CRT_needed?")
for name, speeds in SETS.items():
    sq = smallest_lonely_q_anytype(speeds, LIMIT)
    spp = smallest_lonely_prime_power(speeds, LIMIT)
    if spp is None:
        ppstr = "NONE<=%d"%LIMIT
        crt = "YES (composite only)" if sq else "?"
    else:
        ppstr = f"{spp[0]} (p={spp[1]})"
        crt = "no (single prime)" if (sq and sq==spp[0]) else (f"maybe (any q={sq})")
    print(f"{name:22s} {str(sq):16s} {ppstr:18s} {crt}")
