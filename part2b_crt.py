"""Part 2b: systematically find CRT-only sets (no prime-power lonely time, but composite works).
Scan random/structured m=3,4,5 sets."""
from lrc import lonely_count_all_a
from itertools import combinations

def sieve_primes(limit):
    s=[True]*(limit+1); s[0]=s[1]=False
    for i in range(2,int(limit**0.5)+1):
        if s[i]:
            for j in range(i*i,limit+1,i): s[j]=False
    return [i for i in range(2,limit+1) if s[i]]

def prime_powers_upto(limit):
    pps=[]
    for p in sieve_primes(limit):
        q=p
        while q<=limit: pps.append(q); q*=p
    return set(pps)

LIMIT=2000
PPS = prime_powers_upto(LIMIT)

def classify(speeds):
    """return (smallest_q_any, smallest_q_primepower) within LIMIT."""
    sq=None; spp=None
    for q in range(2,LIMIT+1):
        if lonely_count_all_a(speeds,q)>0:
            if sq is None: sq=q
            if q in PPS and spp is None: spp=q
            if sq is not None and spp is not None: break
    return sq,spp

# enumerate sets {1=..} with small entries to find CRT-only cases
print("CRT-ONLY sets (no prime power works, composite does):")
found=0
for m,maxv in [(3,9),(4,9),(5,10)]:
    for combo in combinations(range(1,maxv+1), m):
        speeds=list(combo)
        sq,spp=classify(speeds)
        if sq is not None and spp is None:
            # composite-only
            print(f"  m={m} {speeds}: smallest_q={sq} (composite), no prime-power<= {LIMIT}")
            found+=1
            if found>=25: break
    if found>=25: break
print(f"total found (capped): {found}")
