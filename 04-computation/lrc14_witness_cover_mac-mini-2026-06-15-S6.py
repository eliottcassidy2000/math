from fractions import Fraction as F
from math import gcd
import random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r
def g(S, t):
    return min(nrm(v*t) for v in S)

# WITNESS-COVER (PROVE-side idea): is there a small finite set of rational
# witnesses tau* s.t. every primitive 13-set has some tau* with g(tau*)>=1/14?
# Witness pool a/N, N<=30.
witnesses=[]
for N in range(8,31):
    for a in range(1,N//2+1):
        if gcd(a,N)==1:
            witnesses.append(F(a,N))
witnesses=sorted(set(witnesses))

random.seed(5)
pool=list(range(1,40))
sample=[]
for _ in range(2000):
    S=sorted(random.sample(pool,13))
    gg=0
    for x in S: gg=gcd(gg,x)
    sample.append([x//gg for x in S])

# greedy set cover
uncovered=set(range(len(sample)))
cover=[]
while uncovered and len(cover)<12:
    best_w=None; best_cov=set()
    for w in witnesses:
        cov={i for i in uncovered if g(sample[i],w)>=F(1,14)}
        if len(cov)>len(best_cov):
            best_cov=cov; best_w=w
    if not best_cov: break
    cover.append((str(best_w),len(best_cov)))
    uncovered-=best_cov
print("Greedy witness cover over 2000 random primitive 13-sets:")
print("  witnesses (tau, #covered):", cover)
print("  remaining uncovered:", len(uncovered))
print("  -> a SMALL finite witness set certifies M>=1/14 for the vast majority;")
print("     tau=1/14 alone covers all sets with no multiple of 14.")
# count how many sample sets have a multiple of 14
withm14=sum(1 for S in sample if any(x%14==0 for x in S))
covered_by_1_14=sum(1 for S in sample if g(S,F(1,14))>=F(1,14))
print(f"  sample sets with a multiple of 14: {withm14}/2000")
print(f"  sample sets certified by tau=1/14 alone: {covered_by_1_14}/2000")
