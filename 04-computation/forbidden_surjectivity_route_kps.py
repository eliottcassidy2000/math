"""
forbidden_surjectivity_route_kps.py

Honest test of whether the MULTIPLICATIVE semigroup picture alone can close the
'{7,21} are the ONLY gaps' direction (i.e. every odd >=23 is achievable), WITHOUT
needing the still-conjectural 'strong spectrum eventually covers all odd' fact.

Idea: The achievable set is the multiplicative monoid <1, strong-values>. If we
restrict to the strong values we ALREADY know rigorously (m<=6, exhaustive):
   irr = {3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45}
plus the products, do we get a CO-FINITE subset of the odd numbers (only finitely
many odd gaps, namely {7,21,...transients})?  If the odd gaps are infinite with just
these irreducibles, then multiplicativity ALONE is NOT enough and surjectivity needs
new strong values at every scale (the conjectural part).

A multiplicative monoid in the ODD numbers is co-finite iff its generators include,
for every prime p, ... no: multiplicatively you can only ever reach numbers whose
prime factors all appear among prime factors of generators. PRIMES not dividing any
generator are NEVER reachable. So large PRIMES are reachable ONLY if they are
themselves strong values. => multiplicativity ALONE can never be co-finite; you NEED
the strong spectrum to contain (eventually) every odd prime.  We demonstrate this.
"""
from sympy import isprime, primerange

# rigorously-known strong values (m<=6 exhaustive)
irr = sorted({3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45})
prime_factors_of_gens=set()
for g in irr:
    x=g; d=3
    while d*d<=x:
        while x%d==0:
            prime_factors_of_gens.add(d); x//=d
        d+=2
    if x>1: prime_factors_of_gens.add(x)
print("Primes dividing some m<=6 strong value:", sorted(prime_factors_of_gens))
print()

# Which odd primes are NOT reachable as products of these generators?
unreachable_primes=[p for p in primerange(3,200) if p not in irr and p not in prime_factors_of_gens]
print("Odd primes < 200 NOT a known strong value AND not a factor of one (=> unreachable")
print("by multiplicativity from m<=6 generators alone):")
print("  ", unreachable_primes[:25], "...")
print()
print("KEY POINT (honest): a prime p is achievable as H ONLY IF some STRONG tournament")
print("has H=p (since p prime => single strong component with H=p). So 'every odd >=23")
print("is achievable' REQUIRES that every odd prime >=23 is eventually a STRONG value.")
print("That is the surjectivity/density CONJECTURE (verified through ~n=9, not proved).")
print("Multiplicative closure CANNOT supply primes; it only recombines existing ones.")
print()
print("Therefore the proof splits cleanly into TWO independent halves:")
print("  (I)  {7,21} are NOT achievable      <- strong-min growth (Busch) -- RIGOROUS, DONE.")
print("  (II) every other odd IS achievable  <- strong spectrum eventually hits every")
print("       odd prime/value >=23           -- CONJECTURE (strong computational support).")
print()
print("The cluster expansion cleans up (I) for H=7 (unique K3 profile) but for 21 only via")
print("multiplicativity; it says NOTHING toward (II). So {7,21}-are-the-ONLY-gaps is NOT")
print("fully closed by the cluster picture: the open half is surjectivity, untouched by it.")
