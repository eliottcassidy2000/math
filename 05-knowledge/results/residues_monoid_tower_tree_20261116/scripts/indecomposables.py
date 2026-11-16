#!/usr/bin/env python3
"""indecomposables.py N: in the monoid M = (odd numbers other than 7 and 21, x), list the indecomposables below N
(elements > 1 that are not a product of two elements > 1 of M), and compare with the prediction: primes other than 7,
7p for primes p >= 5, and 49, 63, 343.  Also the 'first 32 odd numbers' split and the counts per block of 32."""
import sys
N = int(sys.argv[1])
M = set(x for x in range(3, N, 2) if x not in (7, 21))
def is_prime(n): return n > 1 and all(n % d for d in range(2, int(n**0.5)+1))
indec = []
for x in sorted(M):
    dec = False
    d = 3
    while d * d <= x:
        if x % d == 0 and d in M and (x // d) in M: dec = True; break
        d += 2
    if not dec: indec.append(x)
pred = sorted(set([p for p in range(3, N) if is_prime(p) and p != 7] + [7*p for p in range(5, N//7) if is_prime(p)] + [49, 63, 343]))
pred = [x for x in pred if x < N]
print(f"N={N}: indecomposables = prediction: {indec == pred}; count {len(indec)}; extra: {sorted(set(indec)-set(pred))[:10]}; missing: {sorted(set(pred)-set(indec))[:10]}")
hidden = [x for x in indec if not is_prime(x)]
print("hidden generators (indecomposables that are not prime) below N:", hidden[:25], "...", f"({len(hidden)} of them; primes p >= 5 below N/7: {sum(1 for p in range(5, N//7) if is_prime(p))})")
# first 32 odd numbers 1..63: positions of 7 (4th) and 21 (11th), the split 3 | 28 and generators
odds = list(range(1, 64, 2))
print("first 32 odd numbers: 7 is the", odds.index(7)+1, "th, 21 the", odds.index(21)+1, "th; before 7:", [x for x in odds if x < 7], "; after 7:", len([x for x in odds if x > 7]), "numbers")
print("indecomposables among 1..63:", [x for x in indec if x < 64], "; of these hidden:", [x for x in hidden if x < 64])
for k in range(5, 13):
    top = 2**(k+1)
    blk = [x for x in indec if x < top]; hid = [x for x in hidden if x < top]
    print(f"  first 2^{k} odd numbers (< {top}): indecomposables {len(blk)} (primes {len(blk)-len(hid)}, hidden {len(hid)})")
