# AUDIT a4b: reproduce the session's verification counts (472 odd squarefree composites <= 2500 with all prime
# factors <= 199, of which 56 pass; 142 even n <= 400 with squarefree odd part and odd primes <= 100).
from sympy import factorint, isprime
direct = {}
for f in ["a2_direct_4_2000.out", "a2_direct_2001_2500.out"]:
    for line in open(f):
        n, ok, *_ = line.split(); direct[int(n)] = int(ok)
odd = [n for n in range(9, 2501, 2) if not isprime(n) and all(e == 1 for e in factorint(n).values()) and max(factorint(n)) <= 199]
print("odd squarefree composites <= 2500 with prime factors <= 199:", len(odd), " passing (direct):", sum(direct[n] for n in odd))
even = [n for n in range(4, 401, 2) if all(e == 1 for p, e in factorint(n).items() if p > 2) and all(p <= 100 for p in factorint(n) if p > 2)]
print("even n <= 400 with squarefree odd part, odd primes <= 100:", len(even))
print("even n <= 400 with squarefree odd part (no prime bound):", sum(1 for n in range(4, 401, 2) if all(e == 1 for p, e in factorint(n).items() if p > 2)))
