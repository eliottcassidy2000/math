"""Terras-transition reading of Korselt: for odd n, T(x) = (3x+1)/2 is an affine bijection of Z/n with fixed point -1,
conjugate to multiplication by 3/2.  T^(n-1) = id on Z/n  iff  (3/2)^(n-1) = 1 mod n  iff  n | 3^(n-1) - 2^(n-1).
Composite such n are OEIS A073631 (nonprimes k with k | 3^(k-1) - 2^(k-1)); every Carmichael number coprime to 6 is one
(R. Israel's comment there); 561 is excluded by its prime factor 3 (the multiplier).  Census below 2*10^6 + direct Terras checks."""
from sympy import isprime, factorint
from sympy.ntheory import n_order
N = 2_000_000
ps = [n for n in range(9, N, 2) if not isprime(n) and pow(3, n-1, n) == pow(2, n-1, n)]
print("composite (3/2)-pseudoprimes below", N, ":", len(ps)); print("first 40:", ps[:40])
car = [n for n in range(561, N, 2) if not isprime(n) and all(e == 1 for e in factorint(n).values()) and all((n-1) % (p-1) == 0 for p in factorint(n))]
S = set(ps)
print("Carmichael below N:", len(car), "; coprime to 6:", sum(1 for n in car if n % 3), "all in the list:", all(n in S for n in car if n % 3))
print("Carmichael divisible by 3 (excluded):", [n for n in car if n % 3 == 0], "none in the list:", not any(n in S for n in car if n % 3 == 0))
nsf = [n for n in ps if any(e > 1 for e in factorint(n).values())]
print("non-squarefree members:", nsf, "(all multiples of 23^2 = 529? ", all(n % 529 == 0 for n in nsf), ")")
sq = [n for n in ps if n not in nsf]
print("squarefree members satisfy  ord_p(3/2) | n-1 for every p | n:", all(all((n-1) % n_order(3*pow(2,-1,p) % p, p) == 0 for p in factorint(n)) for n in sq))
# direct Terras check of T^(n-1) = id for the first members and for 561, 91
def terras_identity(n):
    inv2 = pow(2, -1, n)
    for x in range(n):
        y = x
        for _ in range(n-1): y = ((3*y+1)*inv2) % n
        if y != x: return False
    return True
for n in (65, 133, 529, 1105, 1729, 91, 561):
    print(f"T^(n-1) = id mod {n}: {terras_identity(n)}  ({'in list' if n in S else 'not in list'})")
# first member 65 = 5 * 13: orders of 3/2
for p in (5, 13): print(f"ord_{p}(3/2) =", n_order(3*pow(2,-1,p) % p, p))
