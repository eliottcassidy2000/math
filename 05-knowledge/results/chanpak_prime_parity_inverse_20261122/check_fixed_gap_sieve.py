"""Exact finite replay for the literal pair sieve; no floating-point checks."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import isqrt


def check(condition, context):
    if not condition:
        raise RuntimeError(context)


def prime_divisors(n):
    primes = []
    d = 2
    while d <= isqrt(n):
        if n % d == 0:
            primes.append(d)
            while n % d == 0:
                n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        primes.append(n)
    return primes


def then(first, last):
    a, b = first
    c, d = last
    return c*a, c*b+d


def evaluate(f, v):
    a, b = f
    return a*v+b


A = (82, 60)
B = (5976, 4396)
C = (96264, 71856)
P = (292, -51)
blocks = (then(A, B), then(B, A))
check(blocks == ((490032, 362956), (490032, 360532)), blocks)
u = (evaluate(P, 1)-1)//2
finish = 20311*u+8654, 14572*u+6207
check((u, finish) == (120, (2445974, 1754847)), (u, finish))
S = (1, 0)
for letter in (A,A,C,A,A,B):
    S = then(S, letter)
F = then(S, finish)
check(F == (63618252914778017089536,47124685853114258307047), F)
E = sorted(set(prime_divisors(2*2424*blocks[0][0])) | set(prime_divisors(F[0])))
check(E == [2,3,7,41,47,83,101,191,26021], E)
fixed_pairs = []
for p in E:
    if p != 101:
        check(F[0] % p == 0, p)
        fixed_pairs.append((F[1] % p, (F[1]+2) % p))
check(fixed_pairs == [(1,1),(2,1),(3,5),(33,35),(8,10),(3,5),(2,4),(11440,11442)], fixed_pairs)
check(all(a and b for a,b in fixed_pairs), fixed_pairs)
check(all((c%101,d%101)==(81,63) for c,d in blocks), blocks)
cycle = set()
v = 0
while v not in cycle:
    cycle.add(v)
    check(all((evaluate(F,v)+j)%101 for j in (0,2)), v)
    v = (81*v+63)%101
check((len(cycle),v) == (25,0), (len(cycle),v))

samples = 0
for n in range(1,13):
    vals = []
    for word in product((0,1), repeat=n):
        v = 0
        for bit in word:
            v = evaluate(blocks[bit], v)
        vals.append(v)
    check(len(set(vals))==2**n,n)
    for q in (5,11,121,17,55,95,275):
        convolution = Counter({0:1})
        for power in range(n):
            new = Counter()
            for r, weight in convolution.items():
                for _, digit in blocks:
                    new[(r+digit*pow(blocks[0][0],power,q))%q] += weight
            convolution = new
        check(convolution==Counter(v%q for v in vals), (n,q))
    samples += len(vals)


def parity_capped(k,p,e):
    sign = 1
    for _ in range(e):
        if k%p:
            break
        sign = -sign
        k //= p
    return sign


for p,e in ((5,1),(5,5),(11,3),(13,3)):
    sums = Counter()
    for v in range(p**e):
        a,b = (parity_capped(evaluate(F,v)+j,p,e) for j in (0,2))
        sums[1] += a
        sums[2] += b
        sums[3] += a*b
    odd_valuation = sum(Fraction(p-1,p**(k+1)) for k in range(1,e,2))
    if e%2:
        odd_valuation += Fraction(1,p**e)
    for mask,count in sums.items():
        check(Fraction(count,p**e)==1-2*mask.bit_count()*odd_valuation, (p,e,mask,count))

print("exceptional primes, 101 orbit and physical chart constants: OK")
print(f"{samples} words: injectivity and exact modular convolution OK")
print("exact prime-power joint capped-sign factors: OK")
