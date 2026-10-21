# AUDIT a13 (independent auditor, 2026-10-21): squarefree Camion pseudoprimes n <= 20000 by the (proved) criterion,
# with D(r) mod p from EXACT integers (r_p <= 141 always holds for n <= 20000; asserted).  Compares with the session's
# camion_squarefree_list.txt and with HYP-9358(c) counts; lists the Chernick-type family 3b(3b-2) and the exotic ones.
from sympy import factorint, isprime, primerange
D = {}
for line in open("a1_exact_D.txt"):
    r, d = line.split(); D[int(r)] = int(d)
N = 20000
lst = []; maxr = 0
for n in range(4, N + 1):
    if isprime(n): continue
    f = factorint(n)
    if any(e > 1 for e in f.values()): continue
    ok = True
    for p in f:
        if p == 2: continue
        r = (n - 1) % (p - 1) + 1; maxr = max(maxr, r)
        assert r in D
        if D[r] % p: ok = False; break
    if ok: lst.append(n)
print("max r_p needed:", maxr)
print("squarefree Camion pseudoprimes <= 20000:", len(lst))
sess = [int(l.split()[0]) for l in open("../camion_squarefree_list.txt")]
print("identical to session list camion_squarefree_list.txt:", sess == lst)
two_p = [n for n in lst if n % 2 == 0 and isprime(n // 2) and n // 2 > 2]
three_p = [n for n in lst if n % 3 == 0 and n % 2 == 1 and isprime(n // 3) and n // 3 > 3]
exotic = [n for n in lst if n not in set(two_p) | set(three_p)]
print("2p:", len(two_p), " (pi(10000)-1 =", len(list(primerange(3, 10001))), ")")
print("3p (p odd prime != 3):", len(three_p), " (#primes 5..6666 =", len(list(primerange(5, N // 3 + 1))), ")")
print("exotic:", len(exotic), "odd", sum(1 for n in exotic if n % 2), "even", sum(1 for n in exotic if n % 2 == 0))
print("exotic odd:", [n for n in exotic if n % 2])
print("exotic even:", [n for n in exotic if n % 2 == 0])
carm = [n for n in range(3, N + 1, 2) if not isprime(n) and all(e == 1 for e in factorint(n).values())
        and all((n - 1) % (p - 1) == 0 for p in factorint(n))]
print("Carmichael numbers <= 20000:", carm, " all in list:", all(n in lst for n in carm))
fam = []
for b in primerange(5, 200):
    a = 3 * b - 2
    if isprime(a) and 3 * a * b <= N: fam.append(3 * a * b)
print("3b(3b-2) family <= 20000:", fam, " all in list:", all(n in lst for n in fam))
# the session's 'not of the form 3p' list (odd, squarefree):
print("odd squarefree passes not of form 3p, first 11:", [n for n in lst if n % 2 and n not in set(three_p)][:11])
