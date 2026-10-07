"""Spot checks of the {2,3,11} catalogue and the -1/2 thread (all exact integer arithmetic)."""
from sympy import factorint, n_order, legendre_symbol, isprime, primerange
from sympy import sqrt, nsimplify, GoldenRatio, expand, simplify
print("2^11-3^7 =", 2**11 - 3**7, " 2^11+3^7 =", 2**11 + 3**7, factorint(2**11 + 3**7))
print("ord_11(2), ord_11(3):", n_order(2, 11), n_order(3, 11), " 3^5-1 =", 3**5 - 1, factorint(3**5 - 1))
print("Fermat quotients mod 11: q(2) =", ((2**10 - 1)//11) % 11, " q(3) =", ((3**10 - 1)//11) % 11)
print("ternary Golay sphere 1+2*11+4*C(11,2) =", 1 + 2*11 + 4*55, "= 3^5?", 1 + 2*11 + 4*55 == 3**5)
print("binary Golay sphere 1+23+C(23,2)+C(23,3) =", 1 + 23 + 253 + 1771, "= 2^11?", 1 + 23 + 253 + 1771 == 2**11)
print("Hamming sphere 1+7 =", 8)
print("base-3 repunits (3^j-1)/2:", [(j, (3**j - 1)//2, factorint((3**j - 1)//2)) for j in range(1, 9)])
print("1093 Wieferich base 2:", pow(2, 1092, 1093**2) == 1, " 11 Wieferich base 3:", pow(3, 10, 121) == 1)
print("-1/2 mod 3^5 =", (-pow(2, -1, 3**5)) % 3**5)
phi = GoldenRatio
print("phi^10 - 1 - 11 phi^5 =", simplify(expand(phi**10 - 1 - 11*phi**5)))
print("(4-phi)(3+phi) =", simplify(expand((4 - phi)*(3 + phi))))
print("phi mod 11 roots of x^2-x-1:", [x for x in range(11) if (x*x - x - 1) % 11 == 0], " orders:", [n_order(x, 11) for x in [4, 8]])
print("Legendre (3|11),(5|11),(-1|11):", legendre_symbol(3, 11), legendre_symbol(5, 11), legendre_symbol(10, 11))
# Borel regular on unordered pairs iff p = 3 mod 4
for p in [5, 7, 11, 13, 23]:
    sq = {pow(a, 2, p) for a in range(1, p)}
    B = [(a, b) for a in sq for b in range(p)]
    pairs = set(); ok = True
    for (a, b) in B:
        img = frozenset(((a*0 + b) % p, (a*1 + b) % p))
        if img in pairs: ok = False
        pairs.add(img)
    print("p=%d |Borel|=%d C(p,2)=%d regular on pairs (orbit of {0,1} free & full): %s ; <2>=squares? %s <3>=squares? %s" %
          (p, len(B), p*(p-1)//2, ok and len(pairs) == p*(p-1)//2, n_order(2, p) == (p-1)//2, n_order(3, p) == (p-1)//2))
# the -17 cycle relative to -1/2 and -1
cyc = [-17]
while True:
    x = cyc[-1]; y = x // 2 if x % 2 == 0 else (3*x + 1)//2
    if y == -17: break
    cyc.append(y)
print("-17 cycle:", cyc)
print("  z = 2x+1:", [2*x + 1 for x in cyc])
print("  y = x+1 :", [x + 1 for x in cyc])
runstarts = [cyc[i] for i in range(len(cyc)) if cyc[i] % 2 != 0 and cyc[i-1] % 2 == 0]
print("  odd-run starts:", runstarts, " 4|m|-1:", [4*abs(m) - 1 for m in runstarts])
print("  -5 cycle odd-run start -5: 4*5-1 =", 19)
lucky = [m for m in range(2, 100) if all(isprime(n*n + n + m) for n in range(0, m - 1))]
print("Euler lucky numbers (n^2+n+m prime for 0<=n<=m-2), m<100:", lucky)
# genus facts on the 101 discriminants
exec(open('idoneal.py').read().split("# (3) Borwein-Choi")[0].replace("DMAX = 30000", "DMAX = 8000"))
print("idoneal n>9 with n = 2 mod 3:", [n for n in idon if n > 9 and n % 3 == 2])
print("odd |D| = 7 mod 8 (2 split):", [N for N in exp2 if N % 8 == 7])
viol = []
for N in exp2:
    for p in primerange(3, 200):
        if p * p < N // 4 and N % p and legendre_symbol((-N) % p, p) == 1:
            viol.append((N, p))
print("exponent-2 |D| with an odd prime p, p^2 < |D|/4, split:", viol)
