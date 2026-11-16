#!/usr/bin/env python3
"""leech_numerics.py: exact arithmetic for the comparison 'Riesel wheel vs {4,2,5,7,7,13} vs Leech' (klein-2026-11-16-two3)."""
from sympy import factorint, n_order, divisors
from math import prod
riesel = [3, 5, 7, 13, 17, 241]
ords = [n_order(2, p) for p in riesel]
phases = [(-pow(2, 0, p)) for p in riesel]  # placeholder
# covering phases of 509203*2^e - 1: e mod ord with p | 509203*2^e - 1
ph = []
for p, d in zip(riesel, ords):
    ph.append([e for e in range(d) if (509203 * pow(2, e, p) - 1) % p == 0])
print("Riesel wheel primes:", riesel, "orders:", ords, "covering phases of 509203*2^e-1:", ph)
print("phase multiset (as in CLOCKS.md):", [x[0] for x in ph], "sum", sum(x[0] for x in ph), "; order sum", sum(ords), "; prime sum", sum(riesel), "; product/(2^24-1):", prod(riesel), (2**24 - 1) // 3)
S_owner_old = [4, 2, 5, 3, 7, 13]; S_owner_new = [4, 2, 5, 7, 7, 13]
print("owner sets: old {4,2,5,3,7,13} sum", sum(S_owner_old), "product", prod(S_owner_old), "; new {4,2,5,7,7,13} sum", sum(S_owner_new), "product", prod(S_owner_new))
# pure-A Niemeier lattices A_n^k: n k = 24, Coxeter number h = n + 1
pureA = [(n, 24 // n, n + 1) for n in divisors(24)]
print("pure-A Niemeier (n, k, h):", pureA, "; Coxeter numbers", [h for _, _, h in pureA])
print("h with h-1 | 12:", [h for _, _, h in pureA if 12 % (h - 1) == 0])
# which Niemeier Coxeter numbers coincide with Riesel orders / primes / phases
print("Riesel orders that are Niemeier Coxeter numbers:", [d for d in ords if d in [h for _, _, h in pureA]])
print("Riesel primes that are Niemeier Coxeter numbers:", [p for p in riesel if p in [h for _, _, h in pureA]])
# Leech / Golay / Conway orders
leech_min = 196560; co0 = 2**22 * 3**9 * 5**4 * 7**2 * 11 * 13 * 23; m24 = 244823040
print("196560 =", factorint(leech_min), "; |Co_0| =", factorint(co0), "; |M24| =", factorint(m24))
print("Riesel primes dividing 196560:", [p for p in riesel if leech_min % p == 0], "; dividing |Co_0|:", [p for p in riesel if co0 % p == 0], "; dividing |M24|:", [p for p in riesel if m24 % p == 0])
# the user's multiset {4,2,5,7,7,13}: sum 38 = max Z/H; as residues?  as orders of 2? (ord_p(2)=4: p=5; 2: p=3; 5: p=31; 7: p=127; 13: p=8191)
print("primes with ord_p(2) in {2,4,5,7,13}:", {d: [p for p in factorint(2**d - 1) if n_order(2, p) == d] for d in (2, 4, 5, 7, 13)})
print("2^d - 1 for d in {2,4,5,7,13}:", {d: 2**d - 1 for d in (2, 4, 5, 7, 13)}, "(Mersenne primes at 2,5,7,13: 3, 31, 127, 8191)")
# density sum of the Riesel wheel and the 'twelve' wheel
print("sum 1/ord Riesel:", sum(1 / d for d in ords), "; divisors of 24:", divisors(24), "; sum 1/d over d|24, d>1:", sum(1 / d for d in divisors(24) if d > 1))
