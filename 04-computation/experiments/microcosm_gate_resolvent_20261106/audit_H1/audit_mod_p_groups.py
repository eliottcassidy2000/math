# audit_H1: group-theory checks for note section 3: <x+2, 2x> mod 7 has order 21 and equals {x -> a x + b : a in QR7}; primes p < 200 with <2> = QR_p split by p mod 4; |{x -> a x + b, a in QR_p}| = p(p-1)/2 and these maps are automorphisms of QR_p.
from sympy import isprime, n_order

def generated_group(p, gens):
    # maps as pairs (a, b): x -> a x + b mod p
    ident = (1, 0)
    G = {ident}
    frontier = [ident]
    while frontier:
        new = []
        for (a, b) in frontier:
            for (c, d) in gens:
                h = ((c * a) % p, (c * b + d) % p)   # apply (a,b) then (c,d)
                if h not in G:
                    G.add(h)
                    new.append(h)
        frontier = new
    return G

G7 = generated_group(7, [(1, 2), (2, 0)])
QR7 = {(x * x) % 7 for x in range(1, 7)}
print("|<x+2, 2x> mod 7| =", len(G7), "; multipliers =", sorted({a for a, b in G7}), "; QR7 =", sorted(QR7))
# automorphisms of QR7 (brute force over S7 is 5040 perms)
import itertools
arcs = {(i, j) for i in range(7) for j in range(7) if (j - i) % 7 in QR7}
aut = [p for p in itertools.permutations(range(7)) if all((p[i], p[j]) in arcs for (i, j) in arcs)]
print("|Aut(QR7)| by brute force =", len(aut), "; equals <x+2,2x>:", {tuple((a * x + b) % 7 for x in range(7)) for a, b in G7} == set(aut))
# Fano lines: x + {1,2,4}
lines = [frozenset((x + d) % 7 for d in (1, 2, 4)) for x in range(7)]
ok = all(len(l1 & l2) == 1 for l1, l2 in itertools.combinations(lines, 2))
print("x + {1,2,4} pairwise meet in one point (Fano plane):", ok)
both = [p for p in range(3, 200) if isprime(p) and n_order(2, p) == (p - 1) // 2]
print("primes p < 200 with ord_p(2) = (p-1)/2:", both)
print("   of these p = 3 mod 4 (QR_p a tournament):", [p for p in both if p % 4 == 3])
print("   of these p = 1 mod 4 (QR_p a graph, not a tournament):", [p for p in both if p % 4 == 1])
p = 23
QR = {(x * x) % p for x in range(1, p)}
arcs = {(i, j) for i in range(p) for j in range(p) if (j - i) % p in QR}
affine = [(a, b) for a in range(1, p) for b in range(p) if all((((a * i + b) % p), ((a * j + b) % p)) in arcs for (i, j) in arcs)]
print("p = 23: affine automorphisms of QR23:", len(affine), "= 23*11 =", 23 * 11)
