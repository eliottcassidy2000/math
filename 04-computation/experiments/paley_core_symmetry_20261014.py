"""Odd-reversal automorphisms of the Paley run cores (THM-4652 part 4): a core whose undirected graph has an automorphism
reversing an odd number of the low->high oriented edges has pattern sum zero for q = 3 mod 4.
Session opus-2026-10-14-S35.
"""
import itertools
CORES = {
 'a_K6-(P3+K2)': [(0,1),(0,2),(0,3),(0,4),(0,5),(1,2),(1,3),(1,4),(2,3),(2,4),(3,5),(4,5)],
 'b_octahedron': [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,5),(2,4),(2,5),(3,4),(3,5),(4,5)],
 'c_K24+2K2':    [(0,1),(0,2),(0,3),(0,4),(1,2),(1,5),(2,5),(3,4),(3,5),(4,5)],
 'd_K33+e':      [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,5),(3,5),(4,5)],
 'W4':           [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(2,4),(3,4)],
 'K6-P6':        [(0,1),(0,3),(0,4),(0,5),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4)],
}
for name, E in CORES.items():
    n = 1 + max(max(e) for e in E)
    S = set(E)
    odd = None
    nauto = 0
    for perm in itertools.permutations(range(n)):
        img = [(perm[u], perm[v]) for u, v in E]
        und = set(tuple(sorted(e)) for e in img)
        if und != S:
            continue
        nauto += 1
        rev = sum(1 for (a, b) in img if a > b)
        if rev % 2 == 1:
            odd = perm
            break
    print(f'{name}: automorphisms scanned {nauto}; odd-reversal automorphism: {odd}  ->  {"VANISHES (THM-4652 part 4)" if odd else "no symmetry obstruction"}')
