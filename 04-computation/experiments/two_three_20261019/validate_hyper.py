"""Brute-force 3-hypertournament counts for n = 4, 5 (orientations of all triples modulo S_n) to validate the Burnside count."""
from itertools import combinations, permutations, product
def count(n):
    triples = list(combinations(range(n), 3)); idx = {t: i for i, t in enumerate(triples)}
    perms = list(permutations(range(n)))
    def canon(orient):
        best = None
        for p in perms:
            img = [0]*len(triples)
            for t, o in zip(triples, orient):
                im = [p[v] for v in t]; s = tuple(sorted(im))
                # cyclic orientation: o = 0 means the cyclic order (t0 t1 t2), 1 means (t0 t2 t1); image orientation = o xor sign of the permutation sorting im
                inv = sum(1 for a in range(3) for b in range(a+1, 3) if im[a] > im[b])
                img[idx[s]] = o ^ (inv % 2)
            tup = tuple(img)
            if best is None or tup < best: best = tup
        return best
    return len({canon(o) for o in product((0, 1), repeat=len(triples))})
for n in (4, 5): print(f"n={n}: brute force 3-hypertournament classes = {count(n)}")
import json, math
seqs = json.load(open('k_tournaments.json'))
print("leading rigidity corrections: even k ~ 2 C(n,3) k^-(2n-4); odd k ~ C(n,2) k^-(n-1):")
for k in (3, 4, 5, 6):
    for n in (10, 12):
        T = seqs[str(k)][n-1]; ratio = T*math.factorial(n)/k**(n*(n-1)//2) - 1
        pred = (math.comb(n,2)*k**(-(n-1))) if k % 2 else (2*math.comb(n,3)*k**(-(2*n-4)))
        print(f"  k={k}, n={n}: ratio-1 = {ratio:.3e}, predicted leading term {pred:.3e}")
