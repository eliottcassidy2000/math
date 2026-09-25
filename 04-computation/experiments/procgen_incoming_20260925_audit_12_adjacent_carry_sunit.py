"""Seed pre-screen: THM-4469 adjacent carry pairs (R_B' = R_B + 1, 3^a > 2^(L+1)) as
{2,3}-unit equations sum_i 3^(a-i) (2^(j'_i) - 2^(j_i)) = 1 (exchange law is q-independent).
Recompute the census for L <= 24 and record, for each pair, the number of non-cancelled
terms (positions i with j_i != j'_i) and whether a proper subsum vanishes."""
from itertools import combinations
import sys
Lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 24
tot = {}
for L in range(1, Lmax + 1):
    for a in range(1, L + 1):
        if not (3 ** a > 2 ** (L + 1)):
            continue
        vals = {}
        for pos in combinations(range(L), a):
            R = sum(2 ** j * 3 ** (a - 1 - i) for i, j in enumerate(pos))
            vals[R] = pos
        pairs = []
        for R, pos in vals.items():
            if R + 1 in vals:
                pos2 = vals[R + 1]
                terms = [3 ** (a - 1 - i) * (2 ** pos2[i] - 2 ** pos[i]) for i in range(a)]
                nz = [t for t in terms if t != 0]
                # vanishing proper subsums among nonzero terms (small check)
                van = False
                n = len(nz)
                if n <= 16:
                    for r in range(1, n):
                        for S in combinations(range(n), r):
                            if sum(nz[s] for s in S) == 0:
                                van = True; break
                        if van: break
                pairs.append((len(nz), van))
        if pairs:
            tot[(L, a)] = pairs
for k, v in sorted(tot.items()):
    print("(L,a)=%s adjacent pairs=%d  nonzero-term counts=%s  degenerate(vanishing subsum)=%s"
          % (k, len(v), sorted(x[0] for x in v), sum(1 for x in v if x[1])))
