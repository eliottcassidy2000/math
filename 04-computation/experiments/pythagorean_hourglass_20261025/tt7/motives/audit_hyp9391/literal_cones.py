# Audit of HYP-9391 (B): literal cone reductions, O(q^3), no Legendre table.
# For a seven-vertex graph H (order orientation: arc i -> j for i < j carries chi(x_j - x_i)) and an independent set C
# of degree-3 vertices, Y_H(C) = sum over injective x on the base B = V - C, with two base vertices fixed to 0 and 1,
# of (base arcs) * prod_{c in C} sum_{t in F_q} prod_{u in N(c)} (arc factor of the pair {c, u} at x_c = t).
# The fixed pair is the LAST two base vertices (the lane fixes the first two; Y_H is AGL(1)-invariant when |E| is
# even, so the value must not change).  Compared with the auditor's C closed forms (ysums.c).
# Usage: literal_cones.py qmax ; writes literal_cones.out
import json, sys, subprocess, os
import numpy as np
from sympy import primerange

HERE = os.path.dirname(os.path.abspath(__file__))
G = json.load(open(os.path.join(HERE, '..', 'genuine37.json')))
YS = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.environ.get('AUDIT_SCRATCH', '.'), 'ysums')

def chi_tab(q):
    c = np.full(q, -1, dtype=np.int64); c[0] = 0
    c[[(x * x) % q for x in range(1, q)]] = 1
    return c

def literal_Y(q, edges, C):
    chi = chi_tab(q)
    E = {tuple(sorted(e)) for e in edges}
    B = [v for v in range(7) if v not in C]
    nb = {c: sorted(u for u in range(7) if tuple(sorted((u, c))) in E) for c in C}
    for c in C:
        assert len(nb[c]) == 3 and all(u in B for u in nb[c])
    fa, fb = B[-2], B[-1]           # fixed to 0, 1
    free = B[:-2]
    t = np.arange(q)
    tot = 0
    for xu in range(q):
        for xv in range(q):
            x = {fa: 0, fb: 1, free[0]: xu, free[1]: xv}
            vals = [x[v] for v in B]
            if len(set(vals)) < 4: continue
            pr = 1
            for i, j in E:
                if i in B and j in B:
                    pr *= int(chi[(x[j] - x[i]) % q])
            if pr == 0: continue
            for c in C:
                f = np.ones(q, dtype=np.int64)
                for u in nb[c]:
                    f = f * (chi[(x[u] - t) % q] if c < u else chi[(t - x[u]) % q])
                pr *= int(f.sum())
                if pr == 0: break
            tot += pr
    return tot

cones = {2: (3, 5, 6), 4: (4, 5, 6), 5: (3, 5, 6), 6: (4, 5, 6), 9: (4, 5, 6), 18: (4, 5, 6)}
order = {2: 1, 4: 2, 5: 3, 6: 4, 9: 5, 18: 6}           # columns of ysums output
qmax = int(sys.argv[1]) if len(sys.argv) > 1 else 100
qs = [q for q in primerange(7, qmax) if q % 4 == 3]
cl = {}
for line in subprocess.run([YS] + [str(q) for q in qs], capture_output=True, text=True).stdout.split('\n'):
    if line.strip():
        r = list(map(int, line.split())); cl[r[0]] = r
lines = []
allok = True
for q in qs:
    lit = {i: literal_Y(q, G[i]['edges'], C) for i, C in cones.items()}
    ok = all(lit[i] == cl[q][order[i]] for i in cones)
    allok &= ok
    lines.append(f"q={q} literal " + ' '.join(f"Y{i:02d}={lit[i]}" for i in cones) + f"  closed-form match: {ok}")
lines.append(f"literal cone sums (fixed pair = last two base vertices) equal the closed forms at all {len(qs)} primes q = 3 mod 4 below {qmax}: {allok}")
lines.append("identities among literal cone sums: Y04 = -Y02 at all: " + str(all(int(l.split('Y04=')[1].split()[0]) == -int(l.split('Y02=')[1].split()[0]) for l in lines[:-1]))
             + "; Y06 = Y05 at all: " + str(all(int(l.split('Y06=')[1].split()[0]) == int(l.split('Y05=')[1].split()[0]) for l in lines[:-1])))
print('\n'.join(lines))
open(os.path.join(HERE, 'literal_cones.out'), 'w').write('\n'.join(lines) + '\n')
