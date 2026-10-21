"""How mixed gluings of small chains fail to be exact (session opus-2026-10-21-S50; THM-4745, remaining step of
HYP-9356).  Note: 05-knowledge/results/h_spectrum_multiplicities_20261021.md.

For A, B in {C_3 = IC(2), R_4 = IC(3)}, every gluing vertex (i in A, j in B) and every orientation of the cross arcs,
list every Hamiltonian path of T (by permutations) and split the surplus H(T) - H(A)H(B) into
  * R-violations: paths whose restriction to A or to B is not a path;
  * extra merges: several paths restricting to the same pair.
Reports, per position, the exact gluings, the least surplus among mixed patterns and the failure types.
Checks: every exact gluing is one-way; every mixed pattern has positive surplus; at the end-to-start positions
(v last in A and first in B, or a 3-cycle factor) with at least one R_4 factor every mixed pattern has an R-violation,
while for C_3 u C_3 four of the fourteen mixed patterns fail only through extra merges (so a proof of the one-way lemma
must use both conditions (R) and (U) of THM-4745 (1)).
Prints ALL CHECKS PASSED.  Runtime about 1 minute."""
import itertools, time
from collections import Counter

FAIL = []
T0 = time.time()


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def chain(c):
    n = 1 + sum(c)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = True
    l = 0
    for p in c:
        r = l + p
        adj[l][r] = False; adj[r][l] = True
        l = r
    return adj


def ham_paths(arc, verts):
    return [p for p in itertools.permutations(verts) if all(arc(p[k], p[k + 1]) for k in range(len(p) - 1))]


def analyse(cA, cB, i, j):
    A = chain(cA); B = chain(cB)
    a, b = len(A), len(B)
    Amap = {}; k = 0
    for x in range(a):
        Amap[x] = a - 1 if x == i else k
        if x != i:
            k += 1
    Bmap = {}; k = a
    for y in range(b):
        Bmap[y] = a - 1 if y == j else k
        if y != j:
            k += 1
    n = a + b - 1
    base = {}
    for x in range(a):
        for x2 in range(a):
            if A[x][x2]:
                base[(Amap[x], Amap[x2])] = True
    for y in range(b):
        for y2 in range(b):
            if B[y][y2]:
                base[(Bmap[y], Bmap[y2])] = True
    v = a - 1
    Avert = list(range(a - 1)) + [v]; Bvert = list(range(a, n)) + [v]
    res = []
    full = (1 << ((a - 1) * (b - 1))) - 1
    for mask in range(full + 1):
        arcs = dict(base)
        for x in range(a - 1):
            for y in range(b - 1):
                if (mask >> (x * (b - 1) + y)) & 1:
                    arcs[(x, a + y)] = True
                else:
                    arcs[(a + y, x)] = True
        arc = lambda u, w: arcs.get((u, w), False)
        HT = ham_paths(arc, list(range(n)))
        nA = len(ham_paths(arc, Avert)); nB = len(ham_paths(arc, Bvert))
        viol = 0; pairs = Counter()
        for P in HT:
            pa = tuple(x for x in P if x in Avert); pb = tuple(x for x in P if x in Bvert)
            if all(arc(pa[t], pa[t + 1]) for t in range(len(pa) - 1)) and all(arc(pb[t], pb[t + 1]) for t in range(len(pb) - 1)):
                pairs[(pa, pb)] += 1
            else:
                viol += 1
        res.append((mask in (0, full), len(HT) - nA * nB, viol, sum(c - 1 for c in pairs.values()), mask == full))
    return res


ok_exact_oneway = True; ok_pos = True; ok_end = True; c33 = []
for cA, cB in (((2,), (2,)), ((2,), (3,)), ((3,), (2,)), ((3,), (3,))):
    a = 1 + sum(cA); b = 1 + sum(cB)
    for i in range(a):
        for j in range(b):
            res = analyse(cA, cB, i, j)
            mixed = [r for r in res if not r[0]]
            exact = [r for r in res if r[1] == 0]
            ok_exact_oneway &= all(r[0] for r in exact)
            ok_pos &= all(r[1] > 0 for r in mixed)
            end_pos = ((i == a - 1 or cA == (2,)) and (j == 0 or cB == (2,))) or ((j == b - 1 or cB == (2,)) and (i == 0 or cA == (2,)))
            if end_pos and (cA, cB) != ((2,), (2,)):
                ok_end &= all(r[2] > 0 for r in mixed)
            if (cA, cB) == ((2,), (2,)):
                c33.append(sum(1 for r in mixed if r[3] and not r[2]))
            print(f'     A = {cA}, i = {i}; B = {cB}, j = {j}: exact {len(exact)}; mixed {len(mixed)}, least surplus '
                  f'{min(r[1] for r in mixed)}; R-violations only {sum(1 for r in mixed if r[2] and not r[3])}, '
                  f'extra merges only {sum(1 for r in mixed if r[3] and not r[2])}, both {sum(1 for r in mixed if r[2] and r[3])}'
                  + ('  [end-to-start position]' if end_pos else ''), flush=True)
check(ok_exact_oneway, 'every exact gluing among these is one-way')
check(ok_pos, 'every mixed pattern has positive surplus')
check(ok_end, 'at the end-to-start positions with an R_4 factor every mixed pattern has a Hamiltonian path of T violating (R)')
check(c33 and all(c == 4 for c in c33), 'C_3 u C_3: at every gluing position exactly 4 of the 14 mixed patterns fail only through extra merges (no R-violation)')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
