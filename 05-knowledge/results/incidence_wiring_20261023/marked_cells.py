"""Source/endpoint refinement on incoming DRTs; exact cyclic-three-cell decoder."""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import argparse, json
import check as core

HERE = Path(__file__).resolve().parent
need = core.need

def paley(n):
    S = {i*i%n for i in range(1, n)}
    return [[int((j-i)%n in S) for j in range(n)] for i in range(n)]

def double(A):
    n = len(A); B = [[0]*(2*n+1) for _ in range(2*n+1)]
    for i in range(n):
        B[n][i] = 1; B[n+1+i][n] = 1
        for j in range(n):
            B[i][j] = A[i][j]; B[n+1+i][n+1+j] = A[j][i]
            B[i][n+1+j] = int(i == j or A[i][j]); B[n+1+j][i] = 1-B[i][n+1+j]
    return B

def refine(A, marks):
    n = len(A); col = [marks.index(v)+1 if v in marks else 0 for v in range(n)]
    history = []
    for _ in range(n):
        history.append(sorted(Counter(col).values()))
        if len(set(col)) == n: return history
        rows = [(col[u], tuple(sorted(col[v] for v in range(n) if A[u][v]))) for u in range(n)]
        ids = {r: i for i, r in enumerate(sorted(set(rows)))}; nxt = [ids[r] for r in rows]
        if len(set(nxt)) == len(set(col)): return history
        col = nxt
    raise ValueError('refinement failed to stabilize')

def triangle_decode(A, D):
    n = len(A); D = list(D); outside = [v for v in range(n) if v not in D]
    need(all(sum(A[u][v] for v in D) == 1 for u in D), 'core is not a cyclic triangle')
    margins = {v: sum(A[u][v] for u in D) for v in outside}
    unit = {v: tuple(A[u][v] if margins[v] == 1 else 1-A[u][v] for u in D)
            for v in outside if margins[v] in (1, 2)}
    # Only the equality matrix is supplied to this clustering step.
    equal = {(v, w): sum(a*b for a, b in zip(unit[v], unit[w])) for v in unit for w in unit}
    groups = []; labels = {}
    for v in unit:
        hit = next((i for i, r in enumerate(groups) if equal[v, r]), None)
        if hit is None: hit = len(groups); groups.append(v)
        labels[v] = hit
    need(len(groups) <= 3, 'more than three core-row labels')
    B = [[0]*n for _ in range(n)]
    for u in outside:
        for v in outside: B[u][v] = A[u][v]
    for v in outside:
        for i, u in enumerate(D):
            if margins[v] in (0, 3): B[u][v] = margins[v]//3
            elif margins[v] == 1: B[u][v] = int(i == labels[v])
            else: B[u][v] = int(i != labels[v])
            B[v][u] = 1-B[u][v]
    chirality = 1
    if len(groups) >= 2:
        a = D[unit[groups[0]].index(1)]; b = D[unit[groups[1]].index(1)]
        chirality = A[a][b]
    for i, j in ((0, 1), (1, 2), (2, 0)):
        B[D[i]][D[j]] = chirality; B[D[j]][D[i]] = 1-chirality
    matches = 0
    for p in permutations(D):
        perm = list(range(n))
        for i, u in enumerate(D): perm[u] = p[i]
        matches += all(A[i][j] == B[perm[i]][perm[j]] for i in range(n) for j in range(n))
    need(matches > 0, 'cyclic core decoder changed a boundary attachment')
    return len(groups)

def main():
    p = argparse.ArgumentParser(); p.add_argument('--write', action='store_true'); args = p.parse_args()
    counts = Counter()
    for n in range(3, 6):
        for mask in range(1 << (n*(n-1)//2)):
            A = [[0]*n for _ in range(n)]
            for k, (i, j) in enumerate(combinations(range(n), 2)):
                A[i][j] = mask>>k&1; A[j][i] = 1-A[i][j]
            for D in combinations(range(n), 3):
                if all(sum(A[u][v] for v in D) == 1 for u in D): counts[triangle_decode(A, D)] += 1
    # Equality tests on a tree no longer suffice with three possible row labels.
    x, y = (0, 1, 0), (0, 1, 2)
    need((x[0] == x[1], x[1] == x[2]) == (y[0] == y[1], y[1] == y[2]) and (x[0] == x[2]) != (y[0] == y[2]), 'ternary tree hostile')
    # Two named boundary columns pin two distinct rows; their common-source
    # overlaps do not remember the orientation of the internal triangle.
    A = [[0]*5 for _ in range(5)]
    for i, j in combinations(range(5), 2): A[i][j] = 1
    A[0][2], A[2][0] = 0, 1
    for u in range(3):
        for v, selected in ((3, 0), (4, 1)):
            A[u][v] = int(u == selected); A[v][u] = 1-A[u][v]
    B = [r[:] for r in A]
    for i, j in combinations(range(3), 2): B[i][j], B[j][i] = A[j][i], A[i][j]
    need(all(A[u][v] == B[u][v] for u in range(5) for v in range(5) if u >= 3 or v >= 3), 'chirality hostile changes boundary')
    need(not any(all(A[i][j] == B[p[i]][p[j]] for i in range(5) for j in range(5)) for q in permutations(range(3)) for p in [q+(3, 4)]), 'chirality hostile has boundary-fixed isomorphism')
    chirality_rows = [[sum(x << j for j, x in enumerate(r)) for r in M] for M in (A, B)]
    drts = {}
    for name, A in (('P31', paley(31)), ('double_double_P7', double(double(paley(7))))):
        need(all(sum(r) == 15 for r in A), 'DRT score')
        need(all(sum(A[u][w]*A[v][w] for w in range(31)) == 7 for u, v in combinations(range(31), 2)), 'DRT pair incidence')
        one = Counter(tuple(refine(A, (u,))[-1]) for u in range(31))
        two = Counter(tuple(refine(A, (u, v))[-1]) for u in range(31) for v in range(31) if u != v)
        need(one == Counter({(1, 15, 15): 31}), 'one-mark DRT rank-three boundary')
        drts[name] = {'single_mark': [[list(k), v] for k, v in sorted(one.items())],
                      'ordered_two_marks': [[list(k), v] for k, v in sorted(two.items())],
                      'example_marks_0_1': refine(A, (0, 1))}
    out = {'status': 'FINITE-EXACT', 'triangle_decodings_by_used_core_rows': dict(counts),
           'chirality_hostile_outmasks_boundary_3_4_fixed': chirality_rows,
           'order31_all_ordered_vertex_pairs_each': 930, 'DRT_mark_tests': drts, 'checks': core.CHECKS}
    out = json.loads(json.dumps(out)); dest = HERE/'marked_cells.json'
    if args.write: dest.write_text(json.dumps(out, indent=2)+'\n')
    else: need(json.loads(dest.read_text()) == out, 'marked-cell receipt differs')
    print('triangle decodings', sum(counts.values()), 'DRT ordered mark pairs', 1860, 'checks', core.CHECKS)

if __name__ == '__main__': main()
