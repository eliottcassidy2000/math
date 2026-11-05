"""aud_seeds.py -- audit S11: the realisable SEAM states (H, S, R) of rooted tournaments of order <= 5, by permutation
brute force (R = number of (permutation, cut) pairs whose prefix is a path ending at the root and whose suffix is a path).
Prints the set of (H, R) with S = 0 (claimed: the twelve plane seeds of plane_bfs.c) and writes aud_plane_seeds.txt."""
import itertools

def tours(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for code in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def is_path(A, seq):
    return all(A[seq[t]][seq[t + 1]] for t in range(len(seq) - 1))

def state(A, n, a):
    H = S = R = 0
    for p in itertools.permutations(range(n)):
        if is_path(A, p):
            H += 1
            S += p[0] == a
        for c in range(1, n + 1):
            if p[c - 1] == a and is_path(A, p[:c]) and is_path(A, p[c:]):
                R += 1
    return H, S, R

allstates = set(); plane = set()
for n in range(1, 6):
    for A in tours(n):
        for a in range(n):
            v = state(A, n, a)
            allstates.add(v)
            if v[1] == 0: plane.add((v[0], v[2]))
print('realisable states (H,S,R) of order <= 5:', len(allstates))
print('S = 0 states (H, R) of order <= 5:', sorted(plane), len(plane))
claimed = {(1, 2), (1, 4), (1, 8), (1, 16), (3, 6), (3, 12), (3, 24), (5, 10), (5, 12), (5, 14), (5, 18), (5, 32)}
print('equal to the twelve seeds of plane_bfs.c:', plane == claimed)
print('3D seeds realisable:', all(s in allstates for s in [(1, 1, 1), (1, 0, 2), (3, 1, 3), (1, 0, 4), (1, 0, 8), (1, 0, 16)]))
print('R == S mod 2 on all of them:', all((r - s) % 2 == 0 for h, s, r in allstates))
with open('aud_plane_seeds.txt', 'w') as o:
    for h, r in sorted(plane): o.write(f'{h} {r}\n')
