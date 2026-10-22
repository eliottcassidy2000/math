"""Computational checks of the configurations used in the proof of the one-way lemma (session opus-2026-10-21-S50;
THM-4745, upgrade of HYP-9356).  Note: 05-knowledge/results/h_spectrum_multiplicities_20261021.md.

The proof needs no computation beyond exact_gluings.py (single blocks NT_3, NT_4, NT_5 at end positions).  These
checks confirm, by exact Hamiltonian-path counts, that the configurations which the proof excludes by explicit
witnesses are indeed never exact, and that the witnesses are valid:
  (C1) the junction configuration F(s, r): A_left = NT_s ending at v, A_right = NT_r starting at v, B = C_3 at v, and
       cross arcs A_left' -> B' -> A_right', for s, r = 3..6: H(T) > 3 H(A); and the path
       L_X b1 R_X v L_Y b2 R_Y (L_X v L_Y a back-arc path of A_left, R_X v R_Y one of A_right) is a Hamiltonian
       path of T whose restriction to B is not a path;
  (C2) the mirrored configuration F'(3,3) (A_right' -> B' -> A_left') and its 7-vertex witness;
  (C3) the interior configurations: A = NT_s with v interior, cross arcs forced by (F1), (B1) (a -> b iff a beats v),
       B = NT_r (r = 3..5) with v at its start or end: never exact, and the pair of paths of Step 2 has two merges.
Prints ALL CHECKS PASSED.  Runtime a few seconds."""
import sys, os, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from h_multiplicities import ham_one                                     # noqa: E402

FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def H(adj):
    return int(ham_one(len(adj), np.array(adj, dtype=np.int64)))


def nt(s):
    adj = [0] * s
    for i in range(s):
        for j in range(i + 1, s):
            adj[i] |= 1 << j
    adj[0] &= ~(1 << (s - 1)); adj[s - 1] |= 1
    return adj


def arcs_of(adj, mp):
    return [(mp[u], mp[w]) for u in range(len(adj)) for w in range(len(adj)) if (adj[u] >> w) & 1]


def build(n, arcs):
    adj = [0] * n
    for u, w in arcs:
        adj[u] |= 1 << w
    return adj


def is_path(adj, seq):
    return all((adj[seq[i]] >> seq[i + 1]) & 1 for i in range(len(seq) - 1))


# (C1)
okH = True; okW = True
for s in range(3, 7):
    for r in range(3, 7):
        v = s - 1; L = list(range(s - 1)); R = list(range(s, s + r - 1)); b1, b2 = s + r - 1, s + r
        n = s + r + 1
        arcs = arcs_of(nt(s), {i: i for i in range(s)})
        mR = {0: v}; mR.update({i: s + i - 1 for i in range(1, r)})
        arcs += arcs_of(nt(r), mR)
        arcs += [(x, z) for x in L for z in R] + [(v, b1), (b1, b2), (b2, v)]
        arcs += [(x, b) for x in L for b in (b1, b2)] + [(b, z) for b in (b1, b2) for z in R]
        T = build(n, arcs)
        okH &= H(T) > 3 * H(nt(s)) * H(nt(r))
        # witness: A_left back-arc path (interior X'' , v, x_0, rest) with X'' = all interior -> L_X = 1..s-2, L_Y = (0)
        LX = list(range(1, s - 1)); LY = [0]
        # A_right back-arc path (X', z_e, v, Y') with X' = all interior: R_X = interior..., z_e; R_Y = ()
        RX = [s + i - 1 for i in range(1, r)]; RY = []
        W = LX + [b1] + RX + [v] + LY + [b2] + RY
        okW &= sorted(W) == list(range(n)) and is_path(T, W) and not is_path(T, [b1, v, b2])
check(okH, '(C1) F(s,r), s,r = 3..6: H(T) > 3 H(A) (never exact)')
check(okW, '(C1) the witness L_X b1 R_X v L_Y b2 R_Y is a Hamiltonian path of T restricting to b1 v b2 (not a path of B)')

# (C2)
arcs = arcs_of(nt(3), {0: 0, 1: 1, 2: 2}) + arcs_of(nt(3), {0: 2, 1: 3, 2: 4})
arcs += [(x, z) for x in (0, 1) for z in (3, 4)] + [(2, 5), (5, 6), (6, 2)]
arcs += [(z, b) for z in (3, 4) for b in (5, 6)] + [(b, x) for b in (5, 6) for x in (0, 1)]
T = build(7, arcs)
W = [5, 0, 1, 2, 3, 4, 6]
check(H(T) > 27 and is_path(T, W) and not is_path(T, [5, 2, 6]), "(C2) F'(3,3): not exact; witness b1 x0 x1 v z1 z2 b2 restricts to b1 v b2")

# (C3)
ok_ne = True; ok_two = True
for s in range(4, 7):
    for t in range(1, s - 1):
        for r in range(3, 6):
            for end in (0, r - 1):
                n = s + r - 1
                arcs = arcs_of(nt(s), {i: i for i in range(s)})
                mB = {}; k = s
                for y in range(r):
                    if y == end:
                        mB[y] = t
                    else:
                        mB[y] = k; k += 1
                arcs += arcs_of(nt(r), mB)
                Ap = [i for i in range(s) if i != t]; Bp = list(range(s, n))
                arcs += [((a, b) if a < t else (b, a)) for a in Ap for b in Bp]
                T = build(n, arcs)
                ok_ne &= H(T) != H(nt(s)) * H(nt(r))
                # the pair: P = x_1..x_t x_N x_0 x_{t+1}..x_{N-1}, Q = (Q1, v, b): count Hamiltonian paths of T restricting to it
                N = s - 1
                P = list(range(1, t + 1)) + [N, 0] + list(range(t + 1, N))
                # a B-path with exactly one vertex after v
                Bpaths = [q for q in itertools.permutations([t] + Bp) if is_path(T, list(q))]
                Qs = [q for q in Bpaths if q.index(t) == len(q) - 2]
                cnt = 0
                Aset = set(range(s)); Bset = set([t] + Bp)
                Q = Qs[0]
                for perm in itertools.permutations(range(n)):
                    if is_path(T, list(perm)):
                        if [x for x in perm if x in Aset] == P and [x for x in perm if x in Bset] == list(Q):
                            cnt += 1
                ok_two &= cnt >= 2
check(ok_ne, '(C3) interior positions with the forced cross arcs: never exact (s = 4..6, all interior t, r = 3..5, both ends)')
check(ok_two, '(C3) the pair (x_1..x_t x_N x_0 x_(t+1)..x_(N-1), Q_1 v b) has at least two merges')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
