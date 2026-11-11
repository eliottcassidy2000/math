"""domination_ray.py - the domination ray of THM-5130 (H2) on every pair (B, O) over the tournament classes of order <= 6
(kind-pasteur-2026-11-10-S22; added after the audit, which asked for the programme behind "3,978 pairs").

Needs the class lists cls_1.txt ... cls_6.txt (first field of each line = upper-triangle code), produced by
classgen.c of 05-knowledge/results/first_creators_genealogy_20261107/ :  ./classgen 6 cls
usage: python domination_ray.py DIR_WITH_CLASS_LISTS

For every class B of order <= 6 and every nonempty O:
  M = U_B(O) is compared with delta_dom = sum over nonempty S in O of 2^(number of vertices outside O beaten by S);
for orders <= 4 the dominating sets of the graph G_t of the closed tournament are counted literally for t = 1, 2, 3 and
compared with rho + 2^t delta_dom and with the inequality H(T_t) - d(G_t) >= H(B) - 1 of Session261."""
import itertools, collections, os, sys
from towers import state


def adj(m, code):
    A = [[0] * m for _ in range(m)]
    b = 0
    for i in range(m):
        for j in range(i + 1, m):
            if code >> b & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            b += 1
    return A


def dom_count(n, edges):
    nb = [1 << v for v in range(n)]
    for a, b in edges:
        nb[a] |= 1 << b
        nb[b] |= 1 << a
    cnt = 0
    for S in range(1 << n):
        cov = 0
        for v in range(n):
            if S >> v & 1:
                cov |= nb[v]
        cnt += (cov == (1 << n) - 1)
    return cnt


d = sys.argv[1] if len(sys.argv) > 1 else "."
bad_ray = bad_ineq = bad_slope = 0
npairs = 0
eq = collections.Counter()
for m in range(1, 7):
    for l in open(os.path.join(d, "cls_%d.txt" % m)):
        code = int(l.split()[0])
        A = adj(m, code)
        for Om in range(1, 1 << m):
            O = [i for i in range(m) if Om >> i & 1]
            I = [i for i in range(m) if not Om >> i & 1]
            r, M, h, E = state(A, O)
            npairs += 1
            dl = 0
            for k in range(1, len(O) + 1):
                for S in itertools.combinations(O, k):
                    dl += 2 ** sum(1 for i in I if any(A[s][i] for s in S))
            bad_slope += (M < dl)
            eq[(M == dl, h == 1)] += 1
            if m <= 4:
                ds = []
                for t in (1, 2, 3):
                    idx = {v: k for k, v in enumerate(O)}
                    ridx = {v: len(O) + k for k, v in enumerate(I)}
                    edges = [(idx[o], ridx[i]) for o in O for i in I if A[o][i]]
                    edges += [(idx[o], len(O) + len(I) + z) for o in O for z in range(t)]
                    ds.append(dom_count(len(O) + len(I) + t, edges))
                rho = ds[0] - 2 * dl
                bad_ray += not all(ds[t - 1] == rho + 2 ** t * dl for t in (1, 2, 3))
                bad_ineq += sum(1 for t in (1, 2, 3) if (r + 2 ** t * M) - ds[t - 1] < h - 1)
print("pairs (B, O) over classes of order <= 6, O nonempty:", npairs)
print("slope below the domination slope:", bad_slope)
print("ray form d(G_t) = rho + 2^t delta_dom fails (orders <= 4, t <= 3):", bad_ray,
      "| inequality H(T_t) - d(G_t) >= H(B) - 1 fails:", bad_ineq)
print("pairs by (M == delta_dom, B transitive):", dict(eq))
