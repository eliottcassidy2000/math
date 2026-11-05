"""aud_F_proof_steps.py -- audit S11 of the proposed addendum (F) to THM-4879: the two EXACT counting identities used in
its proof, checked by brute force (explicit path enumeration) on every labelled tournament G of order 2..5 and every
ordered pair z != u, with T = G + w, w -> u only, every other vertex of G -> w:
    (i)  S_z(T)  = s_z(G) - P_(z->u)(G) + L_(z,u)(G)
    (ii) pi_z(T) = pi_z(G) + sum over paths C of G from u to z of (1 + Pi(G - V(C)))
where pi_z counts nonempty paths (any vertex set) ending at z and Pi counts all nonempty paths.
Also checks the reductions R_z == pi_z (mod 2) on T and G, and the final (F2) congruences."""
import itertools

def tours(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for code in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def paths(A, verts):
    """all nonempty directed paths using only vertices in verts (any vertex subset), as tuples"""
    out = []
    vs = list(verts)
    for r in range(1, len(vs) + 1):
        for sub in itertools.combinations(vs, r):
            for p in itertools.permutations(sub):
                if all(A[p[t]][p[t + 1]] for t in range(len(p) - 1)): out.append(p)
    return out

def ham(A, verts):
    vs = list(verts)
    return [p for p in itertools.permutations(vs) if all(A[p[t]][p[t + 1]] for t in range(len(p) - 1))] if vs else [()]

def two_covers(A, n, first_pred, second_pred):
    """ordered pairs (P, Q) of nonempty disjoint paths covering 0..n-1 with predicates on P and Q"""
    cnt = 0
    for mask in range(1, (1 << n) - 1):
        U = [v for v in range(n) if mask >> v & 1]; W = [v for v in range(n) if not mask >> v & 1]
        cnt += sum(1 for p in ham(A, U) if first_pred(p)) * sum(1 for q in ham(A, W) if second_pred(q))
    return cnt

def Rz(A, n, z):
    tot = 0
    for mask in range(1, 1 << n):
        if not mask >> z & 1: continue
        U = [v for v in range(n) if mask >> v & 1]; W = [v for v in range(n) if not mask >> v & 1]
        tot += sum(1 for p in ham(A, U) if p[-1] == z) * len(ham(A, W))
    return tot

bad = {k: 0 for k in ('i', 'ii', 'redT', 'redG', 'F2a', 'F2b')}; rows = 0
for n in range(2, 6):
    for G in tours(n):
        allp = paths(G, range(n))
        for z in range(n):
            for u in range(n):
                if z == u: continue
                rows += 1
                N = n + 1; w = n
                T = [row[:] + [0] for row in G] + [[0] * N]
                for x in range(n):
                    if x == u: T[w][x] = 1
                    else: T[x][w] = 1
                HT = ham(T, range(N)); HG = ham(G, range(n))
                SzT = sum(1 for p in HT if p[0] == z)
                sz = sum(1 for p in HG if p[0] == z); Pzu = sum(1 for p in HG if p[0] == z and p[-1] == u)
                Puz = sum(1 for p in HG if p[0] == u and p[-1] == z)
                L = two_covers(G, n, lambda p: p[0] == z, lambda q: q[0] == u)
                M = two_covers(G, n, lambda p: p[-1] == z, lambda q: q[-1] == u)
                if SzT != sz - Pzu + L: bad['i'] += 1
                piT = sum(1 for p in paths(T, range(N)) if p[-1] == z)
                piG = sum(1 for p in allp if p[-1] == z)
                rhs = piG
                for C in allp:
                    if C[0] == u and C[-1] == z:
                        rest = [v for v in range(n) if v not in C]
                        rhs += 1 + len(paths(G, rest))
                if piT != rhs: bad['ii'] += 1
                if (Rz(T, N, z) - piT) % 2: bad['redT'] += 1
                if (Rz(G, n, z) - piG) % 2: bad['redG'] += 1
                if (L - Pzu - Puz) % 2: bad['F2a'] += 1
                if (M - L) % 2: bad['F2b'] += 1
print(f"rows (G of order 2-5, ordered z != u): {rows}; failures: {bad}")
