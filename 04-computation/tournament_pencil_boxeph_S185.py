#!/usr/bin/env python3
"""
tournament_pencil_boxeph_S185.py  (HYP-8560, THM-1760)

THE TOURNAMENT PENCIL:  M_T(t,u,v) = det( t I + u D_in - v A ).

(P1) EXPANSION over linear subdigraphs (proof = permutation expansion;
     verified exactly here): M = sum over vertex-disjoint directed-cycle
     sets L of (-1)^{#cycles(L)} v^{|V(L)|} prod_{w not in L} (t + u indeg w).
(P2) THE THREE FACES:
     v=0: SCORE/CUT face  M = prod_j (t + u indeg(j))     [in-score product]
     u=0: SPECTRAL/CYCLE face  M = char poly of A          [THM-506 signed face]
     u=v: FOREST ray  M(t,1,1) = det(tI + L_in) = sum_F t^{comps}
          [matrix-forest; Kirchhoff at d/dt|_0 = sum_r a_r]
     (out-score face = v=0 face of T^op, consistent via P3.)
(P3) COMPLEMENT FUNCTIONAL EQUATION (matrix determinant lemma; A^T = J-I-A,
     D_in^op = (n-1)I - D_in):
       M_{T^op}(t,u,v) = M_T(t~, -u, -v) - v * S_T(t~, -u, -v),
     t~ = t + u(n-1) + v,  S_T(a,b,c) = 1^T adj(aI + bD_in - cA) 1.
     The Z2/complement symmetry forces the DOUBLY-ROOTED (adjugate-sum)
     object into the picture — except on the forest diagonal, where column
     sums collapse it: S(a,b,b) = n M(a,b,b)/a, giving the classic
     complement-spectrum law  det(tI + L_in(T^op)) =
     det((t+n)I - L_in(T)) * t/(t+n)  as the diagonal specialization.
(P4) FINGERPRINT at n=7: the pencil evaluated on the 8x8x8 grid separates
     how many of the 456 classes? In particular klein-THM-1750's 4 resistant
     pairs (same spec A, same sum a, same {a_r} vector, same H): does the
     pencil split them (score face or cross terms)?

boxeph-2026-07-20-S185. Pure python.
"""
import itertools
import math

def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def bits_to_adj(bits, n):
    adj = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs_of(n)):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def build_perm_maps(n):
    prs = pairs_of(n)
    idx = {pr: k for k, pr in enumerate(prs)}
    maps = []
    for p in itertools.permutations(range(n)):
        pm = []
        for (i, j) in prs:
            a, b = p[i], p[j]
            pm.append((idx[(a, b)], 0) if a < b else (idx[(b, a)], 1))
        maps.append(pm)
    return maps


def classes_for_n(n):
    if n == 7:
        with open('05-knowledge/results/n7_class_reps_boxeph_S152.txt') as f:
            return [int(t) for t in f.read().split()]
    perm_maps = build_perm_maps(n)
    total = 1 << (n * (n - 1) // 2)
    visited = bytearray(total)
    reps = []
    for bits in range(total):
        if visited[bits]:
            continue
        members = set()
        for pm in perm_maps:
            nb = 0
            for k in range(len(pm)):
                src, flip = pm[k]
                b = (bits >> src) & 1
                if flip:
                    b ^= 1
                nb |= b << k
            members.add(nb)
        for m in members:
            visited[m] = 1
        reps.append(min(members))
    return reps


def det_int(mat):
    m = [row[:] for row in mat]
    n = len(m)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if m[k][k] == 0:
            piv = next((r for r in range(k + 1, n) if m[r][k] != 0), None)
            if piv is None:
                return 0
            m[k], m[piv] = m[piv], m[k]
            sign = -sign
        mkk = m[k][k]
        for i in range(k + 1, n):
            mik = m[i][k]
            for j in range(k + 1, n):
                m[i][j] = (m[i][j] * mkk - mik * m[k][j]) // prev
            m[i][k] = 0
        prev = mkk
    return sign * m[n - 1][n - 1]


def Mpencil(adj, n, t, u, v):
    indeg = [sum(adj[i][j] for i in range(n)) for j in range(n)]
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        M[i][i] = t + u * indeg[i]
        for j in range(n):
            if i != j and adj[i][j]:
                M[i][j] -= v
    return det_int(M)


def adjsum(adj, n, a, b, c):
    # 1^T adj(aI + b D_in - c A) 1 = sum of all cofactors C_ij
    indeg = [sum(adj[i][j] for i in range(n)) for j in range(n)]
    X = [[0] * n for _ in range(n)]
    for i in range(n):
        X[i][i] = a + b * indeg[i]
        for j in range(n):
            if i != j and adj[i][j]:
                X[i][j] -= c
    tot = 0
    for i in range(n):
        for j in range(n):
            minor = [[X[r][s] for s in range(n) if s != j] for r in range(n) if r != i]
            tot += (-1) ** (i + j) * det_int(minor)
    return tot


def linear_subdigraph_sum(adj, n, t, u, v):
    indeg = [sum(adj[i][j] for i in range(n)) for j in range(n)]
    # enumerate all sets of vertex-disjoint directed cycles (len >= 3)
    cycles = []
    for size in range(3, n + 1):
        for sub in itertools.combinations(range(n), size):
            first = sub[0]
            rest = sub[1:]
            for perm in itertools.permutations(rest):
                cyc = (first,) + perm
                if all(adj[cyc[k]][cyc[(k + 1) % size]] for k in range(size)):
                    cycles.append((frozenset(sub), cyc))
    total = 0
    # sum over disjoint packings, by recursion
    def rec(idx, used, ncyc, nvert):
        nonlocal total
        # count this packing
        term = (-1) ** ncyc * v ** nvert
        for w in range(n):
            if not (used >> w) & 1:
                term *= (t + u * indeg[w])
        total += term
        for k in range(idx, len(cycles)):
            S, _ = cycles[k]
            Sm = sum(1 << x for x in S)
            if Sm & used:
                continue
            rec(k + 1, used | Sm, ncyc + 1, nvert + len(S))
    rec(0, 0, 0, 0)
    return total


def forests_by_comps(adj, n):
    # spanning out-forest counts by #components (enumeration, n <= 5)
    res = {}
    for roots_mask in range(1, 1 << n):
        roots = [v for v in range(n) if (roots_mask >> v) & 1]
        others = [v for v in range(n) if not (roots_mask >> v) & 1]
        choices = [[u for u in range(n) if adj[u][v]] for v in others]
        cnt = 0
        for combo in itertools.product(*choices):
            par = dict(zip(others, combo))
            ok = True
            for v0 in others:
                seen = set()
                x = v0
                while x not in roots:
                    if x in seen:
                        ok = False
                        break
                    seen.add(x)
                    x = par[x]
                if not ok:
                    break
            if ok:
                cnt += 1
        if cnt:
            res[len(roots)] = res.get(len(roots), 0) + cnt
    return res


def charpoly_A(adj, n):
    # Faddeev-LeVerrier, integer char poly of A: x^n + c1 x^{n-1} + ... + cn
    A = [[adj[i][j] for j in range(n)] for i in range(n)]
    Mk = [[0] * n for _ in range(n)]
    cs = [1]
    for k in range(1, n + 1):
        # M_k = A M_{k-1} + c_{k-1} I
        if k == 1:
            Mk = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
        else:
            prod = [[sum(A[i][x] * Mk[x][j] for x in range(n)) for j in range(n)]
                    for i in range(n)]
            for i in range(n):
                prod[i][i] += cs[-1]
            Mk = prod
        AM = [[sum(A[i][x] * Mk[x][j] for x in range(n)) for j in range(n)]
              for i in range(n)]
        tr = sum(AM[i][i] for i in range(n))
        assert tr % k == 0
        cs.append(-tr // k)
    return tuple(cs)


print("=" * 78)
print("THE TOURNAMENT PENCIL  M_T(t,u,v) = det(tI + u D_in - v A)")
print("=" * 78)

print("\n(P1) expansion over linear subdigraphs — exact, all classes n=3..5,")
print("     grid (t,u,v) in {0..3}^3:")
bad = tot = 0
for n in (3, 4, 5):
    for bits in classes_for_n(n):
        adj = bits_to_adj(bits, n)
        for t in range(4):
            for u in range(4):
                for v in range(4):
                    tot += 1
                    if Mpencil(adj, n, t, u, v) != linear_subdigraph_sum(adj, n, t, u, v):
                        bad += 1
print("     %d evaluations, %d mismatches" % (tot, bad))

print("\n(P2) faces (n=3..5, all classes):")
okS = okC = okF = okK = tot = 0
for n in (3, 4, 5):
    for bits in classes_for_n(n):
        adj = bits_to_adj(bits, n)
        indeg = [sum(adj[i][j] for i in range(n)) for j in range(n)]
        tot += 1
        # score face
        okS += all(Mpencil(adj, n, t, u, 0) ==
                   math.prod(t + u * indeg[j] for j in range(n))
                   for t in range(3) for u in range(3))
        # spectral face: M(t,0,1) = charpoly_A evaluated with sign x->t
        cs = charpoly_A(adj, n)
        okC += all(Mpencil(adj, n, t, 0, 1) ==
                   sum(cs[k] * t ** (n - k) for k in range(n + 1))
                   for t in range(4))
        # forest face
        fb = forests_by_comps(adj, n)
        okF += all(Mpencil(adj, n, t, 1, 1) ==
                   sum(cnt * t ** c for c, cnt in fb.items())
                   for t in range(4))
        # Kirchhoff: d/dt at 0 = sum_r a_r = coefficient of t^1
        m0 = Mpencil(adj, n, 0, 1, 1)
        m1 = Mpencil(adj, n, 1, 1, 1)
        m_1 = Mpencil(adj, n, -1, 1, 1)
        # f(t) = sum c_k t^k: c1 = (f(1)-f(-1))/2 for cubic-free... use exact:
        # f(0)=0 always (column sums), c1 = fb.get(1,0)
        okK += (m0 == 0 and fb.get(1, 0) == (m1 - m_1) // 2 - 0 * m0
                if n <= 3 else m0 == 0)
print("     classes checked: %d | score-face ok %d | spectral-face ok %d | forest-face ok %d | det(L)=0 ok %d"
      % (tot, okS, okC, okF, okK))

print("\n(P3) complement functional equation, exact, all classes n=3..5,")
print("     grid (t,u,v) in {1..3}^3 (and the diagonal collapse):")
bad = tot = 0
badd = totd = 0
for n in (3, 4, 5):
    for bits in classes_for_n(n):
        adj = bits_to_adj(bits, n)
        adjop = [[adj[j][i] for j in range(n)] for i in range(n)]
        for t in range(1, 4):
            for u in range(1, 4):
                for v in range(1, 4):
                    tt = t + u * (n - 1) + v
                    lhs = Mpencil(adjop, n, t, u, v)
                    rhs = Mpencil(adj, n, tt, -u, -v) - v * adjsum(adj, n, tt, -u, -v)
                    tot += 1
                    if lhs != rhs:
                        bad += 1
        # diagonal collapse: det(tI + L^op) * (t + n u) = det((t+nu)I - u L) * t
        for t in range(1, 5):
            for u in range(1, 3):
                lhs = Mpencil(adjop, n, t, u, u) * (t + n * u)
                rhs = Mpencil(adj, n, t + n * u, -u, -u) * t
                totd += 1
                if lhs != rhs:
                    badd += 1
print("     general: %d evaluations, %d mismatches | diagonal law: %d evals, %d mismatches"
      % (tot, bad, totd, badd))

print("\n(P4) fingerprint strength at n=7 (grid {0..7}^3 = full coefficient info):")
reps = classes_for_n(7)
n = 7
prints = {}
extra = {}
for bits in reps:
    adj = bits_to_adj(bits, n)
    vec = []
    for t in range(8):
        for u in range(8):
            for v in range(8):
                vec.append(Mpencil(adj, n, t, u, v))
    key = tuple(vec)
    prints.setdefault(key, []).append(bits)
groups = [v for v in prints.values() if len(v) > 1]
print("     456 classes -> %d distinct pencil polynomials; unresolved groups: %d"
      % (len(prints), len(groups)))
if groups:
    print("     unresolved: %s" % groups)
    # do the unresolved share (specA, sum a, vector, H)? (klein's residue check)
    for g in groups:
        for bits in g:
            adj = bits_to_adj(bits, n)
            cs = charpoly_A(adj, n)
            print("       bits=%d charpolyA=%s" % (bits, cs))
print("\nDONE.")
