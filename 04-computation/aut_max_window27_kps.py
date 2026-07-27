# aut_max_window27_kps.py -- kind-pasteur 2026-07-26
# The NEXT crossover after m=9 is m=27 (the tower point). Close its window:
#   f(27) = 1171875 = 3*5^8,  o(27) = 3^13 = 1594323.
# Claim: the only |Aut| value achievable in [f(27), o(27)] is 3^13, and the
# only tournaments attaining it are the Sylow-3-invariant ones, all
# isomorphic to the tower T27 = C3[T9,T9,T9] (strong, non-rigid by SS5).
#
# Step 1 (intransitive cap): for any odd G <= S27 with >= 2 orbits,
#     |G| <= max over partitions 27 = n_1+...+n_k (k>=2, n_i odd) of prod t(n_i).
#     Compute this max; it must be < f(27).
# Step 2 (transitive => 3-group): every odd-order transitive G on 27 points is
#     a 3-group of order <= 3^13 [proved in the deliverable: 27 = 3^3, odd
#     primitive <= odd part of |AGL(3,3)| = 9477; imprimitive embeds in
#     H wr K with H,K odd transitive of degree 9 or 3, all 3-groups by the
#     degree-9 window argument]. Only 3-power in [f,o]: 3^13 (3^12 = 531441 < f).
# Step 3 (uniqueness): enumerate ALL tournaments invariant under the standard
#     Sylow-3 subgroup W = Z3 wr Z3 wr Z3 of S27 (orbit analysis on ordered
#     pairs), check each is strong and isomorphic to T27.

from itertools import product

# ---- t-table (from theory layer) ----
t = {1: 1, 3: 3, 5: 5, 7: 21, 9: 81, 11: 55, 13: 39, 15: 1215,
     17: 17, 19: 171, 21: 45927, 23: 253, 25: 15625, 27: 3 ** 13}

f27 = 3 * 5 ** 8
o27 = 3 ** 13
print(f"f(27) = {f27}, o(27) = 3^13 = {o27}")

# Step 1: best multi-part partition product
best = {0: 1}
def bestprod(m):
    if m in best:
        return best[m]
    b = 0
    for part in range(1, m + 1, 2):
        b = max(b, t[part] * bestprod(m - part))
    best[m] = b
    return b

multi = max(t[part] * bestprod(27 - part) for part in range(1, 27, 2))  # part < 27 => >=2 parts
print(f"Step 1: max prod t(n_i) over partitions of 27 with >= 2 parts = {multi}"
      f"  (< f(27)? {multi < f27})")
assert multi < f27

# 3-powers in window
p3 = [3 ** k for k in range(14) if f27 <= 3 ** k <= o27]
print(f"Step 2: 3-powers in [f(27), o(27)]: {p3} (3^12 = {3**12} is below f)")
assert p3 == [3 ** 13]

# Step 3: W-invariant tournaments
def perm_from_cycles(cycles, n=27):
    p = list(range(n))
    for c in cycles:
        for i in range(len(c)):
            p[c[i]] = c[(i + 1) % len(c)]
    return tuple(p)

gens = []
# level-1 rotations: rotate each inner triple independently (9 gens)
for b in range(9):
    gens.append(perm_from_cycles([[3 * b, 3 * b + 1, 3 * b + 2]]))
# level-2 rotations: within each middle block (of 3 inner triples), rotate the
# 3 inner triples positionwise (3 gens)
for mb in range(3):
    base = 9 * mb
    cyc = [[base + i, base + 3 + i, base + 6 + i] for i in range(3)]
    gens.append(perm_from_cycles(cyc))
# level-3 rotation: rotate the 3 middle blocks positionwise
gens.append(perm_from_cycles([[i, 9 + i, 18 + i] for i in range(9)]))

# orbits of W on ordered pairs (only generators needed: BFS closure)
orbits = []
seen = set()
for u in range(27):
    for v in range(27):
        if u == v or (u, v) in seen:
            continue
        orb = set()
        st = [(u, v)]
        while st:
            (a, b) = st.pop()
            if (a, b) in orb:
                continue
            orb.add((a, b))
            for g in gens:
                st.append((g[a], g[b]))
        orbits.append(orb)
        seen |= orb
sizes = sorted(len(o) for o in orbits)
print(f"Step 3: orbits of W on ordered pairs: {len(orbits)}, sizes {sizes}")

rev = {}
for i, o1 in enumerate(orbits):
    a, b = next(iter(o1))
    for j, o2 in enumerate(orbits):
        if (b, a) in o2:
            rev[i] = j
for i in rev:
    assert rev[i] != i, "reversal inside an orbit (impossible: W has odd order)"
classes = sorted(set(frozenset((i, rev[i])) for i in rev), key=lambda fs: min(fs))
print(f"        reversal-paired classes: {len(classes)} => 2^{len(classes)}"
      f" = {2 ** len(classes)} W-invariant tournaments")

def is_strong(A):
    n = len(A)
    for adj in (A, [[A[j][i] for j in range(n)] for i in range(n)]):
        seenv = {0}; st = [0]
        while st:
            x = st.pop()
            for w in range(n):
                if adj[x][w] and w not in seenv:
                    seenv.add(w); st.append(w)
        if len(seenv) != n:
            return False
    return True

def iso(A, B):
    n = len(A)
    oA = [sum(r) for r in A]; oB = [sum(r) for r in B]
    if sorted(oA) != sorted(oB):
        return False
    perm = [-1] * n; used = [False] * n
    def bt(k):
        if k == n:
            return True
        for w in range(n):
            if used[w] or oA[k] != oB[w]:
                continue
            ok = True
            for j in range(k):
                if A[k][j] != B[w][perm[j]] or A[j][k] != B[perm[j]][w]:
                    ok = False; break
            if ok:
                perm[k] = w; used[w] = True
                if bt(k + 1):
                    return True
                used[w] = False
        return False
    return bt(0)

# reference tower T27
def circulant(n, S):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i + s) % n] = 1
    return A
def lex(Q, parts):
    sizes = [len(p) for p in parts]
    off = [0]
    for s in sizes:
        off.append(off[-1] + s)
    n = off[-1]
    A = [[0] * n for _ in range(n)]
    for b, p in enumerate(parts):
        for i in range(sizes[b]):
            for j in range(sizes[b]):
                A[off[b] + i][off[b] + j] = p[i][j]
    for b1 in range(len(parts)):
        for b2 in range(len(parts)):
            if b1 != b2 and Q[b1][b2]:
                for i in range(sizes[b1]):
                    for j in range(sizes[b2]):
                        A[off[b1] + i][off[b2] + j] = 1
    return A
C3 = circulant(3, {1})
T9 = lex(C3, [C3, C3, C3])
T27 = lex(C3, [T9, T9, T9])

results = []
for choice in product(*[sorted(fs) for fs in classes]):
    A = [[0] * 27 for _ in range(27)]
    for oi in choice:
        for (a, b) in orbits[oi]:
            A[a][b] = 1
    for i in range(27):
        for j in range(i + 1, 27):
            assert A[i][j] + A[j][i] == 1
    results.append((is_strong(A), iso(A, T27)))
print(f"        invariant tournaments: {len(results)}; "
      f"all strong: {all(r[0] for r in results)}; "
      f"all iso to T27: {all(r[1] for r in results)}")

print()
print("CONCLUSION (m=27): only |Aut| value in [f(27), o(27)] is 3^13, attained")
print("exactly by the tower T27 = C3[T9,T9,T9] (strong). THM-2453 SS5")
print("(C3[A,B,C] rigid iff all parts singletons) closes this window: tc(T27) > 1.")
