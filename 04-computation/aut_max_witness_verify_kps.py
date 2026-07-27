# aut_max_witness_verify_kps.py -- kind-pasteur 2026-07-26
# Verification layer:
#  A. Build the o(m) witness tournaments for m = 3..21 and count |Aut| EXACTLY
#     (backtracking), check strongness, compare with the theory table.
#  B. m = 9 crossover window [f(9), o(9)] = [75, 81]:
#     - Lagrange: which odd integers in [75,81] divide oddpart(9!)?
#     - |Aut| = 81 = full Sylow-3 of S9  =>  enumerate ALL tournaments invariant
#       under a fixed Sylow-3 subgroup P = Z3 wr Z3, check they are all
#       isomorphic to T9 = C3[C3,C3,C3], strong, |Aut| = 81.
#     - H(T9) by subset-DP; tc = H/81; independent SS5-style run-transfer
#       cross-check of H(T9).
#  C. Controls from canon THM-2453: C3[C3,1,1] must give H=15, |Aut|=3, tc=5.

import sys
from itertools import permutations, product
from math import factorial

# ---------- tournament constructors (adjacency: A[i][j] = 1 iff i -> j) ----------

def circulant(n, S):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i + s) % n] = 1
    return A

def C3():  return circulant(3, {1})
def R5():  return circulant(5, {1, 2})
def P7():  return circulant(7, {1, 2, 4})
def PT():  return [[0]]

def lex(Q, parts):
    """lexicographic composition Q[parts]: |parts| == len(Q)"""
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

def stack(parts):
    """transitive join: every vertex of an earlier part beats every later part"""
    k = len(parts)
    Q = [[1 if i < j else 0 for j in range(k)] for i in range(k)]
    return lex(Q, parts)

def T9():  return lex(C3(), [C3(), C3(), C3()])

# ---------- exact automorphism count (backtracking) ----------

def aut_count(A):
    n = len(A)
    outd = [sum(A[i]) for i in range(n)]
    cnt = 0
    perm = [-1] * n
    used = [False] * n
    def bt(k):
        nonlocal cnt
        if k == n:
            cnt += 1
            return
        for w in range(n):
            if used[w] or outd[w] != outd[k]:
                continue
            ok = True
            Ak = A[k]; Aw = A[w]
            for j in range(k):
                pj = perm[j]
                if Ak[j] != Aw[pj] or A[j][k] != A[pj][w]:
                    ok = False
                    break
            if ok:
                perm[k] = w
                used[w] = True
                bt(k + 1)
                used[w] = False
        perm[k] = -1
    bt(0)
    return cnt

def is_strong(A):
    n = len(A)
    for adj in (A, [[A[j][i] for j in range(n)] for i in range(n)]):
        seen = {0}
        st = [0]
        while st:
            v = st.pop()
            for w in range(n):
                if adj[v][w] and w not in seen:
                    seen.add(w)
                    st.append(w)
        if len(seen) != n:
            return False
    return True

def ham_paths(A):
    n = len(A)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            for w in range(n):
                if A[v][w] and not (mask >> w) & 1:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[full])

def iso(A, B):
    """isomorphism test by backtracking (small n)"""
    n = len(A)
    if n != len(B):
        return False
    oA = [sum(r) for r in A]; oB = [sum(r) for r in B]
    if sorted(oA) != sorted(oB):
        return False
    perm = [-1] * n
    used = [False] * n
    def bt(k):
        if k == n:
            return True
        for w in range(n):
            if used[w] or oA[k] != oB[w]:
                continue
            ok = True
            for j in range(k):
                if A[k][j] != B[w][perm[j]] or A[j][k] != B[perm[j]][w]:
                    ok = False
                    break
            if ok:
                perm[k] = w; used[w] = True
                if bt(k + 1):
                    return True
                used[w] = False
        perm[k] = -1
        return False
    return bt(0)

# ---------- A. witness verification ----------

print("=" * 78)
print("A. WITNESS VERIFICATION: exact |Aut| by backtracking, strongness")
print("=" * 78)

P7C3 = lambda: lex(P7(), [C3() for _ in range(7)])
R5C3 = lambda: lex(R5(), [C3() for _ in range(5)])

witnesses = [
    (3,  "C3",                 C3(),                              3),
    (4,  "stack[C3,pt]",       stack([C3(), PT()]),               3),
    (5,  "R5",                 R5(),                              5),
    (6,  "stack[C3,C3]",       stack([C3(), C3()]),               9),
    (7,  "P7 (Paley)",         P7(),                             21),
    (8,  "stack[P7,pt]",       stack([P7(), PT()]),              21),
    (9,  "T9=C3[C3,C3,C3]",    T9(),                             81),
    (10, "stack[T9,pt]",       stack([T9(), PT()]),              81),
    (11, "stack[T9,pt,pt]",    stack([T9(), PT(), PT()]),        81),
    (12, "stack[T9,C3]",       stack([T9(), C3()]),             243),
    (13, "stack[T9,C3,pt]",    stack([T9(), C3(), PT()]),       243),
    (14, "stack[P7,P7]",       stack([P7(), P7()]),             441),
    (15, "R5[C3 x5]",          R5C3(),                         1215),
    (16, "stack[T9,P7]",       stack([T9(), P7()]),            1701),
    (17, "stack[T9,P7,pt]",    stack([T9(), P7(), PT()]),      1701),
    (18, "stack[T9,T9]",       stack([T9(), T9()]),            6561),
    (19, "stack[T9,T9,pt]",    stack([T9(), T9(), PT()]),      6561),
    (20, "stack[T9,T9,pt,pt]", stack([T9(), T9(), PT(), PT()]), 6561),
    (21, "P7[C3 x7]",          P7C3(),                        45927),
]

allok = True
for m, name, A, expect in witnesses:
    a = aut_count(A)
    s = is_strong(A)
    ok = (a == expect)
    allok &= ok
    print(f"  m={m:2d}  {name:<18} |Aut| = {a:>6d}  (theory lower bound {expect:>6d})"
          f"  strong={str(s):5}  {'OK' if ok else '*** MISMATCH ***'}")
print(f"  all witnesses match theory: {allok}")

# ---------- C. canon controls (THM-2453) ----------

print()
print("=" * 78)
print("C. CANON CONTROLS (THM-2453 SS1/SS2 values)")
print("=" * 78)
ctrl = lex(C3(), [C3(), PT(), PT()])   # C3[C3,1,1], strong tower, tc = 5
h, a = ham_paths(ctrl), aut_count(ctrl)
print(f"  C3[C3,1,1]: H = {h} (canon 15), |Aut| = {a} (canon 3), tc = {h//a} (canon 5),"
      f" strong = {is_strong(ctrl)}")
assert (h, a) == (15, 3)
hc3, ac3 = ham_paths(C3()), aut_count(C3())
print(f"  C3:         H = {hc3}, |Aut| = {ac3}, tc = {hc3//ac3}  (the GOAL equality case)")
assert (hc3, ac3) == (3, 3)

# ---------- B. the m = 9 window ----------

print()
print("=" * 78)
print("B. THE m=9 CROSSOVER WINDOW [f(9), o(9)] = [75, 81]")
print("=" * 78)

# B1. Lagrange sieve
odd9 = factorial(9)
while odd9 % 2 == 0:
    odd9 //= 2
print(f"  oddpart(9!) = {odd9} = 3^4 * 5 * 7")
cand = [k for k in range(75, 82) if k % 2 == 1]
achievable = [k for k in cand if odd9 % k == 0]
print(f"  odd integers in [75,81]: {cand}")
print(f"  ... dividing oddpart(9!): {achievable}   "
      f"(75=3*5^2 fails: 25 does not divide 9!; 77=7*11 fails: 11 > 9; 79 prime > 9)")
assert achievable == [81]

# B2. |Aut| = 81 = 3^4 = FULL Sylow-3 of S9 (v_3(9!) = 4) => Aut(T) is a Sylow-3
#     subgroup; all are conjugate, so wlog Aut(T) contains the standard
#     P = <(012),(345),(678),(036)(147)(258)> = Z3 wr Z3.
#     Enumerate ALL tournaments invariant under P.

def perm_cycle(c, n=9):
    p = list(range(n))
    for i in range(len(c)):
        p[c[i]] = c[(i + 1) % len(c)]
    return tuple(p)

def compose(p, q):  # (p*q)(x) = p(q(x))
    return tuple(p[q[x]] for x in range(9))

gens = [perm_cycle([0, 1, 2]), perm_cycle([3, 4, 5]), perm_cycle([6, 7, 8])]
s = list(range(9))
for i in range(3):
    for b in range(3):
        s[3 * b + i] = 3 * ((b + 1) % 3) + i
gens.append(tuple(s))

# generate the group
P = {tuple(range(9))}
frontier = list(P)
while frontier:
    new = []
    for g in frontier:
        for h in gens:
            x = compose(h, g)
            if x not in P:
                P.add(x)
                new.append(x)
    frontier = new
print(f"  |P| = {len(P)} (expected 81 = Z3 wr Z3, the Sylow-3 of S9)")
assert len(P) == 81

# orbits of P on ordered pairs; check no orbit contains a reversal; enumerate
pairs = [(i, j) for i in range(9) for j in range(9) if i != j]
seen = set()
orbits = []
for pr in pairs:
    if pr in seen:
        continue
    orb = set()
    st = [pr]
    while st:
        (u, v) = st.pop()
        if (u, v) in orb:
            continue
        orb.add((u, v))
        for g in gens:
            st.append((g[u], g[v]))
    orbits.append(orb)
    seen |= orb
print(f"  orbits of P on ordered pairs: {len(orbits)} with sizes {sorted(len(o) for o in orbits)}")

# pair up orbits with their reversals
rev = {}
for i, o1 in enumerate(orbits):
    u, v = next(iter(o1))
    for j, o2 in enumerate(orbits):
        if (v, u) in o2:
            rev[i] = j
for i in rev:
    assert rev[i] != i, "an orbit contains a reversal -- impossible for odd order"
free = sorted(set(frozenset((i, rev[i])) for i in rev))
print(f"  reversal-paired orbit classes: {len(free)}  =>  2^{len(free)} = {2**len(free)} invariant tournaments")

inv_tournaments = []
for choice in product(*[tuple(fs) for fs in free]):
    A = [[0] * 9 for _ in range(9)]
    for oi in choice:
        for (u, v) in orbits[oi]:
            A[u][v] = 1
    # verify tournament
    for i in range(9):
        for j in range(i + 1, 9):
            assert A[i][j] + A[j][i] == 1
    inv_tournaments.append(A)

T9ref = T9()
print(f"  P-invariant tournaments: {len(inv_tournaments)}")
for k, A in enumerate(inv_tournaments):
    print(f"    #{k}: strong={str(is_strong(A)):5}  |Aut|={aut_count(A):3d}  "
          f"H={ham_paths(A):5d}  iso to T9? {iso(A, T9ref)}")

# B3. H(T9) and tc
H9 = ham_paths(T9ref)
A9 = aut_count(T9ref)
print(f"\n  T9 = C3[C3,C3,C3]:  H = {H9},  |Aut| = {A9},  tc = H/|Aut| = {H9 // A9}"
      f"  (integer: {H9 % A9 == 0}, odd: {(H9 // A9) % 2 == 1})")

# B4. independent SS5-style run-transfer cross-check of H(T9):
#     block walks on the C3 quotient are forced cyclic; a walk of k runs uses
#     run-counts (r_X,r_Y,r_Z) determined by k and the start block; each block
#     contributes p(r) = #ordered r-part path systems of C3: p(1)=3 (Ham paths),
#     p(2)=6 (arc+point, ordered), p(3)=6 (3 singletons ordered).
p = {1: 3, 2: 6, 3: 6}
total = 0
detail = []
for k in range(3, 10):          # k runs, 3 starts each
    q, r = divmod(k, 3)
    counts = [q + (1 if i < r else 0) for i in range(3)]   # run counts along the walk
    term = 3 * p[counts[0]] * p[counts[1]] * p[counts[2]]
    detail.append(f"k={k}:3*{p[counts[0]]}*{p[counts[1]]}*{p[counts[2]]}={term}")
    total += term
print(f"  SS5 run-transfer cross-check: H(T9) = {total}   [{'; '.join(detail)}]")
assert total == H9

print()
print("CONCLUSION (m=9): the only strong tournament on 9 vertices with")
print(f"|Aut| >= f(9) = 75 is T9 (unique with |Aut| = 81), and tc(T9) = {H9 // A9} != 1.")
print("Hence no rigid (H = |Aut|) strong 9-tournament exists: the window closes.")
