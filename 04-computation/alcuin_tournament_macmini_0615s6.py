#!/usr/bin/env python3
"""
alcuin_tournament_macmini_0615s6.py  (mac-mini-2026-06-15-S6)

Independent verification of the Alcuin <-> vertex-cover <-> tournament bridge.
Conflict graph G on ordered [n] -> tournament T_G: for i<j, arc i->j (FORWARD) iff
{i,j} in E(G), else j->i (BACKWARD). forward-arc set = E(G).

Checks:
 (A) tau(G) <= Alcuin(G) <= tau(G)+1 over all non-iso graphs (n<=6, +n=7 if fast); the
     +1 fraction and a structural characterization.
 (B) correspondence: Redei parity (#ham-paths odd); independent set <-> largest
     reverse-transitive monotone set (alpha) and largest transitive subtournament;
     #3-cycles(T_G) == #ordered triples i<j<k that are induced-P3 in G or complement.
 (C) minor-monotonicity: is tau minor-monotone? is Alcuin? (delete-vertex/edge, contract)
 (D) conjectures: Alcuin=tau+1 vs T_G non-strong (natural & all orders);
     G-ham-cycle vs exists-order-T_G-strong; G-ham-path vs exists-order-spine-forward.
"""
import sys, itertools
from itertools import combinations, permutations
sys.stdout.reconfigure(line_buffering=True)

# ---------- graph enumeration (non-iso) ----------
def all_noniso(n):
    """yield one representative (frozenset of edges) per iso class on n labeled vertices."""
    pairs = list(combinations(range(n), 2))
    seen = set()
    reps = []
    perms = list(permutations(range(n)))
    for mask in range(1 << len(pairs)):
        edges = frozenset(pairs[i] for i in range(len(pairs)) if (mask >> i) & 1)
        # canonical key = min over relabelings
        best = None
        for p in perms:
            relab = frozenset((min(p[a], p[b]), max(p[a], p[b])) for (a, b) in edges)
            key = tuple(sorted(relab))
            if best is None or key < best:
                best = key
        if best not in seen:
            seen.add(best)
            reps.append((n, edges))
    return reps

def neighbors(n, E):
    adj = {v: set() for v in range(n)}
    for a, b in E:
        adj[a].add(b); adj[b].add(a)
    return adj

def is_independent(S, Eset):
    for a, b in combinations(sorted(S), 2):
        if (a, b) in Eset or (b, a) in Eset:
            return False
    return True

def alpha_number(n, E):
    Eset = set(E) | set((b, a) for a, b in E)
    best = 0
    for r in range(n, -1, -1):
        if r <= best: break
        found = False
        for S in combinations(range(n), r):
            if is_independent(S, Eset):
                best = r; found = True; break
        if found: break
    return best

def tau_number(n, E):
    return n - alpha_number(n, E)

# ---------- Alcuin number via state-space BFS ----------
def alcuin_number(n, E):
    if n == 0: return 0
    Eset = set(E) | set((b, a) for a, b in E)
    def indep(frozen):
        for a, b in combinations(sorted(frozen), 2):
            if (a, b) in Eset: return False
        return True
    full = frozenset(range(n))
    for b in range(1, n + 1):
        # state: (left_frozenset, boat_side) boat_side: 0=L,1=R
        start = (full, 0); goal = (frozenset(), 1)
        from collections import deque
        seen = {start}; dq = deque([start])
        ok = False
        while dq:
            left, side = dq.popleft()
            if (left, side) == goal: ok = True; break
            # ferryman side items
            bankitems = left if side == 0 else (full - left)
            # choose subset M of bankitems size 0..b
            items = sorted(bankitems)
            for k in range(0, b + 1):
                for M in combinations(items, k):
                    Mset = frozenset(M)
                    if side == 0:
                        newleft = left - Mset; newside = 1
                        leftless = newleft  # side without ferryman = L
                    else:
                        newleft = left | Mset; newside = 0
                        leftless = full - newleft  # side without ferryman = R
                    if indep(leftless):
                        st = (newleft, newside)
                        if st not in seen:
                            seen.add(st); dq.append(st)
        if ok:
            return b
    return n

# ---------- tournament T_G ----------
def build_TG(n, E, order=None):
    """arc matrix A[i][j]=1 if i->j. order: permutation giving the labeling (default identity)."""
    if order is None: order = list(range(n))
    Eset = set(E) | set((b, a) for a, b in E)
    pos = {order[i]: i for i in range(n)}  # vertex -> position
    A = [[0] * n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u == v: continue
            # by POSITION: lower position -> higher position is forward iff edge
            lo, hi = (u, v) if pos[u] < pos[v] else (v, u)
            edge = (lo, hi) in Eset
            # forward arc lo->hi iff edge; else hi->lo
            if edge:
                A[lo][hi] = 1
            else:
                A[hi][lo] = 1
    return A

def out_scores(A, n):
    return [sum(A[i]) for i in range(n)]

def num_3cycles(A, n):
    c = 0
    for i, j, k in combinations(range(n), 3):
        # cyclic if there is a directed triangle among the 3
        verts = [i, j, k]
        # check both orientations of cycle
        def cyc(a, b, c): return A[a][b] and A[b][c] and A[c][a]
        if cyc(i, j, k) or cyc(i, k, j):
            c += 1
    return c

def is_strong(A, n):
    if n <= 1: return True
    def reach(s):
        seen = {s}; st = [s]
        while st:
            x = st.pop()
            for y in range(n):
                if A[x][y] and y not in seen:
                    seen.add(y); st.append(y)
        return seen
    for s in range(n):
        if len(reach(s)) != n: return False
    return True

def count_ham_paths(A, n):
    cnt = 0
    for perm in permutations(range(n)):
        if all(A[perm[t]][perm[t + 1]] for t in range(n - 1)):
            cnt += 1
    return cnt

def largest_transitive_subtournament(A, n):
    def acyclic(S):
        # transitive subtournament <=> induced is acyclic <=> topo order exists
        S = list(S)
        # try: is there a linear order with all arcs forward? equivalently no directed cycle
        # check via: repeatedly remove a source
        nodes = set(S)
        adj = {u: set(v for v in S if v != u and A[u][v]) for u in S}
        indeg = {u: sum(1 for v in S if v != u and A[v][u]) for u in S}
        from collections import deque
        q = deque([u for u in S if indeg[u] == 0])
        removed = 0
        indeg = dict(indeg)
        while q:
            u = q.popleft(); removed += 1
            for v in adj[u]:
                indeg[v] -= 1
                if indeg[v] == 0: q.append(v)
        return removed == len(S)
    for r in range(n, 0, -1):
        for S in combinations(range(n), r):
            if acyclic(S): return r
    return 0

# ---------- graph predicates ----------
def has_ham_path(n, E):
    Eset = set(E) | set((b, a) for a, b in E)
    for perm in permutations(range(n)):
        if all((perm[t], perm[t + 1]) in Eset for t in range(n - 1)):
            return True
    return False

def has_ham_cycle(n, E):
    if n < 3: return False
    Eset = set(E) | set((b, a) for a, b in E)
    for perm in permutations(range(n)):
        if perm[0] != 0: continue
        if all((perm[t], perm[(t + 1) % n]) in Eset for t in range(n)):
            return True
    return False

def num_induced_ordered_P3(n, E):
    """#ordered triples i<j<k forming a 3-cycle in T_G (natural order):
       pattern (edge ij, edge jk, non-edge ik) OR (non-edge ij, non-edge jk, edge ik)."""
    Eset = set(E) | set((b, a) for a, b in E)
    def e(a, b): return (min(a, b), max(a, b)) in Eset
    c = 0
    for i, j, k in combinations(range(n), 3):  # i<j<k
        p = (e(i, j), e(j, k), e(i, k))
        if p == (True, True, False) or p == (False, False, True):
            c += 1
    return c

# ---------- minors ----------
def delete_vertex(n, E, v):
    relabel = {u: (u if u < v else u - 1) for u in range(n) if u != v}
    Enew = frozenset((relabel[a], relabel[b]) for a, b in E if a != v and b != v)
    return (n - 1, Enew)

def delete_edge(n, E, edge):
    return (n, frozenset(x for x in E if x != edge))

def contract_edge(n, E, edge):
    u, w = edge  # merge w into u
    # new vertex u absorbs w; remove w
    Enew = set()
    for a, b in E:
        if (a, b) == (u, w) or (a, b) == (w, u): continue
        a2 = u if a == w else a
        b2 = u if b == w else b
        if a2 == b2: continue
        Enew.add((min(a2, b2), max(a2, b2)))
    # relabel to remove vertex w
    relabel = {x: (x if x < w else x - 1) for x in range(n) if x != w}
    Enew = frozenset((relabel[a], relabel[b]) for a, b in Enew)
    return (n - 1, Enew)

def all_minors_one_step(n, E):
    outs = []
    for v in range(n):
        outs.append(('del-v', delete_vertex(n, E, v)))
    for x in E:
        outs.append(('del-e', delete_edge(n, E, x)))
    for x in E:
        outs.append(('contract', contract_edge(n, E, x)))
    return outs

# ============================ MAIN ============================
print("=" * 78)
print("Alcuin <-> vertex cover <-> tournament  (mac-mini-2026-06-15-S6)")
print("=" * 78)

MAXN = 6
allreps = []
for n in range(1, MAXN + 1):
    reps = all_noniso(n)
    allreps += reps
    print(f"  n={n}: {len(reps)} non-iso graphs")

# (A) bounds + +1 characterization
print("\n--- (A) bound tau <= Alcuin <= tau+1, and the +1 cases ---")
viol = 0
plus1 = []
by_n = {}
for (n, E) in allreps:
    t = tau_number(n, E); a = alcuin_number(n, E)
    by_n.setdefault(n, [0, 0])
    if a == t: by_n[n][0] += 1
    elif a == t + 1: by_n[n][1] += 1
    if not (t <= a <= t + 1): viol += 1
    if a == t + 1:
        plus1.append((n, E, t, a))
print(f"  bound violations: {viol}")
for n in sorted(by_n):
    eq, p1 = by_n[n]
    print(f"  n={n}: Alcuin=tau in {eq}, =tau+1 in {p1}")
# characterize +1: test hypotheses
print(f"  total +1 graphs (n<=in {MAXN}): {len(plus1)}")
# hypothesis: +1 <=> graph has NO edges OR is a 'star-forest'-like / max degree forces it
import collections
def degs(n, E):
    d = [0]*n
    for a,b in E: d[a]+=1; d[b]+=1
    return d
print("  smallest +1 examples (n, edges, tau):")
for (n,E,t,a) in plus1[:8]:
    print(f"     n={n} E={sorted(E)} tau={t}")
# test: +1 <=> tau==0 OR exists a min-vertex-cover vertex adjacent to all of some structure
# empirical: check correlation with 'has isolated-from-cover indep set too big'
n_empty = sum(1 for (n,E,t,a) in plus1 if len(E)==0)
print(f"  of +1 graphs, how many are edgeless (tau=0): {n_empty}")

# (B) correspondence checks
print("\n--- (B) correspondence: Redei parity, indep<->reverse-transitive, 3cycle<->P3 ---")
redei_fail = 0; p3_fail = 0; alpha_vs_trans = collections.Counter()
for (n, E) in allreps:
    A = build_TG(n, E)
    hp = count_ham_paths(A, n)
    if hp % 2 == 0: redei_fail += 1
    c3 = num_3cycles(A, n); p3 = num_induced_ordered_P3(n, E)
    if c3 != p3: p3_fail += 1
    al = alpha_number(n, E)
    lt = largest_transitive_subtournament(A, n)
    alpha_vs_trans[('alpha', al, 'trans', lt, 'ge' if lt >= al else 'LT')] += 1
print(f"  Redei parity (#ham-paths odd) failures: {redei_fail}")
print(f"  3-cycles == ordered-induced-P3 count failures: {p3_fail}")
# is largest transitive subtournament always >= alpha (and >= omega)? and is it ever > both?
ge = sum(v for k, v in alpha_vs_trans.items() if k[4] == 'ge')
lt = sum(v for k, v in alpha_vs_trans.items() if k[4] == 'LT')
print(f"  largest-transitive >= alpha: {ge} cases; < alpha: {lt} cases (should be 0)")

# (C) minor-monotonicity
print("\n--- (C) minor-monotonicity of tau and Alcuin ---")
tau_mono_fail = []; alc_mono_fail = []
alc_del_fail = []
for (n, E) in allreps:
    if n < 2: continue
    tG = tau_number(n, E); aG = alcuin_number(n, E)
    for (op, (nh, Eh)) in all_minors_one_step(n, E):
        if nh < 1: continue
        tH = tau_number(nh, Eh); aH = alcuin_number(nh, Eh)
        if tH > tG: tau_mono_fail.append((op, n, sorted(E), tG, tH))
        if aH > aG:
            alc_mono_fail.append((op, n, sorted(E), aG, nh, sorted(Eh), aH))
            if op in ('del-v', 'del-e'):
                alc_del_fail.append((op, n, sorted(E), aG, aH))
print(f"  tau minor-monotonicity failures (tau(H)>tau(G)): {len(tau_mono_fail)}")
print(f"  Alcuin minor-monotonicity failures (Alcuin(H)>Alcuin(G)): {len(alc_mono_fail)}")
if alc_mono_fail:
    print("    SMALLEST Alcuin minor-increase counterexamples:")
    for row in sorted(alc_mono_fail, key=lambda r: (r[1], len(r[2])))[:6]:
        op, n, E, aG, nh, Eh, aH = row
        print(f"      via {op}: G(n={n},E={E},Alcuin={aG}) -> H(n={nh},E={Eh},Alcuin={aH})")
    print(f"    of those, via deletion (subgraph, not contraction): {len(alc_del_fail)}")
else:
    print("    -> Alcuin appears MINOR-MONOTONE up to n<=", MAXN)

# (D) conjectures
print("\n--- (D) conjectures linking Alcuin/Hamiltonicity to T_G strong-connectivity ---")
def exists_order_strong(n, E):
    for perm in permutations(range(n)):
        A = build_TG(n, E, list(perm))
        if is_strong(A, n): return True
    return False
# tabulate
c_plus1_strongNat = collections.Counter()
c_hamcyc_existstrong = collections.Counter()
c_hampath = 0; c_hampath_consistent = 0
for (n, E) in allreps:
    if n < 1: continue
    t = tau_number(n, E); a = alcuin_number(n, E)
    A = build_TG(n, E)
    strongNat = is_strong(A, n)
    plus1flag = (a == t + 1)
    c_plus1_strongNat[(plus1flag, strongNat)] += 1
    if n >= 3:
        hc = has_ham_cycle(n, E)
        es = exists_order_strong(n, E)
        c_hamcyc_existstrong[(hc, es)] += 1
print("  [conj1] Alcuin=tau+1  vs  T_G strong (natural order):")
for k in sorted(c_plus1_strongNat):
    print(f"     plus1={k[0]}, strongNat={k[1]}: {c_plus1_strongNat[k]}")
# conj1': does Alcuin=tau+1 imply T_G NON-strong under ALL orders?
viol_conj1 = 0
for (n, E) in allreps:
    t = tau_number(n, E); a = alcuin_number(n, E)
    if a == t + 1 and exists_order_strong(n, E):
        viol_conj1 += 1
print(f"  [conj1'] #(Alcuin=tau+1 AND exists-order-strong) = {viol_conj1}  (0 => +1 implies never-strong)")
print("  [conj2] G has Ham-cycle  vs  exists-order T_G strong (n>=3):")
for k in sorted(c_hamcyc_existstrong):
    print(f"     hamcycle={k[0]}, exists-strong={k[1]}: {c_hamcyc_existstrong[k]}")

print("\nDONE.")
