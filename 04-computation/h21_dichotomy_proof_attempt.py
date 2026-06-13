"""
h21_dichotomy_proof_attempt.py — Attempt to prove the dichotomy:
"Every cycle-rich tournament on n >= 9 vertices either has 3 disjoint 3-cycles
OR has a good vertex deletion."

Theoretical approach:
- A vertex v is "bad" if T-v is not cycle-rich.
- T-v fails cycle-rich if: (A) source/sink created, or (B) some u loses all 3-cycles.

(A) source in T-v: exists u with score(u)=1 and u→v (v is u's only out-target)
    sink in T-v: exists u with score(u)=n-2 and v→u

(B) u's only 3-cycles go through v: every 3-cycle containing u also contains v.
    This means: for every pair (a,b) with a,b != v,u: the triple {u,a,b} is NOT
    a 3-cycle. So all C(n-2,2) triples {u,a,b} (a,b not v) are transitive.
    Equivalently: in T-{v}, the vertex u has NO 3-cycle. By Key Lemma (Part J),
    u is acyclic in T-{v}.

    When can u be acyclic in T-{v}? By Key Lemma: all arcs from N-(u,T-v) go to
    N+(u,T-v). Let s = score(u) in T. If v→u, then score(u,T-v) = s, and u
    needs all C(n-2,2) triples to be transitive. If u→v, score(u,T-v) = s-1.

    For u acyclic in T-v: N-(u) → u → N+(u) and N-(u) → N+(u) (in T-v).
    This means u "separates" T-v into layers.

KEY OBSERVATION: If u is acyclic in T-v, then u's neighborhood has layered
structure. Score(u,T-v) = s' means s' vertices beaten, n-2-s' vertices beating u.
All arcs between the two groups go one way.

COUNTING: How many vertices can be "cyclically dependent" on a single vertex v?
Each such vertex u creates a partition of V\{u,v} into N+(u)\{v} and N-(u)\{v}
with all arcs going N-(u) → N+(u). This is a very rigid structure!

If two vertices u1, u2 are both cyclically dependent on v, their layered structures
must be compatible. This severely limits how many can exist simultaneously.

Author: opus-2026-03-07-S43
"""
import random
from itertools import combinations

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def is_cycle_rich(A, n):
    sc = [sum(A[i]) for i in range(n)]
    for s in sc:
        if s == 0 or s == n-1:
            return False
    for v in range(n):
        found = False
        for a in range(n):
            if a == v: continue
            for b in range(n):
                if b == v or b == a: continue
                if A[v][a] and A[a][b] and A[b][v]:
                    found = True
                    break
            if found: break
        if not found:
            return False
    return True

def cyclically_dependent_on(A, n, u, v):
    """Check if u's only 3-cycles all go through v (u is acyclic in T-v)."""
    verts = [w for w in range(n) if w != v]
    for a in verts:
        if a == u: continue
        for b in verts:
            if b == u or b == a: continue
            if A[u][a] and A[a][b] and A[b][u]:
                return False  # u has a 3-cycle not through v
    return True  # All u's 3-cycles go through v

def score_creates_source_sink(A, n, v):
    """Which vertices become source/sink when v is deleted?"""
    result = []
    for u in range(n):
        if u == v: continue
        s = sum(A[u][w] for w in range(n) if w != v and w != u)
        if s == 0 or s == n-2:
            result.append(u)
    return result

# Explore: for cycle-rich n=9 with no good deletion,
# what's the structure of cyclic dependencies?
random.seed(123)
n = 9
trials = 10000000
found = 0

print(f"=== Cyclic Dependency Analysis at n={n} ===")

for trial in range(trials):
    A = random_tournament(n)
    if not is_cycle_rich(A, n):
        continue

    # Quick check: any good deletion?
    has_good = False
    for v in range(n):
        # Check T-v is cycle-rich
        verts = [u for u in range(n) if u != v]
        n1 = n - 1
        sc_Tv = [sum(A[u][w] for w in verts if w != u) for u in verts]

        # Source/sink check
        if any(s == 0 or s == n1-1 for s in sc_Tv):
            continue

        # All-in-3-cycle check
        all_ok = True
        for u in verts:
            cyc_found = False
            for a in verts:
                if a == u: continue
                for b in verts:
                    if b == u or b == a: continue
                    if A[u][a] and A[a][b] and A[b][u]:
                        cyc_found = True
                        break
                if cyc_found: break
            if not cyc_found:
                all_ok = False
                break
        if all_ok:
            has_good = True
            break

    if has_good:
        continue

    found += 1
    # Analyze cyclic dependencies
    dep_matrix = [[False]*n for _ in range(n)]
    for v in range(n):
        for u in range(n):
            if u == v: continue
            dep_matrix[u][v] = cyclically_dependent_on(A, n, u, v)

    ss_matrix = [score_creates_source_sink(A, n, v) for v in range(n)]

    # For each v, which vertices are cyclically dependent on v?
    dep_counts = [sum(dep_matrix[u][v] for u in range(n) if u != v) for v in range(n)]

    # 3-cycle sets
    cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    cycles.append((a,b,c))

    # Max matching (greedy)
    used = set()
    mm = 0
    for c in cycles:
        if not any(x in used for x in c):
            used.update(c)
            mm += 1

    sc = [sum(A[i]) for i in range(n)]

    if found <= 10:
        print(f"\n--- No-good-deletion example #{found} ---")
        print(f"  Scores: {sc}")
        print(f"  t3={len(cycles)}, mm={mm}")
        print(f"  Cyclic dependency counts per vertex: {dep_counts}")
        for v in range(n):
            deps = [u for u in range(n) if u != v and dep_matrix[u][v]]
            ss = ss_matrix[v]
            if deps or ss:
                print(f"  Del v={v}: deps={deps}, source/sink={ss}")

    if (trial + 1) % 2000000 == 0:
        print(f"  Progress: {trial+1}, found {found} no-good-deletion")

print(f"\n=== Found {found} no-good-deletion cycle-rich tournaments in {trials} trials ===")
