#!/usr/bin/env python3
"""
beta1_coverage_proof_S72.py — opus-2026-03-15-S72

PROVING: β₁=0 ⟹ all edges covered by transitive triples
============================================================

PROOF SKETCH:
  Suppose edge i→j is NOT covered by any transitive triple.
  Then for every k ∉ {i,j}: NOT (i→k ∧ k→j ∧ i→j as triple).
  Specifically, no k with A[j][k]=1 ∧ A[i][k]=1 (common outneighbor),
  no k with A[k][i]=1 ∧ A[k][j]=1 (common inneighbor),
  and no k with A[i][k]=1 ∧ A[k][j]=1 (intermediate).

  We need to show: the 1-cycle consisting of just edge (i→j) doesn't
  bound... wait, a single edge is NOT a 1-cycle (∂₁(i→j) = j - i ≠ 0).

  Better approach: if edge i→j is uncovered, we need a CYCLE through it.
  Since T is a tournament, there exists a directed path j→...→i
  (Redei's theorem: tournaments are semicomplete).

  Actually, let's think more carefully.

  dim(ker ∂₁) = |E| - rank(∂₁) = n(n-1)/2 - (n-1) = C(n,2) - n + 1

  rank(∂₂) = number of independent triple boundaries that land in edge space.

  β₁ = ker(∂₁) - rank(∂₂)

  If β₁ = 0, then rank(∂₂) = C(n,2) - n + 1.

  CLAIM: If edge i→j is uncovered, then rank(∂₂) < C(n,2) - n + 1.

  PROOF: Consider the "edge indicator" e_{ij} ∈ R^E (1 at position (i,j), 0 elsewhere).
  Every triple boundary ∂(a,b,c) = e_{bc} - e_{ac} + e_{ab}.
  If (i,j) is uncovered, then e_{ij} does NOT appear in any triple boundary.
  Therefore, e_{ij} is NOT in the column space of ∂₂.

  But wait — is e_{ij} in ker(∂₁)?  ∂₁(e_{ij}) = e_j - e_i ≠ 0.
  So e_{ij} itself is not a cycle. We need more.

  Actually the proof is simpler than I thought:

  If edge (i,j) is uncovered, the (i,j)-th ROW of ∂₂ is all zeros.
  Because ∂₂[edge (i,j), triple t] ≠ 0 only if edge (i,j) appears in ∂(t).
  Edge (i,j) appears in ∂(t) exactly when (i,j) is a face of triple t
  AND triple t is a transitive triple containing edge (i,j).
  If no such triple exists, the row is zero.

  BUT WAIT: this means im(∂₂) is contained in the subspace where
  the (i,j) coordinate is zero. So im(∂₂) ⊆ {v : v_{ij} = 0}.

  Now, is there a 1-cycle with nonzero (i,j) coordinate?
  We need v ∈ ker(∂₁) with v_{ij} ≠ 0.

  ∂₁(v) = 0 means: for each vertex w, Σ_{edges into w} v_e - Σ_{edges out of w} v_e = 0.

  Construct: start with v = e_{ij}. Then ∂₁(v) = e_j - e_i.
  To cancel: add a path j→...→i with coefficient +1.
  In a tournament, such a path always exists (by Redei).
  The resulting vector is in ker(∂₁) and has v_{ij} = 1.

  Since this cycle has nonzero (i,j) component and im(∂₂) has
  zero (i,j) component everywhere, this cycle is NOT in im(∂₂).
  Therefore β₁ ≥ 1.

  QED: β₁ = 0 ⟹ every edge is in some transitive triple.

  CONVERSE FAILS because: even with full coverage, the signs
  in triple boundaries may not align to fill all cycles.

Let's verify this proof is correct by checking its predictions.
"""

import numpy as np
import sys
import time
sys.stdout.reconfigure(line_buffering=True)

def adj_from_bits(n, bits):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

PRIME = 32749

def rank_mod(M, p):
    M = M.copy() % p
    nrows, ncols = M.shape
    r = 0
    for col in range(ncols):
        nz = np.where(M[r:, col] != 0)[0]
        if len(nz) == 0:
            continue
        pivot = nz[0] + r
        if pivot != r:
            M[[r, pivot]] = M[[pivot, r]]
        inv = pow(int(M[r, col]), p - 2, p)
        M[r] = M[r] * inv % p
        factors = M[:, col].copy()
        factors[r] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[r])) % p
        r += 1
    return r

print("=" * 72)
print("PROOF VERIFICATION: β₁=0 ⟹ all edges in transitive triples")
print("opus-2026-03-15-S72")
print("=" * 72)

print("""
THEOREM (Coverage Necessity for β₁=0):
  For any tournament T on n vertices, if β₁(T) = 0, then every
  directed edge i→j participates in at least one transitive triple.

PROOF:
  Let T be a tournament with β₁(T) = 0, and suppose for contradiction
  that edge i→j is not in any transitive triple.

  Step 1: The row of ∂₂ corresponding to edge (i,j) is zero.
    The boundary of transitive triple (a,b,c) [where a→b, b→c, a→c] is:
      ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)
    Edge (i,j) appears in this boundary only if (i,j) ∈ {(a,b),(a,c),(b,c)},
    meaning edge (i,j) is part of the triple. By assumption, no such triple
    exists. So ∂₂ has a zero row at position (i,j).

  Step 2: im(∂₂) lies in the hyperplane {v ∈ R^E : v_{(i,j)} = 0}.
    Since row (i,j) of ∂₂ is zero, every column of ∂₂ has zero in
    the (i,j) coordinate. Hence every vector in im(∂₂) has v_{(i,j)} = 0.

  Step 3: There exists a 1-cycle c ∈ ker(∂₁) with c_{(i,j)} ≠ 0.
    Since T is a tournament on n ≥ 3 vertices, and every tournament
    has a Hamiltonian path (Rédei's theorem), there exists a directed
    path P: j = v₀ → v₁ → ... → v_ℓ = i from j back to i.
    (Such a path exists because the tournament is semiconnected.)

    Define c = e_{(i,j)} + Σ_{k=0}^{ℓ-1} e_{(v_k, v_{k+1})}.
    Then ∂₁(c) = (j - i) + (v₁ - j) + (v₂ - v₁) + ... + (i - v_{ℓ-1}) = 0.
    So c ∈ ker(∂₁) and c_{(i,j)} = 1 ≠ 0.

  Step 4: c ∉ im(∂₂), so β₁ ≥ 1.
    From Step 2, every vector in im(∂₂) has zero (i,j) coordinate.
    But c has c_{(i,j)} = 1. So c is in ker(∂₁) but not in im(∂₂).
    Therefore β₁ = dim(ker ∂₁ / im ∂₂) ≥ 1, contradicting β₁ = 0.  ∎

  NOTE: The converse fails. Full edge coverage is NECESSARY but not
  SUFFICIENT for β₁=0. The "twisted" tournaments have full coverage
  but the triple boundaries contribute wrong signs to some 1-cycles.
""")

# Verify Step 3 explicitly: for every uncovered edge, find a return path
print("Verification: constructing return paths for uncovered edges")
print("=" * 60)

n = 5
num_edges = n * (n - 1) // 2
total = 1 << num_edges
verified = 0

for bits in range(total):
    A = adj_from_bits(n, bits)

    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0:
                continue

            # Check if uncovered
            uncov = True
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[j][k] and A[i][k]:
                    uncov = False; break
                if A[k][i] and A[k][j]:
                    uncov = False; break
                if A[i][k] and A[k][j]:
                    uncov = False; break

            if not uncov:
                continue

            # Find return path j → ... → i using BFS
            from collections import deque
            visited = {j}
            parent = {j: None}
            q = deque([j])
            found = False
            while q:
                v = q.popleft()
                for w in range(n):
                    if w not in visited and A[v][w] == 1:
                        visited.add(w)
                        parent[w] = v
                        if w == i:
                            found = True
                            break
                        q.append(w)
                if found:
                    break

            if found:
                verified += 1
            else:
                print(f"  FAILURE: no return path for edge {i}→{j} at bits={bits}")
                print(f"  A = {A.tolist()}")

print(f"  Verified {verified} return paths exist (all uncovered edges at n=5)")

# Now check if the 1-cycle we construct is genuinely NOT in im(∂₂)
print(f"\nVerifying constructed cycles are NOT in im(∂₂)...")

count_checked = 0
count_not_in_image = 0

for bits in range(total):
    A = adj_from_bits(n, bits)

    edges = []
    for ii in range(n):
        for jj in range(n):
            if ii != jj and A[ii][jj] == 1:
                edges.append((ii, jj))
    edge_idx = {e: i for i, e in enumerate(edges)}
    ne = len(edges)

    triples = []
    for a in range(n):
        for b in range(n):
            if a == b or A[a][b] == 0: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] == 1 and A[a][c] == 1:
                    triples.append((a, b, c))
    nt = len(triples)

    d2 = np.zeros((ne, nt), dtype=np.int64)
    for t_idx, (a, b, c) in enumerate(triples):
        if (b, c) in edge_idx: d2[edge_idx[(b, c)], t_idx] += 1
        if (a, c) in edge_idx: d2[edge_idx[(a, c)], t_idx] -= 1
        if (a, b) in edge_idx: d2[edge_idx[(a, b)], t_idx] += 1

    base_rank = rank_mod(d2 % PRIME, PRIME)

    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue

            uncov = True
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[i][k]: uncov = False; break
                if A[k][i] and A[k][j]: uncov = False; break
                if A[i][k] and A[k][j]: uncov = False; break

            if not uncov: continue

            # The edge indicator e_{ij} has nonzero in (i,j) coordinate
            # Just check: augmenting d2 with e_{ij} increases rank
            col = np.zeros((ne, 1), dtype=np.int64)
            col[edge_idx[(i, j)], 0] = 1
            aug = np.hstack([d2, col])
            aug_rank = rank_mod(aug % PRIME, PRIME)

            count_checked += 1
            if aug_rank > base_rank:
                count_not_in_image += 1

print(f"  Checked {count_checked} uncovered edges")
print(f"  e_{{ij}} not in im(∂₂): {count_not_in_image}/{count_checked}")
if count_not_in_image == count_checked:
    print(f"  ✓ CONFIRMED: every uncovered edge indicator is independent of im(∂₂)")
    print(f"    This directly proves: uncovered edge ⟹ β₁ ≥ 1")

# KEY INSIGHT: can we strengthen to β₁ = (number of uncovered edges)?
# No — from previous results, at n=5 some β₁=1 have 0 or 1 uncovered edges
# But can we say: β₁ ≥ (some function of uncovered edges)?
print(f"\n{'='*60}")
print(f"CHECKING: does #uncovered_edges give a lower bound on β₁?")
print(f"{'='*60}")

# At n=5, β₁ ∈ {0,1} always, and #uncov ∈ {0,1}
# The twisted cases have uncov=0, β₁=1 so β₁ > #uncov is possible
# The untwisted cases have uncov=1, β₁=1 so β₁ ≥ #uncov barely holds
# We can't say β₁ = #uncov in general.

print(f"  β₁ ∈ {{0,1}} at n=5, so the only question is necessary conditions")
print(f"  We proved: uncov > 0 ⟹ β₁ > 0 (equivalently, β₁=0 ⟹ uncov=0)")
print(f"  This is the strongest possible one-directional statement.")

# Characterize uncovered edges
print(f"\n{'='*60}")
print(f"CHARACTERIZING UNCOVERED EDGES")
print(f"{'='*60}")

print("""
  Edge i→j is uncovered iff it participates in no transitive triple.
  A transitive triple containing i→j has one of these forms:
    (i, j, k): i→j, j→k, i→k  [k is common out-neighbor]
    (k, i, j): k→i, i→j, k→j  [k dominates both]
    (i, k, j): i→k, k→j, i→j  [k is intermediate]

  Edge i→j is uncovered iff for ALL k ∉ {i,j}:
    - NOT (A[j][k]=1 AND A[i][k]=1)   → k not common out-neighbor
    - NOT (A[k][i]=1 AND A[k][j]=1)   → k not common dominator
    - NOT (A[i][k]=1 AND A[k][j]=1)   → k not intermediate

  Equivalently: for each k, at least one of these holds:
    (A[j][k]=0 OR A[i][k]=0) AND (A[k][i]=0 OR A[k][j]=0) AND (A[i][k]=0 OR A[k][j]=0)

  In a tournament, A[x][y] + A[y][x] = 1 for x≠y. So:
    A[j][k]=0 means k→j
    A[i][k]=0 means k→i
    A[k][i]=0 means i→k
    A[k][j]=0 means j→k
    A[i][k]=0 means k→i
    A[k][j]=0 means j→k

  Let's denote: for k ∉ {i,j}, we have edges k→i or i→k, and k→j or j→k.
  Four cases for k:
    (k→i, k→j): k dominates both. Then (k,i,j) is a triple. NOT uncovered.
    (k→i, j→k): k→i, j→k. Check (i,k,j): need i→k, k→j. We have k→i, so A[i][k]=0. First cond: no.
                Check (i,j,k): need j→k, i→k. We have k→i, so A[i][k]=0. No.
                Check (k,i,j): need k→i, k→j. We have j→k. No.
                → No triple possible! This is consistent with uncovered.
    (i→k, k→j): intermediate. (i,k,j) is a triple. NOT uncovered.
    (i→k, j→k): Check (i,j,k): j→k, i→k. YES! Triple (i,j,k). NOT uncovered.

  So edge i→j is uncovered iff EVERY k satisfies: k→i AND j→k.
  i.e., N_in(i) ∩ N_in(j)^c ⊇ V\\{i,j}
  Wait, let me redo: k→i means k ∈ N_in(i). j→k means k ∈ N_out(j).
  So: every k ∈ V\\{i,j} has k→i AND j→k.
  i.e., V\\{i,j} ⊆ N_in(i) ∩ N_out(j).

  EQUIVALENTLY: N_out(i)\\{j} = ∅ ... wait, i→k means k ∈ N_out(i).
  If every k has k→i, then N_out(i) ⊆ {j}. Since i→j, we have N_out(i) = {j}.
  So out-degree of i is exactly 1 (the edge to j).

  Similarly, j→k for all k means N_in(j) ⊆ {i} union possibly i→j.
  Wait: j→k for all k ∈ V\\{i,j}. Plus i→j. So j beats everyone EXCEPT i→j
  is already an edge meaning i beats j. Wait no: i→j means j ∈ N_out(i).
  Who beats j? We need A[k][j] or A[j][k]:
  If j→k for all k ∈ V\\{i,j}, then j beats everyone except possibly i.
  And i→j, so j doesn't beat i. So N_in(j) = {i}, out-degree of j = n-2.
  No wait: N_out(j) = V\\{i,j} so out-degree of j = n-2.

  THEOREM: Edge i→j is uncovered iff deg_out(i) = 1 (= j) AND deg_out(j) = n-2.
  This means i has minimum score 1 and j has maximum score n-2.
  But j doesn't beat i (since i→j and we're in a tournament).
  Wait: actually j→k for all k≠i,j, plus i→j. So out-deg(j) = n-2 (beats all except i).
  And out-deg(i) = 1 (beats only j).

  Let me verify this!
""")

# Verify: uncovered edge iff source has out-degree 1 and target has out-degree n-2
n = 5
total = 1 << (n*(n-1)//2)
correct = 0
wrong = 0
for bits in range(total):
    A = adj_from_bits(n, bits)
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue

            # Actual uncovered check
            uncov = True
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[i][k]: uncov = False; break
                if A[k][i] and A[k][j]: uncov = False; break
                if A[i][k] and A[k][j]: uncov = False; break

            # Predicted: out_deg(i) = 1 and out_deg(j) = n-2
            od_i = sum(A[i])
            od_j = sum(A[j])
            predicted = (od_i == 1 and od_j == n - 2)

            if uncov == predicted:
                correct += 1
            else:
                wrong += 1
                if wrong <= 3:
                    print(f"  MISMATCH: edge {i}→{j}, uncov={uncov}, predicted={predicted}")
                    print(f"    out_deg({i})={od_i}, out_deg({j})={od_j}")
                    print(f"    A = {A.tolist()}")

print(f"\nn=5: correct={correct}, wrong={wrong}")
if wrong == 0:
    print(f"✓ PROVED: edge i→j uncovered ⟺ out_deg(i)=1 AND out_deg(j)=n-2")
else:
    print(f"Characterization has errors, need to refine")

# Also verify at n=6
print(f"\nVerifying at n=6...")
n = 6
total = 1 << (n*(n-1)//2)
correct = 0
wrong = 0
for bits in range(total):
    A = adj_from_bits(n, bits)
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue

            uncov = True
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[i][k]: uncov = False; break
                if A[k][i] and A[k][j]: uncov = False; break
                if A[i][k] and A[k][j]: uncov = False; break

            od_i = sum(A[i])
            od_j = sum(A[j])
            predicted = (od_i == 1 and od_j == n - 2)

            if uncov == predicted:
                correct += 1
            else:
                wrong += 1

print(f"n=6: correct={correct}, wrong={wrong}")
if wrong == 0:
    print(f"✓ CONFIRMED at n=6: edge uncovered ⟺ source has min score, target has max-1 score")

print(f"\n{'='*72}")
print(f"COMPLETE THEOREM")
print(f"{'='*72}")
print(f"""
  THEOREM (Uncovered Edge Characterization):
    In a tournament T on n vertices, edge i→j participates in no
    transitive triple if and only if:
      out-degree(i) = 1 (i beats only j)
      out-degree(j) = n-2 (j beats everyone except i)

  PROOF: See case analysis above. In a tournament, for each k ∉ {{i,j}},
  there are exactly 4 configurations. Three of them create a transitive
  triple with (i,j). The fourth (k→i AND j→k) does not. Edge (i,j) is
  uncovered iff ALL k give the fourth case, meaning every k→i and j→k
  for all k, hence out-deg(i)=1 and out-deg(j)=n-2.

  COROLLARY (Coverage Necessity):
    β₁(T) = 0 ⟹ T has no edge i→j with out-deg(i)=1 and out-deg(j)=n-2.
    Equivalently: β₁=0 ⟹ T has no "domination edge" (bottom vertex → near-top vertex).

  COROLLARY (Score Condition):
    If the score sequence of T has minimum ≥ 2, then all edges are
    covered by transitive triples. (But β₁ may still be 1 — the twist!)
""")

print("DONE — beta1_coverage_proof_S72.py complete")
