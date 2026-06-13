# THM-017: B(Li, Rj) = B(Lj, Ri) — Full Proof

**Type:** Theorem
**Certainty:** 5 — PROVED
**Status:** PROVED (induction on |W|, using THM-016)
**Added by:** opus-2026-03-05-S4
**Tags:** #even-odd-split #signed-adjacency #ocf #proof #induction

---

## Statement

For any weighted tournament T on n vertices with arc i→j, define the alternating subset convolution:

  B(Li, Rj) = sum_{S ⊆ W} (-1)^|S| Li(S ∪ {i}) · Rj({j} ∪ (W\S))

where W = V\{i,j}, Li(A) = h_end(T, A, i) (Ham paths on A ending at i),
Rj(A) = h_start(T, A, j) (Ham paths on A starting at j).

**Theorem:** B(Li, Rj) = B(Lj, Ri) for all tournaments T and all arcs i→j.

This holds as a **polynomial identity** in the arc variables.

---

## Proof

**By induction on m = |W| = n − 2.**

**Base cases:** m = 0 (n=2): B = 1 on both sides. m = 1 (n=3): follows from T[a][b]+T[b][a]=1.

**Inductive step:** Assume B(Li,Rj) = B(Lj,Ri) for all tournaments with fewer than m internal vertices.

### Setup

Reparametrize the interface arcs: for each w ∈ W, define
- s_w = 1 − p_w − q_w (where p_w = T[w][i], q_w = T[j][w])
- d_w = p_w − q_w

The difference D = B(Li,Rj) − B(Lj,Ri) is the change under σ: p_w ↦ 1−q_w, q_w ↦ 1−p_w.
This negates all s_w while preserving all d_w and internal arcs.

### D is linear in s

Each monomial in B uses at most one p_w (from the left factor Li) and one q_u (from the right factor Rj), contributing at most degree 2 in s. The bracket identity shows:

  p_u · q_w − (1−q_u)(1−p_w) = −[s_u(1−d_w) + s_w(1+d_u)] / 2

Since σ negates s, evenness of B in s is equivalent to D = 0. The s-degree-2 terms are automatically even. So D is LINEAR in s:

  D = sum_{v ∈ W} α_v · s_v

### α_v = 0 for each v

The coefficient α_v = ∂D/∂s_v|_{s=0} decomposes into contributions from each subset S:

**Part (a) — d-dependent:** For each u ≠ v, the d_u-contribution pairs subsets where u ∈ S (giving factor −(1+d_u)/2) with subsets where u ∈ R (giving factor −(1−d_u)/2). The difference left(u,v) − right(u,v) equals B(Li,Rj) − B(Lj,Ri) on the sub-tournament W\{u,v} of size m−2. By the **induction hypothesis**, this vanishes.

**Part (b) — d-independent:** The remaining terms (boundary + constant part of intermediate terms) reduce to:

  sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v) + (-1)^m h_end(W, v) = 0

This is exactly **THM-016** (Hamiltonian Path Alternating Sum Identity), proved by induction on m.

### Conclusion

Since D is linear in s and all coefficients vanish, D ≡ 0. Therefore B(Li,Rj) = B(Lj,Ri). QED.

---

## Verification Record

| n | m=|W| | Trials | B(Li,Rj)=B(Lj,Ri) | α_v = 0 | left=right | THM-016 |
|---|-------|--------|-------------------|---------|------------|---------|
| 3 | 1 | 50 | 2e-16 | 4e-10 | (trivial) | 0 |
| 4 | 2 | 50 | 2e-15 | 4e-09 | 0 | 2e-16 |
| 5 | 3 | 50 | 1e-14 | 3e-08 | 2e-16 | 2e-15 |
| 6 | 4 | 50 | 6e-14 | 6e-08 | 9e-16 | 4e-15 |
| 7 | 5 | 20 | 3e-13 | — | — | 3e-14 |

Code: `04-computation/q009_full_proof.py`

---

## Significance

This is a NECESSARY CONDITION for OCF but NOT known to be sufficient on its own.
The even-odd split says the alternating sum vanishes:
  sum_S (-1)^|S| Delta(S,R) = 0
while OCF requires the unsigned sum to equal delta_I:
  sum_S Delta(S,R) = delta_I

**CAVEAT (kind-pasteur-S8):** The even-odd split is a CONSEQUENCE of OCF, not
equivalent to it. Proving OCF for all n requires additionally proving delta_H = delta_I
(currently verified for n ≤ 8 via THM-015).

### What this proves and what remains

**Proved for all n:**
1. THM-016 (Claim B path identity): by induction on m
2. THM-017 (this theorem): by induction on m using THM-016
3. The signed adjacency identity / even-odd split

**Proved for n ≤ 8:**
4. THM-015: delta_H = delta_I (polynomial identity, exhaustive verification)
5. OCF: H(T) = I(Ω(T), 2)

**Still open:**
6. OCF for all n (requires proving delta_H = delta_I for all n, or finding a new route)
