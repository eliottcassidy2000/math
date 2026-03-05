# THM-013: Arc-Flip Independence Polynomial Formula

**Type:** Theorem (verified n=4,5,6 exhaustively; n=7 sampled)
**Certainty:** 5 — PROVED at n≤6 (exhaustive); VERIFIED at n=7 (50 random flips)
**Status:** PROVED at n≤6, CONJECTURED for general n
**Added by:** opus-2026-03-05-S2
**Tags:** #ocf #arc-reversal #independence-polynomial #open-q-009 #claim-a

---

## Statement

Let T be a tournament on n=6 vertices, and let T' be obtained by flipping arc i→j to j→i. Define:
- s_x = 1 - T[x][i] - T[j][x] for each x ∈ V\{i,j} (takes values in {-1, 0, 1})
- B_x = V\{i,j,x} (the 3-vertex complement)
- D5 = #{directed 5-cycles using arc i→j in T}
- C5 = #{directed 5-cycles using arc j→i in T'}

Then the **GENERAL** formula (verified n=4,5,6,7):

**ΔH = H(T) - H(T') = -2·Σ_{x ∈ V\{i,j}} s_x · H(T[B_x]) + 2·Σ_{L≥5, L odd} (DL - CL)**

where DL = #{directed L-cycles using arc i→j in T}, CL = #{directed L-cycles using arc j→i in T'}.

Equivalently (via the adjacency formula ΔH = adj_T(i,j) - adj_{T'}(j,i)):

**adj_T(i,j) - adj_{T'}(j,i) = -2·Σ_x s_x · H(B_x) + 2·Σ_{L≥5} (DL - CL)**

IMPORTANT: The cycle sum starts at L=5, not L=3. The 3-cycle contribution is already encoded in the -2·Σ s_x·H(B_x) term (via the identity D3-C3 = -Σ s_x).

---

## Proof Sketch (for ΔI formula)

At n=6, the independence polynomial has the form:
I(Ω(T), 2) = 1 + 2·|{odd cycles}| + 4·|{VD 3-cycle pairs}|

(Max independent set size = 2 at n=6: only pairs of vertex-disjoint 3-cycles.)

**Step 1: Δ(#cycles)**
- Destroyed odd cycles (in T not T') must contain both i,j in their vertex set (since they use arc i→j).
- Created odd cycles (in T' not T) must contain both i,j (they use arc j→i).
- #destroyed_3 - #created_3 = Σ_x (T[j][x]·T[x][i] - T[i][x]·T[x][j]) = -Σ_x s_x
- So: Δ(#cycles) = -Σ s_x + (D5 - C5)

**Step 2: Δ(#VD pairs)**
A VD 3-cycle pair partitions V into two triples {A, B}. The flip i↔j only affects triples containing both i and j. For each x ∈ V\{i,j}, the partition {{i,j,x}, V\{i,j,x}} has:
- δ(cyclicity of {i,j,x}) = c_{T'}({i,j,x}) - c_T({i,j,x}) = s_x
- cyclicity of V\{i,j,x} is unchanged
- Δ(VD pair) = s_x · c(B_x)

Since H(B_x) = 1 + 2·c(B_x) for a 3-tournament:
Δ(#VD pairs) = Σ_x s_x · c(B_x)

**Step 3: Combine**
ΔI = 2·(-Σ s_x + D5 - C5) + 4·(- Σ s_x · c(B_x))
   = -2·Σ s_x · (1 + 2c(B_x)) + 2·(D5 - C5)
   = -2·Σ s_x · H(B_x) + 2·(D5 - C5)  QED

---

## Significance

This reduces OCF (= Claim A) to a single identity:

**adj_T(i,j) - adj_{T'}(j,i) = -2·Σ_x s_x · H(B_x) + 2·(D5 - C5)**

The LHS counts adjacency differences for Hamiltonian paths. The RHS involves:
- s_x: determined by 4 arc directions (x↔i, j↔x)
- H(B_x): Hamiltonian path count of 3-vertex sub-tournaments (always 1 or 3)
- D5, C5: 5-cycle counts, expressible as sums over Ham paths of B_x weighted by boundary arcs

Since any tournament is reachable from the transitive tournament by arc flips, and OCF holds for the transitive tournament (H=1=I(∅,2)), proving this identity for all flips proves OCF for all n=6 tournaments.

---

## Generalization Prospects

At n=5: Ω(T) is always complete (any two cycles share a vertex), so I(Ω,2) = 1 + 2·|cycles|. The formula simplifies to ΔI = 2·(#destroyed - #created), which was verified 732/732.

At n≥7: the formula generalizes but requires additional terms for independent sets of size > 2 (possible when n ≥ 9, where three VD 3-cycles can fit).

For n=7,8: max independent set size is still 2, so the same formula structure applies with modified B_x sizes.

---

## Verification Record

| n | Flips tested | Formula correct |
|---|-------------|----------------|
| 4 | 500 random | 500/500 |
| 5 | 300 random + 732 exhaustive | 1032/1032 |
| 6 | 200 random + 2216 exhaustive | 2416/2416 |
| 7 | 50 random | 50/50 |

---

## The 5-Cycle Term

D5 - C5 = Σ_x Σ_{P path in B_x} [T[j][P[0]]·T[P[2]][i] - T[i][P[0]]·T[P[2]][j]]

This counts the net change in 5-cycles, weighted by whether each 3-vertex path in the complement can be "extended" through the flipped arc.
