# THM-1760: the tournament pencil — one determinant whose faces are the repo's worlds, with a complement functional equation, an ordinal-sum law, and joint completeness at n = 7

**Status:** PROVED (expansion, faces, functional equation, diagonal
collapse, ⊕-law — one-line proofs below) + VERIFIED-EXACT (1152 + 486 +
144 + 6400 evaluations, zero mismatches; census n = 7 complete). The
log-concavity law is a CONJECTURE (530/530, no proof). No completion claims.
**Author:** boxeph-2026-07-20-S185 (HYP-8560)
**Owner:** "find more unifying laws, explore around thoroughly."

## 0. The object

  M_T(t,u,v) := det( tI + u·D_in − v·A ).

**Expansion (permutation expansion, proved; verified 1152/1152):**
  M_T = Σ_{L: vertex-disjoint directed cycle sets} (−1)^{#cycles(L)}
        v^{|V(L)|} Π_{w∉L} (t + u·indeg(w)).

## 1. The three faces (all verified exactly)

- **v = 0 — the SCORE/CUT face:** M = Π_j (t + u·indeg(j)). Pure in-score
  product; the out-score face is the v = 0 face of T^op (consistent via §2).
- **u = 0 — the SPECTRAL/CYCLE face:** M = char poly of A — THM-506's
  signed cycle-packing face.
- **u = v — the FOREST ray:** M(t,1,1) = det(tI + L_in) = Σ_F t^{comps(F)}
  (matrix-forest); Kirchhoff at the t¹ coefficient gives Σ_r a_r; t = 0
  gives det(L_in) = 0 (column sums vanish).

The triangle's cut ⊕ cycle decomposition is realized as the two AXES of one
determinant, with the tree world on the diagonal. H is NOT a face — the
unsigned odd-cycle world stays outside every determinant (THM-505/506
boundary) — consistent with H's #P status vs the pencil's poly-time.

## 2. The complement functional equation (proved; 486/486 + 144/144)

With A^T = J − I − A and D_in(T^op) = (n−1)I − D_in, the matrix determinant
lemma gives, with t̃ = t + u(n−1) + v:
  **M_{T^op}(t,u,v) = M_T(t̃, −u, −v) − v · S_T(t̃, −u, −v),**
  S_T(α,β,γ) := 𝟙ᵀ adj(αI + βD_in − γA) 𝟙  (the all-pairs adjugate sum —
  Chaiken's doubly-rooted forest object).
The Z₂/complement symmetry FORCES the 2-rooted object into the picture —
except on the forest diagonal, where column sums collapse it
(S(α,β,β) = n·M(α,β,β)/α), yielding the classic complement-spectrum law
  det(tI + L_in(T^op)) · (t + nu) = det((t+nu)I − uL_in(T)) · t
as the diagonal specialization (verified 144/144). Off-diagonal, the
correction is genuinely new data: the pencil pairs (M, S) under complement,
resonating with the repo's SC/NS pair structure.

## 3. The ordinal-sum law (proved in one line; 6400/6400)

S-vertices gain n_T in-degree and the pencil matrix is block-triangular on
T ⊕ S, so:
  **M_{T⊕S}(t,u,v) = M_T(t,u,v) · M_S(t + u·n_T, u, v).**
The pencil is multiplicative under ⊕ WITH AN ARGUMENT SHIFT — one law
subsuming: score-face multiplicativity (shift = in-degree gain), spectral
multiplicativity (u = 0: shift invisible), and mac-mini THM-1460(D)'s
det-shift (u = v = 1). H's multiplicativity is NOT a specialization (H is
not a face) — the two monoids are genuinely different unifications, meeting
only in THM-1745's graded interpolation.
**Coning corollary:** M_{1⊕T'}(t,u,v) = t · M_{T'}(t + u, u, v).

## 4. Fingerprint strength and TRANSVERSALITY at n = 7 (census, exact)

- The pencil (full coefficient information, 8×8×8 grid) separates
  **443 of 456** classes: 11 resistant groups (9 pairs, 2 triples), all
  adjacency-cospectral within group.
- **It splits ALL FOUR of klein-THM-1750's resistant pairs** (same spec A,
  Σa, arborescence vector {a_r}, H — their deepest wall; one is actually a
  triple, also fully split).
- klein's vector {a_r} splits 9 of the pencil's 11 groups.
- **JOINT (pencil, {a_r}) residue: exactly 2 pairs** — classes [280,284]
  and [537,538], both with {a_r} = (0,…,0,Σa): all arborescence mass at one
  root, i.e. T = 1 ⊕ T′ — the walls are CONED from n = 6 pencil collisions
  (coning corollary: M_{1⊕T'} determines only M_{T'}(t+u,u,v)).
- Within both pairs H differs (37 vs 33):
  **(pencil, {a_r}, H) is a COMPLETE invariant at n = 7 (456/456).**
  The poly-time pair (pencil, vector) leaves exactly the two coned walls;
  the #P invariant H finishes.

## 5. Census laws of the leaf polynomials (S184 objects)

- **LOG-CONCAVITY (conjecture):** the leaf-stratum sequence (c₁,…,c_{n−1})
  is log-concave — hence unimodal — for ALL 530 iso classes n = 3..7
  (0 violations). Open: proof; real-rootedness status.
- **DEPTH TIGHTNESS (THM-1745 L1 is sharp):** witnesses at n = 5 with
  |Aut| = 3 and 3 ∤ c₃ (c = [3,12,10,1] and [15,27,7,0]) — the free-action
  divisibility |Aut| | c_ℓ genuinely stops at ℓ < p.

## 6. Files

- `04-computation/tournament_pencil_boxeph_S185.py` + `.out` (P1–P4).
- `04-computation/pencil_vector_joint_residue_boxeph_S185.py` + `.out`
  (transversality, joint residue, completeness).
- `04-computation/pencil_sum_law_leafpoly_analytics_boxeph_S185.py` +
  `.out` (⊕-law 6400/6400; log-concavity 530/530; tightness witnesses).
