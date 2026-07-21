# THM-1745: the leaf-graded arborescence filtration — where {7,21} lives, what transfers, and how each hole dies

**Status:** PROVED (L1, L3 — short proofs below) + VERIFIED-EXACT (census
n ≤ 7, all 530 iso classes; L1 on all 80 classes with |Aut| > 1; L3 on 840
evaluations; L2 classic, re-verified n ≤ 7). No completion claims.
**Author:** boxeph-2026-07-20-S184 (HYP-8555)
**Owner:** "see how 7,21 forbiddenness shows up in arborescences instead of
hamiltonian paths, there should be some subtle effect on an aspect of its
structure."
**Builds on (confirm, not compete):** death-star S70 (h-spectrum =
odds ∖ {7,21}, monoid view), S71 (raw arborescence count does NOT inherit
the hole), mac-mini THM-1460 (arborescences = determinantal relaxation;
ordinal-sum det-shift law).

## 0. The object and the one-line answer

Grade the spanning out-arborescences of T rooted at r by LEAF COUNT
(leaf = tree-out-degree 0):
  A_{T,r}(x) = Σ_ℓ c_{r,ℓ} x^ℓ,   c_ℓ := Σ_r c_{r,ℓ}.
Then c₁ = h (a 1-leaf spanning out-tree IS a Hamiltonian path) and
A(1) = a_r (the matrix-tree determinant). The filtration interpolates
death-star's two endpoints, and the answer to the owner's question is
threefold:

1. **Where the hole lives:** {7,21} is confined to the ℓ = 1 stratum
   because that is the UNIQUE stratum that is simultaneously n-stable
   (small values recur at every n) and multiplicative under ordinal sum
   (L3: the (x−1)^j correction terms vanish only at the bottom
   coefficient). Every stratum above has growing band minima
   (c₂-min = 2^{n−1} − n for n ≥ 5, by L2), so no finite forbidden set can
   exist there — death-star's S71 contrast, now resolved per-stratum.
2. **What transfers (the subtle effect):** the hole's MECHANISM — the free
   Aut-action underlying |Aut| | h — extends exactly one stratum up:
   **|Aut(T)| divides c₂ for every tournament** (L1; new). The freeness
   depth is p_min(|Aut|); tournament groups are odd, so strata 1 AND 2 are
   always free, and no further in general.
3. **How each hole dies:** 21 EVAPORATES at the first grading step — it is
   attained within-band by both c₂ and B(2) at n = 5 (L4) — consistent
   with 21 = 3·7 being merely 7's composite casualty (death-star S70).
   7 never gets a window above stratum 1: it falls in the inter-band gaps
   (c₂ bands [3,6] then [11,35]…; B(2) bands [5,10] then [16,55]…), so its
   absence upstairs is arithmetic, not structural. "7 is forbidden" is a
   strictly bottom-stratum phenomenon.

## 1. (L1) The free-action depth theorem  [PROVED + verified 80/80]

**Theorem.** Let σ ∈ Aut(T) have prime order p. Every σ-fixed spanning
out-arborescence (any root) has at least p leaves. Consequently Aut(T)
acts freely on the set of ℓ-leaf rooted arborescences whenever
ℓ < p_min(|Aut(T)|), and
  **|Aut(T)| divides c_ℓ for every ℓ < p_min(|Aut(T)|).**
Since tournament automorphism groups have odd order, p_min ≥ 3:
  **|Aut(T)| divides both c₁ = h (known) and c₂ (new), always.**

*Proof.* Let A be σ-fixed with root r. σ maps A to an arborescence rooted
at σ(r); A fixed forces σ(r) = r, and σ then preserves depth in A. If
σ ≠ id, choose v moved by σ at maximal depth. The orbit of v has size p.
The subtrees below v, σv, …, σ^{p−1}v are pairwise disjoint: if
σ^i v ∈ Subtree(σ^j v) with i ≠ j then depth(σ^i v) > depth(σ^j v),
contradicting depth preservation. Each subtree contains a leaf of A
(its own deepest vertex), giving ≥ p leaves. For freeness: if σ ≠ id
fixes an ℓ-leaf arborescence, take any prime q | ord(σ); the power of
order q also fixes it, so ℓ ≥ q ≥ p_min. Orbits of the free action have
size |Aut|, giving the divisibility (per-root version with Aut_r
verbatim). ∎

Verified: all 80 iso classes with |Aut| > 1, n = 3..7, zero violations.
Visible consequence: at n = 4 the h = 3 class (C₃ block ⊕ 1) has
c₂ = 3 < 4 = c₂(TT₄) — the divisibility law undercuts the transitive
tournament, which is otherwise the c₂-minimizer (n = 5, 6 verified).

## 2. (L2) The Eulerian pole  [classic; re-verified n ≤ 7]

A_{TT_n, source}(x) = Σ_ℓ A(n−1, ℓ−1) x^ℓ (Eulerian numbers), via the
increasing-tree ↔ permutation bijection (arborescences of the transitive
tournament rooted at the source = increasing trees on [n]). In particular
c₂(TT_n) = 2^{n−1} − n and A(1) = (n−1)! (mac-mini's transitive value).
The filtration's two poles: transitive ↔ Eulerian polynomial,
bottom stratum ↔ the {7,21}-monoid.

## 3. (L3) The graded ordinal-sum law  [PROVED + verified 840/840]

Let G_S(x,t) = Σ_{spanning out-forests F of S} x^{leaves(F)} t^{comps(F)}.
For r ∈ T:
  **A_{T⊕S, r}(x) = Σ_ℓ c_{r,ℓ}(T) Σ_{j≤ℓ} C(ℓ,j) (x−1)^j G_S(x, n_T − j).**

*Proof.* In T ⊕ S every S-vertex's parent is either in S or ANY T-vertex,
so an arborescence of T ⊕ S rooted at r ∈ T = (an arborescence A of T
rooted at r) + (a spanning out-forest F of S) + (an independent choice of
attachment T-vertex for each component of F). A leaf of A stays a leaf
iff no component attaches to it; leaves of F stay leaves. For fixed A
with leaf set L and fixed F with c components,
Σ_{attachments} x^{#unhit leaves of A} = Σ_{J⊆L} (x−1)^{|J|}(n_T−|J|)^c
(inclusion–exclusion on "all of J unhit"). Summing against x^{leaves(F)}
over F gives the law. ∎

Specializations: at [x¹] the only surviving configuration is
path × path × attach-at-the-unique-leaf: h-multiplicativity, the {7,21}
monoid (death-star S70). At x = 1 only j = 0 survives:
A(1) = a_r(T)·G_S(1, n_T) = a_r(T)·det(n_T I + L_in(S)) by the
matrix-forest theorem — exactly mac-mini THM-1460(D). The two known
composition laws are the bottom and top faces of ONE graded law, and the
(x−1)-adic corrections measure precisely how the monoid structure (hence
the possibility of a finite forbidden set) is destroyed as the grading
leaves the bottom stratum.

## 4. (L4) The window facts  [census]

n = 4: c₂ ∈ {3,4,5,6}, B(2) ∈ {5,6,9,10}. n = 5: c₂ ∈ {11,…,35} with
21 ATTAINED; B(2) ∈ {16,…,55} with 21 ATTAINED. 7 appears in neither
filtration at any computed level above ℓ = 1 — always inter-band.
(B(k) = arbs with all out-degrees ≤ k: B(1) = h, B(n−1) = a; the maxdeg
ladder shows the same banding.)

## 5. Complexity remark

The filtration also interpolates complexity: c₁ = h is #P-hard, A(1) is
a determinant (poly); the childless inclusion–exclusion computes any c_ℓ
in O(2^n poly) and the TOP strata in poly time. Where hardness sets in
along ℓ is open (c₂'s complexity unknown).

## 6. Files

- `04-computation/leaf_graded_arborescences_boxeph_S184.py` + frozen
  `.out` (census n ≤ 7: inclusion–exclusion cross-checked against direct
  enumeration on ALL classes at n = 4,5; bottom stratum ≡ h asserted on
  all 530 classes; c₂ spectra; (h,c₂) joint table; signed A(−1) — no law
  found, recorded negative; maxdeg ladder n ≤ 6).
- `04-computation/leaf_filtration_laws_boxeph_S184.py` + frozen `.out`
  (L1 80/80; L2 exact rows n ≤ 7 + minimizer anomaly; L3 840/840;
  L4 windows).
