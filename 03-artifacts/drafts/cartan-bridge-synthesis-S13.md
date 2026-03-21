# The Cartan Bridge: Deep Investigation
## Session kind-pasteur-2026-03-21-S13

### Executive Summary

This session investigated the Cartan decomposition gl(n,R) = so(n) + p + R as a bridge between tournament theory and transformer attention analysis. Six investigations produced five theorems and three conjectures.

---

## Theorem 1: S_T^2 Characterizes Double Regularity

**Statement.** For a tournament T on n vertices with signed adjacency S_T:
- S_T^2 = -nI + J if and only if T is doubly regular (DRT).
- For DRT: all off-diagonal entries of S_T^2 equal 1.
- For non-DRT regular: off-diagonal entries of S_T^2 take values from {5-4B : B in {0,1,2,...}} where B is the common out-neighbor count per arc.

**Verified:** n=3 (all regular are DRT), n=5 (no DRT exists since d=0.5), n=7 (only Paley among regulars is DRT, 240/2640).

**Proof sketch:** (S_T^2)[i,j] = 2(p(i,j)+p(j,i)) - (n-2) for i != j, where p(i,j) = #{2-paths i->k->j}. For regular tournaments with score (n-1)/2, this equals 5-4B where B is the mutual out-neighbor count. Uniformity (B=const) is the DRT condition.

---

## Theorem 2: Paley Minimal Polynomial

**Statement.** For the Paley tournament T_p (p prime, p = 3 mod 4):
- The minimal polynomial of S_T is x(x^2 + p) = 0.
- Eigenvalues: {0} (multiplicity 1), {+/-i*sqrt(p)} (multiplicity (p-1)/2 each).
- dim(Alg(S_T)) = 3.

**Verified:** p = 3, 7, 11.

For non-Paley regular tournaments at n=7: minimal polynomial has degree 7, dim(Alg) = 7.

---

## Theorem 3: Paley Commutant Formula

**Statement.** dim(Comm(S_{T_p})) = (p^2 - 2p + 3)/2 for Paley tournament T_p.

**Values:** p=3: 3, p=7: 19, p=11: 51, p=23: 243=3^5.

**Proof:** Comm dim = sum d_k^2 where d_k are eigenspace dimensions. For Paley: d = (1, (p-1)/2, (p-1)/2). So comm = 1 + 2*((p-1)/2)^2 = (p^2-2p+3)/2.

**Verified:** p = 3, 7, 11.

---

## Theorem 4: Power Alternation

**Statement.** For any tournament T:
- S_T^k is antisymmetric (in so(n)) for k odd.
- S_T^k is symmetric (in p) for k even.

**Consequence:** Powers of S_T alternate perfectly between Cartan sectors. The "directed" and "undirected" information exchange with each matrix multiplication.

**Proof:** (S^k)^T = (S^T)^k = (-S)^k = (-1)^k S^k.

---

## Theorem 5: First Distinguishing Power

**Statement.** Among the three isomorphism classes of regular tournaments on 7 vertices:
- tr(S^k) = 0 for all odd k (from antisymmetry).
- tr(S^2) = -42 for all three classes (universal second moment).
- tr(S^4) FIRST distinguishes: H=189 -> 294, H=175 -> 742, H=171 -> 486.

The off-diagonal entries of S^2 also distinguish Paley (all entries 1) from non-Paley ({-3, 1, 5}).

---

## Conjecture 1: Spectral Flatness Principle

**Statement.** Among regular tournaments on n vertices, the one with MINIMUM tr(S^4) (flattest eigenvalue distribution) has MAXIMUM H.

**Verified:** n=3, 5, 7.

**Interpretation:** Spectral flatness = all eigenvalue magnitudes equal = maximum symmetry = maximum Hamiltonian path count. The Paley tournament is the flattest.

---

## Conjecture 2: Commutant Maximality

**Statement.** Among all tournaments on n vertices, the Paley tournament T_p has the LARGEST commutant dimension dim(Comm(S_T)).

**Partial evidence:** At n=7, Paley comm=19 vs generic comm=7. At n=5, all comm=5 (no distinction). Needs checking at n=9 and beyond.

---

## Conjecture 3: The Cartan Bridge for Attention

**Statement.** For transformer attention matrices, the spectral diversity of the antisymmetric part A_anti = (A-A^T)/2 predicts computational richness:
- Paley-like attention (minimal spectral diversity = 2 distinct magnitudes): maximally symmetric directed attention.
- Random attention (maximal spectral diversity = ceil(n/2) magnitudes): generic, unstructured.
- Trained attention: intermediate — structured but not maximally symmetric.

**Verified empirically:** Paley has spectral diversity 2 vs random's ceil(n/2) for n=5,7,11.

---

## Key Quantities

| Invariant | What it measures | Determined by? |
|-----------|-----------------|---------------|
| tr(S^2) | Total "angular momentum" | n only (universal) |
| tr(S^4) | Eigenvalue kurtosis | Spectrum |
| Off-diag(S^2) | DRT deviation | 2-path structure |
| dim(Alg(S_T)) | Complexity of tournament algebra | Minimal polynomial |
| dim(Comm(S_T)) | Symmetry of eigenspace structure | Eigenspace dimensions |
| H(T) | Hamiltonian path count | Full matrix (NOT spectrum alone) |

**The hierarchy:** tr(S^2) -> spectrum -> Alg(S_T) -> H(T) -> S_T

Each level provides more information. H requires more than the spectrum but less than the full matrix.

---

## Scripts and Output Files

- `cartan_bridge_deep.py` → `cartan_bridge_deep.out` (Investigations 1-6)
- `cartan_bridge_deep2.py` → `cartan_bridge_deep2.out` (Investigations A-E)
- `paley_commutant_theorem.py` → `paley_commutant_fast.out`

## New Hypotheses

- HYP-1710: Spectral Flatness Principle (CONFIRMED n=3,5,7)
- HYP-1711: Commutant Maximality for Paley (partial evidence)
- HYP-1712: Trained attention has lower spectral diversity than random (OPEN)

## Open Questions

1. Does the Spectral Flatness Principle hold at n=9? (Need to check regular tournaments at n=9 — computationally feasible with sampling.)
2. Is there a spectral formula for H in terms of eigenvalues + eigenvector invariants?
3. What is the minimal polynomial of non-DRT regular tournaments at n=11?
4. Does trained GPT-2 attention have intermediate spectral diversity?
