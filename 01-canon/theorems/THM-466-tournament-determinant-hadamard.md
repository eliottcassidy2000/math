# THM-466 — The Tournament Determinant Floor: det(I+S) = 2^(n-1) iff local order

**Status:** PROVED (mac-mini-2026-06-10-S2) — adversarial verification PENDING (this session).
**Provenance:** mac-mini-2026-06-10-S2, session theme "Hadamard matrices x tournaments x
odd functions x simplicial geometry". Companions: THM-467 (ceiling), THM-468 (average).
Related: THM-174 (det(I+2A) = Pf^2 — a DIFFERENT determinant: J+S, not I+S), THM-442,
T777, HYP-2384.

## Setup

T a tournament on n vertices, A its 0/1 adjacency, **S = A - A^T** the skew {0,±1}
matrix (S_ij = +1 iff i→j), **M = I + S** the ±1 tournament matrix.
**Switching** at W ⊆ [n]: reverse all arcs across the cut (W, W^c); equivalently
S ↦ DSD with D = diag(±1), D_ii = -1 iff i ∈ W. Switching is an equivalence on
tournaments (the oriented analog of Seidel switching).
A tournament is **locally transitive** (a **local order**) iff every out-neighborhood
and every in-neighborhood induces a transitive subtournament; equivalently iff it
contains no **vortex** = 4-set {v} ∪ C with C a 3-cycle and v→C or C→v.

## Statement

For every tournament T on n ≥ 1 vertices:

1. **(Pfaffian expansion)** det(I+S) = Σ_{K ⊆ [n], |K| even} Pf(S[K])², where each
   Pf(S[K]) is an ODD integer (K = ∅ contributes 1). Hence det(I+S) is a positive
   integer ≥ 2^(n-1) (there are 2^(n-1) even subsets), and det(I+S) is divisible by
   2^(n-1) (classical ±1-matrix divisibility), so **d(T) := det(I+S)/2^(n-1)** is a
   positive integer.
2. **(Invariance)** d(T) is invariant under isomorphism, switching, and reversal
   (T ↦ T^op): it is an invariant of the *switching class* (oriented two-graph) and
   descends to the merged metagraph G_n/Z_2.
3. **(Floor)** The following are equivalent:
   (a) det(I+S) = 2^(n-1) (i.e. d(T) = 1);
   (b) every even principal Pfaffian minor of S equals ±1;
   (c) T is vortex-free (locally transitive, a local order);
   (d) T is a switching of a transitive tournament.

**Corollary.** The number of labeled local orders on n vertices is (2n-2)!! =
2^(n-1)(n-1)! (re-derived; matches the classical local-order count), and the d = 1
stratum of the metagraph is exactly the switching class of the principal-line origin
(the transitive class).

## Proof

**(1).** det(I+X) = Σ_K det(X[K]) over all principal minors (standard). For skew X,
det(X[K]) = 0 when |K| is odd and det(X[K]) = Pf(X[K])² when |K| is even. Pf(S[K])
= Σ_μ sgn(μ) Π s_e over perfect matchings μ of K; with all s_e = ±1 this is a sum of
(|K|-1)!! terms each ±1, so Pf(S[K]) ≡ (|K|-1)!! ≡ 1 (mod 2): odd, hence Pf² ≥ 1.
There are 2^(n-1) even subsets of [n]. Divisibility of the determinant of any n×n
±1 matrix by 2^(n-1): subtract the first row from all others; rows 2..n become
{0,±2}-rows; extract 2^(n-1)... (extract 2 from each of the n-1 rows). ∎

**(2).** det(I + DSD) = det(D(I+S)D) = (det D)² det(I+S) = det(I+S). Reversal:
det(I-S) = det((I-S)^T) = det(I+S). Isomorphism: conjugation by a permutation
matrix. Complement = reversal for tournaments, so d descends to G_n/Z_2.
Note also |Pf((DSD)[K])| = |det(D[K])|·|Pf(S[K])| = |Pf(S[K])| and switching
commutes with restriction to K: ((DSD)[K] = D[K] S[K] D[K]). ∎

**(3).** (a) ⟺ (b): by (1), det = 2^(n-1) iff every one of the 2^(n-1) odd squares
equals 1.

(b) ⟹ (c): the Pfaffian of a 4-set is Pf = s12·s34 - s13·s24 + s14·s23 ∈ {±1, ±3},
and a finite check of the four 4-tournament iso classes gives |Pf| = 3 exactly for the
two vortex classes (vertex → 3-cycle and 3-cycle → vertex), |Pf| = 1 for the
transitive and the strong 4-tournament. So a vortex 4-subset would force a term
Pf² = 9 > 1. [4-set check verified computationally; it is a 4-class finite case
analysis.]

(c) ⟹ (d): induction on n; n ≤ 3 trivial (every tournament on ≤ 3 vertices is a
switching of a transitive: verified directly). Let T be vortex-free on n vertices.
Note vortex-freeness ⟺ all 4-set |Pf| = 1 (the same finite check), which by (2) is
switching-invariant and obviously hereditary. Pick any vertex v. T - v is
vortex-free, so by induction T - v = σ_W(L) for a transitive L; replacing T by
σ_W(T) (still vortex-free) we may assume **T - v is transitive**, with order
u_1 ≻ u_2 ≻ ⋯ ≻ u_{n-1} (u_i → u_j ⟺ i < j). Let ε_i = 1 if v → u_i, else 0.

Every 3-cycle of T contains v, and equals {v, u_j, u_k} with j < k, ε_j = 1,
ε_k = 0 (a "10-pattern"). For such a cycle:
  - any i < j with ε_i = 0 makes u_i an out-vortex apex (u_i → u_j, u_i → u_k,
    u_i → v), and
  - any m > k with ε_m = 1 makes u_m an in-vortex apex (u_j → u_m, u_k → u_m,
    v → u_m).
So vortex-freeness forces: for every 10-pattern (j,k): all ε before j are 1 and all
ε after k are 0.

Claim: ε is monotone non-increasing (1^a 0^b) or monotone non-decreasing (0^a 1^b).
Suppose some 10-pattern exists AND some 01-pattern (p < q, ε_p = 0, ε_q = 1) exists.
If some r > q has ε_r = 0, then (q, r) is a 10-pattern with ε_p = 0, p < q:
contradiction. Otherwise all positions after q are 1; the known 10-pattern (j,k) has
ε_k = 0, so k < q, but then m = q > k has ε_m = 1: contradiction. So if a 10-pattern
exists there is no 01-pattern, i.e. ε = 1^a 0^b; if no 10-pattern exists, ε is
non-decreasing: ε = 0^a 1^b.

Case ε = 0^a 1^b: T itself is transitive (insert v below u_a, above u_{a+1}).
Case ε = 1^a 0^b (a, b ≥ 1): switch at W = {u_1, …, u_a}. The result is transitive
with order u_{a+1} ≻ ⋯ ≻ u_{n-1} ≻ u_1 ≻ ⋯ ≻ u_a ≻ v: arcs inside the two blocks
are unchanged (transitive chains); arcs between blocks are reversed (now lower block
→ upper-block... i.e. u_m → u_i for m > a ≥ i, matching the new order); arcs u_i → v
(i ≤ a) are the reversals of v → u_i (ε_i = 1); arcs u_m → v (m > a) were already
present (ε_m = 0). So T is a switching of a transitive tournament. ∎ (c)⟹(d)

(d) ⟹ (b): induced subtournaments of transitive are transitive; for the transitive
S (all +1 above the diagonal) the Pfaffian expansion along the last vertex gives
Pf_{2k} = Σ_{j=1}^{2k-1} (-1)^{j+1} · 1 · Pf_{2k-2} = Pf_{2k-2} (alternating sum of
an odd number of equal terms), so Pf_{2k} = Pf_2 = 1 for all k. By (2), switching
preserves |Pf| of every even minor. ∎

## Computational verification (this session)

- d = 1 ⟺ locally transitive ⟺ switch-of-transitive: ALL classes agree at
  n = 3..8 exhaustively (456 classes at n=7; 6880 at n=8); labeled equality counts
  = (2n-2)!! at n = 3..7 (8, 48, 384, 3840, 46080).
- Local-order iso-class counts: 2, 2, 4, 6, 10, 16 at n = 3..8.
- Script: 04-computation/hadamard_det_census_macmini_s2.py
  Output: 05-knowledge/results/hadamard_det_census_macmini_s2.out

## Notes

- The "vortex ⟺ |Pf| = 3" 4-set dichotomy is the determinant-lens reason local
  orders are the *circular geometry* stratum: d(T) = 1 is the floor, attained exactly
  on the switching class of the principal line's origin.
- Priority: the lemma "local orders = switching closure of transitive tournaments"
  may be known (Babai–Cameron school, switching classes of tournaments); literature
  check logged in the session's research notes. The determinant characterization
  (a) ⟺ (c) appears to be new to this project.
