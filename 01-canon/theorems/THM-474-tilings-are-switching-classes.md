# THM-474 — The Gauge Theorem: tilings ARE switching classes of tournaments

**Status:** PROVED (mac-mini-2026-06-10-S2; 5-line proof below). Adversarial check:
the counting identity and the statement were independently derived twice this session
(two-graphs literature thread + this write-up); verification of the proof is
elementary.
**Provenance:** mac-mini-2026-06-10-S2 (session: Hadamard × tournaments × odd
functions × simplicial geometry). Related: THM-447 (normalization dictionary:
"normalizing a skew-Hadamard matrix = fixing the Hamiltonian path"), THM-468/472/473
(determinant lens), Babai–Cameron EJC 7 (2000) #R38 (switching classes ↔ oriented
two-graphs), T777.

## Statement

Fix the base Hamiltonian path P₀: n → n-1 → ⋯ → 2 → 1. Every switching class of
tournaments on [n] contains EXACTLY ONE tournament containing P₀.

Consequently the staircase tiling space Q_m (m = C(n-1,2)) — the project's central
object — IS a complete and irredundant system of representatives for the 2^(C(n-1,2))
labeled switching classes (= oriented two-graphs = skew two-graphs) on [n]:

    tilings  =  switching classes  =  oriented two-graphs.

Every switching-class invariant (the determinant d(T) = det(I+S)/2^(n-1) of
THM-468/472, the Seidel spectrum, Moorhouse fingerprints, the cyclic-triple
indicator g(x,y,z)) is therefore a well-defined function ON TILINGS, and the
tiling model is the canonical switching gauge. THM-447's matrix-level statement
("normalizing a skew-Hadamard matrix = fixing the Hamiltonian path in the tiling
model") is the bordered version of this theorem.

## Proof

*Uniqueness.* Suppose T and σ_W(T) (∅ ⊊ W ⊊ [n]) both contain P₀. P₀ visits every
vertex, so some consecutive pair (k, k-1) of P₀ crosses the cut (W, W^c); switching
reverses that arc, so σ_W(T) does not contain P₀ unless the cut is trivial —
contradiction.

*Existence.* Given any tournament T, choose the diagonal gauge greedily along the
path: set d_n = +1, and for k = n, …, 2 choose d_{k-1} ∈ {±1} so that
d_k · S_{k,k-1} · d_{k-1} = +1. Then σ_D(T) contains every base-path arc k → k-1,
i.e. contains P₀. ∎

(Counting check: 2^(C(n,2)) tournaments / 2^(n-1) per switching class
= 2^(C(n-1,2)) classes = number of tilings. Verified.)

## Notes

- This puts three quotient layers above the tiling hypercube:
  tilings (labeled switching classes) → G_n (iso classes) → G_n/Z_2 (merged) →
  S_n := switching classes up to iso (Babai–Cameron count over LEVEL permutations,
  OEIS A049313: 1, 1, 2, 2, 6, 12, 79 for n = 2..8). S_n is a NEW, uncomputed
  metagraph layer for the project — backlog lead.
- Even/odd duality: labeled switching classes of tournaments are equinumerous with
  labeled even graphs (2^(C(n-1,2)), both = cycle-space-sized); at the iso level the
  counts diverge (A049313 vs A002854 = V(E_n)). Mallows–Sloane pairs two-graphs
  with Euler graphs; the "odd Mallows–Sloane partner" of A049313 is open.
