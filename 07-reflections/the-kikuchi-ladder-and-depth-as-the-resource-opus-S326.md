# The Kikuchi ladder and depth-as-the-resource — arXiv:2607.14068 woven into the repo

*opus-2026-07-16-S326; owner: deeply consider the even-uniform hypergraph
Moore bound paper (Bandeira–Kunisky–Nizić-Nikolac–Pesenti–Wang) against the
repo's work; extend usefully.*

## What the paper does

The hypergraph Moore bound problem: how many hyperedges force a short EVEN
COVER (a set of hyperedges covering every vertex an even number of times —
the hypergraph girth analogue, Feige's conjecture lineage). The technique:
lift the k-uniform hypergraph to a KIKUCHI MATRIX indexed by ℓ-subsets of
vertices; spectral norm bounds at level ℓ certify the existence of even
covers that no lower level can see; the even-uniformity assumption makes the
lift's combinatorics clean. The tension the bound quantifies: density vs the
LEVEL needed to certify structure.

## The weave (five threads already in the repo)

1. **Even covers ARE the repo's F₂ world.** A subset of hyperedges with
   every vertex evenly covered is exactly a cycle-space element — the
   repo's EVEN GRAPHS (E_n, the first-class dual of G_n: "tiling → XOR of
   fundamental cycles → even graph"). The paper studies, in hypergraph
   generality, the object the repo has as the metagraph's dual; the locker
   theme (S316: ζ mod 2, Rédei's ⊕P taming, triangular-parity stripes) is
   the same parity functor.
2. **Kikuchi lifts = the repo's body/level constructions.** Indexing by
   ℓ-subsets to gain spectral power is precisely the C(14,4) = 1001
   body-tree move (kps THM-738: the j = 3 level certified the ≤3-far
   branch; j = 4 = 2002 bodies in progress) — the LRC deck IS a Kikuchi
   ladder over the runner hypergraph, discovered independently.
3. **Depth-as-the-resource is one phenomenon.** The repo now has four
   instances: the α-depth cracking the fingerprint twins (THM-890: level-2
   spectra fail, depth ≈ #components succeeds); the OCF 2-adic tower (each
   digit of H needs one more cycle-collection layer, S316); the Vandermonde
   truncation ladder (S317); and the Hunter crossing (second-order
   Bonferroni seeing what first-order provably cannot, THM-856/863). The
   paper's Moore bound is the QUANTITATIVE version: it prices the level ℓ
   needed as a function of density — the missing meta-theorem shape for the
   repo's ladder phenomena. CONJECTURE (the useful extension): the
   two-sheet selector residual has a Kikuchi-type price — the level of
   subset-lift needed to separate survivor packets grows with the return
   complexity N_R (THM-817's Θ(B) satellites = the high-price rows), and
   THM-891/893's coincidence says the price at level "folded-shadow" is
   INFINITE — separation needs the frame-carrying lift, exactly as even
   covers need the right ℓ.
4. **Even-uniformity = the repo's even-n half-coset world.** The paper's
   structural assumption (all hyperedges even) is the same parity door
   through which E₈ entered at n = 8 (all d_v odd ⟺ the half-coset;
   THM-868): evenness assumptions buy exact lattice/spectral structure on
   both sides.
5. **The spectral certification style.** Kikuchi norms certify combinatorial
   existence; THM-863/864 certify non-coverage from pairwise overlap norms;
   THM-894's spectral excess ranks resonance feedback. Same proof-shape:
   moments of a lifted operator standing in for an exhaustive search.

## The extension program (concrete, fleet-ready)

- **The LRC Kikuchi ladder, formalized:** level-1 = single-comb densities
  (the dead union bound); level-2 = the resonance matrix M_ρ (THM-894) and
  the Hunter tree (THM-856); level-3 = the triple-overlap tensor
  μ(D_i ∩ D_j ∩ D_k ∩ E) — the next rung, whose top singular value should
  cross the m′ = 8 wall that the tree bound cannot (22m′ ≤ 165 dies at 8;
  the triple analogue's coefficient arithmetic is the first computation).
- **The even-cover face of the deck:** the c = 7 token criterion (X⁷ − X |
  Π(X − k_a), THM-773) asks for a FULL cover of F₇; the even-cover relaxation
  (cover every residue evenly) is the deck's cycle space — is there a
  THM-779-style walk criterion for even coverage, and does the Moore-bound
  machinery price it?
- **The n = 24 Niemeier question (S319) meets the paper:** the Leech's
  root-freeness is a girth-type statement (no norm-2 vectors = no short
  relations); the score-slice's Niemeier selection is a Moore-bound-flavored
  question about which lattice the tournament hypergraph's even covers can
  see.

## The honest boundary

The paper's regime (random-like dense hypergraphs, asymptotic ℓ) is far from
the repo's exact small-parameter world; what transfers is the LADDER
ARCHITECTURE and the price-of-level question, not the theorems. The
extension program above is stated so each rung is an exact finite
computation in the repo's style.

Cross-refs: THM-738/856/863/864/868/890/891(893)/894; S316 (parity tamer),
S317 (Vandermonde), S319 (Niemeier); the E_n canon.
