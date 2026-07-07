---
source: opus-2026-07-07-S136
status: frame + dictionary (owner directive) with proved lemma-level facts + a new exact engine;
  the (A') evidence it produces is exact-exhaustive on boxes (see companion .out files)
tags:
  - lonely-runner
  - LRC14
  - tiling-model
  - order-cells
  - movie
  - base-path
  - affine-invariance
  - density-floor
---

# The runner movie is a tiling: the sorted-path frame and the order-cell engine

**opus-2026-07-07-S136.** Owner: *a tree on 8 events has 7 edges and that is just the
Hamiltonian path connecting each element in an 8 player tournament; the tiling model leverages
viewing a tournament from one of its Hamiltonian paths to reveal the symmetry in the isomorphism
class structure forced by the nature of intersubjective binary relation itself — apply similar
analysis to lonely runners.* This note is that analysis, carried out to the point where it
produces proofs and a new exact computational engine.

> **CONVERGENCE NOTE (same hours):** the owner fanned this directive fleet-wide; kps-S62
> ("step-gauge": sorted speeds incl 0 = base path, steps = edge weights, reversal = complement,
> palindromic-extremizer test) and mac-mini-S43 (Hamiltonian-path/step-sequence frame:
> palindrome/wall-count/lattice-class theorems) reserved the same frame concurrently with this
> note's checkpoint. Three independent instantiations = the hint was precise. Division: I cede
> theorem-naming on the frame facts to mac-mini-S43's write-ups where they overlap my L1/L2
> (wall-count, lattice-class — convergent derivations, cite both); kps-S62 owns the step-gauge
> strata closure; **this session's differentiated deliverable is the EXACT ORDER-CELL ENGINE +
> the exhaustive exact (A′) verification with certified tight loci and isolation gaps** (§3),
> which neither reserve mentions.

## 1. The dictionary

A tournament is a complete binary relation on n players; the tiling model fixes a Hamiltonian
base path (a spanning tree, n−1 arcs on n events) and codes the remaining C(n−1,2) pairs as free
bits — cut space (the path) ⊕ cycle space (the tiles), with iso classes appearing as orbits of a
hypercube quotient.

The lonely-runner instance carries the same structure, with **time as the walk parameter**:

| tournament side | runner side (14 events: observer 0 + speeds v₁<…<v₁₃) |
|---|---|
| base Hamiltonian path P₀ (n−1 arcs) | the **sorted-speed path**: consecutive gaps `d_j = v_j − v_{j−1}` (13 edges; `d₁ = v₁` is the **observer's edge**) |
| tiles = non-path pairs, free bits | pairwise differences `v_j − v_i` = **interval sums of d**; their *temporal flips* are the free structure |
| a tournament = a tiling ∈ Q_m | the **state at time x** = the circular order of the config (equivalently the half-turn tournament of THM-373) |
| flipping one tile (wiggly move) | crossing one **wall** `x = m/(v_j − v_i)`: one adjacent transposition |
| the tiling hypercube | the pair-flip hypercube; the x-circle traverses a **closed walk** in it (THM-373's "closed finite tournament walk" — this canon fact IS the runner tiling model) |
| transitive tournament = all-zeros tiling | **the AP = all-gaps-one** `(d = 1,…,1)` — the pure base path |
| S_n relabeling; iso classes | the **affine group**: dilation scales `d`; translation shifts **only the observer's edge `d₁`** — affine normalization = gauge-fixing the frame |
| complement symmetry, SC classes, G_n/ℤ₂ | **reversal** `v ↦ (v_max+v_min) − v` = reading the gap vector backwards; μ and E[maxgap] are reversal-invariant, so they live on the **palindrome-quotient**; palindromic-gap families (AP!) are the SC analog |
| cut ⊕ cycle (GF(2)) | in the path basis `K_j = Σ_{i≥j} m_i` (unitriangular), the Fourier **deviation lattice** of every gap functional is `{K : K₁ = 0, Σ_j d_j K_j = 0}` — **gauge-fixed at the observer's edge, orthogonal to the gap vector**. Translation symmetry = K₁ = 0; the frame `d` spans the "cut" direction; deviations are the "cycles" |

Three of this week's hard-won lessons become one sentence in the frame: *the dilation artifact
(kps-S56), the co-offset rotation (boxeph), and the anchored-tail affine failures (opus-S134)
are all "you compared tilings without fixing the base path first."*

## 2. Lemma-level facts the frame yields immediately

- **L1 (movie complexity; proved).** The number of order cells per unit time, counted with
  multiplicity, is `Σ_{i<j} (v_j − v_i)` — and consecutive integers minimize every pairwise
  difference simultaneously, so **the AP is the simplest movie**: 364 cells at k=13
  (`Σ_{d=1}^{12} d(13−d)`), versus 506 (GW), 564 (parity record), 1208 (primes), 1332 (the CE
  adversary). The cell-count functional is even (the closed walk must return: the order at
  `x→1⁻` is the reversal of the order at `x→0⁺`, forcing `Σ(δ_ij − 1) ≡ C(k,2) (mod 2)`).
- **L2 (the deviation lattice in the path basis; proved change of variables).** kps-S59's
  deficit frame ("all deviation flows through zero-sum weight-≥3 relations") is, in gap
  coordinates, the statement that deviations live in `d⊥ ∩ {K₁ = 0}`. The AP (`d ∝ 1`) makes
  `d⊥` the full zero-sum lattice — the maximal deviation lattice (= HYP-4532's Cohn–Elkies
  extremality, now visible as "the base path with no excess"), which is *why* the AP is
  simultaneously max-energy, max-relation, min-distinct-differences, min-cells.
- **L3 (what (A′) says in the frame).** The load-bearing lemma "the AP minimizes the tail
  μ_{1/7}" reads: *the pure base path is the tail-minimal movie* — the exact analog of
  transitive-extremality statements on the tournament side (H(transitive) = 1 minimal). The
  M-ladder rigidity ("detuning jumps to the next rung") is the analog of the metagraph's
  H-gradient. One frame, both monotonicities.

## 3. The engine (new capability): exact μ_θ for arbitrary E

Within an order cell every gap is affine in x; μ_θ(E) is the sum over cells of the superlevel
measure of a max of k affines — **exact rational arithmetic for any integer set**. Implemented
(`lrc_exact_mu_ordercells_opus_S136.py`, adapting death-star-S1's corrected mean integrator —
the mean and the tail are two functionals of one cell decomposition):

- Validation: reproduces all six roof values `691/735 … 477/1078` exactly.
- **First exact general-E μ values in the repo**: `μ_{1/7}(GW) = 4609589/7834365`,
  `μ_{1/7}(parity record) = 17159/32340`, `μ_{1/7}(2AP∪{13}) = 48481/97020`, primes
  `= 4624167679/4888643760`, … (every earlier general-E μ was a ~4-digit grid estimate).
- **Exhaustive exact (A′) with certified tight loci**: on the k=8 height-14 box (1716
  normalized classes): 0 below the bar and **exactly one equality class — the AP** — the first
  *certified* uniqueness of the tight locus (a grid can never distinguish equality from
  4-digit proximity). The k=8..13 boxed sweep with exact **isolation gaps γ_k** (runner-up minus
  bar) is running as `lrc_exact_Aprime_exhaustive_opus_S136`.

The isolation gaps matter beyond evidence: an exact γ_k > 0 on a box turns "(A′) on the box"
into a *finite certificate* (the box is checked; outside the box the diameter/intersected
ledgers of kps-S59/S60 + monad-S2 take over at the DAG bars). The frame thus splits (A′) the
same way the tiling model splits tournament counting: a finite marked-frame core plus a
symmetry-controlled tail.

## 4. Honest status and what's next

- The dictionary is a frame, not a theorem; L1/L2 are proved but elementary; the engine and its
  exhaustive outputs are exact.
- Not yet done: the walk-topology angle (the closed walk's homotopy/winding data as an affine
  invariant — the metagraph analog: which order-state sequences are realizable, and does μ
  factor through walk statistics?); the palindrome/SC-quotient census (are near-palindromic gap
  vectors the μ-extremal families the way SC classes organize G_n?); using L2 to prune the
  lattice sums in monad's CE targets (the S1 bound is a short-vector count in `d⊥ ∩ {K₁=0}`).
- Owner priority honored: no Lean this session; the AP₇₆ ledger certificate remains the named
  formalization task for whoever picks it up after the math settles.
