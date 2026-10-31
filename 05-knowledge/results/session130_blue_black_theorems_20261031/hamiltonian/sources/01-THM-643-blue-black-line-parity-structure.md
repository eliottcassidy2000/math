---
id: THM-643
title: The BLUE/BLACK LINE PARITY STRUCTURE of the tiling fibration (strict explorer definitions) — pure-black = non-self-converse EXACTLY (all n); blue lines live on SC nodes; per-node allocation parities (SC: blue odd / black even; NS: black odd); every node emits an odd number of cross-class line-endpoints; global closed forms; the blue count per SC node is the SELF-CONVERSE HAMILTONIAN PATH COUNT H_sym — a new odd tournament invariant with mass formula Σ = 2^{(m+f)/2}
status: PROVED (T1–T5 below, elementary parity/linear-algebra arguments over the (σ,τ) frame; verified exhaustively n=3..7, all checks TRUE). C1 REFUTED at n=8 and n=9 (blue fibre 39 > 27, so even the odd-order-only cap fails); the other conjecture labels below retain their historical scope and are not promoted by this update. The unrestricted pure-blue-leaf prose is REFUTED at n=5; every even-order pure-blue node is a leaf attached to a mixed node for even n>=4 (2026-10-20 correction below).
source: mac-mini-2026-07-07-S46 (HYP-4977; owner directive: the even/odd duality of blue/black lines and node types, toward complete structural determination of the fibration)
depends_on:
  - LEM-003    # the tiling fibration is free: fiber(C)·|Aut(C)| = H(C) (canon, 2026-06-10)
related:
  - MISTAKE-033  # blue/black lines DO cross classes
  - MISTAKE-065  # grid transpose (no negation) = converse along the reversed path
  - THM-639      # the runner-world sibling frame (reversal = the same involution downstairs)
external: Rédei's theorem (H odd); odd order of tournament automorphism groups (classical).
---

# THM-643 — The blue/black line parity structure

## Current correction — 2026-10-20

The proved parity statements T1–T5 are unchanged. Two extrapolations in the
historical census discussion below are **REFUTED**, not open conjectures:

- Pure-blue does not imply leaf. The regular cyclic tournament on five
  vertices has `H=15`, `|Aut|=5`, and three blue tilings. Their complements
  reach three distinct merged nodes, so its literal and support degrees are
  both three. The [exact witness and independent audit](../../05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020/README.md#4-a-correction-pure-blue-does-not-mean-pendant)
  identify the missing coordinate: purity fixes the colour of the fibre,
  not its size.
- C1's proposed cap `3^floor((n-2)/2)` fails at `n=8`: the explorer-gauge
  masks `40077` and `40109` have blue fibre count `39>27`. Their
  `(H,|Aut|,fibre,blue)` values are `(637,1,637,39)` and `(621,3,207,39)`.
  [Direct independent permutation checks](../../05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020/audit/audit_receipt.md#independent-order-eight-cap-refutation)
  certify these witnesses without assuming the census's global maximum.
  The odd-order-only restriction also fails at `n=9`: a direct `9!`
  permutation witness has `(H,|Aut|,fibre,blue,c3)=(2659,1,2659,39,28)`,
  again `39>27`. The [independent order-nine certificate](../../05-knowledge/results/odd_blue_spectrum_hollow_cubes_20261020/audit_independent/independent_n9_cap_witness.json)
  records its explicit adjacency matrix; the new odd-order pairing theorem
  does not imply the old cap.

The strongest repaired all-order leaf statement is **PROVED**: for every
**even `n>=4`**, a pure-blue node has singleton fibre and its unique blue
line reaches a mixed node. This follows from the
[even-purity/rigidity equivalence and strong-partner proof](../../05-knowledge/results/merged_metagraph_purity_and_strong_core_20261020/PURE_BLUE_STRONG_CYCLE_REDUCTION.md#3-every-even-pure-blue-node-is-rigid).
At odd order a singleton fibre still suffices for a leaf, but purity alone
does not. This result does not require the later rigid-atom classification.

The old leaf sentence, C1 wording, and request for proofs of C1–C4 below are
retained as dated provenance, not active assertions overriding this correction.
The separate C2–C4 questions are not decided by the cap counterexample.


## The (σ, τ) frame

Tiling space = GF(2)^m, m = C(n−1,2) (tiles = non-base-path arcs over the fixed
Hamiltonian base path). Two commuting structural maps:
- **σ** = grid transpose (x,y) ↦ (n−y+1, n−x+1): a LINEAR involution permuting tile
  coordinates, with f = ⌊(n−1)/2⌋ fixed tiles (the anti-diagonal). Fix(σ) = the
  grid-symmetric ("blue") tilings, a linear subspace of dimension f + (m−f)/2 = (m+f)/2.
  Downstairs σ sends the tournament to its CONVERSE (MISTAKE-065).
- **τ** = flip = translation by the all-ones vector 𝟙 (the complement tiling, the d=m
  waggly layer): fixed-point-free; LINES = τ-orbits (2^{m−1} of them). σ𝟙 = 𝟙, so σ and
  τ commute and τ preserves Fix(σ): blue lines = τ-orbits inside Fix(σ).
- The S_n fibration sits on top: tiling ↦ iso class; fiber(C)·|Aut(C)| = H(C) (LEM-003).

## Theorems (all proved; all verified exhaustively n=3..7)

**T1 (odd fibers).** Every fiber size H(C)/|Aut(C)| is ODD. *Proof.* H is odd (Rédei);
|Aut| is odd (an automorphism of even order would contain an involution swapping some
pair u,v, sending the arc u→v to v→u — impossible in a tournament); an odd number
divided by an odd divisor is odd. ∎

**T2 (pure-black = non-SC, all n).** A class is pure black iff it is NOT self-converse;
every SC class has an ODD number ≥ 1 of blue tilings. *Proof.* A grid-symmetric tiling
represents a tournament isomorphic to its converse (σ downstairs = converse), so
non-SC ⟹ zero blue tilings. For SC classes, σ maps fiber(C) to itself (the converse
class is C); an involution on a set of odd size (T1) has a fixed-point set of the same
parity: #blue(C) ≡ |fiber(C)| ≡ 1 (mod 2). ∎
*(This upgrades the project's "verified n=3..7: transpose-self classes are never pure
black; non-transpose-self classes are always pure black" to a theorem for ALL n.)*

**T3 (blue lines live on the SC world).** Both endpoints of a blue line are blue
tilings (τ preserves Fix(σ)), hence both endpoint classes are SC. ∎

**T4 (allocation parities).** Per node: SC: (#blue tilings odd, #black tilings even);
NS: (0, odd). Cross-class line-ENDPOINT counts: every node emits an ODD number of
cross-class endpoints (≥ 1 — no node is line-isolated); SC nodes emit an odd number of
BLUE cross-endpoints (the blue graph on SC nodes has min degree ≥ 1) and an even number
of black cross-endpoints; NS nodes an odd number of black cross-endpoints. *Proof.*
τ is a fixed-point-free involution preserving both the fiber-partition's color classes;
internal pairings consume endpoints two at a time, so each cross-endpoint count is
congruent to the corresponding tiling count mod 2; apply T1/T2. ∎

**T5 (global closed forms).** #blue tilings = 2^{(m+f)/2} (so the exponent identity
(m+f)/2 = ⌊(n−1)/2⌋ + ⌊(n−2)²/4⌋... equivalently the corrected grid-sym fraction of
kind-pasteur-S20ex); #blue lines = 2^{(m+f)/2 − 1}; #black lines = 2^{m−1} −
2^{(m+f)/2 − 1}. *Proof.* Fix(σ) is a linear subspace of dimension (m+f)/2; τ acts
freely on it and on its complement. ∎

**T6 (the new invariant).** For an SC class C, the blue-tiling count is
> **H_sym(C) := #{Hamiltonian paths P of T fixed by (converse-isomorphism ∘ reversal)} / (orbit bookkeeping)**
— concretely, the number of fiber tilings fixed by σ: a "self-converse Rédei count".
It is odd (T2), and satisfies the **blue mass formula** Σ_{SC classes} H_sym-fiber-count
= 2^{(m+f)/2}, the blue companion of the classical Σ_C H/|Aut| = 2^m. Spectrum observed
(exact, n≤7): all odd values {1, 3, 5, 7, 9}, per-class max = 1, 3, 3, 9, 9 at n=3..7.

## Census (exact, n = 3..7) and conjectures

| n | classes | SC | PURE-BLUE | MIXED | PURE-BLACK | gridsym | blue lines | black lines | blue self-loops |
|---|---------|----|-----------|-------|------------|---------|------------|-------------|-----------------|
| 3 | 2   | 2  | 2 | 0  | 0   | 2   | 1   | 0     | 0 |
| 4 | 4   | 2  | 1 | 1  | 2   | 4   | 2   | 2     | 1 |
| 5 | 12  | 8  | 3 | 5  | 4   | 16  | 8   | 24    | 0 |
| 6 | 56  | 12 | 2 | 10 | 44  | 64  | 32  | 480   | 2 |
| 7 | 456 | 88 | 4 | 84 | 368 | 512 | 256 | 16128 | 0 |

Line allocation by endpoint node types (each line once; PB/MX/PK = pure-blue/mixed/
pure-black): blue lines are (MX,MX)-dominated with a few (MX,PB); black lines spread
over (MX,MX)/(MX,PK)/(PK,PK); PB nodes never touch black lines (they have no black
tilings) and appear only as blue leaves. Full tables in the `.out`.

**C1 (the 3-power cap).** max over SC classes of H_sym = 3^{⌊(n−2)/2⌋} (data: 1,3,3,9,9).
A Rédei-refinement: the self-converse path count is not just odd but 3-smooth-capped.
**C2 (blue self-loops only at even n).** Blue self-loop lines (gridsym t with flip(t)
in the SAME class) number 0,1,0,2,0 at n=3..7: they exist only at even n. (Mechanism
candidate: the score-sequence shift under 𝟙-translation obstructs class-preservation
at odd n — the shifted score multiset cannot match.)
**C3 (pure-blue = maximal symmetry).** Pure-blue classes have fiber ∈ {1,3}, i.e.
H ∈ {|Aut|, 3|Aut|}: (H,|Aut|) observed = (1,1),(3,3),(15,5),(9,9) — the transitive
class and the rotationally-symmetric (circulant-like) classes. Equivalently: a class is
pure blue iff its fiber is a single σ-orbit union of fixed points.
**C4 (odd-spectrum fullness).** Every odd value in [1, 3^{⌊(n−2)/2⌋}] is attained by
H_sym at every n ≥ 6 (at n=6,7: {1,3,5,7,9} all attained).

## What this determines (the owner's program)

Known completely: fiber sizes (H/|Aut|), the blue/black split of every fiber (H_sym and
its complement), all global masses (T5), all pairing parities (T4), and the node-type
census (pure-black = non-SC exactly; PB+MX = SC). Remaining for full structural
determination: the explicit cross-class MATCHING (which class each cross line's partner
lands in) — the line-metagraph — now constrained by T3/T4's parities and the allocation
tables; and proofs of C1–C4. The (σ,τ) frame reduces all of it to how two commuting
involutions interact with the S_n fibration — the same two involutions that appear
downstairs in the runner world as step-reversal and window-complement (THM-639).

## Files
`04-computation/gn_lines_parity_census_macmini_S46.py` (+ `.out`): the exhaustive
census, all theorem checks TRUE at n=3..7.
