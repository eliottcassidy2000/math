---
source: opus-2026-07-11-S235
status: INTEGRATION of paths 1 (quantify the detuning) + 3 (inverse-theorem bridge), with a fleet
  course-correction. NEW: the band-edge margin lemma (elementary, verified 0/19999) — clearing at a modulus
  not divisible by 14 forces M > 1/14 — makes path 1's strict margin a FREE COROLLARY of bounded-clearing.
  Path 3 ("covering ⟹ near-AP") is REFUTED as backwards (the hard core is loose; energy floor runs the other
  way); the repo's inverse invariant is E₃/Schur, not E₂/BSG. Both paths reduce to the S230/S231
  anti-concentration.
tags:
  - lrc14
  - band-edge-lemma
  - divisor-complete
  - the-margin
  - additive-energy
  - course-correction
  - E3-schur
---

# The band-edge margin lemma — and why both paths reduce to bounded-clearing

**opus-2026-07-11-S235.** Owner: work paths 1 and 3 simultaneously, pull from the fleet, integrate. Doing so
produces one clean new lemma, collapses path 1 into a corollary, and refutes path 3's framing.

## The band-edge margin lemma (new, elementary, verified 0/19999)

> **Lemma.** If a family `S` has a lonely time at a modulus `q` with `14 ∤ q` (i.e. `bandCount(S,q,p)=0` for
> some `p`), then `M(S) ≥ ⌈q/14⌉/q > 1/14`.

*Proof.* Clearing at `q` means every runner lands in the safe band `[⌈q/14⌉, ⌊13q/14⌋]`, so
`‖v_i p/q‖ ≥ ⌈q/14⌉/q`. When `14 ∤ q`, `⌈q/14⌉ > q/14`, hence `⌈q/14⌉/q > 1/14` strictly. (At `14 | q` the
edge is exactly `1/14`.) ∎

**Corollary (tight-locus characterization).** `M(S) = 1/14` (tight) ⟹ `S` clears **only at multiples of 14**.
Verified: the AP `{1..13}` and the sporadic `V* = {1..11,13,24}` clear at `q ∈ {14,28,42,56}` and at **no**
`14∤q`. (This is the multiplier-level face of THM-610: a tight/covering family hides at a denominator
divisible by 14.)

## Path 1 collapses into a corollary

S234 reduced LRC(14) to **divisor-complete ⟹ M > 1/14** (via THM-366: covering ⟹ divisor-complete).
Divisor-complete families have a multiple of 14, so they do **not** clear at `q=14`; but they **do** clear at a
non-multiple-of-14 modulus `q ∈ [15,41]` (verified; adversarial worst `q=31`). By the band-edge lemma this
gives `M ≥ 2/27 > 1/14` — so **every divisor-complete family is loose**. Therefore:

> **Path 1's strict margin `M > 1/14` is a *free corollary* of "divisor-complete clears at a bounded
> non-14 modulus."** It is not a separate difficulty. Proving the bounded-clearing (the S230/S231
> anti-concentration, verified `q ≤ 60` diameter-free) yields *both* the loneliness certificate *and* the
> strict margin, simultaneously.

There is no independent "detuning bound" to prove — the detuning *is* the clearing at a non-14 modulus, and
the margin follows from the arithmetic of the band edge.

## Path 3 is backwards — the fleet course-correction

The intended path 3 was "covering ⟹ near-maximal additive energy ⟹ (BSG/Freiman) near-AP." Pulling the
fleet's state (THM-656/660, opus-S181/S182/S225, kps cont.36) refutes the framing:

- **The energy floor runs the opposite way.** The AP is the **maximum**-energy configuration and the
  **minimizer** of every floor (`μ ≥ λ²/(R2·V1+λ²)`, `μ ≥ E[W]²/E[W²]`): *high* `R2` ⟹ *weak* floor. The
  theorems prove "even the AP (worst case) clears." "Covering ⟹ near-AP" is a heuristic ordering axis, **not
  a theorem**, and points the wrong direction.
- **kps's decoupling (cont.36):** the window-hard covering cores are **loose, not near-tight** — confirmed
  here (all divisor-complete families are loose, `M > 1/14`). So an energy-floor + "covering ⟹ near-AP"
  strategy pushes against the fleet's own finding that the covering-hard and energy-tight loci are *different
  families*.
- **The right inverse invariant is E₃, not E₂.** opus-S182 (with LEM-015 proved, Lean kernel-pure) shows the
  loneliness deficit tracks the **Schur-triple count `E₃`** (dilation-invariant, matching loneliness's
  symmetry), not the additive energy `E₂`/doubling (translation-invariant — "cannot distinguish the tight AP
  from its loose translate `{2..14}`", which is exactly a divisor-complete loose family). BSG/Freiman-3k−4
  remain parked; if the inverse-structure bridge is pursued, target `E₃` via the S182/LEM-015 machinery.

## Net

Both paths, correctly integrated, reduce to the **same** open statement — *divisor-complete families clear at
a bounded non-14 modulus* (the S230/S231 anti-concentration) — now equipped with:

1. a **free strict margin** `M > 1/14` (the band-edge lemma), so bounded-clearing gives the whole thing; and
2. a clean **tight-locus characterization** (tight ⟺ clears only at multiples of 14), meeting THM-610.

And path 3 is **refuted as framed**: the hard core is loose, the energy floor is minimized at the AP, and the
honest inverse invariant is `E₃` (Schur), not `E₂` (BSG). No new open object was created; the residual is the
anti-concentration, verified but unproved. The genuinely new, proved piece is the band-edge lemma and its
characterization of the tight locus by the divisibility of the clearing denominator.

→ THM-366 (covering ⟹ divisor-complete), THM-610 (covering hides at `q* ≡ 0 mod 14`), opus-S234 (the
divisor-complete reduction), opus-S230/S231 (bounded-clearing anti-concentration — the residual), THM-656/660
(energy floor — AP is the minimizer, the path-3 correction), opus-S182 + LEM-015 (E₃/Schur — the right
inverse invariant), kps cont.36 (decoupling: window-hard = loose). Files:
`lrc14_band_edge_margin_opus_S235.py` (+`.out`).
