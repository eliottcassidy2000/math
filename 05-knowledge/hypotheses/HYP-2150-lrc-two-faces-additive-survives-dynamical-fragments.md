---
id: HYP-2150
status: EXPLORATION + VERIFICATION — the worry-set has two block-diagonalizations; the apex is the
  rank-1 2-block of the DYNAMICAL (mod n) face; the ADDITIVE (mod 2n-1, odd) face is apex-free and
  carries the static rigidity (64-flip-lattice) that survives. Verified; the "multi-sieve = the
  additive face" claim is a conjecture.
source: claudebox-2026-06-03-S602
related:
  - HYP-2145
  - HYP-2140
  - HYP-2063
  - HYP-2075
  - HYP-2097
  - HYP-2117
---

# HYP-2150: the worry-set's two faces — additive (mod 2n-1) survives, dynamical (mod n) fragments

Seeing everything through block-diagonalization (HYP-2145): the worry-set decomposes two ways, by
two group actions, and the LRC obstruction is a block that exists in only one of them.

## The two faces

- **Dynamical face = mod `n`** (the doubling `x↦2x`, the phase dynamics; HYP-2117/S596). For even
  `n` the doubling 2-adically **collapses**; the gcd-divisor decomposition (S599) carries the rank-1
  **2-block** = the **apex** (HYP-2063). The single-corrector obstruction lives here.
- **Additive face = mod `2n−1`** (the antipodal-transversal / residue structure realizing the
  worry-set; the `64` self-converse classes, HYP-2097/oracle-S570). `2n−1` is **always odd**, so the
  doubling is a *permutation* (no collapse) and the gcd-divisor decomposition has **no even divisor
  ⇒ no 2-block ⇒ APEX-FREE**. The **static modular rigidity** lives here and **survives**.

## Verified facts (`lrc_two_faces_additive_dynamical_s602.py`)

1. **The doubling duality.** Mod `n` (even) collapses; mod `2n−1` (odd) permutes. Across `n=6..22`,
   the additive face is often maximal-mixing. **At `n=14`: dynamical mod 14 is dead, but additive
   mod `27 = 3³` is MAXIMALLY mixing** — `ord_27(2) = 18 = φ(27)`, so 2 is a primitive root mod `3³`.
   The static rigidity is *most* alive exactly where the dynamical one collapses.
2. **The additive face is apex-free.** The apex is the 2-block of the mod-`n` divisor decomposition
   (`{2,…}`); mod `2n−1` the divisors are all odd (`n=14`: `{3,9,27}`), so there is **no 2-block**.
   The apex obstruction has no counterpart in the additive face.
3. **The static rigidity skeleton.** The `64 = 2^{(n−2)/2}` self-converse classes are the
   orientations of the `(n−2)/2` antipodal pairs `{k, n−k}` (plus the apex `n/2`) — a Boolean
   flip-lattice with the AP at the tight bottom, all flips loosening. This additive combinatorial
   skeleton is rigid (the mod-`27` transversal census: `8191/8191` lonely).
4. **A 3-adic tower in the additive face.** At `n=14`, mod `2n−1 = 27 = 3³` has the divisor tower
   `{3,9,27}` — a 3-adic height-3 structure (the repo's "curvature cubed / Cassini atom"),
   paralleling the dynamical 2-adic height `v₂(n)` (H6).

## The synthesis: the apex is an even-modulus artifact

The same `2·7` seam, both faces: the prime `2` kills the dynamical face (mod `14` doubling collapse,
the apex 2-block), while the additive face (mod `2n−1` odd) has no `2` to kill it and is even maximal-
mixing. **The apex obstruction is not intrinsic to loneliness — it is the 2-block of the wrong
(even) modulus.** This is exactly why HYP-2075's pair-sum multi-sieve "has no apex": pair-sums
`v+w` push the structure into the odd additive face where there is no 2-block.

> **Conjecture (the multi-sieve is the additive face).** The pair-sum multi-sieve dissolves the
> apex because it computes in the additive face (mod `2n−1` / the pair-sum residues), which is odd
> and apex-free. The single corrector fails because it computes in the dynamical face (mod `n`),
> whose even modulus carries the 2-block. Resolving LRC = transporting the obstruction from the
> dynamical face (where it is the rank-1 2-block) to the additive face (where it vanishes), and the
> static rigidity (the 64-flip-lattice, all transversals lonely) is what guarantees no residue
> survives the transport.

## Open / next

- Verify the "multi-sieve = additive face" claim directly: is the pair-sum sieve's modulus the
  additive face's `2n−1` (or its odd part)? Map HYP-2075's pair moduli to the mod-`2n−1` residues.
- The additive face's own height (the 3-adic `{3,9,27}` tower at n=14) — when `2n−1` is a prime
  power `p^a`, is the static rigidity graded by a `p`-adic height, dual to the dynamical `v₂(n)`?
- The cleanest cases are even `n` with `2n−1` a max-mixing prime power (n=6 mod 11, n=10 mod 19,
  n=14 mod 27): predict these have the cleanest transversal censuses.

**Artifacts:** `04-computation/lrc_two_faces_additive_dynamical_s602.py` (+`.out`),
`07-reflections/lrc-two-faces-block-diagonalization-s602.md`. Builds on HYP-2145 (divisor blocks),
HYP-2140 (rigidity=orbit-type), HYP-2063 (apex/2-block), HYP-2075 (pair sieve), HYP-2097 (64
classes mod 2n-1), HYP-2117 (doubling).
