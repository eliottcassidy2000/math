---
id: HYP-2288
title: Reversible runners and the unit-distance trienement — both are "tournaments with ties" where symmetry controls resonance
status: OPEN (synthesis); Part-1 directional/2n-1 facts VERIFIED, Part-2 tie-degree reframe VERIFIED
source: claudebox-2026-06-03-S631
related:
  - THM-389   # the LRC trienerment (tie conventions)
  - THM-401   # pair-sum sieve modulus 2n-1
  - HYP-2130  # rigidity = orbit-type / the perspective key
  - HYP-2230  # unit-distance ↔ grid-disproof bridge (CM modulus-1 elements)
  - THM-418   # prime-shell dodge ((Z/m)* multipliers)
---

# HYP-2288 — directions, ties, and symmetry across LRC and unit distances

Two reframes (user, S631), one structure: **both problems are trienements (tournaments with ties),
and in both the tie/resonance density is controlled by a symmetry group — the orientation for LRC,
the automorphism group for unit distances — which is the perspective/rigidity key ([[HYP-2130]]).**

## Part 1 — LRC with reversible runners (opposite directions)

A single sign flip is **invisible** (`‖v t‖ = ‖−v t‖`, verified). The content lives in the
**pairwise** relative speeds: runner `i` (dir `ε_i`) vs `j` (dir `ε_j`) has relative speed
`ε_i v_i − ε_j v_j` — **same direction = a difference, opposite = a sum**. So a directional
assignment `ε ∈ {±1}` is a **tournament orientation**, and the reversible mutual-loneliness gap
`max_t min_{i<j} ‖(ε_i v_i − ε_j v_j) t‖` depends on `ε`.

**Verified (`lrc_directional_s631.out`, AP `{1,…,n-1}`):**
- BEST orientation = all-same-direction (pure differences): gap `= 1/(n−1)` (the classic AP value),
  relative speeds `{1,…,n−2}`, modulus `~n`.
- Worse orientations bring in **sums**, and the **witness denominators across orientations contain
  exactly `2n−1`** (7, 9, 11 for n=4,5,6). So **opposite directions ⟹ the `2n−1` modulus.**

> **The pair-sum sieve / shell-`2n−1` framework (THM-401/412/415/418) IS the reversible-runner LRC.**
> The `±` involution (complex conjugation on `(ℤ/(2n-1))*`, S625) is literally direction reversal;
> the orientation `ε` is the tournament; "same vs opposite direction" is "difference vs sum" is
> "modulus `n` vs `2n−1`". The "trick" the user flagged: single flips invisible, but the collective
> orientation selects sum-or-difference per pair = a tournament = sets the modulus.

## Part 2 — the unit-distance trienement

For each pair: `dist < 1` → arrow one way, `dist = 1` → **TIE** (unit edge), `dist > 1` → arrow the
other way. The unit-distance graph is the **tie subgraph**; Erdős' problem = maximize ties. Verified
(`unit_trienement_s631.out`):
- For a lattice the tie-graph is a **Cayley graph**; **interior tie-degree = #norm-1 vectors = #units
  / roots of unity**: triangular `ℤ[ω]=ℚ(√−3)` gives **6**, square `ℤ[i]=ℚ(i)` gives **4** — the
  triangular lattice wins at small `n` because it has more ties per point.
- The norm-form representation count `r(m)` for `x²+xy+y²` (Eisenstein) is the **tie-multiplicity at
  radius √m**: `r(1)=6, r(7)=r(13)=…=12, r(49)=18, r(91)=24` — large at `m` with many primes `≡1 (3)`.

> **The disproof restated:** "maximize unit distances" = "maximize tie-degree" = "find a CM lattice
> whose supply of modulus-1 elements (norm-1 vectors / bounded-height unit-modulus algebraic numbers)
> is superlinear." The tie-degree **is** the symmetry: the tie-graph's automorphism/Cayley structure.
> So **CM beats grid = a tie-graph with a larger symmetry/unit supply** — exactly [[HYP-2230]]'s
> bridge, now as a clean trienement statement.
>
> **Does it simplify the disproof?** The *statement/objective* — yes, cleanly (maximize tie-degree =
> maximize CM modulus-1 supply = maximize symmetry). The *proof's hard core* — no: that class field
> towers deliver superlinear modulus-1 supply (Golod–Shafarevich) is irreducible number theory the
> combinatorial trienement does not touch. Honest.

## The unification

Both are trienements; the tie = the **resonance** (loneliness-tight at `1/n` / unit-distance at `1`),
the arrows = the metric/velocity order, and the resonance density is set by a **symmetry group**:
- LRC: the orientation group `{±1}` (directions) → modulus `n` vs `2n−1`;
- unit distance: the automorphism group of the tie-graph (CM/Cayley) → tie-degree.

And in **both** the controlling number-theoretic object is the **modulus-1 elements of a CM field**:
literal norm-1 vectors (unit distance) and the `(ℤ/(2n-1))*` multipliers of the `2n−1` shell (LRC,
S625) — the same CM-rotation supply ([[HYP-2230]]/[[HYP-2245]]). **Trienement + symmetry-controls-
resonance is the shared skeleton; CM modulus-1 elements are the shared engine.**

## To do
1. The reversible-LRC formalization: is `min_ε max_t min ‖(ε_i v_i−ε_j v_j)t‖` (adversary orients) a
   cleaner invariant than the observer gap, and is its modulus always `2n−1`? Tie to THM-401's proof.
2. Unit-distance tie-degree as an automorphism count: make "max ties = max tie-graph |Aut|" precise;
   does it recover the `n^{1+ε}` lower bound from a symmetry/Cayley statement?
3. The shared CM-modulus-1 object: one generating/partition function (HYP-2245) whose specializations
   are the LRC `2n−1`-multipliers and the unit-distance norm-1 counts.
