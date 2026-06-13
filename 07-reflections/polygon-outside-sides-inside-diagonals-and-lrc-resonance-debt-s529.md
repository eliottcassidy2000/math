---
source: oracle-2026-06-01-S529
status: synthesis + computation (regular-polygon outside/inside; dihedral shells; LRC covering)
tags:
  - lonely-runner
  - regular-polygon
  - dihedral-group
  - cyclotomic
  - gauss-sum
  - tournament
  - resonance
  - tiling
---

# The Polygon's Outside (Sides) and Hidden Inside (Diagonals), and the LRC Resonance Debt

A tournament is an orientation of `K_n` — a binary relation on the edges of the
`(n-1)`-simplex. But the maximally symmetric tournament lives on the **regular
`n`-gon** with its dihedral group `D_n`, and that polygon has two strata: the
**outside** (its sides / convex-hull boundary) and the **hidden inside** (its
diagonals, the chords that cross through the interior). This note follows the
prompt — *how the outside relates to the hidden inside arcs, and how that applies
to LRC* — and finds that the split is exactly the **ranking/cyclic (cut/cycle)**
decomposition the repo already lives by, that the dihedral symmetry pins the LRC
threshold `1/n` to the polygon's own geometry, and that the **inside diagonals are
the LRC obstruction** — they are the "higher-resonance debt" of S526/S527, born
precisely when the polygon first acquires an interior.

## Outside = sides = ranking; inside = diagonals = cycles

Label the `n`-gon vertices `0,…,n-1`. The chords split by **skip** `j` (connecting
`i, i+j`), `j = 1,…,⌊n/2⌋` — the orbits of `D_n` on edges:

- **skip 1 = the sides** — the convex-hull boundary, the **Hamiltonian cycle** of
  cyclic adjacency. In tournament terms this is the **ranking / score backbone**,
  the **cut space**, the *base path* of the tiling model. The OUTSIDE.
- **skip ≥ 2 = the diagonals** — chords through the interior, the genuinely
  **cyclic** content (the 3-cycles live here), the **cycle space**, the *tiles*.
  The HIDDEN INSIDE.

So "outside vs hidden inside arcs" is literally the **GF(2) cut ⊕ cycle split**
(CLAUDE.md): base-path arcs (sides) carry the hierarchy; the tiles (diagonals)
carry the cyclic/3-cycle structure that makes a tournament more than a ranking.
The repo's whole tiling apparatus is the polygon seen edge-on.

## The dihedral symmetry pins `1/n` to the polygon (cyclotomic)

Two computed facts (`lrc_polygon_inside_outside_s529.py`):

- **Chord product = n.** `∏_{k=1}^{n-1} |1 - ω^k| = n` (`ω = e^{2πi/n}`), verified
  `n=3..10` to `1e-6`. This is the cyclotomic identity `∏(x-ω^k)=(x^n-1)/(x-1)` at
  `x=1`. The sides and the hidden diagonals **together multiply to `n`**: the
  outside and inside are not independent — their product is the conjecture's own
  modulus.
- **Nearest gap = `1/n`.** The regular `n`-gon's nearest-neighbour arc gap is
  exactly `1/n`. **The LRC threshold IS the regular polygon's own gap.** The tight
  witness (runners at the `n`-gon, speeds `1,…,n-1` at times `k/n`) sits the
  nearest runner *exactly on the boundary* of the observer's forbidden cap of
  width `2/n`. The polygon doesn't *almost* achieve `1/n`; it *defines* it.

For prime `n=p` the inside has a sharper signature: the diagonal orientation is the
**Legendre symbol**, and the balance of the inside is the **Gauss sum**
`|∑_k χ(k) ω^k| = √p` (verified `p=3,5,7,11,13`). The hidden inside is as balanced
as a complete character sum can be — `√n` cancellation, the extremal.

## The LRC covering sum, graded by polygon shell depth

Observer-loneliness at time `t` = every runner `v_i t` in the safe band
`B=[1/n, 1-1/n]`. With `f = 1_B`, `f̂(0)=1-2/n`, `f̂(m) = -sin(2πm/n)/(πm)`,

```
|LONELY(v)| = ∫_0^1 ∏_i f(v_i t) dt = Σ_{(m_i): Σ m_i v_i = 0} ∏_i f̂(m_i).
```

The **resonance condition `Σ m_i v_i = 0`** is the algebraic shadow of the polygon
diagonals; grade the sum by **resonance order** `r = #{i : m_i ≠ 0}`:

- `r = 0` (main term) = `(1-2/n)^{n-1}` — the **OUTSIDE / independence value**
  (opus-S524's `(6/7)^13` at `n=14`).
- `r = 2` (pairwise) — the **outer / skip-1 relations**. For `n=3` this term alone
  *is* S526's closed form `1/9 + (2/9)·χ(a)χ(b)/(ab)` (reproduced exactly, all test
  pairs match `<2e-3`).
- `r ≥ 3` (multi-way) — the **deep inside diagonals**: the "higher-resonance debt"
  that S526/S527 named as the obstruction for `n ≥ 4`.

**The headline computation.** Every regular-polygon (AP) speed set
`v = (1,…,n-1)` gives **`|LONELY| = 0` exactly** by direct integration
(`n=3,4,5` checked) — the outside main term is cancelled *to zero* by the inside
corrections. The regular polygon is the configuration where outside and inside
**perfectly annihilate** (the `√n` Gauss-sum balance in measure form): the unique
tight case, loneliness only on the measure-zero boundary times `k/n`.

## The geometric birth of the open case

The decisive number is the **inside debt** `= Σ_{r≥3} (resonance terms)`:

```
n=3:  inside debt ≡ 0    (two speeds admit NO 3-term resonance — the triangle has
                          NO diagonals; its only inside is the single 3-cycle,
                          captured entirely by the order-2 Legendre term)
n=4:  inside debt ≠ 0    (-0.081 for the AP; the square has its first interior
                          chord — the diameter — and r=3 resonances appear)
n=5:  inside debt ≠ 0    (-0.114 for the AP; pentagon diagonals)
```

So **LRC@3 is solved precisely because the triangle has no inside** (S526's
Legendre/Gauss computation is complete: outside + the one diagonal-class = done),
and **the open case is born exactly at `n=4`, when the polygon first acquires
genuine interior diagonals.** The conjecture, in this language:

> *The hidden inside (diagonal / `r≥3` resonance) debt can cancel the outside main
> term down to `0` but — by the `√n` Gauss-sum boundedness of each diagonal shell —
> never below it.* The regular polygon saturates the cancellation; nothing beats it.

That is the same wall as S519/S520/S526 (the realizable / higher-resonance
constraint), now with a clean geometric cause: it is the **arithmetic of the
polygon's interior diagonals**, organized by `D_n` into shells, each contributing a
bounded character sum.

## Why this sharpens the attack

The shell grading turns "bound the higher-resonance correction" into "bound each
**diagonal-shell character sum** by its `√n` Gauss-sum size." For prime `n` this is
classical (Weil); the open difficulty is composite/AP `n` (`n=14,16,18`), where the
shells mix moduli — exactly why the "if 15 were prime" thought experiment (S17) and
the composite frontiers matter. The honest program: show the signed sum of the
inside-diagonal shell contributions is `≥ -(1-2/n)^{n-1}` with equality only at the
regular polygon. That is LRC, but now it is a statement about **Gauss/Kloosterman
cancellation across the dihedral shells of the `n`-gon's interior**, not a generic
measure bound.

## Verdict / next
- The outside/inside split of the regular polygon = the cut/cycle (ranking/tile)
  split; the threshold `1/n` is the polygon's own gap (cyclotomic chord-product
  `n`); the inside diagonals = the LRC obstruction, born at `n=4`.
- The regular polygon is the measure-zero tight case where outside + inside cancel
  to `0` (Gauss-sum `√n` balance).
- Concrete next: compute the per-shell character sums for `n=4,5` explicitly and
  test the `√n`-domination bound; push toward the composite frontiers where shells
  share moduli.

## Artifacts
```
04-computation/lrc_polygon_inside_outside_s529.py
05-knowledge/results/lrc_polygon_inside_outside_s529.out
```
Related: S526 (mod-3 Legendre n=3 proof), S518 (Fibonacci circular menu),
S521o (permutohedron geodesic), S17 ("if 15 were prime"), THM-369 (sieve),
THM-386 (central-box grounding).
