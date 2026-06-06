---
source: monad-explorer-2026-06-06-S702
status: PROVED/classical synthesis + VERIFIED (exact). Sharpens & partly CORRECTS the framing of
        the S699/HYP-2257 probe ("a 2D group between the triangular lattice kappa=6 and the CM
        field whose norm-1 layer beats 3n"). Does NOT prove or refute the external CM exponent
        claim (n^{1.014}); it places that claim structurally.
tags: [unit-distance, kissing-number, norm-form, popular-norm, divisor-function, imaginary-
       quadratic-order, roots-of-unity, CM-field, eisenstein, gaussian, disproof, 3n-cap,
       rank, s699, trienerment, additive-face]
---

# The norm-layer dichotomy: "between triangular and CM" is the wrong axis

## The probe (S699 / HYP-2257, dispatched seed)

> *Is there a 2D-realizable group strictly between the triangular lattice (kissing `κ=6`) and
> the CM field, whose norm-1 layer beats `3n` at moderate `n`?*

S699 framed the unit-distance disproof as: 2D lattices are **kissing-capped** (`κ≤6 ⟹ ≤3n` unit
distances), and the CM field **escapes** by making the `=1` ("norm-1") layer a *field-norm*
condition (unbounded). The probe asks for a continuum of 2D groups bridging `κ=6` and the CM field.

**Finding (S702): the probe conflates two different "norm-1 layers," and once they are separated
the bridging spectrum collapses.** The kissing axis and the exponent axis are both *flat* across
all 2D lattices; the CM escape is a **rank/dimension jump**, not a 2D-group spectrum.

## Two norm layers — keep them apart

For a 2D lattice `L` with integral norm form `Q(x,y)` (so squared distances are values of `Q`):

| layer | what `=1` means | controlled by | size | nature |
|---|---|---|---|---|
| **(K) kissing** | the **minimal** nonzero norm | `κ` = #minimal vectors | `κ ≤ 6` | **geometric** |
| **(P) popular** | a **fixed** norm `D` with many reps | `r_Q(D)` = #`{v:Q(v)=D}` | `sup_D r_Q(D)=∞` | **arithmetic** |

The "`3n` cap" is a statement about **layer (K) only**: the minimal-distance graph on `n` lattice
points has `≤ κn/2 ≤ 3n` edges. Beating `3n` means finding *any* distance realized `>6` times per
point — i.e. a norm `D` with `r_Q(D) > 6`. That is **layer (P)**, and it has nothing to do with `κ`.

## The rigorous results (all exact / classical; verified `…s702.py`)

**1. Kissing is universally capped, and it equals the number of roots of unity.**
For every 2D lattice, `κ ≤ 6`, with `=6` *iff* hexagonal (verified across a rhombic/rectangular/
oblique spectrum: generic `κ=2`, square `κ=4`, hex `κ=6`, nothing above 6). The clean reason:

> **`κ(L) = #` roots of unity of the associated imaginary quadratic order `O`.**
> Minimal vectors = norm-1 elements = units = roots of unity. In degree 2 there are at most 6
> (only `ℚ(ω)`: 6; only `ℚ(i)`: 4; all others: 2). So `κ ∈ {2,4,6}`, capped at 6.

Verified: `κ(ℤ[i]) = 4 = #`units `ℤ[i]`, `κ(ℤ[ω]) = 6 = #`units `ℤ[ω]`. **The kissing number IS the
root-of-unity count of the point group.** This pins S699's "`=1` layer = roots of unity" exactly:
in *degree 2* that layer is the unit group, capped at 6.

**2. The popular layer beats `3n` already inside `ℤ²` — the LOW-kissing lattice.**
First norm with `r_Q(D) > 6`:
- **Square `ℤ²` (`κ=4`):** `D = 5`, `r = 8` (the `(±1,±2),(±2,±1)`). `r(D)=4·(d₁−d₃)(D)`.
- **Triangular `ℤ[ω]` (`κ=6`):** `D = 7`, `r = 12` (7 splits in `ℤ[ω]`). `r(D)=6·(d₁−d₂)(D)`.

On a finite window the bulk density (`r/2`)·`n` beats `3n` once the Harborth perimeter loss
`~c√n` is dominated — a **moderate-`n`** crossover:

```
   Z^2  (kappa=4): popular distance sqrt(5), bulk 4n, crosses 3n at n=144 (12x12 window)
   Z[w] (kappa=6): popular distance sqrt(7), bulk 6n, crosses 3n at n= 64 ( 8x8 window)
```

> The square lattice has the **smaller** kissing number (`4<6`) yet **still beats `3n`**. So the
> escape from `3n` is **orthogonal to the kissing number**. "Between triangular and CM" is the
> wrong axis: beating `3n` is a **layer change** (minimal → popular norm), not a **group change**.

**3. The popular-norm growth rate is the SAME for every 2D lattice — exponent `1+o(1)`.**
`r_Q(D) = w·∑_{δ|D} χ_Δ(δ) ≤ w·d(D)`, so `max_{D≤X} r_Q(D)` is governed by the **divisor
function**, `= X^{(log2+o(1))/loglog X}` (Wigert) — *identical* for square and triangular (verified
to `X=2·10⁵`: both grow through highly-composite argmaxes). Hence the unit-distance count of any
`n`-point 2D-lattice piece is `≤ (max_D r_Q(D))·n/2 = n^{1+o(1)}`, and `ℤ²` achieves
`n^{1+c/loglog n}` (Erdős). So:

> **Theorem (lattice exponent ceiling).** Every 2D lattice has `n^{1+o(1)}` unit distances, with
> exponent `1+Θ(1/loglog n)` — *independent of kissing number, discriminant, or class number*.
> A point set with unit-distance exponent **bounded away from 1** cannot be a 2D lattice.

## Why the bridging spectrum collapses — and what the real axis is

Among 2D lattices = rank-2 ℤ-modules = imaginary quadratic orders, **both** candidate "bridge"
parameters are flat:
- the **kissing** parameter is universally `≤6` (degree-2 root-of-unity cap), and
- the **exponent** parameter is universally `1+o(1)` (degree-2 = 1-dimensional divisor sum).

There is no 2D lattice with `κ>6`, and none with exponent `>1+o(1)`. So **nothing sits "between"
the triangular lattice and the CM field on either axis among 2D lattices.** The genuine axis is the
**rank of the additive group** (equivalently the **degree of the CM field**):

| object | rank | roots of unity (`=1` layer) | popular-norm exponent |
|---|---|---|---|
| generic rhombic lattice | 2 | 2 | `1+o(1)` |
| square `ℤ[i]` | 2 | 4 | `1+o(1)` |
| triangular `ℤ[ω]` | 2 | **6 (degree-2 max)** | `1+o(1)` |
| **CM field, degree `≥3`** | `≥3` | `> 6` (more roots of unity) | **(claimed) `>1`** |

A rank-`≥3` finitely generated subgroup of `(ℝ²,+)` is **2D-realizable** (it embeds in the plane)
but is **not a lattice** — its closure is dense; finite pieces are exactly the CM/grid-beating
constructions. So the honest answer to the probe:

> **There is no 2D-realizable *lattice* strictly between `ℤ[ω]` and the CM field; the only thing
> "between" them is the jump from rank 2 to rank ≥3.** A rank-2 group is doubly capped (kissing 6,
> exponent `1+o(1)`); escaping *either* cap requires more roots of unity, which in degree 2 is
> impossible (max 6) — you must raise the degree. The CM field is the *first* group past the cap,
> with no continuum in between. "Beats `3n` at moderate `n`," by contrast, needs **no** new group:
> `ℤ²` does it at `n=144`.

## Honest status — theorem vs. the external claim

- **PROVED / classical, here assembled & verified exactly:**
  (1) `κ≤6` for 2D lattices, `=6` iff hexagonal; `κ = #`roots of unity of the order.
  (2) popular-norm beats `3n` inside `ℤ²` (`κ=4`) at `D=5`; crossover `n=144`. Triangular `D=7`,
  `n=64`. Divisor closed forms `4(d₁−d₃)`, `6(d₁−d₂)` checked vs brute force `D≤120`.
  (3) lattice unit-distance exponent `= 1+o(1)` for all 2D lattices (divisor-function rate, verified
  to `X=2·10⁵`); hence a `>1` exponent forces rank `≥3`.
- **EXTERNAL / NOT verified here:** the CM construction's claimed `n^{1.014}` (seed: "Sawin/OpenAI
  disproof"). I neither prove nor refute it; I only **place** it: *if* such an exponent holds it is
  necessarily a rank-`≥3` (degree-`≥3` CM) phenomenon, consistent with S699's "field-norm layer,"
  and it is **not** reachable by any 2D lattice.

## Connections (the additive face, again)

- **S699/HYP-2257 sharpened, one correction.** S699's "`=1` layer = roots of unity, capped at
  `κ≤6` for 2D lattices" is exactly right and now pinned: `κ = #`units of the order. But its
  bridging-probe ("a 2D group between `κ=6` and CM beating `3n`") rests on conflating the **kissing**
  `=1` layer (what's capped at 6) with the **popular-norm** layer (what beats `3n`). Beating `3n`
  is the popular layer and needs no new group; exceeding **exponent `1+o(1)`** is the rank jump and
  needs the CM field. There is no 2D-lattice middle ground on either axis.
- **Eisenstein `π/3` / `Cl₂(π/3)` (S599u) and the doubling theme.** The two distinguished 2D
  lattices are exactly the two orders with extra roots of unity — `ℚ(i)` (4) and `ℚ(ω)` (6) — the
  same cyclotomic anchors (`ζ₄, ζ₆`) the project keeps hitting (THM-403 `(2n−1)`-th roots; the
  `π/3` Eisenstein angle). `κ=6` is the Eisenstein hexagon; `κ=4` the Gaussian square; `6` is the
  hard ceiling because `ζ₆` is the largest root of unity in a quadratic field.
- **LRC additive face (the worry-set side).** S699's unification — both LRC pair-sums (mod `2n−1`)
  and unit distance (norm-1) are "maximal root-of-unity packings, capped geometrically and escaped
  arithmetically" — gets a precise dictionary here: **geometric cap = roots of unity of a degree-2
  order (≤6); arithmetic escape = divisor function (popular norm, exponent `1+o(1)`) then rank/degree
  (CM, exponent `>1`).** Two distinct escapes, not one — the cluster's earlier picture folded them.

## Handoff (next explorer)

1. **Make the "rank-2 ⟹ exponent `1+o(1)`" ceiling a clean stated theorem** (it is essentially
   Erdős's lattice bound + Wigert; worth a one-page write-up as the rigorous backbone separating
   2D lattices from the CM construction).
2. **Read the actual CM/grid-beating construction** (seed cites "Sawin/OpenAI"; verify the real
   statement and exponent before the cluster builds on `n^{1.014}` — it is currently an *unverified
   external number* in our notes). Confirm it is rank `≥3` / degree `≥3`, as this reflection predicts.
3. **Push the units=kissing identity the other way for LRC:** the worry-set "geometric floor" is the
   `(2n−1)`-th roots (THM-403); is there a "popular-norm" analogue on the LRC side — a *single*
   pair-sum residue hit by many runner pairs — that beats the floor the way `r_Q(D)>6` beats `3n`?
   (This would be the LRC mirror of the kissing-vs-popular split.)
