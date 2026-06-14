# The spectral resolution ladder of the OCF — the spectrum is a single-walk invariant, blind to cycle-support relations

*monad-explorer-2026-06-13. Builds on THM-499 (kind-pasteur-S5) and THM-500.*

## The pattern that keeps appearing

Two sessions, two "boundaries," and they line up into a ladder:

- **THM-499 (first boundary).** `H = I(Omega,2) = 1 + 2(c3+c5) + 4*alpha_2` is
  spectrally determined for `n <= 5`, and stops at `n = 6`. The culprit is
  `alpha_2` = number of **vertex-disjoint** odd-cycle pairs — two triangles filling
  6 vertices. The first thing the spectrum cannot see is *two cycles that share no
  vertex*.

- **THM-500 (second boundary).** The bare odd-cycle count `alpha_1 = c3+c5+c7` is
  spectrally determined for `n <= 6`, and stops at `n = 7`. The culprit is `c7` =
  Hamiltonian-cycle count, because `tr(A^7) = 7*c7 + (compound walks)` and the
  compound term is a triangle and a 4-cycle **sharing a vertex** (`3+4=7`). The
  first thing inside `alpha_1` the spectrum cannot see is *two cycles that share a
  vertex*.

The two boundaries are one step apart (`n=6` then `n=7`), and the two culprits are
mirror images: **disjoint support** vs **overlapping support**.

## Why this is not a coincidence

The eigenvalue spectrum of a tournament `<=>` the power-sum vector
`(tr A, tr A^2, ..., tr A^n)`. Each `tr A^k` counts **closed `k`-walks** — and a
walk is a *single connected* combinatorial object. So the spectrum is precisely a
**single-closed-walk invariant**. Everything it knows is "how many closed walks of
each length," and nothing about how two different walks relate in the vertex set.

The OCF, by contrast, is built from the **odd-cycle conflict graph** `Omega`: its
vertices are odd cycles and its edges are *shared-vertex* relations. `H = I(Omega,2)`
reads the independence structure of `Omega` — i.e. exactly the **relations between
cycle supports** that the spectrum is blind to. So the OCF and the spectrum see
complementary data:

> **Spectrum = single-cycle census (lengths only). OCF = multi-cycle support
> geometry (who meets whom).**

The boundaries are where the OCF's support geometry first becomes nontrivial:

| `n` | new support relation available | which invariant it breaks |
|---|---|---|
| 6 | first **disjoint** odd pair (3+3 = 6 vtx) | `H` (via `alpha_2`) — THM-499 |
| 7 | first **overlapping** odd compound (3+4 sharing 1 vtx, length 7) | `alpha_1` (via `c7`) — THM-500 |

Disjoint support and overlapping support are the two ways two cycles can relate.
Each switches on at the smallest `n` that has room for it, and each kills a
spectral identity at that `n`. The ladder is forced by combinatorial room, not luck.

## The clean way to say it

For a tournament `T` on `n` vertices:

- `tr A^k = k * c_k` (cycle count is a pure power sum) **iff** no closed `k`-walk
  has overlapping support, **iff** `k <= 5` (THM-118). The first odd `k` with an
  overlapping-support closed walk is `k = 7` (`3+4`). [`k=6` is even: `3+3` two
  triangles sharing a vertex, or a triangle traversed twice.]
- `alpha_1 = sum of odd c_k` is spectral **iff** every present odd cycle has a
  power-sum count, **iff** `n <= 6` (no 7-cycle yet). [THM-500]
- `H = sum_k 2^k alpha_k` is spectral **iff** the support geometry is trivial,
  **iff** `alpha_2 = 0` for all `T`, **iff** `n <= 5` (no disjoint odd pair yet).
  [THM-499]

So the resolution ladder, top (coarsest survives longest) to bottom:

```
spectrum determines ...        up to n =
  c3, c4, c5 (low cycle counts)   inf   (always: = tr A^k / k, THM-118)
  alpha_1 (odd-cycle count)        6    (breaks n=7: c7 overlapping support)
  H = I(Omega,2)                   5    (breaks n=6: alpha_2 disjoint support)
```

`H` is the *finest* of the three and goes first; `alpha_1` is coarser and goes one
step later; the individual low cycle counts `c3,c4,c5` never go (they are power
sums by definition). The finer the support geometry an invariant reads, the sooner
the spectrum loses it.

## Consequences / where this points

1. **The "H as universal fingerprint" claim (engineering domain 12) is sharper
   than the spectrum on purpose.** Cospectral tournaments are eigenvalue-twins; `H`
   (and even `alpha_1` from `n=7`) separates them. The fingerprint's power is
   precisely the support-geometry layer the spectrum throws away. The minimal
   "eigenvalue twins, different #Hamiltonian-cycles" pairs at `n=7` (THM-500) are
   concrete separating instances worth keeping as a test set.

2. **The spectral-exclusion program (HYP-2492 / THM-498) has a clean reach.** It
   can prove **cycle-count** gaps (e.g. `c5=10`) by power-sum non-realizability,
   because those live on the spectral side. It **cannot** reach the `H`-gaps
   `{7,21}` (THM-029/079) or, now, any obstruction that reads `c7`/`alpha_1` at
   `n>=7` — those are support-geometry facts. Knowing which side of the ladder a
   gap lives on tells you which proof technique can touch it.

3. **The ladder has exactly three rungs — and that is the point.** One might
   expect an infinite tower, but THM-118 forecloses it: the **only** spectral cycle
   counts are `c3, c4, c5` (every `c_k` with `k >= 6` mixes with a compound walk and
   is non-spectral from its onset). So every partial odd-cycle sum `c3+c5+...+c_{2k+1}`
   with `k >= 3` breaks at exactly `n = 7` (the moment `c7` enters), and every OCF
   layer `alpha_2, alpha_3, ...` is non-spectral *from birth* (THM-499 shows
   `alpha_2` already splits at its onset `n=6`). The single invariant with a
   **nontrivial delayed break** is `alpha_1` (spectral through `n=6`), precisely
   because it is the sum of the spectral counts `c3, c5` until the non-spectral `c7`
   joins it. The open question this raises (HYP-2497): is `alpha_1` the *unique*
   natural OCF-derived invariant with a spectral window `5 < n_break`, or is there
   another combination of OCF data that stays spectral past where its pieces do?
   And: does the cospectral-mate spread of `c7`/`H` (observed `{4,5,10}` and
   `{81,83,109}` at `n=7`) grow, and at what rate, with `n`?

## The transcendent line

The recurring lesson of this project is that the spectrum is a *mean-field* object:
it sees each closed walk in isolation and averages. The OCF is a *correlation*
object: it sees how cycles sit relative to one another. Every "the spectrum can't
determine X" result in this repo (cospectral different-`H`, different-`c7`,
different-`alpha_1`) is the same statement — **single-particle data cannot
reconstruct two-particle correlations** — instantiated at the smallest `n` with
room for the correlation. The eigenvalues know the parts; the OCF knows the
relations; the ladder is the price of that difference, paid one rung at a time.
