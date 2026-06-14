---
id: THM-500
title: The second spectral boundary — the total odd-cycle count alpha_1 is spectrally determined iff n <= 6; it breaks at n=7 because c7 (Hamiltonian-cycle count) is not a power sum
status: PROVED (n<=6 direction: rigorous via THM-118; n=7 direction: explicit cospectral witnesses, brute-force double-checked)
source: monad-explorer-2026-06-13
depends_on:
  - THM-118   # c_k = tr(A^k)/k for k <= 5 (low cycle counts ARE power sums)
  - THM-499   # the (first) spectral-reframe boundary: H non-spectral from n=6 via alpha_2
related:
  - THM-498   # cycle-spectrum gaps = power-sum exclusions
  - THM-029   # H-forbidden 7
  - THM-079   # H-forbidden 21
---

# THM-500 — the second spectral boundary

THM-499 located the **first** spectral boundary: the OCF Hamiltonian-path count
`H = I(Omega,2)` is spectrally determined for `n <= 5` and breaks at `n = 6`,
because `H = 1 + 2(c3+c5) + 4*alpha_2` and the disjoint-odd-cycle-pair count
`alpha_2` is the first non-spectral ingredient. Crucially, **at n=6 the bare
odd-cycle count `alpha_1 = c3+c5` is still spectral** (both `c3=tr A^3/3`,
`c5=tr A^5/5`, THM-118). This file finds the **second** boundary, one step higher:
where the spectrum loses the ability to even *count* the odd cycles.

## Statement

Let `sig(T) = (tr A, tr A^2, ..., tr A^n)` (`<=>` characteristic polynomial of the
adjacency matrix `A` `<=>` the eigenvalue spectrum). Let `alpha_1(T)` = total number
of directed odd cycles.

> **THM-500.** `alpha_1` is a function of `sig(T)` **iff `n <= 6`.**
>
> 1. **(n <= 6, spectral — rigorous.)** For `n <= 6` the only odd cycles have length
>    3 or 5, so `alpha_1 = c3 + c5 = tr(A^3)/3 + tr(A^5)/5` (THM-118). Hence
>    `alpha_1` is a fixed function of the spectrum.
> 2. **(n = 7, non-spectral — explicit certificate.)** At `n=7` the odd cycles are
>    3-, 5-, and 7-cycles, so `alpha_1 = c3 + c5 + c7` with `c7 =` #directed
>    Hamiltonian cycles. **`c7` is NOT spectrally determined:** there exist
>    cospectral tournaments with different `c7`, hence different `alpha_1`.
>    Witness (brute-force verified): two valid 7-tournaments with identical
>    `sig = (0,0,30,68,90,360,1204)` and `c3=10, c5=18`, but `c7 = 4` vs `c7 = 5`,
>    giving `alpha_1 = 32` vs `33` and `H = 81` vs `83`. A third tournament in the
>    same spectral class has `c7 = 10` (`alpha_1 = 38`, `H = 109`).

Explicit witnesses (upper-triangle arc orientations) are in
`05-knowledge/results/second_spectral_boundary_n7_monad.out`; the `c7=4` / `c7=5`
pair is reproduced and independently re-checked (all `7!` permutations) at the
bottom of this file.

## The mechanism: tr(A^7) splits between c7 and a compound-walk term

`tr(A^7)` counts closed 7-walks `= 7*c7 + (non-simple closed 7-walks)`. The
shortest non-simple closed 7-walk uses a triangle and a 4-cycle **sharing a
vertex** (`3 + 4 = 7`), so the correction is a *compound* (overlapping-support)
quantity, not a single-cycle count. For the witness pair `tr(A^7) = 1204` is fixed
by the spectrum, but it partitions as

| witness | `c7` | `7*c7` | correction `tr(A^7) - 7*c7` |
|---|---|---|---|
| A | 4 | 28 | 1176 |
| B | 5 | 35 | 1169 |

The spectrum fixes the **sum** `7*c7 + corr = 1204`; it does **not** fix the split.
That is exactly why `c7` (and hence `alpha_1`) escapes the spectrum at `n=7`,
whereas `c3 = tr A^3/3` and `c5 = tr A^5/5` cannot (no shorter closed walk pads to
length 3 or 5 in a tournament — no loops, no 2-cycles, THM-118).

## The exact correction identity (bridges codex HYP-2498 / OPEN-Q-093)

The compound term is not just bounded below — it has an exact closed form, the odd
analog of codex's `tr(A^6) = 6*c6 + 3*c3 + 6*p33_meet`:

> **`tr(A^7) = 7*(c7 + TQ)`**, equivalently **`c7 = tr(A^7)/7 - TQ`**,
>
> where `TQ` = number of (directed-triangle, directed-4-cycle) pairs with
> **overlapping** vertex support (`|V(tri) ∩ V(4cyc)| ∈ {1,2,3}`).

Verified exact on `600/600` and `120/120` random `n=7` tournaments, 0 exceptions
(`04-computation/trace7_overlap_correction_monad.py`, `.out`). The fit is striking:
splitting `TQ` by overlap size `(tq1,tq2,tq3)` gives the SAME coefficient `7` for
each (least-squares residual `0`), so all overlapping triangle–4-cycle pairs
contribute equally — each such figure (a triangle and a 4-cycle glued on a shared
vertex/edge) is a closed 7-walk traversed in 7 rotations, exactly like a simple
7-cycle.

This **pinpoints the non-spectral carrier**: `tr(A^7)` is spectral (a power sum),
so within a cospectral class `c7 = tr(A^7)/7 - TQ` varies **iff `TQ` varies**.
THM-500's `c7`-split is precisely a `TQ`-split — the overlap count `TQ` is the
n=7 odd-cycle analog of codex's intersecting-triangle-pair count `p33_meet`. So the
trace-correction engine (HYP-2498) and the spectral boundary (THM-500) are two
readings of one object: the **support-overlap geometry** the power sums cannot
resolve.

## The two boundaries form a one-step-offset ladder

The OCF `H = sum_k 2^k alpha_k` and its first layer `alpha_1` lose spectrality at
**different** `n`, for **different** reasons:

| invariant | spectral for | breaks at | why it breaks |
|---|---|---|---|
| `H = I(Omega,2)` | `n <= 5` | `n = 6` | `alpha_2` onset — two **vertex-disjoint** odd cycles (THM-499) |
| `alpha_1` (odd-cycle count) | `n <= 6` | `n = 7` | `c7` onset — a 7-cycle whose count `tr(A^7)` mixes with a **compound** walk |

So `H` goes non-spectral **one step before** the bare odd-cycle count does. The
offset has a unifying cause: the eigenvalue spectrum is a **single-closed-walk**
invariant (the power sums `tr A^k`); it is blind to *relations between cycle
supports*. Both failures are such relations — `alpha_2` is **disjoint** support
(two cycles sharing no vertex), the `c7` correction is **overlapping** support (a
triangle and a 4-cycle sharing a vertex). Disjoint and overlapping support are the
two ways multi-cycle information hides from the spectrum, and they switch on at
`n=6` and `n=7` respectively. (Reflection: `07-reflections/the-spectral-resolution-ladder-of-the-ocf.md`.)

## Corollary (formula extension, verified)

The seed formula `H = 1 + 2(c3+c5) + 4D` (THM-499, `n<=6`) extends exactly to `n=7`:

> `H = 1 + 2(c3 + c5 + c7) + 4*DTP`,  `DTP = #vertex-disjoint triangle pairs`,

verified `4000/4000` random `n=7` tournaments, 0 exceptions
(`04-computation/second_spectral_boundary_n7_monad.py`). At `n=7` the only
disjoint odd-cycle pair is triangle+triangle (6 of 7 vertices; `alpha_3` needs 9),
so this is the exact OCF truncation, with **two** non-spectral ingredients now —
`c7` (inside `alpha_1`) and `DTP = alpha_2`.

## Status / reproduction

- n<=6 direction: rigorous (THM-118 + cycle-length bound).
- n=7 direction: explicit cospectral witnesses; both confirmed valid tournaments,
  cospectral, with `c7 = 4` vs `5`, `H = 81` vs `83`, by an independent brute-force
  over all `7!` orderings.
- Sampling: `168` distinct spectral classes hit in `80000` uniform random
  n=7 tournaments; `46` of them split `c7` (and `alpha_1`, and `DTP`). Stable.

**Artifacts:** `04-computation/second_spectral_boundary_n7_monad.py`,
`05-knowledge/results/second_spectral_boundary_n7_monad.out`.

### Explicit witness pair (arcs as oriented pairs i->j)

```
sig = (0,0,30,68,90,360,1204),  c3=10, c5=18  (both)

c7=4 (H=81):
 0->1 0->2 0->3 4->0 5->0 6->0 2->1 3->1 4->1 1->5 6->1 3->2
 4->2 2->5 6->2 3->4 3->5 6->3 5->4 4->6 5->6
c7=5 (H=83):
 1->0 0->2 3->0 4->0 5->0 6->0 2->1 1->3 1->4 1->5 6->1 3->2
 4->2 5->2 2->6 4->3 5->3 3->6 4->5 6->4 5->6
```
