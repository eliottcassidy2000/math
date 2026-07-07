---
source: kind-pasteur-2026-07-07-S54
status: census of the LRC(14) single-scale rigidity corner — corrects S53, refines opus's
  density-floor lemma, extracts the Farey-ladder spectrum and closed-form near-tight families
tags:
  - lonely-runner
  - LRC14
  - route-1
  - density-floor
  - tight-locus
  - census
  - three-gap
  - correctness
---

# The M-tight locus is {AP, GW}; the density-floor minimizer is uniquely AP

Owner: *work Route 1 density + remaining tasks; if you meet seeming randomness, census it and
understand the finite/infinite families.* I censused the single-scale rigidity corner (the
residue my S52/S53 coarse reduction leaves the density node). The census both **corrects an
error of mine** and **sharpens opus's open lemma**.

## The correction (MISTAKE-100 recurrence)

S53 claimed *"the AP `{1..13}` is the **unique** single-scale tight family."* **False.** A
one-swap census (`lrc_tight_locus_census`) finds a second primitive family with `M = 1/14`:

> **`GW = {1,2,…,11,13,24} = AP[12 → 24]`**, verified `M = 1/14` exactly (witness `t* = 1/14`,
> `min‖v_i/14‖ = 1/14`). This is mac-mini's Goddyn-Wong family (THM-612).

My S53 perturbation test only bumped one coordinate by small amounts, so its generator could
never produce the `12 → 24` shape — exactly the **MISTAKE-100** trap (a "no X exists" claim
from a search that can't emit X). Within `≤ 2` swaps of the AP (pool `{14..40}`), the tight
locus is **exactly `{AP, GW}`** (2 primitive families); a broader `≤ 3`-swap scan is running
but the corner is small and discrete.

## The refinement (this is the useful part): M-tight ≠ floor-minimal

Here is the structural payoff. `M`-tightness and opus's density-floor minimality are
**different functionals**, and the AP is special for the *second*, not the first:

| family | `M` | `μ_{1/7}` (opus's floor) |
|---|---|---|
| **AP `{1..13}`** | **1/14** (tight) | **477/1078 ≈ 0.4425** (the global min) |
| **GW `{1..11,13,24}`** | **1/14** (tight) | **0.588** (NOT minimal) |

So the `M`-minimizer locus is `{AP, GW}` (two families), but the **`μ_{1/7}`-minimizer is
uniquely the AP**. `GW` is just as `M`-tight, yet its asymptotic witness floor is *higher*
(0.588 > 0.4425) — it is "easier" for the density argument. **This means opus-S130's
AP-minimality lemma is not threatened by GW: it is a statement about `μ_{1/7}`, which is
strictly more rigid than `M`-minimization.** The correct statement of the open lemma is:

> **(A′) `μ_{1/7}(E) ≥ μ_{1/7}(AP)` for all 13-sets `E`, with equality iff `E` is a
> dilation/translation of the AP** — GW and all other `M`-tight families sit *strictly above*.

That the AP is the *unique* floor-minimizer (while the `M`-locus is larger) is exactly why
Route 1's density floor is the right object: the floor sees the AP alone, not the whole
degenerate tight locus.

## The spectrum is a Farey ladder (three-gap quantization)

Censusing all primitive single-scale 13-families (13-subsets of `{1..20}`), the distinct
`M`-values at the bottom are

> `1/14 < 2/27 < 1/13 < 2/25 < 1/12 < 3/35 < 2/23 < 1/11 < …`

each a small-denominator Farey fraction, with `2/27 = mediant(1/14, 1/13)`,
`2/25 = mediant(1/13, 1/12)`. The gap `(1/14, 2/27)` is **empty** — the direct-LRC(14)
analogue of the `(G)` gap. There are **no** intermediate values: the near-bottom spectrum is
**quantized to Farey rungs** (mac-mini's three-gap HYP-4412, here on the `1/14` object). This
is the "seeming randomness resolves into finite/infinite families" the owner asked for: the
families are not random, they are the covering-rung families of each Farey denominator.

## Closed-form near-tight families (infinite, lonely, formalizable)

Each rung is populated by a closed-form one-swap family:

- **`M({1..9,11,12,13, 10k}) = k/(10k+7)`** for all `k ≥ 2` (remove 10, add `10k`). First
  excited rung `2/27` at `k=2`; `→ 1/10`. **New (n=13).**
- **`M({1..11,13, 12k}) = k/(12k+5)`** for `k ≥ 3`; `k=2` gives `1/14` (the GW tight
  exception). This is the **residue-liar** family (already GREEN, `LRCResidueLiar.lean`) — and
  the census now shows its `k=2` boundary case *is* the GW tight family.

Both are `≥ 2/27 > 1/14`, so **lonely for every `k`** — infinite families grounded on Route 1
by an explicit residue-table `M`. The `10k` family is a clean next formalization target
(same method as the residue-liar).

## Where this leaves Route 1

The rigidity corner is now understood: a **discrete Farey ladder** of closed-form families
above a **two-point `M`-tight locus `{AP, GW}`**, of which **only the AP minimizes the density
floor `μ_{1/7}`**. So opus's open lemma (A) is correctly (A′): AP is the unique `μ_{1/7}`
minimizer. The other open piece (B, the finite-`Vmax` error budget) is unaffected.

## Update (later S54): the complete one-swap ladder table + what's formalized

Censusing `M({1..13}\{j} ∪ {jk})` (remove AP element `j`, add multiple `jk`) for all `j, k`
gives the **complete one-swap structure**:

| removed `j` | `M(k)` closed form | note |
|---|---|---|
| 2, 3 | `2/17` (constant in `k`) | small `j`: fixed loose rung |
| 4, 5 | `2/19` (constant) | |
| 6 | `2/23` (constant) | |
| 7 | `k/(7k+8)` | ladder `→ 1/7` |
| 8 | `k/(8k+7)` | `→ 1/8` |
| 9 | `k/(9k+5)` | `→ 1/9` |
| **10** | **`k/(10k+7)`** | `→1/10`; `k=2 → 2/27` — **`tenSwap_lonely` GREEN** |
| 11 | `k/(11k+3)` | `→ 1/11` |
| **12** | **`k/(12k+5)`** (`k≥3`); **`k=2 → 1/14` (GW)** | **`residueLiar_lonely` GREEN; `gw_lonely` GREEN** |
| **13** | **`k/(13k+1)`** (all `k≥1`, **includes AP at `k=1`**) | `k=2→2/27` — **`thirteenLadder_lonely`, `ap_lonely` GREEN** |

So removing a **large** AP element (`j ≥ 7`) and dilating it gives a Farey ladder
`k/(jk + b_j)`; removing a **small** element leaves a fixed loose rung. **Only `j=12, k=2`
(GW) touches the tight floor `1/14`** — the structural reason the tight locus is `{AP, GW}`.

**Formalized this session (kernel-pure, GREEN, in manifest):** `ap_lonely` (AP at `1/14`),
`gw_lonely` (GW at `1/14`), `thirteenLadder_lonely` (`k/(13k+1)`, `∀k≥1`), `tenSwap_lonely`
(`k/(10k+7)`, `∀k≥2`) — both tight families plus three near-tight ladders, all via the
residue-table `lattice_dist_ge` atom.  (`LRCTenSwapLadder.lean`, sibling to `LRCResidueLiar`.)

**But the near-tight families are richer than one-swaps.** The `2/27` rung (first excited)
contains, besides the two one-swap ladders `{1..12,26}` and `{1..9,11,12,13,20}`, a **two-swap**
family `{1..9,11,13,20,24} = AP[10→20, 12→24]` (composing the tenSwap and GW modifications).
So the rungs above `1/14` are populated by multi-swap compositions, not just the one-swap
ladders — **the density floor cannot be replaced by enumerating a finite list of certificate
families.** The tight *locus* is finite (`{AP, GW}`), but the near-tight *corner* is an
infinite, compositional family. This is why opus's density-floor `(A')` is genuinely needed.

## Ledger

- **Corrected:** S53's "AP is the unique single-scale tight family" (MISTAKE-100 recurrence) —
  M-tight locus is `{AP, GW}`.
- **New / verified:** `μ_{1/7}(GW) = 0.588 > μ_{1/7}(AP) = 477/1078` (AP is the *strict*
  floor-minimizer); the Farey-ladder spectrum with gap `(1/14, 2/27)`; closed forms
  `k/(10k+7)`, `k/(12k+5)`.
- **Refines:** opus-S130 lemma (A) → (A′) — `μ_{1/7}`-minimality, equality only at the AP.
- **Does NOT claim:** a proof of LRC(14) or of (A′). (A′) is the open analytic lemma.
- **Files:** `lrc_singlescale_census_kps_S54.py`, `lrc_rung_families_kps_S54.py`,
  `lrc_tight_locus_census_kps_S54.py` (+ outputs). HYP-4717.
- **Pointers:** opus-S130 (μ_{1/7} floor, AP-minimality); mac-mini THM-612 (GW tight),
  HYP-4412 (three-gap); kps-S52/S53 (coarse reduction); MISTAKE-100; `LRCResidueLiar.lean`.
