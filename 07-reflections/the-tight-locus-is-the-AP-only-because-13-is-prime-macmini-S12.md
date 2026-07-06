# The tight locus is the AP *only because 13 is prime* — and the induction must not descend through composite levels

*mac-mini-2026-07-06-S12 (HYP-4382). Owner: make LRC(14) a FULL theorem, not a
conditional submission; integrate the fleet during downtime. This note is the
session's structural finding — it does not close the theorem, but it corrects the
route the fleet is formalizing and pins the tight-locus half of
`TightLooseDichotomy` as clean at the prime, so the remaining crux is isolated.*

## The finding (verified)

Enumerate every gcd-1 covering configuration and find the **tight locus**
`{S : M(S) = 1/n}` (the extremal loneliness value LRC(n) attains):

| n | n prime? | tight families | non-AP tight families |
|---|----------|----------------|------------------------|
| 6 | composite (2·3) | **2** | **`{1,3,4,5,9}`** (M = 1/6, covering, gcd 1, *not* a dilated AP) |
| 7 | **prime** | 1 | none — **the AP is unique** |

`M({1,3,4,5,9}) = 1/6` is confirmed two independent ways: the exact profile
solver AND a brute-force fine-grid maximum (`Q = 4620`). It is genuinely a second
tight family, not a solver artifact
(`lrc_peel_lemma_macmini_S12.py`, `lrc_prime_tight_locus_macmini_S12.py`).

**The mechanism.** At `t = 1/6` the residues `{1,3,4,5}` all sit at safe
distances (`1/6, 1/2, 1/3, 1/6`), and `9 ≡ 3 (mod 6)` lands on the *same* safe
distance as `3`. Composite `6` has non-units `2,3,4`, so a tight family may
**skip a residue class (2) and repeat another (3)** yet stay lonely. At a
**prime** modulus every nonzero residue is a **unit**, so tightness forces the
*full* nonzero residue system `{1,…,n-1}` — the AP. This is exactly the content
of the corpus's own `residue_pinning_13` (GREEN): "at the **PRIME** modulus every
nonzero residue is a unit ⇒ all twelve classes are forced." The prime hypothesis
there is not a convenience; **it is the whole reason the pinning is true.**

## Why this corrects the route

The fleet is formalizing an **AP-tower / rigidity-ladder induction** for the
tight side (opus-S103 AP-11 protection; opus-S104 divisor protection; kps-S13
recursive AP tower; my own S59 AP-base rigidity). The inductive step is:

> a 12-family with a **tight 11-subfamily** — which, *by the level-below rigidity*,
> is a **dilated AP `{1,…,11}`** — completes to the AP or is loose (AP-completion).

The hidden assumption is **"tight ⟹ dilated AP" one level down**. That level is
`n = 12` (the 11-subfamily is an 11-speed = 12-runner configuration), and **12 is
composite**. The `n = 6` witness `{1,3,4,5,9}` is the small-scale proof that
composite-level rigidity **fails**: a tight 11-subfamily at composite 12 need not
be a dilated AP, so the AP-completion step has nothing to complete. **The
descending induction is not justified through a composite intermediate level.**

The escape is that **we never needed the descent.** `n = 13` is **prime**, and
`residue_pinning_13` acts **directly** on the 12-set: `M = 1/13` (tight) ⇒ no
runner is `≡ 0 (mod 13)` (a 13-multiple would sit on `0` at the `a/13` witness) ⇒
residues are the full system `{1,…,12}`. No 11-subfamily, no descent through
composite 12. The tight-locus half of `TightLooseDichotomy` is **clean at the
prime** — pinning (formal) plus the lift-rigidity M-minimizer (sibling S11's
HYP-4362, empirical) — and the AP-tower induction is a heuristic scaffold, not the
proof.

## The synthesis: why the AP is the *unique* scale-flow fixed point

opus-S48 (OPEN-Q-108 R2, HYP-4013) established to evidence standard that the
gap-emptiness (G) reduces by a **scale flow**: unbounded-height bounded-ratio
clusters **renormalize down to their difference core**, and the flow **contracts
to the compact AP fixed point**, where the density floor is minimized
(`F_7 = 0.0507` at range 6, rising to `0.099` at range 150, all `≥ 1/36`). The AP
is the fixed point because *differences of an AP are an AP*.

This finding says **why that fixed point is unique**: at the prime the tight locus
is a single point (the AP), so the contraction has nowhere else to land. The
composite artifacts `{1,3,4,5,9}` are fixed points of the flow **at composite
scales** — they are exactly the spurious attractors the descent would hit, and
exactly the ones the prime-13 pinning forbids at the top level. **Primality is
what makes the fixed point unique and the flow's target the AP.**

## What "full theorem" still needs (the honest map)

`TightLooseDichotomy` = tight-locus rigidity + gap-emptiness (G). After this
session:

1. **Tight-locus rigidity** (`M = 1/13 ⇒ dilated AP`): `residue_pinning_13`
   (FORMAL) + lift-rigidity **M-minimizer** (sibling HYP-4362, empirical — the
   canonical lift minimizes M over each mod-25 transversal profile). *Open piece:
   the M-minimizer proof.* CLEAN modulo that — no descent, no composite hazard.
2. **Gap-emptiness (G)** — the general spectral-gap conjecture HYP-2052 at
   `n = 13`, restricted to covering-compressed families. TRUE empirically;
   **not** finite-modulus-decidable (MISTAKE-110). The live route is opus-S48's
   scale-flow contraction to the AP + a **positive density floor** (`≥ 1/36`),
   both at evidence standard. *Open pieces: (a) rigorous scale-flow contraction
   rate, (b) rigorous density floor (the AP-is-min-star-discrepancy statement).*

`CornerLonely` (the loose-branch corner) is the third named hypothesis, handled by
THM-619/620's band sweep — a finite structured sweep needing a uniform argument.

**The bottom line for the owner:** a *full* theorem is not a formalization task —
it requires closing an **open conjecture** (the spectral gap at n = 13). The
tight side is clean at the prime; the open crux is (G), and the fleet's most
advanced route is the scale-flow/density-floor argument, which needs two pieces of
rigor. This session removed a false lead (the composite-descent induction) and
pinned the AP fixed point's uniqueness to primality.

## Pointers

- `lrc_peel_lemma_macmini_S12.py` / `.out` — the peel-lemma test that surfaced the
  n=6 obstruction (3 near-tight configs with no AP-tight peel) vs n=7 (none).
- `lrc_prime_tight_locus_macmini_S12.py` / `.out` — the prime/composite census.
- `residue_pinning_13` (LRCResiduePinning.lean) — the formal embodiment: prime ⇒
  full residue system ⇒ AP.
- OPEN-Q-108 (R2), HYP-4013 (opus-S48) — the scale-flow / density-floor route to (G).
- HYP-2052 (oracle-S552) — the general spectral gap; my S59 AP-base rigidity;
  opus-S103/S104 AP/divisor protection; sibling HYP-4362 lift M-minimizer.
