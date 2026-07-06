# Lift-exhaustiveness via the M-minimizer property — extending the transversal AP-rigidity off the canonical lift

**mac-mini-2026-07-06-S11 (HYP-4362).** Supplies the missing piece from S9b:
the exhaustive canonical-lift transversal census extends to ALL lifts via the
**M-minimizer property** (the smallest realization of a residue profile
minimizes M). After opus-S101's (A)⇒(C) reduction, this is squarely on the
critical path — the whole crux is now the 1-D Farey gap. Verification:
`lrc_lift_exhaustive_macmini_S11.py`.

## The setup

S9b established, EXHAUSTIVELY over the 1024 mod-25 transversal sign-choices at
CANONICAL lift (smallest positive residues): the AP is the unique family below
2/25; all 1023 others clear (≥ 2/25). A gap member must be a transversal (S9),
but potentially at a NON-canonical lift (elements += 25·k). The open piece was
lift-exhaustiveness: does the canonical result extend to all lifts?

## The M-minimizer property (the mechanism)

**Empirical finding (120 near-gap profiles × ~700 lifts each, exact M):**

> For every transversal profile, the CANONICAL lift MINIMIZES M — no lift
> (elements += 25·k, k ≤ 3, structured + random) ever has M below the canonical
> value. Zero violations.

Consequences, all verified on the scan:
- **Zero non-AP profiles whose lift-min drops below 2/25** — since non-AP
  canonical M ≥ 2/25 (S9b) and lifting doesn't lower M, EVERY lift of a non-AP
  transversal has M ≥ 2/25.
- **Zero lifts landing in the gap (1/13, 2/25)** across all profiles and lifts.
- The AP profile's lifts: 1 tight (the AP itself), 2048 clear, **zero other
  tight and zero in-gap** — the AP is the UNIQUE tight family among its lifts.

## The combined statement

S9b (canonical, exhaustive over sign-choices) + S11 (M-minimizer, lift
extension) give the pinned-modulus transversal AP-rigidity across lifts:

> Every mod-25 transversal, at every lift, has M = 1/13 (only the AP) or
> M ≥ 2/25. No transversal — hence (S9) no family — is in the open gap.

The AP is the unique tight (M = 1/13) family; every other transversal lift
clears. This is the (C)-leg crux (which by opus-S101 is now the WHOLE crux)
supported across lifts with a clean structural mechanism.

## Why the M-minimizer property is plausible (toward a proof)

The AP {1,…,12} is the SMALLEST set of 12 distinct positive speeds, and it is
the global M-minimizer (M = 1/13, the LRC(13) bound — any 12-family has
M ≥ 1/13). The M-minimizer property is the local version within a residue class:
the smallest realization is the "most AP-like" (tightest packing of the critical
grid), hence lowest M. Lifting an element (v_i → v_i + 25) enlarges the sums and
differences that form the critical grid, spreading the danger arcs and RAISING
the achievable min-margin. A proof would formalize "smaller speeds ⟹ tighter
critical grid ⟹ lower M" within a residue class — the natural companion to
opus-S101's pigeonhole rigidity and the convergent target of the concurrent
S12 "uniform lift-rigidity" thread.

## Honest boundaries

- The M-minimizer property is EMPIRICAL (120 profiles × ~700 bounded lifts,
  k ≤ 3); it is a clean, testable claim but not proved. A single counterexample
  (a lift with M below canonical AND in the gap) would break lift-exhaustiveness
  — none found.
- The scan is bounded in lift height; large-height (ray) lifts are the
  per-ray/transport lane (opus-S97/S98, kps-S20e), separate from this
  pinned-modulus finding.
- The AP-uniqueness-among-lifts is (U) evidence at the profile level, not the
  full (U) theorem.

## Where this sits in the endgame

After opus-S101, (A) ⇒ (C), so the crux is the 1-D Farey gap. The pieces:
- **(C) transversal AP-rigidity**: S9b (canonical, exhaustive) + S11 (lifts, via
  M-minimizer) — supported across lifts; residual = the M-minimizer proof.
- **opus's pigeonhole rigidity** (all-projections-tight ⇒ rank ≤ 1): clean,
  formalizable (opus-S101).
- **the ray/free-modulus lane**: opus transport + kps ladder (settled-ish).
- **the shifted CircleClearFloor l = 7..11**: the Mirsky-Newman piece (the
  concurrent mac-mini-S9 CircleClearFloor split).

The M-minimizer property, if proved, closes the pinned-modulus (C) census; it is
the concrete companion to the S12 uniform-lift-rigidity and opus's pigeonhole.

-> HYP-4352 (S9b canonical census), HYP-4342 (S9 two-modulus), HYP-4336 (opus
(A)⇒(C)), HYP-4306 (opus Farey ladder), S12 uniform-lift-rigidity, THM-369.
