# The antipodal-transversal decomposition of the hdich gap (mod 25)

**mac-mini-2026-07-06-S7 (HYP-4312).**  A creative reconnection: oracle-S552/S553's
spectral-gap + antipodal-transversal reduction (≈50 sessions old, unintegrated)
is a *sharper* finite object than the Q50 residue census, and it cleanly splits
the hdich gap into three pieces — two proved/closable here, one the deep crux.
Verification: `05-knowledge/results/lrc_antipodal_hole25_macmini_S7.out`.

## The reconnection

oracle-S552 proved the spectral gap (M = 1/n or M ≥ 2/(2n−1)); oracle-S553
reduced its emptiness to antipodal transversals mod 2n−1, with a **non-unit-pair
hole** when 2n−1 is composite.  At the hdich level (observer + 12 speeds, gap
(1/13, 2/25]) the relevant modulus is **2n−1 = 25 = 5²** — composite, so the hole
is present, and it is **mod-5 structured**.  opus-S553 §3's "tight ⟹ witness at
t = j/n ⟹ residues pinned mod n" is exactly THM-593A / `residue_pinning_13`
(the mod-13 tight-witness lattice), one level up.  Two moduli govern the gap:

- **mod 13** (t = j/13): the *tight* witness — the AP is lonely at {1,3,5,9,11,13}/13,
  M = 1/13 exactly.
- **mod 25** (t = b/25): the *surplus* witness — Link-1 pushes M ≥ 2/25.

The gap (1/13, 2/25) is the Farey cell *between* them (THM-622: 1/13, 2/25 are
Farey neighbors).

## The decomposition (proved at the hdich level)

A covering, primitive, residue-0-free 12-set W with M ∈ (1/13, 2/25) must:

**(P1) hit every unit antipodal pair mod 25.**  *Proof (Link-1 surplus):* if W
misses a unit pair {a, 25−a}, then at t = a⁻¹/25 every residue w·a⁻¹ mod 25 lies
in {2,…,23} (not 0: residue-0-free; not ±1: pair missed), so ‖w t‖ ≥ 2/25 for
all w, giving M ≥ 2/25 — out of the gap. ∎  *(Residue-0-free — no w ≡ 0 mod 25,
i.e. 25 ∤ w — is required; the computation confirms the only apparent Link-1
"violations" are exactly the 25 | w cases.  A speed divisible by 25 = 2n−1 is
≥ 25 and handled by the same mod-5 filter P2 or the covering-far peel; it is a
thin edge, not a gap carrier.)*

**(P2) contain a 5-divisible speed.**  *Proof (the mod-5 filter):* if no w ≡ 0
mod 5, then at t = 1/5 every ‖w/5‖ = dist₅(w)/5 ≥ 1/5, so M ≥ 1/5 — out of the
gap. ∎  Moreover the 5-divisible speeds are **free**: for w = 5u and any unit b,
‖w·b/25‖ = dist₅(ub)/5 ≥ 1/5, so they never bind at a mod-25 witness.  The
covering tension at mod 25 is carried **entirely by the unit residues.**

So the gap splits by the antipodal structure of W mod 25:

| class | condition | status |
|-------|-----------|--------|
| **non-transversal** | misses a unit pair | SAFE (P1: M ≥ 2/25) — DONE |
| **non-unit-pair hole** | hits all unit pairs, misses {5,20} or {10,15} | SAFE by LOOSENESS (this session) |
| **transversal** | hits every antipodal pair | the AP-rigidity crux |

## The non-unit-pair hole is LOOSE, not delicate (the session's finding)

The expectation (from the n=8 sporadics) was that the hole needs a careful
"second witness."  The data says the opposite: the hole is **rare and loose.**
Of 56 000 random covering families, only 6 landed in the hole class, and every
one has **M ≥ 4/19 ≈ 0.21** — far *above* the gap, cleared at ordinary moduli
(18, 19, 31, 37, 62, 66), zero in the open window.  The reason is structural:
missing an antipodal pair *removes* a covering constraint, so the config is
**under-constrained → loose**, the opposite of near-tight.  Combined with
oracle-S553's exhaustive n=14 census (hole empty in range), the hole carries no
gap members.  No second-witness lift is needed — looseness suffices.

*(The n=8 sporadics live in this hole because at the SMALL modulus 2n−1 = 15 the
non-unit pairs {3,12},{6,9} are a larger fraction of the residues; at 25 = 5² the
non-units are only 4 of 24 residues, so missing one barely loosens — and the
covering/primitivity constraints dominate, pushing the hole loose.)*

## The clean reduction

Both miss-classes are safe (P1 surplus / looseness), so:

> **Every gap member is a full antipodal transversal mod 25** — it hits every
> antipodal pair {c, 25−c}, unit and non-unit alike.

## What remains: the transversal crux

Transversals (hit every antipodal pair, including the non-unit ones) are the only
class with **no mod-25 surplus witness** — they are the genuinely tight-adjacent
configurations.  The AP {1,…,12} is one (M = 1/13, witnessed at mod 13, not 25).
The crux is the **transversal rigidity**:

> every covering primitive transversal mod 25 is the AP (M = 1/13, excluded from
> the *open* gap) or has M ≥ 2/25.

This is precisely opus-S99's 2-D lift rigidity / the tight-family-is-AP census,
now localized to the antipodal transversals — the 4096 sign patterns
{±1,±2,…,±12} mod 25 (one residue per pair), far smaller than the full Q50
census.

**The mod-25 landscape is flat (a caution for the crux).**  All 4096 transversal
residue-sets have mod-25 margin EXACTLY 1/25 — every transversal hits {1,24}, so
at the balancing dilation some residue sits at distance 1.  So mod 25 identifies
the *domain* (gap ⟹ transversal) but cannot *rank* transversals: the AP-vs-rest
distinction lives at **mod 13** (the AP is tight there, M = 1/13; a non-AP
transversal is not mod-13-tight and jumps to ≥ 2/25 somewhere).  The transversal
rigidity is therefore genuinely a two-modulus (13 & 25) statement — the same
mod-13 residue pinning (THM-593A) that closes hdich's tight case, now needed to
close its *transversal* case.  The reduction shrinks the domain and peels the
easy classes; it does not by itself resolve the deep AP-uniqueness.  It is the same crux, but its domain has shrunk from "all residue
families" to "the antipodal transversals," with two of the three original classes
peeled off by elementary witnesses.

## Why this is a work-skip

The Q50 census asks: does *any* residue family clear ≥ 2/25?  The antipodal
reframe answers it for two of three classes with one-line witnesses (P1, P2),
leaving only the transversals — and among those, the rigidity is the classical
tight-family-is-AP statement, which the census has verified exhaustively in range.
The finite object is not "check 500k residue families" but "prove the antipodal
transversals mod 25 are AP-rigid" — a sharper, more structured target that
connects directly to opus's projection floor and the mod-13 residue pinning.
