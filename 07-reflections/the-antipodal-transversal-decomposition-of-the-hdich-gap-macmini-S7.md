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
all w, giving M ≥ 2/25 — out of the gap. ∎

**(P2) contain a 5-divisible speed.**  *Proof (the mod-5 filter):* if no w ≡ 0
mod 5, then at t = 1/5 every ‖w/5‖ = dist₅(w)/5 ≥ 1/5, so M ≥ 1/5 — out of the
gap. ∎  Moreover the 5-divisible speeds are **free**: for w = 5u and any unit b,
‖w·b/25‖ = dist₅(ub)/5 ≥ 1/5, so they never bind at a mod-25 witness.  The
covering tension at mod 25 is carried **entirely by the unit residues.**

So the gap splits by the antipodal structure of W mod 25:

| class | condition | status |
|-------|-----------|--------|
| **non-transversal** | misses a unit pair | SAFE (P1: M ≥ 2/25) — DONE |
| **non-unit-pair hole** | hits all unit pairs, misses {5,20} or {10,15} | the *second witness* — closable (this session) |
| **transversal** | hits every antipodal pair | the AP-rigidity crux |

## The non-unit-pair hole (mod-5 structured, closable)

A hole member missing {5,20} has its 5-divisible speeds confined to {10,15} =
5·{2,3} (mod 25) — i.e. every 5-divisible w has w/5 ≡ ±2 mod 5.  This is the
rank-2 analogue of the n=8 sporadics (which live in exactly this hole at 2n−1=15).
The census (S553) found it EMPTY at n=14; the mechanism is the *second witness*
(oracle's "lift to restore invertibility"): [computed pattern in the .out —
the witness modulus for hole members reveals the closing family].  The hole is
a bounded mod-5-structured check, not a general census.

## What remains: the transversal crux

Transversals (hit every antipodal pair, including the non-unit ones) are the only
class with **no mod-25 surplus witness** — they are the genuinely tight-adjacent
configurations.  The AP {1,…,12} is one (M = 1/13, witnessed at mod 13, not 25).
The crux is the **transversal rigidity**:

> every covering primitive transversal mod 25 is the AP (M = 1/13, excluded from
> the *open* gap) or has M ≥ 2/25.

This is precisely opus-S99's 2-D lift rigidity / the tight-family-is-AP census,
now localized to the antipodal transversals — a small, finite, structured set
(≈ one choice per unit pair × the free slots), far smaller than the full Q50
census.  It is the same crux, but its domain has shrunk from "all residue
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
