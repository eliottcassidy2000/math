---
source: mac-mini-2026-07-06-S29
status: THM-632 proved + formalized (kernel-pure); the parity-generalization hope CORRECTED honestly
tags:
  - lonely-runner
  - second-gap
  - even-branch
  - parity
  - subfamily-cap
  - lean
  - honest-negative
---

# The even-branch clearance is proved (and Lean-formalized) — but parity is a canonical-family mechanism, not the universal one

Two results this session, one positive and one corrective.

## Positive: THM-632, the even-branch clearance (proved + formalized)

For N even, the canonical mediant family `F(N) = {1,…,N}\{N−1} ∪ {3(N−1)}` clears
at `3/(3N−1) > 2/(2N+1)`, so it is NOT a gap member. The proof is a clean residue
argument at `t = 3·2⁻¹/(3N−1)`:

- `‖v t‖ < 3/Q ⟺ 3v ∈ {0,±2,±4} mod Q` (`Q = 3N−1`);
- using `3⁻¹ ≡ N` and `3(N−1) ≡ −2 (mod Q)`, the collisions are `{2N, N+1, N−1,
  2N−2}`, of which only `N−1` lies in the family's range — and `N−1` is exactly
  the removed element. So every speed clears; the binders `2` and `3(N−1)` sit at
  distance 3. Hence `M(F(N)) ≥ 3/(3N−1)`.

The **N=12 instance is formalized sorry-free and kernel-pure** in
`LRCEvenBranchWitness.lean` (`F12_margin : 3/35 ≤ margin F12 (19/35)`,
`F12_reach_above_gap : 2/25 < reach`), axioms `[propext, Classical.choice,
Quot.sound]`. So the canonical construction's failure at N=12 is now a *machine-
checked* fact: **by parity (12 even), not by 38 = 2·19 compositeness.**

## Corrective: the parity mechanism does NOT generalize to the whole class

I hoped THM-632's parity ("N even ⇒ an odd-denominator competitor `3/(3N−1)`
beats the mediant") would explain gap-emptiness for ALL bordered-AP families at
N=12, lifting the canonical result to the full crux. A sweep of the bordered-AP
candidate class (13,257 of 119,308 families tested) says **no**:

| outcome | count |
|---|---|
| IN gap `(1/13, 2/25)` | **0** (class is gap-empty — confirms opus-S118) |
| tight (`M ≤ 1/13`) | 1 (the AP `{1,…,12}`, at `1/13`) |
| loose (`M ≥ 2/25`) | 13,256 |

But the escapes are **dominated by the plateau `M = 1/12`** (denominator 12,
*even*), i.e. opus-S115's subfamily cap: a retained sub-AP pins `M = 1/12 > 2/25`.
The parity competitor `3/35` (denominator 35, odd) appears for essentially **only
the canonical family**. So:

> The general protector of the gap at N=12 is the **plateau / spectrum structure**
> (opus-S115 cap + the Farey rungs `1/13, 2/25` with nothing between), NOT an
> odd-denominator parity competitor. THM-632's parity is the mechanism for the
> *canonical mediant-aimer* specifically.

This is worth stating plainly because it kills a tempting overgeneralization:
"N even ⇒ gap empty by parity" is **false as a universal statement**. The gap
`(1/13, 2/25)` sits *below* the plateau `1/12`, so a would-be gap member must be
nearly as tight as the AP (near `1/13`) — the plateau doesn't even reach down to
it. What actually excludes such families is the spectrum gap itself (opus-S116),
of which the canonical mediant family is the one natural "aimer," and THM-632
shows it overshoots (to `3/35`, above even the plateau).

## Where this leaves the crux

- **Proved & formalized:** the canonical mediant family is not a gap member for
  even N (THM-632); N=12 machine-checked.
- **Confirmed (with structure):** the bordered-AP class at N=12 is gap-empty; the
  spectrum jumps `1/13 → 2/25` with the AP the unique tight family.
- **Still open (the real crux):** a structural proof that NO 12-speed family lands
  in `(1/13, 2/25)`. The plateau (opus-S115) covers the sub-AP-retaining bulk; the
  residual is families that break the sub-AP enough to approach `1/13` without
  landing on a rung — the genuine n-specific difficulty (opus-S114 "irreducible").

→ THM-632, HYP-4572, HYP-4506/opus-S118 (empty sweep), HYP-4476/opus-S115
(subfamily cap plateau), HYP-4486/opus-S116 (spectrum gap).
