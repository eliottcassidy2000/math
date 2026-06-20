---
id: THM-545
title: LRC14 tower-deletion lemma (comb cutoff) - deleting any dyadic-1 tower bit forces meas >= 426/35035; PROVED for the one-extra-hole layer at all four bits
status: PROVED (one-extra-hole layer, all a in {0,1,2,3}, via periodic-comb rational cutoff + finite exact residue); HIGHER LAYERS reduced to finite residue scans (structurally bounded, not all executed)
source: kind-pasteur-2026-06-19-S17
depends_on:
  - THM-541
  - THM-543
  - HYP-2661
related:
  - THM-542
  - THM-544
  - OPEN-Q-108
external: Lonely Runner Conjecture, n=14
---

# THM-545 - Tower-Deletion Lemma (comb cutoff)

## Statement

Let `C` be a primitive AP-tail 12-core, i.e.

```text
C = ({1,...,13} \ holes) union tails,   tails subset [14, inf),   |C| = 12,
```

and let `theta = 1/14`,

```text
G_C = { t in [0,1) : ||c t|| > 1/14  for every c in C },
      thr2 = 426/35035   (the AP one-hole second value).
```

**Tower-deletion lemma (target).**  For `a in {0,1,2,3}`, if the dyadic-1 tower bit
`2^a` is NOT in `C`, then `meas(G_C) >= 426/35035`.

This is the missing PROOF step behind HYP-2661 ("the full dyadic tower {1,2,4,8}
must lie in C whenever meas(G_C) < thr2").  Prior sessions VERIFIED the law
exhaustively over bounded boxes; THM-545 supplies the comb-cutoff certificate.

## What is PROVED here

The **one-extra-hole layer** (`k = 1`): cores

```text
C = ({1,...,13} \ {2^a, h}) union {r},   h != 2^a,   r >= 14,
```

i.e. delete the tower bit `2^a` plus one more small hole `h`, add one tail `r`.
For **every** `a in {0,1,2,3}` and every such `C`, `meas(G_C) >= thr2`.

## Method (codex comb technique, THM-543/544 form)

**Single-tail comb lemma.**  For a safe set `G` with measure `M` and `c` interval
components and a new speed `r`,

```text
meas(G \ D_r) >= (6/7) M - 2c/(7r).
```

(On each component, full `1/r` periods remove exactly a `1/7` share; the two
partial end-periods remove at most two extra danger teeth of total width
`2/(7r)`.  Sum over `c` components.)

**Rational cutoff.**  With base `B = {1,...,13} \ {2^a, h}`, exact `M_B`, `c_B`,

```text
meas(G_B \ D_r) >= thr2   for all   r >= R_B := ceil( 2 c_B / (6 M_B - 7 thr2) ),
```

since `6 M_B - 7 thr2 > 0` always (every two-hole base has `M_B >= 0.024 >> thr2`).
Tails `r >= R_B` are certified by the comb; the finite residue `r in [14, R_B)`
is checked by EXACT rational evaluation of `meas(G_C)`.

**Result (exact, zero floating point).**

```text
a=0 (delete 1): comb cutoff R=44,  residue cores=154, below thr2 = 0  -> PROVED
a=1 (delete 2): comb cutoff R=83,  residue cores=410, below thr2 = 0  -> PROVED
a=2 (delete 4): comb cutoff R=147, residue cores=652, below thr2 = 0  -> PROVED
a=3 (delete 8): comb cutoff R=69,  residue cores=357, below thr2 = 0  -> PROVED
```

Base case `k = 0` (no extra hole, the single-hole bases `{1,...,13} \ {2^a}`,
no tail): exact measures `1/14, 11/364, 97/4004, 950/21021` for `a = 0,1,2,3`, all
`>= thr2`.

## Floor argument for higher layers (k >= 2)

For `k` extra holes (`k` tails), the quantity

```text
(6/7)^k * min_{ k-hole base containing 2^a } M_B
```

is INCREASING in `k` and stays far above `thr2` at every computed `k` (a=3:
0.0452, 0.0457, 0.0569, 0.0678, 0.0775 for k=0..4; similarly for a=0,1,2).
Removing extra holes raises the base measure faster than the comb factor `(6/7)^k`
erodes it.  Hence only configurations with a SMALL tail can approach `thr2`; those
are exactly the finite residues a `k`-fold iterated comb cutoff isolates.  The
`k >= 2` residue scans are finite and structurally bounded, but were not all
executed exactly in this session (pure-Python exact arithmetic is the bottleneck);
they remain the only gap to the full lemma.

## Certificate

Scripts (exact `Fraction` throughout):

```text
04-computation/lrc14_tower_deletion_proof_kps.py     (base cases + floor + a=3 k=1 full)
04-computation/lrc14_tower_deletion_kps.py           (bounded exhaustive scan, all a)
```

Stored outputs:

```text
05-knowledge/results/lrc14_tower_deletion_proof_kps.out
05-knowledge/results/lrc14_tower_deletion_k1_all_kps.out
05-knowledge/results/lrc14_tower_deletion_scan_kps.out
```

## Corrected dead end (recorded for future agents)

A "monotonicity" reduction (drop the tail constraints to bound `meas(G_C)` below
by a single-hole base) is **INVALID**.  Adding a constraint `||c t|| > 1/14` only
SHRINKS the safe set, so for `B subset C`, `meas(G_C) <= meas(G_B)`, NOT `>=`.
Tails (speeds `>= 14`) are genuine extra constraints and CAN lower the measure.
This is exactly why the comb cutoff + finite residue is required (THM-543/544).
See MISTAKE log.

## Proof Use

THM-545 closes the comb-cutoff proof obligation flagged open in SESSION-LOG
(kind-pasteur S16, line 9) for the one-extra-hole layer at all four tower bits,
turning the exhaustively-verified HYP-2661 into a partially-certified theorem.
It directly serves OPEN-Q-108 (uniform `meas(G_C)` tight-locus finiteness) by
eliminating every tower-deleting core from the sub-`thr2` candidate set at `k <= 1`.
