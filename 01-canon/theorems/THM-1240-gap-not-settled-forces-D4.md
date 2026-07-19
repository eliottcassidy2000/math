---
id: THM-1240
title: (1/14, 3/41) NOT SETTLED — BUT THE INTERVAL FORCES D ≥ 4, 4/55 IS THE CANONICAL TARGET, AND THE OBSTRUCTION TO SETTLING IS NAMED — M = D/s lies strictly inside iff 41D/3 < s < 14D, and for D = 1, 2, 3 that range contains NO INTEGER (13.67<s<14, 27.33<s<28, 41<s<42), so the interval **forces D ≥ 4** by exact arithmetic — one full step beyond boxeph-S123's D ≥ 3, which was for the wider Farey interval. Moreover |1·41 − 3·14| = 1, so 1/14 and 3/41 are FAREY NEIGHBOURS and 4/55 = (1+3)/(14+41) is the mediant: the unique fraction of least denominator inside, hence the canonical first target. ~12,400 families found nothing inside — including a residue-band construction at s = 55 built specifically to realise 4/55. But the candidate (D,s) list is INFINITE (every D ≥ 4 admits an s), so no finite check can settle emptiness: a proof requires a BOUND ON D for a realising family, which I do not have
status: the D ≥ 4 forcing and the Farey/mediant structure are exact arithmetic. The emptiness evidence is ~12,400 families over two structurally-motivated constructions — SEARCH, not proof. The question posed ("settle whether (1/14,3/41) is empty") is NOT settled, and this file says so and names precisely what a settlement needs
source: opus-2026-07-19-S397 (owner: settle whether (1/14, 3/41) is empty)
depends_on: [THM-1235 (which posed the question and corrected the gap edge to 3/41), THM-1230 (the 3/41 witness), boxeph-S123 (the D ≥ 3 stratification this strengthens), THM-1205 (the D-coordinate)]
scripts: 04-computation/settle_gap_opus_S397.py -> 05-knowledge/results/settle_gap_opus_S397.out
---

# THM-1240 — what the interval forces, and why it is not yet settled

## The interval forces D ≥ 4

M = D/s lies strictly inside (1/14, 3/41) exactly when **41D/3 < s < 14D**.
For small D that range holds no integer:

| D | range for s | admissible s |
|---|---|---|
| 1 | (13.67, 14) | **none** |
| 2 | (27.33, 28) | **none** |
| 3 | (41.00, 42) | **none** |
| 4 | (54.67, 56) | 55 → **4/55** |
| 5 | (68.33, 70) | 69 → 5/69 |
| 6 | (82.00, 84) | 83 → 6/83 |
| 7 | (95.67, 98) | 96, 97 |

So the attained interval **forces D ≥ 4**, by arithmetic alone. That is one
full step beyond boxeph-S123's D ≥ 3, which applied to the wider Farey
interval (1/14, 2/27); narrowing the interval to its attained edge sharpens
the determinant floor.

## 4/55 is the canonical target

|1·41 − 3·14| = |41 − 42| = 1, so **1/14 and 3/41 are Farey neighbours**.
Hence every fraction strictly inside has denominator ≥ 14 + 41 = 55, and the
unique one achieving 55 is the mediant

> **4/55 = (1+3)/(14+41)**

So 4/55 is not merely the first candidate by size — it is the unique
least-denominator fraction in the interval, and any realisation of the
interval at all makes 4/55 the natural first place to look.

## The evidence, and its limit

Nothing was found inside, over two structurally-motivated constructions:

- a **residue-band construction at s = 55** — pick t* = p/55, force all
  thirteen speeds into the band [4, 51] mod 55 so that M ≥ 4/55 holds *at
  that point by design*, then test whether any other point beats it;
- the **{1,…,11, x, y}** two-free-slot shapes, which is where the
  near-extremal families live.

Combined with THM-1235's 1552-family scan, roughly **12,400 families** have
now been tested with **zero** hits inside (1/14, 3/41), including zero
realisations of the mediant.

**This is not a settlement.** The candidate (D,s) list is **infinite** —
every D ≥ 4 admits an integer s ≈ 14D — so no finite enumeration closes the
question. In the band construction the residues are never the obstruction
(at D = 4, s = 55 the band holds 48 of 55 residues, 87%); what defeats every
attempt is that some *other* point beats the intended maximiser. That is a
global condition, and it is exactly what a proof would have to control.

## What settling requires

> **Bound D for a family whose maximiser lies in (1/14, 3/41).**

With such a bound the candidate list becomes finite and each (D,s) is a
finite residue problem. Without it, emptiness cannot be decided by search at
all. I record this as the precise obstruction rather than as more evidence:
the honest state is that the interval is *very likely* empty and *not
proved* so, and the missing ingredient has a name.
