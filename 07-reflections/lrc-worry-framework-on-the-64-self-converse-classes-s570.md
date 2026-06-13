---
source: opus-2026-06-03-S570 (remote-control)
status: SYNTHESIS — the n=14 worry-set = 64 self-converse round classes; every realizing transversal is lonely via an unblocked small pair (paired/anchored); AP unique tight
tags: [LRC, n14, worry-set, 64-classes, self-converse, round, transversal, paired-anchored, flip-set, S569, S576o]
---

# The worry framework on the 64 self-converse round classes (n=14)

**Prompt (user):** continue the worry framework with the 64 tournament classes.

The 64 = `2^((n-2)/2) = 2^6` **self-converse round tournament classes** that the
n=14 worry-set collapses to (oracle-S576o). This merges that even-ladder reduction
with my **paired/anchored split** (S569): every speed set realizing one of the 64
worry-classes is lonely via an **unblocked small pair**, and the AP is the unique
floor-tight one.

## 1. The 64 count, grounded

Self-converse round tournaments on `m` (odd) vertices number `2^((m-1)/2)`.
Verified by the circular-points generator (`lrc_selfconverse_round_count_s570.py`):

| m | 3 | 5 | 7 | 9 | … | 13 |
|---|---|---|---|---|---|---|
| self-converse round | **2** | **4** | **8** | **16** | … | **64** |

(`= 2^((m-1)/2)`, exact.) So at `n=14` (`m=13`) the worry-set lives inside **64**
self-converse round classes — a finite, tiny residual vs `A000568(13)≈4.85·10¹³`.

## 2. The realizing speed sets: census of the transversals

The 64 worry-classes are realized by **antipodal transversals mod `2n-1=27`** (one
residue from each pair `{a, 27-a}`; the AP = all-lower). Exhaustive census of all
**8191** gcd-1 transversals (`lrc_64_worry_classes_s570.py`):

- **LONELY: 8191/8191 — 0 counterexamples.**
- **Floor-tight (measure 0, the worry-set proper): exactly 1 — the AP `{1,…,13}`.**
  (V* and the other sporadics are *non-transversal*, the composite-27 cousins; S553.)
- **Unblocked small pair: ALL 8191 have one.** So my S569 mechanism (a pair `(a,b)`,
  `(a+b)/gcd ≤ 14`, with a pair-safe pinch clear of all runners) certifies **every**
  transversal lonely — uniformly across the whole 64-class container.

## 3. The merge (oracle S576o + my S569)

```text
worry-set  ⊆  64 self-converse round classes        (oracle S576o, even-ladder)
each realized by transversals, each of which has an  (this census)
   UNBLOCKED small pair  ⟹  lonely                   (my S569 paired/anchored split)
AP = the unique floor-tight class = the regular       (S557/S566/S568)
   rotational = the maximally-symmetric self-converse
```

So the paired/anchored mechanism is **uniform over all 64 worry-classes**: each
worry-config's private obstruction is its **straddle pair**, which is neither
**paired-away** (no multiple of the sum) nor **anchored** (sum `=n`, empty window),
hence open → witness. The 64 classes don't each need a bespoke argument — the *same*
unblocked-straddle mechanism handles them all.

## 4. The flip-set structure: AP is the unique tight, the 63 flips loosen

The 64 self-converse classes are the **flip-sets** `F` (which of the 6 free
antipodal pairs are flipped to their upper residue; `F=∅` is the AP). The census
shows:

> **`F=∅` (the AP) is the UNIQUE floor-tight class (`M=1/14`); all 63 non-empty
> flips LOOSEN it (`M>1/14`, positive measure).**

So flipping any antipodal pair strictly relaxes loneliness — the AP/regular
rotational sits alone at the floor, and every deviation along the 64-class lattice
moves *up*, away from the wall. This is the dual-Burnside fix-side (S565) seen as a
`2^6` Boolean lattice with the AP at its tight bottom.

## 5. What this leaves for a proof of LRC(14)

The reduction is now: **prove every speed set whose optimal-time tournament is one
of the 64 self-converse round classes is lonely.** The census handles the canonical
**transversal** realizations (bounded speeds) — all lonely, via unblocked small
pairs. The honest gap is the **structural lift**: a *general* speed set realizing
one of the 64 classes (unbounded speeds, the non-transversal composite-27 cousins
like V*) — show it too has an unblocked small pair. By S569 this is the single open
lemma *measure-0 ⟹ unblocked small pair*, now localised to the 64-class container.

## 6. Honest status

Grounded/verified: the 64 = `2^((m-1)/2)` count; the 8191-transversal census
(0 counterexamples, AP unique floor-tight, every transversal has an unblocked small
pair); the flip-lattice (AP unique tight, 63 flips loosen). New: the **merge** of
oracle's 64-class even-ladder reduction with my paired/anchored split — one uniform
mechanism over all 64 worry-classes — and the flip-lattice picture. Not a proof of
LRC(14): the remaining lemma is the structural lift (measure-0 ⟹ unblocked small
pair) for the non-transversal cousins.

**Artifacts:** `04-computation/lrc_64_worry_classes_s570.py` (+`.out`),
`lrc_selfconverse_round_count_s570.py` (+`.out`). Builds on oracle-S576o (even-
ladder/64), HYP-2095/S569 (paired-anchored split), S565 (dual Burnside), S568
(floor-tight target), S553 (transversals/sporadics). New: **HYP-2097**.
