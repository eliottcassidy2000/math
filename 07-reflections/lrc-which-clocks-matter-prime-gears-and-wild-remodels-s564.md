---
source: opus-2026-06-02-S564 (remote-control)
status: SPECULATION + grounded core — worry/ignore separation, the prime-gear clock, and a menu of drastic remodels for "does the orbit hit B even once"
tags: [LRC, n14, clocks, worry-ignore, prime-gears, CRT, covering-systems, remodel, speculation, S557, S559, S563]
---

# Which clocks matter, what to ignore, and wild ways to ask "does the orbit hit B once"

**Prompt (user):** which clocks matter? separate the cases we must worry about
from those we can ignore. Find a more elegant determination of "does this orbit
hit the safe box even once," remodelling the scenario if needed. Be wild.

## A. The grounded core: separate WORRY from IGNORE

The orbit `γ(t)=(v_i t)` either spends a positive-length *interval* in the box `B`
or only *grazes* it. That single distinction is the whole game:

- **IGNORE — positive safe-measure.** The orbit sits in `B` for an interval, so it
  hits trivially. This is **almost everything**: every incommensurate set (Weyl-
  dense, S563), every low-resonance set, and — measured — **300/302** random n=14
  sets. A handful of random clocks `a/M` certify these (each clock is safe with
  prob ≈ `(1-2/n)^{13} ≈ 0.135`, so ~30 random clocks certify >99%). *Never compute
  anything hard for these.*
- **WORRY — measure-zero.** The orbit only touches `B` at isolated instants (tight)
  or misses (the hypothetical counterexample). Measured: **2/302** (exactly the AP
  and V*). These are the **resonance-maximal** sets (S563) — vanishingly rare and
  highly structured.

> **You only ever have to worry about the measure-zero, resonance-maximal sets.
> Everything else is certified by a few random clocks.**

## B. Which clocks matter

- **Complete clock family:** the `O(k²)` **pair-sum clocks** `t=m/(v_a+v_b)`. The
  optimal witness is always one of these (S557 pinch; S562 pinch sieve is complete).
  Not `{1,…,M}` (incomplete at any finite `M`, S551), not the `c=(k+1)` lift.
- **For the WORRY set, ONE clock:** the **n-clock** `t=j/n`. Tight ⟹ the binding
  pair sums to `n` ⟹ witnessed at `j/n` (S553/S556). So the entire hard question is
  the n-clock on the resonance-maximal family.

## C. The drastic remodel that I think is *right*: TIME IS A STACK OF PRIME GEARS

Factor the n-clock by CRT: `ℤ/14 = ℤ/2 × ℤ/7`. Then **time at resolution `1/14` is
two gears turning together** — a parity gear (mod 2) and a 7-gear (mod 7). At gear
position `j`:

> a runner `v` is **at the observer iff it aligns on EVERY gear** — i.e. `14 | vj`
> ⇔ (`2 | vj`) **and** (`7 | vj`).

So a runner is dangerous only at a **full alignment**, the coincidence of all prime
gears. The AP `{1,…,13}` has **no runner that is a multiple of 14**, so **no runner
ever fully aligns — every gear setting is CLEAR** (verified: AP clear at all
`j∈{1,3,5,9,11,13}`). Danger can only come from runners that are multiples of `n`,
and whether they can be dodged is the **gear-coupling** problem (the apex/parity
interplay, S559). This reframes LRC@14 as:

> **Can the parity gear and the 7-gear always be co-set so that no runner is fully
> aligned?** The conjecture says yes; the worry is exactly the resonant configs
> where the two gears are forced to conspire.

This is elegant because it makes the prime structure of `n` *the model of time*: the
runners live on the gear-coincidences, and the difficulty is a small CRT
constraint, not a continuum.

## D. Wilder remodels (speculation menu)

1. **LRC ⇄ COVERING SYSTEMS (Erdős).** A counterexample = the danger arcs *exactly
   cover* the circle. At the n-clock these are residue classes; "the arcs cover" =
   a **covering system of congruences**. The divisibility sieve (THM-369 / Rosenfeld)
   *is* a covering system. So LRC@14 ⇄ "no AP-arc covering system tiles the circle
   with these widths," and the sieve's unbounded modulus (S551) is the Erdős
   minimum-modulus phenomenon. **This is the most promising bridge** — covering-
   system technology (Hough's bound on minimum modulus, etc.) may import directly.

2. **THE TWO DUAL RHYTHMS (beats).** Pair-**differences** `|v_i−v_j|` are the
   *collision beats* (when two runners coincide, the tournament edge flips);
   pair-**sums** `v_i+v_j` are the *pinch clocks* (when loneliness is optimal). Time
   has two interleaved frequency combs; loneliness is a sum-beat threading between
   difference-beats. Categorise speed sets by the alignment of these two combs.

3. **FOLD 13 RUNNERS INTO 1.** Where speeds are pairwise coprime, CRT folds them
   into a single super-runner on a circle of circumference `∏v_i` facing a product
   forbidden set. For general speeds, fold the 2-adic part (even→half, S554) and the
   coprime parts separately — a recursive collapse toward a *single* clock.

4. **THE ORBIT AS A CYCLIC CODE.** The closed orbit visits a cyclic sequence of
   cells (tournaments); this sequence is a length-`?` cyclic word determined by the
   speeds (the H-spectrum fingerprint, S26). Loneliness = a designated symbol occurs.
   LRC = "this speed-indexed cyclic code is never `lonely`-free." Coding-theory
   weight bounds?

5. **REVERSE THE ARROW — the box radiates.** Put sources in `B`; the speed-
   directions it fails to illuminate are the counterexamples. By Combinatorial
   Nullstellensatz the non-illuminated integer directions are the zero set of an
   explicit polynomial; show it has no integer points (the polynomial-method, S559,
   is exactly this for `k+1` prime — the `2q` apex is its zero-divisor).

6. **REPLACE CONTINUOUS TIME BY THE FAREY DISSECTION.** The critical times are
   Farey-like fractions (`j/v_i`, pair-sums); the three-gap theorem governs the
   spacing. Walk the Stern–Brocot tree of speed ratios; loneliness is a node
   condition. (Repo: Zeckendorf/Stern-Brocot threads, S460.)

## E. The efficient determination, assembled

```
hits_box(V):
   # IGNORE tier: certify the easy 99.x%
   if any of ~30 random clocks a/M is safe:        return HITS        # positive measure
   # WORRY tier: only the resonance-maximal remain
   if the n-gear t=j/n has a clear setting:         return HITS (tight)
   else:                                            return CANDIDATE COUNTEREXAMPLE
```

The first line disposes of all incommensurate and low-resonance sets in O(1)
clocks; the second is the single CRT-gear check on the rare resonant residual. The
*only* object a proof must control is the n-gear on the resonance-maximal family —
and remodel C/D-1 say that is a **covering-system** question about coupled prime
gears.

## F. Honest status

The worry/ignore split, the pair-sum/n-clock identification, and the prime-gear
model are grounded (proved or measured). The covering-system bridge (D-1) and the
other remodels are *speculation* — but D-1 in particular connects LRC@14 to a
mature body of technique (Erdős covering systems, Hough's theorem) that the repo
has not yet mined, and it is the same object every thread converged on: the
resonance-maximal closed orbit, now seen as a covering system of prime gears.

**Artifacts:** `04-computation/lrc_which_clocks_gears_s564.py` (+`.out`).
Builds on S557 (pinch clocks), S559 (CRT/apex), S562 (multi-sieve), S563 (orbit/
resonance), THM-369. New: **HYP-2081**.
