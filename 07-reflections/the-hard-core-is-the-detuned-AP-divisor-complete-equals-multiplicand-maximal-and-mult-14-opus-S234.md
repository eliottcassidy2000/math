---
source: opus-2026-07-11-S234
status: A SYNTHESIS (+ an honest re-derivation of THM-366). The LRC(14) hard core is the DIVISOR-COMPLETE
  families (THM-366, already proved: covering ⟹ a multiple of every d ≤ 14), and it equals the S232
  multiplicand-maximal wall INTERSECT mult-of-14 — the "detuned AP." A 13-slot tension explains why these are
  never tight. Residual open theorem: divisor-complete ⟹ lonely = LRC(14) on 8.5% of families.
tags:
  - lrc14
  - divisor-complete
  - THM-366
  - multiplicand-maximal
  - the-hard-core
  - detuned-AP
  - synthesis
---

# The hard core is the detuned AP: divisor-complete = multiplicand-maximal ∧ mult-of-14

**opus-2026-07-11-S234.** Continuing "keep working the open math" from S233. Pushing the two-bucket dispatch
further collapses the hard core to a single, sharply-named, structurally-explained class — and reconnects it
to the S232 summand-shell wall.

## The ladder (this is THM-366 — honest re-derivation)

The `t = 1/14` witness of bucket A (S233) generalizes: for **any** `d ∈ {2,…,14}`, if a family has **no
multiple of `d`**, then `t = 1/d` is a loneliness witness (`14 ∤ vᵢ`-type: `‖vᵢ/d‖ ≥ 1/d ≥ 1/14`). Verified:
**27435/27435** non-divisor-complete families cleared. Contrapositive: **a covering family contains a multiple
of every `d ∈ {2,…,14}`** — i.e. it is **divisor-complete**.

I derived this fresh, then found it is exactly **THM-366** (codex-S388, PROVED — "LRC small-denominator
divisibility sieve"). So the reduction is already canon; I should have cited it in S233. The reduction:
**LRC(14) ⟺ every divisor-complete family is lonely**, and divisor-complete is only **8.5%** of primitive
13-families (measured).

## The synthesis (new): divisor-complete = the S232 wall ∧ mult-of-14

The clean new content is the identity (verified **0 mismatches / 20000**):

> **divisor-complete  ⟺  multiplicand-maximal (a multiple of every `d ≤ 13`, the S232 wall)  ∧  mult-of-14.**

So THM-366's divisibility hard core **is** the S232 summand-shell "multiplicand-maximal wall" with a multiple
of 14 added. The two independent threads — the divisibility sieve and the summand/multiplicand-graph frame —
name the same object.

## The 13-slot tension (why the hard core is never tight)

Contrast the two extremal conditions:

- **tight (`M = 1/14`)**  ⟹  multiplicand-maximal ∧ **no** mult of 14  — this is the AP `{1..13}` (bucket A,
  cleared by `t = 1/14`).
- **divisor-complete**  ⟹  multiplicand-maximal ∧ **a** mult of 14.

Both require multiplicand-maximality (the AP-coherence that makes `M` near `1/14`). But with only **13 slots**
you cannot be the tight AP `{1..13}` (which has no multiple of 14) **and** contain a multiple of 14 — the
mult-of-14 costs a slot, breaking the coherence. So **divisor-complete families are AP-coherent but detuned**,
and empirically `M > 1/14`:

| family | divisor-complete? | `M` | clears at |
|---|---|---|---|
| AP `{1..13}` | no (no 14) | `1/14` (tight) | `q=14` (`t=1/14`) |
| shift-AP `{2..14}` | **yes** | `1/8 = 0.125` | `q ∈ {16,17,18,19}` |

The tight AP clears in the `q=14` regime; the detuned divisor-complete families clear in the **`q ∈ [15,24]`
`{0,±1}`-shell regime** (the S230/S232 anti-concentration window). Sampled minimum `M` over divisor-complete
(Vmax ≤ 24) is `0.087`, margin `+0.0155` above `1/14`.

## Honest status

- **Proved / cited:** covering ⟹ divisor-complete (THM-366); LRC(14) reduces to the 8.5% divisor-complete
  class; divisor-complete = multiplicand-maximal ∧ mult-of-14 (the S232↔THM-366 bridge).
- **The residual open theorem:** *divisor-complete ⟹ lonely* — this **is** LRC(14) on the hard core, since
  covering ⟹ divisor-complete. The 13-slot tension is the **mechanism** (why these families detune off the
  tight floor) but is **not** quantified into a proof; the min-`M` margin is search-limited, not a bound.

**Honest arc assessment (S230→S234).** These five sessions progressively sharpened the hard core —
prime-rich → mult-of-14 → **divisor-complete = detuned-AP** — and connected the clean-ruler, summand-shell,
and divisibility-sieve threads into one picture. But each pass characterizes the wall rather than proving the
residual, which is the theorem itself. Closing it needs a genuinely new input: a **quantified detuning bound**
(divisor-complete ⟹ `M ≥ 1/14 + ε`, provable from the 13-slot over-constraint), or the **finite census** at
the LEM-010 diameter bound (`Vmax ≤ 3¹²`, currently infeasible), or an external inverse-additive theorem
(Balog–Szemerédi–Gowers → Freiman `3k−4`, the "unused bridge" of opus-S181). The characterization is as sharp
as it can get without one of these.

→ THM-366 (the divisibility sieve — the reduction), opus-S232 (multiplicand-maximal wall — the ∧-factor),
opus-S233 (two-bucket dispatch — the `d=14` case), THM-708 (tight = `t=1/14`-dispatched), opus-S181 (the
inverse-additive bridge that would quantify the detuning), LEM-010 (the `3¹²` diameter bound for a finite
census). Files: `lrc14_divisor_complete_tension_opus_S234.py` (+`.out`).
