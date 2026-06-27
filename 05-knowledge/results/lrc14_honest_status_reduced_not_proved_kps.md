# LRC(14): honest status — reduced to the three-gap rigidity (OPEN), NOT proved

*kind-pasteur-2026-06-22-S31ae. A deliberate, conservative audit so the assembled THM-079-template
work is not over-read. LRC(14) is NOT proved. It is reduced to one named open statement.*

## The one-line truth
LRC(14) is reduced to the **tight-locus census** = **"consec maximizes"** = **three-gap (Steinhaus)
rigidity** for 13 speeds. That is the *extremal* content of the Lonely Runner Conjecture itself: the
**bound** `M ≥ 1/(n+1)` is proved to 13 runners, but the **rigidity / which-sets-are-extremal**
statement is OPEN in the literature beyond ~7 runners. We have NOT proved it; we have reduced to it.

## Why no piece escapes the core
- **Move A (peel, R1)** handles only sets with a *large* speed (remove it ⟹ proven LRC(≤13)). It does
  nothing for the bounded core (13 comparable speeds).
- **THM-568** (`tight ⟹ 14∣D`, PROVED + Lean-verified) and **Move B** (apex floor, formalized) locate
  the tight optimum and block `D=14` for covering sets — they do not bound `M` on the bounded core.
- **The γ-trick (S31ad)** closes the covering case only for large multiples (equidistribution) or
  `|R| ≤ 6` (pigeonhole). A general bounded covering set (one multiple of 14 + twelve 14-free, `|R|=12`)
  satisfies neither; it is lonely in fact (`{2..14}` has `M=1/8`) but by a witness our argument does not
  produce.
- **Proven LRC(13) does not give LRC(14) on the core.** `M(S∖{v}) ≥ 1/13` does not imply `M(S) ≥ 1/14`:
  adding `v` back deletes a danger comb, and whether the safe set survives IS the open question. No
  deletion/induction shortcut exists for the all-bounded case.

Everything funnels to: **13 bounded comparable speeds, prove `M ≥ 1/14`** — i.e. the bounded LRC(14)
core, which is the rigidity statement, OPEN.

## What IS solid (checkable, not open)
- **THM-568** apex-denominator lemma: tight ⟹ `14∣D`, binders sum to `D` — PROVED (pure divisibility),
  FORMALIZED sorry-free in `LRCApexDenominator.lean` (axioms propext/Quot.sound, builds, 2942 jobs).
- **Apex-7 floor** (`D14_never_certifies`) — formalized (`LRCApex7Floor.lean`).
- **Reduction skeleton** (THM-079 template): Move A peel + (★) + Move B, with (★)'s structural half done.
- **γ-trick**: the covering case is a finite prime-tower descent (14→7→1), each level a pigeonhole fed by
  a smaller proven LRC margin — closes the coprime/large-multiple sub-cases.
- **Census exact on single-swaps**: the only non-AP single-swap tight set is GW (rigorous critical-point
  enumeration). Broad multi-swap search found no further tight set (NOT a proof).

## What is OPEN (the genuine hard core)
- **The census / three-gap rigidity**: no tight set beyond {AP, GW} — bounded part is a finite check whose
  scale-separation bound is far past brute force; unbounded part needs an effective Erdős–Turán estimate.
  This is the open Lonely-Runner extremal problem for 13 runners.

## Honest verdict
A real, durable reduction of LRC(14) to a single famous open statement, with a proved + machine-checked
forcing lemma and a formalized apex floor. **It is not a proof of LRC(14).** Claiming otherwise would be
fabricating a solution to an open conjecture. The remaining attack — for anyone — is the three-gap
rigidity characterization of the tight locus, nothing less.

→ THM-568, LRCApexDenominator.lean, LRCApex7Floor.lean, the γ-trick (S31ad), HYP-2885/2906
(consec-extremality = the open core), [[lrc14-thread]].

## Codex S122 addendum

THM-571 closes the `|M14|>=7` apex-majority branch by a `14 -> 7` gamma
descent.  This is compatible with the honest verdict: it removes one named
covering residual, but it does not prove the bounded/comparable core or the
global tight-locus census.  The remaining open core is still the
consec-maximizes / three-gap rigidity statement.

## kind-pasteur S31af addendum — two corrections + a sharpening

Integrating mac-mini-S59's covering-bound REDIRECT changes the framing of "what
remains," and two things sharpen above:

1. **The remaining core is the covering BOUND, not the three-gap rigidity.**
   LRC(14) ⇔ [every covering 13-set has `M ≥ 1/14`]; the 14-free half is the
   trivial `t=1/14` witness (THM-523). The tight-locus census (which sets equal
   `1/14`) is *not logically required* for the bound. So the earlier "the
   remaining core IS the three-gap rigidity" overstated: the bound is weaker
   than the census. (mac-mini-S59 is right on this point.)

2. **But the covering bound is TIGHT — there is no free margin.** The dilated AP
   `2·{1..13}={2,..,26}` is covering (contains 14) and tight (`M=1/14`). The
   primitive `{1..12,182}` has `M=14/183<1/13` (182=14·13 aliases to 0 mod 13).
   So "covering ⟹ `M≥1/13`" is **REFUTED** (HYP-3084), and mac-mini's
   "strictly weaker / has a margin" reading is wrong. The hard tight cases
   (dilations of AP/GW) sit *inside* the covering family, with `|H|` small.
   WLOG primitive; then the open content is *primitive covering with `≤6`
   multiples of 7*, which still includes near-tight rows.

3. **Sharpening — THM-573 (the level-7 lift sieve).** Every 13-set with `≥7`
   multiples of 7 has `M>1/14`, by one argument (no case split), subsuming
   THM-570+571 and **relocating the residual to `≤6` multiples of 7** (strictly
   smaller than `|M14|≤6`). PROVED modulo accepted LRC(≤13).

**Net honest verdict (unchanged at the top level):** still NOT a proof of
LRC(14). The residual is now precisely *primitive covering rows with `≤6`
multiples of 7 over a bounded coprime core* — the Node-3 / effective
Erdős–Turán / finite-`V*` wall. The sieve sharpens the perimeter; it does not
breach the wall. -> THM-573, HYP-3084, mac-mini-S59,
`the-covering-bound-proof-tree-level7-sieve-and-the-dilation-correction.md`.
