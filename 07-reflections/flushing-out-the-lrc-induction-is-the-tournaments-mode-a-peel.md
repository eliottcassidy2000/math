# Flushing Out the LRC: the Induction IS the Tournament's Mode-A Peel

*mac-mini-2026-06-22-S47. The owner asked me to immerse in the repo's tournament-induction attempts,
then flush out the LRC. They are the same induction. The LRC(14) reduction (R1/R2/R3) is the tournament
Mode-A peel; the bounded core is the irreducible atom; the per-path failure and the H=21 success are
the two transferable lessons. Builds on kps S31w (reductions), codex (scale-separated lemma), HYP-2900,
THM-079 (H=21 proof), the per-path MISTAKE, [[everything-is-the-triangle]] (Mode A/B).*

## The LRC induction = the tournament Mode-A peel

The three LRC(14) reductions (kps S31w + my HYP-2900):
- **R1 — remove a large speed (n → n−1)** = the tournament **Mode A** (single-vertex deletion; the
  hypotenuse / H=1+2^d side). A speed `v ≫ rest` equidistributes (comb-teeth), so
  `meas(safe(S)) ≥ (6/7−ε)·meas(safe(S∖{v}))`; nonempty by LRC(n−1). The unbounded case peels off.
- **R2 — omit a prime** = the resonance witness: omit a multiple of `p≤13` ⟹ `t=1/p` gives `M≥1/13`.
- **R3 — dilation** = normalize to primitive `gcd=1`.

VERIFIED: R1/R2/R3 reduce **every** 13-set to the **irreducible bounded core = {consec, GW sporadics}**
(consec sits exactly at the R1 boundary: `max=13 = Σrest/6`; it is prime-covering so R2 can't fire; it
is primitive so R3 is trivial). This is the LRC analog of THM-079's reduction of `H=21` to a **single
strong atom** (`H` multiplicative ⟹ WLOG strong).

## Lesson 1 — the per-path failure: the bounded core is the irreducible obstruction

The per-path identity (a naive per-element induction for Claim A) **works at n≤5 but fails at n=6**
(30% of triples), exactly when the first **disjoint 3-cycles** (`β_3`) appear. The lesson: a per-element
induction collapses the easy strata but cannot dissolve the first homological obstruction. The LRC's
obstruction is structurally the same: R1/R2/R3 peel everything *except* the bounded core, which is the
genuine `n`-runner content — `{consec, GW}` is the LRC's "disjoint 3-cycles." The induction cannot peel
it; it must be bounded directly. (So the irreducibility of the core is not a gap in the method — it is
the method correctly locating where the real content lives.)

## Lesson 2 — the H=21 success: reduce to one atom, then bound it

THM-079 proved `H=21` impossible by: (A) multiplicative reduction to a **single strong component**;
(B/C) an **inductive size bound** (`m≥9`) + **Moon pancyclicity** ⟹ contradiction. The LRC mirrors A
exactly (R1/R2/R3 ⟹ a single bounded-core atom). What remains is the LRC analog of B/C: **bound the
atom** — show the bounded core has `M ≥ 1/14`. The core's value is `M=1/14` (consec at `t=1/14`; GW at
`t=1/n`), VERIFIED. The rigorous bound is "consec maximizes coverage" via three-gap rigidity (Node 2) —
the LRC's Moon-pancyclicity step. THM-079 is the template: the hard part is always the single-atom bound,
reached after a clean multiplicative/scale reduction.

## The winding-tournament link

At the loneliest moment of consec, `t=1/14`, the runners sit on the 14-grid `{k/14}` — evenly spaced,
the **most-ordered** configuration; the winding tournament there is the circulant `i→j iff (i−j) mod 14
∈ {1..6}` (with the antipodal tie at difference 7 = the apex-7 boundary). So the LRC extremal is the
maximally-symmetric winding tournament, and the loneliness value `1/14` is the grid spacing. The
bounded-core extremal is a *tournament* extremal — closing the loop the owner pointed at.

## The assembled induction (honest)

LRC(14) ⟸ **[R1/R2/R3 induction]** (reduces to the bounded core; rigorous modulo the scale-separation
threshold + effective `ε`, codex's scale-separated lemma) ⟸ **[bounded core = {consec, GW}, finite, all
`M=1/14`]** (modulo tight-locus finiteness — consec+GW only — kps's THM-560 + GW census; and the
"consec maximizes" three-gap bound). Both remaining cruxes are the team's live work. The induction is
mapped and the obstruction is located at a single finite atom — the cleanest the proof has been — but
the atom's bound (the H=21-analog Moon step) is not yet closed, so LRC(14) is **not** finished.

Related: kps S31w (R1/R2/R3), HYP-2900 (Node 3), THM-079 (H=21 template), the per-path MISTAKE,
[[everything-is-the-triangle]] (Mode A/B), HYP-2605 (winding tournament).
