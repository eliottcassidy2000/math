# Cutting sequences and the aggregate extremiser (LRC(14), Angle B)

*mac-mini-2026-06-18-S7 — a reflection on what the symbolic-dynamics view of the seven-sector
cover taught us, and why three natural proofs of AP-extremality all die.*

## The coordinate that was hiding in plain sight

The LRC(14) crux had been written as `meas(S7(E)) ≤ cap_k`, with `S7` the set of `x` where the
phases `frac(e_i x)` hit all seven sectors `[j/7,(j+1)/7)`. The sector index `σ_e(x)=⌊7ex⌋ mod 7`
is *literally* the cutting sequence of a slope-`7e` line through the integer grid — a billiard
word. For the consecutive cluster `{0,…,k−1}` (the conjectured worst case), a single substitution
`θ=7x` collapses the whole thing: `σ_e = ⌊eθ⌋ mod 7`, and the residues hit are the **partial sums
of a mechanical word**. The seven-sector cover of the AP is a *Sturmian cover-time*. That is a
genuinely clean re-coordinatization — the same object the three-distance theorem governs. It made
two facts immediate that had only been observed numerically: the cover *vanishes* for `k≤6` (you
visit at most `k<7` partial sums) and it is *scale-invariant but not translation-invariant* (the
`e=0` runner pins residue `0` at the origin).

## The lesson: the extremiser is aggregate, not structural

I went in expecting the AP to win *structurally* — per inclusion–exclusion block, or by a
monotone "spreading loses coverage" argument, or by an AP-orbit symmetry. **All three are false.**
- Per-block: the AP is neither the max nor the min of the individual `meas(B_M)` terms. About half
  the signed differences against the AP are *negative*. The AP only wins after all 64 signed terms
  cancel into the aggregate.
- Monotone: pushing the top offset out by one both helps and hurts; `{0..6,10}` beats `{0..6,9}`.
- Symmetric: `meas(S7)` has no translation symmetry to exploit; only the scaling orbit is free.

This is the real content of the session, and it is a *negative* result that sharpens the target.
The reason THM-534's moment-LP dual and THM-535's cap-split had to work so hard — a degree-4
Bonferroni majorant for k=8, three genuinely-tight rational rows — is that there is **no cheaper
certificate**. AP-extremality of a covering measure is an irreducibly aggregate rearrangement
inequality, the kind that resists term-by-term or symmetry proofs and yields only to a global
convexity/moment argument. The cutting-sequence view doesn't make it easier; it makes *why it is
hard* legible.

## A pattern across the whole LRC(14) effort

This rhymes with what the project keeps finding on this problem. Every time a clean "the AP is
extremal" statement is dispatched, the honest version turns out to be: *true on the dangerous rows,
false (but harmless) on the slack rows, and never reducible to a local move.* HYP-2604 (meas(S7)),
THM-534 (the moment functional `L_y`), HYP-2602 (`μ_{1/7}`): same shape every time. The AP is the
most *aligned* cutting sequence — its partial sums sweep the residues most thoroughly — but
"alignment" is a global property of the whole word, not of any factor or any single flip. The
three-distance theorem tells you the *gaps* of one mechanical word; the cover needs the *joint*
word over all of `E`, and the joint word's surjectivity is exactly where the local structure
dissolves into an aggregate.

## What I'd tell the next session

The two rigorous lemmas (k≤6 vanishing, pointwise subset domination) are worth keeping — the
domination lemma in particular certifies a real band of bounded-span shapes for free, and it
composes with THM-535's cap-split. But the live gap is unchanged and it is now *characterized*:
the only route left for the three tight rows (k=8,9,10) is a global moment/convexity argument
(THM-534's `L_y(E) ≤ cap_k`), because every structural shortcut has been explicitly closed. Spend
no more time looking for a per-block, monotone, or symmetry proof of AP-extremality of `meas(S7)` —
those doors are now known to be walls. The Sturmian partial-sum cover-time is the right language to
*state* the moment inequality cleanly; whether it helps *prove* it (via a transfer operator on the
mechanical word, say) is the open question worth the next session's energy.

→ THM-536, THM-534, THM-535, HYP-2604, HYP-2602, THM-531, OPEN-Q-108.
