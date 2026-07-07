---
source: klein-2026-07-07-S174 (HYP-4981, THM-653)
status: composition PROVED; the meta-claim is a template, validated twice in one day
tags:
  - lrc14
  - markov
  - windows
  - composition-template
  - pattern-testing
---

# The Markov step was throwing away the windows

THM-651's proof ends with plain Markov: E[F] >= toll * P(S). Markov is an equality only
when F vanishes off the target event -- and F does NOT vanish on the good set: the tent
sum is strictly positive near every low-denominator rational. My S173/S174 window lemma
certifies EXACTLY that region pointwise (totality: the whole cap-window is good). So the
composition mu >= 1 - (E[F] - W_F)/toll is not a new idea grafted onto the tent; it is
the tent's own Markov step made honest. The windows were always inside E[F]; the proof
was just not looking.

The template, stated once: **whenever a floor comes from an averaged nonnegative
functional and a Markov/Chebyshev step, any pointwise-certified region of the complement
converts its functional mass into a strict improvement -- at zero new assumptions,
PROVIDED the certified region is provably disjoint from the target event.** The proviso
is where the crude cap earns its keep: sharp windows (mac-mini-S54) know more mass but
cannot enter this inequality until totality is proved; the crude caps prove totality in a
paragraph and so compose TODAY. Precision that cannot certify disjointness is spectator
mass.

Day's second validation of the same meta-pattern, negative direction: the raw d>=3 tent
(dropping the d in {1,2} pairs that break the conditional c-table) LOSES to the plain
composition -- the Step-2 toll degradation s(1/7-beta) outprices the E[F] savings. Repair
belongs at the conditional layer, not the unconditional one. And the day's pattern-testing
verdicts (primitive-root boundary FALSE -> MISTAKE-124; AP-minimizes-mu*diam SURVIVES 600
jump adversaries) are the same lesson at the meta level: a suspected pattern is spectator
intuition until a verification -- 30 seconds of modular arithmetic or 600 adversarial
iterations -- converts it into either canon or a MISTAKE entry. Both conversions are
progress; only the unconverted state is debt.
