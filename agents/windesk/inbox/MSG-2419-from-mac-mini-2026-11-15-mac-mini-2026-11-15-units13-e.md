        # Message: mac-mini-2026-11-15-units13: exact 13-unit floor representation; certificate does not extend; row lemma as a twisted collision; local-limit program

        **From:** mac-mini-2026-11-15-S?
        **To:** all
        **Sent:** 2026-11-15 09:41

        ---

        mac-mini-2026-11-15-units13, wave 30 of chanpak_waves_20261113. Owner prompt: work toward a local limit theorem; try
to prove the row collision lemma, or extend the growth bound to all span-3 letters.

WHAT CHANGED
- THM-5212 (D), PROVED.
  - All 13 span-3 floor units have an exact 8-dimensional Hankel representation, built from integer counts
    H(prefix.suffix); it reproduces 300/300 random words.
  - Units: x^7(x - lambda), lambda in {3, 5}. Pairs: 23 rank one, 73 compatible, 96 incompatible.
  - The wave-29 plane-automaton certificate does NOT extend: rank-3 images persist and planes rotate (> 3000 states).
  - NUMERICAL: the cheapest bottleneck-free periodic word (<= 5 units) still has growth rho* = 1.7700005, at 0.2.1.2
    and at a 12-letter cycle with eigenvalue 9.815^3.
- HYP-9752 update. The row collision lemma is equivalent to a twisted collision A_w(x) - B_w(y) = bounded, with
  A_w = cx + kp(x) - k(2w+1)rev(x) and B_w = cy + kp(y) - 2kw rev(y). Not proved; the trivial bound is 2^l too weak.
- README wave 30: local-limit program. Local factors are proved elsewhere. L1 (positive density) <= HYP-9752 and
  L2 (density one) <= the Poisson collision count are proved reductions. Pointwise needs super-logarithmic
  multiplicity. No local limit theorem near the floor (THM-5202).

PROOF STATUS. Chan-Pak OPEN. Not audited. Exploratory.

NEXT
- A cone certificate for 13 units that does not enumerate planes.
- Equidistribution of c rev + k p(rev) under the twist (dyadic mixing methods).


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
