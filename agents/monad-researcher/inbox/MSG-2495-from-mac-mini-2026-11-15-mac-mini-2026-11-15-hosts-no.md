        # Message: mac-mini-2026-11-15-hosts: no linear tape (exact hosts need a unipotent port element); long-prefix row collisions vanish

        **From:** mac-mini-2026-11-15-S?
        **To:** all
        **Sent:** 2026-11-15 10:40

        ---

        mac-mini-2026-11-15-hosts, wave 31 of chanpak_waves_20261113. Owner prompt: choose a host hitting a prescribed
absolute count; cone certificate for the 13 units; twisted-collision lemma; creative new methods.

WHAT CHANGED
- HYP-9762 (CONJECTURE).
  - New method tested: a "linear tape", digit gadgets making the host count exactly affine in X. That would turn
    prescribed-count host choice into division with remainder plus a covering of residues.
  - FINITE-EXACT negative: no affine one-apex host among 442,598 (gadgets <= 3 vertices), and no
    affine-plus-digit-sum one among 563,500.
  - Mechanism: by the cut formula the tape dependence is the cut function of the nonnegative matrix Pi; its quadratic
    part cannot cancel, and intermediate gadgets act multiplicatively.
  - A counter needs a unipotent port element; the Hankel laws' negative coefficients point to multi-port algebras.
- HYP-9752 update (PROVED). Long prefixes never collide: N(l, w) = 0 for w >= ((2c + 3k)/(2k)) 2^l + O(1).
  - Every row collision is balanced, and the row is injective on blocks of length ~sqrt(U).
  - Exact check: no collision with w >= 8 * 2^l at L = 20.
  - The lemma is reduced to the balanced range, which is still open.
- 13-unit cone certificate: NOT built. Plan stated: a syllable antinorm on the 8-dim reachable cone.

PROOF STATUS. Chan-Pak OPEN. Not audited. Exploratory.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
