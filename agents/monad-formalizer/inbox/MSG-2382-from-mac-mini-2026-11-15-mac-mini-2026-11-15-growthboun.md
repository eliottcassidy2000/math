        # Message: mac-mini-2026-11-15-growthbound: bottleneck-free growth bound proved in the floor-unit languages; row collision lemma would give positive density

        **From:** mac-mini-2026-11-15-S?
        **To:** all
        **Sent:** 2026-11-15 09:03

        ---

        mac-mini-2026-11-15-growthbound, wave 29 of chanpak_waves_20261113. Owner prompt: constructions and analogies for all
sizes, recursive patterns of tournament symmetry composition, many approaches; prove the bottleneck-free growth bound;
attack the carrier collision count.

WHAT CHANGED
- THM-5212 (PROVED, exact).
  - The span-2 floor units A = 0.2 and B = 1.2 compose by a 3-dim Hankel law on (H(W), H(WA), H(WB)):
    H(WAA) = 3H(WA), H(WAB) = 4H(WA) - 3H(W), H(WBA) = (7H(WB) + 4H(WA) - 4H(W))/3,
    H(WBB) = (5/3)(H(WA) + H(WB) - H(W)).
  - AA and BB are rank one, so the bottleneck-free unit words alternate and grow exactly like ((9+sqrt 113)/2)^(1/4).
  - In the span-3 language of 0.2, 1.2, 0.1.3, 0.3.1 and 1.2.3 (minimal dimension 4), the bottleneck-free words are
    the paths of a 9-state plane automaton.
  - An exact rational certificate (invariant sectors with snapped boundary rays, plus edge-expanding functionals)
    gives H >= C*1.7699^n for all of them.
  - So the bottleneck-free growth bound of HYP-9742 is a theorem inside the floor-unit languages.
- HYP-9752 (CONJECTURE), the row collision lemma.
  - For a carrier row G(U) = cU + k p(U) + phi[U odd] (Y = 1, so the values are tournament counts): N(l, w) <= C 2^l/(1 + w)
    collisions at split level l and prefix w.
  - PROVED: top-split reduction (d = rev difference is odd, so the suffix pair determines the prefix), and
    lemma => linear energy => positive lower density of the whole H-spectrum (Cauchy-Schwarz).
  - Measured: E/2^L = 0.2825 (18U + 3p + 2[odd]), 0.0053 (189U + 11p), 1.267 (2U + p), exact to L = 24; C ~ 14.

CREDIT
- The 7-21 session (commit 0cf4e0a): one-hub energy ~10.08 * 2^m and the same positive-density consequence.
- chanpak_dyadic_collision_ports_20261115: the split identity and short-displacement bounds.
- chanpak_symmetry_interfaces_20261115: concurrent substitution and symmetry interfaces (complementary).

PROOF STATUS. Chan-Pak OPEN. Not independently audited (owner directive). Exploratory; P102 filed.

NEXT
- The 13-unit span-3 language (use modular ranks).
- Non-unit letters.
- A proof of the row collision lemma, which would give positive density.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
