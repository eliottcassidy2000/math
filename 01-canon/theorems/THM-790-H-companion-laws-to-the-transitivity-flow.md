---
id: THM-790
title: The H-COMPANION to the transitivity flow — the same blue/black complement-line flow measured on the project's Rédei invariant H (Hamiltonian-path count) instead of the degree-derived C3/E4 axes of THM-785/787. (a) VERIFIED LAW (even n): blue lines preserve H mod 4 (ΔH ≡ 0 mod 4; n=4,6 exhaustive), while at odd n ΔH ≡ 0 mod 2 only (Rédei-trivial; ΔH=2 occurs at n=5 — the mod-4 law is genuinely a parity-of-n phenomenon, the H-analogue of THM-785's ΔE4 ≡ 0/8 mod 16 dichotomy); (b) the blue |ΔH| distribution is near-flat vs the decaying black one; (c) score majorization-comparability: s(t) and s(t̄) are majorization-COMPARABLE for every tiling (zero incomparable pairs, exhaustive n ≤ 6) — finer than the 1-D C3 flux, following the exact affine law s(t̄) = C − s(t), C = (n−2, n−1, …, n−1, n) (= THM-785 §2's degree transformation); (d) the pure-blue census formula pure_blue(n) = ⌊(n+1)/2⌋ − [n even] (kps-S66 conjecture) is CONFIRMED through n=7 (2,1,3,2,4; n=7 count 4 = codex THM-785 §6's independent atlas)
status: CLAIMED+VERIFIED (kind-pasteur-2026-07-14-S128 cont.5): (a),(b) exhaustive n=4..6, n=7 job running (single-threaded canonicalization behind the THM-741 overnight run); (c) exhaustive n≤6, n=7 pending; (d) n≤7 via two independent atlases (mine n≤6 + codex n=7). No general-n proofs yet — the mod-4 blue law and the comparability law are the two proof targets; note H is NOT a function of the degree sequence, so none of this follows from THM-785/787
source: kind-pasteur-2026-07-14-S128 (cont.5; owner directive: trace the flow of transitivity, quantify symmetry/imbalance, build the ordering)
depends_on: []
related:
  - THM-785 (codex-S9: C3 flux, closed blue binomial law, categorical interface (13)/(14), oriented quotient imbalance, the flow address (15)) — the frame this companions; my independent n≤6 run reproduced its interface law and node censuses (PB,M,PK) = (3,5,2), (2,10,22) exactly before reading it
  - THM-787 (opus-S304: E4 axis, blue-avoids-pure-black proved, parity/max conjectures proved by 785)
  - HYP-6855 (this session's log entry); kps-S66/HYP-4997 (blue/black line counts; the pure-blue census conjecture now confirmed)
  - MINIMAL-INVARIANT ORDERING remark: the 4-stage order (score-seq via (Σs² desc, lex), phase PB<MX<PBk, H, canonical word) resolves 3/10/28 of 3/10/34 classes at n=4/5/6 before the single arbitrary stage (canonical word breaks 0/1/6 residual ties) — a low-coordinate alternative to THM-785's 8-coordinate flow address when minimal arbitrariness is the goal; H does real separating work in it (stage gains +0/+2/+11)
---

# THM-790 — the H-companion laws

THM-785/787 put the transitivity flow on degree-derived axes (C3, E4). H — the Hamiltonian-path
count, the invariant this project is built on (Rédei; tilings·|Aut| = H) — is not determined by
the degree sequence, so the flow's H-behaviour is independent structure. Measured exhaustively
(all complement-tiling lines, strict explorer conventions):

**(a) The blue mod-4 law (even n).** Along every blue line at n=4, 6: H(cls(t̄)) ≡ H(cls(t)) mod 4.
At odd n only the Rédei-trivial mod 2 holds (n=5 has blue ΔH = 2, 10). The even/odd-n dichotomy
mirrors THM-785's ΔE4 ≡ 0 vs 8 (mod 16) — but for H, where no degree-sequence argument applies.
n=4: |ΔH| ∈ {0,4}; n=6: {0,4,8,12,16,20,24,28} (all ≡ 0 mod 4, distribution 5,6,5,6,6,1,2,1).

**(b) Blue is flat, black decays.** The blue |ΔH| distribution is near-uniform across its support
(n=6 above); black's decays (n=6: 48,38,48,66,50,42,36,44,12,30,24,18,10,10,4 at |ΔH|=0..28 even).
The H-axis sees the same "blue = long symmetric strides, black = short diffusive steps" picture
as the C3 axis, on an invariant the C3 laws cannot reach.

**(c) Majorization comparability.** s(t̄) = C − s(t) exactly (C above; = THM-785's §2 transform
summed). Over all lines at n ≤ 6: the score multisets of the two ends are ALWAYS majorization-
comparable (no incomparable pair exists; blue n=6: 16 toward-regular / 12 level / 4 toward-
transitive; black: 240/108/132). So the flow is not merely variance-monotone (C3) — it moves
along CHAINS of the majorization lattice, with the blue layer 4:1 regular-directed at n=6.

**(d) Pure-blue census.** pure_blue(n) = ⌊(n+1)/2⌋ − [n even] = 2,1,3,2,4 for n=3..7 (confirmed;
n=7 via THM-785 §6's (4,84,184)). The pure-blue nodes are the flow's sources; at odd n they
include BOTH endpoints of the spectrum (transitive AND regular — n=5: {transitive, (0,2,2,2,4),
regular}), at even n only the transitive end plus interior nodes (n=6: {transitive, (1,1,1,4,4,4)}).

## Evidence log

- [x] n=4,5,6 exhaustive (transitivity_flow_merged_metagraph_kps_S128c5.py + .out): all laws as stated;
      independent reproduction of THM-785's interface law + node censuses before reading it
- [ ] n=7 exhaustive (job running): (a) prediction ΔH ≡ 0 mod 2 with mod-4 violations present (odd n);
      (c) comparability test over 2^15 lines; (d) identify the four pure-blue classes' score seqs
- [ ] proof targets: (a) via the grid-sym anti-automorphism's action on path count parity classes;
      (c) via s̄ = C − s and the near-constancy of C
