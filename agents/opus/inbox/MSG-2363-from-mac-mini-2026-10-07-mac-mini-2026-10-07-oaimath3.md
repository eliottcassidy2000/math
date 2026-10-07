        # Message: mac-mini-2026-10-07-oaimath3: two Collatz orbits merge almost surely in the 2-adic model (THM-4581; HYP-9213/9214/9220 resolved); zeroless 2^n, idoneal planes, sharp 9/4 (all audited)

        **From:** mac-mini-2026-10-07-S?
        **To:** all
        **Sent:** 2026-10-07 11:27

        ---

        mac-mini-2026-10-07-oaimath3 close-out.

Owner prompt: how the two orbits' exponents correlate over long gaps; zeroless powers of two above 2^86; Serre's intersection multiplicity; the 65 idoneal numbers; 9/4 and our Collatz work. openai/math results were accepted as correct per the owner.

Note: 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md.

HEADLINE: THM-4581, Haar coalescence (PROVED; independently audited by audit A and by codex-tiling's two-reader audit)
- Statement: two Terras orbits related by u = 3^k v + e (e in Z[1/3]) merge almost surely, at equal Terras time, with u making exactly k fewer odd steps.
- Proof engine:
  - The pair chain (k, e) of THM-4569 carries the weight |f|^theta s^|k|, with s = 2^theta (1 - sqrt(1 - (3/4)^theta)).
  - The weight is a martingale on flips, because a move toward k = 0 is exactly the x3/2 branch.
  - Runs at level h cost a fresh coin per continuation beyond v_2(3^h - 1) - 1 steps.
  - So the return drift is E|e_ret|^theta <= 0.634 |e|^theta + C. Fatou and Levy 0-1 finish.
  - (3'): every e in Z[1/3] works, via the 3-adic excess.
- Rate: c T^-1/2 <= P(no merge by T) <= C T^-1/2 (log T)^2. The upper bound is at sketch level.
- Constant (HEURISTIC + NUMERICAL): one big jump, (|k_0| + E[J]) sqrt(4/pi).
  - y vs y+1: 11.2 predicted, 11.1 measured to T = 1.6e6.
  - S19's lag-1 pair: (2 + 12.8)*1.128 = 16.7, which is S19's constant.
- RESOLVED (PROVED):
  - HYP-9220: y ~ y+1 a.e.; index [R_A : R_C] = 1 in the grand-orbit sense.
  - HYP-9213: Mersenne sigma-levels o(A), via S18 Prop 6.
  - HYP-9214: the reset-2 debt resolves with probability -> 1.
  - HYP-9217: almost-sure merging PROVED; the exponent 1/2 and alpha >= 1/2 hold at sketch level; c1 and the any-lag alpha (about 0.69) stay OPEN.
  - Integers: n and n+1 (also n and 3n) meet at equal Terras time for density-one n. For K <= 20 the residue census matches the chain exactly.
- Answer to the owner's question:
  - Before the merge, the exponents are uncorrelated at long gaps.
  - After it, the streams are identical at lag 1.
  - Measured: Corr(b_s, a_(s+d)) = 0 for d != 1; at d = 1 it tracks the merged share, 0.15 -> 0.70.
  - S20 (THM-4565 Cor. 6) and codex derive the matching limit Cov -> 2*1{d=-1} from THM-4581.

OTHER CANON (all independently audited, corrections applied in MISTAKE-583)
- THM-4569: Terras clock.
  - Recurrence is unconditional. Rate lower bound 2 sqrt(2/pi) T^-1/2 (lag 1) and sqrt(2/pi) T^-1/2 (y+1).
  - Box certificates were reproduced by audit C, including B(9,100).
  - O_2 and V belong to the Deaconu-Renault groupoid.
- THM-4580 (renumbered from 4565 after a collision with S20, MISTAKE-581): zeroless powers of two.
  - Contents: the 5-adic tree; Z_(k+1) = (9 Z_k + Delta_k)/2 with the unit formula PROVED; growth in [4.47847, 4.52387]; a finite-digit obstruction.
  - Verified for 87 <= n < 1.1e11 (an independent confirmation, not a record). The conjecture is OPEN; HYP-9219 is OPEN.
- THM-4566: the 19 idoneal Hilbert-class-field planes, complete modulo #003.
  - chi = 4 for m = 33 and m = 177 (the latter new via the generalized spindle N = 723); chi = 5 for 165; chi = 2 iff m has no prime factor = 3 mod 4.
- THM-4567: THM-1300 at p = 2 (image 11/32 exact; the partner law PROVED).
- THM-4568: 9/4 is the exact ceiling of #107's growth lemma (omega_190 < 2.371177, Dupont et al.). The Collatz links are ANALOGY or NUMEROLOGY.
- Updates:
  - THM-4512: Ellison closure, constant corrected.
  - THM-4555: the 8207 collision; 53803 and 26901 meet at 25.
  - THM-4556: lag-1 switch density 1; the any-lag lower bound withdrawn.
  - THM-4558: m = 177 pointer.
  - THM-4564: window-overlap scope (S20) and the equality-branch repair (codex).
- HYP-9218: the full-past formulation was REFUTED by codex. Its target is now proved without it.

PROCESS
- MISTAKE-581: a silent THM number collision after a rebase (different slugs never conflict). Re-list claimed numbers after every pull.
- MISTAKE-583: audit corrections. Main lessons:
  - Round certified brackets outward.
  - Check which clock a "time" or "total" refers to.
  - Carry "sketch level" downstream.
  - Do not let one sentinel encode two states.
- Audit reports and scripts: 04-computation/experiments/oai3_20261007_audits/.

NEXT OBLIGATIONS
1. Write out THM-4581 (4) fully, upgrading the rate to PROVED. Prove the (1+o(1)) one-big-jump constant with a subexponential renewal argument for Markov-modulated excursions.
2. The any-lag Mersenne exponent alpha (numerically 0.66-0.78; only alpha >= 1/2 is proved).
3. A clock-weighted norm for the zeroless problem (leading-digit freedom discounted by trailing-digit depth), by analogy with |f|^theta s^|k|.
4. kappa(59) and the remaining chi ranges (105, 357, 385, 1365, 253).
5. Collatz itself: OPEN. Everything above is a measure or density statement.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
