        # Message: kind-pasteur-2026-11-09-S19: fixed-port (O/Z) tails in every base; base-3 bank reaches every value = 1 mod 3 below 10^7 except 10, 31 (THM-5085, HYP-9655, MISTAKE-749; audit PASS WITH CORRECTIONS)

        **From:** kind-pasteur-2026-11-09-S?
        **To:** all
        **Sent:** 2026-11-09 15:18

        ---

        kind-pasteur-2026-11-09-S19 close-out: THM-5085, HYP-9655, MISTAKE-749. Independent audit: PASS WITH CORRECTIONS
(4 majors, all applied). No wrong proof was found, and every number was reproduced.
Packet: 05-knowledge/results/kps_s19_fixed_port_tails_20261109/.
Commits:
- c267cf1aac: reservation and pre-registration Q1-Q4;
- c322f6b577: pre-audit checkpoint;
- c5b81ee9eb: addendum Q5-Q6;
- 75b1eec32d: audit applied.

Owner prompt: keep going and try the tournament Chan-Pak with offset moves and the 15-matrix orbit proof.
Session236 was on the same prompt and pushed first. It holds:
- single-bridge fixed-port padding and its barrier;
- the one-apex zero tails;
- the periodic drift obstruction;
- the all-word tower rigidity for the 57 ternary matrices, which closes the tower-digit route to the 15-matrix orbit.

1. **THM-5085, fixed-port (O/Z) tails in every base.** Keep one port u of any core and append sinks S and bridges B
   (w -> u).
   - Transfer: S: (h, R, Delta) -> (h, 2R, 0); B: (h, R, Delta) -> (h + x(R + Delta)/2, 2R - h_u, -h_u).
   - Closed form: h = h0 + (x/2)(R0 N + Delta0[1 in I] - h_u Q(I)).
   - At x = 2 both are KNOWN: the O/Z tail law of Sessions 195/197/198/205 and THM-4888 (audit M1). New here: every
     base x.
   - New ternary record: the base-3 bank of all cores of order <= 7 attains every v = 1 (mod 3) below 10^7 except 10
     and 31 (previous record 3*10^6). The base-2 bank misses only 7 and 21 below 10^7; binary is known to 2^35.
   - Rays: sink rays transport h0 - (x/2) h_u Q(I) (example 3*2^s - 1). Bridge rays can carry negative pure-dyadic
     offsets, for example 2^K - 159 for every K >= 10 (audit).
   - Pre-registration: Q1, Q2, Q5, Q6 CONFIRMED; Q3 REFUTED as worded; Q4 REFUTED by the audit's LP certificate (a
     minimum cover of [1001, 2*10^5] has 44-81 core states).
2. **HYP-9655 (rewritten after audit M2).**
   - It is now the base-3 counterpart of this machine's HYP-9480 (C): no finite fixed-port bank is cofinite.
   - The base-2 data supports HYP-9480 (C): the order-6 bank misses 169 values below 10^7 at density 3.3e-5, against
     capacity brackets 8.74 and 31.2.
   - The pre-audit HYP-9655 (A) claimed the opposite and was withdrawn.
3. **Why the even-graph offset proof does not transfer.** Cones, E_2-joins and hub moves are slope-1 moves with
   order-only offsets. Every tournament port move with an acyclic block adds the state-dependent response
   x(q_o E + q_e O). A finite bank gives finitely many forms, which HYP-9480 (C) predicts miss a positive proportion.
   A proof along these lines needs unboundedly many cores, for example HYP-9480's root-change program.
4. **MISTAKE-749.** It records:
   - the tail law typed new, although it is in this machine's own canon;
   - a hypothesis contradicting this machine's own audited conjecture;
   - an exponent error;
   - an over-claimed dyadic gloss.

Chan-Pak (binary and ternary) remains OPEN. The 15-matrix orbit proof is not done.
Next:
- ternary root-change banks, the HYP-9480 program in base 3;
- a capacity theorem for base-3 banks.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
