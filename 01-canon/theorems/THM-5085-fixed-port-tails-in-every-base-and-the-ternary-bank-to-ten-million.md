---
id: THM-5085
title: "Fixed-port (O/Z) tails in every base, and the ternary bank to 10^7. Keep a port u of any core tournament T0 and append letters S (a sink) and B (a vertex beaten by all except u, with w -> u). (1) The port state (h, R, Delta) = (I(Omega,x), E+O, E-O) moves by S: (h, 2R, 0) and B: (h + x(R+Delta)/2, 2R - h_u, -h_u), h_u = I(Omega(T0 - u), x) invariant. (2) For bridge positions I, h = h0 + (x/2)[R0 N(I) + Delta0 [1 in I] - h_u Q(I)]. Both are KNOWN at x = 2 (the O/Z tail law: Sessions 195/197/198/205, THM-4888); the base-x versions are new. (3) Sink rays transport the offset h0 - (x/2) h_u Q(I); bridge rays B^k transport h - (x/2)(R - h_u) and can carry negative offsets on pure dyadic leading parts (2^K - 159 for every K >= 10). (4) The banks of all cores of order <= 7 miss exactly {7, 21} below 10^7 in base 2 and exactly {10, 31} below 10^7 in base 3 (ternary record; previously 3*10^6); the order-<= 6 bank misses 169 odd values below 10^7."
status: >
  KNOWN at x = 2:
  - (1) is the K/D (= Z/O) law of session195/combined-growth section 2 (also Session197 ONE_CORE_BINARY_INVERSE (1)).
  - (2) is the O/Z tail formula of Session205 finite-oz-tail-bank/ARITHMETIC_REDUCTION. Also Session197 ONE_CORE (6),
    ROOTED_ORDERED_AMALGAMATION (7)-(8), and Session198 DIGIT_BLOCK (1), whose pair sum C(A) is Q(I). It is restated in
    THM-4888 (kind-pasteur S13).
  - The dictionary at x = 2 is e = -Delta (e = Hamiltonian paths ending at the port) and d = h_u.
  PROVED (new): the base-x versions of (1) and (2), with remainders counted by I(Omega(T - V(P)), x) and
  Delta = E - O. The proof is the x = 2 proof verbatim.
  - Checked: 800 random tails of order <= 15 in bases 2 and 3 (fp_check.py). The audit adds 19,811 exhaustive and
    2,100 random instances in bases 2, 3 and 4, with 0 mismatches.
  (3) is a corollary. The sink-ray part follows from Q(2m) = Q(m). The bridge-ray formula and the 2^K - 159 example
  are the audit's, brute-forced at orders 9-15.
  FINITE-EXACT: (4), by fp_cover.c; the audit reproduced it with two independent engines.
  Pre-registration:
  - Q1, Q2 (10^6) CONFIRMED; Q5, Q6 (10^7) CONFIRMED.
  - Q3 REFUTED as worded. Its only predictive clause, at most 500 misses, holds.
  - Q4 REFUTED by the audit's LP certificate: every cover of [1001, 2*10^5] needs at least 44 core states in base 2
    (at least 50 in base 3); greedy covers use 81 and 86.
  AUDIT: PASS WITH CORRECTIONS (05-knowledge/results/kps_s19_fixed_port_tails_20261109/audit/AUDIT.md). The M1 credit,
  the M2 conflict with HYP-9480, the M3 example and the M4 dyadic gloss are applied; MISTAKE-749.
  NOT CLAIMED: tournament or ternary Chan-Pak (both OPEN). HYP-9480 (C) conjectures that every finite bank, these
  included, misses a positive proportion beyond any census.
session: kind-pasteur-2026-11-09-S19
source: 05-knowledge/results/kps_s19_fixed_port_tails_20261109/README.md
proof: 05-knowledge/results/kps_s19_fixed_port_tails_20261109/PROOF.md
related:
  - HYP-9655
  - HYP-9480
  - THM-4888
  - THM-5070
  - HYP-9640
  - THM-4940
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - THM-4875
  - 05-knowledge/results/session205/finite-oz-tail-bank/ARITHMETIC_REDUCTION.md
  - 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/offsets/PROOF.md
verification:
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/port_stats.c
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/fp_cover.c
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/fp_greedy.c
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/fp_check.py
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/out/
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/audit/
---

# THM-5085 - fixed-port (O/Z) tails in every base

## Setting

- **Core.** A core is a tournament T0 with a port u.
- **Port state.** For the current tournament T, let:
  - h = I(Omega(T), x);
  - E, O = the sums of I(Omega(T - V(P)), x) over the paths P of T from u of even and odd order;
  - R = E + O and Delta = E - O.

  Also h_u = I(Omega(T0 - u), x).
- **Letters.**
  - **S** (O in the session-195-206 notation) appends a vertex beaten by every existing vertex.
  - **B** (Z or K there) appends a vertex w beaten by every existing vertex except u, with w -> u.
  - The port is never moved. "t transitive vertices, then a bridge" is the word S^t B.
- **x = 2 dictionary.** The session-195-206 notation is related by e = -Delta (Hamiltonian paths ending at u) and
  d = h_u.

## Statement

**(1) Transfer.**
- S: (h, R, Delta) -> (h, 2R, 0).
- B: (h, R, Delta) -> (h + x(R + Delta)/2, 2R - h_u, -h_u).
- h_u is unchanged by both letters.

KNOWN at x = 2; new for general x.

**(2) Closed form.** For a word with bridge set I, from (h0, R0, Delta0):

    h = h0 + (x/2) [ R0 N(I) + Delta0 [1 in I] - h_u Q(I) ],

with N(I) = Σ_{i in I} 2^(i-1) and Q(I) = Σ_{j < i in I} c(i - j), where c(1) = 2 and c(d) = 2^(d-1) for d >= 2.
Equivalently, Q(2m) = Q(m) and Q(2m+1) = Q(m) + m + (m mod 2). KNOWN at x = 2; new for general x.

**(3) Rays.**
- **Sink rays.** For s >= 1, the tail S^s (pattern I) has value o(I) + 2^s (x/2) R0 N(I), with
  o(I) = h0 - (x/2) h_u Q(I).
  - Example: T0 = {u} with I = {1, 2} gives 3 * 2^s - 1 (Session236's one-apex negative offsets, offsets/PROOF.md
    (3), give the same values).
  - For even x, a power-of-2 leading coefficient on a sink ray forces |I| = 1 and the offset h0. For odd x the leading
    coefficient is never a power of 2.
- **Bridge rays** (the audit's formula). From a state (h, R, -h_u), the word B^k gives

      [h - (x/2)(R - h_u)] + 2^k (x/2)(R - h_u).

  Bridge rays can carry negative offsets on pure dyadic leading parts. The order-7 core with rows
  [94, 84, 32, 6, 12, 27, 60] and port 0 has base-2 state (71, 287, -5, 31). After one B, R - h_u = 512, and the value
  of B^(k+1) is exactly 2^(k+8) - 159; this was brute-forced for k = 2..8. So it reaches 2^K - 159 for every K >= 10.
  Session236's barrier concerns its single-final-bridge mechanism and remains correct there.

**(4) Bank coverage.** The bank of order n0 consists of every class of order <= n0 with every port.

| n0 | x | misses (complete list up to the bound) | bound | distinct core states | search states |
|---|---|---|---|---|---|
| 7 | 2 | 7, 21 | 10^7 | 2,831 | 573,476,037 |
| 7 | 3 | 10, 31 | 10^7 | 2,901 | 340,144,258 |
| 6 | 2 | 7, 21 and 167 others (169 in all) | 10^7 | 273 | (audit) |
| 6 | 3 | 10, 31, 5401, 40597, 97282, ... (15 in all) | 10^6 | 278 | (audit) |

**What this adds.**
- **Ternary.** Every v = 1 (mod 3) below 10^7 other than 10 and 31 is I(Omega(T), 3) of a fixed-port tail over a core
  of order <= 7. The previous ternary bound was 3*10^6 (the S18 audit's extension of THM-5070).
- **Binary.** The base-2 bank lies far inside the known range: 2^35, for a one-apex x order-<= 7 bank (HYP-9478 P5-P6,
  HYP-9480 (C)).
- **Minimum covers.** A minimum core cover of the odd targets in [1001, 2*10^5] has between 44 (the audit's LP
  certificate) and 81 (greedy) states. In base 3 the range is [50, 85].

## Proof

The full proof is in PROOF.md.

**(1).**
- The bridge w beats only u. So every new odd cycle is a path P from u of even order, then w. A packing holds at most
  one new cycle, and its remainder is T - V(P). This gives h' = h + xE.
- Paths from u gain an optional last vertex, except (u, w). So f'(z) = (1 + z) f(z) - [B] z^2 h_u, where
  f(z) = Σ_P z^|P| I(Omega(T - V(P)), x).
- This is the session-195 proof with Hamiltonian-path counts of remainders replaced by I(Omega(.), x).

**(2).** Unroll R_i = 2 R_{i-1} - h_u [i in I] and add the bridge increments, as in Session205.

**(3).** Sinks double R and clear Delta. After a bridge, R - h_u doubles at each further bridge.

**(4).** `fp_cover.c` takes the closure with h <= LIM. This is sound because h is nondecreasing. The audit reproduced
every row with a state-closure engine and a closed-form enumeration engine.

## Remarks

- **Relation to the even-graph proofs (THM-5040, THM-5055).** Those are driven by slope-1 moves whose offsets depend
  only on the order: cones, E_2-joins and hub moves. Here S is the transport and B a digit whose weight is the current
  port response R, with drift h_u. No port move with an acyclic block is a cone: its increment x(q_o E + q_e O)
  depends on the port response.
- **Finite banks.** HYP-9480 (C) conjectures that every finite O/Z bank misses a positive proportion of the odd
  numbers, of Poisson order exp(-capacity). The audit's capacity brackets are 8.74 (order <= 6) and 31.2 (order <= 7).
  - The order-<= 6 bank's 169 misses below 10^7 have a stable density of about 3.3*10^-5. This is consistent with a
    positive miss proportion.
  - The order-<= 7 bank's empty miss set below 10^7 is what a proportion near exp(-31) predicts.
  - HYP-9655 records the base-3 counterpart.
