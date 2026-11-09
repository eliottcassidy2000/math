---
id: THM-5085
title: "Fixed-port tails: tournament offset moves with a closed form. Keep a port u of any core tournament T0 and append letters S (a sink) and B (a vertex beaten by all except u, with w -> u). (1) The state (h, R, Delta) = (I(Omega,x), E+O, E-O) of the port moves by S: (h, 2R, 0) and B: (h + x(R+Delta)/2, 2R - h_u, -h_u), with h_u = I(Omega(T0 - u), x) invariant. (2) For bridge positions I, h = h0 + (x/2)[R0 N(I) + Delta0 [1 in I] - h_u Q(I)], N(I) = sum 2^(i-1), Q(I) = sum over pairs j < i in I of c(i-j), c(1) = 2, c(d) = 2^(d-1). (3) Prepending sinks transports the signed offset h0 - (x/2) h_u Q(I) along the ray 2^s (x/2) R0 N(I); pure dyadic rays force one bridge and a positive H-valued offset. (4) Fixed-port tails over all cores of order <= 7 realise every odd value below 10^6 except 7, 21 (base 2) and every value = 1 (mod 3) below 10^6 except 10, 31 (base 3); cores of order <= 6 leave 18 further sparse misses in base 2."
status: >
  PROVED: (1), (2) and (3).
  - (1) applies THM-5070 Lemma B (one-arc gluing) and the path-response bookkeeping of THM-5070 Lemma C, here with the
    port kept fixed.
  - (2) and (3) follow by induction on the word.
  - Brute-force check: 800 random tails of order <= 15 in bases 2 and 3, 0 mismatches (fp_check.py).
  FINITE-EXACT: (4), by fp_cover.c over every class and port of order <= 7 (classgen + port_stats).
  - Pre-registered Q1 and Q2 CONFIRMED.
  - Q3 REFUTED as worded (its "all misses above 30,000" omitted the permanent gaps 7 and 21); the intended statement
    is CONFIRMED.
  - Q4 UNDECIDED: a minimum core cover of [1001, 2*10^5] has between 17 and 81 core states.
  KNOWN and credited:
  - The one-apex recurrence (session25; Session191 one_apex/PROOF.md (1)) is the core T0 = {u}.
  - The single-bridge rays, their positive-offset barrier, the one-apex zero tails, the periodic-suffix drift
    obstruction and the all-word tower-ray rigidity are Session236's
    (05-knowledge/results/session236_chanpak_offsets_orbits_20261109, concurrent, pushed first). Here (3) extends the
    barrier to multi-bridge pure-dyadic rays, and shows that non-dyadic multi-bridge rays carry negative offsets.
  - The even-graph offset template is THM-5040 (opus-S89) and THM-5055 (mac-mini).
  - The port-chain transfer is THM-5070.
  NOT CLAIMED: tournament or ternary Chan-Pak (both OPEN; HYP-9655 (A) is the reduction to finitely many tail forms).
session: kind-pasteur-2026-11-09-S19
source: 05-knowledge/results/kps_s19_fixed_port_tails_20261109/README.md
proof: 05-knowledge/results/kps_s19_fixed_port_tails_20261109/PROOF.md
related:
  - HYP-9655
  - THM-5070
  - HYP-9640
  - THM-5055
  - THM-5040
  - THM-4875
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/offsets/PROOF.md
verification:
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/port_stats.c
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/fp_cover.c
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/fp_greedy.c
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/fp_check.py
  - 05-knowledge/results/kps_s19_fixed_port_tails_20261109/out/
---

# THM-5085 - fixed-port tails

## Setting

- **Core.** A core is a tournament T0 with a port u.
- **Port state.** For the current tournament T (T0 plus a tail), let:
  - h = I(Omega(T), x);
  - E, O = the sums of I(Omega(T - V(P)), x) over the paths P of T from u of even and odd order;
  - R = E + O and Delta = E - O.

  Also h_u = I(Omega(T0 - u), x).
- **Letters.**
  - **S** appends a vertex beaten by every existing vertex.
  - **B** appends a vertex w beaten by every existing vertex except u, with w -> u.
  - The port is never moved. The letter "t transitive vertices, then a bridge" is the word S^t B.

## Statement

**(1) Transfer.**
- S: (h, R, Delta) -> (h, 2R, 0).
- B: (h, R, Delta) -> (h + x(R + Delta)/2, 2R - h_u, -h_u).
- h_u is unchanged by both letters.

**(2) Closed form.** For a word over {S, B} with bridge set I, from (h0, R0, Delta0):

    h = h0 + (x/2) [ R0 N(I) + Delta0 [1 in I] - h_u Q(I) ],

with N(I) = Σ_{i in I} 2^(i-1) and Q(I) = Σ_{j < i in I} c(i - j), where c(1) = 2 and c(d) = 2^(d-1) for d >= 2.

**(3) Signed offset transport.** For every bridge pattern I and s >= 1, the tail S^s (pattern I) has value

    o(I) + 2^s (x/2) R0 N(I),    o(I) = h0 - (x/2) h_u Q(I).

- o(I) can be negative. For example, T0 = {u} with I = {1, 2} gives 3 * 2^(s+1) - 1.
- If the leading coefficient 2^s (x/2) R0 N(I) is a power of 2, then |I| = 1 and the offset is h0, the value of T0.

**(4) Coverage.** Run over every tournament class of order <= n0 and every port, with letters S, I_0 = B, I_1 = SB and
I_2 = SSB:

| n0 | x | misses below 10^6 |
|---|---|---|
| 7 | 2 | 7, 21 |
| 7 | 3 | 10, 31 |
| 6 | 2 | 7, 21, 34869, 85885, 103061, 116587, 153451, 245923, 535405, 594793, 595543, 609873, 618723, 767135, 774073, 782103, 824119, 826373, 868367, 949777 |

## Proof

The full proof is in PROOF.md.

**(1).**
- The only new arc into the old vertices is w -> u. Every new odd cycle is a path P from u of even order, followed by
  w. Packings contain at most one new cycle, and its remainder T - V(P) has no extra cycles, because w's other arcs
  point into deleted vertices. This gives h' = h + xE.
- Paths from u in the new tournament are old paths, optionally extended by the new vertex. The exception is (u, w),
  forbidden since w -> u. So f'(z) = (1 + z) f(z) - [B] z^2 h_u, where f(z) = Σ_P z^|P| I(Omega(T - V(P)), x).
- Evaluate at z = +1 and z = -1.

**(2).** Unroll R_i = 2 R_{i-1} - h_u [i in I] and add the bridge increments. Consecutive bridges contribute the extra
-h_u from Delta, which turns c(1) from 1 into 2.

**(3).** Sinks double R and clear Delta.

**(4).** `fp_cover.c` takes the closure with h <= LIM. This is sound because h is nondecreasing. A state whose
cheapest bridge overflows is a leaf, since Delta <= R. The pre-registered outcomes and state counts are in
PREREGISTRATION.md.

## Remarks

- **Relation to the even-graph proofs.** Those are driven by slope-1 moves whose offsets depend only on the order:
  cones, E_2-joins and hub moves. Tournament port moves add x(q_o E + q_e O), which depends on the port response. That
  response is a pure path count only at an apex.
  - S is the exact analogue of the order-shifting transport (THM-5055 Lemma 6).
  - B is a digit whose weight is the current response R, with drift h_u.
- **Status of Chan-Pak.** It remains OPEN in both bases. HYP-9655 (A) reduces each to the cofiniteness of finitely many
  forms h0 + (x/2)(R0 N - h_u Q(N) + Delta0 [N odd]).
