---
id: THM-5070
title: "Port chains: a nonnegative 3-dimensional transfer for odd-cycle packing counts in every base. (1) Gluing a block (B, w, u') after a ported tournament (T, u) - every old vertex beats every block vertex except w -> u, new port u' - maps the state (I(Omega(T), x), E, O) by a nonnegative integer 3 x 3 matrix N_x(B, w, u'), where E, O sum I(Omega(T - V(P)), x) over paths P from the port of even and odd order. (2) Apex bridges: if the odd cycles of each side all pass through its glued vertex, I = I_1 I_2 + x c; for one-apex sides I = (1 + x a(W))(1 + x a(W')) + x[ceil(W/2)(1 + floor(W'/2)) + (1 + floor(W/2)) ceil(W'/2)]. (3) The one-vertex block is the tower, characteristic polynomial t^3 - t^2 - t - (x - 1), = (t - 2)(t^2 + t + 1) at x = 3. (4) Port chains with blocks of order <= 5 attain every v = 1 (mod 3), v <= 10^6, except 10 and 31, and every odd v <= 2*10^5 except 7 and 21; with blocks of order <= 4 the base-3 misses below 10^6 are exactly 10, 31, 61, 106, 217, 274, 337, 346, 607, 865 (already for 15 explicit generators) and the base-2 misses below 2*10^5 are exactly 7, 21, 47, 67, 127, 173, 491. (5) Among all products of <= 3 order-<= 4 base-3 blocks, the only matrices preserving the tower's dominant direction are the tower powers. (6) In bases 4-8 the alpha-sum gaps below 30,000 are exactly the known gaps below 81: {3,14,27}, {3,10,27}, {3,10,14,21,46}, {3,10,14,29}, {3,35,52}."
status: >
  PROVED: (1), (2) and (3).
  - (1) is the base-x packing form of the seam transfer. It is KNOWN at x = 2 through OCF (SEAM.md,
    seven_return_native_recursions_20261029; THM-4875 (B); two-block count Session29 ADDENDUM section 2). At x = 2
    SEAM.md's coordinates are (H, S, R) = (h, O - E, E + O), observed by the audit on 60 chains.
  - THM-4952 (A) is the general-base precedent: a finite-state transfer for I(Omega, x) on span-w ribbon words. A port
    chain with blocks of order <= k is a ribbon word of span <= 2k - 1. The new content is the 3-dimensional
    nonnegative port state.
  - S86's packet asked for such a ternary coverage engine (packing_bases_repunits_20261108/README.md, section 8).
  - (2) specialises the one-arc join law (THM-4866 (B)). Its hypothesis holds exactly for one-apex sides (T - u
    transitive). Its one-apex statistics use the inherited one-apex recurrence
    (session25, Session191 one_apex/PROOF.md (1)).
  - The recurrence in (3) is KNOWN (THM-4952 (C), THM-4865 (c)); the matrix form is new.
  FINITE-EXACT: (4), (5) and (6).
  - (6): the gap sets of bases 4-8 to alpha-sum 30,000.
    - m <= 80: THM-4940 (6). The audit independently confirmed its attained side.
    - 80 < m <= 30,000: port chains, plus 61 hill-climbed witnesses of orders 8-11. The witnesses were checked by a
      second script of this session, and then independently by the audit (all 61 pass).
    - The order-join closure settled nothing in (80, 30000].
  - (4) by saturated breadth-first closure; saturation is exact for nonnegative matrices.
  - Pre-registered P1-P3 all CONFIRMED (PREREGISTRATION.md). The audit ran both order-5 rows directly to 2*10^5, and
    the 15 generators to 3*10^6 (293,669,502 states, the same ten misses).
  - The transfer is checked against an independent brute-force packing count on 600 random chains of orders <= 14
    (0 failures).
  NOT CLAIMED: ternary or binary Chan-Pak, both OPEN; see HYP-9640 (A), (B) for the reductions.
  AUDIT: PASS WITH CORRECTIONS (05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/audit/AUDIT.md;
  M2 over-claim in the gloss of (5); m3, m4, m6, m10, m11 applied; MISTAKE-748).
session: kind-pasteur-2026-11-09-S18
source: 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/README.md
proof: 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/PROOF.md
related:
  - HYP-9640
  - HYP-9530
  - THM-4940
  - THM-4950
  - THM-4952
  - THM-4875
  - THM-4866
  - THM-4865
  - THM-1370 (h-spectrum-omits-7-21-all-n)
verification:
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/chain_blocks.py
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/chain_cover.c
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/chain_check.py
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/packing_count.c
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/tower_direction.py
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/bases_gaps.py
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/hunt.c
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/verify_witnesses.py
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/out/
---

# THM-5070 - port chains

## Setting

- **Ported tournament.** A ported tournament (T, u) has state s(T, u) = (h, E, O):
  - h = I(Omega(T), x);
  - E = sum of I(Omega(T - V(P)), x) over the paths P of T starting at u with |V(P)| even;
  - O = the same over paths with |V(P)| odd. The path (u) is included.
- **Blocks.** A block is a tournament B with marked vertices w, u'.
- **Gluing.** T (.) B is the ordinal sum T => B with the single arc between u and w reversed (w -> u). Its port is u'.
- **Port chains.** A port chain is any result of repeated gluing from ({u}), whose state is (1, 0, 1).

Up to converse, this is SEAM.md's operation X (.) Y.

## Statement

**(1) Transfer.** s(T (.) B, u') = N_x(B, w, u') s(T, u). The matrix N_x has nonnegative integer entries for integer
x >= 0. Its first row is (eta, x q_o, x q_e), where:
- eta = I(Omega(B), x);
- q_o (resp. q_e) sums I(Omega(B - V(Q)), x) over the paths Q of B ending at w with |V(Q)| odd (resp. even).

The E- and O-rows are the explicit path sums of PROOF.md, Lemma C.

**(2) Apex bridges.** Suppose every odd cycle of T passes through u and every odd cycle of B passes through w. In a
tournament this is equivalent to T - u and B - w being transitive (every cycle implies a 3-cycle), i.e. to one-apex
sides. Then

    I(Omega(T (.) B), x) = I(Omega(T), x) I(Omega(B), x) + x c,

where c is the number of odd cycles through w -> u.

For one-apex sides:
- T = TT_k + u, where u beats the chain vertices marked by the binary word W (most significant bit first);
- B = TT_l + w, where the chain vertices beating w are marked by W' (least significant bit first).

Then E = ceil(W/2), O = 1 + floor(W/2), and the value is

    (1 + x a(W))(1 + x a(W')) + x [ceil(W/2)(1 + floor(W'/2)) + (1 + floor(W/2)) ceil(W'/2)],

with a(2W+1) = a(W), a(2W) = a(W) + ceil(W/2), a(0) = 0.

**(3) Tower.** The one-vertex block (w = u') has matrix [[1, x, 0], [0, 0, 1], [1, 1, 0]]. Iterating it builds the
tower 1-^n. Its characteristic polynomial is t^3 - t^2 - t - (x - 1). At x = 3 this is (t - 2)(t^2 + t + 1), with
dominant eigenvector (3, 1, 2).

**(4) Coverage (FINITE-EXACT).** The audit reproduced every row with its own code: every state count matched, and
every reached value has a certificate word. It also ran the order-5 rows directly to 2*10^5, and the 15 generators to
3*10^6 with the same ten misses.

| family | base | misses (complete list up to the bound) | bound |
|---|---|---|---|
| blocks of order <= 5 (271 matrices) | 3 | 10, 31 | 10^6 |
| blocks of order <= 4 (57 matrices); also the 15 generators of `out/generators_base3_15.txt` | 3 | 10, 31, 61, 106, 217, 274, 337, 346, 607, 865 | 10^6 for the 15 (hence for the 57); runs: 2*10^5 (57), 10^6 (15) |
| blocks of order <= 5 (264 matrices) | 2 | 7, 21 | 2*10^5 |
| blocks of order <= 4 (56 matrices) | 2 | 7, 21, 47, 67, 127, 173, 491 | 2*10^5 |
| 13 generators (`out/generators_base2_13.txt`) | 2 | 7, 21, 47, 67, 127, 163, 173, 331, 491, 511, 799, 959 | 2*10^5 |

**Previous records and placement.**
- The base-3 range was 10^5 (span-2 ribbon words, THM-4952 (D)) and 29,035 (order joins of strong tournaments of
  order <= 10, THM-4950 (3)).
- The base-2 rows are a control: THM-4875 (D) reaches 2^32 with one-arc joins of one-apex pieces.

**(6) Gap sets in bases 4-8 to alpha-sum 30,000 (FINITE-EXACT).** For b = 4, 5, 6, 7, 8 the set of gaps
Gap_b = {m : 1 + b m is not I(Omega(T), b) for any tournament T} satisfies Gap_b ∩ [0, 30000] = Gap_b ∩ [0, 80], namely:

| b | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|
| gaps below 30,000 | 3, 14, 27 | 3, 10, 27 | 3, 10, 14, 21, 46 | 3, 10, 14, 29 | 3, 35, 52 |

Proof sketch:
- **m <= 80.** THM-4940 (6) gives the exact gap sets here.
- **80 < m <= 30,000.** The values are covered by:
  - the order-<= 5 port chains in base b;
  - 61 explicit tournaments of orders 8-11 for the alpha-sums the chains miss: 1 in base 4, 7 in base 5, 6 in base 6,
    7 in base 7, 40 in base 8.

  Closing the chain values under order joins was tried and added nothing in this range.
- **Witness check.** The 61 witnesses were found by hill-climbing. A second script of this session verified them from
  scratch (cycle DP, packing DP, OCF; `out/bases_4to8_witnesses_verify.out`). The audit verified all 61 again
  independently.

This is evidence for HYP-9530 (B) (finite gaps in every base); see HYP-9640 (E).

**(5) Tower directions (FINITE-EXACT).**
- Take all products of length <= 3 of the 57 base-3 blocks of order <= 4. The distinct ones that map (3, 1, 2) to a
  multiple of itself are exactly M_T^k, k = 1..12.
- So, among products of length <= 3 (and, by the audit, of length exactly 4), no block product other than a tower
  power preserves the tower direction. Longer products are not covered by this check.

## Proof

The full proof is in PROOF.md (Lemmas A-C, sections (2)-(5)).

**(1).**
- The only backward arc is w -> u, so every new cycle uses it once, and all new cycles share u. This gives the first
  row.
- A path from the new port either stays in B, or runs Q_1 (to w), then P_T (in T, from u), then Q_2 (back in B).
- Each remainder is either an ordinal sum, whose count is a product, or a smaller instance of the first row. Every
  coefficient is a packing count of a subtournament of B.

**(2).** Each side's remainders are acyclic once the glued vertex is removed. The one-apex statistics count increasing
chain paths by parity.

**(3).** The case analysis of Lemma C for B = {w}, then a 3 x 3 determinant.

**(4).** `chain_cover.c` takes the breadth-first closure of (1, 0, 1), keeping h <= LIM and saturating E, O at
floor(LIM/x) + 1.
- **Pruning.** Pruning at h > LIM is sound because h never decreases: h' >= eta h >= h, since eta = I(Omega(B), x) >= 1.
- **Saturation is exact.** Each h-row is (eta, x q_o, x q_e), so a saturated coordinate with a positive coefficient
  forces h > LIM. For nonnegative integer rows, cap(row . s) = cap(row . cap(s)), where only E and O are capped.
- **Audit check.** The audit obtained identical miss sets with this cap, with cap = LIM, and with no cap.
- The order-5 family contains the order-4 family. So the order-5 runs to 3000 settle all misses below the
  order-4 bounds.
- The 10^6 base-3 run used 64,420,757 states.
- `chain_check.py` compares the transfer with `packing_count.c` (Hamiltonian-cycle DP plus set-partition DP) and with
  OCF on 600 random chains of orders <= 14: 0 failures.

**(5).** `tower_direction.py`, output `out/tower_direction.out`.

## Remarks

- **What a proof would need.** HYP-9640 (A), (B) restate ternary and binary Chan-Pak as the cofiniteness of the orbit
  of (1, 0, 1) under 15 (resp. 56) explicit nonnegative 3 x 3 integer matrices.
  - (5) shows that the cheapest route, a common eigen-direction with affine digits, does not appear among short
    products of order-<= 4 blocks.
  - Exploratory runs (not claimed) show that the apex-bridge family alone, and the 2-dimensional subsystem made by
    appending a sink after each block, both leave gaps that keep growing.
- **Base 2 and base 3 together.** The same order-<= 5 blocks (387 classes) realise both Chan-Pak spectra to the
  stated bounds. The matrices depend on the base: 264 distinct at x = 2, 271 at x = 3. The only misses are the permanent gaps: alpha-sums 3 and 10 in both bases (THM-4940 (3), (4); THM-1370).
