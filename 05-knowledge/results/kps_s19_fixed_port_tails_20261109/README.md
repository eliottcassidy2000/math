# Tournament offset moves as fixed-port (O/Z) tails, in every base

**Session:** kind-pasteur-2026-11-09-S19. **Canon:** THM-5085, HYP-9655.

**Status.**
- KNOWN at x = 2: the transfer and the closed form (Sessions 195/197/198/205; THM-4888).
- PROVED (new): the base-x versions.
- FINITE-EXACT: the bank coverage to 10^7 in bases 2 and 3.
- OPEN: tournament and ternary Chan-Pak.

**Independent audit:** PASS WITH CORRECTIONS ([audit/AUDIT.md](audit/AUDIT.md)). Its four majors are all applied:
- M1: x = 2 is KNOWN;
- M2: the pre-audit HYP-9655 contradicted HYP-9480 (C) and was withdrawn;
- M3: an exponent in an example;
- M4: an over-claimed dyadic gloss.

It also decided Q4 (REFUTED) and confirmed Q5 and Q6. MISTAKE-749.

Proofs are in [PROOF.md](PROOF.md); the predictions are in [PREREGISTRATION.md](PREREGISTRATION.md).

## Owner prompt (2026-11-09)

> keep going and try the tournament Chan-Pak with offset moves and the 15-matrix orbit proof

## Inheritance (read the audit's M1 first)

**The fixed-port tail family is the "O/Z tail bank" of Sessions 195-206.**
- session195/combined-growth: the K/D laws.
- Session197 path-polynomials: ONE_CORE_BINARY_INVERSE (1) and (6), ROOTED_ORDERED_AMALGAMATION (7)-(8).
- Session198 DIGIT_BLOCK (1): its pair sum is Q(I).
- Session205 finite-oz-tail-bank ARITHMETIC_REDUCTION: the tail formula H_s(m) = h0 + R0 m - d Q(m) - e0 (m mod 2).
- Session206: the capacity theorems and the root obstruction.
- This machine's S13: THM-4888 and HYP-9480.

At x = 2 the transfer and the closed form of THM-5085 are these results, under the dictionary e = -Delta and d = h_u.
This session missed that priority before its audit (MISTAKE-749).

**Session236** (same prompt, pushed first) holds:
- the single-final-bridge padding and its barrier;
- the one-apex zero tails, including negative one-apex offsets such as 6 * 2^m - 1;
- the periodic-drift obstruction;
- the all-word tower-ray rigidity for the 57 ternary matrices. This upgrades S18's finite THM-5070 (5) and closes the
  tower-digit route to HYP-9640 (A).

**Even-graph offset proofs:** THM-5040 (opus-S89) and THM-5055 (mac-mini).

## What is new here

**(a) Every base.** The transfer S: (h, R, Delta) -> (h, 2R, 0) and B: (h, R, Delta) -> (h + x(R + Delta)/2, 2R - h_u,
-h_u) holds for the base-x packing count. So does the closed form

    h = h0 + (x/2) [ R0 N(I) + Delta0 [1 in I] - h_u Q(I) ].

- Remainders are counted by I(Omega(T - V(P)), x).
- Delta = E - O replaces -e, which has no base-x meaning.
- Checked: 800 random tails in bases 2 and 3; the audit adds 19,811 exhaustive and 2,100 random instances in bases 2,
  3 and 4. All with 0 mismatches.

**(b) The ternary bank to 10^7.** Over all cores of order <= 7 (2,901 distinct core states), the base-3 fixed-port
tails attain every v = 1 (mod 3) below 10^7 except 10 and 31 (pre-registered Q2, Q6 CONFIRMED).
- This is a new ternary bound: the previous one was 3*10^6.
- Every reached value is the packing count I(Omega(T), 3) of an explicit tournament.

**(c) Base-2 bank data** (pre-registered Q1, Q5 CONFIRMED):
- The order-<= 7 bank misses only 7 and 21 below 10^7. This is far inside the known binary range 2^35, which is
  reached by a one-apex x order-<= 7 bank (HYP-9478, HYP-9480 (C)).
- The order-<= 6 bank misses 169 odd values below 10^7, at a stable density of about 3.3*10^-5.
- Both are consistent with HYP-9480 (C): every finite bank misses a positive proportion, of Poisson order
  exp(-capacity), and the capacity brackets are 8.74 and 31.2. HYP-9655 states the base-3 counterpart.

**(d) Minimum bank sizes.** A minimum core cover of the odd targets in [1001, 2*10^5] has between 44 and 81 states:
44 from the audit's LP certificate, 81 from greedy. In base 3 the range is [50, 85]. So pre-registered Q4 ("at most
40") is REFUTED.

**(e) Rays.**
- **Sink rays.** S^s (pattern I) transports o(I) = h0 - (x/2) h_u Q(I).
  - Example: the core {u} with I = {1, 2} gives 3 * 2^s - 1.
  - On sink rays a pure dyadic leading part forces one bridge and the offset h0.
- **Bridge rays.** B^k from a state (h, R, -h_u) transports h - (x/2)(R - h_u) (the audit's formula).
  - These can carry negative offsets on pure dyadic leading parts: an order-7 core gives 2^K - 159 for every K >= 10.
  - Session236's barrier is scoped to its single-final-bridge mechanism and remains correct there.

## Why the even-graph proof does not transfer

Both even-graph proofs build every odd count by writing leading bits on top of an arbitrary smaller count. The moves are
the cone (d -> 2^n + d), the E_2-join (d -> 3*2^n - 2 + d) and the hub moves. Each has slope 1, and its offset depends
only on the order.

A tournament port move with an acyclic block adds x(q_o E + q_e O) (THM-5070 (1)). That depends on the port response.
So no port move of this type is a cone. Within fixed-port tails:
- S is the exact counterpart of the order transport (THM-5055 Lemma 6);
- B is a digit whose weight is the current response R, with drift h_u.

A finite set of cores gives only finitely many such forms, and HYP-9480 (C) predicts that they miss a positive
proportion. A proof of Chan-Pak along these lines needs unboundedly many cores. One candidate is HYP-9480's
root-change program; another is a different mechanism.

## The 15-matrix orbit (HYP-9640 (A))

**Not proved.** Session236's all-word tower rigidity rules out the digit route along the tower at every word length.

THM-5085 (b) supplies an independent base-3 family, with a closed form, to 10^7. By HYP-9655 (A), no finite bank of
this type is expected to be cofinite, so it is evidence, not a proof route.

## Files

| file | purpose |
|---|---|
| `port_stats.c` | (h, E, O, h_u) in bases 2 and 3 for every class and port (input: classgen output) |
| `fp_cover.c` | bank coverage (letters S and I_0, I_1, I_2; I_t = S^t B) |
| `fp_greedy.c` | greedy core cover and the disjoint-cover-set lower bound (Q4) |
| `fp_check.py` | brute-force check of the transfer on random tails, using the S18 `packing_count.c` |
| `fixedport.py` | exploratory Python version of `fp_cover.c` |
| `out/` | outputs of every run quoted above, and the SHA-256 of the regenerable port-statistics file |
| `audit/` | the independent audit: report, programs, outputs, `run_all.sh` |

**Reproduce.**
1. Generate the classes: `gcc -O2 -o classgen ../first_creators_genealogy_20261107/classgen.c && ./classgen 8 cls`.
2. Compute the port statistics:
   `for n in 1 2 3 4 5 6 7 8; do ./port_stats cls_$n.txt $n | sed "s/^/$n /"; done > ports_le8_n.txt`.
3. Run the coverage: `./fp_cover ports_le8_n.txt 3 10000000 2 1 7`.
