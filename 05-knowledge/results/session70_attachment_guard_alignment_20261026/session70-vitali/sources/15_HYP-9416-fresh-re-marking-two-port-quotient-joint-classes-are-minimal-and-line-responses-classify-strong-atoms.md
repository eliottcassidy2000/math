---
id: HYP-9416
title: "Fresh re-marking two-port quotient. Rooted two-port states [C; s, t] evolve by fresh (unnamed, multiset) re-marking along Hamilton paths from s to t followed by complement. (a) The free classes of two-sided gluings k*h*k' with rooted contexts of order <= 2 already separate all strong joint two-port classes, so no closed compression of the joint class exists for strong states. (b) One-sided gluings never separate the one-sided collision pairs, but closure under the dynamics separates all of them. (c) The port-forgotten line-response histogram of a strong atom determines the atom. (d) If only the unrooted class must be determined at every step, the fresh-re-marking quotient is a proper compression. This is the stochastic quotient left open by the vitali-response packet's C5."
status: >
  OPEN for all orders. FINITE-EXACT, exhaustively over all tournaments, for n <= 7 (17824 rooted states and 14162 strong
  at n = 7). Validated against the published 289/253 hostile. Not yet independently audited.
session: mac-mini-2026-10-25-pyth (2026-10-26)
source: 05-knowledge/results/fresh_remarking_quotient_and_divisor_arrangements_20261026.md
scripts:
  - 04-computation/experiments/joint_port_bisimulation_20261026/rooted.py, canon_nauty.py: rooted canonical forms (exact over the n-2 unmarked vertices) and nauty labelg free forms
  - 04-computation/experiments/joint_port_bisimulation_20261026/bisim3.py (+ bisim3_*.out, bisim3_n*.json): quotients for both dynamics and four observations, n = 3..7
  - 04-computation/experiments/joint_port_bisimulation_20261026/bisim_detail.py, bisim_detail2.py (+ .out): static versus closed separation, one-sided collisions, two-sided contexts
  - 04-computation/experiments/joint_port_bisimulation_20261026/atom_responses.py (+ .out): strong-atom response histograms, n = 3..7
  - 04-computation/experiments/joint_port_bisimulation_20261026/control_81_27.py (+ .out): the 81/27 control with entrance phases
  - 04-computation/experiments/joint_port_bisimulation_20261026/validate_289_253.out: the sandwich code reproduces H = 289, 253
related:
  - 05-knowledge/results/strong_atom_continuation_vitali_20261026/README.md   # C1-C5: named, guarded re-marking; this HYP is its fresh-re-marking companion
  - 05-knowledge/results/rooted_lines_frobenius_sheets_20261026/README.md    # R1-R5 joint two-port response, 289/253
  - 05-knowledge/results/session69_prime_ports_power_aliquot_20261026/README.md
  - 05-knowledge/hypotheses/HYP-9414-module-shattering-consequences-injective-fibres-pair-module-partners-and-the-orbit-size-law.md
---

# HYP-9416 — The fresh-re-marking two-port quotient

## Setting

Conventions follow the rooted-lines note.
- **States.** A two-port state `h = [C; s, t]` is an isomorphism class of a tournament with an ordered, jointly marked entrance `s` and exit `t`, `s ≠ t`.
- **Fixed-port dynamics V1 (fresh re-marking).** The successor multiset of `h` is `{[C^op xor E(P); s, t] : P a Hamilton path of C from s to t}`. Each path is fresh: it is unnamed, and only the multiset is kept.
- **Free re-marking V2.** The successors range over all Hamilton paths `P`, with ports reset to the ends of `P`.
- **Quotient.** For an observation `O`, take the coarsest partition refining `O` in which equivalent states have equal multisets of successor cells. This is the bisimulation (operation-kernel) quotient of the vitali-response C5, with multisets in place of named guarded operations.
- **Observations used:**
  - `free`: the unrooted class;
  - `sand`: `H(e*h*e)`;
  - `glue1`: the free classes of the one-sided gluings `h*k` and `k*h`, over all rooted `k` of order `≤ 3`;
  - `glue2`: the free classes of the two-sided gluings `k*h*k'`, over all rooted `k, k'` of order `≤ 2`.

## Statement

1. **Joint classes are minimal (two-sided).** For every `n ≥ 5`, the static observation `glue2` separates any two distinct strong joint two-port classes. So no closed compression of the joint class exists for strong states, even without the dynamics.
2. **Closure upgrades one-sided data.**
   - A one-sided collision is a pair of distinct strong joint classes with equal entrance-rooted and equal exit-rooted classes.
   - `glue1` depends only on one-sided classes, so it separates none of these pairs.
   - The V1 quotient refining `glue1` separates all of them.
3. **Line responses classify strong atoms.** Let `F_A` be the joint response histogram of rooted-lines R1. Distinct strong atoms have distinct histograms, and this survives forgetting the ports, i.e. keeping only the multiset of target free classes.
4. **Compression for the class observable.** The V1 quotient refining `free` is a proper compression on strong states for `n ≥ 5`. The V2 quotient refining `free` is the class partition itself.

## Evidence (FINITE-EXACT, all tournaments)

| n | rooted states | strong | V1/free (strong cells) | V1/sand | V1/glue1 (strong) | glue2 static (strong) | V2/free |
|---|---|---|---|---|---|---|---|
| 3 | 8 | 2 | 4 (2) | 5 | 8 (2) | – | 2 |
| 4 | 32 | 12 | 11 (5) | 14 | 30 (12) | – | 4 |
| 5 | 176 | 92 | 60 (42) | 64 | 173 (92) | 92 | 12 |
| 6 | 1408 | 946 | 691 (607) | 445 | 1394 (946) | 946 | 56 |
| 7 | 17824 | 14162 | 11605 (11033) | 6165 | 17767 (14162) | 14162 | 456 |

**One-sided collisions** (strong states):
- at `n = 5, 6, 7` there are 7, 20, 70 pairs, in 2, 12, 39 fibres;
- `glue1` separates 0 of them statically and all of them after closure;
- `glue2` separates all of them statically;
- the singleton sandwich `H(e*h*e)` alone separates 6, 17, 64.

**Strong atoms** at `n = 3..7`:
- there are 1, 1, 6, 35, 353 atoms, with 1, 5, 44, 547, 9985 distinct joint targets;
- the histograms `F_A` are pairwise distinct in every order, also after forgetting the ports;
- 1, 1, 6, 9, 79 of the atoms have blue presentations.

**Validation.** The rooted-lines R3 targets `C = (26,12,41,48,38,3)` and `D = (38,20,56,19,33,10)` have equal entrance-rooted and exit-rooted classes but different joint classes, and give `H(e*h*e) = 289` and `253`.

## Relation to prior results

- **vitali-response C5** (named guarded re-marking): legal-guard signatures alone separate every presentation through `n = 7`. That packet flags the fresh-re-marking quotient as a different, possibly smaller, question.
- **Answer here.**
  - For two-sided attachment observations there is no smaller quotient on strong states (1).
  - For the bare class observable there is a smaller one (4).
  - Closure transfers joint-port information into one-sided observations (2).
- **81/27 control** (`control_81_27.out`):
  - both 12-vertex tournaments have `H = 462915` and 67 cyclic triples;
  - the partner map is injective on presentations in both (5715 and 17145 partner classes);
  - only 74 partner classes are shared;
  - the period-3 decoration's three entrance phases carry equal path mass, 154305 each.
  - This confirms HYP-9414 (1) at `m = 3` and the control of the vitali-response packet. No claim beyond those is made.

## Open

- All-order proofs of (1)–(3).
- A structural reason why order-2 two-sided contexts suffice.
- Whether the class-observable compression ratio (strong cells 0.46, 0.64, 0.78 at `n = 5, 6, 7`) tends to 1.
- The same questions for multi-run (arity `k > 1`) attachments.
