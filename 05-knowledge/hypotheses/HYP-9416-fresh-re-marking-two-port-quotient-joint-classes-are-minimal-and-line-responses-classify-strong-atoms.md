---
id: HYP-9416
title: "Fresh re-marking two-port quotient. Rooted two-port states [C; s, t] evolve by fresh (unnamed, multiset-counted) re-marking along Hamilton paths from s to t followed by complement (V1), or along all Hamilton paths with ports reset (V2). (a) PROVED (audit reconstruction lemma): already the single sandwich e*h*e determines every strong joint two-port class for n >= 3, so the static observer of the two-sided gluings k*h*k' with contexts of order <= 2 trivially admits no closed compression on strong states (scalar two-sided observers such as H(e*h*e) do compress). (b) One-sided gluings depend only on the one-sided classes (PROVED) and never separate the one-sided collision pairs; one V1 step, whose successor multiset is a joint-port statistic, separates all of them even from the bare class observation (FINITE-EXACT n <= 7); V2 separates none. (c) The port-forgotten line-response histogram of a strong atom determines the atom (FINITE-EXACT n <= 8; even the set of partner classes does). (d) If only the unrooted class must be determined at every step, the V1 quotient is a proper compression (FINITE-EXACT n = 4..7); the V2 quotient is the class partition (PROVED, trivial). This is an analogue, on two-port classes and with counted multisets, of the stochastic quotient mentioned in vitali-response C5, not that quotient."
status: >
  (1) PROVED for every n >= 3 by the independent audit (reconstruction lemma: e*h*e alone determines a strong two-port state);
  the one-sidedness of glue1 and the V2/free statement are PROVED (trivial). (3) is FINITE-EXACT through n <= 8 (6008 strong atoms); all-order OPEN. (2) closure part and the V1 part of (4):
  FINITE-EXACT, exhaustively over all tournaments, for n <= 7 (17824 rooted states and 14162 strong at n = 7), all numbers
  independently reproduced; all-order versions OPEN. Validated against the published 289/253 hostile. Independently audited
  2026-10-26: 04-computation/experiments/joint_port_bisimulation_20261026/audit_hyp9416/.
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
- **Quotient.** For an observation `O`, take the coarsest partition refining `O` in which equivalent states have equal multisets of successor cells. It has the form of the operation-kernel quotient of vitali-response C5, but on a different carrier (two-port classes of all tournaments, not the presentations X_n), with V1 restricted to the current ports and with successor multisets counted rather than normalized: the normalized-distribution quotient is slightly coarser (V1 refining `free`: 11031 instead of 11033 strong cells at n = 7; refining `sand`: 5888 instead of 5889).
- **Observations used:**
  - `free`: the unrooted class;
  - `sand`: `H(e*h*e)`;
  - `glue1`: the free classes of the one-sided gluings `h*k` and `k*h`, over all rooted `k` of order `≤ 3`;
  - `glue2`: the free classes of the two-sided gluings `k*h*k'`, over all rooted `k, k'` of order `≤ 2`.

## Statement

1. **Joint classes are minimal for the two-sided gluing observer: PROVED, trivially.** For every `n ≥ 3` the single two-sided gluing `e*h*e` determines a strong two-port state. In `e*h*e` the added entrance vertex `a` has out-degree 1 and the added exit vertex `b` has in-degree 1. Since `C` is strong, the only other vertex of out-degree 1 can be `s` (when `d⁺_C(s) = 1`), whose out-neighbour has out-degree at least 2, and the only other vertex of in-degree 1 can be `t`, whose in-neighbour has in-degree at least 2. As `a → s` and `t → b`, the vertices `a` and `b` are determined, and `s = N⁺(a)`, `t = N⁻(b)`, `C = e*h*e − {a, b}` (audit `sandwich_reconstruction.py`). Hence `glue2` separates any two distinct strong joint two-port classes in every order (through `n = 7` each of its nine contexts does so alone), and no closed compression refining it exists. This concerns this observer only: the two-sided scalar observer `sand = H(e*h*e)` does compress (5889 strong V1 cells of 14162 at `n = 7`).
2. **Fixed-port re-marking separates one-sided collisions.**
   - A one-sided collision is a pair of distinct strong joint classes with equal entrance-rooted and equal exit-rooted classes.
   - `glue1` depends only on one-sided classes (PROVED: an isomorphism of `(C, t)` extends by the identity to `h*k`, one of `(C, s)` to `k*h`), so it separates none of these pairs.
   - The V1 quotient refining `glue1` separates all of them (FINITE-EXACT, `n ≤ 7`). This is not specific to `glue1`: the V1 quotient refining the bare class observation `free` separates every pair at its first refinement round, i.e. by the multiset of unrooted partner classes over the Hamilton paths from `s` to `t`, a joint-port statistic by construction. The number of `s → t` Hamilton paths alone (the size of that multiset) separates 1, 4, 16, 65 of the 1, 7, 20, 70 pairs at `n = 3, 5, 6, 7`. The V2 quotient separates none of them. One-sided observers that see only one side do not suffice: V1 refining only the free class of `h*e` leaves 12371 strong cells of 14162 at `n = 7`.
3. **Line responses classify strong atoms.** Let `F_A` be the joint response histogram of rooted-lines R1. Distinct strong atoms have distinct histograms, and this survives forgetting the ports, i.e. keeping only the multiset of target free classes.
4. **Compression for the class observable.** The V1 quotient refining `free` is a proper compression on strong states for `n ≥ 4` (5 strong cells of 12 at `n = 4`; FINITE-EXACT). The V2 quotient refining `free` is the class partition itself (PROVED, trivially: V2 successor multisets do not depend on the ports).

## Evidence (FINITE-EXACT, all tournaments)

| n | rooted states | strong | V1/free (strong cells) | V1/sand | V1/glue1 (strong) | glue2 static (strong) | V2/free |
|---|---|---|---|---|---|---|---|
| 3 | 8 | 2 | 4 (2) | 5 | 8 (2) | 2 | 2 |
| 4 | 32 | 12 | 11 (5) | 14 | 30 (12) | 12 | 4 |
| 5 | 176 | 92 | 60 (42) | 64 | 173 (92) | 92 | 12 |
| 6 | 1408 | 946 | 691 (607) | 445 | 1394 (946) | 946 | 56 |
| 7 | 17824 | 14162 | 11605 (11033) | 6165 | 17767 (14162) | 14162 | 456 |

**One-sided collisions** (strong states):
- at `n = 3, 4, 5, 6, 7` there are 1, 0, 7, 20, 70 pairs, in 1, 0, 2, 12, 39 fibres (the `n = 3` pair is the cyclic triangle with `s → t` versus `t → s`);
- `glue1` separates 0 of them statically and all of them after closure;
- `glue2` separates all of them statically;
- the singleton sandwich `H(e*h*e)` alone separates 6, 17, 64;
- the number of `s → t` Hamilton paths alone separates 4, 16, 65 (together with `H(e*h*e)`: 7, 17, 68), and every pair splits at the first V1 refinement round of `free` (audit).

**Strong atoms** at `n = 3..7`:
- there are 1, 1, 6, 35, 353 atoms, with 1, 5, 44, 547, 9985 distinct joint targets (these counts, the blue counts and every coefficient of `F_A` and `B_A` are already in vitali-response `strong-atlas.json`; the audit matched all 396 sources coefficient by coefficient);
- the histograms `F_A` are pairwise distinct in every order (vacuously at `n = 3, 4`), also after forgetting the ports; even the set of partner free classes, without multiplicities, separates all atoms through `n = 7`, while the multiset of partner score sequences gives 351 classes for the 353 atoms at `n = 7`;
- 1, 1, 6, 9, 79 of the atoms have blue presentations.

**Order-eight extension (FINITE-EXACT, 2026-10-26):** the [unnamed-carrier packet](../results/unnamed_carrier_dissection_20261026/README.md) exhausts all 2,097,152 fixed-path presentations. Its 6,008 strong classes have 2,032,504 presentations, 269,860 joint targets, and 162 atoms with blue presentations. All 6,008 joint and free response histograms are distinct; both unweighted support signatures also remain distinct. The C++ extension independently matches the preceding counts through seven. No order-eight V1 closure is computed. The original `atom_responses.py` dictionary-overwrite pushforward is repaired there; the independent audit already used the correct sum.

**Validation.** The rooted-lines R3 targets `C = (26,12,41,48,38,3)` and `D = (38,20,56,19,33,10)` have equal entrance-rooted and exit-rooted classes but different joint classes, and give `H(e*h*e) = 289` and `253`; they also have 4 and 2 Hamilton paths from entrance to exit, so their V1 successor multisets already differ in size.

## Relation to prior results

- **vitali-response C5** (named guarded re-marking): legal-guard signatures alone separate every presentation through `n = 7`. That packet mentions a quotient under distributions of fresh uniform re-marking on its presentation carrier as a different, possibly smaller, question that must not stand in for the named protocol; the quotients here are a related construction on two-port classes (Session70 likewise records them as a different operation menu).
- **Answer here.**
  - For the observer `glue2` there is no smaller quotient on strong states (1), trivially, since `e*h*e` alone is a complete invariant of strong two-port states; the two-sided scalar observer `sand` compresses (5889 strong V1 cells of 14162 at `n = 7`).
  - For the bare class observable there is a smaller one (4).
  - V1 (fixed-port) successor multisets are joint-port statistics, so one V1 step separates every one-sided collision pair (2); V2 adds nothing to one-sided data.
- **81/27 control** (`control_81_27.out`):
  - both 12-vertex tournaments have `H = 462915` and 67 cyclic triples;
  - the partner map is injective on presentations in both (5715 and 17145 partner classes);
  - only 74 partner classes are shared (already recorded in HYP-9414 (3), 74 between `T_000` and `T_100`); all 74 come from the phase-0 entrance family of the period-3 decoration;
  - the period-3 decoration's three entrance phases carry equal path mass, 154305 each. This is forced: A- = A+^op has the same ordered path-cover counts as A+ (fibre F4 of the tournament-composition packet), the T+ blocks carry the same masses, and the joint (entrance, exit) phase mass depends only on the phase difference (42930, 45036, 66339).
  - This re-checks HYP-9414 (1) and (3) at `m = 3` and the control of the vitali-response packet. No claim beyond those is made.

**Separate literal protocol, now PROVED:** for presentations X_n with complement as its own operation and fresh reset, the minimal exactly lumpable quotient preserving q is the ordered pair `(q(x),q(tau x))`, when the whole reset row depends only on q. The [U2 proof and hostile controls](../results/unnamed_carrier_dissection_20261026/README.md) cover free and joint observations. Blue joins distinguish states identified by that pair, so it is not a minimum for the full attachment grammar. This does not replace the arbitrary-port V1 state space above.

## Open

- All-order proofs of (2) (closure part) and (3); (1) is proved (audit lemma).
- Answered: order-1 contexts suffice, by the degree reconstruction of `e*h*e`. Open: whether `A2(s→t)*h*A2(s→t)` (`k = k'` the rooted arc entrance → exit), which is injective on all 8, 32, 176, 1408, 17824 states at `n = 3..7`, is a complete invariant of all two-port states in every order.
- Whether the class-observable compression ratio (strong cells 0.46, 0.64, 0.78 at `n = 5, 6, 7`) tends to 1.
- The same questions for multi-run (arity `k > 1`) attachments.
