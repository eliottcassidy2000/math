# Ternary Chan-Pak lanes, and lines between packing-number classes of the tiling metagraph

opus-2026-11-09-S87. Canon: [THM-4950](../../../01-canon/theorems/THM-4950-ternary-chan-pak-lane-reduction-and-order-ten-census.md),
[THM-4951](../../../01-canon/theorems/THM-4951-packing-number-lines-tight-boundary-staircases.md),
[HYP-9540](../../hypotheses/HYP-9540-packing-number-line-laws-tight-words-allowed-pairs.md); updates
[HYP-9530](../../hypotheses/HYP-9530-ternary-chan-pak-and-finite-gaps-in-every-base.md) (A). Independent audit:
[audit/AUDIT.md](audit/AUDIT.md), PASS WITH CORRECTIONS, all applied (MISTAKE-742).

**Owner directive (2026-11-08).** "Start another session proving ternary Chan-Pak and other conjectures as you
propose. Think about the metagraph and black and blue lines and conjecture about the number of blue or black lines
that can be between new classifications of the nodes you can describe."

## Part I - ternary Chan-Pak (base-3 packing counts)

**Reduction (THM-4950 (1), (2), PROVED).**
- H_3(T) = I(Omega(T), 3). It is multiplicative over order joins and never takes the values 10 or 31 (THM-4940).
- So HYP-9530 (A) is equivalent to: every irreducible of M3 = {v = 1 (mod 3)} minus {10, 31} is H_3 of a strong
  tournament.
- The irreducibles are:
  - the primes = 1 (mod 3) other than 31;
  - the products of two primes = 2 (mod 3) other than 10;
  - the products 10p, 31p and 31pq;
  - five sporadic values: 40, 250, 310, 9610, 29791.
- This is the base-3 analogue of THM-4094's binary lanes: p != 7, 7p with p not 3 or 7, and the carry atoms 49, 63,
  343. It has more infinite lanes because the monoid of residues 1 (mod 3) is not factorial.

**Coverage (THM-4950 (3), (4), FINITE-EXACT; pre-registered Q1-Q3 CONFIRMED).**
- The complete census of the strong tournaments of order 10 (`ternary_order10.c`) gives 10,070 distinct H_3 values,
  from 199 to 31,549.
- Each value has a stored witness, re-checked by independent Python (`verify_witness10.py`: 10,070 / 10,070). The
  audit re-derived the whole set without canonical forms.
- Every v = 1 (mod 3) up to 29,035 except 10 and 31 is an order join of strong tournaments of order <= 10. The
  previous range was 5,398.
- The threshold is m* = 9679, at the value 29,038 = 2 * 14,519, a lane-Q irreducible.

**Strong spectrum (THM-4950 (5)).**
- Below 28,981 the values v = 1 (mod 3), v > 1, with no strong realiser are exactly 10, 31 and 64.
- This is exact below 337: H_3 - 1 >= (3/2)(H - 1), and Busch's floor is f(11) = 225. It is constructive above.
- The method is the concurrent lane's THM-4952 (E). The result confirms the strong-spectrum sentence of their
  HYP-9541 (4) on that range.

Ternary Chan-Pak stays OPEN. Its infinite lanes P and Q are as hard as the binary lanes.

## Part II - lines between packing-number classes

**The new node classification.** nu(T) is the packing number: the maximum number of vertex-disjoint cycles, which
equals the maximum number of disjoint 3-cycles.
- It is an isomorphism and converse invariant, so it classifies metagraph nodes.
- Session26 proved the line bound tau(t) + tau(t + 1) >= floor(n/3) for the feedback number tau. Its
  consecutive-triples proof gives the same bound for nu, and its window certificates nu_c are the P(A), P(B) used
  here.
- Session26 did not study the equality case or count lines between classes. That is the content of THM-4951.

**Results.**

| item | statement | status |
|---|---|---|
| census | lines between nu-classes, black and blue, n <= 9; blue n <= 12 by the audit | FINITE-EXACT |
| (1) | on a tight line (nu(t) + nu(t + 1) = floor(n/3)) the packing numbers are read off the triple word: nu = P(A) | PROVED |
| (2a) | prefix staircase: every split {j, k - j} of the boundary occurs, with tau = nu at both ends | PROVED |
| (2b) | double staircase: a blue tight line for every n >= 10, and the pair {2j, k - 2j} in blue for every n >= 6j + 4 | PROVED |
| (3a, b, d) | at n = 3k + 2 the tight words are exactly the 2^k block words, each realised by a multi-staircase; blue = palindromic block words, so blue pairs are {j, k - j} with j even or k odd | PROVED |
| (3e) | at n = 3k + 1 the tight words are the block words with one block of 2, (k + 1) 2^(k-2) of them up to complement | PROVED |
| (3f) | run-window theorem: every word whose runs all have length >= 2 is realised with nu = (P(A), P(B)); so the word law and the blue word law hold at n != 0 (mod 3) | PROVED (follow-up audit) |
| counts | tight colourings (k^2 + 9k - 2) 2^(k-4) at n = 3k | DP-fitted to k = 13 |
| tight lines | 1, 1, 1, 57, 31, 13, 3673, 1802, 541, 563,609 for n = 3..12 (blue 1, 1, 1, 1, 1, 1, 39, 10, 19, 297) | FINITE-EXACT |
| mod 4 | blue lines (H mod 4 pairs) at odd n: (1,1) = (3,3) always; 1:2:1 at n = 5, 7, 9; 4,014,080 : 8,749,056 : 4,014,080 at n = 11; (1,3)-share minus 1/2 = -0.0009, +0.0017, -0.0001 at n = 13, 15, 17 | KNOWN mechanism (THM-4922 (4)) / FINITE-EXACT / SAMPLED |

**Conjectures (HYP-9540).**
- (A) Word law at n = 3k. It is proved at n != 0 (mod 3) and verified for n <= 12.
- (B) Allowed pairs:
  - black lines join nu-classes i <= j iff i >= 1 and i + j >= k;
  - blue lines join them iff i >= 1 and i + j > k, or {i, j} is the P-pair of a palindromic tight colouring (for
    n <= 26: j even or k odd).
  The first blue clause, "black pairs minus (1, 1)", is refuted at n = 12 (audit, M1).
- (C) Blue word law at n = 3k. It is proved at n != 0 (mod 3).
- (D) Sign law: the blue (1,3) excess has sign (-1)^(c+1), n = 2c + 1 >= 11. The first version ("the excess
  persists") is refuted at n = 13 (U1).

**Pre-registrations.**
- [PREREGISTRATION.md](PREREGISTRATION.md): Q1-Q3, ternary coverage. All CONFIRMED.
- [PREREGISTRATION-2.md](PREREGISTRATION-2.md): R1-R3, the blue-line mod-4 law at n = 11. R1 CONFIRMED, R2 and R3
  REFUTED.
- [PREREGISTRATION-3.md](PREREGISTRATION-3.md): T1-T5, tight lines at n = 12-14. All CONFIRMED.
- [PREREGISTRATION-4.md](PREREGISTRATION-4.md): U1-U4, the (1,3) excess at n = 13, 15, 17. U1 REFUTED; U2-U4
  CONFIRMED.

**Lessons (MISTAKE-742).**
- The mod-4 mechanism (mirror-path decomposition, quotient tournament, Redei) was re-derived here before its source
  was found: the audit proof of THM-4922 (4). It is credited, not claimed. The packing bound was likewise Session26's.
- The first blue clause of HYP-9540 (B) contradicted THM-4951 (3d) and T4, both proved or checked in this session.
- The first (D) stayed in the canon text after U1 refuted it.

## Prior art and concurrent work

**Credited.**
- Session26 (`session26_random_bridges_20261022/session26-blue/RESULT.md`): the line bound, the window certificates
  and the independence corollary.
- THM-554: consecutive triples.
- THM-4731 (6) and THM-337: the full staircase T*_n and its cycle structure.
- THM-643, THM-646, THM-4729: the tiling frame, closed ends.
- THM-791: blue lines at even n have H = 1 (mod 4), and the odd-n contrast.
- THM-4922 (4) and its audit: the odd-order mod-4 balance and (P2).
- THM-4094: binary lanes.
- THM-4940 and HYP-9530: base-b packing counts.
- THM-4744: Busch's formula.
- kind-pasteur's `classgen.c`.

**Concurrent.** mac-mini-2026-11-08-ternary (THM-4952/4953, HYP-9541) works the same directive with different node
classes: number of strong components, eternal colour, first creator, self-converse.
- It also works span-w ribbon words in base 3.
- THM-4952 (E) supplied the method of THM-4950 (5).
- THM-4952 (C) also records the staircase's interval cycles.
- Session229's "ternary closure" is a construction grammar for binary Chan-Pak, unrelated to base 3.

## Files

| file | content |
|---|---|
| `ternary_order10.c` | strong order-10 census of H_3 values; witness mode |
| `verify_witness10.py`, `verify_witness10.out`, `witnesses_order10.txt` | independent re-check of the 10,070 order-10 witnesses (value, canonical code) |
| `ternary_coverage10.py`, `coverage10.out` | multiplicative closure, m*, Q1-Q3, missing irreducibles by lane |
| `ternary_lanes.py`, `ternary_lanes.out` | irreducibles of M3 and the shape theorem below 200,000 |
| `strong_spectrum10.py`, `strong_spectrum10.out` | THM-4950 (5): the strong base-3 spectrum below 28,981 |
| `nu_lines.c`, `nu_pairs.out` | every line by nu-pair and colour, n <= 9 |
| `tight_lines.c`, `tight_lines.out`, `t5_first_hits.out` | tight lines n <= 12 (pairs, colours, triple-word histogram); first-hit realizations at n = 13, 14 |
| `staircase_check.py`, `staircase_check.out` | the staircase constructions of THM-4951 (2), (3b) and the tight-colouring DP |
| `run_window_check.py`, `run_window_check.out` | run-window realisation of all tight colourings, n <= 16 |
| `run_window_theorem_check.py`, `run_window_theorem_check.out` | THM-4951 (3e), (3f) by brute force |
| `blue_lines.c`, `blue_lines.out` | every blue line by (H mod 8) pair and nu-pair, n <= 9 |
| `mirror_paths.py`, `mirror_check.c`, `mirror_check_n11.out` | the mod-4 decomposition and the exact n = 11 blue-line split |
| `mirror_sample.c`, `mirror_sample.out` | sampled (1,3)-shares at n = 13, 15, 17 (PREREGISTRATION-4) |
| `quotient_tournament_check.py` | re-check of the quotient-tournament parity (KNOWN, THM-4922 (4)) |
| `blue_flip_probe.py`, `line_census.c`, `line_classes.py` | exploratory tools: blue orbit flips, merged line census by class, aggregation by class invariants |
| `audit/` | the independent audit: report and all its own programs and outputs |

**Reproduce.**
- Order-10 census: compile `ternary_order10.c` with `gcc -O3`. Run it on cls_9.txt from kind-pasteur's `classgen.c`
  (191,536 classes) as `ternary_order10 cls_9.txt s 12 t10` for s = 0..11. That takes about 20 minutes on 12 cores.
- Coverage: `python ternary_coverage10.py CENSUS_DIR few_cycle_strong_80.txt t10 40000`, where CENSUS_DIR holds S86's
  alpha_n.txt.
- Tight lines: `tight_lines n` for n <= 11 runs in seconds. For n = 12, use `tight_lines 12 s 12 8` for s = 0..11
  (about 25 minutes).
- Sampler: `mirror_sample n 1 SAMPLES SEED s 12`.
