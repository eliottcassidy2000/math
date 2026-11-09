# Ternary Chan-Pak lanes, and lines between packing-number classes of the tiling metagraph

opus-2026-11-09-S87. Canon: [THM-4950](../../../01-canon/theorems/THM-4950-ternary-chan-pak-lane-reduction-and-order-ten-census.md),
[THM-4951](../../../01-canon/theorems/THM-4951-packing-number-lines-tight-boundary-staircases.md),
[HYP-9540](../../hypotheses/HYP-9540-packing-number-line-laws-tight-words-allowed-pairs.md); updates
[HYP-9530](../../hypotheses/HYP-9530-ternary-chan-pak-and-finite-gaps-in-every-base.md) (A).

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
- This is the base-3 analogue of THM-4094's binary lanes. It has more lanes, because the monoid of residues 1 (mod 3)
  is not factorial.

**Coverage (THM-4950 (3), (4); pre-registered Q1-Q3 CONFIRMED).**
- The complete census of the strong tournaments of order 10 (`ternary_order10.c`) gives 10,070 distinct H_3 values,
  the largest 31,549.
- Each value has a stored witness, re-checked by independent Python (`verify_witness10.py`: 10,070 / 10,070).
- With S86's pieces, every v = 1 (mod 3) up to 29,035 except 10 and 31 is a base-3 packing count. The previous range
  was 5,398.
- The construction threshold is m* = 9679, at the value 29,038 = 2 * 14,519, a lane-Q irreducible.
- **Strong spectrum (THM-4950 (5)).** Below 28,981 the values v = 1 (mod 3) with no strong realiser are exactly 10,
  31 and 64. This is exact below 225, by H <= H_3 and Busch's floor f(11) = 225, and constructive above. It confirms
  the concurrent lane's HYP-9541 (4) on that range; mac-mini's THM-4952 (E) proved the 64.
- Ternary Chan-Pak stays OPEN. Its infinite lanes P and Q are as hard as the binary lanes.

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
| census | lines between nu-classes, black and blue, n <= 9 (`nu_pairs.out`) | FINITE-EXACT |
| (1) | on a tight line (nu(t) + nu(t + 1) = floor(n/3)) the packing numbers are read off the triple word: nu = P(A) | PROVED |
| (2a) | prefix staircase: every split {j, k - j} of the boundary occurs, with tau = nu at both ends | PROVED |
| (2b) | double staircase: blue tight lines with pair {2j, k - 2j}, for every n >= 10 | PROVED |
| (3) | at n = 3k + 2 the tight words are exactly the 2^k block words, each realised by a multi-staircase; blue = palindromic block words, so blue pairs are {j, k - j} with j even or k odd | PROVED |
| counts | tight colourings 2^(k-1), (k+1) 2^(k-2), (k^2 + 9k - 2) 2^(k-4) at n = 3k + 2, 3k + 1, 3k | PROVED / DP-fitted to k = 13 |
| tight lines | 1, 1, 1, 57, 31, 13, 3673, 1802, 541, 563,609 for n = 3..12 (blue 1, 1, 1, 1, 1, 1, 39, 10, 19, 297) | FINITE-EXACT |
| mod 4 | blue lines (H mod 4 pairs) at odd n: (1,1) = (3,3) always; 1:2:1 at n = 5, 7, 9 but 4,014,080 : 8,749,056 : 4,014,080 at n = 11 | KNOWN mechanism (THM-4922 (4)) / FINITE-EXACT |

**Conjectures (HYP-9540).**
- (A) Word law at all n. Proved at n = 3k + 2; verified n <= 12, and at n = 13 by first-hit search.
- (B) Allowed pairs:
  - black lines join nu-classes i <= j iff i >= 1 and i + j >= k;
  - blue lines join the same pairs except (1, 1) from n = 6 on, plus the transitive line (0, k).
- (C) Blue word law: blue tight words = palindromic tight colourings.
- (D) The blue-line (1,3) excess persists at odd n >= 11.

**Pre-registrations.**
- [PREREGISTRATION.md](PREREGISTRATION.md): Q1-Q3, ternary coverage. All CONFIRMED.
- [PREREGISTRATION-2.md](PREREGISTRATION-2.md): R1-R3, the blue-line mod-4 law at n = 11. R1 CONFIRMED, R2 and R3
  REFUTED.
- [PREREGISTRATION-3.md](PREREGISTRATION-3.md): T1-T5, tight lines at n = 12-14. All CONFIRMED.

**Lesson recorded.**
- The mod-4 mechanism (mirror-path decomposition, quotient tournament, Redei) was re-derived here before its source
  was found: the audit proof of THM-4922 (4). It is credited, not claimed.
- The packing bound was likewise Session26's.
- Both were found by the background survey before the canon was written.

## Prior art and concurrent work

**Credited.**
- Session26 (`session26_random_bridges_20261022/session26-blue/RESULT.md`): the line bound and window certificates.
- THM-554: consecutive triples.
- THM-643, THM-646, THM-4729: the tiling frame, closed ends.
- THM-791: blue lines at even n have H = 1 (mod 4).
- THM-4922 (4) and its audit: the odd-order mod-4 balance.
- THM-4094: binary lanes.
- THM-4940 and HYP-9530: base-b packing counts.
- kind-pasteur's `classgen.c`.

**Concurrent.**
- mac-mini-2026-11-08-ternary (HYP-9541, THM-4952/4953) works the same directive with different node classes: number
  of strong components, eternal colour, first creator, self-converse. It also works span-w ribbon words in base 3.
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
| `staircase_check.py`, `staircase_check.out` | the staircase constructions of THM-4951 (2), (3) and the tight-colouring DP |
| `blue_lines.c`, `blue_lines.out` | every blue line by (H mod 8) pair and nu-pair, n <= 9 |
| `mirror_paths.py`, `mirror_check.c`, `mirror_check_n11.out` | the mod-4 decomposition and the exact n = 11 blue-line split |
| `quotient_tournament_check.py` | re-check of the quotient-tournament parity (KNOWN, THM-4922 (4)) |
| `blue_flip_probe.py`, `line_census.c`, `line_classes.py` | exploratory tools: blue orbit flips, merged line census by class, aggregation by class invariants |

**Reproduce.**
- Order-10 census: compile `ternary_order10.c` with `gcc -O3`. Run it on cls_9.txt from kind-pasteur's `classgen.c`
  (191,536 classes) as `ternary_order10 cls_9.txt s 12 t10` for s = 0..11. That takes about 20 minutes on 12 cores.
- Coverage: `python ternary_coverage10.py CENSUS_DIR few_cycle_strong_80.txt t10 40000`, where CENSUS_DIR holds S86's
  alpha_n.txt.
- Tight lines: `tight_lines n` for n <= 11 runs in seconds. For n = 12, use `tight_lines 12 s 12 8` for s = 0..11
  (about 25 minutes).
