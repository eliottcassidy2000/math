# Audit P1: mac-mini-2026-11-09-paleyeven

VERDICT: PASS WITH CORRECTIONS

Independent adversarial audit, 2026-11-09. I read PROOF.md, README.md, THM-4995, THM-4996, HYP-9575 and
PREREGISTRATION.md in their 10:25 versions, which include the coordinator's two mid-audit additions: the Royle-odd
spectrum {1, 5, 37} and the certified order-17 witness. I also read the session's code and outputs. Every number
was recomputed with programs written for this audit (in `audit/`); none of the session's code was reused.

Overall:
- Every number checked is right.
- Every proof checked is correct. Theorem A2's case analysis is complete, and Lemma A1's ranges are right apart
  from one rounding.
- D(C_9{1,2,3,5}) is a certified escaper of order 18 with decomposability defect 5. Both witnesses were recounted
  by two different algorithms.
- Theorem B4 was reproduced by an independent census of all 474,101 classes, with every per-value count agreeing.
- The corrections are one missing priority (MAJOR: Lemma A1 is a 1974 theorem, and md is the classical
  "simplicity index") and thirteen minor items (Findings 2-14): typing, citations and wording. Finding 13 is
  already resolved.

## Findings

### 1. MAJOR: missing priority for Lemma A1 and for the notion md

**What is already known.**
- md(T) is the **simplicity index** s(T): the least number of arc reversals that make T non-simple, that is, give
  it a module M with 2 <= |M| <= n-1.
- Müller and Pelant (*On strongly homogeneous tournaments*, Czechoslovak Math. J. 24 (1974) 378-391) proved
  s(T) <= (n-1)/2 for every n-tournament, with equality if and only if T is doubly regular. This is cited in
  Boussaïri, Lakhlifi, Talbaoui, *The simplicity index of tournaments*, Contrib. Discrete Math. 17(2) (2022) 172-186,
  arXiv:1907.11777, whose abstract states that "equality holds if and only if T is doubly regular".
- So **Lemma A1 is the equality case of Müller-Pelant**. The S S^T = qI - J / Cauchy-Schwarz argument is a new proof
  of a known theorem.
- BLT 2022, Proposition 3.1, gives s(T) <= min(delta0, sigma), where sigma is the twin distance. That covers the
  md <= delta0 and md <= twin parts of Lemma A0.
- lambda is Kirkland's *reversal index* (Linear Multilinear Algebra 34 (1993) 343-351, as summarised in BLT).
- So the "module defect of Paley-type tournaments" is classical.

**Where the canon presents this without credit.**
- THM-4995 line 3 (title: "The input is md(B) = (q - 1)/2, from S S^T = qI - J by a Parseval bound").
- THM-4995 line 6 and line 64.
- THM-4995 KNOWN list, lines 18-27.
- PROOF.md lines 25-66 (Lemmas A0, A1).
- README lines 24-25.

**What I found no prior art for.** The double construction and its results:
- md(D(B)) = lambda = delta0 = (q+1)/2;
- the conditional escapers;
- the order-18 certificate.

My searches covered simplicity and decomposability indices and doubly regular doubling.

**Fix.**
- Credit Müller-Pelant 1974 for Lemma A1. Present the Parseval proof as an alternative proof.
- Credit BLT 2022 for s(T) = md and for Proposition 3.1.
- Credit Kirkland 1993 for lambda as a reversal index.
- Reword the THM-4995 title as "the input is Müller-Pelant's md(B) = (q - 1)/2 (re-proved here via S S^T = qI - J)".
- THM-4985's "decomposability defect" is the same notion and is also uncredited. I recommend a follow-up correction
  there.
- Note: Belkhechine-Ben Salha's "decomposability index" (arXiv:2003.06503) is the opposite notion (reversals to make
  T *in*decomposable). Do not cite it for md.

### 2. MINOR (typing): T_7 and the R_5 chain are stated as excluded or failed, which is not established

T_7 is an escaper iff 4,730,397 is in S(13) and 686,457 is in S(12). Neither membership is decided:
- max H(12) and max H(13) are unknown;
- THM-329's 531,205 is only "strongly believed" to be a(12).

So the escaper half of Q1 is **undecided** (not certifiable now), not refuted. It is false iff 4,730,397 is not in S(13)
or 686,457 is not in S(12). Sufficient conditions for that are a(12) = 531,205 or max H(13) < 4,730,397.

Places that overstate it:
- README line 115 ("Q1 REFUTED");
- HYP-9575 lines 24-25 ("T_7 is excluded");
- HYP-9575 line 47 ("Smaller candidates fail. T_7 at order 14 fails by its count");
- PROOF.md line 166 ("T_7 is not:", which reads as "T_7 is not an escaper").

Likewise for the order-15 R_5 chain:
- H = 29,678,025 = 2.79 E[H(14)] is not known to lie outside S(14);
- its largest deletion, 5,081,917, exceeds only the best *known* order-13 count.

Fix:
- "Q1: md = lambda = delta0 = 4 CONFIRMED; escaper status UNDECIDED (not certifiable; false if THM-329's a(12) is
  exact)";
- "not certifiable" in HYP-9575 (A) and (C);
- "the argument does not cover T_7" at PROOF line 166.

THM-4995 (2) ("not a certifiable escaper") is correctly worded.

### 3. MINOR: a wrong "iff" (THM-4995 line 29)

"it is one iff max H(13) < 4,730,397" is wrong: only "if" holds. T_7 is a first creator iff 4,730,397 is not in S(13).
S(13) can have gaps below its maximum, as S(11) has 3,635 gaps in [80,453, 95,093] (THM-4889).

Fix: "(it is one iff 4,730,397 is not in S(13), which holds if max H(13) < 4,730,397; max H(13) is unknown)".

### 4. MINOR (proof detail): the large-q step of Theorem A4 has no explicit constant (PROOF lines 160-161; THM-4995 line 13)

The step says "q log(3/4) + O(log q) ... hence negative beyond the verified range" but gives no constant. The
conclusion is true. From Robbins' bounds I derived, for both log-ratios (natural logs) and every q >= 7,

    R(q) <= q log(3/4) + 2.85 log q + 4.8.

- The right side is decreasing for q >= 10 and negative for q >= 57.
- Numerically the inequality holds for all q = 3 mod 4 in [7, 2*10^6], with minimum slack 0.96 at q = 1975
  (`bregman_tail.out`).

Suggested replacement: state this inequality. Together with the exact table for 55 <= q <= 20,000 it gives all
q >= 55.

Also, "exact log-gamma evaluation" (PROOF line 157, THM-4995 line 12) is double-precision lgamma. It agrees with
exact-integer factorial logs to 9e-13 (`bregman_check.out`). Say "lgamma (checked against exact factorials)".

### 5. MINOR: Lemma A1, case s = 3 (PROOF line 58)

"iff q >= 10.5" should read "iff q >= 10.464", so q >= 11. An exact integer check over all q = 3 mod 4 up to 20,000 and
all 3 <= s <= q - 2 finds the bound f(s) >= (q - 1)/2 failing only at q = 7 (s = 3 and s = 5). That case is covered by
the computation, as the proof says (`lemmaA1_bound.out`).

### 6. MINOR (pre-registration): the outcomes section is empty

PREREGISTRATION.md lines 37-39 say outcomes are "appended after the runs", but the section is empty. The outcomes
appear only in README section 4.

Fix: append them there, with Finding 2's correction for Q1, and record the sub-predictions:
- Q2: in-session certification of T_11 (p = 0.35) not achieved;
- Q4: "21 unattained" (p = 0.60) refuted; "some other value <= 255 unattained" (p = 0.55) confirmed, by 13.

Q2, Q3 and Q4 as stated in README match the predictions. Q3 should mention that its input lemma is Müller-Pelant
(Finding 1).

Timing: the pre-registration commit (719e10f41d, 09:14:09) precedes the binaries (09:14:50) and the first census
output (09:25:50). That is consistent with "committed before the runs".

### 7. MINOR (citation): Brouwer's parity theorem

THM-4996 line 19 and PROOF line 194 cite "Brouwer 2009". The standard citation is A. E. Brouwer, P. Csorba, A.
Schrijver, *The number of dominating sets of a finite graph is odd*, preprint (2009), as cited by Heinrich-Tittmann,
arXiv:1701.03453.

### 8. MINOR (attribution): the even-graph notion and the name "Royle-even"

Places: PROOF line 197, README lines 52-53, THM-4996.
- Even graphs, and the conjecture that they are equinumerous with tournaments, are due to Pontus von Brömssen
  (né Andersson); see OEIS A334335 and the RPGFD abstract.
- RPGFD (J. Algebraic Combin. 57 (2023) 515-524, arXiv:2204.01947) proved the conjecture. This citation is correct.
- "Royle-even" is a coinage that credits the wrong person. Suggest "even (von Brömssen; RPGFD)".
- THM-4996 line 13 / PROOF line 258 ("reproduced class by class"): these are counts, and von Brömssen already
  verified them to n = 10 (OEIS comment). Say "count by count, a re-check".

### 9. MINOR (credit): the cone and join formula

d(A v G) = (2^|A| - 1)(2^|G| - 1) + d(A) + d(G) is the x = 1 case of the domination-polynomial join theorem
(Alikhani-Peng, arXiv:0905.2251, Theorem 2; the domination polynomial is from Arocha-Llano 2000).

It is used in Lemma B1 (PROOF lines 204-210, "very possibly folklore") and in dom_cover / B5. Cite it there. The
completeness and conciseness of d on the odd numbers itself: no prior art found.

### 10. MINOR (credit): Brégman for Hamiltonian paths (PROOF lines 136-142; THM-4995 line 24)

Using Brégman's permanent bound to bound Hamiltonian paths and cycles of tournaments is Alon's method (*The maximum
number of Hamiltonian paths in tournaments*, Combinatorica 10 (1990) 319-324). Credit it next to Brégman 1973.

### 11. MINOR (citation detail): Reid-Brown (THM-4995 line 23)

The attribution is correct. Add the details: K. B. Reid, E. Brown, *Doubly regular tournaments are equivalent to skew
Hadamard matrices*, J. Combin. Theory Ser. A 12 (1972) 332-338.

### 12. MINOR (evidence wording): the description of the unreached values

Places: PROOF line 282, README line 74, HYP-9575 lines 4 and 38-39, THM-4996 line 17.

The claim is that the unreached values sit "at 2^m + {3, 7, 13} and just below 2^m" ("near powers of 2"). The printed
list in `dom_cover_K511_22_genjoin.out` shows clusters near 2^m - 2^k as well:
- 761-773 around 768 = 2^10 - 2^8;
- 889 near 896;
- 949-953 near 960;
- 1531 near 1536;
- 1783-1799 near 1792;
- 1913 near 1920;
- 3571-3581 near 3584.

Describe them as "near 2^m - 2^k (including 2^m)". B5 is evidence only, and I did not re-run it.

### 13. MINOR (wording, already resolved)

README's earlier "twin-free realiser" for 21 was wrong: K_1 v (P_3 u K_1) has non-adjacent twins. Line 89 now reads
"a realiser whose only twins are non-adjacent". I confirm the new wording; the census shows |Aut| = 2.

### 14. MINOR (wording): PROOF.md Theorem A5 "Housekeeping" bullet (lines 182-183)

The note about the lost seed-1709 output is process history, not proof. Move it to README section 4 or a log.

## Checked and found correct (no finding)

- **Chan-Pak citation.** Exact, checked against arXiv:2308.10214v2:
  - Section 5.1 (complete, almost complete);
  - Definition 5.7 (concise);
  - Proposition 5.14 (exponential growth plus almost complete implies concise);
  - Conjecture 6.2 (h takes all odd values except 7 and 21);
  - Theoret. Comput. Sci. 1015 (2024) 114776.

  THM-4996 (1) correctly says "complete on the odd numbers" (Chan-Pak's "complete" means all of N).

- **Lemma A0.** Both cut cases are correct.
- **Lemma A1.**
  - The identity S S^T = qI - J and its equivalence with double regularity hold. S^T S = S S^T because S is skew.
  - The Parseval step sum_w D_w^2 = qs - s^2 and the Cauchy-Schwarz bound f(s) are correct.
  - The ranges s = q - 2 (q >= 7.83), s = q - 3 (equality at q = 7) and 4 <= s <= q - 4 (concavity) are correct.
  - The q = 7 base case holds: all 240 labelled doubly regular 7-tournaments (= 5040/|Aut(QR_7)|) have md 3.
- **Theorem A2.**
  - I re-derived every case 1-7 (including (4a) with u != v: d+(u) + d-(v) + 2([v->u] - [u->v])). The cases cover
    every (s, t) with s + t >= 2 and (s, t) != (q, q). Overlaps at (q, 1) and (1, q) are harmless.
  - The upper bound via delta0 = mu + 1 is correct.
  - Strong and prime follow from md >= 1.
  - The diagonal Aut(B) action is correct, and so is the anti-automorphism a_x -> b_phi(x). I also checked it by
    computer for T_7, T_11 and D(C_9).
- **Proposition A3 / Theorem A4.**
  - The reduction H(D) <= per(A(D + w)) and the row sums d+(v) + 1 are correct.
  - The out-degrees (3q-3)/2 on A and (q+1)/2 on B' are correct.
  - The worst deletion is a b-vertex: it removes the smallest factor. I checked that the b-deletion bound is at least
    the a-deletion bound for every q.
  - The asymptotic q log(3/4) is correct (derivation in Finding 4). The conditional logic of A4 follows: Rédei gives
    oddness, H is far from 7 and 21, and Paley gives infinitely many q.
- **Part B.**
  - Lemma B1 holds.
  - All six items of Lemma B2 hold, including:
    - the swap sign (-1)^s;
    - the cone condition (no universal vertex; outside it, every cone over G with a universal vertex is odd: 52 of 52);
    - the join condition (complements connected, A not isomorphic to G);
    - the tree criterion.
  - Lemma B3 holds: the strict monotonicity, the BFS prefix argument and the child formula are all correct.

## Independently reproduced

All programs are in `audit/`, written for this audit. Each final run used under 30 s and well under 1 GB. An early
self-test variant, with random graphs up to n = 14, was stopped after 3 minutes and replaced by the n <= 9 test. The
wit jobs were untouched.

**Part A**

| claim | program | output | result |
|---|---|---|---|
| QR_7, QR_11, T_7, T_11, D(C_9{1,2,3,5}): md/twin/lambda/delta0 = 3/3/3/3, 5/5/5/5, 4/4/4/4, 6/7/6/6, 5/5/5/5 | `build_tournaments.py` (own builders), `tstruct.c` | `tournaments.txt`, `tstruct.out` | all confirmed. T_7 minimisers: 35 (21 cross pairs (a_x, b_y) with y - x a square, plus 14 complements). T_11: 22 (complements). D(C_9): 54. D(QR_3) md 0 |
| md(QR_p) = lambda = delta0 = (p-1)/2 for p = 19, 23 | `tstruct.c` | `tstruct_qr19_23.out` | 9/9/9/9 and 11/11/11/11 |
| QR_p doubly regular, S S^T = qI - J (p = 7, 11, 19, 23); translations and anti-automorphisms of the doubles | `build_tournaments.py` | `build_tournaments.err` | all true |
| q = 7 base case of Lemma A1 | `drt7.c` | `drt7.out` | 240 labelled DRTs, md = 3 for all |
| Lemma A1 inequality ranges | `lemmaA1_bound.py` | `lemmaA1_bound.out` | fails only at q = 7 (s = 3, 5) |
| H(T_7) = 4,730,397; all H(T_7 - v) = 686,457 | `ham_dp.c` (subset DP) and `ham_ie.c` (inclusion-exclusion over walks, mod 2^64) | `ham_dp_small.out`, `ham_ie_small.out` | confirmed by both |
| H(D(C_9)) = 6,450,103,269; all 18 deletions 773,858,997 | same | same | confirmed by both |
| H(T_11) = 17,667,825,987,845; all 22 deletions 1,822,800,617,515 | `ham_ie.c` (every deletion computed directly, no symmetry used) | `ham_ie_T11.out` | confirmed (4.3 s) |
| H(QR_7) = 189, H(QR_11) = 95,095; ratios 132.43, 1953.74; 0.3626 / 0.3928; 3.1116 / 2.9350; 1.1884 / 1.2120 | `bregman_check.py` (4) | `bregman_check.out` | confirmed |
| 16-witness H = 773,858,997; 17-witness H = 6,450,103,269 | `witness_to_masks.py`, `ham_dp.c`, `ham_ie.c` | `wit16.txt`, `wit17.txt`, `wit_recount.out` | both confirmed by both algorithms |
| THM-4995 (5): 2-regular patterns md 2, Fano md 0; Paley chain md = lambda = 4, H = 110,624,292,765; R_5 chain md 4, H = 29,678,025, max deletion 5,081,917; D(C_9{1,2,3,4}) md 3 | `build_variants.py`, `tstruct.c`, `ham_dp.c`, `ham_ie.c` | `tstruct_variants.out`, `ham_variants_*.out` | all confirmed |
| previous defect record 3: the THM-4985 bands have md 1 | `tstruct.c` | `tstruct_bands.out` | band_3(13), band_4(15), band_4(16): md = 1; H(band_4(16)) = 64,053,917 |
| Brégman table | `bregman_check.py`, `bregman_tail.py` | `bregman_check.out`, `bregman_tail.out` | first q with both ratios negative = 55; every row of `bregman_table.out` matches to 3 decimals; both ratios decrease from 55 on; negative to q = 10^6 |

The Brégman checks used exact-integer factorial logs for q <= 303 and lgamma beyond. The Robbins upper bound is
negative for 59 <= q <= 10^7.

**Part B**

| claim | program | output | result |
|---|---|---|---|
| census of connected graphs with d <= 511: 1, 1, 2, 6, 21, 112, 853, 11,117, 261,080, 199,053, 1,855, none at 12 (474,101 classes); even 1, 0, 1, 2, 8, 43, 398, 6,413, 184,596, 159,079, 1,434 | `dom_census.c census 511 16` (own canonical forms by individualisation-refinement; parity by enumerating **all** of Aut(G) from the least-form leaves; own child formula; direct recount of every d) | `census_K511.err` (class list regenerable in 9 s) | identical, level by level |
| K = 255 census (9,140 and 94 at n = 9, 10) | same | `census_K255.err` | identical |
| per-value class and even counts | `analyze_census.py` | `analyze_census_K511.out` | all 256 rows of `dom_K511_summary.txt` match |
| self-tests | `dom_census test` | `dom_census_test.out` | child formula 620,124 tests, 0 mismatches; 20,000 random relabellings, 0 failures |
| even values below 512 miss exactly {3, 7, 13}; connected even miss exactly {3, 7, 13, 25} | `analyze_census.py` | `analyze_census_K511.out` | confirmed |
| Royle-odd values below 512 miss exactly {1, 5, 37}; the 4 connected realisers of 37 are all even | `analyze_census.py` | same | confirmed (odd values = {a·k : a a connected-odd value, k odd}) |
| d = 3, 7, 13 realised only by K_2, K_3, K_4 - e (all odd); d = 25: 5 connected classes, all odd; 21 = d(K_1 v (P_3 u K_1)), even, \|Aut\| = 2 | `analyze_census.py`, `census_K37.txt` | same | confirmed |
| Euler(connected-even) = A000568 through n = 9; inverse Euler(A000568) = A334335; Euler(connected) = A000088 | `analyze_census.py` (A000568 from my own Burnside sum over odd-cycle permutations) | same | confirmed |
| Lemma B1 for every odd k <= 1023 on ceil(log2(k+1)) vertices | `partB_small.py` | `partB_small.out` | 0 failures |
| brute-force parity (all n! permutations) = census parity, 996 connected classes with n <= 7 | same | same | 0 disagreements |
| B2.1: 98,221 census classes with adjacent twins, all odd | same | same | confirmed |
| B2.2: 188 unions | same | same | 0 failures |
| B2.3: 155 cones | same | same | 0 failures |
| B2.4: 184 joins, plus 4 A v A cases | same | same | 0 failures |
| B2.5: 2,000 random (graph, permutation) pairs | same | same | 0 failures |
| B2.6: all 436 trees up to 11 vertices (`dom_census trees 11`), criterion by centres and AHU encoding | same | same, `trees11.txt` | 0 failures |
| B6: on even graphs below 32, I(G,2) takes exactly 1, 3, 9, 11, 17, 19, 25, 27, 29 and misses 5, 7, 13, 15, 21, 23, 31 | `analyze_census.py` | `analyze_census_K511.out` | confirmed |

Why the B6 check only needs n <= 7: a connected even graph has at most one universal vertex, so its complement has at
least (n-1)/2 edges. Then I(G,2) >= 1 + 2n + 4·#nonedges >= 4n - 1, and n = 8 would need I >= 33. So the n <= 7 part of
the census suffices.

## Literature and priority checks

| question | answer |
|---|---|
| Module defect of Paley / doubly regular tournaments | **Known**: Müller-Pelant 1974, the simplicity index (Finding 1) |
| Doubles D(B) of doubly regular tournaments | Not found in the literature |
| Reid-Brown 1972 | Correct (Finding 11) |
| RPGFD 2023 | Correct |
| Chan-Pak 2024 | Correct (see "Checked and found correct") |
| Brouwer 2009 | Should be Brouwer-Csorba-Schrijver 2009 (Finding 7) |
| Value set of d, or d restricted to automorphism-even graphs | No prior study found |
| arXiv 2610.28149 (Zhang, de Rham-Witt / Nygaard stacks) and 2610.31289 (Shao, integral structure sheaves on open varieties) | Exist and are described fairly; README section 3 is typed ANALOGY |

## Not audited

- Proposition B5's counts (12,326 and 18,766) were not re-run. They are evidence only.
- New witness runs started after 10:25 (`wit_d8_*`) are outside the claims under audit.

## Sources consulted

- arXiv:1907.11777 (abstract and ar5iv text)
- eudml.org/doc/12801 (Müller-Pelant metadata)
- arXiv:2003.06503
- arXiv:2308.10214v2 (HTML)
- arXiv:2204.01947
- arXiv:1701.03453
- arXiv:0905.2251
- OEIS A334335, A000568
- arXiv:2610.28149, arXiv:2610.31289
- Alon 1990 bibliographic data

## Addendum (order-16 example)

Follow-up audit of THM-4995 (4b), as worded at THM-4995 lines 115-121 and PROOF.md Theorem A5 lines 202-215, and of
the order-15 negative data. Same rules: own code (`audit/addendum_d8.py` plus the audit's `tstruct`, `ham_dp` and
`ham_ie`), writes only in `audit/` and this report.

ADDENDUM VERDICT: (4b) CONFIRMED in full. One problem: a wrong number in the order-15 negative data (Finding A1).

### Checks

**1. d8_144 = D(B8): confirmed.**
- `d8_144_block.txt` is exactly the stated code 1010100100111001101001010101.
- Rebuilding D(B8) with the audit's own `double()` (A = 0..7, B' = 8..15, b_x -> a_x, all other cross arcs A -> B')
  reproduces d8_144 arc for arc.
- B8 has scores 3, 4, 3, 4, 4, 3, 3, 4 and twin distance 2.

**2. md = lambda = delta0 = 4: confirmed.**
- Twin distance is also 4.
- md is attained by 20 sets: 12 pairs and the 8 complements of size 15.
- md >= 1 gives strong and prime.

**3. H(D8_144) = 144,227,835: confirmed** by subset DP and by inclusion-exclusion.

**4. Deletion counts: confirmed.** The 16 deletion counts take exactly the eight claimed values, with multiplicities:

| value | count |
|---|---|
| 16,750,389 | 1 |
| 17,527,823 | 3 |
| 18,351,387 | 3 |
| 18,908,079 | 1 |
| 19,218,519 | 1 |
| 19,720,143 | 3 |
| 20,574,141 | 3 |
| 21,469,761 | 1 |

Ratios: 144,227,835 = 1.81 E[H(15)], and the largest deletion is 2.02 E[H(14)].

**5. Witnesses: all 15 recounted correctly by both algorithms.**
- 1 order-15 witness, with H = 144,227,835.
- 14 order-14 witnesses: 1 in part1, 6 in part2, 7 in part3. Together they cover all eight deletion values;
  21,469,761 appears only in part3.
- So D(B8) is strong, prime, has lambda = 4 >= 3, H in S(15), and every deletion in S(14). It is a certified escaper of
  order 16 with simplicity index 4.
- The upper bound 16 in HYP-9575 (C) is therefore established.

**6. Order-15 negative data (optional item).**
- All 400 tournaments in `d87_list.txt` have the stated "unequal double" structure:
  - A is near-regular;
  - B' = QR_7;
  - there are 7 back arcs b_y -> a_f(y), with f injective and covering the four in-degree-3 vertices of A.
- md distribution: 27 / 151 / 208 / 14 for md = 1 / 2 / 3 / 4. The 14 with md 4 are exactly the session's
  `d87_md4_tags.txt`.
- Each of the 14 has exactly one deletion above 3,719,831 (THM-329's best order-13 count). So "not certifiable" is
  correct.
- The stated threshold is wrong (Finding A1).
- The extra sample claim in PROOF.md line 212 is also confirmed:
  - 300 doubles of near-regular 8-tournaments, all of the form D(B);
  - 61 with md 4, exactly the session's set;
  - H in [144,227,835, 145,280,921].

### Finding A1. MAJOR (wrong number; the conclusion is unaffected): the order-15 threshold "4.33·10^6"

Places:
- THM-4995 lines 128-129: "all have a deletion above 4.33·10^6";
- HYP-9575 line 58: "Their deletions (>= 4.33·10^6)".

What the data show (the session's own `d87_md4_tprop.out` gives the same counts):
- The largest deletion of the 14 tournaments ranges from 4,152,889 (D87_122) to 4,345,085 (D87_109).
- Only 7 of the 14 exceed 4.33·10^6.
- Every other deletion lies between 2,709,519 and 3,719,831. HYP-9575's plural "deletions (>= 4.33·10^6)" therefore
  misdescribes them twice.

Fix: "each has exactly one deletion above 3,719,831, the best known order-13 count (largest deletions
4.15-4.35·10^6), so none is certifiable". PROOF.md lines 214-215 already have this correct wording.

### Addendum files (in `audit/`)

| file | content |
|---|---|
| `addendum_d8.py` | structure checks; writes `d8_144_masks.txt` and `d87_masks.txt` |
| `addendum_d8_structure.out` | items 1 and 6 (structure) |
| `addendum_d8_tstruct.out` | item 2 |
| `addendum_d8_ham.out` | items 3-4 |
| `addendum_d8_wit_masks.txt`, `addendum_d8_wit_recount.raw`, `addendum_d8_wit_recount.out` | item 5 |
| `addendum_d87_tstruct.out`, `addendum_d87_md4_tags.txt`, `addendum_d87_md4_masks.txt`, `addendum_d87_md4_ham.out`, `addendum_d87_summary.out` | item 6 |
| `d8_sample_masks.txt`, `addendum_d8_sample_tstruct.out`, `addendum_d8_sample_md4_tags.txt`, `addendum_d8_sample_md4_ham.out` | the 300-sample check |
