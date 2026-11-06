# Independent audit of THM-4888 and HYP-9480 (kind-pasteur-2026-11-06-S13)

Audited version: the packet, `01-canon/theorems/THM-4888-*.md` and `05-knowledge/hypotheses/HYP-9480-*.md` as committed
in bf97753f19 (14:11). Earlier drafts, read at 13:48 and 13:57, were also checked. Their differences are noted where they
matter. All code in this directory is the auditor's own (see `REPRODUCE.md`). The packet's programs were not run; their
raw outputs were only read and compared.

## Verdict: PASS WITH CORRECTIONS

- **Proofs.** THM-4888 (A), (B), (D)2 and (D)3 are correct as proved. (C) and (D)1 are correct and are now, rightly,
  re-typed KNOWN. Attribution still needs precision: issue 5.
- **Computations.** Every raw number that I recomputed with independent code agrees exactly. This covers:
  - P1, all 109 seeds;
  - P2;
  - P3, both registered windows, per class;
  - P6a, both windows, all 63 classes;
  - P6b, all 6,426 (seed, class) pairs;
  - the section-4 table, every row;
  - the bank-A intensities and the weighted-capacity ratio;
  - P4: the complete `root_changes_1e8.txt`, recomputed to 10^8 by a move-only enumerator with no closed-form tail.
- **Pre-registered decisions.** I confirm all of them: P1, P2, P3, P4, P5, P6a and P6b are each CONFIRMED (PASS for P5).
- **Pre-registration integrity.** `PREREGISTRATION.md` at HEAD is identical to d36f57c903. `PREREGISTRATION_P6.md` and
  `residue_law.c` are identical to e3fdff4ceb, and `p6_seeds.txt` is identical up to CR bytes. All confirmatory outputs
  carry mtimes after their registration commits. Push times cannot be checked from a clone; only commit times can.
- **Problems.** The four MAJOR issues all concern what HYP-9480 states, or what the evidence is said to show:
  1. the kappa = 1 law, as stated, is false for seeds the packet itself names;
  2. P1 cannot distinguish the Poisson law from a no-collision model;
  3. the "Scale" arithmetic is wrong by about five orders of magnitude;
  4. the "Not explained" bank correlation is in fact fully explained by the identities (D)2-(D)3.

  None of them touches a proof or a raw number. There are also 12 MINOR issues, numbered 5-16.

## Checks

### 1. THM-4888 (A), the mod-3 law

**Proof, step by step.**
- floor(m/2^(j+1)) = sum_(i>j) b_i 2^(i-j-1), and 2 = -1 (mod 3) turns it into sum b_i (-1)^(i-j-1).
- With eps_j = (-1)^j b_j: (-1)^(i-j-1) = -(-1)^(i+j), so each cross term equals -eps_i eps_j.
- sum_(j<i) eps_i eps_j = (S^2 - sum eps_j^2)/2 is an exact integer identity.
- Since 2^(-1) = 2 (mod 3), -1/2 = 1. S = m (mod 3) gives S^2 = m^2, and sum eps_j^2 = s_2(m).
- A run of L ones contributes L to s_2 and L - 1 to a, so s_2 - a = r.

Every step is valid. The digit identity it uses is Session 206's (PARAMETER_CAPACITY.md), quoted correctly. (One
notation nit, issue 16: the canon writes "S = m" for "S = m (mod 3)".)

**Computation.** `a_identities.c`, output `out_identities_2p28.txt`. Q was built by its recurrence for every m < 2^28,
ahead of the packet's 2^26. The run found 0 failures of:
- the digit identity;
- Q = m^2 - r(m) (mod 3);
- Q even, and Q <= m - 1.

A deliberately wrong law (+r in place of -r) fails on 2,366 of the m below 4096, so the test is not vacuous. A
targeted search of the Session 205 and 206 packets found no prior statement of (A).

### 2. THM-4888 (B), involution rigidity

**T and the identity.**
- Session 206's T(m) = 2^n - 1 - rev_(n-1)(m - 2^(n-1)) gives T(9) = 11, T(4) = 7, T(11) = 9 and T(7) = 4.
- For 1 <= m < 2^28: T is an involution of each dyadic block, and Q(T(m)) - Q(m) = 2(T(m) - m) holds with no exception.
- The number of fixed points is 16,383 = 2^14 - 1, as Session 206's count predicts.
- The identity H_s(T(m)) - H_s(m) = (R - 2d)(T(m) - m) - e((T(m) mod 2) - (m mod 2)) is a one-line consequence.

**The "iff".** At m = 9, T(m) - m = 2 with equal parity, forcing R = 2d. At m = 4, T(m) - m = 3 and the parity term is
1, so 3(R - 2d) - e = 0 forces e = 0. The converse is immediate.

**Classification.** The predecessor formulas in the proof match Session 198's exact inverse (bilinear-descent PROOF.md
section 3) letter for letter, including the side conditions and the ranges of q:
- O needs e = 0 and f - a even, with 0 <= q <= min(a, p).
- Z needs p = 0 and a - f - e even and nonnegative, with 0 <= q <= e.
- J needs f = 0, with 0 <= q <= min(a, e - p).

At y = (0, 0, 0, d):
- O gives (0, -d, 0, d), which is negative.
- Z is impossible because p = d >= 1.
- J has an empty range, min(0, -d) < 0.

So "y has no physical predecessor" is right. With Session 198 Theorem 1 (every nontrivial response other than t has a
physical predecessor), this forces s = t.

**Independent BFS** (`b_states.py`, output `out_states_U4096.json`).
- I re-derived the inverse by brute-force solving the moves. On all 4,353 states with U <= 512 it equals the printed
  formulas, and every BFS parent is among the printed predecessors.
- BFS with U <= 4096 finds 94,010 states (194 with U <= 64). Only (1, 0, 2, 1) has U = 2d and E = 0; 5,678 states
  have U = 2d with E > 0.
- y = (0, 0, 0, d) has no nonnegative predecessor for d = 1..50.

**Completeness of a U-bounded BFS.** Session 198's inequalities (1) give U >= D + d - E >= 2d. The U-increments of O, Z
and J are therefore U, U - d >= d and D, all at least 1. U strictly increases along every move, so every ancestor of a
state with U <= UMAX lies inside the search. All enumerated states satisfy (1), U >= 2d and D odd, and no move fails to
raise U.

### 3. THM-4888 (C), confinement and parity-split capacity (KNOWN)

**Proof.** The congruence is immediate from the tail formula. In Session 206's capacity proof, u_L(m) = b_0/2 + b_1/4 +
..., so j >= K/2 holds exactly when b_0 = 1: the parity is the leading digit of the truncated bit-reversal. Restricting
to a parity is legitimate:
- the sandwich between the H-count and the F-counts consists of pointwise inclusions in m, which survive any
  restriction of m;
- each residue-class count still errs by O(1);
- a Riemann sum over a half-interval still errs by O(1/K).

The parity-split integrals are correct.

**Numerics** (`h_parity.c`, output `out_parity_capacity.txt`). Ten seeds, among them the root, a d = 1 seed (3, 3, 1, 1),
JJ, (81, 99, 27, 45) and (159, 160, 78, 79), at X = 2^24 and 2^26. Both parity counts agree with
(X/d) ln(R/(R - d/2)) and (X/d) ln((R - d/2)/(R - d)) to within 0.006 sqrt(X).

**Prior art.**
- Both halves are in Session 206's weighted-parameter-capacity theorem (de303ea26d, pushed 13:38): section 1 is the
  parity split, section 2 the zero outside the coset. That packet's audit (`audit/REVIEW.md`) checks the 1/K
  normalisation independently.
- The congruence is also, and earlier, the necessity half of Session 206's exact modular image
  (`exact-modular-image/DYADIC_SURJECTIVITY.md`, "Exact modular image for every seed", with M = 2p and g = p). That
  packet was in this checkout before the session started and is cited in PROOF.md section 1. See issue 5.

### 4. THM-4888 (D), and the literal tournaments

**Algebra.**
- J then Z gives (D + 2U, D, D + 2U, D), and O then J gives the same.
- (D)2 reads off the tail formula; Z depends on (D - E, U, d) and O on (D, U, d).
- (D)3: (h - d) + d(1 + 2m - Q) - e(m mod 2) = h + 2dm - dQ - e(m mod 2).

**Numerics** (`out_D_identities.txt`). (D)3 holds on all 16 R = 2d bank-A seeds for m < 2^20. The 11 bank-A pairs that
share (D - E, U, d) have identical odd-m halves and identical Z-children. No pair shares (D, U, d).

**Literal construction** (Session 198 section 2; Session 206 ROOT_OBSTRUCTION section 1), as built in `c_literal.c`:
- t is root c plus one vertex, with c -> v0.
- O appends a vertex beaten by everything.
- Z appends a vertex w with w -> root, and every other vertex beats w.
- J does the same as Z and then makes w the root.

D, E, U, d and the number of root-start Hamiltonian paths are counted by subset DP over actual adjacency matrices.
U = sum over S containing the root of G(S) H(V - S), with H(empty) = 1. No formulas are used. Results:
- **Moves.** All 29,524 words of length <= 9 (up to 11 vertices): the literal (D, E, U, d) equals the move arithmetic,
  and exactly one Hamiltonian path starts at the root. In all 3,280 JZ / OJ word pairs the two literal responses are
  equal.
- **Tails.** 7,744 tail tournaments (prefix of length <= 4, tail m < 64, up to 12 vertices): the literal H equals
  h + Rm - dQ(m) - e(m mod 2), with Q from Session 205's pair-sum definition.
- **Seeds.** The 85 bank-A seeds and 10 P1 seeds, built literally from their BFS words, have exactly the listed
  (h, R, e, d). All their tails with m < 16 match the formula (1,520 tournaments).

### 5. HYP-9480 (A) evidence

**Engine** (`d_image.c`). A different algorithm from the packet's pruned DFS: Q is a full recurrence array, and every m
from 0 to the complete horizon (H >= h + (R - d) m) is scanned linearly into a window bitmap.

**Raw counts.**
- P1: all 109 hit counts are identical to `p1_density.txt`.
- P2: 33,456,317 hits, density 0.498538, identical.
- P6a, root at [2^25, 2^26) and [2^26, 2^27): all 63 (n_c, p_c, v_c) identical to `p6_root_w2*.txt`.
- P6b: all 102 seeds x 63 classes identical to `p6_seeds_w24.txt`.

**Predictions** (`e_decide_p1_p6.py`). I implemented the registered local-factor model from the text of
PREREGISTRATION.md. It equals `verify.predict` on all 109 seeds; 21 of them carry a local factor. The registered text
leaves the piece intensities implicit. The parity-split capacities are the only natural reading, and they are what
verify.py uses (issue 16).

**Decisions** (`out_decisions_p1_p2_p6.json`).

| prediction | registered criterion | result | outcome |
|---|---|---|---|
| P1 | >= 95% within 0.004, all within 0.010 | 109/109 within 0.004; 108 within 2 x 10^-4; maximum 0.0019986 at (81, 99, 27, 45) | CONFIRMED |
| P2 | density in [0.4978, 0.5005] | 0.498538 | CONFIRMED |
| P6a | every class within 0.006; mod-7 span >= 0.05 | maximum deviation 0.002175 (mod 63, 2^25); mod-7 spans 0.0902 and 0.0803 | CONFIRMED |
| P6b | >= 95% within 0.004, all within 0.010 | 5,908 live pairs, all within 0.000841; 95th percentile 0.000118; mean 4.23 x 10^-5; class-blind mean 0.00247; the 518 zero-intensity classes are empty | CONFIRMED |

**Alternative models** (issue 2). The no-collision model miss_c = 1 - kappa lambda_c:
- passes P1 (109/109 within 0.004; maximum 0.0031);
- fails P6b only through its "all within 0.010" clause: 5,892 of 5,908 pairs are within 0.004, and the maximum is 0.032;
- is off by more than 0.2 for the root (P6a).

The class-blind model fails P6b (82% of pairs within 0.004).

**THM-4878 (5).** The root image below 2^21 contains exactly 521,915 odd values. Its mod-7 miss rates are 0.455, 0.487,
0.551, 0.581, 0.533, 0.466 and 0.443; mod 3 they are 0.501, 0.490, 0.516; mod 5, 0.501-0.503. These reproduce THM-4878
(5) (0.487 against 0.488 is rounding), and confirm that the one-apex image equals the root image there.

**Mod-7 spans over seven scales.** Miss-rate spans 0.133, 0.105, 0.091, 0.090, 0.080 and 0.074. Intensity spans 0.526,
0.410, 0.360, 0.356, 0.319 and 0.293, against a uniform limit of 2 ln 2 = 1.386. This is consistent with Session 206's
proof that the root's per-class intensity becomes uniform. The coordinator's new paragraph "The bias decays" is
correct.

**Exceptions** (`out_exceptions_w22_summary.txt`; [2^22, 2^23), mod 63). Mean signed deviation of miss_c from
exp(-lambda_c):

| kind | seeds | mean signed deviation |
|---|---|---|
| d = 1 | (1, 4, 0, 1), (3, 3, 1, 1), (5, 7, 1, 1), (1, 8, 0, 1) | +0.024 to +0.124 |
| R = 2d | (5, 6, 2, 3), (9, 10, 4, 5), (17, 18, 8, 9), (45, 50, 20, 25) | +0.021 to +0.105 |
| balanced, nested in an R = 2d seed | (9, 9, 3, 3), (15, 15, 5, 5), (25, 25, 5, 5) | +0.024 to +0.114 |
| balanced, nested in a generic seed | (53, 53, 17, 17), in (27, 35, 9, 17) | +0.0001 |
| generic, in-sample | 5 seeds with lambda 0.06-0.24 | within 3 x 10^-4; no-collision model 0.002-0.028 |

(Issue 1.)

### 6. Section 4 table (`f_section4.py`, `out_section4_w22.json`; [2^22, 2^23))

**Table rows.** Every row reproduces.

| bank | union miss | product of the members' miss rates |
|---|---|---|
| bank A, 85 seeds (mod 3) | 0.02187 | 0.00355 |
| without the 35 balanced seeds | 0.02187 | 0.01573 (mod 3), 0.01692 (mod 945), 0.01574 (mod 3 x 4096) |
| the 34 generic seeds (mod 945) | 0.09086 | 0.07316 |
| the 16 R = 2d seeds (mod 945) | 0.20915 | 0.19680 |

**Nesting.** Bank A has exactly 35 balanced seeds, 16 R = 2d seeds (root included) and 34 generic seeds, and no
non-balanced seed with h = R. Each balanced seed equals J(O(y)) for a physical y. Following J(y), which may itself be
balanced, every chain ends at a non-balanced bank-A seed: 12 end at R = 2d seeds and 23 at generic seeds. Removing the
35 therefore cannot change the union (issue 13).

**Shared halves.** (25, 30, 10, 15) and (27, 30, 12, 15) meet in exactly 60,807 values: their two odd-m halves, which are
identical sets. Their even halves are disjoint, in classes 10 and 12 mod 15. The packet's baseline of 9,714 is resolved
mod 3 only; at the seeds' confinement modulus 15 the baseline is 26,447 (issue 12).

Among all 1,225 pairs of the 50 non-balanced seeds, the 11 largest excesses are exactly the 11 (D)2 pairs (excess at
least 24,996). The next is 532.

**Decomposition** (`out_w22_halves_product.txt`; issue 4).
- The product over the 62 distinct half-images of the 34 generic seeds is 0.0908 (mod 945), against their union 0.0909.
- The 16 R = 2d seeds give a union of 0.2092, against 0.1530 for the product of their distinct halves.
- The product of the exact generic-block union and the exact R = 2d-block union is 0.02189, against the bank's 0.02187.

**New section-4 numbers** (`out_bankA_intensity_and_weighted_test.txt`).
- 2 sum c_s = 8.288. It is 7.595 with the root's involution, 5.153 without the 35 nested seeds, and -ln(0.02187) =
  3.823. All are correct.
- The weighted-capacity ratio has its minimum 4.699 at (p, a) = (75, 1), so "at least 4.7" is right. It counts the 35
  nested seeds; without them the minimum is 3.39. The conclusion is unchanged.

### 7. P3 (`out_p3_bankA_w2*_M3.txt`)

- The per-class miss counts are identical to the packet: 734,628 misses (0.021894) at [2^26, 2^27) and 1,457,666
  (0.021721) at [2^27, 2^28).
- At [2^24, 2^25): 0.022086, with class rates 0.0060, 0.0411 and 0.0191 mod 3. This fits the slow downward drift.

### 8. Root-change families F_k (`g_rootchanges.c`)

**Enumerator.** It walks every word t . w_1 (J w_2) ... with the moves only, using no closed-form tail (the moves are
the ones verified literally in check 4). A word is recorded at its last Z or J letter, or at t itself; trailing O letters
are not recorded. The pruning is complete:
- D never decreases along a word;
- after an O the next Z or J gives D + 2U, so an O-child is explored only when D + 2U <= X.

This is the same convention as `root_changes.c`: the final block is canonical without trailing O, and intermediate
blocks keep all their letters.

**Results.**
- **X = 4,000,000** (3 s). F_0, F_(<=1) and F_(<=2) miss 1,005,150, 26,330 and 140 odd N. The largest F_(<=2) miss is
  3,780,967. The 15 misses above 10^5 are exactly the list in the brief, 14 of them = 1 (mod 3); the exception is
  136,607.
- **X = 10^8** (264 s, `out_rootchanges_X1e8.txt`). Misses 25,095,920, 224,380 and 140. All 27 dyadic-window lines are
  identical to `root_changes_1e8.txt`. So P4 is independently confirmed:
  - (a) no miss in (4 x 10^6, 10^8];
  - (b) a single miss above 10^5 not = 1 (mod 3), namely 136,607;
  - (c) nine ratios from 1.5300 to 1.6311.
- **Cross-check.** F_0 hits 24,904,080 odd values up to 10^8, matching Session 205.

**Masses.**
- N_0(X)/X - ln 2 is at most 5 x 10^-5.
- For exactly one root change, (N_1 - (pi^2/12) X ln X)/X = -1.1693, -1.1710, -1.1708, -1.1709 and -1.1709 at
  X = 10^5, 10^6, 4 x 10^6, 10^7 and 10^8. This confirms -1.17. For "at most one" the constant would be -0.478 (issue 9).
- The exact three-point fit for exactly two root changes is 0.5002 X ln^2 X - 2.0296 X ln X + 4.1317 X. It extrapolates
  to 10^8 within 0.12%. The rounded (0.50, -2.0, 4.2) is 0.5-0.8% off (issue 10).

**Class 1 (mod 3)** (`i_classmass.c`). In [X/2, X] with X = 4 x 10^6, the representation mass of exactly-2 words splits
0.494 / 0.245 / 0.261 over the classes mod 3. Final blocks with 3 | gcd(d, R) carry 24% of the mass, and about 80% of
the class-1 deficit. That supports the "mechanism" remark for F_2 as a heuristic (issue 6).

### 9. Citations, status labels, wording

**Verified against the sources.**
- Session 205: the tail formula, the Q recurrence, the digit identity, JZ = OJ, and 24,904,080.
- Session 206:
  - the capacity c_s and its residue-class proof;
  - T and Q(T(m)) - Q(m) = 2(T(m) - m);
  - the d = 1 collapse;
  - root residue capacity;
  - the exact modular image, with (5, 6, 2, 3) omitting 1 (mod 3);
  - the weighted capacity: coefficients, coset, example and the X/252 bound.
- Session 198: the moves, the positive coordinates and the inverse.
- Session 199 BALANCED_SOURCES.md: "Z J = J O = M_2" with right-to-left composition, which is JZ = OJ in word order.
- THM-4878 (5): opus-2026-11-04-S79, 521,915 / 1,048,576, the mod-7 range 0.443-0.581, flat mod 16 and mod 5.
- THM-4885: opus-S82, (pi^2/24) X ln X for one-apex x one-apex pairs u <= v.
- S82's source note: F1 miss density N^-0.35, ratios 1.5879 / 1.5878 / 1.5876, class means (60.443, 41.411, 42.738).

**Status labels.** THM-4888 and HYP-9480 now separate PROVED, KNOWN, CONJECTURE, EXPLORATORY and FINITE-EXACT
correctly in their status lines. Five places in the body still blend them: issues 2, 3, 6, 7 and 16.

## Issues

### MAJOR

**1. The kappa = 1 law is false as stated: its domain must be restricted.**

*Where.*
- HYP-9480 title: "kappa = 1/2 for the root (converse involution) and 1 otherwise".
- HYP-9480 (A) and PROOF.md section 3, "The law": "and kappa = 1 for every other seed".
- HYP-9480 (A): "Exceptions are exactly the extra symmetries of THM-4888 (B) and (D): ... balanced seeds, which are
  nested". PROOF.md section 3 has the same list under "**The exceptions are exactly the extra symmetries.**".

*Problem.* The packet itself lists d = 1 seeds, R = 2d seeds and balanced seeds as exceptions. Measured at
[2^22, 2^23) mod 63 (check 5), they violate kappa = 1 by +0.02 to +0.14 per class.
- Nesting is not the cause for balanced seeds. (53, 53, 17, 17), nested in the generic (27, 35, 9, 17), follows the law
  to 10^-4. (9, 9, 3, 3), (15, 15, 5, 5) and (25, 25, 5, 5), nested in the R = 2d seeds (5, 6, 2, 3) and (9, 10, 4, 5),
  inherit their collisions.
- That the list is "exactly" the exceptions is untested.

*Correction.*
- In the HYP-9480 title, replace "and 1 otherwise" with:
  > and kappa = 1 for generic seeds (d > 1, R != 2d, and not nested by JZ = OJ in a seed with d = 1 or R = 2d)
- In HYP-9480 (A) and PROOF.md section 3, replace "and kappa = 1 for every other seed." with:
  > and kappa = 1 for generic seeds: d > 1, R != 2d, and not of the form W O J with W J a d = 1 or R = 2d seed.
  > Seeds with d = 1, seeds with R = 2d, and balanced seeds nested in them violate kappa = 1 by 0.02 to 0.14 per class
  > (audit).
- Replace "Exceptions are exactly the extra symmetries of THM-4888 (B) and (D):" (and the bold PROOF.md heading) with:
  > Known exceptions, each from THM-4888 (B) or (D) (the list is not shown to be complete):
- Replace the bullet "balanced seeds, which are nested" with:
  > balanced seeds W O J whose parent W J is itself exceptional, which inherit its collisions. A balanced seed nested
  > in a generic seed, such as (53, 53, 17, 17) in (27, 35, 9, 17), follows kappa = 1.

**2. P1 does not test the Poisson form, and P6b tests it only weakly.**

*Where.*
- HYP-9480 status: "(A) CONFIRMED out of sample (P1: 109/109 unseen seeds within 0.004; ...)".
- HYP-9480 (A), Evidence.
- PROOF.md section 0 item 2 ("The law predicts the miss rate ... within 0.00084") and section 3, Evidence.
- README, Status.

*Problem.* P1 seeds have lambda <= 0.04. There the Poisson law and a no-collision model (each parameter its own value,
miss = 1 - lambda) differ by at most lambda^2/2 = 0.0008, far inside the 0.004 tolerance.
- The no-collision model passes P1 (109/109; maximum 0.0031). P1 therefore confirms the intensities, which are Session
  206 theorems, and the confinement; it does not confirm the collision law.
- P6b separates the two models only through its "all within 0.010" clause. The no-collision model has 5,892 of 5,908
  pairs within 0.004 and fails on 16 high-intensity classes, mostly (81, 99, 27, 45).
- The Poisson form's tenfold smaller mean error (4.2 x 10^-5 against 4.3 x 10^-4) is real, but it was reported, not
  registered.
- The decisive out-of-sample tests are the root's, P2 and P6a, where the no-collision model is off by more than 0.2.

*Correction.*
- In the HYP status, replace the (A) clause with:
  > (A) CONFIRMED out of sample. For the root (kappa = 1/2), P2 and P6a are decisive: a no-collision model misses by
  > more than 0.2. For generic seeds (kappa = 1), P1 and P6b are CONFIRMED but discriminate weakly. The no-collision
  > model miss_c = 1 - lambda_c also passes P1 (109/109 within 0.004) and fails P6b only on 16 of 5,908
  > high-intensity classes. The Poisson form's mean error, 4.2 x 10^-5 against 4.3 x 10^-4, is reported, not
  > pre-registered.
- In PROOF.md section 0 item 2, after "within 0.00084.", add:
  > At these intensities, a model with no collisions at all is within 0.004 for 99.7% of the classes. The decisive
  > out-of-sample tests of the Poisson form are the root's (P2, P6a).
- In PROOF.md section 3 Evidence, replace "P1: 109 unseen seeds, overall density." with:
  > P1: 109 unseen seeds, overall density. This tests the intensities and the confinement. It does not separate the
  > Poisson law from a no-collision model, which also passes.

**3. The "Scale" estimate is wrong by about five orders of magnitude, and S82 shows F2 is not Poisson.**

*Where.*
- PROOF.md section 4, "Scale": "This is what (C) predicts: even with bank A's loss of half the intensity, its first
  misses would be expected only beyond about 10^15."
- HYP-9480 (C): "That is consistent with (C): its first misses would be expected only far beyond 10^15."

*Problem.* Take S82's hardest class, class 1 mod 3, with mean 41.411.
- With half the intensity, lambda = 20.7. The class has 2^35/6 odd N up to 2^35, so a Poisson model expects
  2^35/6 x e^-20.7 = 5.8 misses there; the first would be expected near 6 x 10^9 < 2^35.
- The census found none, so under the stated assumption it would contradict, not confirm, the transfer.
- A first miss beyond 10^15 needs a loss of at most about 20%. Bank A's after-nesting ratio 3.82/5.15 gives 1.3 x 10^14;
  with no loss the first miss is near 6 x 10^18.
- S82 also records that F2 is far from Poisson: P(r2 <= 10) = 4.6 x 10^-6, about 10^5 times the Poisson value. A Poisson
  extrapolation of the first-miss location is unsupported either way.

*Correction.*
- In PROOF.md, replace the two sentences from "This is what (C) predicts" through "either way." with:
  > Under a Poisson model with no loss, the first misses of F2 in its hardest class (mean 41.4) would be expected only
  > near 6 x 10^18. With bank A's loss of half the intensity they would be expected near 6 x 10^9, inside the census,
  > which found none. S82 also found F2's lower tail about 10^5 times heavier than Poisson
  > (P(r2 <= 10) = 4.6 x 10^-6). So the census neither confirms nor refutes (C). It shows only that a positive miss
  > proportion can lie beyond any feasible census.
- In HYP-9480 (C), replace "That is consistent with (C): its first misses would be expected only far beyond 10^15." with:
  > A Poisson model would place its first misses anywhere from about 6 x 10^9 (with bank A's halving of the
  > intensity) to about 6 x 10^18 (no loss). S82's F2 is far from Poisson, so this census does not test (C).

**4. The "Not explained" correlation is fully explained by (D)2 and (D)3.**

*Where.*
- PROOF.md section 4, bullets "*Dilated roots.*" and "*Not explained.*".
- HYP-9480 (B): "The 34 generic seeds still miss 24% more than their product."

*Problem* (check 6).
- The 34 generic seeds contain six (D)2 shared-half pairs. Merging identical halves, the product over the 62 distinct
  half-images is 0.0908 (mod 945), against the union miss 0.0909. The 24% excess is entirely the shared halves.
- For the 50 non-balanced seeds, the generic block and the 16-seed R = 2d block are independent: the product of the two
  exact block unions is 0.02189, against 0.02187.
- All remaining correlation sits inside the dilated-root block: union 0.2092, against 0.1530 for its independent
  halves.
- The 11 largest pairwise excesses are exactly the 11 (D)2 pairs; the next is about 50 times smaller.
- No correlation beyond the identities is visible.

*Correction.*
- In PROOF.md, replace the "*Not explained.*" bullet with:
  > *Nothing left over.* The 34 generic seeds behave as independent half-images. The product of the miss rates of their
  > 62 distinct half-images, per class mod 945, is 0.0908, against their union miss 0.0909. Their 24% excess over the
  > seed product is therefore the six shared-half pairs of (D)2. The 16 R = 2d seeds are strongly correlated among
  > themselves (union 0.2092; product of distinct halves 0.1530). The generic and dilated-root blocks are independent
  > of each other: the product of the two block unions is 0.02189, against the bank's 0.02187. Residue or position
  > resolution cannot close the seed-level gap, because the gap comes from the identities (D)1-(D)3.
- Replace the "*Dilated roots.*" bullet with:
  > *Dilated roots.* The R = 2d seeds of (D)3 are translates and dilates of the two root halves. As a block they miss
  > 0.2092, against 0.1530 for independent halves. This is the only correlation left after (D)1-(D)2.
- In HYP-9480 (B), replace "The 34 generic seeds still miss 24% more than their product." with:
  > The rest is accounted for. Once their six shared-half pairs (D)2 are merged, the 34 generic seeds behave as
  > independent half-images (0.0908 against 0.0909). The remaining gap lies inside the dilated-root block (D)3 (audit).

### MINOR

**5. Attribution: (C), (D)1 and (D)3.**
- (C): the congruence is also the necessity half of Session 206's exact modular image. That packet predates this
  session and is cited in PROOF.md section 1, so "found independently" should be confined to the parity split.
- (D)1: JZ = OJ is stated in Session 205's ARITHMETIC_REDUCTION.md itself (line 44), not only in its audit.
- (D)3: generalises Session 206's identity H_s = 3H + 2 - 2(m mod 2) for (5, 6, 2, 3) (ROOT_RESIDUE_CAPACITY section 4).

Suggested THM-4888 status text:
> (C) is KNOWN: the congruence is the necessity half of Session 206's exact modular image (DYADIC_SURJECTIVITY.md); the
> parity split and the exact weights (g/p) c_b are Session 206's weighted-parameter-capacity (de303ea26d), the parity
> split found independently here. (D)1 JZ = OJ is KNOWN (Session 205 ARITHMETIC_REDUCTION.md and its audit; Session 199,
> ZJ = JO in operator order). (D)3 generalises Session 206's (5, 6, 2, 3) identity (ROOT_RESIDUE_CAPACITY section 4).

Apply the same changes to PROOF.md section 7.

**6. A heuristic inside a theorem.**

*Where.*
- THM-4888 (C) corollary: "which is the mechanism behind the hard class 1 (mod 3) of root-change and join families".
- PROOF.md section 2, Consequence 3.
- PROOF.md section 4 and HYP-9480 (C): "The class-1 deficit is THM-4888 (C)".

*Problem.* This is an interpretation, not part of the proof. For F_2 it is supported: confined final blocks carry about
80% of the class-1 mass deficit. For S82's families, S82 derived its own analogous local factor (3 | y_h and 3 | y_R).

*Correction.* Replace with:
> Heuristic remark (not part of the theorem): this confinement carries most of the class-1 (mod 3) deficit of F_2
> (about 80% of the deficit of representation mass in [X/2, X], X = 4 x 10^6, audit). S82 found the analogous
> 3 | y_h, 3 | y_R local factor for its F2.

**7. An exploratory number inside the theorem file.**

*Where.* THM-4888, final paragraph: "In bank A, the 35 balanced seeds are redundant, and the union misses six times the
product of the members' miss rates."

*Problem.* Most of the "six times" is the trivial nesting redundancy.

*Correction.* Replace with:
> Exploration (HYP-9480 (B)): after the 35 nested balanced seeds of bank A are removed, the union misses 1.39 times the
> product of the members' miss rates. That excess is accounted for by (D)2 and (D)3.

**8. The S82 constant is off by a factor of 2.**

*Where.*
- PROOF.md section 5 table: "the pi^2/12 matches S82's THM-4885 (B)".
- PROOF.md section 5, Misses: "with the same mass constant pi^2/12".
- HYP-9480 (D), in both places.
- PROOF.md section 0: "The constant for one root change is the one opus-S82 proved for its one-apex joins."

*Problem.* THM-4885 (B) is (pi^2/24) X ln X for pairs u <= v, with mean multiplicity (pi^2/12) ln N. F_1 has
(pi^2/12) X ln X, with mean (pi^2/6) ln N: twice S82's F1. And F_1's constant here is numerical, not proved.

*Correction.* Replace with:
> the constant pi^2/12 = integral_0^1 ln(1 + t)/t dt of THM-4885 (B), whose unordered pairs u <= v have mass
> (pi^2/24) X ln X. F_1 here has twice that mass (mean (pi^2/6) ln N against (pi^2/12) ln N); its constant is numerical,
> not proved.

**9. The mass convention is not stated.**

*Where.* PROOF.md section 5, table header; HYP-9480 (D), "Masses (words with H <= X)".

*Problem.* F_k is defined with at most k root changes, but the masses count exactly k. The constant -1.17 holds only for
"exactly 1"; "at most 1" gives -0.478.

*Correction.* Use the header:
> N_k(X) = #(words with exactly k root changes and H <= X; trailing O letters of the final block deleted, intermediate
> blocks with all letters)

**10. The k = 2 fit and multiplicity.**

*Where.* PROOF.md section 5.

*Problem.* The exact three-point fit for exactly two root changes is 0.5002 X ln^2 X - 2.0296 X ln X + 4.1317 X; the
rounded coefficients are 0.5-0.8% off. "Roughly 230 at N = 4 x 10^6" is the leading term ln^2 N alone; the fitted mean
multiplicity is about 204.

*Correction.* Change the formula to:
> about 0.500 X ln^2 X - 2.03 X ln X + 4.13 X

and change "roughly 230" to:
> about 204 (231 from the leading term ln^2 N)

**11. "Unseen" P6b seeds.**

*Where.* PROOF.md sections 0 and 3, and HYP-9480 status and (A): "102 unseen (generic) seeds".

*Problem.* P6 was registered at 13:43:24 on the P1 seeds and the P1 window. Their overall densities there had been known
since 13:16.

*Correction.* Replace with:
> 102 generic P1 seeds (not used in the P6 exploration; their overall densities in this window were known from P1)

**12. The shared-halves example.**

*Where.* PROOF.md section 4: "They meet in 60,807 values, against 9,714 expected for independent sets".

*Problem.* The baseline is resolved mod 3 only, although both seeds are confined mod 15. And the overlap is an exact set
identity, not a statistical excess.

*Correction.* Replace with:
> They meet in exactly their common odd-m half, 60,807 values. Independent sets with the same mod-15 confinement would
> meet in about 26,400.

**13. The nesting chain.**

*Where.* PROOF.md section 4: "J(y) is itself in bank A, or has d = 1 and lies inside the root image."

*Problem.* For bank-A seeds d > 1, so the second branch never occurs. J(y) can itself be balanced, for example
(15, 15, 3, 3) -> (9, 9, 3, 3) -> (5, 6, 2, 3).

*Correction.* Replace with:
> J(y) has the same d > 1 and a smaller U, so it is in bank A; if it is balanced, repeat. All 35 chains end at
> non-balanced bank-A seeds (12 at R = 2d seeds, 23 at generic seeds; audit).

**14. Forced collisions.**

*Where.* PROOF.md section 4, "What a proof needs": "by (B) the involution that halves the root's image exists for no
other seed". HYP-9480 (C): "Collision structure cannot suffice either."

*Problem.* (B) excludes only full T-invariance.

*Correction.* Replace with:
> (B) excludes the full root involution for every other seed. But R = 2d seeds keep a partial involution (T preserves
> H_s when m and T(m) have the same parity), d = 1 seeds lie inside the root image, and other collision structures are
> not excluded.

**15. The mod-7 wording, and its credit.**

*Where.*
- README: "which also explains the mod-7 bias of the one-apex image".
- PROOF.md section 3: "The parameter intensity of the root is not uniform mod 7."
- PROOF.md section 3, "The bias decays": credits only the weighted capacity.

*Problem.* The README still uses the superseded wording. And Session 206's ROOT_RESIDUE_CAPACITY (16e5d8ee81, 12:57,
before the P1-P5 registration) already proved the uniform limit for the root.

*Correction.*
- In the README, use:
  > which accounts, at each finite scale, for the slowly decaying mod-7 bias of the one-apex image
- In PROOF.md section 3, use:
  > is not uniform mod 7 at the tested scales
- In "The bias decays", credit:
  > Session 206 root-residue-capacity (root) and weighted-parameter-capacity (all seeds)

**16. Provenance, wording and small numeric nits.**

- **Provenance.**
  - `root_changes.c` says it is "byte-identical in logic to the exploratory kj.c that produced root_changes_1e8.txt",
    while the README says `./root_changes` produced it. Make the two consistent.
  - `density.c` was last modified at 13:28:22, after `p1_density.txt` (13:16:37), `p2_root.txt` (13:17:16) and
    `p3_bankA_w26.txt` (13:24:01). The outputs are reproduced exactly by the audit's independent engine, so no number
    is in doubt. Still, say which binary produced which file, or rerun (each run takes at most 4 minutes).
- **Pre-registration gap.** The registered P1 model leaves the piece intensities implicit; (C) was not among the
  registered claims. Add to section 6:
  > The registered model leaves the piece intensities implicit; verify.py uses the parity-split capacities.
- **Session 199.**
  - Current: "Session 199 found a similar census for the balanced readout with at most three root changes".
  - Problem: Session 199 counts J letters of the source word s, and the readout H(s O J) adds one more J.
  - Replacement:
    > Session 199 found the balanced readout H(s O J) complete from 27,783 to 10^6 for source words s with at most three
    > J letters. Two J letters still miss 293 values, the largest 879,001.
- **Overconfident phrasing.**
  - "**Answer:** almost certainly not" (PROOF.md and README): replace with "Conjecturally no (HYP-9480 (C))".
  - "One free block can never suffice because its multiplicity is bounded (that is (C))" and "The trichotomy says why
    one free block is never enough": replace with "Conjecturally (HYP-9480 (C)), one free block never suffices: its
    mean multiplicity is bounded (Session 206)".
  - "(D)(2) plus a finite check would prove Chan-Pak": replace with "an effective form of (D)(2), with an explicit N_0,
    plus a finite check below N_0 would prove Chan-Pak".
- **Small nits.**
  - Post-hoc note: 0.02393 should be 0.02392 (the value is 0.0239208).
  - "Exploration at 2^20, 2^22 and 2^24 found spans ... with the same agreement": at 2^20 the mod-63 maximum deviation is
    0.0091, above P6a's 0.006, though mod 7 it is 0.0010. Say "the same agreement mod 7".
  - "the root sits about 0.001 above the law in every class": true mod 7 (+0.0005 to +0.0010) and mod 9 (+0.0012 to
    +0.0021), but 2 of 63 (2^25) and 10 of 63 (2^26) classes mod 63 sit slightly below. Say "in every class mod 7 and
    mod 9".
  - THM-4888 title, "the tail image of a physical seed is T-symmetric": say "the tail map H_s is T-invariant". In
    PROOF.md section 2 (B), "every other seed's image lacks it" should read "no other seed's tail map is T-invariant".
  - THM-4888 (A): "S = sum (-1)^j b_j = m" should read "= m (mod 3)".
  - The weighted-capacity factor "at least 4.7" counts the 35 nested seeds; without them it is 3.39. The conclusion is
    unchanged.

## What was not checked

- The "first-generation children" and "16 position bands" products of section 4.
- The 0.022 in every dyadic window from 2^10.
- The P6 exploration's in-sample figures (33 seeds, 1,890 pairs).
- Whether a push happened before each run. Only commit times are available in a clone.
