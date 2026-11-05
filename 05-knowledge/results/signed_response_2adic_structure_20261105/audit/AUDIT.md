# Independent adversarial audit of THM-4883 / HYP-9476

Auditor: independent agent, 2026-11-05, own code only. Nothing was imported from the packet: no `verify.py`, no
`dependencies/bilinear.py` or `core.py`, no `witness_length.c`. g is evaluated three independent ways; the one-apex
sources are built from the O/Z/J birth rules; the C2 parent/child graphs of T4 are built directly from the text of
BILINEAR_C2 (Section 2 graft, Section 3 pair forms N and M, Section 4 converse), and every Hamilton path, cycle,
endpoint and cover count is a subset DP in `audit_census.c`.

Files: `audit_census.c` (census to 2^24/2^26, subset DP, witness-length DFS), `independent.py` (driver),
`INDEPENDENT.json` (deterministic output; `python independent.py --exe <binary> --out INDEPENDENT.json`, about 5
minutes; the binary is compiled outside the repository). The `-O` rerun is byte-identical.

## Verdict: PASS WITH CORRECTIONS

The new mathematics is correct. T1 (one full class modulo 2^(s+1) iff 1/2 < a/b < 1 and b' > 1, every other class
one-side bounded), T2 (the digit formula for C*), T3 (exact half-lines and thresholds), Lemmas 1-6, Theorem A0 and the
T4 offset identity survived a line-by-line reading and every computation below; every number in the receipts that I
recomputed (629 / 1,902 / 1,143 / 759 / 186 = 69 + 73 + 44 / 340 / 544,115 / 320,748 / P6 lengths) is exact. Two
statements typed PROVED are false as written and need one-line repairs (M1, M2); neither touches T1-T3.

## MAJOR

**M1. A false example inside T4: "every J-born source over a strong source" does not have d_v - E_v = 1.**
Location: `PROOF.md` section 7, lines 271-272 ("Every left source with d_v = E_v + 1 (for example t, O t, and every
J-born source over a strong source) yields all odd offsets"); THM-4883, T4, line 64 ("every odd offset when
d_v - E_v = 1, e.g. v = t, O t, any J-born source over a strong source").
Problem: J(D, E, U, d) = (D - E + U, D - E, D + U, D) (BILINEAR_C2 (9), DIGIT_BLOCK), so a J-born source has
d - E = E(parent), which is 1 only when the parent has E = 1. Strong sources have E >= 1, not E = 1.
Counterexample: the strong source t.JJZ = (9, 3, 9, 3) (literal 5-vertex tournament, strongly connected, one
root-start path, profile by subset DP) has the J-child v = t.JJZJ = (15, 6, 18, 9) (literal 6-vertex profile, one
root-start path), d_v - E_v = 3, alpha = 69, beta = 66, c = 24, gamma = gcd(69, 3) = 3. Every offset
M(v, Jq) - 2N(v, Oq) is = 3 (mod 6) (all q with odd n < 2^20; this is T4's own class), so neither -1 nor +1 occurs.
Literal C2 graphs built from BILINEAR_C2's text confirm it: q = one-apex 3 and 5 give parents of orders 11, 12 with
H = 417, 711, children of orders 12, 13 with H = 867, 1461, exactly two Hamiltonian cycles each, offsets 33 and 39.
Similarly t.JZJZ = (15, 5, 15, 5) (strong, literal) has a J-child with gamma = 5.
Repair: "e.g. v = t, O t, any J-born source over a source with E = 1 (for instance over a strong unit source), any
O-born source over a source with d = 1". The main statement of T4 (one class modulo 2 gcd(alpha, d_v - E_v)) is
correct and unaffected.

**M2. T5 is stated for every reduced slope, but its "iff" holds only for 1/2 < a/b < 1.**
Location: `PROOF.md` section 8, T5 (lines 277-283: "Fix a reduced slope a/b ... independent of the source iff b is
odd, ord_b(2) divides k and w = W^(k/ord_b(2)) (Session 198's clock)"); README line 21-22 ("the only words that move
aU - bD by a source-independent amount are the powers of Session 198's clock"); THM-4883 title clause (T5) and T5
(lines 67-70).
Problem: Session 198's clock W (Lemma 6) exists only for 1/2 < a/b < 1. For other slopes with odd b there are
source-independent words that are not Session 198 clocks, and for slopes above 1 the condition "b odd, ord_b(2) | k"
is not sufficient. The proof's own computation (lambda_w = a(2^k - 1) - b A_w = 0) needs 0 <= A_w <= 2^k - 1, i.e.
0 <= a/b <= 1, which the statement omits.
Counterexamples (constant change on every strong head; checked on heads 3-13): a/b = 1/3, w = ZO changes aU - bD by
+1 (for odd m, F(4m - 2) = 4m - 2 - 3(g(m) + m - 1) = F(m) + 1); 2/5: OZZO (-2); 1/5: ZZOO (+3); 3/7: ZZO (+3);
1/7: ZOO (+3); 2/9: OZZZOO (-2); 4/9: OOZZZO (-2); degenerate 0/1: O and 1/1: Z (change 0). For 4/3 and 6/5 (b odd,
ord_b(2) | k possible) no word of length <= 10 works.
Repair: either restrict T5 to 1/2 < a/b < 1, or state: the change is source-independent iff b A_w = a(2^k - 1), i.e.
iff b is odd, ord_b(2) | k, 0 <= a/b <= 1 and w is the (k/ord)-th power of the ord-letter word of a(2^ord - 1)/b,
which is Session 198's clock exactly when 1/2 < a/b < 1. Adjust the README line and the THM title accordingly.
(The interior statement is correct: 153 slopes with b <= 31, all 5,013,198 words of length <= 14, no exception.)

## MINOR

1. **HYP-9476 sub-labels disagree.** The HYP file uses (A) predictions, (B) density, (C) witness size, but
   `PROOF.md` line 288 says "HYP-9476 (B)" for witness size, line 296 "HYP-9476 (A) records their densities",
   README lines 30-31 "(A) density, (B) witness size", and the frozen `PREREGISTRATION.md` lines 61/67 "(A)", "(B)".
   Repair: relabel PROOF/README to (B)/(C) and add one sentence in HYP-9476 mapping the pre-registration's (A)/(B).
2. **Redundant maximum in c_+.** For s >= 2, max(5(a - b), F(2^s - 1)) = F(2^s - 1) always: F(2^(j+1) - 1) =
   2^j(2a - b) - a - b increases with j and F(3) = 3(a - b) > 5(a - b). (All 18 positive-exact cases with s >= 2 in
   b <= 64.) Locations: `PROOF.md` line 223, THM-4883 line 54 (the frozen pre-registration may keep it). Repair:
   c_+ = 5(a - b) (s = 1), c_+ = F(2^s - 1) = 2^(s-1)(2a - b) - a - b (s >= 2).
3. **THM-4883 proof summary (lines 43-46).** "exactly one branch has an interior slope" is false in the dyadic case
   b' = 1, a = b/2 + 1 (3/4, 5/8, 9/16, ...), where no branch is interior (one has slope exactly 1/2; `PROOF.md`
   handles this as the no-j* case). "one-signed by the cone 2 g(n) >= n + 3" omits g(n) <= n, which gives the
   slope >= 1 side. Repair: "at most one branch is interior (none iff b' = 1 and a = b/2 + 1); the others are
   one-signed by (n + 3)/2 <= g(n) <= n".
4. **T4 uses E_v <= d_v without saying so.** The dichotomy "E_v < d_v / E_v = d_v" needs E <= d, which BILINEAR_C2
   does not list (it proves 1 <= d <= D <= U and 0 <= E <= D). One line: a Hamilton path ending at c restricts
   injectively to a Hamilton path of V - c. (0 of 1,950 bank sources violate it.)
5. **Missing credit for Lemma 3 and Lemma 2.** Lemma 3 is THM-4882 (2) (the branch formulas for n = 2m - 1 with m odd
   and m = 2m' even) iterated with g(2^(j-1) m') = g(m'); cite it. The sign statements of Lemma 2 (a)/(b) are Session
   198 section 4 and THM-4882 (3) "outside the interval"; the attained extrema (alpha - beta at n = 1, 3(alpha - beta)
   at n = 3) are the new part.
6. **`PROOF.md` lines 132-133:** "for even b it is one full class modulo 2^(s+1) plus one-side-bounded remainders" is
   false for b = 2^s (no full class). Repair: "for even b that is not a power of two".
7. **Stale status of THM-4882 (`PROOF.md` line 10, "audit pending at the time of writing").** THM-4882's audit (PASS
   WITH CORRECTIONS, MISTAKE-722, commit 5c029326f0, 15:08:18) precedes this packet's receipts (c4e9654f72,
   15:24:44), and section 6 of the same file already cites it. Update the sentence.
8. **HYP-9476 (B), lines 33-34:** "If true, THM-4883 is the complete list of full residue classes and S(a, b) is not
   eventually periodic for even b." The first clause does not need (B): by T1 every class other than C* is one-side
   bounded, so every two-sided infinite arithmetic progression (any modulus) contained in S(a, b) already lies in C*.
   The second clause does not follow from (B): a class can meet S in density strictly between 0 and 1 and still be
   eventually periodic with a finer period (e.g. half of its subclasses modulo 2^(s+2)). Repair: replace by the T1
   corollary (Upgrade 4) and "(B) rules out that a non-full class is cofinite on its unbounded side".
9. **P6 (b) at j = 9 for 8/9 was not decided by the depth-34 search.** `PROOF.md` lines 326-327 list "at j = 9 four
   targets need more than 34 letters" under "(b) fails", and the receipt's `b_failing_j` is [4, 9]; but 3j + 8 = 35
   at j = 9, so a target first found at length 35 would satisfy (b) (`verify.py` counts any target missing at depth
   34 as a failure). My own depth-35 search settles it: targets -427, -425, -405 and 277 still have no witness of
   length <= 35, so (b) does fail at j = 9 as well. The REFUTED typing stands (it rests on j = 4: maxlen 21 > 20 at
   r = -15, decided, reproduced); fix the stated reason ("need at least 36 letters, depth-35 search").
10. **T5 remark on HYP-9475 (C) (`PROOF.md` lines 285-289, THM-4883 line 69).** T5 excludes a second
    source-independent move for the same functional; another slope's clock still changes aU - bD by
    lambda_w U + const, which is target-dependent but usable in a search, so "unavailable there" overstates. "Every
    witness of a target r is a prefix P followed by t >= 0 clock periods" is vacuous (t = 0). Soften both.
11. **"Constructive" witnesses (`PROOF.md` lines 214-216).** "tail T after a strong head chosen by Theorem A0's recipe"
    describes the proof; `verify.py` finds heads by direct search over all odd m < 2^20 in the right class modulo K,
    not by A0's branch heads 2^(j+1) L - 1 with gadget-built L. A0's heads are explicit but impractically long when K
    is large (the gadget word has O(K^2) bits, hence about 2^(O(K^2))/K clock periods). With L restricted to < 2^20 my
    implementation of A0's branch-head route stayed under a 40,000-letter cap for 923 of 974 targets (51 failures,
    K = 50, 62, 74, 98) and 467 of 506 T2 heads; every target was then witnessed by the direct-head recipe. Say "a
    head with the right residue modulo K (in practice by direct search), then clock repetitions".
12. **Naming.** "The automaton" (`PROOF.md` line 185) is the proof's recursion, not a finite automaton; it invites
    confusion with THM-4878's cut-state automaton. Suggest "classification algorithm".
13. **Post-registration tie-break in P6 (c).** Clause (c) depends on which minimal witness is used when several have
    the same length; the rule (first in DFS order, O before Z) is fixed only in `witness_length.c` (written after the
    pre-registration commit). Not outcome-changing (184/522 and 160/433 against 1/2), but record it as a choice made
    after registration.
14. **Informational.** Commit c4e9654f72 also folds the unrelated HYP-9467 and HYP-9464 lines of
    `05-knowledge/hypotheses/INDEX.md` into one line (content preserved). Mention it in the commit message or leave
    it; it is not an error.

## UPGRADES

1. **T4, explicit class and a parity lemma.** When E_v < d_v the offsets fill exactly the class
   c_v + d_v - E_v = U_v - D_v + E_v + 2 d_v (mod 2 gamma) (from S(beta', alpha') = all integers of parity
   (E_v - d_v)/gamma). U - D + E is odd for every physical source (t: 1; O gives 2U - D, Z leaves it unchanged, J gives
   D; D is always odd), so the class is odd, consistent with Redei for m - 2n. In the case E_v = d_v the bound c_v is
   attained (q with U = D, i.e. n = 2^k + 1; e.g. v = J t, q = 3: offset 3 = c_v, literal graphs H = 77 and 157).
   Checked on all 1,950 bank sources and 15,972,450 (v, q) pairs.
2. **The nine "lower-bound-only" classes of P5 are decided.** For the slope-1 piece {n - g(n)}, every value <= X is
   attained at some odd n <= 4X + 3: by Lemma 3 with alpha = beta, odd m gives F(2m - 1) = F(m) (a smaller n), and
   m = 2^j m' gives F = 2^j m' - g(m') - 1 >= (2^j - 1) m' - 1. For the slope-1/2 piece {n - 2 g(n)} the values are
   1 - 2g(m) (m odd) and -1 - 2g(m') (the same for every j), so every value with |F| <= X is attained with
   m <= X + 1 and j = 1, i.e. at n <= 4X + 3. The window is therefore decided by n of order 4X times the scale of
   the piece (far below 2^22 for X = 20000). Hence the classes reported
   as lower bounds (0.526-0.546) are exact window fractions, and P5 / HYP-9476 (B) cover all 59 classes with range
   0.097-0.617. Data: all 59 classes identical at n < 2^22 and n < 2^26; the raw pieces give 0.510102 and 0.517538
   already at n <= 4X + 8 and unchanged to 2^27.
3. **HYP-9476 (B) is stated for every slope; add exterior evidence.** For 15 exterior slopes (2/1, 3/1, 3/2, 5/2,
   4/3, 5/3, 5/4, 1/1, 0/1, 1/3, 1/4, 2/5, 1/2, 3/8, 9/8) every class has window fraction in [0.190, 0.601] on the
   unbounded side (500 <= |r| <= 20000), unchanged from 2^22 to 2^26.
4. **A corollary worth stating in T1:** every two-sided infinite arithmetic progression contained in S(a, b) lies in
   C* (none exists when b is a power of two or the slope is exterior).
5. **Dyadic example.** For b = 2^s the proof gives every class's sign and bound; e.g. 3/4 modulo 8: classes 1 and 5
   bounded below by 1 and -3, classes 3 and 7 bounded above by -5 (attained in class 3; class 7 has maximum -9). This
   is the modulo-8 refinement of THM-4882's two-sided modulo-4 confinement at (3, 4); worth one line in T1.
6. **Dictionary, exhaustive small-order check.** All labelled tournaments of orders 2-6 rooted at a fixed vertex: 2,141
   have s = d = 1, every one with profile (g(U), U mod 2, U, 1); the 1,070 with E = 1 are exactly the strongly
   connected ones and have odd U >= 3. This can be cited next to the KNOWN dictionary.

## What I verified and could not break

Proof reading (line by line): Lemma 1; Lemma 2 (a)-(c) including g(2^k - 1) = 2^(k-1) + 1, g(2^k + 1) = 2^k + 1;
Lemma 3 (both branches, exhaustive pairing of (m odd >= 3) and (j >= 1, m' odd >= 1)); Lemma 5 (peeling gives
A_T = sum eps_i 2^i, kappa independent of the odd head); Lemma 6 (U_i = 2^i(U - u*) + u_i, E_i = eps_(i-1) cyclically on
a strong input because the last letter is Z, F' = F - beta Delta, run contribution v + (1 - v) 2^(-j) > 0 with the
wrap-around run, K integral and even); Theorem A0 against Session 197 Theorem 1 (even modulus, any u, odd z, arbitrarily
large L: exactly what A0 uses; gcd(sigma, beta) = gcd(2 alpha - beta, beta) = 1; sigma y - beta x reaches every residue
modulo K/2); the routing table, every row for s = 1 and s >= 2 (gcd computations, parities, v_2(r + a)); the group
sizes 2^(s - nu); the sign analysis (j < j* below 1/2, j > j* above 1, the no-j* case exactly b' = 1, a_1 = 1,
j_0 = s - 2); cases (i)-(iv), including the bijection of P_(s'') onto the classes of G_(j*+1) under phi and the
harmlessness of the S* point; T2 (A_T = floor(2^(L-1) a/b), the head slope has denominator exactly b' and lies in
(1/2, 1)); T3 (level maps r -> 2^i r + (b - a)(2^i - 1) because b_i - a_i = b - a at every level, minima
(2^(i+1) + 1)(a - b) at n = 2^(i+1) + 1, maxima 5(a - b) and F(2^(j+1) - 1), both converses); T4 (expansion of
BILINEAR_C2 (11)-(12), alpha odd, alpha >= 3D_v >= 3d_v > 2(d_v - E_v) from the cone proved in BILINEAR_C2 section 1);
T5 inside the interval.

Computations (all in `INDEPENDENT.json`):

- g: recurrence, set-bit expansion and O/Z transfer agree on n < 2^14 and on 24 integers of up to 1,500 bits;
  Lemma 1 holds for all n < 2^24.
- Dictionary: 358 literal one-apex tournaments (orders <= 14): (D, E, U, d, s, cycles) = (g(n), n mod 2, n, 1, 1,
  n mod 2); exhaustive converse at orders 2-6 (Upgrade 6). Lemma 3 multiset identity on 8 slopes including
  alpha <= 0.
- Lemma 6 on all 1,053 reduced slopes with odd b <= 101: last letter Z, mixed, first letter O iff alpha even,
  K = beta Delta exactly (rationals), run decomposition exact, K > 0 even, the same drop on 12 strong heads up to
  3^50; K ranges over [2, 1250].
- Theorem A0 by its own route: my implementation of the gadget construction gives the prescribed joint residues for
  every even M <= 40, every u and every odd z (5,740 cases); 974 targets (|r| <= 12,345) on 106 slopes with odd
  b <= 31 witnessed, 923 through the branch heads n0 = 2^(j+1) L - 1 (80 with the gadget's own L) and the other 51 by
  the direct-head recipe (Minor 11), all verified by three g evaluations.
- Census of all 629 interior slopes with b <= 64 over odd n < 2^24: exactly one class with values beyond +-2^18 on
  both sides iff b' > 1, and it is T2's class (598/598; none for the 31 dyadic slopes); my own recursion marks exactly
  that class FULL; 0 bound violations in the 1,902 non-full classes; 1,143 classes miss their group bound (first:
  3/4, class 7, bound -5, maximum -9) and all 759 groups attain it (the same 759 whether groups are formed by
  (sign, bound) label or by identical piece sets); T2 suffix: 398,982,982 odd n < 2^24 carry the suffix, 0 lie
  outside C*; T3: 186 slopes, 69 positive-exact, 73 negative-exact, 44 neither, every threshold attained exactly and
  every "neither" case unbounded beyond +-2^18 on both sides.
- P2 window: 340 slopes, 544,115 targets of C* with |r| <= 2000, 320,748 hit directly below 2^24 (identical to the
  receipt); 3,377 sampled unhit targets on 325 slopes (234 with odd b) witnessed by the T2 tail plus a head (467 heads
  by A0's branch route) and verified, including the low-bit condition n - 1 = sum (1 - d_(i+1)) 2^i (mod 2^L); up to
  38,578 bits.
- 66 further interior slopes with b in {96, ..., 256} (including 129/256, 255/256, 65/128, 127/128, 97/192, 191/192)
  and 20 exterior slopes (including a <= 0) over n < 2^22: structure, bounds and T3 exactly as stated; exterior S is
  one-signed with extremum 3(a - b) and no two-sided class.
- Worked examples: 5/6 -> 3 (mod 4), suffix ZZ; 7/12 -> 5 (mod 8), ZOOZ; 11/12 -> 7 (mod 8), ZZZ; outside C*:
  max S(7, 10) = -15, min S(5, 6) = -3, max S(7, 12) = -15.
- T4: symbolic identity (sympy); literal Jq profiles; 1,950 bank sources (U <= 300, all literal); 362 of them have
  gamma > 1 (gamma in {3, 5, 7, 9, ..., 53}); literal C2 pairs of orders 8-19 from v = t (157/313, 187/375), O t
  (711/1421, 77309/154619), OJO = (5, 0, 10, 1) (6683/13367), J t (E = d: offsets -23 and 3 <= c = 3) and JJZJ
  (gamma = 3: offsets 33, 39), each with H equal to the bilinear form and exactly two Hamiltonian cycles.
- P5: the 59 classes of the 14 slopes are identical at n < 2^22 and n < 2^26 (Upgrade 2).
- P6 (own DFS): 3/5 maxlen 9, 17, 22, 22, 32, 32 at j = 4..9, all found by depth 32; 8/9 21, 21, 21, 29, 32 at
  j = 4..8; depth 35: four targets of |r| <= 512 need >= 36 letters (Minor 9).
- Integrity: `dependencies/*.py` byte-identical to Session 198's (SHA256 33cfe6ce..., dc32dc1e...);
  `VALIDATION_NORMAL.json` = `VALIDATION_OPTIMIZED.json`.

Pre-registration: reservation 6cb451411f (14:57:04); `PREREGISTRATION.md` alone in 2c11d95d56 (15:04:07) and
unchanged since; `verify.py`, `witness_length.c` and all receipts first appear in c4e9654f72 (15:24:44). The
exploratory P1 run and the 2/3, 4/5, 6/7 witness searches are disclosed. Outcomes are reported honestly (P1 per-class
attainment and P6 for 8/9 REFUTED). The post-registration choices (P4's search, MQ = 2^12, words of length <= 5 ending
in O, child order <= 19; P5's decided-class rule; P6's tie-break) cannot manufacture a confirmation: P4's pairs are
literal graphs, P5's excluded classes are decided anyway and lie in the same range (Upgrade 2), and P6 (c) has
margin (Minor 13).

Novelty and credit: no packet in session197-199, `05-knowledge/results` or `01-canon/theorems` contains the 2-adic
classification for even b, the digit formula for C*, the half-line thresholds, the offset formula for a general left
source or clock uniqueness (searched for "2-adic", "v_2", "2^(s+1)", "offset", "aU-bD", "signed response"; read
Session 198 rational-response-coverage, source-digit-blocks, upper-replacement-fibers, primitive-branching-barrier;
Session 199 hidden-structure, physical-interpolation; Session 197 ONE_CORE_BINARY_INVERSE; THM-4882 and its note as
amended after its audit). The credits for the dictionary, the block transfer (Lemma 5), the clock, the gadget theorem,
BILINEAR_C2, THM-4882 (3)(iii) and its -3(b - a) threshold are present and accurate; Minor 5 lists the two missing ones.
