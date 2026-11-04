# Audit of THM-4872 / HYP-9465 (opus-2026-11-04-S76)

**Scope.** Under audit:
- the note `05-knowledge/results/unit_start_pm2_port_pairs_symmetry_split_mrb_20261104.md` (cited as N);
- the canon `01-canon/theorems/THM-4872-unit-start-plus-minus-two-port-pairs-and-the-z3-trivial-symmetry-split-of-three-cycle-tournaments.md` (T);
- the pre-registration `05-knowledge/hypotheses/HYP-9465-...md` (H);
- the four session programs and their outputs.

All code in this folder was written for the audit; none of the session's programs was reused.

## Verdict: PASS WITH CORRECTIONS

Every number in THM-4872 (A) and (B), in HYP-9465 P1-P3 and in the MRB section was reproduced. The methods differ
from the session's: Burnside counting instead of canonical forms, a root-prepending generator instead of the F1
component grammar, and an own canonical form with a separate automorphism counter.

The proof of (B) is correct, and no mathematical claim fails. There are 2 MAJOR findings, both about provenance and
overlap (no mathematical error), and 13 MINOR ones. One minor finding is a factual slip: 47 should be 91.

The pre-registration is genuine. HYP-9465 was committed in 9036e4dcb9 at 15:23:12, and the order-13/14 outputs first
appear in 3192fab75f at 15:26:48. The programs did not change between the two commits.

## 1. What was reproduced

**(A) at order 7** (`aud_order7.c` → `aud_order7.out`, 4 s).
- **Burnside count.** Over the 225 odd-order elements of S_6 that fix the root there are 2,188,800 fixed labelled
  tournaments, and every class count is integral: 3,040 rooted classes, 2,395 strong, 55 strong with s = 1, and 89
  with s = 1 in all. F by the subset formula agrees with brute force over all 7!·7 cut orderings on 3,000 random
  tournaments.
- **Pairs** of strong classes with h' = h + 2 and F' = F:
  - **689** in all;
  - **1** with s = s' = 1, namely (h,s,q,F,e) = (29,1,15,157,10) and (31,1,9,157,9);
  - **47** with s' = s and q' = q, none with s = 1 (23 of them also have e' = e);
  - **91** with q' = q alone (m1).
- **A and B** have the stated profiles and 720 labelled copies each, so their rooted automorphism groups are trivial.
- **Routes.** BFS over the 2^21 labelled tournaments, root fixed, strong intermediates:
  - distance **4** in both directions;
  - **38** shortest routes into 5 copies, by a dynamic programme and by explicit enumeration;
  - **0** routes keep s = 1;
  - **2** routes per direction keep the wrappers below the start, in the closed form and literally for m = 1..60
    (128-bit); these are the note's arc orders;
  - distance **7** if every intermediate keeps s = 1, and still 4 if strongness is dropped.
- **Arc sets.** `aud_route_sets.py` (pure Python, relabelling B over 720 permutations) shows the 38 routes are the
  strong orderings (9 + 4 + 12 + 9 + 4) of the 5 four-arc sets at Hamming distance 4. Both below-start routes order
  ONE set: {01,03,12,34} from A to B, and {01,04,12,34} from B to A.
- **Wrapper.** The literal subset-DP count equals h + 2^m F on 4,500 random cases (294 non-strong, 1,453 with s ≠ 1)
  and on all 456 route intermediates. A and B give 343/345, 657/659, 1285/1287, 2541/2543 and 5053/5055 for m = 1..5.

**(A) at orders 3-14** (`aud_unitstart.c` → `aud_unitstart_14.out`, 13 s).
- **Normal form and generator.** A rooted class with s = 1 is exactly one labelled tournament whose unique path from 0
  is 0, 1, ..., n-1. Since s_T(0) = Σ over j in N+(0) of s_(T-0)(j), the new root must beat the old root and may beat
  any subset of the vertices with s = 0. This does not use F1's component decomposition. At n <= 8 the generated set
  is identical to a brute-force scan.
- **Counts.** All unit-start classes number F(2n-3), up to 75,025; the strong ones number F(2n-4), up to 46,368 (P2).
- **Pairs.**
  - The ±2 pairs at orders 7-14 are **1, 2, 5, 1, 6, 10, 4, 8**.
  - The distinct F counts agree: 48, 126, 323, 868, 2242, 5880, 15415, 40639.
  - q-equal pairs first appear at order 14 (2 of them).
  - The smallest pairs run from (29,31) to (1959,1961), and the first order-9 pair in F order is 105/107 (F = 573).
  - All pair rows are printed.

**(B)** (`aud_fewcycle.c` → `aud_fewcycle_14.out`, 23 s).
- **Method.** Generation is complete by Moon's theorem and the insertion injection. Each parent's Hamiltonian-path
  matrix comes from a pruned DFS, and a Gray code runs over the in-sets.
- **Canonical form and |Aut|.** An own individualisation-refinement canonical form, with |Aut| = the number of optimal
  leaves. A separate backtracking counter recounts |Aut| on every sole-exit-free class with hc = 3, 4, 5:
  0 mismatches. hc is recounted by DP on every class to n = 12 and on all u_3-u_5 classes to n = 14: 0 mismatches.
- **Class totals.** The classes with hc <= 5 number 1, 1, 6, 35, 140, 468, 1506, 4708, 14390, 43225, 128077, 375339
  at n = 3-14. At n = 5-10, t_1-t_5 equal the S75 audit's exhaustive census of all classes (`aud_census_exh_10.out`).
- **The u_3 split**, as (asymmetric, Z_3), is (F(2n-11), F(2n-10)) at every n = 6-14:
  (1,1), (2,3), (5,8), (13,21), (34,55), (89,144), (233,377), (610,987), (1597,2584).
- **hc = 4 and 5.** Every u_4 and u_5 class is asymmetric at n = 7-14. u_4 is already rigid at n = 6, but there 2 of
  the 4 hc = 5 classes have |Aut| = 3.
- **The proof** was checked step by step.
  - B1 fixes D, u and l.
  - The prefix stabiliser is trivial, by THM-4868 (E). More directly: an automorphism of a uniquely Hamiltonian
    tournament rotates its unique cycle, so fixing one vertex fixes all.
  - The chain is fixed.
  - Aut(T) is isomorphic to Aut(W, A): the prefix and chain vertices dominate W, and u sees W only through A.
  - A C3-block core of order d has 2^(d-5) choices: u, a chain of d - 4 vertices, and d - 5 free bits.
  - x^5/(1-2x) · x(1-2x)/(1-3x+x^2) = x^6/(1-3x+x^2), with coefficients F(2n-10).

**MRB** (`aud_mrb.py` → `aud_mrb.out`).
- **The constant and the integral.** C_MRB agrees to 60 digits by two methods, and the Abel-Plana integral agrees to
  9e-48.
- **Theorem 1.** M_0 = 1/4 and M_1 = 1/8. h has exact Taylor coefficients 1, -1/2, -3/4, 9/8, -271/720, -523/630.
  R = 3.596424620671202485179340542732e-4, equal to the remainder integral. The residues are -1/16 and 1/64.
- **Lambert W.** N(a_n) = n for n = 3..20 (n = 2 returns 4). The mirrors are r_3 = 2.47805, r_4 = 2, r_5 = 1.76492
  and r_8 = 1.46250.
- **Euler table.** The transform reproduces the attachment's table digit for digit at M = 51, 101, 151, 201.
- **Numerology.** 2 - 29/16 = 3/16, and -1/4 → -7/4 → 5/4 → -1/4 under x^2 - 29/16.

**The §4 table.** Rows 1, 2, 3, 4 and 6 match the THM-4870 note (§0, §1/(A), §5, (C), (B)). Rows 5 and 7 need
corrections (m7, m9).

## 2. MAJOR findings (provenance and overlap; no mathematical error)

**M1. A concurrent same-day answer to the same prompt is missed.**
- **Where.** N lines 29-40 credit only kind-pasteur. The same gap is in §4 (lines 187-197), in N lines 155-156 and in
  `05-knowledge/results/INDEX.md:149`.
- **What was missed.** Session190 had already published two checkpoints; S76's first commit 9036e4dcb9 has CP2 as its
  parent.
  - CP1 (3f218bec77, 15:08): the synthesis of the four-colour, hypohamiltonian and Thomassen papers, with a
    boundary-state patch transfer.
  - CP2 (a8108cfdb2, 15:22), `session190_constructive_links_20261104/`:
    - `planar/`: planar interfaces and 124 square-sum deletion cards;
    - `bijection/`: tournament/even-graph Walsh coupling, with no natural bijection;
    - `groups/session190-groups/PROOF.md`: SO(3) against planar groups and Maschke; P_7 with F_21 and AGL(1,7), and
      7/21/42; Theorem 3, the signed metagraph-edge stabiliser is 1 or C_2.
- **Fix.**
  - Credit Session190 in "What was already done" and in the INDEX line.
  - Add Session190 rows to §4.
  - At lines 155-156 add "...; Session190 (groups, Theorem 3) also obtains {1, C_2} as signed edge stabilisers".
  - Record a MISTAKE entry: this is the same failure type as MISTAKE-677/694/699.

**M2. The (A) payoff, "a reusable ±2 adjustment across a clocked family", has uncited prior art.**
- **Where.** N lines 90-98 and T line 54.
- **Session176 §1** (`session176_tournament_information_20261104/spectrum/PROOF.md`, 4fde8ff005, 2026-11-03)
  proves a persistent single-edge ±2 gate across a dyadic clock: H = 6·2^r ∓ 1, both strong at every r. Its §2 shows
  why the gate does not recurse into ternary closure.
- **Session185 §5** (`session185_proof_frontiers_20261104/alternative/PROOF.md`, 2026-11-03) has a three-state
  version.
- **Without s = 1, D_m-wrapper ±2 pairs already exist at order 5:** (h, s, F) = (13, 3, 15) and (15, 3, 15), with
  wrappers 43/45, 73/75, ..., 3 arc reversals apart; and 27 profile pairs at order 6 (`aud_small_pairs.out`,
  `aud_order5_pair.out`). So unit start is exactly what makes order 7 first.
- **B is already in the repository.** It is the order-7 nested core, at mark 4, of the named source of 170807:
  `zero_ready_component_returns_20261103/component/COMPONENT.md:141` (2026-11-03), with
  (h, e, P(2), F, S_w(2)) = (31, 9, 644, 157, 31). `aud_lineage_core.py` rebuilds it from the stored rows, and its
  path-order rows equal B's.
- **What is new.** A, the pair, the routes, uniqueness at order 7 and the counts at orders 8-14 were not recorded
  before S76.
- **Fix.** Add "Prior art and scope" paragraphs to N (after line 98) and T (after line 54), and add the zero-ready,
  Session176 and Session185 packets to T's `related`.

## 3. MINOR findings

**m1. N line 76: 47 should be 91.** 47 counts pairs with s' = s AND q' = q; with q' = q alone there are 91.
Fix: "Requiring also s = s' and q = q' leaves 47 pairs (91 with q = q' alone); none has s = s' = 1."

**m2. T line 3: the title omits the order ranges.** "Four or five cycles are asymmetric" is false at n = 6, where 2 of
the 4 hc = 5 classes have Aut Z_3. u_3 = F(2n-9) needs n >= 6, since u_3(5) = 0 ≠ F(1). Fix: add "(n >= 6; four or
five cycles: n >= 7)".

**m3. T line 8, T lines 68-75 and N lines 128-140: the dependency list is incomplete.** The proof also uses:
- THM-4869 (C2), the decomposition;
- B2, that Q is uniquely Hamiltonian with sole exit z;
- B3, the correspondence of classes;
- (C4) with THM-4868 (F), the prefix factor.

Also say that the base of (C1) is a computer census. The rotation argument above can replace the two-case
stabiliser step.

**m4. N lines 47 and 62-65, T line 51: the wrapper formula is checked on a narrower class than it is used on.**
- It was re-verified only on s = 1 strong classes (220 cases), but it is applied to intermediates with s = 2.
- It holds for every rooted U and every m >= 1: if w comes first, Z follows, giving h paths; otherwise r → w and Z
  splits freely between the two slots, giving 2^m F paths. At m = 0 the count is h + F - s.
- The formula is the zero-ready clock (`triangle_return_170807_20261103/ZERO_READY.md`), which Session189 calls
  "inherited".

Fix: state the general proof and credit the zero-ready clock.

**m5. N lines 80-88 and T lines 46-49: the 38 are orderings of 5 arc sets.** The 2 below-start routes per direction are
two orderings of one four-arc set. Fix: say so.

**m6. N lines 76 and 90-94: "matching port response" is only F.** A and B differ in e, q, P(2) (580/644) and S_w(2)
(34/31). Other attachments break the ±2: a sole exit x → 0 gives 53/53, and x beating N+(0) gives 49/73
(`aud_scope_funnel.out`). Fix: restrict "reusable" to D_m. The ±2 also holds at m = 0, where h + F - s gives 185/187.

**m7. N line 150 and §4 row 5 (line 195): wrong source for {2,5,13,34,89} = u_3(6..10).** This is klein's HYP-9459 (2),
also recorded in HYP-9464 (C) and in mac-mini's `cpq_defect_family_mirror_folds_20261104.md` §0. THM-4870 §4, which
row 5 cites, does not contain it. Fix: cite these sources.

**m8. MRB, N lines 158-185.**
- (a) No MRB script was committed, so the 30-40 digit checks cannot be reproduced from the repository. `aud_mrb.py`
  can serve as one.
- (b) Lines 173-174, the cautionary remark.
  - 0.19164 was not reproduced in 36 natural set-ups (`aud_mrb_nsum_variants.out`).
  - On the paired series, the alternating, Levin and Richardson methods give 0.192-0.197; on the unpaired series,
    Euler-Maclaurin gives 0.171 or 0.180.
  - So "its alternating-series and Euler-Maclaurin methods give the correct value" is wrong as written. Only
    alternating-on-unpaired and Euler-Maclaurin-on-paired are right. Fix: reword, or record the exact call.
- (c) Lines 176-182, the mirror reading. It is also in mac-mini's same-day note §5 (14:44) and in the S75 note
  (lines 408-412). Fix: credit both next to session31.
- (d) Observation for the owner: the attachment's table is labelled M = 51..201, but its code prints M = 50..200,
  which differ by about 1e-20. "40 digits" at line 171 is a precision limit; the 200-term value is right to 67 digits.

**m9. §4 row 7 (line 197).** {-3}-self-duality is not an item of this prompt. HYP-9464 is no longer pending:
kind-pasteur confirmed it at 15:34 (9aa6e0836f).

**m10. H lines 28-38 and N lines 212-218: pre-registration typing.**
- P2 is the PROVED F2 count.
- P1 follows from the proof sketch in H lines 19-26.
- P3 is a low-risk existence prediction.
- The order 10-13 counts in H line 31 have no output committed in 9036e4dcb9.

Fix: type P1 and P2 as consistency checks.

**m11. `fewcycle_symmetry_and_unit_pairs_20261104.c` has stale comments.**
- Lines 1-23 are the S75 census header, with stale usage lines.
- Lines 206-208 describe a "Part 2" that the program never runs.
- Line 203 says |Aut| comes from refinement leaves, but `aut_count` is plain backtracking.

**m12. N line 151 and T line 54 / N line 97: wording.**
- "The trivial group takes the smaller part" fails at n = 6, where the parts are 1 + 1.
- "No 3H + 2 / 3H + 4 parent" is asserted without proof. A one-line proof: 29 + 157·2^m ≡ 2 + (-1)^m (mod 3), so the
  lower wrapper is never ≡ 2 (mod 3). Session190 states this in its spectrum packet, §5.1.

**m13. `05-knowledge/results/INDEX.md:149` and `05-knowledge/hypotheses/INDEX.md:2`.** Add Session190 (M1), and change
"2 of the 38" to "2 per direction, one arc set".

## 4. Prior art and overlap

**The git log since 13:00** was read in full (`gitlog_since_1300.txt`, 32 commits including those after the
checkpoint). Before S76 it shows Session190 CP1-CP2 (M1), klein's HYP-9459, kind-pasteur's HYP-9464 (m7, m9) and
mac-mini's THM-4871 note (m8c).

**Packets grepped:** session18x and session190; joint_entry_return_20261104; zero_ready_component_returns_20261103
and zero_ready_strong_switch_20261103; the *_returns_20261103 adapters. The search terms were the 157/29/31 tuples,
four-arc/four-edge, ±2/"differ by two", and 512423/512425. The hits are B (M2) and the earlier ±2 gates (M2). The pair
(A, B), its routes and the counts at orders 8-14 are recorded nowhere before S76.

**After S76's checkpoint.** This is not an S76 error, but it must be reconciled at close-out. Session190 CP3
(3cb7969b5f, 15:34), `session190_constructive_links_20261104/spectrum/`:
- builds the same pair and route independently: root 2, A0 = (6,28,8,1,109,15,47), equal to N's second below-start
  route under an exact rooted relabelling;
- credits S76;
- adds the residue argument, an obstruction at orders <= 6 that agrees with this audit, and fixed-parent exclusions
  for S76's order-11 pairs 575/577 (F = 3819) and 605/607 (F = 15909).

S77 reserved THM-4873 (the ternary parent) at 15:56.

**The symmetry split (B)** appears nowhere else. The search covered F(2n-10)/F(2n-11), "C3 block", and the |Aut|
statements in THM-4869, HYP-9459/9461/9464, the S75 note and audit, the THM-4870 note, mac-mini's note and Session190
groups. It is new, as a corollary of THM-4869. HYP-9464 (B), which counts uniquely Hamiltonian deletions, is a
different invariant.

**MRB.** "3/16 = 2 - 29/16" appears nowhere else, so NUMEROLOGY is the right typing. Grogan was considered in S73, as N
says.

## 5. Reproduction (gcc 15.2 -O3 -march=native -fopenmp, 12 threads; Python 3.13, mpmath 1.3, sympy 1.14)

    gcc -O3 -march=native -fopenmp -o aud_order7 aud_order7.c && ./aud_order7 > aud_order7.out
    gcc -O3 -march=native -fopenmp -o aud_unitstart aud_unitstart.c && ./aud_unitstart 14 > aud_unitstart_14.out
    gcc -O3 -march=native -fopenmp -o aud_fewcycle aud_fewcycle.c && ./aud_fewcycle 14 12 > aud_fewcycle_14.out
    python3 aud_route_sets.py; python3 aud_lineage_core.py; python3 aud_small_pairs.py; python3 aud_order5_pair.py
    python3 aud_scope_funnel.py; python3 aud_mrb.py; python3 aud_mrb_nsum_variants.py

## 6. Summary

PASS WITH CORRECTIONS.
- **Theorem A.** Every pasted claim and every number in THM-4872 (A) through order 14 checks out: the pair, the
  distances 4 and 7, the 38 routes, the 2 below-start routes per direction, the wrappers, and the counts
  1, 2, 5, 1, 6, 10, 4, 8.
- **Theorem B** is proved correctly and reproduced at n = 6-14.
- **MRB.** The checks hold, except the cautionary remark.
- **Fix before close-out:**
  - the missed same-day Session190 answer (M1);
  - the uncited ±2 prior art: Session176/185, the order-5/6 wrapper pairs, and B as the 170807 core (M2);
  - the 47/91 slip (m1);
  - the ranges missing from T's title (m2).
- Then reconcile with Session190 CP3.
