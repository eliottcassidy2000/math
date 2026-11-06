# Independent audit of opus-2026-11-04-S78 (THM-4875, HYP-9468)

**Auditor:** independent adversarial audit subagent, 2026-11-04.

**Audited version.**
- Pre-audit checkpoint 61397dd494.
- The two later S78 commits: 94805064e9 (2^32 result filled in) and 78c462fb63 (residue profile).
- Note line numbers refer to the current HEAD file `05-knowledge/results/old_mark_bound_lorentzian_gadgets_one_apex_join_coverage_20261104.md` (sections 0-3 are unchanged since the checkpoint).
- THM line numbers refer to `01-canon/theorems/THM-4875-...md` at HEAD.

**Scope rule.** Nothing outside this folder was modified. Executables were deleted after use.

## Verdict: PASS WITH CORRECTIONS

**All four mathematical claims are correct.** I reproduced each with independent code and a fresh derivation:
- (A) the old-mark bounds;
- (B) the transfer law and the adjoint;
- (C) the descent lemma itself;
- (D) the 2^32 coverage.

The 2^32 coverage was rerun with an independently derived g, an independently generated bank and an independent sieve. The misses are exactly {7, 21}, and all 22 per-block F1 miss counts agree exactly.

**The pre-registration is sound.** Both pre-registered outputs postdate commit 19562af223.

**Three major corrections are needed:**
- (B)'s general adjoint, and the general transfer, were proved on 29 October and are not cited.
- (A) claims as new a bound already in canon (THM-4873 (ii)).
- (C)'s stated consequence for Chan-Pak is not proved, and as worded it concerns a false target.

**Upgrade.** The pre-registered conjecture "non-sink old digits >= 9" (HYP-9468 P3) is PROVED below at every order. Its equality case is exactly (B, a, w) = (TT_n, a sink, w penultimate).

## 1. What I reproduced, and how

All code and outputs are in this folder. My programs share no code with the session's programs.

| claim | auditor's method | result |
|---|---|---|
| order-<=7 rooted-triple bank | `aud_bank.c`: all 2^(n(n-1)/2) labelled tournaments of each order n = 1..7 (2,131,019 in total), rooted at vertex 0, pull-style DP | 2,732 distinct triples, identical as a set to `rooted_triple_bank_order7_20261104.txt` and to the bank the 2^32 run read (`C:/tmp/s78/bank7.txt`); only CRLF differs |
| one-apex triple (g(v), [v even], 1+v) | `aud_witness.c` (1): direct subset-DP counts of T(w) for all 8,191 words of length 0..12, leading zeros included | 0 failures |
| (D) coverage | `aud_cover.c`: own g (closed form below), own 2^16 split, byte map with relaxed atomic stores, (row, v-chunk) work items | misses exactly {7, 21} up to 2^24, 2^28 and **2^32** (`aud_cover_32.out`, 940 s on 12 threads; 47,964,383,242 F1 pairs, 134,010,742,236 F2 evaluations); the F1-only miss counts for k = 10..31 equal the session's `.out` exactly |
| F0 needed? | `aud_cover.c -DNOF0` to 2^28 | F1 + F2 alone miss exactly {7, 21} (upgrade U3) |
| explicit witnesses | `aud_witness.c` (2): find a representation, **build** the join tournament, count its Hamiltonian paths by direct subset DP | 41 witnesses, orders 14-26, N up to 40,376,785 (32 F1, 9 F2): 41 OK, 0 mismatches (`aud_witness.out`, `aud_witness_large.out`) |
| (A1)-(A5), sink formula, digit-4 equality, non-sink >= 9 | `aud_oldmark.c`: own class generation by canonical augmentation (ordered invariant partition, minimum code over cell-respecting permutations); class counts 1, 1, 2, 4, 12, 56, 456, 6,880 = A000568 | every check, at every order 1-8, has 0 failures (`aud_oldmark_8.out`) |
| (B) transfer formula | `aud_gadget.c`: M **fitted** from three probes by exact Cramer (det 2) for all 26,705 marked gadgets of order <= 5 | fitted M integral and equal to the auditor's own implementation of the note's closed formula; linearity confirmed on 456,745 further base checks (89 bases per gadget of order <= 4, namely all 75 labelled rooted bases of order <= 4 plus 14 of order 5-6; the 14 for order 5) |
| (B) adjoint | same, fitted matrices for G and G^op (order <= 5); closed formula at order 6 | 0 failures on 26,705 + 1,179,648 marked gadgets |
| (C) inequality and equality | `aud_descent.c`: every gadget of order <= 5 against every one of the 2,732 realisable triples | 0 violations of h' >= R >= h (72,958,060 pairs); the 13,661 cases with h' = h are exactly the characterisation in M3(v) |
| exploratory claims | the session's `greedy_bank.c` rerun (`greedy_repro_20.out`); `aud_f1res.c` with the auditor's g | greedy needs exactly 111 triples; F1-miss residues on (2^22, 2^24] are mod 3 = 0.0036 / 0.0124 / 0.0039, matching the note |
| session's own (A) program | `old_mark_digit_bound_20261104.c` rebuilt and rerun (`omb_rerun_8.out`) | output identical to the committed `.out`; runtime 0.75 s |

**Pre-registration timing** (file times from `C:/tmp/s78`, all -0700):
- HYP-9468 was committed at 17:45:23 (19562af223).
- P3: `omb8.out` was last written at 17:45:41. The program needs 0.75 s for orders 1-8, so the order-8 run started after the commit. The order-8 line is written only at the end.
- P1 and P2: `cover.exe` was built at 17:46:08, `cover26.out` written at 17:46:22, and `cover32.out` at 18:08:01.
- So both pre-registered outputs postdate the commit.

**Concurrent work.** Session 193 (origin/main 18:44-18:45, after the S78 checkpoint) checked the S78 families to 2^20 and calls the 2^32 result "SOURCE-REPORTED and UNAUDITED". This audit now supplies the independent 2^32 reproduction.

## 2. Verification of the proofs

### (A1) (note lines 50-65)
- The bit-pair walk is correct:
  - A_k is valid iff a_k = 1 and b_(k+1) = 0.
  - The first exit from (1,1) and the last entry into (0,0) satisfy k <= k'.
  - If k = k', both A_k and B_k are valid.
  - Covers are distinguished by {last(P1), second(P2)}.
- Since q = (P1 minus its last vertex)(P2 minus its first two) is also recoverable, the proof gives the stronger count 2·h(T - {v,p,q}) per pair. This stronger form is READY (7) itself (see U1).

### (A2) (lines 67-86)
- The five cases by the start of P2 are exhaustive and each count is right:
  - P2 = z gives h(B') [paths ...y x] + R_a(B') [paths ending at y, x followed by a];
  - P2 = z y gives h(B') - e_a(B');
  - P2 = z y x ... gives 0;
  - P2 = z x gives h(B');
  - P2 = z x Q gives R_a(B').
- The equality analysis is right, including R_i(TT_m) = 2^(m-i).
- Confirmed on all 3,192 sink cases at order 8 and on every lower order.

### (A3)-(A4) (lines 88-92)
Correct. The equality case was checked in both directions: exactly one digit-4 old vertex per order.

### (A5) (lines 94-100)
- Correct, including "never 1" and "= 2 iff U is transitive with w its sink or its penultimate vertex".
- Every case is a special case of RESULT (12): sigma_a(C2(B,a)) = 2R_a(B) - h(B) - h(B - a) (see m10). I verified this on all 58,653 pairs (B, a) of orders 1-8.

### (B) transfer (lines 114-142)
- I re-derived all three rows. A path uses the arc u -> v at most once; the G-runs and T-runs then factor.
- In row 3, P2 either stays in G (one-arc join of T with G - S') or enters T through u -> v (plain join, factor R_v·K, with K summed over ordered triples (S, Y, Z)).
- The special cases (point gadget and C2) check by hand.

### (B) adjoint (lines 150-156)
- The double-join argument is valid. Reading (T1 ⋈ G) ⋈ T3 gives <M x1, y3>_J. The converse (T3^op ⋈ G^op) ⋈ T1^op gives <M' y3, x1>_J, and H(T) = H(T^op).
- (1,1,1), (1,0,2) and (3,1,3) are realisable, with determinant 2, so they span Q^3. Hence M'^T J = J M.

### (C) lemma (lines 168-173)
- The lemma is correct. Its proof uses only that M_i^(-1) is integral and that h is a nonnegative integer on C ∩ Z^3; positivity is not needed.
- The application has gaps; see M3.

### (D)
- **Join law.** I re-derived it: paths avoiding b -> a give h_A h_B - e_a s_b; paths using it split as T1 B1 a T2 B2, giving R_a(A)·R_b(B^op).
- **One-apex counts, derived independently.** A Hamiltonian path of T(w) is Q1 x Q2 with Q1 and Q2 increasing. This gives H(T(w)) = [w_1 = 1] + [w_n = 0] + #01(w) + sum_{b<a, w_b=1, w_a=0} 2^(a-b-1), i.e. g(v) = 1 + #10(v) + sum_{zero bits i} (v >> (i+1)). It also gives e_x = [v even] and R_x = #{Y : Y empty or w_min(Y) = 1} = 1 + v.
- **Suffix law** (lines 220-221). It follows from the closed form: for v = k 2^m + r,

      g(v) = g(k) + k C(r) - [k even][top bit of r = 0] + Q(r),

  with C(r) = sum over zero bits i < m of 2^(m-1-i) and Q(r) = [r even] + #01_m(r) + S_m(r). This is the session's form, with eps(r) = [top bit 0] and K(r) = S_m(r) + #10_m(r) = Q(r) - eps(r). It is a theorem, not just a checked identity (m4).

## 3. MAJOR issues

### M1. (B) is KNOWN: the general transfer and the general Lorentzian adjoint were proved on 29 October and are not cited

**Prior art.** `05-knowledge/results/seven_return_native_recursions_20261029/SEAM.md` (commit a421179977, 2026-10-29; INDEX.md line 229) proves:
- **Every block.** For arbitrary nonempty blocks with arbitrary entry and exit ports, coincident ports included (SEAM lines 6-7, 44-46), there is an integer 3x3 seam matrix with v(X⊙Y) = M_Y v(X) (eq. (1), construction (3)).
- **Composition.** M_(X⊙Y) = M_Y M_X (eq. (5)).
- **Adjoint.** M_(X*) = J M_X^T J with J = diag(1,-1,1) (eq. (6), line 119). The proof is the same: the two-block formula read from both sides, plus spanning vectors (1,0,2), (1,0,4), (3,1,3) (lines 105-119).
- **Point gadget.** The singleton matrix [1 -1 1; 1 -1 0; 1 0 1] (line 140) is S78's tower step.
- **Finite check.** `seam_transfer.py` line 90 checks the "universal mirror adjoint" on all labelled tournaments of order <= 4 with every port pair.

SEAM's v = (H, S, R_end) on X is S78's (h, e, R) on X^op, and X⊙Y is the converse of S78's T => G. So eq. (6) is exactly S78's M(G^op, z, u) = J M(G,u,z)^T J.

The general block transfer is also proved in `golden_hamilton_transfer_20261103/RECURRENCE.md` ((G1)-(G4), lines 35-110, 2 November). Its (G1) uses S78's own join convention.

**S78's claims.** S78 credits only session29's 18-tile case and claims the general identity as new:
- note lines 15-16: "the explicit coefficients and the general adjoint identity are added here";
- note lines 158-159;
- typing table lines 240-241;
- THM-4875 status lines 6-8 and section (B).

**What is new in (B).** Only:
- the closed subset-sum expressions for the nine entries (E_u, E_(u,z), J_1, J_2, K);
- the larger finite checks (order <= 5 against direct counts; order 6 for the adjoint).

**Retype:**
- transfer law KNOWN (SEAM (1)-(3); golden_hamilton (G1)-(G4));
- adjoint KNOWN (SEAM (6));
- closed-form entries PROVED (new presentation);
- FINITE-EXACT as stated.

The S73 audit, which cited only session29 and polarized_path_calculus, is probably where the omission started.

### M2. (A) claims as new a bound already in canon; THM-4873 (ii) is misdescribed

**THM-4873 (ii) already contains it** (canon, opus-S77, same day; lines 55-62; S77 `c2_digit_heredity_20261104_audit_S77/AUDIT.md` section 4, lines 142-177):
- (D) = (A1), with the same walk proof.
- "a small digit can sit only at x, at y, at the sink z of B, or at a with d_B(a) <= 1". So every other old vertex has digit >= 6 (S78's (A3)).
- "**The sink.** sigma_z >= h(B - z) + B_a(B - z) + 2 >= 4" (S78's bound in (A2) and (A4)).
- The whole of (A5), including "never 1" and the digit-2 case.

S78 lists THM-4873 as INHERITED, but then misdescribes it:
- Note lines 108-110: "THM-4873 (ii) ... classifies digits 1 and 2. Here the bound is completed to the sharp old-vertex minimum 4". THM-4873 (ii) already proves the minimum 4.
- THM-4875 status lines 5-6: "(A2) and (A4) are new (exact sink digit; sharp minimum 4 and its equality case)". The minimum-4 bound is not new.
- Note lines 13-14 and the title present "the sharp minimum 4" as the session's result.

**What is new in (A):**
- the exact sink formula 3h(B') + 2R_a(B') - e_a(B') (it strengthens THM-4873's h(B') + R_a(B') + 2);
- the equality case of digit 4 (the "only if" direction already follows from THM-4873's bound; the "if" needs the exact value);
- with this audit, the >= 9 theorem (U1).

**Further credits:**
- (A5) is a corollary of RESULT (12) (`multiport_defect_grammar_20261103/construction/RESULT.md` line 177: S_a' = 4S - H(B-a), i.e. sigma_a = 2R_a(B) - h(B) - h(B - a)).
- ONE_EXIT_AND_PERSISTENCE (lines 118-121, 149-152) already had ">= 6" for every old non-mark vertex when B is strong.

**Retype:** (A3), (A4)'s bound and (A5) KNOWN (THM-4873 (ii)); (A2)'s exact formula and both equality cases PROVED (new).

### M3. (C): the lemma is fine; its stated consequence for Chan-Pak is not proved, and as worded it targets a false statement

The consequence appears in note lines 175-186 and in THM-4875 (C), lines 83-85: "So a finite-base unimodular cone-covering induction for Chan-Pak is impossible".

**(i) The finite base is not handled.**
- The lemma yields one point of C outside the union of the M_i(C). A finite base F can absorb finitely many such points.
- Excluding every finite F needs a further argument. The uncovered set is relatively open in C (the union is closed) and is a cone. When C has nonempty interior, it therefore contains an open subcone, which holds infinitely many lattice points.
- That argument is absent, and the statement is false for lower-dimensional C (for a ray, one base point suffices). See U2 for a corrected statement and proof.

**(ii) The target is false, so the no-go is vacuous as worded.** Line 182 reads: "all primitive triples of such a cone are realisable".
- Every realisable triple has h odd (Rédei). (2, 1, 3) is primitive, lies in {0 <= e <= h, R >= (1+c)h} for c <= 1/2, and is not realisable.
- Even on the coset {h odd}, (3, 3, 4) is primitive and lies in the cone (c <= 1/3) but is not realisable. By THM-4876 (7), or the last-strong-component argument in THM-4873 (i)'s proof, e = h forces the root to be the sink, and then R = h.
- So no cone-intersect-lattice target of this kind is realisable. Ruling out a proof method for it says nothing about Chan-Pak, which is a statement about h alone.
- A meaningful version would have to describe the realisable set itself. That set is not a cone intersected with a coset, and the descent argument gives no obstruction for it: an induction "every realisable triple in C is M_i of a smaller realisable triple, or in F" is consistent with descent.

**(iii) Unstated hypotheses.**
- The inequality h' >= R uses 0 <= e <= h (for h_G h >= s_u e), so the cone must lie inside {0 <= e <= h}.
- {0 <= e <= h, R >= (1+c)h} contains (0, 0, 1), where h = 0. So the lemma's hypothesis "positive integral functional on C - 0" (line 169) fails as written. It is harmless, because nonnegativity suffices, but it should be stated.

**(iv) Few gadgets are unimodular.**
- The lemma needs det M_i = ±1. The application paragraph ("every gadget transfer ... the lemma applies", lines 175-180) never says which gadgets qualify.
- By exact census (`aud_gadget.out`, `aud_descent.out`):

  | gadget order | det ±1 | singular | abs(det) >= 2 | total |
  |---|---|---|---|---|
  | 4 | 96 | 640 | 288 | 1,024 |
  | 5 | 1,200 | 12,240 | 12,160 | 25,600 |
  | 6 | 20,160 | 384,768 | 774,720 | 1,179,648 |

- So the no-go covers only inductions built from a minority of gadgets: 9.4% at order 4, 4.7% at order 5, 1.7% at order 6.

**(v) The equality clause is misdescribed.** Lines 177-179 say "with equality only for neutral insertions (R = h, i.e. digit 0; for example a sole exit into the sink)".
- R = h is necessary, not sufficient.
- **Proof of the exact condition.** h' - h = (h_G h - s_u e) + (E_u - 1)R + (R - h), with three nonnegative terms. The first vanishes iff h_G = s_u and e = h, the second iff E_u = 1, and the third iff R = h. With u the source of G, E_u = h(G - u). So h' = h iff:
  - e = h = R, i.e. the base port is a sink (a non-sink has e <= (h-1)/2 by THM-4876 (7));
  - G is transitive and u is its source.
- Checked: all 13,661 equality cases among the 72,958,060 (gadget of order <= 5, realisable triple) pairs match exactly.

**(vi) Uncited precursors of the "moves that keep h fixed" conclusion** (lines 183-186):
- THM-4865 (f) (line 55): "The missing ingredient is a construction that moves B in small controlled steps at fixed h".
- HYP-9453 section 2 (line 30): "what is needed is a construction that moves B by one while fixing h".
- klein's reflection `07-reflections/the-root-response-one-object-behind-...-klein-20261104.md` line 18: "Local moves cannot do it".

**Retype:** (C) lemma PROVED (elementary). The application is a heuristic remark: delete the THM sentence, or replace it with U2.

## 4. MINOR issues

1. **Class and mark counts mixed** (note line 102; THM line 9). "All 7,412 classes of order 1-8, all 58,652 marks" mixes ranges: 7,412 includes order 1, while 58,652 excludes its single pair. The output's pair counts sum to 58,653 over orders 1-8. Write "7,411 classes of order 2-8 and their 58,652 marks".

2. **Probe bases** (lines 161-163, THM line 10). The "five probe bases" realise only four distinct triples: (1,0,2) occurs twice, as TT_2 at its source and TT_3 at its middle vertex (`gadget_adjoint_20261104.out` line 1). Linearity beyond the three spanning probes was therefore tested at a single extra point per gadget. The auditor's 456,745 checks now cover this.

3. **Neutral insertions** (lines 177-179). See M3 (v).

4. **Suffix law typing** (lines 220-221). The suffix law is a theorem: Session191 `one_apex/PROOF.md` section 2 ("The dependence of the eventual count on an initial state (h,A,p) is h + A C(v) + p epsilon(v), plus a term independent of that initial state"), or section 2 of this audit. Type it PROVED (inherited) rather than "checked on 2,000,000 random v".

5. **One-apex credits** (lines 39-44; THM lines 13-14).
   - The apex triple is proved in `session179_niche_methods_20261104/spectrum/REPORT.md` (1), lines 23-27: "exactly 1-p old path ending at a", "1 + sum_{i:b_i=1} 2^(m-i) = 1+s". It is not merely "checked by brute force".
   - g is session124's digital function f (`session124_digital_constructor_inverse_20261031/PROOF.md` lines 59-61; overnight 20261031d README line 95).
   - The direct count is Session191 `one_apex/PROOF.md` eq. (2).

6. **The new conjecture and the two-apex language** (lines 207-208: "almost every odd number is the count of a join of two one-apex tournaments").
   - F1 lies inside the two-apex (two-feedback-vertex) word language. Deleting x1 and x2 leaves chain1 => chain2^op, which is transitive.
   - That language already has an OPEN completeness conjecture: `square_sum_spectrum_overnight_20261030/README.md` lines 41-47 ("Every positive odd number other than 7 and 21 is H(w) for this single four-letter, two-feedback language").
   - It also has finite coverage: to 99,999 there; to 230,519 in Session191 `two_apex/TWO_APEX.md` lines 165-172.
   - Cite both. The new conjecture is a density statement for a subfamily.

7. **THM-4876 is not cited anywhere in S78's files.**
   - Statement (6), the "one-centre family" TT_m + u, is the one-apex family.
   - Statement (7), the end-count bound e_a <= (h-1)/2, was sent to S78 as MSG-2448 and corrected in MSG-2449. It sharpens the e_a <= h steps in (A2).
   - The audit brief lists THM-4876 among S78's credits, but no S78 file contains it.

8. **(C) precursors.** See M3 (vi): THM-4865 (f), HYP-9453 section 2, klein's reflection.

9. **Home of (A1)** (line 50, table line 236). The canonical home of (A1) is THM-4873 (ii) (D). READY (7) (`two_entry_funnels_20261103/READY.md` lines 18-19) is stronger than d(d-1): it states delta = F - h - q + s = 2 sum b_ab with b_ab >= H(T - {r,a,b}). The full form is what proves U1.

10. **(A5) and RESULT (12).** (A5) is a special case of RESULT (12) (M2). It holds on all 58,653 pairs (B, a) of orders 1-8 (`aud_oldmark_8.out`, column "RESULT(12)").

11. **P2 index.** HYP-9468 P2 pre-registers m_(k+1)/m_k for k = 22..31, which needs m_32. The results report m_k/m_(k-1) for k = 22..31. All values lie in the window either way.

12. **Runtime** (line 201). "About 30 minutes": the file times give 21.7 minutes (17:46:22 to 18:08:01).

13. **Greedy output** (lines 227-228). The 111-triple greedy result had no committed output. The auditor reproduced it with the session's `greedy_bank.c` and bank order (`greedy_repro_20.out`: 10,892 targets, 111 picks).

14. **Row-3 wording** (lines 137-140). The row-3 proof reuses S for two different sets ("S contains T", "T[S] => G[S2]"), while the definitions of J_1 and J_2 use S as a subset of G. Name the G-part of P1 explicitly.

15. **F0 is redundant** (THM (D) line 92; note lines 191-192). See U3. The third bullet can be dropped, which makes the THM title ("one-apex one-arc joins realise every odd number") literally true.

16. **Proved bound** (line 106, table line 239, THM lines 10 and 58-59). "Proved bound 6" is superseded by U1: PROVED >= 9, sharp.

## 5. Upgrades proved by the auditor

### U1. Non-sink old vertices have digit at least 9, with equality only for (TT_n, a sink, w penultimate)

This settles HYP-9468 P3's all-order conjecture.

**Theorem.** Let T1 = C2(B, a), |B| >= 2, and let w != a be an old vertex with d = d_B(w) >= 1. Then sigma_w(T1) >= 9. Equality holds iff B is transitive, a is its sink and w is its penultimate vertex.

**Lemma (READY (7), dualised in full).** For every vertex v,

    sigma_v >= E_v + 2 * sum over pairs {p,q} in N+(v) of h(T - {v,p,q}),

where E_v = h(T - v) - e_v(T) counts the Hamiltonian paths of T - v ending in N+(v).

*Proof of the lemma.*
- The covers with P2 = (v) are exactly the paths counted by E_v.
- For |P2| >= 2, take each pair and each Hamiltonian path Q of T - {v,p,q}. The (A1) walk gives at least 2 covers of the form (Q_<=k p', v p'' Q_>k).
- (pair, Q) is recoverable from the cover.
- This is READY (7) (delta = F - h - q + s, b_ab >= H(T - {r,a,b})) applied to T^op.

*Proof of the theorem.*
- If d >= 2, (A1) gives (d+2)(d+1) >= 12.
- Let d = 1, with w -> c the only out-arc of w in B, and U = B - w. Then N+_T1(w) = {c, x, y} and T1 - w = C2(U, a).
  - E_w = e_c(C2(U,a)) + h(U) + R_a(U). The last two terms are the paths ending at x and at y, as in (A2).
  - The pairs give h(U), h(U - c) (y is a sink there) and h((U - c) + x) >= 1.
  - Hence

        sigma_w >= e_c(C2(U,a)) + R_a(U) + 3h(U) + 4.        (*)

**Case 1: c = a.**
- The paths [U - a] y x a end at a, so e_a(C2(U,a)) >= h(U - a) >= 1, and (*) >= 1 + 1 + 3 + 4 = 9.
- Equality forces h(U) = R_a(U) = 1, i.e. U is transitive with sink a. Since w -> a only and U - a -> w, B is then transitive with a its sink and w penultimate.
- Conversely, for that B, counting by the start of P2 gives exactly 3 + 1 + 1 + 2 + 2 = 9.

**Case 2: c != a and h(U) >= 3.** Then R_a(U) >= h(U), so (*) >= 3 + 9 + 4 = 16.

**Case 3: c != a and U = u_1 -> ... -> u_m transitive, with c = u_j and a = u_k.**
- **(3a) k < j.** The path [U - {a,c}] y x a c gives e_c >= 1. Also R_a(U) = 2^(m-k) >= 2. So (*) >= 10.
- **(3b) k > j.** Then c is not the sink of U. The cover

      P1 = u_1 ... u_(j-1) y x,   P2 = w u_j u_(j+1) ... u_m

  is valid and is not counted in (*). For the pair {c, x}, the only Hamiltonian path of T1 - {w,c,x} is u_1 ... û_j ... u_m y. The prefix u_1 ... u_(j-1) y would need j = m.
  - If a is not the sink of U, then R_a >= 2 and sigma_w >= (*) + 1 >= 10.
  - If a = u_m is the sink, counting by the start of P2 gives the exact value

        sigma_w = 2 + 0 + [j = m-1] + (2 + [j = m-1]) + 3·2^(m-j) >= 12.

    The five terms are: P2 = (w); w y; w y x...; w x...; w c... . In the last term, W(c) = 2^(m-j) + 2^(m-j-1) + 2^(m-j-1) + 2^(m-j): P1 ends at x with y in P1 or in P2, or P1 ends at y with x in P1 or in P2.
    The formula was checked by brute force for every 2 <= m <= 6 and 1 <= j < m (e.g. B = C_3 gives 12).

So sigma_w >= 9 in every case, with equality only in Case 1's configuration. ∎

**FINITE-EXACT.** At every order 2-8 (`aud_oldmark_8.out`, independent class generation):
- the minimum non-sink old digit is 9;
- exactly one (class, mark, vertex) attains it, and it is (TT_n, a sink, w penultimate);
- the lemma's inequality holds at every vertex of every child.

**Remark.** Digit 4 sits at (a penultimate, w sink) and digit 9 at (a sink, w penultimate). On B = TT_n the two minimal configurations are exchanged by swapping the roles of a and w.

### U2. A correct form of (C)

**Statement.** Let C be a closed convex cone with nonempty interior inside {0 <= e <= h, R >= (1+c)h} (c > 0). Let M_1, ..., M_r be gadget matrices with det M_i = ±1. Let Λ ⊂ Z^3 be any coset of a full-rank lattice. Then for every finite F,

    C ∩ Λ is not contained in F ∪ M_1(C ∩ Λ) ∪ ... ∪ M_r(C ∩ Λ).

**Proof.**
- On C - 0, h(M_i q) >= E_u q_R >= q_R > q_h if q_h > 0. If q_h = 0, then q_e = 0 and q_R > 0. In both cases h strictly increases.
- h is a nonnegative integer on C ∩ Z^3. By the session's descent lemma, some p0 in C lies outside the closed set ∪ M_i(C).
- A ball around p0 misses ∪ M_i(C) and meets int C. The open cone it generates misses every M_i(C) and contains infinitely many points of Λ. ∎

This rules out only inductions that use unimodular gadgets exclusively and aim at a whole cone-intersect-coset target. M3 (ii) explains why such targets are never fully realisable.

### U3. F0 is redundant: F1 + F2 alone reach every odd N <= 2^32 except 7 and 21

- **g-values are F1 values.** F1(0, t) = 1 + t + g(t) - e(t) = g(2t) (Session 193 `one_arc_joins/PROOF.md` section 3, post-checkpoint), and g(2t+1) = g(t). The F1 loop includes u = 0 and every t with F1(0,t) <= X.
- **Bank h-values.** The bank h-values are at most 189. The auditor's `-DNOF0` run shows that F1 + F2 miss exactly {7, 21} on [1, 2^28] (`aud_cover_noF0_28.out`).
- **Conclusion.** Combined with the full 2^32 run, every odd N <= 2^32 other than 7 and 21 is the count of a one-arc join of a one-apex tournament with either a one-apex tournament or a rooted tournament of order <= 7.

### U4. Exact digit in the residual family

For B = TT_m + w with w -> u_j only (j < m) and a = u_m:

    sigma_w(C2(B, a)) = 4 + 2[j = m-1] + 3·2^(m-j).

## 6. Files (this folder)

| file | purpose |
|---|---|
| `aud_bank.c`, `aud_bank7.txt` | independent rooted-triple bank (2,732 triples) |
| `aud_cover.c`; `aud_cover_28.out`, `aud_cover_32.out`, `aud_cover_noF0_28.out` | independent coverage (build with `-DNOF0` for U3) |
| `aud_witness.c`; `aud_witness.out`, `aud_witness_large.out` | one-apex triples to length 12; 41 explicit witnesses by direct DP |
| `aud_oldmark.c`; `aud_oldmark_8.out` | (A1)-(A5), sink formula, equality cases, U1, RESULT (12), READY (7) dual, orders 1-8 |
| `aud_gadget.c`; `aud_gadget.out` | (B): probe-fitted M vs closed formula, linearity, adjoint, determinants |
| `aud_descent.c`; `aud_descent.out` | (C): h' >= R >= h, equality characterisation, determinant histogram |
| `aud_f1res.c`; `aud_f1res_24.out` | residue profile of F1 misses with the auditor's g |
| `greedy_repro_20.out`, `omb_rerun_8.out` | reruns of the session's greedy_bank.c and old_mark_digit_bound_20261104.c |

**Reproduction** (from this folder, gcc 15 with OpenMP):

    gcc -O3 -march=native -o aud_bank aud_bank.c && ./aud_bank 7 > aud_bank7.txt
    gcc -O3 -march=native -fopenmp -o aud_cover aud_cover.c && ./aud_cover 32 aud_bank7.txt
    gcc -O3 -march=native -fopenmp -DNOF0 -o aud_cover_noF0 aud_cover.c && ./aud_cover_noF0 28 aud_bank7.txt
    gcc -O3 -march=native -fopenmp -o aud_witness aud_witness.c && ./aud_witness 12 40 7 22 && ./aud_witness 0 8 11 26 23 25
    gcc -O3 -march=native -fopenmp -o aud_oldmark aud_oldmark.c && ./aud_oldmark 8
    gcc -O3 -march=native -fopenmp -o aud_gadget aud_gadget.c && ./aud_gadget 5
    gcc -O3 -march=native -o aud_descent aud_descent.c && ./aud_descent 5 aud_bank7.txt
    gcc -O3 -march=native -o aud_f1res aud_f1res.c && ./aud_f1res 24
