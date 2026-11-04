# Independent audit: HYP-9464 and the note `fermat_self_duality_fewcycle_coincidences_20261104.md`

**Auditor:** an independent general-purpose subagent of kind-pasteur-2026-11-04-S10, briefed only with the files to
audit and instructions to re-derive and break every claim with its own code. It ran on 2026-11-04 against the
pre-audit checkpoint 220666527f (pushed to main as 9aa6e0836f) and the working copy. The harness did not let the
subagent write a `.md` file, so the session transcribed the report from the auditor's final message. The six programs
and their outputs in this folder are the auditor's own.

**Verdict: PASS WITH CORRECTIONS.** 1 MAJOR, 11 MINOR. Every PROVED claim re-derived correctly, and every FINITE-EXACT
number was reproduced by independent code. The checks extend two ranges: order 11 for self-duality, and orders 11-14
for the deletion law. The working copy changed during the audit; that edit already fixed F6 and the rigidity part of
F11.

All corrections were applied in the commit that adds this file (MISTAKE-712).

## Findings

**F1 (MAJOR): the mirror parity reading is canon and was presented as new.**
- The involution iota(P) = reverse(sigma(P)), H ≡ Fix (mod 2) and the oddness of Fix were already canon:
  - THM-582 (mac-mini, 2026-06-29);
  - THM-647 (Anti-Redei, monad, 2026-07-07);
  - THM-024 (every self-converse tournament has an involutory anti-automorphism);
  - THM-583 (the half-system principle behind Part 3).
- Only the explicit AT_n count, Fix(AT_(2k+2)) = H(AT_(k+2)), is new.
- Fix: the note's section 0 rows 4-5, sections 2, 3 and 9, and the HYP Results now cite these theorems.

**F2 (MINOR): the Fermat condition was backwards.** Section 4 read "a power of 2 being one more than a prime". It now
reads "a prime being one more than a power of 2".

**F3 (MINOR): x = 1 also solves x^(x^2) = (x^2)^x.** Fix: x > 0, x != 1.

**F4 (MINOR): the module dichotomy.**
- The scratch-only script reference is replaced by the committed `selfdual3_modules_20261104.py`.
- "TT_2 to TT_8" was false at order 10.
- The claim is now typed FINITE-EXACT at orders 7 and 8 and CITED (Theorem 4) from order 9 on.

**F5 (MINOR): order 8 has a short proof** of the exclusion of C3 modules (two outside vertices on the same side), and
n >= 9 follows from Theorem 4. Both are added.

**F6 (MINOR, fixed during the audit): "is every (AT_m, {m-1, m}) a block?" was listed as OPEN.** The working copy
already contained a proof. The auditor confirmed it for m = 4..12, and S75's block list by brute force at orders 3-6.

**F7 (MINOR): four credits.**
- (a) The rooted count R(m) = F(2m - 5) is THM-4868 (F); Garey 1972 counts the classes, F(2m - 6).
- (b) Hereditary self-converse is KNOWN: Reid-Thomassen 1976, as quoted in arXiv:1204.2513 section 1.2, with the
  converse for strongly self-dual tournaments of order >= 8.
- (c) The witness facts are credited to S75 section 1 and the S74 audit.
- (d) The Ogg level filtration is credited to the repo note, with the lists from Aricheta 2019.

**F8 (MINOR): THM-4748 was overstated.**
- Nothing proves that an apex is necessary on the Mersenne side.
- Two arcs can subtract (15, 23); what fails is reaching 2^k - 1 for k >= 5.

**F9 (MINOR): the nsum caution.**
- 2.09688... is 3 - length(E), not length(E).
- The error depends on precision: 2.09732 (dps 15), 2.09693 (dps 30), 2.09688 (dps 40).

**F10 (MINOR): records cited but not committed.**
- The networkx check is now `selfdual3_networkx_check_20261104.py`.
- The raw selfdual3 and hcdump outputs are now in `selfdual3_fewcycle_20261104_records/`.
- Section 10 says how to regenerate them.

**F11 (MINOR): qualifiers.**
- Span d >= 2 in the single-arc formula.
- n >= 4 for Fact 4.
- "Rotation-rigid" for the circulant sense of "rigid".
- "Every almost transitive tournament" in section 7, in place of "every Fermat tournament".

**F12 (MINOR, upgrade): ranges extended by the audit.**
- (A): none indecomposable at order 11 (complete pruned search, validated at orders 7 and 8).
- (B): the law holds at orders 11-14 (hc = 3: 233, 610, 1,597, 4,181 classes; hc = 4: 178, 466, 1,220, 3,194).

## Independently reproduced

- **`aud_census.c`** (orbit-marking census of every labelled tournament, orders 3-8):
  - class counts 2, 4, 12, 56, 456, 6,880;
  - self-converse classes 2, 2, 8, 12, 88, 176;
  - {-3}-self-dual classes 56, 10, 18, with multisets identical to the explore record; order 8 has 16 decomposable
    classes plus the indecomposable H = 633 pair;
  - all modules at orders 7 and 8 are transitive;
  - the pruned search finds only TT and AT at orders 9-11;
  - THM-4749's n = 8 row, including N(8, 257) = 30;
  - maximum H 3, 5, 15, 45, 189, 661.
- **`aud_fewcycle.c`** (every strong class with hc <= 5, orders 3-14):
  - t_1 = F(2n - 6), and u_3, u_4, u_5 equal S75's values;
  - the deletion law at orders 6-14, with the hc = 3 deletion profiles matching the session's probe;
  - uniquely Hamiltonian classes with two sole exits number F(2n - 7);
  - S76's split, F(2n - 10) classes with |Aut| = 3 and F(2n - 11) rigid.
- **`aud_at_family.py`** (ALL CHECKS PASSED):
  - substitution form, and heredity for n <= 12;
  - H(AT_n) and the mirror law to n = 16, by two path enumerators plus a DP;
  - the parity check on 240 random cases;
  - end counts, the block list, and the generating-function identity to n = 44;
  - |Aut(AT_n)|;
  - the order-7 and order-8 C3 lemmas;
  - the Ogg and Platonic arithmetic.
- **`aud_witnesses.py`**: the four order-9 witnesses have hc = 9, |Aut| = 9, two C3 modules and the forced arc x -> y.
  They match S75, and none is {-3}-self-dual.
- **`aud_mrb.py`**:
  - B = 2.0967911243915261257250306750849, by Euler-Boole and Cohen-Villegas-Zagier; all 28 digits printed in the
    PDF are correct;
  - B = 3 - length(E) to 2e-51;
  - the MRB constant to 50 digits, and the tree-function identity.
- **`aud_block_construction.py`**: the order-(m + 2) construction on the AT_m block, for m = 4..12.
- **Code review of the `selfdual3` mode:** no logic errors. The pre-registration timing is consistent: commit at
  15:12:18, order-10 runs at 15:13 and 15:17.

**Observation (not a finding).** For hc = 5 the deletion split follows S75's decomposition term by term at orders 7-14.
L(2n - 12) classes have exactly 2 deletions and F(2n - 13) exactly 3. This is a candidate extension of HYP-9464 (B).
