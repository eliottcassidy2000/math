# Independent audit of the packet `mod4_race_square_term_20261111` (THM-5131, THM-5132, HYP-9663)

Auditor: independent adversarial subagent, 2026-11-11. Everything below was checked with my own reasoning and my own
programs (this folder); the author's programs were read only to learn file formats. Nothing outside this folder was
edited, nothing was committed.

**Audited version.** Commit `aa0edea9c1` ("checkpoint (pre-audit)") **plus the uncommitted changes that appeared in the
working tree at 09:22-09:25, after that commit and during the audit**: `README.md` (+40/-8 lines: the remark "Expected
kernel size", extra proof text in 6.2 (3) and (5), the paragraph "The period-8 phase (IDENTITY) and the Gaussian
integral (ANALOGY)", two lines in section 10) and the untracked files `check_selfconverse_mod8.py`,
`kernel_size_check.py`, `out/check_selfconverse_mod8.out`, `out/kernel_size_check.out`. Section names are used instead
of line numbers because the note was moving.

## Verdict: PASS WITH CORRECTIONS

- **No mathematical error was found.** Every proof (Theorems T, S, Proposition 5.1, the deck remarks, Theorem P
  (1)-(5), the conditional form, the corollaries, the kernel-size remark) was re-derived by hand, and every statement
  and every table entry I could reach was recomputed with independent code, in most cases beyond the range asked for
  (all tables to order 10, `D_11` and the residues mod 8 at order 11 exactly, an independent order-12 Monte Carlo).
- **Two MAJOR findings, both of the kind "known result typed as new"** (the third session in a row for this author):
  the mod-8 refinement Theorem S (2) is in the repository twice, with the same proof, and was pointed out to this
  author three days ago; the closed form of Theorem T (3) has been in canon since June (THM-555, HYP-2707).
- **22 minor findings**: two wrong last digits in a row labelled exact, a `10^6` typo, wrong or loose sentences in the
  outcome list and in section 8, citation errors, one overstatement of the order-12 result, a wrong start time in the
  pre-registration file.

---

## A. Findings

### MAJOR

**M1. Theorem S (2), "H = 1 + 4s (mod 8)", is not new: identical statement and proof are in the repository twice.**

- *Where typed new.* README section 4, Theorem S (2) and its proof (no credit); THM-5131 title ("more precisely
  H = 1 + 4s mod 8 ..."), status line "NEW here: (T); the mod-8 refinement (S) (2); the count (S) (3) ...", body
  "(S) 2. (PROVED.)".
- *Prior statements.*
  - `05-knowledge/results/session29_framed_profiles_20261022/session29-converse-congruences/REPORT.md` (commit
    `25a4f6ba21`, 2026-10-22), section 3 "The reflection-overlap congruence": "For every self-converse tournament,
    H(T) ≡ 1+2F+4I (mod 8). Indeed OCF gives H≡1+2 alpha_1+4 alpha_2 modulo eight. Here alpha_1=F+2A. ... so
    alpha_2≡D modulo two. ... At even order this specializes to H(T) ≡ 1 (mod 4), H(T) ≡ 1+4I (mod 8)", with `I` the
    number of mirror pairs `{C, rho(C)}` whose supports meet. This is Theorem S (2) word for word (`s = I`), by the
    same argument. Its section 1 also proves what Theorem S (1) uses and THM-791 does not contain (an anti-automorphism
    gives an involutive one fixing at most one vertex) and, beyond the note, that all anti-involutions are conjugate,
    so that `I` is an invariant of the class.
  - `05-knowledge/results/session130_blue_black_theorems_20261031/hamiltonian/PROOFS.md` (commit `41a9d96041`,
    2026-10-31), Theorem 1: "For every finite tournament with an involutory anti-automorphism,
    `H(T) = 1 + 2q + 4i (mod 8)`. ... For even order, `q=0`, hence `H(T)=1+4i (mod 8)`. This refines, but does not
    replace, inherited THM-791."
  - The pointer was given to this author: `05-knowledge/results/eternal_colours_audit_K1_20261108.md` (audit of this
    machine's THM-4930), lines 257-259: "Session29's converse-congruences report (§§1, 3) proves the same parity, with
    a mod-8 refinement `H ≡ 1 + 4I (mod 8)` at even order. The THM-4922 audit observed it on all 8,784 self-converse
    classes of order 10." THM-4930 (same machine, 2026-11-08) cites that report four times.
- *Related, same theorem.* (S) (1): the reduction to an involution is Session29 section 1 / THM-4930, not THM-791.
  (S) (3): the count per involution is `session29-halfpath/REPORT.md` section 5 ("Fix sigma and choose uniformly among
  all labelled tournaments admitting it. There are 2^(m^2) choices at n=2m and 2^(m(m+1)) at n=2m+1"; section 5b: "A
  self-converse class T contributes exactly 2^m m!/|C_G(sigma)| distinct labelled tournaments in the fixed-sigma
  ensemble"), also THM-4743 and THM-4930 ("enumerating the tournaments anti-invariant under one fixed involution covers
  all of them"); the bound of (3) is the union bound over the `(n-1)!!` (resp. `n!!`) involutions, one line. The numbers
  `704, 7040, 380416` (= `SC_5, SC_6, SC_7`) and `960, 7680, 430080` (the bounds) are printed in
  `05-knowledge/results/selfk_type_ledger_opus_S312.out` and `selfk_type_ledger_fixed_opus_S312.out`.
- *My check of the mathematics.* True: `audit_misc.py S` finds an anti-automorphism by backtracking, builds the
  involution, lists all directed odd cycles and verifies `H = 1 + 4s (mod 8)` on all 2 + 12 + 176 self-converse classes
  of orders 4, 6, 8 (and `H = 1 + 2F + 4I` on the 2 + 8 + 88 of orders 3, 5, 7): 0 failures.
- *Correction.* Retype (S) (2) as KNOWN (Session29 converse-congruences section 3; Session130 Theorem 1), credit (S) (1)
  to Session29 section 1 / THM-4930 as well as THM-791, state in (S) (3) that the count per involution is Session29
  half-path section 5 and that the union bound is its corollary; remove both from "NEW here" in THM-5131 (title, status,
  body) and from the outcome list; keep as new only the exact labelled values `SC_n` (n = 3, 4, 8, 9, 10; orders 5-7
  were on file as by-products), the split of `D_n`, the shares, the class split `82/94`, `4448/4336`, and the reading as
  a square term. Record a MISTAKE entry (same lesson as MISTAKE-738, MISTAKE-753, MISTAKE-755: the report was cited
  by the author's own theorem of 2026-11-08).

**M2. Theorem T (3): the closed form `E (-1)^(c_3) = 2^(-floor((n-1)/2))` has been in canon since June; "Not found in
the repository" is false.**

- *Where typed new.* README section 3, Typing: "Elementary. Not found in the repository (THM-478 has the quadratic form,
  THM-552 the self-converse case)"; outcome item 1; THM-5131 status "NEW here: (T)".
- *Prior statements.*
  - `01-canon/theorems/THM-555-tiling-cycle-moment-hierarchy.md` (kind-pasteur-2026-06-20-S20, commit `f727806f8b`),
    line 49: "E[(-1)^c3] = 1 / 2^{floor((n-1)/2)}   (0 at n=3) => Pr[c3 odd] = 1/2 - 1/2^{floor((n-1)/2)+1}  [VERIFIED
    n<=10]", and under "Repo connections": "THM-552 (c3-parity dichotomy): the ensemble even-bias
    `E[(-1)^c3]=1/2^{floor((n-1)/2)}>0` is the statistical SHADOW of the pointwise dichotomy." THM-552 is the file the
    note cites.
  - `05-knowledge/hypotheses/HYP-2707-tournament-clifford-magic-hierarchy.md` (kind-pasteur-2026-06-20-S22), "PROVED
    CORE (kps-S22, exact n=4..8)": "The c3-parity ... is a GF(2) quadratic ... function ... Its associated bilinear
    form ... has GF(2) rank `r = 2*floor((n-1)/2)` ... The standard Boolean Gauss sum then gives `E[(-1)^c3] = ... =
    2^{-floor((n-1)/2)}`". This is the mechanism of the note's proof (a Gauss sum of the triangle form).
  - Both are stated for the tiling ensemble (tournaments containing a fixed Hamiltonian path) and are verified, not
    proved, for general `n`. The two ensembles give the same value for every `n >= 4` for a one-line reason (the
    complement of a Hamiltonian path in `K_n` is connected for `n >= 4`, so the score-parity vector is still uniform on
    its coset; at `n = 3` it is not, whence the 0). Checked exhaustively for `n = 3..8` (`audit_tiling.py`).
- *What is new.* The proof for all `n`; the uniform ensemble (including `n = 3`); the pointwise statement (T) (2)
  (every base has cube sum `2^ceil(m/2)`), which I found nowhere. Theorem T therefore also upgrades THM-555's line from
  VERIFIED `n <= 10` to PROVED for all `n >= 4`; say so.
- *Correction.* Replace "Not found in the repository" by the credit to THM-555 and HYP-2707; in THM-5131 type (T) (3)
  as "closed form first recorded in THM-555 (verified, tiling ensemble; Gauss-sum mechanism HYP-2707); proved here for
  all `n` in the uniform ensemble, with the pointwise form (T) (2)". See also m1 (the result is a two-line corollary of
  a classical formula).

### minor

**m1. Theorem T and the first formula of Theorem P (1) are two-line corollaries of the classical formula
`c_3 = C(n,3) - sum_v C(s_v, 2)`** (Kendall-Babington Smith; in the repository: THM-118, THM-466 "(v)" and its table
row "c3 mod 2 (classical)", THM-559). Since `sum_v s_v^2 = q_1 (mod 4)` and `sum_v s_v = C(n,2)`,
`c_3 = C(n,3) + (q_1 - C(n,2))/2 (mod 2)`, which equals `C(q_1, 2) + floor(n/4)` for every admissible `q_1` (`audit_kbs.py`,
`n <= 200`). (T) (2) then follows because `x -> (s + 1 + x, |x|)` is a bijection from the out-sets onto the admissible
score-parity vectors of `R + v_X`, so the cube sum is `(-1)^floor(n/4) sum_{k = C(n,2) (2)} C(n,k) (-1)^C(k,2) =
2^floor(n/2)`, `n = m + 1`. The arc-reversal proof of the note is correct; the typing should say "elementary corollary
of the classical score formula" and (T) should cite it. No outside source for the parity statement was found either
(`literature_notes.md`).

**m2. Row `n = 11` of the table of residues mod 8 (README section 5; outcome item 5; THM-5131 FINITE-EXACT list) is not
exact in the last digit.** The note has `1.00052, 0.99948, 1.00044, 0.99956`. Exact (own census of all 9,733,056
bases of order 10 with the vertex-cube formula, `audit_cube.c`): labelled counts
`9011754468980736, 9002583394936832, 9011181642270720, 9003277512775680` of `2^55`, i.e. fractions
`0.2501264, 0.2498719, 0.2501105, 0.2498911` and relative masses **`1.00051, 0.99949, 1.00044, 0.99956`**. The note's
digits are the five-digit fractions of the second audit of THM-5127 multiplied by 4. The sentence "rows 10 and 11
are exact" should be kept only after the correction. (The reading "the pull is `5 * 10^(-4)`" is unaffected.)

**m3. Typo of a factor `10^6`** (README section 6.3, first bullet): "the full race is `10^(-3) 4^(-n)` at `n = 12`"
must be `10^3 * 4^(-12)` (`D_12 = 989 * 4^(-12)`).

**m4. Outcome item 2: "otherwise about twice and three times the triangle values"** is wrong for `E (-1)^(c_5)`: it is
about three times the *square* of the triangle value (`E_n[{5}] / E_n[{3}]^2 -> 3`; at `n = 11`: `49/16384` against
`1/32`). Section 6.1 (5) has it right.

**m5. "Size `2^(-n/2)`" is claimed for the pentagon alone.** Outcome item 7 ("Triangles and pentagons keep biases of
size `2^(-n/2)` for ever") and section 6.3 ("The first two columns are permanent laws of size `2^(-n/2)`"): the first
two columns of the table are `{3}` and `{5}`, and `E_n[{5}]` is of size `2^(-n)` (`3 * 2^(-(n-1))`, or `2^(-(n-2))` for
`n = 2 mod 4`). The statement holds for `{3}` and `{3,5}`.

**m6. Section 6.2 (5), case `n = 2 (mod 4)`: "the exponent of `-1` is even in both cases".** The exponent of the
summand, `C(q,2) C(n-q,2) + floor(q/4) + floor((n-q)/4)`, equals `(n-2)/4` in both cases (odd when `n = 6 mod 8`); it is
its sum with the prefactor exponent `floor(n/4) = (n-2)/4` that is even. The conclusion is right. For the other rows
the finite-dimensionality argument is rigorous as written in the working-tree version (I checked that both sides are
annihilated by the recurrence with characteristic polynomial
`(x - 2^-4)(x - 2^-8)(x - 2^-12)(x^2 + (17/512) x + 2^-20)` in `k`, `n = 8k + r`; that the five ratios are distinct;
and the identities for all `n <= 240`, thirty members of each class, and by a direct double sum that does not use the
factorised formula for `n <= 100`). Remark: only the three rational ratios actually occur (Hankel rank 3 for
`E[{3,5}]`, 2 for `E[{5}]`, 1 for `n = 2, 6 mod 8`); "at most five" is correct.

**m7. Citation of THM-1470 for the uniformity of the score-parity vector** (README sections 1, 6.1 (4), 6.2 Typing;
THM-5132 status). THM-1470 is klein-2026-07-20-S338, not "opus S409" (that is THM-1460-babai-cameron); it proves the
statement inside a switching class for odd `n` and says for even `n`: "the parity vector is pinned only up to global
complement. **No analogue.**" The statement for all `n` is
`05-knowledge/results/session72_even_graph_unifying_carriers_20261026/session72-even/SCORE_PARITY_QUOTIENT.md`
("Tour/Z ≅ P_n as affine S_n-sets"), and it is re-proved inside the proof of (4), which is self-contained. Fix the
citation.

**m8. The order-12 result is over-read.** Outcome item 6 ("the lead decays faster than `4^(-n)`"), section 7 ("the
data now say only that it is faster than `4^(-n)`"), HYP-9663 (G) ("The decay is faster than `4^(-n)`; nothing more is
known"). One order shows that the plateau `2 * 10^3` does not continue to `n = 12`; it does not establish a rate
(`D_n 4^n = 1706, 2722, 1972, 2541, 1988` already moves by 25% from one order to the next). Suggested wording: "the
registered plateau failed: `D_12 4^12 = 1008 +- 77` (author's and audit samples pooled); the rate is OPEN".

**m9. Section 8, table row 1.** "7 and 21 are missed ... carried by bands `[f(n), max H(n)]` that do not yet overlap
... ends at value 21, order 6": the bands are `[3,3], [5,5], [9,15], [15,45], [25,189]`; non-overlap explains 7 only.
21 lies inside the band of order 6 and is excluded by `21 notin strong(6)`, `f(7) = 25` and `21 = 3 * 7` (THM-1370
(4)-(5); its crevasses are the floor gaps `(5,9)` and `(15,25)`). The fold paragraph says correctly "THM-1370's reason
for 7".

**m10. Section 8, "Folds and bands": "the three peaks are near 5" holds, loosely, for the first sequence only**; the
second peaks at 3 and the third has no peak. `pi^k/k!` peaks at `k = 3` (dimension 6; the integer-dimension volume at 5); `f(n+1)/E[H_n]`
peaks at `n = 3` (`3, 3.33, 3, 2, ...`); the self-converse share bound `(n-1)!! 2^(-n(n-2)/4) = 1, 3/4, 15/64, ...` is
monotone decreasing and has no interior peak, so row 11 ("fold", "ends at 5") describes a threshold crossing (the
labelled share passes 1/2 between `n = 5` and `6`), not "a value realised twice on each side of a peak". The paragraph
is typed ANALOGY, but the factual sentence should be corrected or dropped (the rule it cites, MISTAKE-559, applies).

**m11. Section 8, same paragraph: "its length grows like `1.86 ln h`", inside a sentence introduced as "an exact
form".** `1.86 ln h` (`= 3 ln h / ln 5`) is the heuristic for the *upper end* `B(h)` of the window, from HYP-9520,
section "Heuristic (HEURISTIC)"; it is not in THM-4889 (3), which the sentence cites, and it is not the length (at
`h = 223` it gives 11 against lengths at most 4). Type it HEURISTIC with its source, or drop it.

**m12. Section 8, "Micro and macro": "a quadratic form whose triangle part `sigma . x + e_2(x)` is the same for every
base (Theorem T)".** The triangle part depends on the base through `sigma`; what is the same for every base is `e_2`
and the cube sum `2^ceil(m/2)`.

**m13. "The tower stops at the pentagon" (section 6.3), "The construction stops at the heptagon" (THM-5132 title;
outcome item 7).** What is established is that the direct continuation
`tr(A^3 + A^5 + A^7) = C(|A^3 1|, 2) + floor(n/4)` is false (confirmed: failure rates
`0, .376, .373, .493, .480, .405, .406, .479, .500, .499, .500, .502` at `n = 3..14`, 100,000 samples each, mean 41.8%)
and that `c_7` is not a trace. HYP-9663 (H) (3) rightly keeps the question open; the titles should say "the direct
continuation fails".

**m14. Outcome item 3: "At `n = 4` every other tournament that is not transitive-like votes the other way"** is
unclear. The fact: at `n = 4` the two self-converse classes (transitive, `H = 1`; strong, `H = 5`) carry 48 labelled
tournaments with `H = 1 (mod 4)` and the 16 others (a 3-cycle with a source or a sink) all have `H = 3`.

**m15. "Self-converse tournaments are the majority exactly for `n <= 5`"** (section 4; section 8 row 11) is true for
the labelled fractions (`1, 1, 1, 3/4, 11/16, 55/256`); among classes there is a tie at `n = 4` (2 of 4). Say
"labelled".

**m16. Theorem P (3) can be stated more sharply.** For `m = 1 (mod 4)`: if `K_1` is empty then `A` is singular, so
`|K_0| >= 2` and `G(R) >= 2^((m+1)/2)`; if `K_1` is not empty, `G(R) = 2^((m+1)/2) n_13`. Hence for every `m`:
**`G(R) = 0` or `G(R) >= 2^ceil(m/2)`, and `G(R) = 0` only if `m = 1 (mod 4)`** (observed: `0, 8, 16` at `m = 5`;
`0, 32, 64, 128` at `m = 9`). Suggestion, not an error.

**m17. R2 was a weak test, and the data show a small deviation.** The reference law quoted in R2 is right (exact,
from the eigenvalues of the association scheme of alternating forms, brute-force checked for orders 4, 5, 6:
`0.838981 / 0.159650 / 0.00137`). Run A has `0.838685` and `0.159917`, with standard error `0.000097`: `-3.1` and
`+2.8` standard errors. By the registered rule (tolerance `0.006`, sixty standard errors) R2 is CONFIRMED, and the
note's "to `3 * 10^(-4)`" is accurate; it should add that this difference is three standard errors, i.e. probably
real at order 11.

**m18. Credit lines missing for items that are on record** (all minor; the note does credit THM-5127 for `D_n`):
`D_4, ..., D_7` (THM-466, "P(α_1 odd) climbs 0.250 → 0.297 → 0.375 → 0.448"); `E_5[{5}] = 7/32` (THM-466, "n = 5:
400/1024 tournaments have c5 odd"); the class race at `n = 6` (THM-1980, "n=6: 32 classes `≡1`, 24 `≡3`"); the
unsigned `T = 1/(1 - S)` behind Proposition 5.1 (THM-4816, THM-4817); the matrix `A` over `F_2` and de Caen's rank
bound (THM-4857, same machine); `A^2 1` over the integers is the Kendall-Wei second-order score (THM-896; THM-555:
"Adding `A^2` row-sums does not rescue the distribution"), so "Not found in the repository: the vector `t = As`"
should read "its reduction mod 2 and its weight mod 4". Outside sources worth one line each: `tr A^5 = 5 c_5`
(closed 5-walks are cycles; Savchenko), `H = 1 + 2 alpha_1 + 4 alpha_2 (mod 8)` (Grinberg-Stanley, cited in THM-466).

**m19. THM-5131 lists the residues mod 8 at `n = 10, 11` both under KNOWN (status) and under FINITE-EXACT / in the
title as results of this packet.** Order 10 was recomputed here; order 11 is the second audit of THM-5127 (and is
wrong in the last digit, m2). Make the two places agree.

**m20. HYP-9663 (H) (2) is not falsifiable as written** ("tends to a limit, near 0.2206"; falsifier "a drift of the
fraction at `m = 41, ..., 61`"). State a window (for instance: the fraction lies in `[0.2195, 0.2215]` for every
`m = 1 (mod 4)`, `m >= 17`). My samples of `10^6` bases: `0.2179, 0.2213, 0.2203, 0.2198, 0.2206` at
`m = 13, 17, 21, 25, 37` (standard error `0.0004`); classes: 2 of 12 at `m = 5`, 43,993 of 191,536 at `m = 9`
(both as stated). (G), (A), (H) (1), (H) (3) are well posed and their falsifiers are real.

**m21. Section 8 row 5: the logarithmic density `0.9959` is conditional** (Rubinstein-Sarnak, under GRH and linear
independence of the zeros); the row is typed CITED, add "conditional".

**m22. PREREGISTRATION.md, Run B: "started at 08:02 MST" is inaccurate.** The six output files of Run B and its script
were created at 08:06:07 (the cycle program itself at 08:04:48); the run started at 08:02:02 is the self-converse
census of order 10 (`run_scrace2_10.sh`). Immaterial for the integrity claim (the outputs were written at 08:14:35 -
08:16:20, after the commit of 08:09:37; section C). The registered text must stay unchanged, so add an erratum line
under "Outcomes".

---

## B. What was verified, and how

All programs are in this folder; outputs in `audit/out/`. Conventions: vertices `0..n-1`; in the files of `gentourng`
the character of the pair `(i,j)`, `i < j`, is read as `1` iff `i -> j` (the lists of classes are closed under converse
and every class-level check is run for the tournament and its converse or is converse-invariant, so the other
reading gives the same tables); `A_ab = 1` iff `a -> b`; `X` is the out-neighbourhood of the new vertex.

| # | claim | how | size | result |
|---|---|---|---|---|
| 1 | `H = 1 + 2 c_odd (mod 4)`; OCF `H = 1 + 2 alpha_1 + 4 alpha_2 + 8 alpha_3` exactly | `audit_small.c` (DFS paths and cycles), `audit_census.c` (subset DP for `H`, DP for Hamiltonian cycles of every vertex subset, lists of disjoint odd sets) | all labelled `n <= 7`; all 9,932,002 classes of orders 3-10 | 0 failures |
| 2 | Theorem T (1) exactly, T (2) | `audit_small.c`: cycles by DFS for every labelled tournament; per-base sums | every labelled base of order `1 <= m <= 6` (33,867 bases, 2,131,018 extensions) | 0 failures |
| 3 | T (2), P (3) by direct cube sums | `audit_rand.c p3class`, `p3`: `c_3`, `c_5` counted directly in `R + v_X` for all `2^m` out-sets; `K_0`, `K_1` by trying all `u` | every class of orders 5-8 and its converse (14,808 bases); random bases: 3000, 2000, 1200, 900, 150 of orders 6-10, half of them forced singular (source, sink, twins) | 0 failures; singular `A` in 16 + 76 + 648 + 10,000 class cases |
| 4 | P (1) | `audit_small.c` (labelled `n <= 7`), `audit_census.c` (all classes `n <= 10`, tournament and converse), `audit_rand.c p1` (5000 random tournaments at each `n = 7..12`, direct counts) | as stated | 0 failures |
| 5 | P (2), both equalities | `audit_small.c` (walk form by integer matrix-vector products; weight form) | every labelled extension, `n <= 7` | 0 failures |
| 6 | P (3) table, weights of `K_0`, `K_1` | `audit_small.c` | all labelled bases of order `m <= 6` | 0 failures; values `{8,4}` (`m=4`), `{0,8,16}` (`m=5`), `{8,16,32,64}` (`m=6`) |
| 7 | P (4) joint law of `(s,t)` | `audit_small.c`: full table of counts | `n <= 7` | exact (each admissible `t` has `2^(C(n,2)-(2n-3))` tournaments) |
| 8 | P (5), conditional form | `audit_closed.py`: direct double sum over `(k_1, k_2)` from (1) and (4), `n <= 100`; the note's factorised sum, `n <= 240`; recurrence; exhaustive `n <= 7`; census `n <= 10`; **exact census at `n = 11`** (`audit_census.c -g`) | see m6 | 0 mismatches; `E_11[{3}], E_11[{5}], E_11[{3,5}]`: see below |
| 9 | Section 4: `SC_n`, bounds, shares | `audit_small.c` (brute force over all `n!` permutations, labelled, `n <= 7`), `audit_census.c` (classes, `n <= 10`; backtracking checked against all permutations for `n <= 8`), `audit_misc.py S3` (count per involution, `n <= 6`) | all classes to order 10 | `SC_n = 8, 48, 704, 7040, 380416, 6737920, 958578688, 31596429312`; all table entries agree |
| 10 | Section 5: three strata, class race, self-converse splits mod 4 and mod 8, strong race (both routes), residues mod 8 | `audit_tables.py` on `audit_census.c` | orders 3-10; signed strong counts to order 11 | every entry agrees (82 comparisons); row 11 of the mod-8 table: m2 |
| 11 | Exact `D_n` | `audit_small.c` (`n <= 7`), `audit_census.c` (`n <= 10`), `audit_cube.c census` (orders 3-11 from the bases one order lower) | `D_11` from all 9,733,056 bases of order 10 | `1/2, 1/2, 13/32, 1/4, 853/8192, 1361/32768, 63103/2^23, 650471/2^28`, **`260546929/2^39`** (exact) |
| 12 | Section 6.3 table, all sets `K` | `audit_tables.py` | `n <= 10` (15 sets at orders 9, 10) | every entry agrees, including `E_10[{3,7,9}] = 533081/2^28`, `E_10[{9}] = 1050641/2^28`, the two product identities and "no other" |
| 13 | "`G` not constant at order 5", tower at `k = 3`, zero fraction | rows 3, 6; `audit_rand.c tower`, `zero` | see m13, m20 | as stated |
| 14 | Theorem S (1), (2); deck remarks; `Y_k` | `audit_misc.py` | classes `n <= 8`; `n <= 7` | 0 failures |
| 15 | Kernel-size remark (working tree) | `audit_rand.c kernel`, `audit_closed.py (g)` | exhaustive `m <= 7` (`2, 4, 14, 104, 1728, 60416, 4161536`) | agrees; second derivation of the row `n = 1 mod 4` checked to `n = 197` |
| 16 | Section 7 from the raw files; verdict rules | `audit_mcraw.py` | six files of Run A, four of the control, eight of Run C | all summary numbers reproduced; verdicts follow the registered rules |
| 17 | Order 12, independent estimator | `audit_cube.c mc`: `H(R + v_X) mod 8` for all `2^11` out-sets from the three-case path decomposition, own DP, xoshiro256**; self-test against a direct DP count on 33,344 (base, out-set) pairs of orders 3-11 | 9,000,000 bases of order 11, seeds 77120001-3 | **`D_12 4^12 = 1039 +- 125`** (runs 1266, 949, 902); author `988.5 +- 98.5`; difference `+0.32` combined s.e.; pooled `1008 +- 77` |
| 18 | Same estimator one order lower | `audit_cube.c mc` | 4,000,000 bases of order 10, seed 77110001 | `D_11 4^11 = 2038.9 +- 66.8`; exact `1987.8`; `+0.76` s.e. |
| 19 | Section 8 numbers | `audit_sec8.py` | - | all quoted numbers agree; m9-m11 |
| 20 | Normal form on 5 vertices | `audit_anf5.py` | `2^10` | 120 monomials: 3, 25, 2, 20, 10, 60 |
| 21 | THM-555 ensemble against Theorem T | `audit_tiling.py` | `n = 3..8` exhaustive | equal for `n >= 4`; `0` against `1/2` at `n = 3` |

Further numbers obtained on the way (not in the note): exact `rho_8, rho_9, rho_10 = 1.2441, 1.1558, 1.0406` and
`E|b_R| = 0.022563` at order 10 (agree with the values quoted from THM-5127); relative masses mod 8 at `n = 12` from my
sample: `1.00007, 0.99996, 1.00006, 0.99992`.

Exact order 11 for the triangle and pentagon signs (`audit_census.c -g`: for every one of the 9,733,056 bases of order
10, parities of `c_3`, `c_5` by direct DP and the increments on all `2^10` out-sets by walk counts): the labelled sums
are `1125899906842624`, `107752139522048`, `2214416418340864` of `2^55`, i.e. **`E_11[{3}] = 1/32`,
`E_11[{5}] = 49/16384`, `E_11[{3,5}] = 1007/16384`**, as in the note. A further sample of `2 * 10^8` tournaments at each
order 11-18 through the score formulas (`audit_q2mc.c`, `out/q2mc_summary.out`): all 24 comparisons with the closed
forms within `2.03` standard errors (Run C had `-2.28` at `n = 11`; mine has `-1.72` there, and the exact census
settles it).

## C. Pre-registration integrity

Established from git (`git log --follow -p`, `git show --stat`, diffs between the three versions):

- `80eaf478f0` (author and commit date 08:09:37): `PREREGISTRATION.md` with R1-R6 (51 lines), the programs
  `biasmc.c`, `cycrace.c`, `scrace2.c`, ...; **no file under `out/`**. It is an ancestor of another agent's commit
  `661a2a7e52` (commit date 08:17:16), so it was on the shared remote by then.
- `6070103767` (08:25:37 / 08:25:44): +26 lines (R7, R8) and `c35mc.c`; no outputs. Ancestor of another agent's
  `c76a6c4d51` (08:49:02).
- `aa0edea9c1` (09:16:11 / 09:16:40): all outputs (`out/mc11_*.out` Run A, `out/cyc_10.txt` Run B, `out/c35_*.out`
  Run C) and the Outcomes table (+30 lines).
- The registered text is unchanged: the first 51 lines of the second version and the first 77 lines of the third are
  byte-identical to the earlier versions; the working-tree file equals the committed one.
- The registered programs are unchanged between registration and the final commit (`biasmc.c`, `cycrace.c`,
  `c35mc.c`: empty diffs).
- So each prediction block was committed, and visible to other agents, before the outputs of its run existed in the
  repository.
- The Outcomes table reports each verdict by the registered rule (recomputed from the raw files, `audit_mcraw.py`):
  R1 point estimate `988.5`, `5.19` standard errors outside `[1500, 2500]`: FAILED by rule; R2 `|0.838685 - 0.83898| =
  0.0003 <= 0.006` and `0.001399 in [0.0012, 0.0030]`: CONFIRMED (see m17); R3 `1.01719 in [1.015, 1.045]` (`5.6`
  run-scatter standard errors above the lower end): CONFIRMED; R4 `0.50399 in [0.500, 0.512]`: CONFIRMED; R5
  `533081/2^28 != 650471/2^28` and R6 `1/256, 1/16` against `23/2048, 263/2048`: FAILED (both recomputed exactly by my
  census); R7 `+1.13`, `-0.52` standard errors; R8 `+0.73`, `+0.60`; Run C control: largest deviation `+2.68` at
  `n = 13`; the 16 comparisons with the exact values: largest `2.28`.

Not verifiable from git, but corroborated by the file system. The session's scratch directory still holds the
author's run scripts and raw outputs; each script writes to a file created at launch and renamed at the end, so the
creation time of an output file is the start of its run and its modification time the end (`stat`; these times are
not tamper-proof, but they are independent of the narrative):

| run | script and output files created | output written | registration commit |
|---|---|---|---|
| in-sample: `scrace2` order 10; cycle table orders 3-9 | 08:02:02; 08:04:48 | 08:03:49; 08:05:01 | - |
| Run B (`run_cyc10.sh`, six parts) | **08:06:07** | 08:14:35 - 08:16:20 | `80eaf478f0`, 08:09:37 |
| Run A (`run_mc11.sh`, six seeds) | 08:09:59 | 08:30:21 - 08:30:28 | `80eaf478f0`, 08:09:37 |
| Run C (`run_c35.sh`, eight orders) | 08:26:05 | 08:31:18 - 08:52:38 | `6070103767`, 08:25:37 / 08:25:44 |

- Run B was running at the first commit and wrote its output five to seven minutes after it: the statement "no output
  existed when this text was written" is supported. (My timing of `cycrace.c`, 0.16 ms per class, about 26
  CPU-minutes in all, is compatible with the observed 8.5 to 10 minutes for six parts on a loaded machine.) The start
  time given in the file ("08:02 MST") is not: see m22.
- Run A was started 22 seconds after the first commit, Run C 21 seconds after the second: R1-R4 and R7-R8 were
  registered before their runs began. Run B had ended (08:16:20) before the second block was written, and Run A was
  still running then, as the file says.
- The in-sample numbers of the first block existed before it (08:03:49 and 08:05:01). I recomputed them: `D_n 4^n`,
  the rigid remainders `1417.5, 1181.2, 1965.7, 1653.3`, `rho_8, rho_9, rho_10 = 1.244, 1.156, 1.041`, the values at
  `n = 6, 7, 8, 9` all agree.
- Not verifiable at all: that no earlier, unrecorded run of the same quantities was looked at. The outcomes speak
  against it: R1, R5 and R6 failed, and R7, R8 are now theorems.

## D. Prior art

Repository: full notes in `prior_art_repo_notes.md` (patterns, quotes, line numbers). Overlaps: M1, M2, m1, m7, m18.
Searched and **not found** in the repository: (T) (2); the formula `c_3 = C(q_1,2) + floor(n/4)` in this form;
Theorem P (1) for `c_3 + c_5` and `c_5`, (2) second equality, (3), (4), (5) and the fractions `13/32` (as a pentagon
law), `67/256`, `263/2048`, `1007/16384`, `11/256`, `23/2048`, `49/16384`; the corollary "`H mod 4` is decided by
`q_2 mod 4` for `n <= 6`" (THM-1980 records the exhaustive fact, THM-500 the trace formula; neither has the score form;
HYP-6935 [ii] has a partial law at `n = 5, 6`, see the notes file); the three strata; the class race beyond `n = 6`; the signed strong counts; the table `E_n[K]`;
`chi(T) = prod_v chi(T - v)`; the "square term" reading. THM-172 / THM-173 (`c_5` is determined by the lambda graph)
are related only. Citations checked and accurate: THM-466 (digit lemma, `alpha_2` = pairs of disjoint odd cycles),
THM-478, THM-552, THM-791 (fixed-point-free involutive anti-automorphism), THM-4922 (4), THM-1980 (`d(6) = 2`),
THM-500, THM-118, THM-5127 (fractions, bound, Corollary 4.6, gate), HYP-9355, THM-1370, THM-4889 (3), THM-4026 /
THM-4028. Inaccurate: THM-1470 (m7), the source of `1.86 ln h` (m11).

Literature and OEIS (generic queries only; `literature_notes.md`): the formula for `c_3` is classical (Moon, *Topics
on Tournaments*, section 5); no statement of the parity law, of `E (-1)^(c_3)`, of the parity of `c_5`, or of the
kernel formula was found; `tr A^5 = 5 c_5` is in the literature; de Caen's rank bound is the classical neighbour of
the kernel statements. OEIS: the class counts are A002785; the labelled self-converse counts
`8, 48, 704, 7040, 380416, 6737920, ...`, the signed sums `1, 2, 4, 32, 416, 8192`, the signed strong counts and the
class race are not there. The typing "no claim about the literature" is appropriate.

## E. Typing and wording

Status lines match what is established, except: THM-5131 "NEW here" and "(S) 2. (PROVED.)" (M1, M2), the label
"exact" of the row `n = 11` (m2), the sentences listed in m3-m5, m8-m15. The pre-registered outcomes are reported
honestly (three failures are stated as failures, in the outcome list, the titles and the hypothesis file). HYP-9663:
see m8, m20. Numbers in the prose agree with the tables and the outputs wherever not flagged above. The placeholders
`AUDITLINE` (README) and `AUDITSTATUS` (THM-5131, THM-5132, HYP-9663) are still to be filled. The uncommitted
additions of 09:22-09:25 (kernel-size remark, the two check programs, the paragraph "The period-8 phase ...") are
correct as far as they state identities (`Re (1+i)^n = (2/n) 2^((n-1)/2)` for odd `n`; `E|ker A|`); the rest of that
paragraph is typed ANALOGY and claims no map.

## F. Not checked

- The 120,000,000 random bases of orders 11-40 of `k35cube.c` and the author's censuses were not re-run with the
  author's programs; the same statements were checked with mine on smaller random samples (orders 6-10 by direct
  sums) and on all classes to order 10.
- `D_12` is NUMERICAL (two independent samples); no exact value exists.
- The three-strata split and the class race at order 11.
- Section 1 (inheritance pass) was checked only through the citations listed in D; the claim that no same-day
  session treats `H mod 4` was spot-checked by pattern only.
- Full texts of some external papers (see `literature_notes.md`).

## G. Reproduction of the audit (from this folder; binaries in a scratch directory `$S`)

    cc -O2 -o $S/audit_small audit_small.c;  $S/audit_small 7 1 > out/audit_small_7.out          # about 2 minutes
    cc -O2 -o $S/audit_census audit_census.c
    for n in 3 4 5 6 7 8; do gentourng -q $n | $S/audit_census -b $n > out/census_$n.txt 2> out/census_$n.err; done
    gentourng -q 9 | $S/audit_census 9 > out/census_9.txt 2> out/census_9.err
    for r in 0 1 2; do gentourng -q 10 $r/3 | $S/audit_census 10 > out/census_10_part$r.txt 2> out/census_10_part$r.err & done; wait
    for r in 0 1 2; do gentourng -q 10 $r/3 | $S/audit_census -g 10 > /dev/null 2> out/cube35_10_part$r.out & done; wait   # exact order 11, c3 / c5
    cc -O2 -o $S/audit_cube audit_cube.c
    for m in 3 5 7 9 11; do $S/audit_cube selftest $m 40 $((700+m)); done                        # (m = 11: 3 bases)
    for m in 2 3 4 5 6 7 8 9; do gentourng -q $m | $S/audit_cube census $m; done                  # exact D_3 .. D_10, residues mod 8
    for r in 0 1 2; do gentourng -q 10 $r/3 | $S/audit_cube census 10 > out/cube_census_10_part$r.out & done; wait        # exact D_11
    for s in 1 2 3; do $S/audit_cube mc 11 3000000 $((77120000+s)) > out/mc_audit_11_$s.out & done; wait                 # order 12, 14 minutes each
    $S/audit_cube mc 10 4000000 77110001 > out/mc_audit_10_1.out
    cc -O2 -o $S/audit_rand audit_rand.c -lm                                                    # modes p1, p3, p3class, tower, zero, kernel
    cc -O2 -o $S/audit_q2mc audit_q2mc.c
    python3 audit_tables.py > out/audit_tables.out;  python3 audit_closed.py 100 240 > out/audit_closed.out
    python3 audit_mcraw.py > out/audit_mcraw.out;    python3 audit_misc.py > out/audit_misc.out
    python3 audit_sec8.py; python3 audit_anf5.py; python3 audit_tiling.py; python3 audit_kbs.py

Outputs of every run quoted above are in `out/` (`audit_small_6.out`, `audit_small_7.out`, `census_*.txt/.err`,
`cube_census_*.out`, `cube35_*.out`, `cube_selftest.out`, `mc_audit_*.out`, `rand_p1.out`, `rand_p3.out`,
`rand_p3class.out`, `rand_tower.out`, `rand_zero.out`, `rand_kernel.out`, `q2mc.out`, `q2mc_summary.out`,
`audit_tables.out`, `audit_closed.out`, `audit_mcraw.out`, `audit_misc.out`, `audit_sec8.out`, `audit_anf5.out`,
`audit_tiling.out`, `audit_kbs.out`). The two notes files `prior_art_repo_notes.md` and `literature_notes.md` hold the
search patterns, quotations and sources of part D.
