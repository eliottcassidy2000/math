# AUDIT - kind-pasteur-2026-11-11-S25 (Euler fields in twisted de Rham classes; one-apex chain semigroup; addendum D)

Independent adversarial audit, own code, exact arithmetic. Auditor: sub-agent of the session, no access to the
author's reasoning beyond the committed files.

- Audited commit for parts A and B: `8b1398bae9` (README, PREREGISTRATION, THM-5147, THM-5148, HYP-9668 parts
  (E'), (Q)).
- Addendum D (added to the scope during the audit): PREREGISTRATION section "Addendum D" (`870ad0c5d7`), its
  outcomes (`d67741ce21`), `scripts/direct_joins.py`, `runs/d1_bruteforce.txt`, `runs/direct_closure_*.txt`, and
  the statements that carry it at `7319c023bc` (README B4, THM-5148 (D) and (E) first item, HYP-9677 (A)).
- The packet moved under the audit (checkpoint 2, `7319c023bc`, 19:02). Locations are given as
  `README@8b13 line n` (audited commit) and `README@7319 line n` (HEAD when this file was written; Part A lines are
  shifted by +4 there).
- Not in the assigned scope: addendum E (README B5, HYP-9677 (B), THM-5148 (E) items 2-3). Its numbers were
  reproduced as a by-product (section 6.6); its statements and credits were not audited.

Python 3.13.14, sympy 1.14.0, numpy 2.3.2, numba 0.64.0. The author's scripts were run with `NUMBA_CACHE_DIR` in a
scratch folder and `PYTHONDONTWRITEBYTECODE=1`; nothing was written outside `audit/`.

## 0. Verdict

**PASS WITH CORRECTIONS.** No proved statement of part A is false except one remark; every computed number of
part B and of addendum D reproduces with independent programs. Five major corrections: one false remark (A), one
mislabelled headline number (B), and three credit / novelty errors (B, D): the family of addendum D and its list
of 29 misses are already in canon, to `10^9`, in an audited packet of the same agent (HYP-9494, 2026-11-07).

| statement | verdict |
|---|---|
| Lemma A1 | PASS |
| Theorem A2 (1), (2), (3), all `l >= 1` | PASS (upgraded: U1, U2) |
| Corollary A3 | PASS; the remark after it: PASS WITH CORRECTIONS (m2, m3, m4) |
| Proposition A4 | PASS |
| "Uses" (i), (iii) | PASS |
| "Uses" (ii) (`UV - c` is never a null element) | **FAIL as written** (M5): needs `[1]_z != 0` |
| THM-5147 | PASS WITH CORRECTIONS (m1, m7; add U1-U3) |
| HYP-9668, parts (E'), (Q) | PASS WITH CORRECTIONS (m1, m4) |
| Proposition B1 (both formulas) | PASS (mathematics); credit: KNOWN (M2) |
| Corollary B2 | PASS WITH CORRECTIONS (m9; credit M2) |
| C1, C2 CONFIRMED; C3, C4 FAILED | PASS: labels and numbers right (independent closures) |
| C5 CONFIRMED, "mean multiplicity 82.9" | **FAIL as a multiplicity** (M3): 14.3 distinct chains per odd number |
| Reading of B3 (`W >= A`, `N >= 3 A_1 A_2`) | PASS; the sentence on the bank: PASS WITH CORRECTIONS (m11) |
| B1 credits | accurate for what is cited, **incomplete** (M2) |
| B4@8b13 (= B6@7319) remarks | PASS WITH CORRECTIONS (m12, M3) |
| Pre-registration integrity (C1-C5, D1-D3) | PASS (section 5.5), with the caveat M1 (iii) for D2 |
| D1 (law of direct joins; Proposition B3@7319) | PASS (mathematics); credit: KNOWN (M1, M2) |
| `direct_joins.py` | PASS: no bug found (section 6.2) |
| 29 misses below `2^13`, `2^16` | PASS (numbers); **novelty FAIL**: KNOWN to `10^9` (M1) |
| HYP-9677 (A) | restates the statement of HYP-9494 (D), with a smaller verified range (M1, M4) |

## 1. MAJOR findings

### M1. The family of addendum D is the O/Z/J root-change grammar; its 29 misses are KNOWN to `10^9`

Location: README@7319 lines 22-32 ("Added here: (1) ... (2) FINITE-EXACT, pre-registered: with no bank ... 29
numbers"), B4 lines 200-238; THM-5148@7319 title (D), (E), status lines 5-15, body line 72; HYP-9677 (A) lines 5-7,
13-15, 36; PREREGISTRATION lines 80-119.

Finding. "Apex-path tournaments" (iterated direct joins of one-apex blocks at their apexes) are exactly the words
of the O/Z/J grammar of Sessions 197-200. THM-4904 (kind-pasteur-2026-11-07-S14), Setting: words over `{O, Z, J}`
act on `(D, E, U, d)` from the rooted edge `(1, 0, 2, 1)` by
`O = (D, 0, 2U, d)`, `Z = (D - E + U, d, 2U - d, d)`, `J = (D - E + U, D - E, D + U, D)`; "Each is a literal vertex
birth, and J is a root change." Dictionary: `J` = a new apex (the block with the empty word), `O` = a spine letter
1, `Z` = a spine letter 0, the rooted edge = the block with word `1`; `(D, E, U) = (A, d, W) = (H, o - e, o + e)`
and the fourth coordinate is `H(T - root)`. (Proof: each letter is the join law for one new vertex; `E' = d`
because a path ending at the root must pass `z -> root` last; `U' = 2U - d` because every path from the root with
at least two vertices can be prolonged by `z`.) Checked: 29,524 words of length at most 9 against the literal
tournaments built from the blocks; 797,161 words of length at most 12 against the products `M_(w_k)...M_(w_0)`;
the closures of both descriptions below `2^11` have the same 28 misses (`d4_ozj_identification.out`, Z1-Z3). The
strong members of the family are uniquely Hamiltonian (9,840 tournaments checked, Z5; in general
`hc(J) = s_a(A) e_b(B)`, THM-4866 (B)), and Session197 (ROOTED.md Theorem 2, C1_COMPLETE.md Theorem 4) proves that
every uniquely Hamiltonian strong tournament arises this way.

The miss set is in canon:
- `session197/lowcycle/C1_COMPLETE.md`, section 3: the missing odd values of the C = 1 spectrum through 5,001 are
  exactly 7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583,
  637, 821, 1067, 1339, 1861, 2647;
- THM-4884 (2b): "the C1 readout misses exactly Session 197's 29 values through 50,001", for the profiles of "the
  physical O/Z/J grammar from (1, 0, 2, 1)";
- HYP-9494 (D) and THM-4904, status (kind-pasteur-2026-11-07-S14, audited): "the C = 1 spectrum misses exactly
  Session197's 29 values below 10^9", with the root-change depth of every value (maximum 7, at 1877 and 1897),
  `F_<=2` (three blocks) missing exactly 140 odd values below `10^9`, and the statement that every odd `N > 2647`
  is the path count of a uniquely Hamiltonian strong tournament, derived there from the conjecture that two root
  changes suffice above 3,780,967, with the remark "that would prove the Chan-Pak conjecture".
(THM-4884 concerns all words; Session197 and HYP-9494 the strong ones. The two miss sets coincide in the ranges
computed, and the strong statement implies the one for all words.)

Consequences.
(i) The result of addendum D (29 misses below `2^13`, `2^16`, `2^18`) is a re-computation of a sub-range of a
KNOWN, audited statement; it is not an addition of this session.
(ii) Law D1 / Proposition B3 is the block form of the KNOWN letter moves (see M2).
(iii) D2 ("misses exactly 7 and 21", confidence 40 per cent) was registered against a fact already in canon, and
D3 (13 is missed, at most 50 misses) was readable from THM-4884. They are not out-of-sample tests. The
registration itself is in order (section 5.5).
(iv) HYP-9677 (A) is implied by, and in every computed range identical to, the statement of HYP-9494 (D); it is
filed again with evidence to `2^18` instead of `10^9`, and its list of "related families in canon" omits HYP-9494,
THM-4904, THM-4888, HYP-9480, THM-4884.
(v) The sweep quoted in the README ("sub-agent digest") missed the author's own S13 and S14 packets. A grep for
the first misses ("7, 13, 21, 41, 43, 55") finds THM-4884 at once.

Proposed text, README B4, in place of "**Result (FINITE-EXACT; ...)**":
"**Result (KNOWN; re-computed here).** This family is the O/Z/J root-change grammar of Sessions 197-200 (THM-4904,
Setting: J = a new apex, O / Z = a spine letter 1 / 0, `(D, E, U) = (A, d, W)`); its strong members are the
uniquely Hamiltonian strong tournaments (Session197). That it misses exactly the 29 odd numbers 7, 13, ..., 2647
is Session197's C1_COMPLETE.md (to 5,001), THM-4884 (2b) (to 50,001) and HYP-9494 (D) (to `10^9`, with the
root-change depth of every value). `scripts/direct_joins.py` re-computes the list below `2^18` in block
coordinates; the audit re-computed it below `2^16` with an independent program, and below `2^11` with two more,
one of them without any pruning. D2 was registered without this knowledge: its outcome was already in canon."
Proposed text, README section 0, part B: replace "Added here: (1) ... (2) ..." by "Re-derived here in block form
(KNOWN: SEAM.md, THM-4875 (B), THM-4904): ...; re-computed (KNOWN: HYP-9494 (D)): ...; new: the coverage data of
the pure chains (B3) and of finitely many blocks (B5)."
HYP-9677: keep (B) only; state (A) as "HYP-9494 (D), FINITE-EXACT to `10^9`".

### M2. Proposition B1, Corollary B2 and Proposition B3 are block forms of KNOWN laws; "not found in this form" is wrong

Location: README@8b13 lines 21-25 ("What is added: ... a product of 2 x 2 integer matrices ... a matrix
semigroup"), 144-169; THM-5148@8b13 status lines 5-8; THM-5148@7319 status lines 5-9 ("(A) and (D) add the
transport of the root data to the new apex, which makes iterated joins a matrix semigroup. Not found in this form
in one sweep of the repo ...; no novelty is claimed beyond the formulation").

Finding. The formulation is the standard one in canon:
- `05-knowledge/results/seven_return_native_recursions_20261029/SEAM.md`, (1)-(7): every block `Y` with arbitrary
  entry and exit ports, "including coincident ports", has an integer 3 x 3 matrix with
  `v(X.Y, exit(Y)) = M_Y v(X, a)`, `M_Y e_0 = v(Y)`, `e_0 = (1, 0, 0)^T`, `M_(X.Y) = M_Y M_X`, and the generating
  function `e_0^T (I - A(z))^(-1) e_0` of a block alphabet. Proposition B3 is the case "one-apex block, both ports
  at the apex".
- THM-4875 (B): "(h, e, R)(N, z) = M(G, u, z)(h, e, R)(T, v)" for a gadget `G` with entry `u` and exit `z`; the
  law is labelled KNOWN there (SEAM.md), and closed entries of `M(G, u, z)` are given.
- THM-4904, Setting: the letter moves O, Z, J for exactly this family (M1); `M_w` is the product of the letter maps
  of the block (797,161 words checked).
- THM-4888 (D) 1 (Session205): `JZ = OJ`. This is the normal form "leading letters 0 of a word are absorbed in the
  exponent" (and, read the other way, "trailing letters 1 of a block are chain vertices", which the note does not
  state; see M3). Checked as maps on states (`d4`, Z4).
- THM-5130 (A) (the author's S22, one day earlier): states transported by explicit integer matrices through
  gadgets, rays `r + 2^t M` from inserted transitive chains.
The pure chains of B2-B3 are the O/Z/J words in which every `J` is immediately preceded by `O`: the closure of
that sub-language computed from the letter moves misses the same 414 odd numbers below `2^12` (`d4`, Z4).
Proposition B1 itself is THM-4866 (B) applied to `(X => C, u)`, whose triple is `(H(X), 2^i W+_X(u), 0)`, plus
the third row of the seam matrix of the block `C => Y`.

So even "no novelty beyond the formulation" claims too much. What the note adds: the closed entries of a whole
one-apex block (`M(p, i)`, `M_w`), and the finite coverage data of the sub-language "every root change is preceded
by a sink" (C2-C4), for which no earlier record was found by this audit (a search, not a clearance).

Proposed text, THM-5148 status: "(A), (D) KNOWN, re-derived in block form: the transport of the rooted triple
through a block by an integer matrix and the semigroup law are the universal seam transfer (SEAM.md (1)-(7);
THM-4875 (B)); for this family the letter form is the O/Z/J grammar (THM-4904, Setting), and the absorption of
leading letters 0 is Session205's JZ = OJ (THM-4888 (D) 1). New: the closed entries of one block. (B) the same.
(C) FINITE-EXACT, pre-registered: coverage of the sub-language in which every root change is preceded by a sink.
(E), first item: KNOWN (HYP-9494 (D), to `10^9`); re-computed to `2^18`."

### M3. "Mean multiplicity 82.9" counts each tournament many times; the number of distinct chains per odd number is 14.3

Location: README@8b13 lines 27, 184, 193, 206 (@7319: 30, 188, 197, 266); THM-5148@8b13 title "(C) ... with mean
multiplicity 82.9" and line 55; THM-5148@7319 line 54 ("10,866,860 products (mean multiplicity 82.9)");
PREREGISTRATION line 66; `scripts/chain_semigroup.py` lines 159-206 (`cnt2`, `cnt3`).

Finding. `depth_cover` counts parameter tuples, and the parametrisation is redundant in three ways: the first and
the middle block run over all values `p` (trailing letters 1 included), and the last block over all `q` (leading
letters 0 included). Trailing letters 1 of a block and leading letters 0 of the next one are both chain vertices
(`M(p', i') M(2p + 1, i) = M(p', i' + 1) M(p, i)`: 9,600 cases, `b2`, K4; canon: `JZ = OJ`). Hence the same
tournament with chain lengths `c_1` (depth 2), resp. `(c_1, c_2)` (depth 3), is counted `c_1 (c_1 + 1)/2`, resp.
`c_1 c_2 (c_2 + 1)/2` times. With these weights the auditor's enumeration of fully stripped chains reproduces the
author's two counts exactly (`b4_depth_and_multiplicity.out`):

| below `2^18` | distinct chains | tuples of the registered matrices `M(p, i)` (last block without trailing 1's) | author's tuples |
|---|---|---|---|
| depth 2 | 671,731 | 1,352,226 | 2,659,887 |
| depth 3 | 1,199,555 | 4,406,181 | 8,206,973 |
| per odd number | **14.3** | 43.9 | 82.9 |

(Depth 1 adds 181,720 cores, 1.39 per odd number, the KNOWN `2 ln 2`. Even 14.3 counts some tournaments twice: a
depth-2 chain whose second block is a single vertex is a one-apex tournament.) Read literally, the registered
count is infinite (trailing letters 1 of the last block do not change the top-left entry); with the last block
stripped it is 43.9. So the label "C5 CONFIRMED" survives only for tuples, and 82.9 is not a number of
"products" nor a multiplicity of the family. THM-4904 counts words; the note should use the same convention.

The qualitative reading stands with the corrected number: a Poisson model with mean 14.3 would miss a fraction
`6 * 10^-7`, the family misses 5.26 per cent.

Proposed text, B3 table: "C5 (at least 40 tuples per odd number): CONFIRMED for tuples only (43.9 tuples of the
registered matrices; 82.9 in the program's parametrisation). The tuples are redundant: a chain with chain lengths
`c_1, c_2` is counted `c_1 c_2 (c_2 + 1)/2` times. Distinct chains per odd number: 14.3 (depth 2: 5.1; depth 3:
9.2), below the registered 40." Replace "mean multiplicity 83" by "14 chains per odd number" in section 0, in the
Reading, in remark 3, and in THM-5148 (title@8b13, (C)).

### M4. Canon side of M1: THM-5148 (E), first item, and HYP-9677 (A) carry labels that hide a KNOWN result

Location: THM-5148@7319 status line 12 ("(C), (E) FINITE-EXACT ..., pre-registered"), line 72; HYP-9677 status
lines 5-7 ("(A) holds for all odd numbers up to 2^18").

Finding. FINITE-EXACT without KNOWN blends a re-computation with a new result, and the stated range (`2^18`) is
below the range in canon (`10^9`). Proposed: "(E) first item: KNOWN (Session197; THM-4884 (2b); HYP-9494 (D), to
`10^9`), re-computed here to `2^18` and by the audit to `2^16`"; in HYP-9677: "(A): HYP-9494 (D); evidence to
`10^9` there."

### M5. "`UV - c` is never a null element" is false without `[1]_z != 0`

Location: README@8b13 lines 111-113 (@7319: 115-117), "Uses" (ii).

Finding. The identity `[(UV - c)^m]_z = m! e_m(-c) [1]_z` is right (THM-5147 (D) states only the identity). The
conclusion needs `[1]_z != 0`: then `e_1(-c) = 1 - c` and `e_2(-1) = 1/2` exclude every `c` (consecutive
truncations differ by `t^m/m!` and `e_m(0) = 1`). Proposition A4 is stated for every etale map and every centre,
and there the sentence fails: for the double cover of S24, Proposition 6.1, at the centres `(0, b)` the Stein space
is zero, so every `UV - c` is a null element.

Proposed sentence: "(ii) Where `[1]_z != 0` (every centre for an automorphism; every centre algebraically
independent over a field of definition, S24 Theorem 5.4), `UV - c` is not a null element for any `c`:
`[(UV - c)^m] = m! e_m(-c) [1]`, `e_1(-c) = 1 - c`, `e_2(-1) = 1/2`. Where `[1]_z = 0` (the double cover of S24 at
the centres `(0, b)`) every `UV - c` is null."

## 2. MINOR findings

- **m1.** README@8b13 lines 16-17; HYP-9668 status line 22. "(`W^m` lies in it for all `m >= 1`, `x W^m = m! x`
  never does); the same for every `W = x + x^(l+1) y`" is inexact for `l >= 2`: `[W^m] = (1-l)...(m-l) [1] != 0`
  for `m < l` (for `l = 2`: `W + 1` lies in `M_W`). Proposed: "for `W = x + x^(l+1) y`: `W^m` lies in the image
  for `m >= l` and `x^l W^n = n! x^l` never does".
- **m2.** README@8b13 line 90. "the obstruction is absent exactly there": "exactly" is not justified. The
  obstruction is absent for every potential without an Euler field of non-positive integer divergence, for example
  `x^3 + y^3` (every Euler field of constant divergence has `delta = 2/3`) and `x (x^2 y + x + 1)` (U3). Proposed:
  "the obstruction is absent there".
- **m3.** README@8b13 lines 87-88. "two necessary conditions that do not involve 1 being in the image only" is
  garbled (the 1-property is about 1 being in the image). Proposed: "two necessary conditions: the 1-property (1
  is not in the image unless the image is everything) and the Euler obstruction, which can fail with 1 outside the
  image".
- **m4.** README@8b13 lines 92-95; HYP-9668 line 115. (a) "without an Euler field": every polynomial without
  critical point has fields with `xi(W) = W` (`xi = W (p d/dx + q d/dy)`, `p W_x + q W_y = 1`); what is meant, and
  true, is "without an Euler field of constant divergence" (proved in U3; the note asserts it without proof).
  (b) The example is not open: the residue form of Theorem A2 settles it (U2): for `W = x (x^2 y + x + 1)`,
  `W^m` lies in `M_W` for `m >= 2` and `x^2 W^n = n! x^2` does not. Proposed: "Not settled: rank-one potentials
  outside the monomials `x^p y^q` (1-property) and the family `x^(l+1) y + r(x)` (audit, U2), i.e.
  `x^m (x^l y + p(x))^n` with `(m, n) != (1, 1)` in the classification of polynomials with general fibre `C^*`"
  (the classification is used there without a reference; add one).
- **m5.** README@8b13 line 72. The direct check is announced for `l = 1` only; the theorem is for all `l`. The
  term-by-term computation for general `l` is two lines (section 4.1) and was checked for `l <= 4`.
- **m6.** README@8b13 line 65; THM-5147 (B) 3. "`dim B/M_W = 1`" rests on S24's machine control and a cited
  theorem; an elementary proof for all `l` is in U1.
- **m7.** THM-5147 status lines 12-15. Accurate, but say which statement: (B) is an instance of the statement of
  Zhao's Conjecture 1.3 (arXiv:0902.0210: any commuting set of operators `h(z) + sum c_i d/dz_i`), not of the
  Image Conjecture as restated by van den Essen, Wright and Zhao (arXiv:1008.3962, Conjecture 2.1: operators
  `d/dz_i - a_i` with `(a_i)` a regular sequence of `A`), whose special case SIC is what THM-2801 refutes. The
  S25 README has no sentence on this; add one after Theorem A2 so that it is not read as a first counterexample.
- **m8.** README@8b13 lines 150-151. The proof of B1 quotes THM-4876 for an arc-reversal formula; THM-4876 (1) is
  the vertex form. One line is missing (apply (1) at `v` to `T` and to `X => C => Y`, or quote (4)). A direct
  proof without the odd-cycle formula is in U4.
- **m9.** README@8b13 lines 157-159; THM-5148@8b13 line 45. "`W- = R(p) + 1`, `tau = 2^(L(p))`" hold only for a
  word without leading letter 0. For `k` leading letters 0: `W- = 2^k (R(p) + 1)`, `tau = 2^(k + L(p))` (511 words
  checked). The remark on absorption comes ten lines later; state the hypothesis in the sentence.
- **m10.** `scripts/chain_semigroup.py` line 211. `depth()` builds the depth-1 table from `g_list(X)`, i.e.
  `p <= X`; 16,328 one-apex values of `[2^17, 2^18]` are attained only for `p` in `(X, 2X]`. Harmless here: each
  of them is a depth-2 value, the miss counts do not change (`b4`). Use `g_list(2 * X)` as `run_closure` does.
- **m11.** README@8b13 lines 27-28, 193-195; THM-5148@8b13 status lines 14-16. "the joins without an inserted
  chain and the bank ... carry the coverage", "this test shows they cannot be dropped": the test drops both
  ingredients at once and cannot separate them, and canon already said that the bank is not needed (HYP-9494 (D)).
  Replaced at `7319c023bc`; keep the replacement and cite HYP-9494 there.
- **m12.** README@8b13 lines 199-204 (B6@7319). "the first two are KNOWN": in remark 1 the bounded mean is PROVED
  (THM-4885 (C)); "so completeness up to `2^35` ... is not evidence" is a Poisson heuristic, the content of
  HYP-9480 (A), (C) (OPEN). Label it HEURISTIC.
- **m13.** README@7319 line 230; THM-5148@7319 line 12. "every state with `A <= X` is generated": the program
  discards states with `o` or `e` above `X/2`. What is exhaustive is the set of values (section 6.2). Proposed:
  "every value `A <= X` of the family is generated".
- **m14.** README@8b13 lines 128-133. THM-4885 (D) at `2^35` is "one run, not independently recomputed" (`2^30`
  audited): say so. The two-apex statement below `2^23` is S93's registered W1 (commit `057f7b27c8`), not part of
  THM-5141; S94 has since reached 226,079,703.
- **m15.** `scripts/direct_joins.py` line 58 (and `finite_blocks.py`). The hash multiplier `0x9E3779B97F4A7C15`
  exceeds `2^63`; the code relies on numba's integer typing and wrap-around. No effect on correctness (any
  deterministic hash works; the three independent closures of the audit use Python and numba sets).
- **m16.** README@7319 line 222 ("Along every orbit `0 <= d <= A <= W`"): true after the first block (the virtual
  state `(1, 0, 0)` has `W < A`), and KNOWN: Session200's inequalities `D - E >= d >= 1`, `D <= U` (quoted in
  THM-4904).

## 3. Upgrades (proved or computed by the audit; the author may use them)

- **U1 (codimension one, elementary, all `l` and `r`).** Let `W = x^(l+1) y + r(x)`, `l >= 1`, `r` any polynomial.
  Then `B = C x^l + M_W`. Proof. `(d/dy - W_y)(x^k y^j)` gives `x^(k+l+1) y^j = j x^k y^(j-1)`, so `B/M_W` is
  spanned by `x^k y^j`, `k <= l`. Let `V_j` be the span of these with `y`-degree at most `j`. For `1 <= k <= l`,
  `(d/dx - W_x)(x^k y^j)` and the previous relation give
  `(k - (l+1)(j+1)) x^(k-1) y^j = r'(x) x^k y^j`, and the coefficient is not zero; by descending induction on `k`,
  `x^(k-1) y^j` lies in `C x^l y^j + V_(j-1)`. For `k = 0`: `x^l y^(j+1) = -r'(x) y^j/(l+1)` lies in `V_j`. So
  `V_j = V_0 = C x^l` modulo `M_W`. With `lambda(x^l) != 0` the codimension is one and `M_W = ker lambda`. No
  cited comparison theorem is needed for Theorem A2.
- **U2 (the family `x^(l+1) y + r(x)`, `r(0) = 0`: no Euler field needed).** The form
  `lambda_r(x^k y^j) = j! [x^((l+1)(j+1)-1-k)] e^(-r(x))` vanishes on `M_W` (same proof as A2 (1)), and
  `lambda_r(x^k W^m) = [x^(l-k)] (m! e^(-r) e_m(r))`, `e_m` the truncated exponential; since
  `m! e^(-s) e_m(s) = m! + O(s^(m+1))`: `W^m` lies in `M_W` for `m >= l`, `x^l W^n = n! x^l` modulo `M_W`, and
  `lambda_r(x^l) = 1`. **So `M_W` is not a Mathieu-Zhao subspace for any `W = x^(l+1) y + r(x)`**, in particular
  for `x (x^2 y + x + 1)`, which the note lists as not settled. `1` lies in `M_W` iff `[x^l] e^(-r) = 0`. (If
  `r(0) != 0`, replace `W` by `W - r(0)`: `M_W` depends on `dW` only.) Exact certificates (vector fields over `Q`)
  for five potentials and the moment formula on 175 cases: `a3`, `a5`.
- **U3.** `x (x^2 y + x + 1)` has no vector field with `xi(W) = W` and constant divergence: it would give
  `[W] = delta [1]`, `[W^2] = delta (delta + 1) [1]`; but `W - 1` and `W^2` lie in `M_W` and `[1] != 0`
  (`lambda(1) = -1/2`), so `delta = 1` and `2 [1] = 0`.
- **U4 (B1 without the odd-cycle formula).** A Hamiltonian path of `T` either avoids `v -> u` (it is a path of
  `X => C => Y`: `H(X) H(Y)`), or is `(A_X)(A_C)(A_Y -> v) -> (u -> B_X)(B_C)(B_Y)` with `B_X` a path of `X` from
  `u`, `A_X` a Hamiltonian path of the rest of `X`, `A_C + B_C = C` (`2^i` splittings), `A_Y` a path of `Y` ending
  at `v`, `B_Y` a Hamiltonian path of the rest: `2^i W+_X(u) W-_Y(v)`.
- **U5 (`W >= A`).** Cutting a Hamiltonian path of `X` at `u` maps it injectively to a pair (path from `u`,
  Hamiltonian path of the complement): `W+_X(u) >= H(X)`, and dually. Checked on all 202,013 labeled rooted
  tournaments with at most 6 vertices; equality is not confined to sink roots (every vertex of the 3-cycle).
- **U6.** Corollary A3 does not need an eigen-polynomial: any `h` with
  `[(xi + delta)(xi + delta + 1)...(xi + delta + N - 1) h] != 0` for infinitely many `N` works. In applications
  `[h] != 0` must come from a linear form (as `lambda` in A2), not from a normal form relative to a truncation.
- **U7.** Explicit witnesses (adjacency rows, literal recount) for the 27 misses other than 7 and 21:
  `d3_witnesses_for_the_29.out`.
- **U8.** `x_e = o - e` (README B1, B4): the alternating sum over the paths from the root equals the number of
  Hamiltonian paths ending at the root, on all 202,013 rooted tournaments with at most 6 vertices (`b6`).

## 4. Part A in detail

### 4.1 Proofs, line by line

Lemma A1. `(d/dx - W_x) p + (d/dy - W_y) q = div(v) - v(W)`; with `v = h xi`: `xi(h) + (div xi - xi(W)) h`.
Powers: `xi(W^(m-1) h) = W^(m-1) (xi + m - 1) h` and induction give the product
`(xi + delta)...(xi + delta + m - 1)`. Correct. On a smooth affine surface with a volume form `omega`:
`(d - dW ^) i_v omega = (div v - v(W)) omega` and every 1-form is `i_v omega`. Correct.

Theorem A2 (1), all `l`. With `N = (l+1)(j+1) - k`:
`lambda((d/dx - W_x) x^k y^j) = (-1)^N (j!/N!) (k + N - (l+1)(j+1)) = 0` for `N >= 1`; for `N = 0` it is
`j! (k - (l+1)(j+1)) = 0`; for `N < 0` all three terms vanish. `lambda((d/dy - W_y) x^k y^j)`: both terms have the
same index `n = (l+1) j - 1 - k` and cancel (`j >= 1`); for `j = 0` the index is negative. Machine: 672 generators
(`l <= 4`, `k < 14`, `j < 6`), 0 failures; `lambda(1) = (-1)^l/l!`, `lambda(x^l) = 1` (`l <= 8`); the formula
against sympy residues, 108 cases. The integral description in the README is correct (the inner integral of a
`y`-derivative is `-h(x, 0) e^(-x)`, entire in `x`).

Theorem A2 (2). `xi = x d/dx - l y d/dy`, `xi(W) = W`, `div xi = 1 - l`: correct; `lambda(W h) = (k - l j + 1 - l)
lambda(h)` on 336 monomials. Independent certificates (`a2`): polynomial vector fields over `Q`, found by exact
Gauss-Jordan elimination and re-expanded, for `l = 1`: `W^m` (`m <= 4`), `x W^n - n! x` (`n <= 4`); `l = 2`:
`W + 1`, `W^2`, `W^3`, `W^4`, `x^2 W^n - n! x^2` (`n <= 4`); `l = 3`: `W + 2`, `W^2 - 2`, `W^3`, `W^4`,
`x^3 W^n - n! x^3` (`n <= 3`). Controls `1`, `x^l`, `x^l W^2`: no solution in the box; non-membership is proved by
`lambda`.

Theorem A2 (3) and the definition. The definition used ("if `P^m` lies in `M` for all `m >= 1`, then for every
`h`, `h P^m` lies in `M` for all large `m`") is Zhao's. `P = W^l`: `P^m = W^(lm)`, `lm >= l`, so every power
lies in `M_W`; `x^l P^m = (lm)! x^l` modulo `M_W`, outside `M_W` for every `m`, not only for infinitely many.
Correct.

Corollary A3. `[W^m] = delta (delta + 1)...(delta + m - 1) [1] = 0` for `m >= 1 - delta`;
`[h W^N] = (e + delta)...(e + delta + N - 1) [h]`. The hypothesis `[h] != 0` is enough; nothing is missing in the
statement (see U6 for how `[h] != 0` has to be certified).

Proposition A4. In the frame `du ^ dv`: `div E = (d_u U + d_v V)/2 = 1`, `E(UV) = UV`; the rotation field has
divergence 0 and kills `UV`. Correct. Checked with explicit functionals (`a4`): standard pair with a symbolic
centre (functional `(e^(d_x d_y) h)(a, b)`), `m <= 3`, 64 identities; double cover of S24 at the centre (4, 3)
with the two residue functionals, `m <= 2`, 126 identities. Use (i): correct (and it is S24's own Proposition
2.2 (e) applied to the charge grading). Use (ii): M5. Use (iii): an analogy, consistent with THM-1645.

`scripts/euler_field.py` re-run: output identical to `runs/euler_field.txt`.

### 4.2 Literature (one search, three abstracts / papers read through a summarising fetch)

- W. Zhao, arXiv:0902.0210, Conjecture 1.3: for any commutative algebra over a field of characteristic zero and
  any commuting set of operators `h(z) + sum c_i d/dz_i`. Section 2 reduces to `d/dz_i - d q/d z_i` and proves
  finite codimension when `A` is a field (Corollary 2.11). Its failure examples have a non-constant leading
  coefficient (Example 2.6) or positive characteristic (Example 2.7); no example over a field of characteristic
  zero with constant leading coefficients.
- van den Essen, Wright, Zhao, arXiv:1004.0521, Example 2.4: the derivation `d/dx - y^2 d/dy`; 1 lies in the
  image and `y` does not, so the image is not a Mathieu subspace, by the 1-property. Derivations only; no twisted
  pair, no `x + x^2 y`.
- van den Essen, Wright, Zhao, arXiv:1008.3962, Conjecture 2.1: the Image Conjecture for `d/dz_i - a_i` with
  `(a_i)` a regular sequence of `A`; counterexample without regularity only in characteristic `p`.
- The search ("Mathieu subspace", "x + x^2 y", image conjecture, commuting operators of order one) returned no
  record of the failure for `W = x + x^2 y`, nor of a failure with 1 outside the image for a pair
  `(d/dx - W_x, d/dy - W_y)` over a field.
Conclusion: not found; this is a search, not a clearance; THM-5147's "no priority is claimed" is the right label.
THM-2801 (canon) refutes SIC(n), `n >= 2`, hence the general conjecture; THM-5147's sentence is accurate (m7).
The public sources of 2026 quoted by S24 were not re-checked.

## 5. Part B in detail

### 5.1 Proposition B1, Corollary B2 (own code, `tt.py`: subset dynamic programme, numba)

- B1, both formulas: all 81,225 ordered pairs of labeled rooted tournaments with at most 4 vertices, `i = 1`;
  3,000 random pairs for each of `i = 2, 3`; 600 random pairs of orders 3 to 5, `i = 1..3`: 0 failures (`b1`).
  With `i = 0` the same formula fails in 5,183 of 14,400 cases, as the note says (parity split needed).
- One-apex data on all 511 words of length at most 8: `H = g(P)`, `W+ = P + 1`, `W- = Q + 1`, `tau = 2^n`,
  `Q + 1 = 2^k (R(p) + 1)`, `n = k + L(p)`, `g(Q) = g(P)` (`b2`, K1). A closed formula for `g` independent of
  the recursion (count of the pairs (set before the apex, set after)) agrees on all words of length at most 13.
- Corollary B2: 400 random chains, depth 2 to 4, words of at most 4 letters, chain lengths 1 to 3, at most 17
  vertices: `H` and `W+` at the last apex equal the entries of the matrix product, 0 failures (K2).

### 5.2 `scripts/chain_semigroup.py`

Re-run: `c1`, `depth 18`, `closure 12`, `closure 14` reproduce the stored output lines character for character
(`rerun_chain_semigroup.out`). Reading of `closure`:
- every product is enumerated. Odd `p` is handled by `M(2p + 1, i) = diag(1, 2) M(p, i)` (same `A'`, `W'`
  doubled): correct. Leading letters 0 are in the exponent: correct.
- stopping rule in `p`: the loop leaves when `g(p) A + 2W > X` and `p > 2X/A + 4`; safe because
  `g(2v) = g(v) + v + (v mod 2) >= v + 1`. The cap `p <= X` is safe because for even `p >= 2`,
  `2 (R(p) + 1) > p + 2`, so a hit needs `p < 2X/3`.
- states with `W > X/2` are dropped: correct (`A' >= A + 2W`). The index of `seen` uses that `A` is odd: correct.
- initial states: `p_0 < 2X + 8` covers every one-apex value (even cores `<= 2X - 2`): correct.
`depth`: complete for the hits of depth 2 and 3; the depth-1 table is incomplete (m10, harmless); the tuple
counts are redundant (M3).

### 5.3 Independent closure (C4)

`b3_independent_closure.py`: no normal form; every binary word of length at most `K + 1` is a block with its own
`(g, P + 1, Q + 1, n)`, every chain length `i >= 1`; all states `(A, W)` with `A <= X`, a state being expanded iff
`A + 2W <= X` (exact: `A' >= A + 2W`). Below `2^12` (16,383 words, no cut): 414 misses; the full list equals the
author's (closure function called directly) and the 80 stored values. Below `2^14` (65,535 words, words sorted by
`Q + 1` and scanned while `2^i W (Q + 1) <= X - A`, an exact cut): 1214 misses, full list identical; its
restriction to `[1, 2^12]` is the list at `2^12`. A third computation from the O/Z/J letter moves gives the same
414 (`d4`, Z4).

### 5.4 Depth 2 and 3 (C2, C3)

Word-level computation, level by level, no normal form (`b5`): below `2^18`, depth at most 2 misses 16,430 odd
numbers (7,250 in `[2^17, 2^18)`, 11.06 per cent), depth at most 3 misses 9,190 (3,446, 5.26 per cent); lists
identical to those of `depth_cover`. Below `2^14`: 1,537 and 1,225. The stripped enumeration (`b4`) gives the same
four numbers at `2^18`. "Depth at most 3 misses 1225 below `2^14`, the closure 1214": confirmed.

### 5.5 Pre-registration integrity

- `git show 7479352f36 --stat`: PREREGISTRATION.md (54 lines), `scripts/euler_field.py`, `runs/euler_field.txt`
  (disclosed as already run) and three stubs. No `chain_semigroup.py`, no part-B output.
- The first 54 lines of PREREGISTRATION.md are identical at `7479352f36`, `8b1398bae9`, `870ad0c5d7`,
  `d67741ce21`, `7319c023bc`: nothing above "Outcomes" was changed.
- `chain_semigroup.py` and `runs/c1_bruteforce.txt`, `closure_12.txt`, `closure_14.txt`, `depth_18.txt` first
  appear in `8b1398bae9` (18:42:48). Reflog of `origin/main`: `7479352f36` pushed at 18:30:18. The numba cache of
  `chain_semigroup.depth_cover` and `closure` in `scripts/__pycache__` (untracked) dates from 18:32:32 and
  18:32:47, after the push.
- Addendum D: the first 104 lines are identical at `870ad0c5d7`, `d67741ce21`, `7319c023bc`; `870ad0c5d7` (pushed
  18:46:52) contains the text only; `direct_joins.py` and its outputs first appear in `d67741ce21` (18:52:01); the
  untracked `direct_closure_16.txt` and `_18.txt` were written at 18:50 and 18:51.
- Limits: git cannot exclude an earlier run from a scratch file (the registration itself discloses one scratch run
  of another family); and for D2 the outcome was in canon before the registration (M1 (iii)).
- Outcome labels: C1, C2 CONFIRMED and C3, C4 FAILED are right. C5: M3. D1, D3 CONFIRMED and D2 FAILED are right
  as arithmetic.

### 5.6 Credits of B1 (README)

| cited | as described? |
|---|---|
| THM-4876 (increment at a vertex) | yes: statement (1), KNOWN there |
| THM-4866 (B) (`H = x_h y_h + x_R y_R - x_e y_e`) | yes; `x_e = o - e` checked (U8) |
| Session263, Theorem RJ2 | yes: the packing polynomial of a return join; at `(t, u) = (2, 1)` it is the direct join law |
| THM-4875, THM-4885, THM-4890 (one-apex `g`, triple `(g(v), e(v), 1 + v)`) | yes |
| THM-4885 (B), (C), (D) | yes; `2^35` is one run, not recomputed (m14) |
| THM-5141, S93, S94, Sessions 275-276 (CONCURRENT) | yes, with m14 |
| **missing** | SEAM.md (1)-(7); THM-4875 (B); THM-4904 and HYP-9494 (O/Z/J words, the 29 misses, depth); THM-4888 (D) 1 (`JZ = OJ`); HYP-9480 (D) (multiplicity by number of root changes); THM-5130 (A); Session197 (C = 1 class) |

Session274 (`session274_rooted_seams_20261111`; there is no `00-navigation/S274.md` in this worktree) proves the
fan law `Q_k(t) = A (t+1)^k + B t^k + C (t-1)^k`, whose value at `t = 1` is the doubling of the excess by an
inserted transitive chain: the same mechanism as the factor `2^i` of B1, not the two-sided formula. THM-4866 (B)
lists "the near-join `h + 2^k B`" as a special case; THM-5126 (A) is the tail law `c_0 + M 2^t` (KNOWN there);
THM-5130 (A) is the 4 x 4 gadget law. None states B1 with a general `Y` after the chain; all of them, with
SEAM.md, contain it as a one-line specialisation.

### 5.7 Status labels

Blends found: M4 (FINITE-EXACT for a KNOWN list); m12 (KNOWN for a heuristic inference); README section 0, part A,
"PROVED" covers the remark M5; THM-5148 "(A), (D) PROVED" should read "KNOWN, re-derived" (M2). No blend of
CONCURRENT or OPEN was found; "Chan-Pak OPEN, no proof" is stated everywhere it should be.

## 6. Addendum D in detail

### 6.1 Law D1

Derivation checked by hand: the cycles through the reversed arc are (path of `X` from `u`, path of `Y` ending at
the apex) with an odd total, and for a one-apex block the paths ending at the apex number `O = 1 + floor(Q/2)`
(odd) and `E = ceil(Q/2)` (even); paths from the new apex inside `Y` number `1 + floor(P/2)` and `ceil(P/2)`;
through `v -> u` they are followed by a subset of the spine (`t_o`, `t_e`). The form `(A, W, d)` follows from
`O + E = Q + 1`, `O - E = [Q even]`, `t_o + t_e = 2^n`, `t_o - t_e = -[n = 0]`.
Brute force (own code): 600 random chains of depth 2 to 4, words of at most 4 letters, chain lengths 0 to 2, at
most 16 vertices, the state `(A, o, e)` compared with the literal tournament at every one of 1,105 stages; and all
3,234 chains of depth 2, 3 with words of length at most 2 and chain lengths 0, 1, 2: 0 failures (`b2`, K3). The
author's `d1` compares only the final `H`; the audit's check of `o`, `e` at each stage is stronger.
`python direct_joins.py d1`, `closure 11`, `closure 13`: outputs identical to the stored ones.

### 6.2 `scripts/direct_joins.py`

- Cores: every word is `0^a . core . 1^b` uniquely, the core being empty or beginning with 1 and ending with 0
  (an even `P`); `P + 1 = 2^b (P_core + 1)`, `Q + 1 = 2^a (Q_core + 1)`, `Q_core = R(P_core)`, `g` unchanged.
  Correct.
- Loop on cores: a core is skipped when `g A + 2e > X` (right: `A' >= g A + 2 e O`, `O >= 1`), and the loop ends
  when moreover `P > 2 (X // A) + 4` (right: `g(P) >= P/2 + 1` on even `P`; the list reaches `2X + 4`).
- Loop on leading letters 0: `A'` is non-decreasing in `a` and unbounded (`o >= 1`): the `break` is right.
- Loop on trailing letters 1: `e'` increases with `b`; `break` at `e' > X/2` is right (`A'' >= A' + 2 e'`).
- States with `o' > X/2 >= e'` are not stored and the value `A' + 2 e'` is recorded: **correct and sufficient.**
  From such a state a block with `E >= 1` gives `A'' >= A' + 2 o' > X`; `E = 0` means a word without letter 0,
  for which `g = 1`, `O = 1`, `A'' = A' + 2 e'` whatever its length; and the state reached then has
  `e'' >= o' > X/2`, so it is dead.
- Initial states: `P < 4X + 8`; unstored ones have `e > X/2` (dead) or `P = X` (`o = X/2 + 1`, `e = X/2`,
  `A + 2e > X`): no value is lost.
- Hash set: open addressing with linear probing and full key comparison; keys `(A << 40) | (o << 20) | e` with
  `o, e <= X/2 <= 2^17`, `A <= 2^18`; load at most 0.25. Correct (m15).

### 6.3 Independent closures (`d2_independent_direct_closure.py`)

No normal form: every binary word of length at most `K + 2` is a block, the law D1 is applied with its own
`(g, P, Q, n)`.
- literal mode (every state with `A <= X` stored and expanded, no pruning at all): `2^8`: 8,815 states, 15 misses;
  **`2^11`: 242,785 states, 28 misses**: 7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299,
  421, 449, 503, 533, 569, 583, 637, 821, 1067, 1339, 1861. Identical to `runs/direct_closure_11.txt` and to the
  author's closure function.
- states expanded iff `A + 2e <= X` (exact: the empty word gives the least successor `A + 2e`): `2^11` (numpy and
  numba): the same 28; `2^13`: 29 misses, the 29th is 2647; `2^16` (2,676,377 states): the same 29, none in
  `(2647, 2^16]`. Full lists identical to the author's closure and to the stored files.
- `2^18` was not re-run (instruction of the coordinator); the stored output (29 misses, 32,860,481 live states)
  agrees with HYP-9494 (D).

### 6.4 The 29 numbers are counts of other tournaments

KNOWN: bof and Royle (MathOverflow 232751; THM-4749): every odd number below 80,557 except 7 and 21; Session270: a
literal witness on at most 16 vertices for every odd number up to `10^6` except 7 and 21; THM-4875 (D) and
THM-4885 (D) to `2^32` and `2^35`. Own witnesses, each recounted on the literal adjacency matrix
(`d3_witnesses_for_the_29.out`): 13 (order 5), 41, 43 (order 6), and the other 24 as joins, with the arc between
the roots reversed, of a rooted tournament with at most 6 vertices with a one-apex block or with a second rooted
tournament with at most 6 vertices (67 and 83 at order 7; order at most 12 in all). So the statement is about
this family only, as the note says.

### 6.5 Name and prior statements of the family

See M1. Names in canon: "O/Z/J words", "root-change words" (THM-4904), "the physical O/Z/J grammar" (Session200),
"the C = 1 class" / uniquely Hamiltonian strong tournaments and their one-seam grammar (Session197). Depth 2 is
the family F1 of THM-4875 (D) (inside G(2,0) of THM-4878). THM-4878 and THM-4880 treat other grammars (spine
words of bounded cutwidth, ribbon words); THM-4885 treats F1 and bank joins. The miss set is recorded in
Session197 C1_COMPLETE.md, THM-4884 (2b), THM-4904 (status), HYP-9494 (D).

### 6.6 By-product: the numbers of addendum E (not audited as statements)

The audit's closure restricted to blocks of word length at most `L` reproduces `runs/finite_blocks_16.txt` and
`runs/finite_blocks_trend_17.txt`: below `2^16`, `L = 0..4`: 32,750, 28,170, 10,078, 1,078, 126 misses; below
`2^17`, `L = 3, 4, 5`: 1,355, 127, 49 misses, largest 130,247, 130,247, 10,391, with the same counts per dyadic
block (`e1_finite_blocks_bonus.out`). HYP-9677 (B) and README B5 need the same prior-art check as M1 (HYP-9480,
HYP-9494, THM-4888) before they are filed; that check was not made here.

## 7. What was recomputed

| item | program | result |
|---|---|---|
| `lambda` kills `M_W`, `l <= 4`; values; residue formula | `a1_lambda_form.py` | 672 + 336 + 108 checks, 0 failures |
| membership certificates over `Q`, `l = 1, 2, 3` | `a2_membership_over_Q.py` | 23 vector fields found and re-expanded; 9 controls without solution |
| `x (x^2 y + x + 1)` and two more potentials | `a3_no_euler_field_example.py` | no constant-divergence Euler field; `W^2`, `W^3`, `W^4`, `x^2 W^n - n! x^2` in `M_W` |
| family `x^(l+1) y + r(x)` | `a5_family_x_l1_y_plus_r.py` | five potentials: U1 on boxes, the moment formula of U2 (175 cases), certificates |
| radial and rotation identities | `a4_radial_identity.py` | 64 + 126 identities; resultants of `e_m`, `e_(m+1)` nonzero |
| `euler_field.py` | re-run | identical to `runs/euler_field.txt` |
| Proposition B1; `W >= A` | `b1_join_law.py` | 87,825 joins; 202,013 rooted tournaments |
| one-apex data, Corollary B2, law D1, redundancy | `b2_one_apex_chains.py` | 511 words; 400 + 600 + 3,234 chains; 9,600 identities |
| `chain_semigroup.py`, `direct_joins.py` | re-run | identical to the stored outputs |
| closure of the pure chains | `b3_independent_closure.py` | 414 below `2^12`, 1214 below `2^14`, full lists identical |
| depth 2, 3 below `2^18`; multiplicities | `b4_depth_and_multiplicity.py`, `b5_depth_word_level.py` | 16,430 / 7,250 / 9,190 / 3,446; 14.3, 43.9, 82.9 |
| `x_e = o - e` | `b6_xe_identity.py` | 202,013 rooted tournaments |
| closure of the direct joins | `d2_independent_direct_closure.py` | 28 below `2^11` (literal), 29 below `2^13` and `2^16` |
| witnesses for the 27 values | `d3_witnesses_for_the_29.py` | 27 of 27, literal recounts |
| identification with O/Z/J | `d4_ozj_identification.py` | 29,524 + 797,161 words; 28 misses; `OJ = JZ`; 414; 9,840 cycle counts |
| addendum E numbers (by-product) | `d2_... blocks` | 8 closures, all equal to the stored ones |

## 8. Not verified

- The run of `direct_joins.py` at `2^18`, and the census of HYP-9494 (D) to `10^9` on which M1 relies for the
  range above `2^16`. (Below `2^16` the audit's own closures give the same list.)
- Session197's classification of the uniquely Hamiltonian strong tournaments (ROOTED.md Theorem 2, C1_COMPLETE.md)
  and the proofs of SEAM.md, THM-4875 (B), THM-4904: quoted, not re-audited. The identification of the families
  was checked directly on words and tournaments.
- The papers of section 4.2 were read through a summarising fetch tool, not line by line; the public sources of
  2026 on the Image Conjecture quoted by S24, and THM-2801's proof, were not re-checked.
- The classification of polynomials with general fibre `C^*` used in m4 and U2 (quoted from memory: H. Saito,
  Osaka J. Math. 14 (1977); the note gives no reference).
- The S24 engine (`stein.py`) and the cited comparison theorem: background, not re-audited; part A no longer needs
  them (U1).
- Whether the coverage data of the pure chains (C2-C4) is recorded elsewhere in the repository. A grep for their
  first misses (7, 13, 21, 29, 35, 41, 43, 49) finds only lists of the one-apex family, which begin in the same
  way but contain 57, a value of the pure chains; THM-4904, THM-4888, HYP-9480 and HYP-9494 do not treat this
  sub-language. A search, not a clearance.
- The number of non-isomorphic tournaments behind the 14.3 chains per odd number.
- Addendum E as a set of statements (section 6.6 reproduces its numbers only); HYP-9677 (B).
- That no run preceded the registrations: git and file times are consistent with the claim and cannot prove it.

## 9. Files of this audit

`AUDIT.md`; `tt.py` (toolkit); `a1`-`a5`, `b1`-`b6`, `d2`-`d4` (`.py` and `.out`);
`b3_misses_closure_12.txt`, `b3_misses_closure_14.txt` (full miss lists of the pure chains);
`rerun_euler_field.out`, `rerun_chain_semigroup.out`, `rerun_direct_joins.out`; `e1_finite_blocks_bonus.out`.
Reproduce: run each script from this folder (`python b3_independent_closure.py 12`, `... 14 cut`;
`python b4_depth_and_multiplicity.py 18`; `python b5_depth_word_level.py 18`;
`python d2_independent_direct_closure.py literal 11`, `numba 16`, `blocks 17 3 4 5`); the longest takes about one
minute.

Sources read for section 4.2: https://arxiv.org/abs/0902.0210 , https://arxiv.org/abs/1004.0521 ,
https://arxiv.org/abs/1008.3962 .
