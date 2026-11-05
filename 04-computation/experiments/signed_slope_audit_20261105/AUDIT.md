# Independent adversarial audit of the opus-2026-11-05-S80 packet (THM-4882 / HYP-9475)

**Auditor:** independent agent (own code only; nothing in `C:/tmp/wt-s80` or the main repository was modified).
**Audit directory:** `C:/Users/Eliott/AppData/Local/Temp/claude/800365894db7d24c/f973fbc2-2adf-4ee3-a3b5-b87fe066eb0a/scratchpad/s80_audit/`
**Worktree state audited:** `C:/tmp/wt-s80`, branch `codex/session-s80`. The audit started at HEAD `f9f1400aee`
(the pre-audit checkpoint). The author's session kept committing while the audit ran: HEAD moved to `1a7bbde9dc`
(14:40:13) and then `d2bab54307`. `git diff f9f1400aee HEAD` on the six audited files touches only the note's
section 7 (lines 190+) and HYP-9475 section (B); sections 0-6 of the note, THM-4882, the pre-registration and the
two scripts are byte-identical across the three HEADs. The revised section 7 / HYP-9475 (B) (two-apex drift
functionals, `two_apex_drift_functionals_20261105.py`) were added after the audit snapshot and are **not audited**.

## VERDICT: PASS WITH CORRECTIONS

The two new mathematical claims, THM-4882 (3)(ii) (both odd: all even integers) and (3)(iii) (even denominator:
never complete, one residue class modulo 4 on one sign side), are **correct and proved**; I could not break them with
small cases, large negative targets, slopes outside the note's list, or an independent physical (subset-DP)
materialisation of the recipe. The barrier (4) is a correct (trivial) corollary of the known monotonicity. The
FINITE-EXACT records reproduce exactly. The corrections concern (a) missed prior art for statement (1) and for the
"boundary term" of the transfer, (b) the typing of the pre-registered outcomes P1/P2 (the literal small-witness
predictions are refuted for most of the larger-denominator slopes, contrary to the note's wording), and (c) several
loose statements (non-negativity of the period matrix, "left eigenvector", slope count, bounds that are far from sharp).

---

## 1. Task 1 - proof check of THM-4882 (3)(ii) (a, b odd: all even integers)

Checked line by line against the note's section 4 and Session 198's packet.

| step | verdict | remark |
|---|---|---|
| A = a(2^k - 1)/b is odd | correct | 2^k - 1 odd, b odd => (2^k-1)/b odd; a odd. |
| B/2 <= A < B - 1, word starts and ends with Z, mixed | correct | A > (B-1)/2 and A integer gives A >= B/2; a/b < 1 gives A < B-1 (not all ones). k >= 2 since b >= 3. |
| transfer D' = D + AU - C_lin d - eps_0 E, U' = BU - Rd, E' = d | correct | derived independently: before letter i, U_i = 2^i U - d sum_{j<i} eps_j 2^(i-1-j), E_0 = E, E_i = eps_(i-1) d (i >= 1); a Z adds U_i - E_i. The term -eps_0 E is the input ending count seen by a first-letter Z. **Prior art:** this exact formula with the term `-E epsilon_0` is already proved in `session198/source-digit-blocks/DIGIT_BLOCK.md`, equation (1) and the displayed derivation `D' = D + UA - E epsilon_0 - d C(A)` (see Major 1). |
| F' = F - K_lin d + b eps_0 E | correct | uses bA = a(B-1). |
| K_cyc = K_lin - b on strong inputs (E = d = 1, eps_0 = 1); K_lin after an O (E = 0) | correct | verified on random inputs and on 77 physical instances (section 7 below). |
| positivity of K_cyc via Session 198's run argument | correct, and the cyclic indexing is genuinely needed and genuinely available | Session 198's (5) is Delta = sum_i eps_i (u_i - eps_(i-1)) with eps_(-1) = eps_(k-1); its run decomposition treats the cyclic word, so the run through the boundary k-1 -> 0 is one run. With eps_0 = 1 the identity is Delta_cyc = AR/(B-1) - C_cyc = K_cyc/b (the cyclic sum of eps_i eps_(i-1) is C_lin's adjacent sum plus eps_(k-1) eps_0). Session 198's own audit (`audit/rational-replay/AUDIT.md`) already notes "since the first bit is zero, the cyclic ending-count convention adds no missing boundary term", i.e. the only difference for a first-letter Z is this boundary term. My `recipe_test.py` recomputes Delta_cyc with exact Fractions and the run decomposition for 13 slopes: Delta_cyc = K_cyc/b every time, all runs positive (e.g. 7/11: runs 6/11, 8/11, 4/11, K_cyc = 18). |
| evenness of K_cyc | correct | F = aU - bD is even on strong orbit points (U, D, a, b odd) and W maps strong to strong (last letter Z, d = 1). K_lin = K_cyc + b is odd. |
| residue equation 2aU_0 - bD_0 = r + b (mod K_cyc) with odd D_0, U_0 | correct and solvable | D_0 = 1+2x, U_0 = 1+2y gives 4ay - 2bx = r + 2b - 2a (mod K_cyc), halved to 2ay - bx = (r+2b-2a)/2 (mod K_cyc/2); gcd(2a, b) = 1. The case K_cyc = 2 (modulus 1 after halving) is fine. |
| use of STRONG_SOURCE_INTERFACE | correct as a black box | Its literal sources are V(L) = T(L)^op with profile (g(L) + L - (L mod 2), 1, 2L - 1, 1) = (g(2L-1), 1, 2L-1, 1): exactly the odd orbit point n = 2L - 1, i.e. the one-apex tournament rooted at the apex whose chain word is the bits of L - 1 followed by a 0 (the vertex a). "Arbitrarily large D" follows from JOINT_MODULAR_GADGET's arbitrarily large L and D >= L. It requires M even and odd residues: K_cyc is even and >= 2. Empirical support for the moduli in play: for every even M <= 50 every pair of odd residues (g(n) mod M, n mod M) occurs among odd n < 2^22, minimum multiplicity 3,065 at M = 50 (`residue_and_misses.py`). |
| h = (S - b - r)/K_cyc >= 1, final strongness | correct | S = 2aU - bD >= (2a-b)D >= r + b + K_cyc gives h >= 1; after the first W (drop K_lin) the input is strong, later drops are K_cyc; the final letter is Z, so E = d = 1, s = 1, transitive root deletion. |

**Independent reproduction (own code, `recipe_test.py`, 29 s):** 13 both-odd slopes - the note's 3/5, 5/7, 5/9, 7/9,
7/11, 9/11, 7/13, 9/13 and the new 11/13, 9/17, 13/15, 11/19, 15/17 - with 23 targets each (0, +-2, +-4, +-6, +-100,
-1000, -10000, 20000, -30000 and 10 random even targets in [-20000, 20000]); every target hit exactly by
seed -> O -> W^h with a strong final point, and the final point checked to lie on the orbit (D_final = g(U_final)
recomputed from the big integer U_final, an independent check of the whole transfer). Physical check: 77 small clock
outputs (orders 8-19, seeds n_0 in {3,...,15}, h = 1, 2, 3, five slopes) were materialised as tournaments and their
(H, E, U, d, s) counted by an independent subset DP (`hamcount.cpp`): 77/77 agree with the O/Z profile, including
the first-period drop K_lin and the later drops K_cyc.

## 2. Task 2 - proof check of (3)(iii) (b even)

- **Lemma 1** (g odd, g(n) <= n, g(n) >= (n+3)/2 for odd n >= 3): correct; the induction as written is sound (the
  sentence "g(2n-1) <= 2n - (n mod 2), and for even n this is at most 2n-1 because g(n) = g(n/2) <= n/2 <= n-1" is
  garbled - for even n the bound 2n - (n mod 2) equals 2n, and the intended argument is g(2n-1) = g(n) + n <= n/2 + n
  <= 2n - 1 - but the conclusion is right). Verified numerically to n = 2^24 (`mod4_check.py`); the lower bound is
  attained exactly at n = 2^j - 1 (23 cases below 2^24).
- **Lemma 2**: correct (direct substitution).
- **Case 4a <= 3b**: odd m >= 3 gives F <= (2a - 3b/2) m - a - b/2 < 0 (correct; the maximum of the odd-m branch
  observed numerically is exactly -a - b/2 for 3/4 and below it otherwise); even m = 2m' gives F = (4a - 2b) m' - a - b g(m')
  with 4a - 2b = 0 (mod 4) and b g(m') = b (mod 4) (b even, g odd), so F = -(a+b) (mod 4). Correct. The positive odd
  integers in the other class are missed: infinitely many.
- **Case 4a > 3b**: even m: F >= (4a-3b) m' - a, non-negative for m' >= a/(4a-3b), and > -a - ab/(4a-3b) = -N_0 for
  smaller m'. Odd m: F = (2a-b) m - a + b - b g(m); 4 | b gives (2a-b)m = 2 (mod 4) and F = 2 - a (mod 4); b = 2 (mod 4)
  gives 2a - b = 0 (mod 4) and F = -a (mod 4). Correct in both subcases; the elements below -N_0 lie in one class, the
  other class below -N_0 is missed: infinitely many.
- "Misses infinitely many odd integers" follows in both cases. The parity statement (all elements odd) is correct.

**Numerical search for counterexamples (`mod4_check.py`, n < 2^24, 1 m 45 s):** all 35 coprime slopes with even
b <= 26 inside (1/2, 1) plus 17/18, 19/20, 11/20, 13/20, 17/24, 23/24, 19/24, 13/24, 25/26, 15/26, 31/32, 17/32:
0 violations of the positive-side confinement (4a <= 3b), 0 violations of the confinement below -N_0 (4a > 3b),
odd-m branch negative in every 4a <= 3b case. The note's seven + eight slopes are reproduced (classes 1/3 as listed).
Observation: N_0 is far from sharp - the most negative even-m value is always -3(b - a) (attained at n = 3), e.g. -3 for
5/6 against N_0 = 20, -9 for 11/14 against N_0 = 88 (Upgrade 2).

## 3. Task 3 - the dictionary (1) and the barrier (4)

**Dictionary, own brute force (`dict_bruteforce.py`):** all 254 chain words of length 1..7 (both leading bits),
apex x beats v_i iff b_i = 1, chain v_1 -> ... -> v_m. Subset DP for H, E (paths ending at x), s (paths starting at
x), d = H(T - x), and U = sum over X subset of chain of H(T[X]) * (paths of T[V - X] starting at x) with the empty
X counted once; permutation enumeration as a second count for orders <= 7. Result: 0 failures of
(D, E, U, d) = (g(1 + bin(b)), 1 - b_m, 1 + bin(b), 1), H odd, s = b_1; all 252 single-bit appends satisfy
bit 1 = O (D, 0, 2U, d) and bit 0 = Z (D - E + U, d, 2U - d, d); the seed word "1" has profile (1, 0, 2, 1) = Session
198's t. The probe's own check (g(n) = H(word of value n-1) for n <= 4096, 0 mismatches) reproduces. The dictionary
is therefore correct - but it is not new (Major 1).

**Barrier (4):** correct. In G(k, 0) (THM-4878: k <= 2 apexes x_1 -> x_2, each beating an arbitrary subset of a
transitive chain, no reversed chain arcs) deleting the chain vertices of one period of W from T(p W^(h+1) q) leaves the
induced subtournament whose vertices are the apexes with their fixed mutual arc and the remaining chain vertices in
chain order with their letters' apex incidences, i.e. T(p W^h q); H(T - v) <= H(T) by Redei insertion (THM-1370),
hence H(p W^(h+1) q) >= H(p W^h q). In G(k, w) the induced subtournament equals T(p W^h q) only when no reversed arc
joins a deleted period vertex to a surviving vertex ("every slot of W closes inside the same period", as hedged); a
slot of p that crosses the whole of W into q survives deletion but changes its letter-distance span, so the hedge
should say "slots whose tail or head lies in W". The rhetorical consequence ("the dominant eigenvalue survives in
every non-negative readout") is interpretation, not part of the proved statement.

## 4. Task 4 - status typing and prior art

| claim | packet label | audit |
|---|---|---|
| (1) dictionary O/Z = G(1,0), g = one-apex count | PROVED, elementary ("merge", "new") | correct but **KNOWN**: see Major 1. |
| (1) clock = left eigenvector with Jordan drift | PROVED | loose: phi = aU - bD is a generalised eigenvector (phi M_W = phi - K_lin psi + b eps_0 chi_E; on the slice E = d this is phi - K_cyc psi), and M_W is **not** non-negative in the (D, E, U, d) coordinates (rows E: (-eps_0, 0, 0, 0), d: (-C_lin, 1, -R, 1)); the Jordan-block reading of an additive drift H(u B^k) = H(u) + kc is already in the ribbon lane (`05-knowledge/results/ribbon_tournaments_bounded_span_transfer_20261104.md` lines 248-252; `ribbon_transfer_audit_20261104.md` lines 48-56) and in Session 194 two-apex (Jordan justification, section "Self-contained common-flag and Jordan justification"). Minor 1. |
| (2) lemmas | PROVED | correct. |
| (3)(i) | KNOWN (Session 198) | correct citation; Session 198's own audit verdict PASS. |
| (3)(ii) | PROVED, new | **correct and new.** No packet on the branch treats both-odd slopes: `git grep` of "both odd / even integers / even denominator / odd denominator" over session194-199, 01-canon, 05-knowledge, 04-computation returns no tournament-lane hit other than Session 198's even-a/odd-b hypothesis (`session198/rational-response-coverage/PROOF.md:4`). |
| (3)(iii) | PROVED, new | **correct and new.** The only prior modulo-4 statement is Session 198's 3/4 control (`PROOF.md:183-185`, "every positive value 3n-4g(n) is 1 modulo 4"), which is the single slope a = 3, b = 4 of the 4a <= 3b case; no other packet contains the even-denominator obstruction (`git grep -i -E "modulo ?4|mod ?4"` over session197-199: 138 lines, all listed in `grep/mod4_sessions.txt`, none about signed slopes except those). |
| (4) barrier | PROVED | correct; elementary corollary of monotonicity that is used verbatim in earlier packets (Minor 3): `session195/full-binary-return/PROOF.md:43` "Counts h never decrease, because every operation preserves the old tournament as an induced subtournament"; `session197/minimum-positive-clock/PROOF.md:57` (Induced-count monotonicity); `session197/opposed-binary/PROOF.md:34`; `session196/spectral-moves/PROOF.md:206`; `session197/capped-ports/PROOF.md:85`. Session 194 two-apex (`session194/two-apex/PROOF.md`, Theorem 1) proves the stronger fact that every mixed two-apex block of length <= 6 carries a non-constant polynomial counter in the repetition number. |
| FINITE-EXACT (probe, recipe) | FINITE-EXACT | reproduced exactly (section 6). |
| HYP-9475 P1-P6 outcomes | P1 supported, P2 supported/PROVED, P3-P6 as stated | P3, P4 (refuted), P5, P6 correctly typed; **P1 and P2 mis-typed** (Major 2). |

**Prior art for the identification of O/Z with the one-apex automaton / HYP-9468's g (Task 4 item):**
- `session197/path-polynomials/ONE_CORE_BINARY_INVERSE.md:13-16`: "Call them O and Z here, to keep the binary
  encoding visible. O appends a universal sink. Z appends a vertex beating only the retained root a. Their exact
  state laws are O(h,e,R,d) = (h,0,2R,d), Z(h,e,R,d) = (h+R-e,d,2R-d,d)"; line 20: "The scalar binary recurrence used
  below is also inherited. It is exactly f(q)=g(q+1) in the reviewed modular-gadget theorem"; line 59-63: the scalar
  output of any O/Z suffix reads the bits of t "from least significant to most, with Z=1 and O=0".
- `session198/source-digit-blocks/DIGIT_BLOCK.md:7-33`: "Let g be the already proved one-apex scalar function";
  equation (1) `D' = D + (U-d)A + d(g(A+1)-1) + (d-E)p` with p = A mod 2 = eps_0, derived from
  `D' = D + UA - E epsilon_0 - d C(A)` - the transfer with the first-letter boundary term for arbitrary words.
- `session198/rational-response-coverage/PROOF.md` section 3 "An equivalent arithmetic image of the one-apex sequence",
  (7) (D,E,U,d) = (g(n), n mod 2, n, 1).
- `session198/rational-response-coverage/dependencies/JOINT_MODULAR_GADGET.md:12-14` f(n) = g(n+1) with the recurrence
  (1) and the gadget "root beats x_i exactly when bit i of L-1 is zero".
- `session199/hidden-structure/PROOF.md:78-85` (concurrent): "every integer U>=2 has one O/Z word from t. Its D is
  the function g(U)".
- HYP-9468 (Setting): the rooted triple at the apex (g(v), e(v), 1 + v) with x -> i iff the i-th bit is 1 (credited
  there to session179 spectrum/REPORT.md (1) and session124's f); R_x = sum_X h(A[X]) s_a(A - X) is literally U.

**Pre-registration (Task 5):** `git log` in the worktree: `0e5ba011be` 2026-11-05 14:22:55 -0700 adds
`signed_slope_prereg_20261105.md` (P1-P5, 25 lines) and `signed_slope_probe_20261105.cpp`, nothing else;
`c6d7e35dc4` 14:26:54 adds `signed_clock_recipe_20261105.py` and appends the two P6 lines to the pre-registration
(diff shown by `git show c6d7e35dc4 -- 04-computation/experiments/signed_slope_prereg_20261105.md`: only `+` lines
26-27); `f9f1400aee` 14:33:00 adds the note, THM-4882, HYP-9475 and index rows. The probe and recipe sources are
byte-identical between their first commit and HEAD (`git diff 0e5ba011be HEAD -- ...cpp` and
`git diff c6d7e35dc4 HEAD -- ...py` are empty). No `.out` file was committed, so the recorded outcomes were checked
by re-running both scripts (section 6): every number in the note's section 5 matches. The probe results quoted in the
c6d7e35dc4 commit message (dictionary exact; P3 confirmed; P4 refuted) are consistent with the re-run. Timeline is
consistent with "P1-P5 before the probe, P6 before the recipe run" (the probe takes 1.4 s, the recipe 6.4 s).

---

## 5. Majors

**Major 1 - statement (1) and the "boundary term" are prior art (status typing).** THM-4882 (1) is presented as the
session's "merge" ("PROVED, elementary", "new" in the index row and ledger). In substance it is already on the branch:
Session 197 `path-polynomials/ONE_CORE_BINARY_INVERSE.md` defines O/Z as the sink / root-beating appends "to keep the
binary encoding visible" and identifies the scalar with f(q) = g(q+1); Session 198 `source-digit-blocks/DIGIT_BLOCK.md`
calls g "the already proved one-apex scalar function" and proves the general block transfer **including the first-letter
term -E epsilon_0** that the note's section 4 introduces as its own addition ("which Session 198 did not need because
its first letter is O"); Session 198's rational-response-coverage (7) and its audit, and the concurrent Session 199
hidden-structure (3), state D = g(U) on the O/Z orbit; HYP-9468 / session179 give the apex triple (g(v), e(v), 1+v).
Correction: type (1) as KNOWN/CITED (Session 197 ONE_CORE_BINARY_INVERSE, Session 198 DIGIT_BLOCK (1) and (7),
HYP-9468/session179), with S80's contribution limited to (a) the explicit bit-order dictionary in THM-4878's G(1,0)
vocabulary (b_1 = most significant bit, value n - 1) and (b) the automaton cross-check to n = 4096; cite DIGIT_BLOCK (1)
for the transfer with the boundary term in section 4; add a MISTAKES ledger entry (missed prior art in the same
session198 tree the owner handed over).

**Major 2 - P1 and P2 outcomes are mis-typed and the note's explanation of the misses is false.** P1 and P2 predict
"every ... target in [-999, 999] is attained with n < 2^24 except possibly a few of the largest". The probe (re-run)
misses, with n < 2^24: 6/11: 157 targets including 7 and -11; 8/11: 125 including -5 and 15; 8/13: 144 including -1;
10/11: 138 including 3 and -15; 10/13: 132 including 13; 4/7: 42 including 131; 8/9: 69 including 79; and in class
(ii) 7/11: 132 including 0; 9/11: 132 including 2; 9/13: 134 including 4; 7/13: 151 including -8 and 10; 5/9: 73
including 74 (`residue_and_misses.py`, part (b)). So the literal small-witness clause of P1 is clearly REFUTED for 7 of the 10
slopes (8/9, 4/7, 6/11, 8/11, 8/13, 10/11, 10/13: misses at |t| <= 131) and that of P2 for 6 of the 8 slopes (5/9, 7/9,
7/11, 9/11, 7/13, 9/13: misses at |t| <= 78, including the targets 0, 2 and 4); it holds cleanly only for 2/3 (0 misses)
and borderline for 4/5, 6/7, 3/5, 5/7 (14-19 misses, all at |t| >= 177, which are not "the largest" either). The note's
sentence "the misses are large targets whose witnesses need n beyond 2^24" is wrong (the witnesses do need n beyond
2^24, but the targets are small). Correction: record P1 as "small-witness clause REFUTED for the larger denominators; structural clause KNOWN
(Session 198)", P2 as "parity clause CONFIRMED, small-witness clause REFUTED for the larger denominators, exact
coverage PROVED"; fix the sentence in section 5 and the HYP-9475 title/status ("P2 CONFIRMED then PROVED" vs
"P2 supported and then PROVED" are also mutually inconsistent). This does not touch the theorem: (3)(ii) is proved
by the recipe, whose witnesses are large by design (h = O(|r|) periods, U with ~kh bits).

## 6. Minors

1. **"Non-negative matrix M_W" and "left eigenvector"** (note section 2 and THM-4882 (1) title): in the profile
   coordinates the letter matrices are not non-negative (Z has the entries -1 on E -> D and d -> U) and M_W has the
   negative entries -eps_0, -C_lin, -R; the clock functional phi is a generalised eigenvector (Jordan chain (phi, psi) at
   the eigenvalue 1), and the drift identity phi M_W = phi - K psi holds only on the slice E = d (K = K_cyc) - on the
   input after an O it is phi - K_lin psi. The non-negativity holds in THM-4878's cut-state basis, not in (D, E, U, d).
   Cite the ribbon lane and Session 194 for the Jordan-block reading of additive drifts.
2. Section 2: "Session 198's orbit of the seed t is exactly G(1,0)" - it is the sub-language of G(1,0) words with
   leading bit 1 (s = 1); words with a leading 0 have the same H and U but s = 0 and are not sources.
3. Barrier (4): cite the earlier verbatim uses of the monotonicity in the lane (session195 full-binary-return line 43,
   session197 minimum-positive-clock line 57, opposed-binary line 34, capped-ports line 85, session196 spectral-moves
   line 206) and Session 194 two-apex Theorem 1 (non-constant counter for mixed blocks), and label (4) "PROVED
   (elementary corollary of KNOWN monotonicity)".
4. "36 slopes" (THM-4882 Verification, note section 5): the lists contain 37 slopes (10 + 8 + 7 + 8 + 4).
5. "all 2^22 positive values in the single class" (section 5): the positive-value counts are 4,194,282 (3/4),
   2,097,131 (5/8), ... - not 2^22 - and there are 2^23 odd n below 2^24; write "all positive values".
6. Lemma 1's proof sentence for the upper bound is garbled (see section 2 above); the conclusion stands.
7. Section 4, proof of 3: the identity Delta_cyc = AR/(B-1) - C_cyc (the cyclic adjacent sum includes eps_(k-1) eps_0)
   is used but not written; one line would make the appeal to Session 198's (5) self-contained.
8. THM-4882 "Verification": the probe prints residues of the **positive** values only; the negative-side residues
   (classes below -200) come from part (3) of the recipe script. Also the threshold -200 in the script is unrelated to
   the theorem's N_0 (it happens to exceed every N_0 of the listed slopes: max N_0 = 88 for 11/14).
9. The recipe script's inline comment on `target_S` is garbled ("S - target_S >= K_cyc... S - b eps_0 - r >= K_cyc");
   the code is right (it tests S - r - b eps_0 = 0 mod K_cyc and >= K_cyc).
10. The probe prints "least n" only for odd targets, so the class-(ii) lines have empty lists; cosmetic.
11. The worktree advanced during the audit (`1a7bbde9dc`, `d2bab54307`); the revised section 7 and HYP-9475 (B)
    (two-apex drift functionals, "FINITE-EXACT data; naive form withdrawn") are outside this audit and should be covered
    by a separate check before the hypothesis is cited.
12. In G(k, w) the hedge "when every slot of W closes inside the same period" should read "when no reversed arc has
    exactly one end in the deleted period" (a slot of p whose head lies in q survives but its span in letters changes).

## 7. Upgrades

1. **(3,4) is confined on both sides.** At 4a = 3b (only (3,4) in lowest terms) the even-m branch is
   F = 4m' - 3 - 4g(m') >= -3, so every element below -3 comes from the odd-m branch and is 2 - a = 3 (mod 4); the
   positive side is 1 (mod 4) as proved. Numerically: negative values below -200 are all 3 mod 4 for 3/4 (recipe part (3)
   and `mod4_check.py`).
2. **Sharp negative threshold for 4a > 3b.** For m' >= 1, F(even m) >= (4a - 3b) m' - a >= 3a - 3b, with equality at
   n = 3 (m' = 1). Hence every element of S(a, b) below -3(b - a) lies in the single class (2 - a or -a mod 4), and the
   bound is attained: -3(b - a) is in S. This replaces N_0 = a + ab/(4a - 3b) (e.g. 20 -> 3 for 5/6, 88 -> 9 for
   11/14, 204 -> 15 for 17/22); verified for all 35 slopes to 2^24 (the most negative even-m value is -3(b-a) every time).
3. **(3)(ii) verified on five further slopes** not in the packet (11/13, 9/17, 13/15, 11/19, 15/17; K_cyc = 18, 2, 8,
   50, 2) and physically on 77 materialised clock outputs (orders 8-19).
4. The modular black box was checked for every even modulus <= 50: all odd residue pairs (g(n), n) mod M occur with
   large multiplicity among odd n < 2^22 (supporting evidence only; the theorem is Session 198's).

## 8. Reproductions (exact commands and outputs)

All commands run from `C:/Users/Eliott/AppData/Local/Temp/claude/800365894db7d24c/f973fbc2-2adf-4ee3-a3b5-b87fe066eb0a/scratchpad/s80_audit/`
(Git Bash; python3 3.13 with numpy; g++).

### 8.1 Git evidence (pre-registration)

```
$ cd C:/tmp/wt-s80 && git log --format='%h %ad %an %s' --date=iso -- 04-computation/experiments/signed_slope_prereg_20261105.md \
    04-computation/experiments/signed_slope_probe_20261105.cpp 04-computation/experiments/signed_clock_recipe_20261105.py \
    05-knowledge/results/signed_responses_slope_classification_clock_barrier_20261105.md
f9f1400aee 2026-11-05 14:33:00 -0700 ... checkpoint (pre-audit) - THM-4882 ...
c6d7e35dc4 2026-11-05 14:26:54 -0700 ... probe results (...); pre-register P6 (...) before the recipe run
0e5ba011be 2026-11-05 14:22:55 -0700 ... pre-register P1-P5 ... before running the probe; probe source
$ git show --stat 0e5ba011be   -> signed_slope_prereg_20261105.md (+25), signed_slope_probe_20261105.cpp (+41)
$ git show --stat c6d7e35dc4   -> signed_clock_recipe_20261105.py (+92), signed_slope_prereg_20261105.md (+2: the P6 lines)
$ git diff 0e5ba011be HEAD -- 04-computation/experiments/signed_slope_probe_20261105.cpp     (empty)
$ git diff c6d7e35dc4 HEAD -- 04-computation/experiments/signed_clock_recipe_20261105.py     (empty)
$ git log --oneline -3   (at the end of the audit)
d2bab54307 opus-2026-11-05-S80: block-12 constant terms ...
8aa0d1ff18 Session199: prove four-block source compiler ...
1a7bbde9dc opus-2026-11-05-S80: two-apex drift functionals ...
$ git diff f9f1400aee HEAD -- <six audited files> | grep -E '^(diff|@@)'
diff ... HYP-9475-...md   @@ -19,14 +19,26 @@      (section (B) only)
diff ... signed_responses_slope_classification_clock_barrier_20261105.md   @@ -190,10 +190,22 @@   (section 7 only)
```

### 8.2 Probe re-run (author's script, 37 slopes)

```
$ g++ -O2 -o probe.exe C:/tmp/wt-s80/04-computation/experiments/signed_slope_probe_20261105.cpp
$ time ./probe.exe 24 999 2 3 4 5 6 7 8 9 4 7 6 11 8 11 8 13 10 11 10 13 3 5 5 7 5 9 7 9 7 11 9 11 7 13 9 13 \
      3 4 5 8 7 10 7 12 9 14 11 16 9 16 7 8 5 6 9 10 11 12 13 14 15 16 13 16 11 14 1 3 2 7 4 3 5 4 > probe_rerun.out
real 0m1.397s
dictionary check g(n) = H(one-apex word of value n-1) for n <= 4096: 0 mismatches
attained (class i): 2/3 1000, 4/5 986, 6/7 981, 8/9 931, 4/7 958, 6/11 843, 8/11 875, 8/13 856, 10/11 862, 10/13 868 of 1000
attained (class ii): 3/5 984, 5/7 980, 5/9 926, 7/9 935, 7/11 867, 9/11 867, 7/13 848, 9/13 865 of 999
positive residues mod 4, a/b <= 3/4: 3/4 {1}, 5/8 {3}, 7/10 {3}, 7/12 {1}, 9/14 {1}, 11/16 {1}, 9/16 {3}
positive residues mod 4, a/b > 3/4: both classes for 7/8, 5/6, 9/10, 11/12, 13/14, 15/16, 13/16, 11/14
outside: 1/3 and 2/7 no positive values; 4/3 and 5/4 no nonpositive values
```
(full output in `probe_rerun.out`; every number equals the note's section 5.)

### 8.3 Recipe re-run (author's script)

```
$ time python3 C:/tmp/wt-s80/04-computation/experiments/signed_clock_recipe_20261105.py > recipe_rerun.out
real 0m6.376s
slope 2/3: k=2 B=4 A=2 W=OZ R=1 K_lin=2 K_cyc=2 | drops verified: True | recipe: 40/40 ... (max repetitions 4966)
slope 4/5: ... W=OOZZ ... K_lin=2 K_cyc=2 | True | 40/40 (4731)      slope 6/7: W=OZZ K_lin=4 K_cyc=4 | 40/40 (2456)
slope 4/7: W=OOZ K_lin=4 K_cyc=4 | 40/40 (2457)                       slope 8/11: W=OOOZOZZZOZ K_lin=18 K_cyc=18 | 40/40 (509)
slope 3/5: W=ZOOZ K_lin=7 K_cyc=2 | 40/40 (4654)                      slope 5/7: W=ZOZ K_lin=11 K_cyc=4 | 40/40 (2436)
slope 5/9: W=ZZOOOZ K_lin=11 K_cyc=2 | 40/40 (4789)                   slope 7/9: W=ZOOOZZ K_lin=11 K_cyc=2 | 40/40 (4734)
slope 7/11: W=ZZOZOOOZOZ K_lin=29 K_cyc=18 | 40/40 (565)               slope 9/13: W=ZZOOZOOOZZOZ K_lin=31 K_cyc=18 | 40/40 (541)
b even: 3/4 positive [1] ... below -200: [3]; 5/8 [3]; 7/10 [3]; 7/12 [1]; 9/14 [1]; 11/16 [1]; 9/16 [3];
        7/8 below -200 [3]; 5/6 [3]; 9/10 [3]; 11/12 [3]; 13/16 [1]; 11/14 [1]; 13/14 [3]; 15/16 [3]
```
(matches the note's constants and section 5; full output in `recipe_rerun.out`.)

### 8.4 Own code

```
$ python3 dict_bruteforce.py 7
words of length 1..7: 254 checked (126 also by permutation enumeration); dictionary failures: 0
transitions checked: 252; O/Z failures: 0
seed word (1): (1, 0, 2, 1, 1)  Session 198 seed t = (1,0,2,1)

$ time python3 recipe_test.py > recipe_test.out            (29 s; 13 both-odd slopes x 23 targets, all hit exactly)
slope 3/5: k=4 A=9 W=ZOOZ R=9 C_lin=4 C_cyc=5 K_lin=7 K_cyc=2 runs=['2/5'] | 23 targets hit exactly; (target, seed n0, h, bits of final U): [(0, 7, 6, 28), (2, 7, 5, 24), (-2, 5, 1, 8), ...]
slope 7/11: ... C_lin=530 C_cyc=531 K_lin=29 K_cyc=18 runs=['6/11', '8/11', '4/11'] | 23 targets hit exactly ...
slope 11/19: k=18 A=151767 W=ZZZOZOZZOOOOZOZOOZ ... K_lin=69 K_cyc=50 runs=[118/171, 16/57, 112/171, 106/171, 22/57] | 23 targets hit exactly ...
77 small words written to small_words.txt

$ g++ -O2 -o hamcount.exe hamcount.cpp && ./hamcount.exe < small_words.txt > hamcount.out
slope 3/5 seed n0=3 h=1 order 8: DP (H,E,U,d,s) = (53,1,87,1,1) recorded (53,1,87,1) response -4 expected -4 -> OK
...  77 lines, 77 verified physically, 0 skipped;  OK lines: 77  MISMATCH lines: 0

$ time python3 mod4_check.py 24 > mod4_check.out          (1 m 45 s)
Lemma 1 verified for n <= 2^24: g odd, g <= n, g(n) >= (n+3)/2 for odd n >= 3; equality at n = [3, 7, 15, 31, ...] (23 cases, all 2^j - 1)
3/4 (<= 3/4): positive classes mod 4 [1] predicted {1} OK; odd-m branch max -5 ...
5/6 (> 3/4): N_0 = 20; classes below -N_0 [3] predicted {3} OK; most negative even-m value -3 ...
... (35 slopes incl. 5/8, 9/10, 15/16, 17/18, 11/20, 13/20, 19/20, 17/24) ...
slopes checked: 35; violations: 0

$ python3 residue_and_misses.py > residue_and_misses.out
(a) M=2..50 even: missing odd residue pairs (g(n) mod M, n mod M) among odd n < 2^22: 0 for every M; min multiplicity 3,065 at M=50
(b) smallest-|t| misses with n < 2^24: 6/11: [7, -11, 27, ...] (157); 8/11: [-5, 15, ...]; 8/13: [-1, 17, ...]; 10/11: [3, -15, ...];
    7/11: [0, -18, 20, ...] (132); 9/11: [2, -16, ...]; 7/13: [-8, 10, ...]; 9/13: [4, 22, ...]; 2/3, 3/5, 4/5, 5/7, 6/7: only |t| >= 177
```

### 8.5 Prior-art listings

`git grep -n` over session194-199, 05-knowledge, 01-canon, 04-computation, 06-writeups, 07-reflections, 00-navigation,
saved untruncated under `grep/`: `classification.txt` (795 lines; 10 in session19x, none about signed slopes except
`session198/rational-response-coverage/PROOF.md:4`), `mod4_sessions.txt` (138), `signed.txt` (204), `barrier.txt`
(605), `oz.txt` (1421), `eigen.txt` (11,926). The decisive hits are quoted in section 4 above.

## 9. Files in the audit directory

`AUDIT.md` (this report), `dict_bruteforce.py`, `recipe_test.py` (+ `recipe_test.out`, `small_words.txt`),
`hamcount.cpp` (+ `hamcount.exe`, `hamcount.out`), `mod4_check.py` (+ `mod4_check.out`), `residue_and_misses.py`
(+ `.out`), `probe.exe`, `probe_rerun.out`, `recipe_rerun.out`, `grep/*.txt`.
