# Independent audit of `collatz_posets_dags_20260927_spine_descent_tree.md` (candidate THM-4514)

**Auditor:** adversarial subagent (Claude), 2026-09-27, blind re-derivation first, then
comparison with the note's proofs; own script
`04-computation/experiments/collatz_posets_dags_20260927_audit.py`
(sha256 `22bcf7c4559ef82ba4fc8bfaedd6c31f1f889e70ee62c9dbab08d62522859c39`), output
`05-knowledge/results/collatz_posets_dags_20260927_audit.out`
(sha256 `f2f8f4b2ce6d46cc6f0037637f89bd3e0498caac84c26409903ae0ad1c85cd9e`), 84 checks,
83 pass; the one deliberate failure is item N4 below. Written without importing or
copying the session script; exact integer arithmetic throughout.

## Verdict

**SOUND WITH CORRECTIONS.** Propositions 1(a),(b), 2, Theorem 3, Proposition 4(a),(c),(d)
and section 5 are correct as mathematics and every printed number reproduces. Four
statements need repair before THM-4514 can be promoted: Proposition 1(c)'s "iff it
diverges" is false in the note's own `T_b` setting and unproven for `3x+1` (C1, also in
THM-4514's title and status); Proposition 4(b) silently assumes the orbit reaches `1`
(C3); the lower end of the `c_D` bracket is the partial sum rounded up (C6, the interval
is true but not by the printed derivation); and the "Collatz iff ... above its threshold"
reading and ranked step (3) overclaim what `N(w) < 2^|w|` buys (C11, C12). The rest are
attribution, path and wording fixes. No conclusion of section 6 is affected.

## 1. Claims checked

### Section 2 (Proposition 1)

* **1(a)** minimal = strict upper records, maximal = strict future minima, dimension `<= 2`.
  **CONFIRMED** (by definition; intersection of two linear orders). Verified from the bare
  definition on the eight record orbits (script G).
* **1(b)** forest order; roots = strict height lower records; infinitely many descendants
  iff strict coefficient future minimum; locally finite when `h -> infinity`; König;
  infinite branches are tails of the spine. **CONFIRMED.** Own derivation: parent(j) is
  the nearest `i < j` with `h_i < h_j`; children of `j` are the running minima of the
  walk on `(j, infinity)` with height above `h_j`; verified against the bare definition
  (transitivity, ancestor chains, roots, children) on the orbits of 27, 703 and 250 steps
  of the `5x+1` orbit of 7. Cosmetic: the chain step "`h > h_(i') > h_i` on `(i', j]`,
  whence `i ⊑ i'`" is a non sequitur as written (`i ⊑ i'` needs `h > h_i` on `(i, i']`,
  which holds because `(i, i'] ⊆ (i, j]`); the conclusion is right. Also: when
  `h -> infinity` the unique infinite branch is the *whole* spine from its first element
  (the global argmin, which is a root), not merely a tail.
* **1(c)** "a strict value future minimum exists iff the orbit diverges, and then there
  are infinitely many". **CORRECTED (C1).** The proof sentence *"If the orbit is
  eventually periodic its global minimum is attained infinitely often, so no time is a
  strict value future minimum"* is false: an eventually periodic orbit can enter its
  cycle from below the cycle minimum. In the note's own control map, the `5x+1` orbit of
  `5` is `5, 13, 33, 83, 208, 104, 52, 26, 13, ...` (cycle minimum 13), so time 0 is a
  strict value future minimum of a non-divergent orbit (script G). For `3x+1` the claim
  is unproven: if a hypothetical cycle had odd minimum `c ≡ 2 (mod 3)`, then
  `(2c-1)/3 < c` is an odd integer with `T((2c-1)/3) = c` (script G), so its orbit would
  be eventually periodic with a strict future minimum. Repaired statement: *a strict
  value future minimum exists iff the orbit's global minimum is attained exactly once
  (i.e. the orbit diverges, or is eventually periodic with its minimum in the
  pre-periodic part); there are infinitely many iff the orbit diverges; an orbit that
  reaches `1` has none.* The other three clauses of 1(c) (coefficient future minima are
  value future minima; value lower records are height lower records; the step out of a
  value future minimum is odd) are **CONFIRMED**. Nothing downstream uses the false
  direction: Proposition 2 uses only the injectivity argument of its own part (i).
* **Remark (Terras along the chain).** **CONFIRMED**; equality of the two record sets
  reproduced for all `n <= 10^4` and the eight record orbits.

### Section 2 (Proposition 2)

* `x_j in E_inf` iff `j` is a strict coefficient future minimum. **CONFIRMED.**
* **(i)** `n in E_inf` implies injective orbit, `x_j -> infinity`, `sum 1/x_j < infinity`
  (THM-4476, hypothesis "pairwise distinct terms" satisfied), `C_infinity < infinity`,
  `h_j -> infinity`, `sum 2^(-h_j) < infinity`, infinite spine. **CONFIRMED.** The cycle
  identity `(3^o(u) - 2^|u|) m + C_u = 0` correctly excludes a positive period height.
* **(ii)** `kappa <= sigma`: **CONFIRMED unconditional** (an actual descent at `k` forces
  `3^(o_k) < 2^k` because `C_(w,k) >= 0`; my scan confirms `kappa <= sigma` to
  `2 * 10^6`). `Z^+ ∩ E_inf = {kappa = infinity} ⊆ {minima of divergent orbits}`:
  **CONFIRMED unconditional** (via (i)). Equality iff Terras's equality at the minima of
  divergent orbits; closure under `s` with `s(m) > m`: **CONFIRMED.** The parenthetical
  characterisation *"(`sigma(m) = infinity` and `m` not a cycle minimum)"* of "minimum
  of a divergent orbit" is **CORRECTED (C2)**: `sigma(m) = infinity` and "not a cycle
  minimum" also hold for a pre-periodic global minimum below a hypothetical cycle
  (`5x+1`: `m = 5`). Repaired: *`sigma(m) = infinity` and the orbit of `m` is not
  eventually periodic.* The scope note is accurate.
* **(iii)** eventually periodic positive orbits meet `E_inf` nowhere. **CONFIRMED**
  (period height negative, heights `-> -infinity`).

### Section 3 (Theorem 3)

* **(i)** unique decomposition into spine blocks, free concatenation, infinite case,
  `sum b_l t^l = 1 - 1/W = 1 - exp(-sum B_n t^n/n)`. **CONFIRMED**; own proof agrees with
  the note's; decomposition and free concatenation verified for all no-descent words of
  length `<= 14` and all pairs of blocks of length `<= 8`; both generating-function forms
  agree to `t^60` (own series inversion and own exponential); `W_k` computed by two
  routes (ballot DP and THM-4495's recurrence) agree to `k = 60`.
* **(ii)** first letter `1`, `0 < H < log_2 3 - 1`, forced ones `a(l) = ceil(l log_3 2)`,
  height `(1 - {l log_3 2}) log_2 3`, existence iff `{l log_3 2} > log_3 2` iff `l - 1`
  not a Beatty number `floor(a log_2 3)` iff `l = 1 + floor(b beta)`, `1^a 0^(l-a)` is a
  block, density `1 - log_3 2`. **CONFIRMED**, every equivalence re-derived (the interval
  `[(l-1) log_3 2, l log_3 2)` argument is right; equality is impossible by
  irrationality; Rayleigh with `1/log_2 3 + 1/beta = 1`) and verified exactly for
  `2 <= l <= 60` in all three forms plus the Rayleigh partition of `1..100`. Brute force
  over all `2^l` words for `l <= 16` and DFS to `l = 20` give
  `b_l = 1, 0, 1, 0, 0, 2, 0, 0, 7, 0, 30, 0, 0, 113, 0, 0, 525, 0, 2652, 0` and every
  block of length `2..20` starts with `1`, has `a(l)` ones and `H < log_2 3 - 1`.
* **(iii)** carry bound `C_w/3^a < a/2` and `s_(i+1) < (3/2) s_i + (3/4) l_i`.
  **CONFIRMED.** The bound needs only that the height after the `i`-th one is `>= H > 0`
  (it equals `H` when the last letter is a one, which never happens for `l >= 2` but the
  argument does not need it); verified `2 C_w < a 3^a` on every block of length `<= 20`.
  The one-letter block gives `(3 s_i + 1)/2 < (3/2) s_i + 3/4`, consistent. The same
  argument with `5` in place of `3` holds on all 510 interior blocks of the `5x+1`
  control: `2^H s_i < s_(i+1) < 2^H (s_i + l/2)` exactly (script G).

### Section 4 (Proposition 4)

* **(a)** `D(m) in [ceil(m/2), m - 1]`. **CONFIRMED** (last step is a halving of a value
  `>= m`); verified for all `2 <= m <= 2 * 10^6`.
* **(b)** shells and record bound, equality iff power of two. **CORRECTED (C3)**: the
  statement and proof assume the descent chain reaches `1` (*"the length of the descent
  chain to 1"*, *"when it passes from a value `>= 2^(j+1)` to a value `< 2^(j+1)`"*). For
  a chain ending at a hypothetical `r` with `sigma(r) = infinity` only the shells from
  `floor(log_2 m)` down to `floor(log_2 r)` are hit. Repaired: *for every `m` whose orbit
  reaches `1` (all `m <= 2 * 10^6`; all `m` if Collatz holds) ...*; unconditionally,
  `#records >= floor(log_2 m) - floor(log_2 r_last) + 1`. The induction for the equality
  case is correct. Verified: bound for all `m <= 2 * 10^6`, equality set exactly the 20
  powers of two, and every shell hit for all `m <= 3 * 10^5` by direct chain walk.
* **(c)** affine bijection `{r + t 2^k > N(w)} -> {m_0 + t 3^o}`, exactly the edges from
  `sigma = kappa` sources, in-degree formula with defect term. **CONFIRMED.** Own
  derivation: the class of `w` is the unique solution of `3^o m ≡ -C (mod 2^k)`; `m_0 =
  T^k(r)`; the threshold `m' > N(w)` is exactly `T^k m' < m'` (`m' = N(w)` is a fixed
  point, only `w = 10, m' = 1`); negative `t` are automatically excluded because such a
  source would be negative. Verified: (1) prediction from the 3-adic classes of the
  190069 words of length `<= 26` equals the direct count restricted to `sigma <= 26` for
  every landing `m <= 10^5` (2071 longer sources excluded on both sides, as in the note);
  (2) **the full formula, with sources of any `sigma` (up to 135)**: for every source
  `m' <= 40000`, its own first-descent word satisfies `m' ≡ r(w) (mod 2^k)`, `m' > N(w)`
  and `D(m') = m_0(w) + t 3^o`, so the in-degree of every landing `m <= 20000` is the
  full 3-adic count (the note's script only checks the `|w| <= 26` truncation, which is
  what its *"exact to `10^5`"* means). Also: `N(w) < 2^|w|` for every first-descent word
  of length `<= 26` (max `N = 24.54`), and the only least residue `r <= N(w)` is `r = 1`
  for `w = 10`, so no Terras-defect source exists for any word of length `<= 26`,
  regardless of the size of the source.
* **(d)** `#{m' : D(m') <= X}/X -> c_D = sum_w 3^(-o(w))`, `1 < c_D < 2`, bracket.
  Limit argument **CONFIRMED** (`R_K(X) <= 2X W_K 2^(-K) + O_K(1)`, `W_K 2^(-K) -> 0`,
  `X -> infinity` then `K -> infinity`). Two **CORRECTIONS**: (C5) *"`3^(-o) < 2^(1-k)`
  for a first-descent word (`3^o > 2^(k-1)` from the no-descent prefix)"* and
  *"`|H(w)| = k - o log_2 3 in (0, 1)`"* fail for `w = 0` (`k = 1, o = 0`: `3^0 = 2^0`,
  `|H| = 1`); true for `k >= 2` (verified on all 190069 words), and `c_D < 2` survives
  because the `k = 1` term is exactly `1` and the `k >= 2` terms are strictly below
  `2^(1-k)`. (C6) the printed derivation gives `c_D in [1.669582, 1.700498]`, and
  *"`c_D in [1.6696, 1.7005]`"* rounds the lower end **up**; the interval is nevertheless
  true: with `|w| <= 28` (502524 words) the partial sum is `1.672002` and the tail bound
  `2 W_28/2^28` gives `c_D in [1.672002, 1.698262]`.

### Section 5

* Partition coordinates, boundary `g_i = min({j - 1 : f(j) >= i} ∪ {b})`, emptiness iff
  `f(b) > a`, interval `[empty, g]`, covering = swap `10 -> 01`, carry increment
  `3^(a-i) 2^(i-1) 2^(e_i)`, strict order embedding, `(5,1)` residues `31, 47, 39, 27`
  along the chain, `{27, 31}` = top and bottom, dominance filter above the mechanical
  word, carry decreasing along dominance. **ALL CONFIRMED** (own derivation of the
  equivalence `e_(f(j)) <= j - 1` for all `j` iff `e <= g`; 870 covering pairs with the
  exact increment; residues recomputed; dominance checked for all words of length
  `<= 16`).

### Section 6 (DIRECTION; factual claims only)

* 6.1 block sequence of `(1, 1, 2)^N` is `(1)(110)` repeated. **CONFIRMED** (heights
  `0.585, 1.170, 1.755, 0.755`; future minima at times `0, 1, 4, 5, ...`).
* 6.2 Higman statistics. **CONFIRMED** (165/300, max wait 141, mean 61.01, none equal
  to 1). Wording in section 0, *"the spine has Higman-good pairs within at most 141
  blocks"*, applies only to the 165 points that found a partner inside the 300-point
  window (C10).
* 6.4 `E_inf` closed, nowhere dense, null, dimension `h*`. **CONFIRMED** (S13 Prop. 1 +
  `W_k 2^(-k) -> 0`).
* 6.5 *"the value is the least positive residue `x_f = (2^(-d_f) S_(f-1) mod 3^f)` (S8)"*
  at a spine **time** `f`. **CORRECTED (C9)**: S8's formula is indexed by the odd-iterate
  count `l` with modulus `3^l`; at a `T`-time `f` the correct statement is
  `x_f ≡ 2^(-f) C_(w,f) (mod 3^(o_f))`, equal to the least residue once `2^f > n C_f`.
  The literal `mod 3^f` is wrong whenever `3^(f-o_f)` does not divide `n` (orbit of 7:
  fails at every `f >= 4`, script G).
* 6.6 *"A dipper `i` of a landing point `j` (THM-4506) is exactly a node whose
  `⊑_D`-subtree ends at `j - 1`"*. **UNVERIFIED as stated**: dippers are defined on
  values, `⊑_D` on heights, and the two differ by the carry factor `C_k/C_i >= 1`; the
  chain statement holds verbatim if `⊑_D` is defined on values (`x_l >= 2^(-D) x_i` on
  `(i, j]`). Heuristic section; no numbers depend on it.
* Ranked step (3): *"a proof that `N(w) < 2^(|w|)` for every first-descent word ... would
  make the in-degree formula unconditional for every `m`"*. **CORRECTED (C12)**:
  `N(w) < 2^|w|` only confines a Terras defect to the least residue `r_w` of each class
  (THM-4512's one-member statement); the defect term vanishes for all `m` iff
  `r_w > N(w)` for every first-descent word, i.e. Terras's equality for all `m`.
* "Reading" of section 4: *"Collatz iff every `m >= 2` lies in some source class above
  its threshold, i.e. `kappa(m) < infinity` ... and `m > N(w_m)`"*. **CORRECTED (C11)**:
  the right side implies Collatz (then `sigma = kappa` everywhere), but Collatz does not
  imply `m > N(w_m)` (a source with `kappa < sigma < infinity` reaches `1` while lying
  below its threshold). Correct: *Collatz together with Terras's equality iff ...*; and
  *Collatz implies `kappa(m) < infinity` for all `m`* (S13 Prop. 3's direction).

### Section 7 (controls)

* Sturmian element: heights in `(0, log_2 3)`, min `0.000063`, max `1.584621`, real
  series `> J/3`, residues `0.980, 0.359, 0.744, 0.788, 0.671, 0.253, 0.005, 0.676`,
  window future-minimum gaps `{1, 19, 84, 569, 1054}` with `12/19, 53/84, 359/569,
  665/1054` above `log_3 2`. **CONFIRMED.** *"has **no** strict future minimum ..., so
  its spine is empty and Proposition 1(b)'s infinite tree does not exist for it"*:
  **CORRECTED (C8)**: time `0` is a strict future minimum of every element of `E_inf`,
  so the spine is `{0}`; the tree rooted at `0` **is** infinite (every time is a
  descendant of `0`); what fails is local finiteness (the children of `0` are the
  one-sided best-approximation times `1, 3, 11, 19, 84, 569, 1054, ...`, infinitely many)
  and hence König's hypothesis, and there is no infinite branch.
* `5x+1` orbit of 7: 512 = 512 future minima, block lengths max 174 mean 5.87, first
  twenty as printed, forced ones and height law on every interior block, **451 bits**
  (note correct; the `.out`'s `451.0` is `log_2` rounded, true value `450.973`).
  **CONFIRMED.**
* Record orbits: all eight lines (records, bound, leaders, roots = records, depth, only
  the final `1` a future minimum). **CONFIRMED.**

### Citations

* THM-4495 Step 1 = minimum decomposition into positive min-ending words with
  `W = 1/(1 - P)`, Step 2 = reversal, identity (A): **CONFIRMED** in the theorem file and
  the S9 note (lines 102-116). **Misnomer (C7)**: the header calls Step 1 *"the ladder
  decomposition"* and section 0 calls spine blocks *"THM-4495's ladder blocks"*;
  THM-4495's "ladder blocks" are its Step 3 objects, the first-passage words cut at
  ascending ladder epochs, which are the **reversals** of the spine blocks (equinumerous,
  not equal). Theorem 3's own text ("positive min-ending words `P_l`, its Step 1") is
  right.
* THM-4476: `sum 1/|x_i| < infinity` for orbits with distinct terms; no bounded strip;
  no narrow log-band; `sum 2^(d_l)/3^l < infinity`. **CONFIRMED** as attributed.
* THM-4512: cited at the path
  `01-canon/theorems/THM-4512-collatz-coefficient-descent-classes-certified-thresholds.md`
  in the note's header and in THM-4514's `related:`; **the file does not exist** (C4).
  The actual file is `01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md`.
  Its content is represented correctly (Terras's equality for odd `n <= 10^7` in the
  `T`-coding, `N(w) < 2^A` for `j <= 5000` odd steps, all-`j` needs an effective cutoff);
  the note's *"THM-4512 has it for `|w| <= 5000`"* is weaker than the source (which is
  for `o(w) <= 5000` ones), hence not wrong.
* S11 correction (synthesis section 2s) and S13 Proposition 4 (no prefix rank; the
  `(1,1,2)^N` witness): **CONFIRMED** as attributed. THM-4503's `(5,1)` example and
  `m_j`: **CONFIRMED.** THM-4506's shell lemma: **CONFIRMED.** S8's
  `m_l = 2^(-d_l) S_(l-1) mod 3^l` (Prop. 3.3, valid once `2^(d_l) > Q_l`): cited
  correctly in the header, section 1 and the section 4 reading; mis-indexed in 6.5 (C9).

## 2. Numbers recomputed (all AGREE unless marked)

* `W_1..W_24`: agree with the `.out` table; `W_26 = 1037374`.
* `b_1..b_24`: `1, 0, 1, 0, 0, 2, 0, 0, 7, 0, 30, 0, 0, 113, 0, 0, 525, 0, 2652, 0, 0, 11433, 0, 0` (agree; `b_1..b_16` also by brute force over all `2^l` words, `b_1..b_20` by DFS).
* block lengths `<= 60`: agree; `beta = 2.709511`; density on `2..60`: `22/59 = 0.373` vs `0.369`.
* `b_l/W_l` on block lengths `40..60`: min `0.070452`, max `0.112406` (agree at 4 dp).
* descent tree to `2 * 10^6`: landing range holds; max `sigma = 224` at `1126015`; no Terras defect; record bound holds; equality set = the 20 powers of two; in-degree distribution over landings `<= 10^6` identical to the note's 16 entries (sum `10^6`); total `1690291`, mean `1.690291`, max `17` at `293501`.
* first-descent words `<= 26`: `190069`; `F_k = 2 W_(k-1) - W_k`.
* partial sums (exact): `sum_(|w|<=22) 3^(-o) = 1.659807329`, `sum_(|w|<=22) 2^(-|w|) = 0.977774143`; at 26: `1.669582336` and `0.984541923 = 1 - W_26/2^26` (agree at 6 dp); tail bound `2 W_26/2^26 = 0.030916154`; bracket from 26: `[1.669582, 1.700498]` (**DISAGREE with the printed lower end `1.6696`**, see C6); from 28: `[1.672002, 1.698262]`.
* example `w = 11100`: `o = 3, C = 19, r = 23, N = 19/5, m_0 = 20`, landings `20, 47, 74, 101, 128, 155` (agree).
* sources with `sigma > 26` landing `<= 10^5`: `2071` (agree).
* cells `a + b <= 12`: `870` covering pairs (agree); `(5,1)` residues `27, 39, 47, 31` (agree).
* Sturmian: min height `0.000063` (at `j = 1054`), max `1.584621`; gaps `{1, 19, 84, 569, 1054}`; residues agree at 3 dp; `sum 2^(-h) = 12138.8`.
* `5x+1`: 512/512; max 174; mean 5.87; 451 bits (agree; `log_2 = 450.973`); Higman 165/300, 141, 61.01, 0 (agree).
* record orbits: all 32 printed integers agree.

## 3. Overclaims and misattributions (the correction list)

* **C1** Prop. 1(c), section 0 ("a positive-integer orbit has a strict future minimum
  **iff it diverges**"), THM-4514 title and (P1): false for `T_b` in general (`5x+1`,
  `n = 5`), unproven for `3x+1`; the proof sentence *"If the orbit is eventually periodic
  its global minimum is attained infinitely often"* is wrong. Repair as in section 1.
* **C2** Prop. 2(ii): *"(`sigma(m) = infinity` and `m` not a cycle minimum)"* -> *"`sigma(m)
  = infinity` and the orbit of `m` is not eventually periodic"*.
* **C3** Prop. 4(b), THM-4514 title and (P4): add the hypothesis "the orbit of `m`
  reaches `1`" (or state the shell count down to the last chain element).
* **C4** THM-4512 cited at a non-existent path (note header, THM-4514 `related:`); the file
  is `THM-4512-coefficient-descent-classes-one-member.md`.
* **C5** Prop. 4(d): `3^o > 2^(k-1)`, `3^(-o) < 2^(1-k)` and `|H(w)| in (0, 1)` need
  `k >= 2` (equality for `w = 0`); `c_D < 2` still follows.
* **C6** Prop. 4(d) and section 0, THM-4514: *"`c_D in [1.6696, 1.7005]`"* is not what the
  printed derivation gives (`[1.66958, 1.70050]`); it is true by the `|w| <= 28` sums
  (`[1.67200, 1.69826]`), which should be cited instead, or the lower end rounded down.
* **C7** header and section 0: *"ladder decomposition"* / *"THM-4495's ladder blocks"* ->
  "minimum decomposition (Step 1)"; ladder blocks (Step 3) are the reversed,
  first-passage words.
* **C8** section 7: the Sturmian point's spine is `{0}`, not empty; the tree at `0` is
  infinite, not locally finite, and has no infinite branch.
* **C9** section 6.5: `x_f = 2^(-d_f) S_(f-1) mod 3^f` at a `T`-time `f` -> `x_f ≡ 2^(-f)
  C_(w,f) (mod 3^(o_f))`, least residue once `2^f > n C_f`.
* **C10** section 0 and concept board: *"the strict lower records (the descent chain `n >
  D(n) > ...`) as its roots"*, *"roots = lower records = descent chain"*: the roots are the
  strict **height** lower records (the coefficient-descent chain); equal to the value
  chain iff Terras's equality holds along the chain (true for `n <= 10^7`). Also section
  0: *"a positive integer lies in `E_inf` iff it is the minimum of a divergent orbit"*
  drops Prop. 2(ii)'s Terras caveat (only "only if" is unconditional). And *"Higman-good
  pairs within at most 141 blocks"* holds for the 165 of 300 points with a partner in the
  window.
* **C11** section 4 "Reading": *"Collatz iff every `m >= 2` lies in some source class above
  its threshold"*: the "only if" direction needs Terras's equality.
* **C12** ranked step (3): `N(w) < 2^|w|` for all `w` does not make the in-degree formula
  unconditional; `r_w > N(w)` for all `w` (Terras's equality) does.
* Minor: the numerical face *"the 3-adic prediction matches the direct in-degree for
  every landing `m <= 10^5`"* is a check of the `|w| <= 26` truncation against `sigma <=
  26` sources only ("reconciled separately" = excluded on both sides); the full formula
  is a theorem given Terras's equality and is now verified here for all landings
  `<= 20000` with sources of every `sigma`.

## 4. What promotion of THM-4514 to PROVED still requires

1. Apply C1-C3 to the note, and to THM-4514's title and status block (the title currently
   asserts "a positive orbit has a strict future minimum iff it diverges" and the
   unconditional record bound).
2. Fix the THM-4512 path (C4) in both files; fix C5-C9 in the note; fix C10-C12 in
   sections 0, 4 and 6 (wording only, but they are the places a reader would quote).
3. Record this audit in section 9 with the script/output hashes above; the `audit:` field
   of THM-4514 should point to it.
4. Optionally replace the truncated in-degree check by the full-formula check (script D,
   all `sigma`) and cite the `|w| <= 28` bracket for `c_D`.
5. No new mathematics is needed: after the restatements every clause of THM-4514 is an
   elementary theorem, with the two conditional clauses (equality in Prop. 2(ii); the
   vanishing of the defect term) explicitly conditional on Terras's equality.
