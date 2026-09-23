# Procedural search for sign-specific order laws of Syracuse orbits: 6,016 candidates, nothing beyond the sign law

**Status: FINITE-EXACT for the grammar search: 6,016 candidate order statements evaluated on every odd `n <= 10^6` (search range) and every odd `n <= 10^7` (extension) on both sheets. The C engine and an independent pure-Python recomputation (exact integers and Fractions) agree exactly at `N = 20001` on every output, and on 6,000 random per-`n` rows at `10^7`. PROVED (elementary): the direction lemma (on the plus sheet only decay windows can contradict their word's order prediction; on the minus sheet only growth windows can); order = word + crossings, and the transfer theorem (a window-local order law that holds on the plus sheet and fails on the minus sheet can fail only in a window that contains a minus gate crossing); the ladder-epoch (Sturmian) lemma; the cross-orbit laws for `R_b(n) = 4n + b` and the directional stopping-time monotonicity along residue classes. PROVED (computer-assisted: census plus an explicit gate bound `G_b(l)`): the complete list of gate crossings with lag `<= 38` on the plus sheet (12 windows at lag 17, 56 at lag 29) and lag `<= 35` on the minus sheet (`165, 309, 549` at lag 12); short-lag genericity (every comparison at lag `< 17` on the plus sheet, `< 12` on the minus sheet, equals its word prediction); `sigma = tau` for `tau(n) <= 38` (plus) and `sigma(n) <= 35` (minus); the minus-sheet Sturmian record law for gaps `<= 36`. OPEN: `sigma = tau` in general on either sheet (the plus sheet case is the accelerated form of Terras's coefficient stopping-time conjecture; the minus mirror is inherited as OPEN from the control lane's S3); completeness of the crossing lists at larger lags. RESULT: no sign-specific order law beyond the sign law was found. Each of the 158 one-sided survivors is a window sign-law instance (16), a Sturmian crossing instance (12), terminal (5), a normalization artifact (80) or a small-`N` artifact (45); 10 more one-sided tail artifacts appear at `10^7`. None is unexplained. No independent-agent audit yet.**

Session `collatz-procgen-20260922` (machine `mac-mini`), lane `order_laws`, 2026-09-22. Script [collatz_procgen_20260922_order_laws.py](../../04-computation/experiments/collatz_procgen_20260922_order_laws.py) with C engine [collatz_procgen_20260922_order_laws.c](../../04-computation/experiments/collatz_procgen_20260922_order_laws.c); output [collatz_procgen_20260922_order_laws.out](collatz_procgen_20260922_order_laws.out). Section numbers `S0`-`S7` refer to the `.out`.

## 0. Inheritance

* **Closest proved mechanisms.** The first is the word-function theorem ([w6 probes S6a](collatz_mod6_20260922_w6_sign_specific_probes.md), inherited from the [control lane S2](collatz_mod6_20260922_minus_sheet_positive_control.md)): an invariant that depends only on `(b, word)` cannot tell the sheets apart. The second is the sign law `sign(B) = b` for the carry in `2^K T^L(n) = 3^L n + B` ([portrait S2](collatz_mod6_20260922_counterexample_portrait.md)). The third is the control lane's S3: "coefficient descent implies actual descent on the minus sheet and the converse holds on the plus sheet". That statement is the start-window case of the direction lemma in section 3.1. The same S3 already records, as OPEN, the minus mirror of Terras's coefficient stopping-time conjecture (verified there to `10^6`, and for T-map first-descent times `<= 20`).
* **Canonical hostile.** The minus seven-cycle `{17, 25, 37, 55, 41, 61, 91}` (clock `11/7`). Here it is joined by the minus *transient* near-cycle `165 -> ... -> 163`, twelve odd steps with clock `19/12`. That orbit segment is the smallest minus-sheet witness against every plus-sheet order law found.
* **Corrected near miss.** The task statement records one sign-specific order law: the minimum of a positive plus cycle is `3 mod 4` and the maximum `1 mod 4`, reversed on the minus sheet. In signed residues `b*x mod 4` the two sheets say the same thing: `b*min = 3` and `b*max = 1 mod 4`, because the minimum ascends (`k = 1`) and the maximum descends (`k >= 2`) on both sheets. So the "reversal" is an artifact of reading residues without the sign; it is not sign content. The grammar rediscovers this 68 times as raw-residue record statements (class NORMALIZATION, section 4.5). For example, the predecessor of every running-max record is `3, 7, 11, 15 mod 16` on the plus sheet and `1, 5, 9, 13` on the minus sheet, but in signed residues it is `3, 7, 11, 15` on both (checked at `10^7`).
* **Least-used sidecar.** The Beatty/Sturmian structure of first-passage lags: `floor(g log2 3) - floor((g-1) log2 3)` in `{1, 2}`. With it, record gaps become a sheet-detecting statistic (section 4.3).

## 1. Objects

* **Sheets.** `U_b(x) = (3x + b)/2^(v_2(3x + b))` on positive odd `x`, with `b = +1` (plus) or `b = -1` (minus).
* **Stopped orbit.** `x_0 = n, ..., x_T` is followed until the first element of `C_+ = {1}`, respectively `C_- = {1, 5, 7, 17, 25, 37, 41, 55, 61, 91}`. Its values are pairwise distinct. A continuation of 16 cycle steps is used only to evaluate `sigma` and `tau`.
* **Word.** `k_i = v_2(3 x_i + b)` and `K_j = k_0 + ... + k_(j-1)`.
* **Window.** A window `(i, i+l)` has lag `l` and exponent sum `K = K_(i+l) - K_i`. It is a GROWTH window if `3^l > 2^K` and a DECAY window if `2^K > 3^l`; `2^K = 3^l` never happens. Write `2^K x_(i+l) = 3^l x_i + B`, where `B = b * sum_(m<l) 3^(l-1-m) 2^(K_(i+m) - K_i)`. The window's GATE is `g = B/(2^K - 3^l)`, the fixed point of its affine map. The gate is positive exactly for decay windows on the plus sheet and growth windows on the minus sheet.
* **Crossing and margin.** A window is a (gate) CROSSING if its actual comparison `x_(i+l)` vs `x_i` differs from the word prediction (up iff growth). The margin of a positive-gate window is `mu = x_i/g - 1`, and a crossing is `mu < 0`.
* **Stopping times.** `sigma(n)` is the first `j` with `x_j < n`. `tau(n)` is the first `j` with `2^(K_j) > 3^j` (the coefficient stopping time, accelerated convention).
* **Records.** A running-max record is an index `j >= 1` with `x_j > x_i` for all `i < j`; a running-min record is defined the same way. `fmr(n)` is the first running-max record, `rmax` and `rmin` count the records, and `rho(n) = #{j : x_j > n}`.
* **Signed residues** `b*x mod M` are used throughout. With them every word function is the same function on the two sheets, because the word bijection sends plus class `r` to minus class `-r`.

## 2. The grammar and the evaluation (S2, S3)

Each candidate has the form "for every positive odd `n <= N`, feature `F` of the stopped orbit has property `P`". It is evaluated on both sheets.

| family | features and properties | candidates | plus-only | minus-only | both | neither |
|---|---|---|---|---|---|---|
| P | ordinal pattern (rank permutation) of windows of `w = 3..6` values; all `w!` patterns; scopes: all windows, initial window, windows not containing `x_T`, windows with all values `>= 1000`; property "never occurs" | 3480 | 0 | 0 | 2640 | 840 |
| L | lag comparisons `x_(i+l)` vs `x_i`, `l = 1..64`, from the start (`i = 0`, i.e. `n` vs `x_l`) or anywhere; "growth window goes up", "decay window goes down", "never up", "never down"; plus the all-lag growth/decay statements | 516 | 6 | 10 | 244 | 256 |
| G | gate margins, equivariant: "every positive-gate window has `mu > 0`, resp. `>= c`", `c = .01 ... 1` | 7 | 0 | 0 | 0 | 7 |
| R | running max/min records: signed and raw residues of the record value and its predecessor mod 16 and mod 27, the step `k_(j-1)`, record gap `1..64`, record index `1..128`; scopes `j >= 1`, `j < T`; "never takes value `v`" | 1344 | 64 | 44 | 242 | 994 |
| S | per-`n` features `sigma, tau, T, argmax, rho, rmax, rmin, fmr`: "never equals `v`", `v = -1..64` (`-1` = infinity) | 528 | 5 | 13 | 99 | 411 |
| X | cross-orbit pairs `(n, m)`: same word through `sigma(n)` (`m = n + 2^(K_sigma + 1)`), same stopped word (`m = n + 2^(K_T + 1)`), `m = R_b(n) = 4n + b` (shared image), raw `m = 4n + 1` (non-equivariant control), `m = n + 2`; nine features; `<=`, `>=`, `=` | 135 | 15 | 1 | 28 | 91 |
| C | `sigma = tau`, `sigma >= tau`, `sigma <= tau`, `sigma < infinity`, `tau < infinity`, `sigma = tau` at every orbit point | 6 | 0 | 0 | 6 | 0 |
| total | | **6016** | **90** | **68** | 3259 | 2599 |

Classification of the 158 one-sided survivors at `10^6`. Every label was verified mechanically on the witnesses (S3, S6).

| label | plus-only | minus-only | meaning |
|---|---|---|---|
| SIGN-LAW | 6 (L) | 10 (L) | "growth windows grow" (plus) or "decay windows decay" (minus): PROVED on its sheet by `sign(B) = b`; the other sheet fails exactly at its gate crossings |
| SIGN-LAW-STURMIAN | 0 | 12 (8 R, 4 S) | record gap or first-record index `g` that no word can produce; the plus sheet realizes it only through a crossing |
| TERMINAL | 5 (1 S, 4 X) | 0 | fails only at `n = 33, 49`, whose minus orbits enter `C_-` at `61` before dropping below their start (names the cycles) |
| NORMALIZATION | 45 (34 R, 11 X) | 35 (34 R, 1 X) | raw residues or the raw map `4n + 1`; the signed (equivariant) versions are sheet-blind |
| SMALL-N | 34 (30 R, 4 S) | 11 (2 R, 9 S) | tail values (record index `92..128`, gap `62`, `k = 19..22`, `rmax = 28..35`, `rmin = 25..28`, `fmr = 57..64`); the other sheet realizes each by `10^7` |
| UNRESOLVED | 0 | 0 | |

At `10^7` a further 10 statements that held on both sheets at `10^6` become one-sided (`k_(j-1) = 22..25` for min records, `rmax = 39, 40`, `rmin = 31..34`). They are tails of the same kind: each first appears between `10^6` and `10^7` on one sheet only. The maxima of these features at `10^7` agree closely: `rmax` 43/43, `rmin` 34/30, `fmr` 91/93 (plus/minus).

## 3. Structural results

### 3.1 Direction lemma (PROVED)

For every window, `x_(i+l) - x_i = (B - (2^K - 3^l) x_i)/2^K` with `sign(B) = b`. Therefore:

* on the **plus** sheet a growth window always goes up, and a decay window goes up iff `x_i < g`;
* on the **minus** sheet a decay window always goes down, and a growth window goes down iff `x_i < g`.

Equality `x_(i+l) = x_i` means `x_i = g`, a periodic point, so it cannot occur in a stopped orbit. Consequently the complete order relation of a stopped orbit is determined by its word and its crossing set. On the plus sheet every crossing is a decay pair that goes up (type U); on the minus sheet every crossing is a growth pair that goes down (type D). The `i = 0` case is the control lane's S3.

### 3.2 Transfer theorem (PROVED)

Let `F` be a set of window types `(v, pi)`: a word `v` of length `l` and a permutation `pi` of `l + 1` points. Suppose no plus stopped orbit, for any `n`, contains a window of a type in `F`, and a minus window has type `(v, pi)` in `F`. Then that minus window contains a crossing pair.

*Proof.* Suppose the minus window has no crossing. Then `pi` is the generic pattern of `v`, the one predicted pairwise by the word. Take `n` in the plus class of `v` (one odd class mod `2^(K+1)`), with `n` larger than `2^K` times the largest gate of any sub-window of `v`. Along the window `x_m(n) > 3^m n/2^(K_m) >= n/2^K` because `B_m > 0`, so no sub-window crosses and no `x_m` equals `1`. So the plus sheet realizes `(v, pi)` inside a transient, contradicting the hypothesis. QED. The statement with the sheets swapped holds by the same argument, choosing `n` large enough for the minus carries.

**Consequence.** Every window-local order law that is true on the plus sheet and false on the minus sheet fails only in minus windows that contain one of the minus crossings of section 3.3. Its sign content is therefore the sign law. Global features (records, `sigma`, pair statements) can also separate the sheets through the terminal structure, because stopped orbits end in `C_b`. The TERMINAL class is exactly that.

### 3.3 Complete crossing lists and short-lag genericity (PROVED, computer-assisted; S4)

**Gate bound.** For a word of length `l` and exponent sum `K`, each `K_m <= K - (l - m)`. Hence `|B| <= B_max(l, K) = 3^(l-1) + 2^(K-l+1)(3^(l-1) - 2^(l-1))`, with equality for the word `(K-l+1, 1, ..., 1)`. A crossing at lag `l` on sheet `b` therefore has `x_i < G_b(l)`, where `G_b(l)` is the maximum of `B_max/|2^K - 3^l|` over the positive-gate side. For the plus sheet, which has infinitely many `K`, the value is a finite maximum plus an explicit monotone tail bound. A brute-force check over all words with `k_i <= 11`, `l <= 6`, confirms that the bound dominates every true gate.

**Completeness argument.** A crossing window from `x_i` is a prefix of the stopped orbit of `x_i`. So if `G_b(l) <= 10^7`, every crossing of lag `l` starts at a census start and has been seen.

* `G_+(l) <= 10^7` for `l <= 38` (`G_+(39) = 4.05e7`). The **plus crossings with lag `<= 38` are exactly the 12 windows of lag 17 and the 56 of lag 29**.
* `G_-(l) <= 10^7` for `l <= 35` (`G_-(36) = 3.51e7`). The **minus crossings with lag `<= 35` are exactly `165 -> 163`, `309 -> 307`, `549 -> 547` at lag 12**.
* **Short-lag genericity.** On transients every comparison at lag `l < 17` (plus) or `l < 12` (minus) equals its word prediction; the largest gate bound used is 2012 (plus) and 162.5 (minus). In particular every window of at most 12 consecutive odd iterates has its generic pattern on both sheets. The first sheet-separating window has 13 values.

**Pattern family, exactly.** The generic (word-realizable) ordinal patterns number 6, 16, 48 and 140 for `w = 3, 4, 5, 6` (of 6, 24, 120, 720). Words with `k_i <= 9` suffice, since any `k >= 8` already makes every window containing it a decay window. At `10^6`, in every scope and on each sheet, the realized set equals the generic set (checked). So the whole P family, 3,480 candidates, is word-level: 2,640 statements are true on both sheets (forbidden patterns) and 840 are false on both.

### 3.4 Ladder-epoch lemma (PROVED; S4b)

Suppose `x_j > x_i` and `x_m < x_i` for `i < m < j` (a first passage above `x_i` at lag `g = j - i`), and no pair `(i, i+m)` is a crossing. Then `2^(K_m) > 3^m` for `0 < m < g`, `2^(K_g) < 3^g`, and the last exponent is `1`. This forces `floor(g log2 3) - floor((g-1) log2 3) = 2`. The excluded lags have density `2 - log2 3 = 0.415`: `3, 5, 8, 10, 13, 15, 17, 20, 22, 25, 27, 29, 32, 34, 37, 39, 41, 44, 46, ...`. A dynamic program over exponent words gives the same set for `g <= 64`. The dual statement for first passage below has no obstruction, since a descent can use any `k >= 2`.

## 4. The survivors, with witnesses and status

### 4.1 Window sign law: "every growth window grows" (L, plus-only, 6 statements)

* **Statements.** At lag 12, at lag 53, and at all lags, each both from the start and anywhere.
* **Plus sheet.** PROVED for all `n`: `2^K x_(i+l) - 3^l x_i = B > 0`.
* **Minus sheet.** False exactly at its crossings, 15 distinct windows at `10^6` and again at `10^7`. Each was re-verified in exact arithmetic: the word, the identity `2^K x_l = 3^l x_0 + B`, and `0 < x_i < g`.
  * Clock `19/12` (convergent `p4/q4`, `2^19/3^12 = 0.98654`): `165 -> 163`, `309 -> 307`, `549 -> 547`. This list is complete (`G_-(12) = 6291`).
  * Clock `84/53` (convergent `p6/q6`, `0.997914`): `2981 -> 2953`, `4471 -> 4429`, `6005 -> 5995`, `6705 -> 6643`, `17537 -> 17497`, `24933 -> 24911`, `26305 -> 26245`, `29593 -> 29525`, `39457 -> 39367`, `44325 -> 44287`, `44389 -> 44287`, `44433 -> 44287`. Complete to `10^7`, not certified beyond (`G_-(53) = 6.9e11`).
* **Smallest minus witness.** `165` with word `(1,2,1,1,2,2,1,1,5,1,1,1)`: `165 -> 247 -> 185 -> 277 -> 415 -> 311 -> 233 -> 349 -> 523 -> 49 -> 73 -> 109 -> 163`.
* **Reduction.** This is the sign law read on windows. What is new is the anatomy (section 5): the minus sheet's transient near-cycles sit at its next growth-side convergents after the cycle clocks `1/1, 3/2, 11/7`.

### 4.2 The dual: "every decay window decays" (L, minus-only, 10 statements)

* **Minus sheet.** PROVED (`B < 0`).
* **Plus sheet.** False at 148 crossing windows, the same set at `2*10^4`, `10^6` and `10^7`:
  * clock `27/17`: 12 windows, `165 -> 167`, `171 -> 175`, `231 -> 233`, `257 -> 263`, `259 -> 263`, `387 -> 395`, `389 -> 395`, `391 -> 395`, `437 -> 445`, `581 -> 593`, `587 -> 593`, `589 -> 593`; complete;
  * clock `46/29`: 56 windows, `x_i` in `[313, 3073]`; complete;
  * clock `65/41`: 62 windows, `x_i` in `[313, 3031]`;
  * clock `73/46`: 18 windows, `x_i` in `[1071, 3075]`.

  `27/17` and `46/29` are the intermediate fractions between `8/5` and `65/41`. `73/46` is not a best approximation: it is the mediant of `27/17` and `46/29`, and 9 of its 18 windows contain both a lag-17 and a lag-29 crossing as sub-windows (S4).
* **Status.** A statement true on the minus sheet and false on the plus sheet: the wrong direction for a Collatz proof, and again the sign law.

### 4.3 Sturmian record law (R and S, minus-only, 12 statements)

* **Statement.** Consecutive running-max records are never `g` apart, and the first record is never at index `g`, for `g = 17, 29, 41, 46`. More generally, for every excluded `g` of section 3.4.
* **Minus sheet.** PROVED for `g <= 36`, for all `n`. An excluded gap needs a crossing `(x_i, m)` with `m < g` whose start is a running maximum. For `m <= 35` the only minus crossings are `165, 309, 549` at lag 12 (section 3.3), and a running maximum `x_i` forces `n <= x_i <= 549`, inside the census. FINITE-EXACT for every `g` to `10^7`: no excluded gap of any size occurs among the 63 realized gap values up to 112.
* **Plus sheet.** False at exactly the excluded gaps `17, 29, 41, 46`, the crossing lags. First witnesses by gap: `171` (record `257` then `263` after 17 steps), `327` (`2369 -> 2429`, 29), `1287` (`2897 -> 3077`, 41), `2035` (`3053 -> 3077`, 46). First witnesses by first-record index: `257, 685, 2897, 3053`. Every witness window is one of the plus crossings of section 4.2 (verified).
* **Status.** The word part is PROVED. The sheet asymmetry is the sign law, since only a decay window that goes up can create an excluded first passage.

### 4.4 Terminal (5)

* **Statements.** `S rmin=0` ("every stopped orbit has a new minimum") and `X SPT sigma/tau >=, =`.
* **Minus failures.** Only at `n = 33` (`33 -> 49 -> 73 -> 109 -> 163 -> 61`) and `n = 49`. Both enter `C_-` at `61` before dropping below their start, so `sigma(n) > T(n)`.
* **Status.** These statements only name the cycles.

### 4.5 Normalization (80)

* **Raw residues (68).** The raw-residue record statements come in mirror pairs `r <-> -r`. They are the trivial residue law of the task statement, in record form (section 0). Every signed-residue record statement is two-sided (both or neither).
* **Raw `4n + 1` (12).** The map `m = 4n + 1` shares images with `n` only on the plus sheet. On the minus sheet `U_-(4n + 1) = 6n + 1`, and the equivariant partner is `R_-(n) = 4n - 1`.

### 4.6 Both-sheet order laws (flagged: order statements that are not word functions by construction)

* **Coefficient stopping time `sigma(n) = tau(n)` (C, 6 statements).**
  * FINITE-EXACT on both sheets for all 5,000,000 plus and 4,999,991 minus non-cycle odd `n <= 10^7`, and at every orbit point where both times occur inside a stopped orbit.
  * PROVED for `tau(n) <= 38` (plus) and `sigma(n) <= 35` (minus). A failure is a crossing from the start at lag `tau` (plus) or `sigma` (minus), hence `n < G_b(l) <= 10^7`.
  * Minus cycle elements: `sigma = tau = infinity` at `1, 5, 17`, and equal and finite at the others.
  * Relation to Terras. The accelerated statement is implied by Terras's T-map conjecture but weaker: the two T-map times need only fall in the same halving run. It already implies that no nontrivial positive `3n+1` cycle exists, because at the minimum `tau <= L` (decay clock) while `sigma = infinity`. The minus cycles satisfy it because their words never turn decay.
  * Inheritance. The plus statement is Terras's conjecture (accelerated form; attribution UNCITED-RECOLLECTION, as in the control lane). The minus statement and its `10^6` verification are inherited from the control lane S3. New here: `10^7`, the orbit-point version, and the certified ranges in odd steps (the control lane has first-descent `<= 20` T-steps).
  * **CST margin** (answers the w6 note's stopping question numerically). Among all first-decay windows `(0, tau(n))` with `n <= 10^7` on the plus sheet, exactly one has `n < 2g`: `n = 63` (`tau = 34`, `K = 54`, `x_34 = 61`, `mu = 0.747`), an orbit that rides the 27-trajectory. On the minus sheet no growth prefix before the first descent has `mu < 1`. Both lists are complete for lags `<= 33` (plus) and `<= 35` (minus), where `2G_b(l) <= 10^7`.
* **Cross-orbit (X, 28).**
  * `R_b(n) = 4n + b`. `U_b(4n+b) = U_b(n)` because `3(4n+b) + b = 4(3n+b)`. Hence `T` is equal, `sigma(m) = 1 <= sigma(n)`, `tau(m) <= tau(n)` (the word is `(k_0 + 2, k_1, ...)`), `M(m) >= M(n)`, `argmax(m) <= argmax(n)`, `rho(m) <= rho(n)`, `rmax(m) <= rmax(n)`, `rmin(m) >= rmin(n)` and `fmr(m) >= fmr(n)`. All PROVED and sheet-blind.
  * Same prefix through `sigma(n)`. Along the shared word `x_j(m) - x_j(n) = 3^j (m-n)/2^(K_j) > 0`, so order is preserved (PROVED, sheet-blind; also checked). The descent set inside a residue class is upward closed on the plus sheet and downward closed on the minus sheet (the gate direction). So `sigma(m) <= sigma(n)` is PROVED on plus and `sigma(m) >= sigma(n)` is PROVED on minus: a sign-specific direction. Equality holds on both sheets over 4.73 million pairs each at `10^7`, and is a consequence of `sigma = tau`.
  * Same stopped word (SPT). `M(m) >= M(n)` is PROVED; `T(m) >= T(n)` is PROVED on plus. The rest are FINITE-EXACT on only 380 (plus) and 3,919 (minus) pairs at `10^7`, because such pairs exist only for short stopped words: weak evidence.
* **L (244).** "Growth windows grow" and "decay windows decay" at every lag with no crossing. These are order statements; they are PROVED for the certified lags and FINITE-EXACT beyond.

## 5. Anatomy of the failures (S4)

* **Dip and return.** Every crossing window, on either sheet, dips strictly below both endpoints and then climbs back to within a few units of its start. Examples: `165 -> 167` on plus (dip `31`), and `165 -> 163`, `309 -> 307`, `549 -> 547` on minus (dips `49, 65, 65`). The dip values are `{31, 47, 71, 91, 103}` on plus and `{49, 65, 229, 481, 961, 1249}` on minus. This is `sigma = tau` in window form: no crossing happens before the first dip.
* **Clocks on the sheet's own side of `log2 3`.**
  * Plus crossings use decay-side clocks `27/17, 46/29, 65/41` (and the composite `73/46`). The only positive plus cycle `{1}` sits at `2/1`; `8/5` hosts nothing, since a lag-5 crossing would need `x_i < 86`.
  * Minus crossings use the growth-side convergents `19/12, 84/53`, the next ones after the cycle clocks `1/1, 3/2, 11/7`.
  * So the near-cycles are the failed cycles at the next admissible clocks. Which side is admissible is decided by the sign law.
* **The 27-highway.** All 148 plus crossing windows share values with the orbit of 27, and every plus crossing ends on it, at `167, 175, ..., 2429, 3077`. Because that highway is heavily travelled, 1,170,879 of the 5,000,001 plus orbits up to `10^7` contain a crossing, against 15,255 minus orbits. This is a frequency asymmetry of the tree near the root, not a law.
* **Near misses.** Windows with `0 < mu < 1`: 170 on plus (closest `3079` at lag 46, `mu = 0.0107`) and 14 on minus (closest `809` at lag 12, `mu = 0.221`). The equivariant margin family G fails on both sheets for every `c`: both sheets have gate crossings, and they differ only in which clocks carry them.

## 6. The three best candidates

1. **Window sign law, plus-only.** "If `3^l > 2^(K_(i+l) - K_i)` then `x_(i+l) > x_i`." PROVED on the plus sheet. False on the minus sheet at transient near-cycles: smallest witness `165`, 12 odd steps to `163` (clock `19/12`); the lag-12 list `{165, 309, 549}` is PROVED complete; 12 more windows at lag 53. Plus sheet has 0 failures to `10^7` (as it must). It is the sign law restated on windows.
2. **Sturmian record law, minus-only.** "No running-max record gap and no first-record index equals an excluded `g` (`floor(g log2 3) - floor((g-1) log2 3) = 1`)." PROVED on the minus sheet for `g <= 36`, FINITE-EXACT to `10^7` for all `g`. False on the plus sheet at `g = 17, 29, 41, 46` (first witnesses `171, 327, 1287, 2035`), each witness window a plus crossing. It is the sign law again, in the minus-true direction.
3. **Coefficient stopping time `sigma(n) = tau(n)`, both sheets.** It is order-based, not a word function. FINITE-EXACT to `10^7` on both sheets and at every orbit point. PROVED for `tau <= 38` (plus) and `sigma <= 35` (minus), with CST margin `>= 0.747` (worst case `n = 63`). OPEN in general. It is sheet-blind as a statement. Its power to exclude plus cycles comes entirely from the sign law: plus cycle words are decay words, while minus cycle words never become decay.

## 7. Skeptical assessment

* **No new sign-specific law.** The search found no sign-specific order law beyond the sign law, and the transfer theorem explains why for window-local laws. A window law that separates the sheets must fail at a minus crossing, and the direction of every crossing is `sign(B) = b`. The only remaining channels are the terminal structure (which names the cycles) and non-equivariant normalizations.
* **Most survivors are artifacts.** 45 small-`N` tails at `10^6` (all realized on the other sheet by `10^7`) and 10 new tails at `10^7`: statistics near the end of their range, where each sheet's first occurrence is a matter of which residue class (`r` or `-r`) is hit first. 80 normalization artifacts. 5 terminal.
* **Nothing at bounded window length.** The ordinal-pattern family, which the task statement suggested as the natural grammar, is completely barren. By short-lag genericity every window of at most 12 odd iterates has its word's generic pattern on both sheets, and this is PROVED, not just observed.
* **The crossing sets look finite, but that is not proved.** They are unchanged from `10^6` to `10^7` (plus: unchanged from `2*10^4`). Heuristic reason: a crossing needs a dip to a value of order `1/(3|2^K/3^l - 1|)`, then a climb back to the start. On the plus sheet, climbs from small values are capped by the 27-trajectory (odd maximum 3077), and every plus crossing ends on it. Completeness is certified only for lags `<= 38` (plus) and `<= 35` (minus). The crude bound `G_b(l)` is very loose: actual crossing starts are `<= 3075` while `G_+(41) = 9.7e8`.
* **Weakly tested pair statements.** The SPT pair statements rest on few pairs (short stopped words only), and the per-lag L statements stop at `l = 64`. The all-lag L statements and the crossing census have no lag limit.

## 8. Reproduction

```bash
cd <worktree>
ORDER_LAWS_RERUN=1 python3 04-computation/experiments/collatz_procgen_20260922_order_laws.py \
    > 05-knowledge/results/collatz_procgen_20260922_order_laws.out
# rerun: reuses the C outputs in scratch/procgen_order/ when present (~30 s), else regenerates them:
python3 04-computation/experiments/collatz_procgen_20260922_order_laws.py | \
    diff - 05-knowledge/results/collatz_procgen_20260922_order_laws.out      # identical (checked, also under -O,
                                                                             # and from a clean scratch dir)
```

* **Scratch files.** The per-`n` binaries (`*_pern.bin`, 200 MB each at `10^7`) and the compiled engine are regenerable. `scratch/procgen_order/.gitignore` excludes them, so they are never committed. They were deleted after the final check, so the first rerun takes the full 2.5 minutes.

* **What runs.** The script compiles the C engine (`cc -O2`, into `scratch/procgen_order/order_laws`) and runs it at `N = 20001`, `10^6 + 1` and `10^7 + 1` on both sheets: about 72 s (plus) and 41 s (minus) at `10^7`, peak engine memory under 10 MB. It then loads the outputs (peak RSS about 1.05 GB, from the two `10^7` per-`n` arrays), recomputes everything at `20001` in pure Python and checks a random sample at `10^7`.
* **Time and checks.** The full run takes about 2.5 minutes. Every check raises explicitly, and timing goes to stderr.
* **Engine usage.** `order_laws b N prefix [lagmax]` writes the pattern, lag, crossing, margin, record and `sigma/tau` tables and a per-`n` binary file.

## 9. Stopping boundary / next question

The lane stops at the classification: every survivor is explained, and the only order-level carrier of sign is the gate-crossing set, whose direction is the sign law. Two concrete next steps:

1. **Sharpen the gate bound.** Use the realizability constraint `x_m >= 1` inside the window, for example by a DP over words that tracks the smallest admissible start. This should certify the crossing lists, and with them `sigma = tau`, at lags 41, 46 and 53. `G_b(l)` overestimates the true crossing starts by five orders of magnitude.
2. **Test the uniform CST margin at the word level.** Is `n >= (1 + c) g` for every stopping word, with some `c > 0`? The census gives `c = 0.747` (attained at `n = 63`, a 27-trajectory rider) to `10^7`. A proof of any uniform `c > 0` would prove the accelerated coefficient stopping-time law, hence no nontrivial positive `3n+1` cycle. It is therefore at least as hard as the no-cycle problem, and it is sign-specific only through the sign law.
