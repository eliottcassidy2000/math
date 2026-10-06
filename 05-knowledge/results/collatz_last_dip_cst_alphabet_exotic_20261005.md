# The last-dip lemma is a corollary of Terras's coefficient-stopping-time conjecture, hence `D = ⋃ rising cones` for every `n ≤ 2.8·10^19`; the hovering-word alphabet and the explicit Diophantine residual; the order-plus-path envelope of the general-`m` window alphabet; exotic spheres (28, 992, 16256), Dynkin trees and the perfect-number factor; the lemniscate / `|x|` / `x^3 sin(1/x)` / `x·1_Q` quartet

**Session:** mac-mini (claude) 2026-10-05, worktree `codex/session-trees-tournaments-20261005`, third wave
(after [trees5_tournaments4_erdos_sos_20261005](trees5_tournaments4_erdos_sos_20261005.md) and
[sumner_t5_petersen_collatz_20261005](sumner_t5_petersen_collatz_20261005.md)).
**Owner's seed:** prove the last-dip lemma for `l ≥ 6` via route (c); work toward the general-`m` alphabet;
consider deeply the Poincaré conjecture, exotic spheres, 992 and 16256, Milnor's `S^7` (28 total); merge the
lemniscate of Bernoulli and `|x|` at 0 against `x^3 sin(1/x)` and `x·1_Q`.

**Status.** PROVED (elementary): Prop. 1.1 (last-dip violation ⟹ failure of Terras's coefficient stopping
time conjecture, CST), Prop. 1.2 (the hovering alphabet: `A = ⌈l log_2 3⌉`, the Beatty gate, `n`-window),
Prop. 2.1 (window classes lie in the order-plus-path family `F_m`). CITED: Rozier–Terracol 2025
(arXiv:2502.00948, Cor. 5.2: CST holds for `2 ≤ n ≤ 2.8·10^19`), Terras 1976, Garner 1981, Lagarias 1985;
Kervaire–Milnor 1963; Milnor 1956; Brieskorn 1966. FINITE-EXACT: the C sweeps below `2^28` and `2^32`
(143,745,904 and 2,300,119,872 last-dip segments, no violation); `|F_m| = 2, 4, 11, 37, 143, 622` for `m = 3..8`; the Cartan
determinants; the Kervaire–Milnor orders `k ≤ 7`. Typed verdicts in §5. Collatz OPEN; the Lemma's residual is
explicit (§1.4) and is a statement about orbits with stopping time above `6.6·10^9`. Audit OWED on Props. 1.1–1.2.
Scripts in `04-computation/experiments/`: `collatz_last_dip_sweep_20261005.c`, `collatz_stopping_time_records_20261006.c`, `collatz_last_dip_residual_20261005.py`,
`collatz_last_dip_residual_cf_20261005.py`, `syracuse_window_alphabet_order_plus_path_macmini_20261005.py`, `order_plus_path_F8_macmini_20261005.py`,
`trees5_dynkin_kervaire_milnor_20261005.py`; outputs `*.out` in this directory (§7).

## 0. Inheritance and concept board

- Last-dip lemma (previous note §4): an excursion above `n` that starts at the orbit's last value `m* < n`
  before `n` has a rising word (`2^A < 3^l`); it implies `D = ⋃ rising cones` (the S6/S7 OPEN question).
  Proved there: `a_1 = 1`, `(2n-1)/3 ≤ m* < n`, all proper suffixes non-rising, `c_w/3^l ≤ n((1+1/(3n))^l - 1)`,
  `l ≤ 5` given a sweep to `2^20`.
- THM-4512 (`coefficient-descent-classes-one-member`): coefficient descent vs actual descent; "universal
  coefficient stopping and universal stopping-time equality remain OPEN"; only the least representative of a
  coarse class can fail the threshold test.
- Terras 1976: stopping time `t(n) = min{j : T^j(n) < n}`, coefficient stopping time
  `τ(n) = min{j : C_j(n) < 1}`, `C_j(n) = 3^{q_j}/2^j` (`q_j` = odd terms among the first `j`), `T(n) = (3n+1)/2`
  or `n/2`; always `τ ≤ t`. **CST conjecture** (named by Lagarias 1985): `t(n) = τ(n)` for all `n ≥ 2`; it implies
  that no nontrivial cycle exists. Verified: Terras to `250,000`, Garner 1981 to `1,150,000`, Rozier–Terracol
  2025 to `2.8·10^19` (Cor. 5.2, from Theorem 5.1: exactly 593 "paradoxical" sequences with `7 ≤ n ≤ 4614`
  and none for `4615 ≤ n ≤ 2.8·10^19`, using the convergence verification and the excursion/delay records).
- S12/S13 (opus, today): the Syracuse `m`-window alphabet (time arcs + numerical order, gauge DESC), Haar law
  `11:5:5:11` at `m = 4`, classes `4, 8, 25, 64, 194` for `m = 4..8`, condensation multiplicativity, Ryser
  twins (`syracuse_window_alphabet_20261005.out`, `syracuse_window_alphabet_general_m_20261005.out`).
- THM-868 (`e8-bridge-score-lattice`): every 8-tournament's score deviation `(2s_v - 7)/2` is an `E_8` vector
  (half-integer coset); THM-067 (`c1-formula-mersenne`): `c_1 = 2^{f+1} - d - 2` vanishes at `n = 2^{f+1} - 1`
  (1, 7, 31, 127, 511); HYP-2559 already typed the perfect-number echo as "same arithmetic engine, distinct loci".
- The only prior exotic-sphere work: `04-computation/surgery_exotic_23.py` (opus S84, 2026-03-14, untyped
  numerology in the "(2,3) lens"; its `bP` table is wrong for odd `k`, §3.3 and the MISTAKE entry).
- Board: {CST; hovering words; convergents of `log_2 3`; order + Hamiltonian path; plumbing along a tree;
  Euclid's factor `2^{p-1}(2^p - 1)`; a corner at 0 vs infinitely many oscillations vs a single point of continuity}.

## 1. The last-dip lemma for `l ≥ 6`

### 1.1 Proposition 1.1 (PROVED). A last-dip violation is a counterexample to Terras's CST conjecture.

Let `m* < n` be odd with Syracuse orbit `m* = n_0 -> … -> n_l = n`, `n_j > n` for `1 ≤ j < l`, and word `w` with
`2^A > 3^l`. In the `T`-orbit of `m*` (one `(3x+1)/2` step and `a_j - 1` halvings per Syracuse step; `A` steps in
all) every value after time 0 is a power-of-two multiple of some `n_j` with `j ≥ 1`, or equals `(3m*+1)/2 = 2^{a_1-1} n_1`,
hence exceeds `n > m*`; the value at time `A` is `n > m*`. So `t(m*) > A`. The coefficient at time `A` is
`C_A(m*) = 3^l/2^A < 1`, so `τ(m*) ≤ A`. Hence `τ(m*) < t(m*)`: CST fails at `m*` (`m* ≥ 3`). ∎

**Corollary 1.1.** CST ⟹ last-dip lemma ⟹ `D = ⋃ rising cones` (previous note, Prop. 4.2). Since CST holds for
`2 ≤ m* ≤ 2.8·10^19` (Rozier–Terracol Cor. 5.2, CITED) and `m* < n`, **every `n ∈ D` with `n ≤ 2.8·10^19` lies in a
rising cone**; the two sets agree below `2.8·10^19`. The independent C sweep (`collatz_last_dip_sweep_20261005.c`)
confirms this directly below `2^28` (143,745,904 last-dip segments) and below `2^32` (2,300,119,872 segments):
none non-rising, extremal word `84/53` (`2^84/3^53 = 0.997914`, first at `93119 -> 93317`), reproducing the
Python sweep at `2^20`.

**Reading.** The non-rising ancestor words that motivated the question (`165 -> 167`, `l = 17`, `A = 27`) are
exactly Rozier–Terracol's *paradoxical sequences* (`C_j(n) < 1` with `T^j(n) ≥ n`): there are 593 of them below
4615 and none up to `2.8·10^19`. A last-dip violation is a paradoxical sequence that is also a CST failure (no
dip below the start before the paradoxical moment); their Theorem 5.1 says every paradoxical sequence below
`2.8·10^19` has already dipped. The lemma is therefore *weaker* than CST: it concerns only paradoxical moments
that are the first return to a value above the start.

### 1.2 Proposition 1.2 (PROVED). The hovering alphabet.

For a violation `(w, m*, n)` with `g = n - m* ≥ 2` and `δ = 2^A/3^l - 1 > 0`:

1. **Exact slope.** `2^A/3^l = (m*/n) ∏_{j=1}^{l} (1 + 1/(3 n_{j-1}))`, so
   `1 < 2^A/3^l < (1 + 1/(3m*))(1 + 1/(3n))^{l-1} ≤ (1 + 2/(3n-3)) e^{(l-1)/(3n)}`. Hence **`A = ⌈l log_2 3⌉`**
   and `δ = δ_l := 2^{1 - {l log_2 3}} - 1`, with the **Beatty gate** `1 - {l log_2 3} < log_2((1 + 2/(3n-3)) e^{(l-1)/(3n)})`
   `≈ (l+2)/(3 n ln 2)`. (In particular `{l log_2 3} > log_2(4/3)` for every hovering word, since `n/m* < 3/2 + o(1)`.)
2. **The gap identity.** `g = c_w/3^l - δ_l n` and `X_w - n = g/δ_l`, `X_w - m* = g(1+δ_l)/δ_l`, where
   `X_w = c_w/(2^A - 3^l)` is the word's positive cycle point: `m* < n < X_w`, and for every proper suffix `w'`
   (non-rising) `n > X_{w'}`. So a violation is an excursion trapped **below the positive cycle point of its own word
   and above those of all its suffixes** — S8's Proposition R ("a non-rising word descends iff `x > x_w`") as an
   obligation. For `w = (1) ∘ w'`: `X_w - X_{w'} = (X_{w'} + 1)/(3 ε)` with `ε = δ/(1+δ)`, so `2 ≤ g < (X_{w'} + 1)/3`.
3. **The `n`-window.** `2 + δ_l n ≤ c_w/3^l ≤ n((1 + 1/(3n))^l - 1)`, so `n ≤ N_l` (the largest solution) and
   `n < l/(3 ln(1 + δ_l))`.
4. **Few letters, one source.** `c_w ≥ 2^{A - a_l}` gives `A - a_l ≤ l log_2 3 + log_2 (l/3) + o(1)`; the source is
   `m* ≡ X_w (mod 2^A)` and, by THM-4512's one-member bound, it must be the least representative of its coarse class
   (a violation is an actual non-descent of a coefficient-descent word).

### 1.3 Route (c): what it gives and where it stops

Route (c) was: "a word of length `l` fixes `m* mod 2^A`, so long near-balanced excursions should force large
sources against `n < l/(3 ln(1+δ_w))`". Items 1–4 make this exact: the word is forced to be the *best* upper
approximation at its length (`A = ⌈l log_2 3⌉`), the source is the least representative of the class, and `n`
is boxed into `[X, N_l]` once sources below `X` are excluded. What route (c) does **not** give is a lower bound on the
least representative `r_w` of a hovering word: `r_w` is a 2-adic residue `X_w mod 2^A` with no known regularity, and
for the words that survive the Beatty gate the window `[X, N_l]` is nonempty. The crude bound of the previous
note (`c_w/3^l ≤ (l/3) e^{l/(3n)}`) proves the Lemma for `l ≤ 3` outright and for `l ≤ 5` given any sweep past
`2^20`; from `l = 6` on the content is Diophantine, exactly as for CST itself.

### 1.4 The explicit residual (FINITE-EXACT via the continued fraction of `log_2 3`)

`log_2 3 = [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 15, …]`. The Beatty gate at `n ≥ X` forces
`1 - {l log_2 3} < (l+2)/(3 X ln 2)`. The residual length set is therefore the **Bohr set**
`B(X) = {l : 1 - {l log_2 3} < (l+2)/(3 X ln 2)}`; since the gate is linear in `l`, sums and small multiples of
qualifying lengths qualify too (at `X = 2^28`: `15601, 31202, 46803, 47468 = 15601 + 31867, 62404, 63069, …`).
Below `√(3X ln 2/2)` the gate is stronger than `1/(2l)`, so there `l` must be an upper convergent denominator; hence the
**smallest** residual length is the first upper convergent denominator past that bound. The table lists the smallest
lengths and the convergent/intermediate-fraction members with nonempty windows `[X, N_l]`:

| after sources `≤ X` are excluded | residual lengths `l` (odd steps of the excursion) | first window |
|---|---|---|
| `X = 2^20` (Python sweep) | 2291 lengths `≤ 10^5`, smallest `l = 2966` (`A = 4701`) | `n ∈ [2^20, 1.16·10^6]` |
| `X = 2^28` (C sweep) | 12 lengths `≤ 10^5` (exact Bohr set); smallest `l = 15601` (`A = 24727`), then 31202, 46803, 47468, …, 79335; convergent members `≤ 10^14`: 15601, 47468, 79335, 190537, 10781274, … | `n ∈ [2^28, 2.86·10^8]` |
| `X = 2^32` (C sweep) | one length `≤ 10^5`: `l = 79335` (`A = 125743`, `δ_l = 3.7·10^{-6}`); then 190537, 10781274, … | `n ∈ [2^32, 7.2·10^9]` |
| `X = 2.8·10^19` (CST, Rozier–Terracol) | smallest **`l = 6,586,818,670`** (`A = 10,439,860,591`, `δ_l = 1.0·10^{-11}`), then its small multiples and 72057431991, 137528045312, … | `n ∈ [2.8·10^19, 2.2·10^20]` |

**Independent stopping-time route (FINITE-EXACT, "wanted entries only").** A violation in the window `(l, A, [X, N_l])`
needs a source `m* < N_l` with `T`-stopping time `σ_T(m*) > A`. The record sweep `collatz_stopping_time_records_20261006.c`
gives `max σ_T(m) = 447` for `m < 2^32` (record holder `2788008987`; `183` below `2^20`, `287` below `2^24`, `395` below
`2^28`), against `A ≥ 4701` for every residual length at `X = 2^20`, `A ≥ 24727` at `2^28`, `A = 125743` at `2^32`. So
every residual window with `N_l ≤ 2^32` is empty without the CST input, and the two routes (direct last-dip sweep,
stopping-time maxima) agree below `2^32`. The record is unchanged below `2^33` (still `447`), which empties the
`l = 79335` window `[2^32, 7.2·10^9]`: the lemma holds for all `n < 2^33` with no CST input
(`collatz_five_papers_synthesis_20261006.md` §5.1).

So any violation of the last-dip lemma (and any non-cone element of `D`) is an orbit that stays above its
target for at least `6.6·10^9` odd steps, i.e. has stopping time above `10^10`, starting from `m* > 2.8·10^19`.
No present theorem excludes orbits with stopping time `≫ log n`; none has been observed (the record stopping
times below `2^68` are in the low thousands). The Lemma is thus proved up to the same Diophantine residue as
CST and the "logarithmic deadline" (S7). Typed: **PROVED for `n ≤ 2.8·10^19` (CITED input), FINITE-EXACT
residual, OPEN beyond** — a statement about the existence of astronomically long excursions along the upper
convergents of `log_2 3`.

## 2. Toward the general-`m` alphabet: the order-plus-path envelope

**Setting (S12).** A window `(x, Sx, …, S^{m-1}x)` of the Syracuse map, read as the tournament with time arcs
`x_i -> x_{i+1}` and every other pair oriented by size (DESC: larger `->` smaller). In the Haar/large-`x` regime
the window depends on the valuation word only through the order pattern of `v_j = j log_2 3 - A_j`.

**Proposition 2.1 (PROVED).** Every window tournament has its time path as a Hamiltonian path whose complement is
the restriction of a linear order (hence acyclic). Define
`F_m := {classes of m-tournaments having a Hamiltonian path whose complement is acyclic}` ("order + path",
the tournaments `T(π)`, `π ∈ S_m`, obtained from `TT_m` by rerouting the arcs along one Hamiltonian path).
Then `W_m^{DESC}, W_m^{ASC} ⊆ F_m`, `F_m` is closed under converse, and `F_m` contains no regular tournament for
`m = 5, 7` (checked; for `m = 3` the 3-cycle is in `F_3`). A letter above `C_m = ⌈(m-1) log_2 3 - (m-2)⌉` cannot
change the pattern, so the Haar law is an exact dyadic rational (reproduces S12: `11/32, 5/32, 5/32, 11/32` at
`m = 4`; `19/64, 13/64, 1/8, 1/8, 5/64, 5/64, 3/64, 3/64` at `m = 5`; `73/512, 33/256, …` at `m = 6`).

| `m` | labeled DESC windows (= S13's "chord patterns") | `|W_m^{DESC}|` | `|W_m^{ASC}|` | `|F_m|` | all classes |
|---|---|---|---|---|---|
| 3 | 2 | 2 | 2 | 2 | 2 |
| 4 | 5 | 4 | 2 | 4 | 4 |
| 5 | 16 | 8 | 7 | 11 | 12 |
| 6 | 49 | 25 | 26 | 37 | 56 |
| 7 | 145 | 64 | 70 | 143 | 456 |
| 8 | 458 (S13) | 194 (S13) | — | 622 | 6880 |

`W_m = F_m` exactly for `m ≤ 4`; from `m = 5` the staircase constraint (`v_j` from a slope-`log_2 3` integer
staircase) bites: `3, 12, 79, 428` order-plus-path classes are not windows at `m = 5, 6, 7, 8`. The two gauges give
different class alphabets (S12 §D: converse of a DESC window is the ASC window of the *reversed* sequence, a
backward orbit). `|F_m| = 2, 4, 11, 37, 143, 622` (`m = 3..8`) is not in OEIS (checked 2026-10-05); the window alphabet is a shrinking
fraction of it (`1, 1, 8/11, 25/37, 64/143, 194/622`). S13's census (`m ≤ 8`:
`4, 8, 25, 64, 194` classes; condensation multiplicativity; Ryser twins) is the finer object; `F_m` is its
combinatorial envelope and the natural place to ask which permutations `π` are staircase patterns: the
excluded `F_m \ W_m` classes are the "non-staircase" permutation tournaments.

## 3. Exotic spheres, 28 / 992 / 16256, and the trees

### 3.1 The exact bridge: trees are plumbing graphs (CITED, classical; FINITE-EXACT census here)

Plumbing `D^4`-bundles over `S^4` of Euler number 2 along a tree `T` gives an 8-manifold with intersection form
the Cartan matrix `C(T) = 2I - Adj(T)`; its boundary is a homotopy 7-sphere iff `det C(T) = ±1`. The three trees on
5 vertices are **`A_5` (path), `D_5` (fork), affine `D̃_4` (star `K_{1,4}`)** with `det = 6, 4, 0`: their plumbings
bound a rational homology sphere with `|H_3| = 6`, one with `|H_3| = 4`, and a manifold with infinite `H_3`. Among the
23 trees on 8 vertices exactly one is unimodular, `E_8`; exactly three are positive definite (`A_8`, `D_8`, `E_8`);
`E_8`'s plumbing boundary is Milnor's generator of `Θ_7 ≅ Z/28` (Milnor 1956 via `S^3`-bundles; Brieskorn 1966
realises all 28 as `Σ(2,2,2,3,6k-1)`). So the only tree whose plumbing bounds an exotic sphere has 8 vertices,
and none of the 5-vertex trees does. This is the same object as THM-868's lattice: the score vectors of
8-tournaments live in `E_8`; here `E_8` is the intersection form of the `E_8`-tree plumbing. Same lattice, two
carriers (a tree on 8 vertices; a tournament on 8 vertices), no map between the carriers is claimed.

### 3.2 The numbers (CITED Kervaire–Milnor 1963; FINITE-EXACT `k ≤ 7`)

`|bP_{4k}| = 2^{2k-2}(2^{2k-1} - 1)·num(4B_k/k)` (`B_k = |B_{2k}|`): `28, 992, 8128, 261632, 1448424448, 67100672` for
`k = 2..7`; `Θ_7 = bP_8 = Z/28`, `Θ_11 = bP_12 = Z/992`, `Θ_15 = bP_16 ⊕ Z/2`, `|Θ_15| = 16256`. The factor
`2^{2k-2}(2^{2k-1} - 1)` is Euclid's `2^{p-1}(2^p - 1)` with `p = 2k - 1`, a perfect number exactly when `2^p - 1` is a
Mersenne prime: `p = 3, 5, 7` give `28, 496, 8128`; `p = 9, 11` do not (`511, 2047` composite); `p = 13` does
(`33550336`, `|bP_28| = 67100672`). So `992 = 2·496 = σ(496)` and `16256 = 2·8128 = σ(8128)`, but the two doublings
have different sources (`num(4B_3/3) = 2` vs the `Z/2` of `coker J` in dimension 15): EXPLAINED COINCIDENCE of the
Euclid factor with perfect numbers, NUMEROLOGY for `σ`.

### 3.3 Relation to our objects (typed)

- THM-067's vanishing loci `n = 2^{f+1} - 1` (1, 7, 31, 127, 511) are the Mersenne factors of `|bP_{4k}|` with
  `f = 2k - 2`: the shared engine is `2^p - 1`; HYP-2559's verdict ("same arithmetic engine, distinct `n`-loci, not
  the same event") stands. NUMEROLOGY.
- `28 = C(8,2)` = the arcs of the Sumner host; `E_8` has 8 vertices; `|Θ_7| = 28`. No map. NUMEROLOGY.
- The genuine inheritance is structural, not numerical: the three 5-vertex trees are `A_5, D_5, D̃_4`; the end-defect
  Klein group of the first note (fold the two ends of `A_5`) moves `A_5 -> D_5 -> D̃_4`, i.e. from the finite Dynkin
  path through the finite `D_5` to the affine star — each fold lowers `det C` from `6` to `4` to `0`. The affine
  diagram is where the plumbing stops being a rational homology sphere, and it is the star, the tree that alone
  sees the Erdős–Sós half and the Sumner degree wall. DICTIONARY (exact), no Collatz content.
- The prior exotic-sphere script `surgery_exotic_23.py` (S84) is wrong at every odd `k` (it applies both `a_k = 2`
  and the numerator; prints `|bP_12| = 1984`, `|bP_20| = 523264`, `|bP_28| = 134201344`), hard-codes `bP_20 = 130816`,
  mislabels Milnor's `λ` as "signature mod 7", and calls `π_1(Σ(2,3,7))` finite; logged as a MISTAKE entry today.

## 4. The quartet at 0: lemniscate, `|x|`, `x^3 sin(1/x)`, `x·1_Q`

The four objects are the textbook archetypes of behaviour at a single point, and each has an exact counterpart
among this session's objects. Everything here is a DICTIONARY (typed per line); no theorem crosses.

| archetype at `0` | what fails / holds at `0` | this session's counterpart | type |
|---|---|---|---|
| `|x|` | continuous, Lipschitz, no derivative: one-sided slopes `±1`, a corner | the two Collatz sheets `(3n ± 1)/2 = (3n + sgn n)/2`: the sheet is the sign, i.e. the derivative of `|n|`, undefined exactly at the fixed point `0` (THM-4474(E): sheet involution = negation) | EXACT restatement |
| lemniscate `(x^2+y^2)^2 = 2a^2(x^2-y^2)` | a node at `0` with tangent lines `y = ±x`: two smooth branches crossing; the upper branch pair *is* the corner of `|x|` to first order | the two lobes exchanged by `x -> -x` are the two sheets exchanged by `n -> -n`, meeting at `0` (the "vanishing" class of the pentagon trichotomy, THM-4521) | DICTIONARY |
| `x^3 sin(1/x)` | `C^1` with `f'(0) = 0`, but `f''(0)` does not exist: infinitely many oscillations of shrinking amplitude; smooth-looking, not smooth | an orbit that hovers: finitely many oscillations above `n` before landing (§1); the archetype says what an *infinite* hover would look like — the residual of §1.4 is exactly the question whether the Collatz walk can oscillate `6.6·10^9` times within a shrinking band | ANALOGY |
| `x·1_Q` (= `x` on `Q`, `0` off `Q`) | continuous at `0` only; squeezed by `|x|`; the difference quotient takes the values `0, 1` | the rational/irrational dichotomy of the parity-vector conjugacy `Φ` (Lagarias's periodicity conjecture, the repo's PC): `Q ∩ Z_2` ↔ eventually periodic words. But `Φ` is a 2-adic homeomorphism (barrier atlas), so nothing is "continuous only at 0" there: the archetype's *discontinuity* has no counterpart, only its *dichotomy* does | ANALOGY (weak) |

Two things are more than analogy. (i) *Bernoulli twice.* The lemniscate is Jacob Bernoulli's curve (1694); the
exotic-sphere orders of §3 use Jacob Bernoulli's numbers (1713). The genuine mathematical relation between the
two is the circle/lemniscate pair of period lattices: `B_{2k}` are the Laurent coefficients of `x/(e^x - 1)` and
give `ζ(2k) ∈ π^{2k} Q`, the periods of the circle; the lemniscatic sine `sl` (arc length `∫ dr/√(1-r^4)`, the CM
curve `y^2 = x^3 - x`, lemniscate constant `ϖ = 2.62205…`, `AGM(1, √2) = π/ϖ`) has as its Laurent data the
**Hurwitz numbers**, the `Z[i]`-analogue of the Bernoulli numbers, giving the Eisenstein values of `Q(i)` in `ϖ^{4k} Q`.
The Kervaire–Milnor formula lives entirely on the circle side (image of `J`, `K`-theory Adams operations);
there is no lemniscatic exotic sphere — that would be a `tmf`-level statement, and none is claimed here. CITED
(Hurwitz 1899; Kervaire–Milnor 1963). (ii) *The corner and the merged node.* `|x|` is `max(x, -x)`: the
non-smooth selection of the two smooth branches of the lemniscate's node. The first note's merged class `[±]`
is the same selection in the tournament world: forgetting which of the two converse diamonds (forward map /
inverse tree) one looks at keeps only `u = 1`, the number of contracting moves — the absolute value of the
sheet sign. EXACT as a restatement of THM-4472 + Prop. 1 of the first note.

### 4.1 Cross-references (repo sweep, 2026-10-05)

The repo had the two halves of the quartet in two threads that were never joined; `x^3 sin(1/x)` and `x·1_Q`
occur nowhere before this note.

- **`|x|` at 0 against the Collatz sign law** — `collatz_procgen_20260924_fixed_points_kawasaki_audit.md` §3
  (mac-mini, 2026-09-24), from the owner's phrases "the sharp corner at 0 introducing orthogonality" and "`|x|` as
  the precursor of polynomials": Proposition L (PROVED) `|T_+(x)| = T_{sgn x}(|x|)`, "the corner of `|x|` at 0 is
  exactly where the selector jumps", "smoothing the corner erases the sign" (even observables are side-blind);
  Theorem I (PROVED): `x ≥ 0` gives the upper equality case, `x ≤ 0` the lower, with the owner's middle point `z = 0`;
  verdicts REAL (sign law as the two equality cases with 0 in the middle) and ANALOGY (orthogonality; "precursor").
  Also THM-4471 (`kawasaki-fixed-point-collatz-proof-refuted`, line 114), the mod-192 note's Theorem 6 / Corollary 7
  ("counted by `|x|` it cannot separate the sheets", PROVED), and `MISTAKES.md:676` ("`x -> -x` preserves `|x|`").
  The `|x|` row above restates Proposition L; it adds only the reading "sheet = derivative of `|n|`".
- **The lemniscate node as the Lonely-Runner pinch** (July 2026, LRC only, never Collatz): kps-S95
  (`the-lemniscate-and-the-second-moment-of-gaps-kps-S95.md`; owner's `(x^2+y^2)^2 = x^2 - y^2` as "a strange source
  of inspiration"; Cue 1 REAL = the second-moment bound `maxgap ≥ Σgap^2/Vmax`, Lean; Cue 3 "thematic": the origin
  crossing is a `Z_2` quotient, "the same complement symmetry as the merged tournament metagraph `G_n/Z_2`" — the
  present "corner and merged node" paragraph is this cue made exact through THM-4472 and the first note's
  Prop. 1(d)); kps-S107 (`grid-invisible-pinches-are-lemniscate-nodes`): the maxgap function has a corner with slopes
  `-12 -> +9` at the pinch, "the two crossing lines `y = ±x` … the same local singularity", smooth surrogate = corner
  (Fourier decay `1/m^2`) vs sharp indicator = jump (`1/m`), Lean `LRCPinch.lean` (pinches are rational `m/d`);
  opus-S177/S169: node = collision = exact resonance = arc edge, the tight point `M = 1/14` ("one object, four views";
  HYP-5547/5660 CONFIRMED). The trigger was MISTAKE-130: the retracted tent `maxgap = 1 - spread·|x|` near `x = 0`.
- **Lemniscatic constants without `|x|` content:** THM-3012 (the spherical/Euclidean wall `4/k + |1 - 4/k|` has its
  corner at `k = 4`, the lemniscatic case; `S(4)` lemniscatic), THM-3560 (the lemniscatic curve "not as decorative
  analogies"), HYP-9075 ("an adjacency of tools, not a link"), kps-S146–S148 (S148 retracts S146/147's `ϖ` claims).
  **Name collision:** "Hurwitz number" in the repo means the automorphism bound `42(g-1)` / `84(g-1)`; the lemniscatic
  Hurwitz numbers of this section are a different object (Hurwitz 1899).
- **Rational vs irrational:** `the-rational-irrational-duality.md` (kps-S31, LRC; `lonely_abs_iff` folds every sign
  choice onto `|v|`); the PC/Φ work (transversality foundry, hard-class Theorem D; HYP-9127); the barrier atlas records
  `Φ` as a 2-adic homeomorphism — hence the weakening of the `x·1_Q` row.
- The owner's SOURCE files contain none of the four; the nearest owner phrasing is the glued-affine blueprint's
  signed reflection law `U_b(-n) = -U_{-b}(n)` across the origin, which is the `|x|` row in the owner's own words.

**Merged verdict.** The quartet is a dictionary of regularity-at-a-point for the sheet structure: corner (`|x|`,
REAL identities already in canon), node (lemniscate, the same `Z_2` quotient as the merged metagraph, EXACT as a
statement about involutions, "inspiration" as a tool), infinite oscillation (`x^3 sin(1/x)`, the open residual of §1
in one picture), and the rational/irrational dichotomy (`x·1_Q`, weak). It changes no status.

## 5. Verdicts

| claim / question | verdict | carrier |
|---|---|---|
| last-dip lemma, `l ≥ 6` | PROVED for all `n ≤ 2.8·10^19` as a corollary of CST (CITED verification); FINITE-EXACT direct sweep to `2^28`; explicit residual (`l ≥ 6.6·10^9`, `n > 2.8·10^19`); OPEN beyond | Prop. 1.1 + Rozier–Terracol Cor. 5.2 |
| `D = ⋃ rising cones` | holds for all `n ≤ 2.8·10^19`; OPEN beyond, same residual | Prop. 4.2 of the previous note |
| route (c) | makes the word, source and `n`-window exact (Prop. 1.2); cannot bound the least representative | 2-adic class vs Beatty gate |
| general-`m` alphabet | `W_m ⊆ F_m` (PROVED); `|F_m| = 2,4,11,37,143,622`; equality iff `m ≤ 4`; no regular tournament in `F_5, F_7` | order + Hamiltonian path |
| 28 / 992 / 16256 vs our work | Euclid factor explained; `σ` doubling NUMEROLOGY; `C(8,2)` NUMEROLOGY; trees = Dynkin plumbing graphs EXACT; `A_5 -> D_5 -> D̃_4` under the end-defect folds EXACT | Kervaire–Milnor; Cartan determinants |
| lemniscate / `|x|` / `x^3 sin(1/x)` / `x·1_Q` | see §4 | — |

## 6. Open items and cheapest next tests

- Prove the last-dip lemma outright: it is strictly weaker than CST, so a proof may exist where CST's does not.
  The cleanest target is Prop. 1.2(2): show no excursion can stay between `max_{w'} X_{w'}` and `X_w` while
  following `w`; equivalently bound the least representative of a hovering class from below by `l/(3 ln(1+δ_l))`.
- Extend the direct sweep to `2^36` in C (hours at the current rate); `2^32` is done and leaves `79335` as the only
  residual length below `10^5` for the sweep-only route, independent of the CITED input.
- The characterisation of `F_m \ W_m` as non-staircase permutation tournaments; whether `|F_m|` (`2, 4, 11, 37, 143, 622`)
  has a formula (permutations modulo the symmetries of the order-plus-path presentation).
- Plumbing dictionary: which trees on `n ≤ 10` vertices are unimodular (`E_8`, `E_10 = T(2,3,7)`, …) and whether the
  end-defect folds have a meaning for the plumbed manifolds (fold = blow-down?).

## 7. Reproduction and hashes

```
cc -O2 -o lastdip 04-computation/experiments/collatz_last_dip_sweep_20261005.c && ./lastdip 268435456 > 05-knowledge/results/collatz_last_dip_sweep_2e28_20261005.out   # ~ 1 min
./lastdip 4294967296 > 05-knowledge/results/collatz_last_dip_sweep_2e32_20261005.out   # ~ 20 min
python3 04-computation/experiments/collatz_last_dip_residual_20261005.py 100000 > 05-knowledge/results/collatz_last_dip_residual_20261005.out            # ~ 1 min
python3 04-computation/experiments/collatz_last_dip_residual_cf_20261005.py > 05-knowledge/results/collatz_last_dip_residual_cf_20261005.out
python3 04-computation/experiments/syracuse_window_alphabet_order_plus_path_macmini_20261005.py 7 > 05-knowledge/results/syracuse_window_alphabet_order_plus_path_macmini_m7_20261005.out   # ~ 8 min; 6 is instant
python3 04-computation/experiments/order_plus_path_F8_macmini_20261005.py > 05-knowledge/results/order_plus_path_F8_macmini_20261005.out   # ~ 1 min
python3 04-computation/experiments/trees5_dynkin_kervaire_milnor_20261005.py > 05-knowledge/results/trees5_dynkin_kervaire_milnor_20261005.out
```

Controls: the C sweep reproduces the Python sweep at `2^20` (561,601 segments, same extremal segment); the
alphabet script reproduces S12's `m = 4, 5` Haar laws and S13's labeled-window counts `5, 16, 49, 145`; the
Kervaire–Milnor orders match the literature values for `k = 2..7` (`|bP_24| = 1448424448`, `|Θ_23| = 48·|bP_24|`);
the Cartan determinants match `n+1, 4, 3, 2, 1, 0` for `A_n, D_n, E_6, E_7, E_8`, affine. Hashes (SHA-256, raw
bytes) are in the commit metadata of this session's checkpoints; the `.out` files are the frozen outputs.

References (CITED). R. Terras, *A stopping time problem on the positive integers*, Acta Arith. 30 (1976).
C. Garner, *On the Collatz 3n+1 algorithm*, Proc. AMS 82 (1981). J. Lagarias, *The 3x+1 problem and its
generalizations*, Amer. Math. Monthly 92 (1985). O. Rozier, C. Terracol, *Paradoxical behavior in Collatz
sequences*, arXiv:2502.00948 (2025), Thm. 1.1/5.1, Cor. 5.2. M. Kervaire, J. Milnor, *Groups of homotopy
spheres I*, Ann. Math. 77 (1963). J. Milnor, *On manifolds homeomorphic to the 7-sphere*, Ann. Math. 64 (1956).
E. Brieskorn, *Beispiele zur Differentialtopologie von Singularitäten*, Invent. Math. 2 (1966). OEIS A001676.
