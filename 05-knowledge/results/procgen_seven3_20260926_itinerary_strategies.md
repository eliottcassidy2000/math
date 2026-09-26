# Itinerary-coded sign strategies for 7n±1: an exact variable-depth evaluator, rigidity of rejoining flips, searches, adversaries, and the limit object

**Status.**
- **Lane question** (is 7n±1 bounded-lookahead provable, i.e. is `lim rho*(7,k) < log_7 2`?): **still OPEN.** Coding strategies and adversaries by the max-halving (MH) itinerary gave new tools and structure, but it did not move either bound.
  - *Upper side.* No rule below `37/100` was found. The best new explicit rules reach the known optima at their levels, down to `7/18`.
  - *Lower side.* No adversary above `1/3` at every level was found. The all-level floor stays `1/3` (Theorem N of THM-4486).
- **PROVED** (hand proofs in §§2–4; the runner re-checks their finite content):
  - *Fact I.* Itinerary cylinders are exactly the residue classes, so finite itinerary trees are finite bit trees. Rules that read the recent *past* of the orbit are not sign strategies. The genuinely new objects are itinerary-*regular* (cyclic) rules, which read unboundedly many symbols and are truncated at a bit depth.
  - *Lemma M (Markov refinement).* The exact `rho_max` of any variable-depth rule is the maximum mean cycle of the leaf graph of its Markov refinement. No uniform level `2^k` is needed: it evaluates rules of depth 46 with 1.1·10^4 nodes.
  - *Lemma R (rigidity of rejoining)* and *Corollary R.* If, after a flip, the orbit of `x` meets the MH orbit of `x` again, it does so after exactly the same numbers of odd steps and of halvings, except for `x` in a countable set of rationals. So the local gains of the flip calculus (Corollary F of THM-4508) are paid back in full on rejoining. Consequence: a rule with `rho_max < 1/2` must move the orbit of every generic point of `S_inf` permanently off its MH orbit.
  - *Lemma S (shadowing).* After the gaining flip `(s,2),(s,2),(-s,·)`, the orbit sits at `E = (x_2 - s)/2`. `E` copies the valuations of `x_3, x_4, ...` for as long as those symbols have sign `-s` and valuation 2 or 3.
  - *Family bound.* The rules "flip `(s,2)(s,2)` and alternating `v = 2` runs of length `>= D`" have `rho_max >= D/(2D+1)` for every `D`, with equality FINITE-EXACT for `D = 3..15` (`D = 2` gives `1/2`).
- **FINITE-EXACT** (runner; certificates checked exactly):
  - rho*(7,k) reproduced with an independent game engine for `k = 8..20`: Stern–Brocot from scratch for `k <= 16`, both certificates at every level.
  - Exhaustive searches: pattern sets; all itinerary automata with at most 2 transient states over `{=2, !2, >=3}`; all with 1 transient state over `{2, 3, >=4}`. Best values: `3/7` (pattern sets and the 2-state automata) and `8/17` (the 1-state automata).
  - Counterexample-guided itinerary search from MH: `2/5 -> 15/38 -> 7/18`, reaching `rho*(7,16) = 7/18` with a sparse rule of depth 46. It stalls there: 350 more moves stay at `7/18`, and re-solving Min's game on the rule's own Markov partition fails below `7/18`.
  - Optimal level-k rules (k = 12..20): their Markov refinements have 16–25% of `2^k` nodes (no sparsity).
  - Adversaries: top lift exactly `1/3` at k = 8..18; d-scaled top lifts, MH-density lifts and hill-climbed itinerary-prefix lifts are all `<= 1/3` at `k >= 12`.
- **EMPIRICAL:**
  - *Forced decisions.* Bit cylinders explain the forced decisions of the canonical optimal potentials better than itinerary cylinders with the same number of cells. The itinerary coding is *not* the natural coordinate system of the optima.
  - *Single flips.* A single flip at a random point of `S_inf` never rejoined (200 + 3000 trials), and the orbit afterwards has typical MH statistics (mean valuation 3.00).
  - *Abstract-arena refinement.* Counterexample-guided refinement of abstract arenas (CEGAR) did not reach `15/38` in 400 rounds.
- **OPEN:** `lim rho*(7,k)` versus `log_7 2`, and an all-level floor above `1/3`. No HYP or THM file was created. Collatz is untouched.

Session `collatz-procgen-20260922`, lane **seven3**, 2026-09-26.
- **Scripts** (`04-computation/experiments/`): `procgen_seven3_20260926_{lib,itree,search,run}.py`, `procgen_seven3_20260926_{markov,game,adv}.c`, data `procgen_seven3_20260926_rules.json`.
- **Reused read-only:** `procgen_seven2_20260926_{rhomax,restrict}.c`.
- **Output:** [procgen_seven3_20260926.out](procgen_seven3_20260926.out).
- **Parents:** [THM-4508](../../01-canon/theorems/THM-4508-structure-of-optimal-sign-strategies-flip-calculus-and-7n-plus-minus-1-at-level-31.md) (Lemmas C, MH, F), [THM-4486](../../01-canon/theorems/THM-4486-min-max-cycle-density-game.md) (the game, Lemmas G1, G2, L, Theorem N), and the notes [seven2](procgen_seven2_20260926_seven_structure.md) and [floor](procgen_floor_20260926_density_floor.md).

## 0. Answers in brief

| task (lane brief) | answer | status |
|---|---|---|
| A. upper side: itinerary rules below 2/5, 37/100, log_7 2 | Below 2/5: yes. A counterexample-guided search over itinerary trees, starting from MH, reaches `2/5`, `15/38` and then `7/18 = rho*(7,16)` with a sparse rule of depth 46. It stalls there: 350 more moves stay at `7/18`, and re-solving Min's game on the rule's own Markov partition fails below `7/18`. Below 37/100: **no**. Small itinerary automata and short pattern sets never go below `3/7`. Optimal rules are not sparse in either coding. | FINITE-EXACT; EMPIRICAL (search) |
| B. lower side: all-level adversaries above 1/3 | **None found.** Top lift is exactly `1/3` at k = 8..18. Negative-rational worlds (`-c/d`), MH-density lifts and hill-climbed itinerary-prefix lifts are `<= 1/3` from k = 12 on; one reaches `14/41` at k = 10 but only there. The all-level floor stays `1/3` (Theorem N). | FINITE-EXACT; OPEN (uniform floor) |
| C. the limit object | `lim_k rho*(7,k) = inf over locally constant sign functions sigma of max over T_sigma-invariant measures of mu(odd)`. In itinerary coordinates: `inf` over clopen flip sets `F` of the densest periodic orbit of "shift off `F`, flip map on `F`". `< log_7 2` is certified by one finite rule plus a potential; `>= log_7 2` needs a uniform family of Lemma G2 certificates. Corollary R shows that typical orbits cannot decide it: the obstruction lives on exceptional rational coincidences. | PROVED (formulation, Corollary R); OPEN (value) |
| D. arena cross-check | rho*(7,k) = 3/7, 3/7, 2/5 (x4), 15/38, 15/38, 7/18, 7/18, 19/49, 13/34, 34/89 for k = 8..20, reproduced by an independent engine | FINITE-EXACT |

## 1. Setting

Setting as in THM-4486 and THM-4508, with `q = 7`.
- **Sign strategies.** A sign strategy (rule) `sigma` gives each odd 2-adic integer a sign. It is *locally constant* (LC, "finite") if it is constant on the classes of a finite partition into residue classes. Its *depth* is the largest modulus exponent of that partition.
- **The map.** `T_sigma(x) = x/2` (x even) and `(7x + sigma(x))/2` (x odd). `rho_max(sigma)` is the largest odd density of a periodic orbit. By Lemma C of THM-4508 this is the densest cycle of the parity graph at any level `>= depth`.
- **Max-halving (MH).** MH is `s(x) = x mod 4` (in `{±1}`). Its accelerated step is `x -> x_1 = (7x + s_1)/2^(v_1)` with `v_1 >= 2`.
- **Itineraries.** The *itinerary* of `x` is `((s_1,v_1),(s_2,v_2),...)`. Lemma MH(ii): a prefix of `n` symbols is one class mod `2^(1 + v_1 + ... + v_n)`. `S_inf` = {all `v_i = 2`}.
- **Normalized words.** Negation maps the itinerary `(s_i,v_i)` to `(-s_i,v_i)`. So a negation-symmetric rule is described by *normalized words*:
  - the first valuation;
  - then tokens `=v` (same sign as the previous symbol) or `!v` (opposite sign);
  - a token `>=V` means "valuation at least V" and is terminal.
- **Examples.** The gaining pattern of Corollary F is `2 =2`. The alternating orbit `±1/11` is `2 !2 !2 ...`.
- A *flip* at `x` uses `-s_1`: `y = (7x - s_1)/2 = 2^(v_1-1) x_1 - s_1`, one odd step of "valuation" 1.

## 2. Fact I: what itinerary coding is, in bits

**Fact I.**
- (a) The cylinder of an exact word `w` is the class mod `2^(1+|w|_v)`, and the cylinder of `w(s, >=m)` is the class mod `2^(|w|_v + m)` (Lemma MH(ii)).
- (b) Conversely, every odd class `c mod 2^d` is one of these. Let `w` be the symbols determined by the class. If `d = 1 + |w|_v`, the class is `[w]`. Otherwise the next sign `s` is known and `7 x_n + s = 0 mod 2^(d - |w|_v)`, since a nonzero residue would determine one more symbol. Then the class is `[w(s, >= d - |w|_v)]`.
- So the binary trie of residue classes *is* the itinerary tree, with `>= V` nodes as binary chains. A finite rule has the same leaves in both codings, and "itinerary tree size" equals "bit tree size".
- (c) A rule that reads the *past* of the orbit is not a function of the point (points have two or three `T_sigma`-preimages). It is therefore not a sign strategy in the sense of THM-4474. At a fixed level, memory cannot beat positional play anyway (Theorem G: positional determinacy).

*Proof.* (a) is Lemma MH(ii) with the last symbol only partly fixed. (b) runs the MH steps on the known bits: a symbol is determined exactly when the residue `7 x_n + s mod 2^m` is nonzero, and `x_(n+1)` is then known mod `2^(m-v) >= 2`. ∎

*Checked* (runner M0b): all 4095 odd classes with `d <= 12` re-parse to the same class.

**What is new** in the itinerary angle is therefore not the finite rules but:
- (i) *regular* (cyclic) itinerary automata, which read unboundedly many symbols and are truncated at a bit depth;
- (ii) evaluation at *variable depth* (Lemma M);
- (iii) the dynamics of flips in itinerary coordinates (§4).

## 3. Lemma M: exact rho_max at variable depth (Markov refinement)

**Lemma M.**
- *Hypotheses.* Let `P` be a finite partition of `Z_2` into residue classes (leaves), with `sigma` constant on the odd leaves. Assume that for every leaf `L = (c mod 2^d)` (`d >= 1`) the image `T(L)` is a union of leaves. `T(L)` is the class `T(c) mod 2^(d-1)`, because `T` is affine on `L` with slope `1/2` or `7/2`.
- *The graph.* Let `G` have an edge `L -> L'` iff `L'` is contained in `T(L)`.
- *Conclusions.*
  - (a) Every periodic orbit of `T_sigma` gives a closed walk of `G` of the same length and odd count.
  - (b) Every closed walk `L_0 -> ... -> L_(n-1) -> L_0` comes from exactly one periodic point.
  - (c) `rho_max(sigma)` is the maximum mean cycle of `G` with weight 1 on odd leaves.
  - (d) Every finite partition on whose odd leaves `sigma` is constant has such a refinement ("Markov refinement"), of depth at most its own.

*Proof.*
- (a) `x_(i+1) = T(x_i)` lies in `T(L_i)`, a union of leaves, so its leaf `L_(i+1)` is contained in `T(L_i)`.
- (b) Put `K_n = {x in L_0 : T^i x in L_i, i <= n}`.
  - By induction `T^n(K_n) = L_n`: `T^n` maps `K_(n-1)` bijectively (affinely) onto `T(L_(n-1))`, which contains `L_n`.
  - For a closed walk, `T^n : K_n -> L_0` is an affine bijection that multiplies 2-adic distances by `2^n`. Its inverse maps the closed ball `L_0` into `K_n` (a subset of `L_0`) and contracts by `2^(-n)`, so it has exactly one fixed point.
  - Every periodic point realizing the walk is such a fixed point.
- (c) This follows from (a) and (b).
- (d) While some image `T(L)` lies strictly inside a coarser leaf `L'`, split `L'` along the path to `T(L)`. Every split is at depth `<= depth(L) - 1`, so depths never exceed the initial maximum and the process stops. ∎

*Remarks.*
- The uniform level-`k` parity graph of Lemma C is the special case `P` = all classes mod `2^k`.
- Without the Markov condition, the over-approximated graph (`L -> L'` whenever `L'` meets `T(L)`) still gives a sound *upper* certificate. This was used by the CEGAR attempt (§5.6).

**The evaluator** (`procgen_seven3_20260926_markov.c`, and an independent Python version in the lib). It computes the refinement and then:
- runs Howard policy iteration for a candidate densest cycle;
- verifies an **exact** integer potential on every edge (Lemma G1 on `G`), with Dinkelbach fallback;
- reports the critical set (leaves on cycles of maximal density);
- re-walks the witness as a rational periodic point with `Fraction` arithmetic.

*Checks* (runner M, F, L; the two implementations always agree):

| rule | depth | Markov leaves | uniform level-depth nodes | rho_max | witness |
|---|---|---|---|---|---|
| MH | 2 | 3 | 4 | 1/2 | 2/3 |
| 18-class 2/5 rule of THM-4508 | 10 | 170 | 1024 | 2/5 (also seven2 engine at k = 10, 12, 14) | -12/17 |
| greedy optimal tree, k = 14 (300 leaves) | 14 | 1151 | 16384 | 15/38 | |
| greedy optimal tree, k = 16 (508 leaves) | 16 | 2595 | 65536 | 7/18 | |
| itinerary-search tree (783 normalized leaves) | **46** | 11181 | 2^46 | 7/18 | period 54 |

## 4. Flips in itinerary coordinates: rigidity, orbit switching, shadowing

**Lemma R (rigidity of rejoining).**
- *Setting.* A branch word `B` is a finite word of steps `x -> x/2` and `x -> (7x + s)/2`, with `a` odd steps out of `p`. It acts by the affine map `g_B(x) = (7^a x + c_B)/2^p`.
- *Statement.* There is a countable set `X_exc` of rationals such that, for `x` outside `X_exc`, `g_B(x) = g_B'(x)` implies `(a, p) = (a', p')` and `c_B = c_B'`.
- *In particular.* For such `x`, suppose the orbit of `x` under any strategy (for example one flip followed by MH) meets the MH orbit of `x`. Then it does so after the same number of odd steps and the same number of halvings as MH.

*Proof.*
- If `(a,p) ≠ (a',p')`, then `7^a/2^p ≠ 7^(a')/2^(p')` by unique factorization. So `g_B - g_B'` is a non-constant affine map with rational coefficients, and it vanishes at one rational point only.
- Let `X_exc` be the set of these points over all pairs of branch words. It is countable.
- If `(a,p) = (a',p')`, then `g_B(x) = g_B'(x)` forces `c_B = c_B'`. ∎

**Consequence for the flip calculus.**
- Corollary F of THM-4508 measures the gain of a flip over two odd steps.
- By Lemma R, a flipped orbit that rejoins pays that gain back exactly. On the lemma's pattern `(s,4),(-s,·)` the rejoin is immediate and the gain is 0 at once.
- A flip changes densities only through orbits that never rejoin, or through exceptional (rational) coincidences. Periodic orbits consist of such coincidences.

**Corollary R (orbit switching is necessary).** Let `sigma` be a finite rule with `rho_max(sigma) < 1/2`. Then every `x` in `S_inf` outside the countable set `X_exc` meets its own MH orbit along its `T_sigma`-orbit at only finitely many times.

*Proof.*
- Suppose `T_sigma^(n_i)(x) = x_(j_i)` for infinitely many `i`.
- By Lemma R the two paths have equal counts. On `S_inf` every MH symbol has `v = 2`, so the first `n_i` steps of the `T_sigma`-orbit contain exactly `n_i/2` odd steps.
- A weak* limit of the empirical measures along `n_i` is then `T_sigma`-invariant, with odd mass `1/2`. The odd set is clopen, so its indicator is continuous.
- For a locally constant `sigma`, the maximum of `mu(odd)` over invariant measures is the maximum cycle mean of the parity graph (Lemma C: the system is a vertex shift). So `rho_max(sigma) >= 1/2`, a contradiction. ∎

The same argument shows the following. If the MH orbit of `x` has a limiting odd density `d_x`, and the `T_sigma`-orbit of `x` rejoins it infinitely often, then `rho_max(sigma) >= d_x`.

**Lemma S (shadowing, q = 7).**
- *Hypothesis.* Let `s` be a sign and `z` odd, with MH itinerary `(-s,v_1), ..., (-s,v_n), (-s,·)`, where all `v_i` lie in `{2,3}`.
- *Conclusion.* `E = (z - s)/2` is odd, its first `n` MH valuations are `v_1, ..., v_n`, and `MH^i(E) = (z_i - s)/2` for `i <= n`.
- *Application.* By Lemma F, a flip on the gaining pattern `(s,2),(s,2),(-s,·)` leads (after the valuation-4 step) to `E = (x_2 - s)/2`. That point shadows the MH orbit `x_3, x_4, ...`, with equal valuations, for as long as the orbit keeps sign `-s` and valuations 2 or 3.

*Proof.* `z - s = -2s mod 4`, so `E` is odd. Write `7z = 2^v z_1 + s`, where `(-s, v)` is the first symbol.
- `v = 2`.
  - `7E = 2 z_1 - 3s`.
  - Since `z_1` is odd, `7E + s' = 2 - 3s + s' (mod 4)`, which forces `s' = s`.
  - Then `7E + s = 2(z_1 - s)`, and `z_1 - s = 2 (mod 4)` because `z_1` has sign `-s`.
  - So `v(E) = 2` and `MH(E) = (z_1 - s)/2`.
- `v = 3`.
  - `7E = 4 z_1 - 3s`, which forces `s' = -s`.
  - Then `7E - s = 4(z_1 - s)`, so `v(E) = 3` and `MH(E) = (z_1 - s)/2`.
- Induct on `n`. ∎

*Checked* (runner R):
- **S1:** Lemma S on 3000 random cylinders.
- **R1:** 1500 random 1000-bit `x` with one flip and MH afterwards. Every rejoin (equality mod `2^200`) used equal odd counts and halvings.
- **R2:** at 200 random `S_inf` points (and 3000 more in an exploratory run), one flip followed by MH **never** met the MH orbit again. The valuations afterwards were typical: mean 3.00, with frequencies `1/2, 1/4, 1/8, ...`.
  - The defect `w_n = 2^(j_n) z_n + beta_n` (exploratory) has halving deficit `j_n -> -infinity` and more than 6·10^5 distinct states `(j, beta)`.
  - So there is no finite "defect transducer" on `S_inf`. A single flip turns an `S_inf` orbit into a typical orbit of density about `1/3 < log_7 2`.
  - The limit is therefore decided on the countable exceptional set, not by typical behaviour.

## 5. Task A: the upper side

### 5.1 The simplest itinerary family

Let `G` be the gaining pattern `2 =2`, and `A_D` the alternating run `2 !2 ... !2` of length `D`. The rule `G + A_D` is MH, flipped on the union of the two.

**Proposition A.**
- For every `D >= 2`, `rho_max(G + A_D) >= D/(2D+1)` (PROVED).
- Equality holds for `D = 3..15` (FINITE-EXACT): both Markov engines, and the seven2 uniform engine for `2D + 1 <= 19`. For `D = 2` the value is `1/2`: flipping all of `2 !2` is too coarse.

*Proof of the bound.*
- Take the MH-periodic orbit whose itinerary alternates signs throughout and has valuations `2^(D-1) 3` repeated. Its period is `D` symbols for even `D` and `2D` for odd `D`, so that the alternation closes. It exists and is unique by Lemma MH(ii).
- No point of it starts with two equal signs of valuation 2, nor with `D` alternating valuation-2 symbols. So it is never flipped.
- Its density is `D/(2(D-1)+3) = D/(2D+1)`.
- The runner builds it as a rational, re-walks it under the rule, and checks the value. ∎

So cutting the alternating orbit alone never gets below `3/7` (`D = 3`). Orbits that mix alternating runs with valuation-3 steps must also be broken. The 18-class 2/5 rule of THM-4508 does that: in normalized form its flips are `2 =2`, `2 !2 !2 =*`, `2 !2 =2 !*`, `2 =3 =2`, `2 !2 !2 !2 !*`, `2 !2 !3 =2`, `2 !2 !2 !3`, `2 !2 =3 !2`, `2 !2 !3 !2` (`*` = any valuation).

### 5.2 Exhaustive pattern sets and automata (FINITE-EXACT, runner P1, U1)

| rule class | size | least rho_max |
|---|---|---|
| complete normalized trees of length 2 over valuations `{2,3,>=4}`, every flip set | 13 leaves, 8192 sets | 1/2 |
| complete trees of length 3 over `{2,>=3}`, every flip set | 11 leaves, 2048 sets | 3/7 |
| automata: start + 1 transient state over `{=2,!2,>=3}`, truncated at bit depth 16, default keep or flip | 216 x 2 | 8/17 (the automaton is `G + A_8`) |
| start + 2 transient states over `{=2,!2,>=3}` | 32768 x 2 | 3/7 |
| start + 1 transient state over `{2,3,>=4}` | 5832 x 2 | 8/17 |
| 2 or 3 transient states over `{2,3,>=4}`, 12000 random each (exploratory) | | 8/17, 1/2 |

- *Automaton model.* The start state reads the first valuation class. Transient states read (relative sign, valuation class). Terminal classes must go to FLIP or KEEP. Classes still undecided at bit depth 16 get the default.
- *Reading.* Small itinerary automata stop at `3/7`. Reaching `2/5` already needs valuation-3 context in specific positions, i.e. a tree, not a small automaton. This mirrors the bit-automaton result of THM-4508 §6 (never below 1/2 there).

### 5.3 Counterexample-guided search over itinerary trees

*Method.*
- The rule is a normalized itinerary tree (`procgen_seven3_20260926_itree.py`).
- Evaluate it (Lemma M) and take its witness periodic orbit. For each odd point of the orbit, try two moves:
  - toggle the leaf containing the point;
  - refine that leaf into its itinerary children (`=2, =3, =>=4, !2, !3, !>=4`) and toggle the child containing the point.
- Move to the best non-tabu candidate, ordered by `rho_max`, then the size of the critical set, then the size of the tree.
- No level-`k` arena is used, so the depth is unbounded.

*Results* (the stored trees are re-evaluated by both Markov engines, runner L):
- From MH (400 moves, deterministic; the runner re-runs it) the best value falls to `2/5` at move 20 (48 normalized leaves). It then passes 19/48, 17/43, `15/38` (move 140, 344 leaves), 13/33, 11/28, 9/23 and 16/41, and reaches `7/18` at move 346 (783 normalized leaves).
- The `7/18` tree has 1566 residue classes of depths 4..46, mostly 13..23. Its Markov refinement has 11181 leaves.
- It equals `rho*(7,16)` exactly, so it is optimal among rules of depth `<= 17`. Being deeper, it is not a new bound.
- Its flips cover 17.8% of the odd numbers, all with `v_1 = 2`: 10.2% on the gaining pattern `2 =2` (82% of that pattern), 4.5% on `2 !2`, 2.2% on `2 =3`, 0.9% on `2 !3` (runner L6). This matches Corollary F and the exception statistics of THM-4508.
- **It stalls at 7/18** (EMPIRICAL):
  - Continuing from the `7/18` tree for 350 more moves stayed at `7/18`. The critical set fluctuated between 1700 and 2700 leaves while the tree grew to 1467 leaves.
  - Re-solving Min's game on the tree's own Markov partition (11181 leaves, free sign per leaf, over-approximated successors) fails at `47/121`, `26/67` and `19/49`: Min loses from every leaf.
  - Below `7/18` the partition must be refined, which is the wall met by CEGAR (§5.6). Single-cycle moves cannot make the coordinated change that the level-18 optimum makes.

### 5.4 Optimal rules are not sparse (FINITE-EXACT, runner O1)

For each level, the optimal strategy is MH wherever the least Min potential admits it, and the other sign elsewhere. Its minimal bit tree (= minimal itinerary tree, Fact I) and its Markov refinement:

| k | rho*(7,k) | leaves | flipped leaves | Markov leaves | share of `2^k` |
|---|---|---|---|---|---|
| 12 | 2/5 | 270 | 126 | 883 | 21.6% |
| 14 | 15/38 | 1268 | 574 | 4085 | 24.9% |
| 16 | 7/18 | 3480 | 1660 | 12593 | 19.2% |
| 18 | 19/49 | 11972 | 5810 | 44879 | 17.1% |
| 19 | 13/34 | 22894 | 11264 | 87775 | 16.7% |
| 20 | 34/89 | 42498 | 20624 | 163471 | 15.6% |

- In every row the Markov engine returns exactly `rho*(7,k)`. This is a second, independent confirmation of those values.
- Greedy trees (THM-4508's construction, own driver; runner G1) are much smaller: 46 / 300 / 508 leaves with Markov refinements 170 / 1151 / 2595 at k = 10 / 14 / 16. But they also grow roughly like `2^(0.6k)`.
- Evaluating is cheap at depth 46 (Lemma M). *Finding* rules below `37/100` needs depth above 30, and the only systematic constructions known are the level-`k` games, which are infeasible beyond `k = 31`.

### 5.5 Is the itinerary the natural coordinate of the optima? (EMPIRICAL, runner I1)

*Setup.*
- Classify each odd residue at level `k` by the least Min potential at `rho*`: forced-MH, forced-flip or free. This classification is canonical.
- Call a cylinder *conflicting* if it contains both forced-MH and forced-flip residues.
- Compare the share of residues in conflict-free cylinders, for itinerary cylinders of `n` symbols and for bit cylinders with about as many cells.

| k | rho* | free / forced-MH / forced-flip | itinerary n = 3 vs bits | n = 4 vs bits | n = 5 vs bits |
|---|---|---|---|---|---|
| 14 | 15/38 | 0.065 / 0.754 / 0.181 | 416 cells: 0.567 vs b = 10 (512 cells): 0.676 | 1392: 0.776 vs b = 12 (2048): 0.870 | 2288: 0.925 vs b = 13 (4096): 0.949 |
| 16 | 7/18 | 0.086 / 0.748 / 0.166 | 516: 0.547 vs b = 11 (1024): 0.727 | 2676: 0.733 vs b = 13 (4096): 0.872 | 6628: 0.894 vs b = 14 (8192): 0.947 |
| 18 | 19/49 | 0.076 / 0.761 / 0.163 | 552: 0.547 vs b = 11 (1024): 0.621 | 4032: 0.648 vs b = 13 (4096): 0.745 | 15024: 0.804 vs b = 15 (16384): 0.876 |

- The table gives the share of odd residues lying in conflict-free cells. Itinerary cells lump valuations `>= 6`; the bit depth `b` is the least with at least as many cells.
- Exploratory runs at k = 12 and k = 20 give the same picture. For example, at k = 20 and n = 5: 26948 itinerary cells 0.761, against 32768 bit cells 0.867.

*Reading.*
- At equal numbers of cells, bit cylinders are **at least as consistent** as itinerary cylinders, and usually more. The optimal decisions are not organized by short itinerary windows.
- The forced-flip set covers 15–18% of the odd residues. The free set is 4–18% (k = 12..20, exploratory for k = 12, 20).

### 5.6 CEGAR on abstract arenas (negative, exploratory)

*Method.* Min's game on a variable partition with over-approximated successors (§3 remark). The partition is refined along counterexample cycles of Max's strategy.

*Result.* Started from the 18-class 2/5 rule with target `15/38`, it ran 400 rounds and reached 2007 leaves. Min still lost from every leaf, while the known 15/38 rule has a Markov partition of 1151 leaves. Refining one counterexample at a time does not find it.

### 5.7 Best upper rules (summary)

| value | rule | depth | description size | status |
|---|---|---|---|---|
| 2/5 | 18 flipped classes (THM-4508), 9 normalized patterns | 10 | 170 Markov leaves | FINITE-EXACT |
| 15/38 | greedy tree k = 14 | 14 | 300 leaves | FINITE-EXACT |
| 7/18 | itinerary search tree | 46 | 783 normalized leaves | FINITE-EXACT |
| 34/89 | optimal rule k = 20 (MH-preferring) | 20 | 42498 leaves | FINITE-EXACT |
| 37/100 | THM-4508 (level 30) | 30 | not explicit here | FINITE-EXACT (THM-4508) |

**Nothing below `37/100`** is claimed. The best rule found in this lane is the k = 20 optimum, at `34/89`.

## 6. Task B: the lower side (adversaries)

*Engine* (`procgen_seven3_20260926_adv.c`). For a fixed lift strategy `tau` it computes, exactly:
- **Value.** The value `F` is the largest `F` with a nonempty set `W`, closed under Min's moves, on which every cycle has density `>= F`. It is certified by Lemma G2: `W` is checked closed, and the least-fixed-point potential is checked on every edge.
- **Maximality.** From *every* node, Min can reach a cycle of density `<= F`. This is certified by the least fixed point of `h(s) = max(0, e(s) + min over t of h(t))`, finite everywhere.
- **Witness.** A cycle of density `F`, re-walked in Python.

| adversary family (uniform in k) | values | status |
|---|---|---|
| top lift = negative integers (Theorem N) | exactly 1/3 at k = 8..18 | FINITE-EXACT |
| d-scaled top lift: choose the lift t with `d t mod 2^k` in the top half (follows `-c/d`) | d = 1, 3, 5: 1/3. d = 11: 9/29, 11/36, 12/43. d = 13: 8/25, 11/39, 3/10. d = 17: 9/28, 20/63, 14/45. d = 23: 1/4. All at k = 10, 14, 18, all `<= 1/3` | FINITE-EXACT |
| MH-density greedy lifts (prefer the lift whose determined MH itinerary is denser) | `<= 1/3` at k = 12..16 (e.g. 7/26, 13/44 at k = 16); 5/17 to 1/3 at k = 8..11 | FINITE-EXACT |
| itinerary-prefix lifts `tau(P) = f(first n MH symbols of P)`, `f` hill-climbed at k = 10 | n = 1, 3: 1/3. n = 2: 14/41 at k = 10, then 1/5, 11/41, 12/49 at k = 12, 14, 16 | FINITE-EXACT |

*Reading.*
- **Low part of the window.** The lift decides the *top* bit, which is the far end of the itinerary window. Rules that read the itinerary prefix (the low bits) react to the present by setting the distant future. They fall back to `<= 1/3`, as the low-bit adversaries of THM-4508 do.
- **Top part of the window.** Rules that read the top part follow real (rational) dynamics. There Min has free sign choices along orbits, and reaches a density-`<= 1/3` cycle: `(1,4,2)` for `d = 1`, `1 -> 5 -> 16 -> 8 -> 4 -> 2 -> 1` of `7c ± 3` for `d = 3`, and cycles below 1/3 for larger `d`. Note that `-1/9` lies on a density-1/4 cycle.
- **What a floor above 1/3 needs.** Any all-level certificate above `1/3` must exploit local constancy: it must pair points that are congruent to high order and need opposite signs (Proposition S of THM-4486). By Corollary R and §4 R2, it cannot come from typical orbits.
- **No uniform proof.** No inductive certificate was found. The all-level floor remains `1/3` (Theorem N, PROVED). The finite-level floors remain THM-4508's (`rho*(7,31) >= 7/19`).

## 7. Task C: the limit object

**7.1 Formulation (PROVED).** Let `LC` be the set of locally constant sign functions on the odd 2-adic integers.
- For `sigma` in `LC`, `T_sigma` is conjugate to the vertex shift of its parity graph at level `depth(sigma)` (Lemma C). So `rho_max(sigma) = max over T_sigma-invariant probability measures mu of mu(odd)`, attained on a periodic orbit.
- Every level-k strategy is LC, and every LC `sigma` is a level-`depth(sigma)` strategy. Since `rho*(7,k)` is non-increasing,
  `rho_inf(7) := lim_k rho*(7,k) = inf over sigma in LC of max over mu in M(T_sigma) of mu(odd)`.
  This is an ergodic-optimization min-max: Min picks a locally constant sign function, Max picks an invariant measure, or equivalently a periodic orbit.
- **Itinerary form.** Off the countable set of points whose MH orbit reaches `±1/7` (and then 0), the MH itinerary identifies the odd 2-adic integers with `Omega = ({±1} x {2,3,...})^N`, and MH becomes the shift.
  - An LC rule is a clopen flip set `F`. On odd points, the accelerated `T_sigma` is the shift off `F` and the flip map `Phi(x) = 2^(v_1-1) x_1 - s_1` on `F`.
  - A periodic orbit has density `(number of accelerated steps)/(sum of weights)`, with weight `v_1` for a shift step and 1 for a flip.
  - So `rho_inf(7) = inf over clopen F of the densest periodic orbit of this shift-or-flip system`.

**7.2 What would decide `rho_inf(7) < log_7 2` versus `>= log_7 2`.**
- **`<`** is certified by **one** finite object: a clopen `F` whose Markov refinement (Lemma M) carries an integer potential at some `a/p` with `2^p > 7^a`. Equivalently, some level `k` at which Min's least fixed point at `log_7 2` is finite.
  - This is semi-decidable (enumerate `k`), and nothing below 37/100 is known.
  - Necessary features (PROVED): by Corollary R, `F` must switch every generic `S_inf` orbit off its MH orbit permanently.
  - Features seen in the data (EMPIRICAL): optimal rules flip on 15–18% of odd numbers, essentially all with `v_1 = 2`.
- **`>=`** needs a Lemma G2 certificate at **every** level for every `theta < log_7 2`: a uniform family of lift strategies with potentials, as in Theorem N, whose proof is inductive in `k`.
  - Proposition S of THM-4486 forbids a finite list of rational cycles as the obstruction.
  - Lemma R and the single-flip data R2 indicate that typical orbits do not see it either: after flips they have density about 1/3.
  - So the certificate must organize infinitely many exact coincidences (residue collisions).
  - Every uniform family tested here stops at `<= 1/3`.

**7.3 The data (EMPIRICAL).**
- `1/rho*(7,k)` rises: 2.500 (k = 10), 2.533 (14), 2.571 (16), 2.579 (18), 2.615 (19), 2.618 (20), 2.625 (21), 2.643 (22), 2.667 (24), 2.686 (27), 2.692 (28), 2.703 (30). The threshold is `log_2 7 = 2.807`.
- The decrement per level falls from about 0.0018 (k = 10..20) to 0.0012 (k = 20..30).
  - A ratio of about 2/3 per 10 levels would cross the threshold.
  - A ratio of about 1/2 per 10 levels would not.
- The itinerary search descends through the same values as the level optima (2/5, 15/38, 7/18) at larger depth. That is consistent with a slowly decreasing limit, and says nothing about which side.
- Optimal rules have Markov refinements of about 16–25% of `2^k` (§5.4). There is no sign of a finite "renormalized" description.
- **Undecided.**

**7.4 Why q = 5 is different, in itinerary language.**
- **q = 5.** The MH fixed point `1` of `S_inf` (symbol `(-,2)`) flips to `y = 3`, and `MH(3) = 1`. The flipped orbit *rejoins at the same point* after the valuation pair `(1,4)`. This is an exceptional coincidence in the sense of Lemma R (at a rational point).
  - It closes the sporadic cycle `(1,3,8,4,2)`, of density `2/5`.
  - As `(-1,-3,-8,-4,-2)` of `5u - 1`, the same cycle lies on the negative integers. There the top-lift adversary forces it at every level (the `u^2` potential of THM-4486).
  - So Min's first rejoining flip and Max's forced cycle coincide. The limit `2/5` is reached as soon as the residue collisions around this rule are separated (k = 15).
- **q = 7.** The fixed-point flip at `±1/3` closes `(1,5)`, of density 1/3.
  - That is *too good to be forced*. The top lift forces only the free cycle `(1,4,2)`, whose point `-1` has MH symbol `(-,3)` and is not in `S_inf`. Its value is 1/3.
  - The binding obstructions are the neighbours. The pattern `(s,2),(s,2),(-s,·)` gains only 1, and it closes the denominator-17 cycles (density 2/5, THM-4508).
  - Below `2/5` every improvement must switch orbits permanently (Corollary R). The switched orbit first shadows the old one (Lemma S), so the rule must read far ahead.
  - The optimal decisions are then global: they are not captured by short itinerary windows (§5.5), and the rules grow without bound (§5.4).
- **In one line.** For q = 5 the first exceptional closure is at once Min's optimum and Max's forced cycle. For q = 7 it is too good to be forced, and the value is set by an unbounded hierarchy of exceptional closures.

## 8. Task D: the arena, reproduced (FINITE-EXACT, runner D)

- **Engine.** This lane's engine (`procgen_seven3_20260926_game.c`) was written without reading the seven or floor engines. It computes the least fixed points of Min's and Max's operators with a FIFO worklist.
- **Search.** A Stern–Brocot search with caps `64·F.denominator` ran from scratch for k = 8..16. Caps can only make Min lose wrongly, and the value is accepted only when both certificates hold.
- **Direct check.** For k = 17..20 the least fixed points were computed directly at the published values.
- **Certificates.** Both certificates are checked with exact integers at every level.

| k | 8 | 9 | 10–13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| rho*(7,k) | 3/7 | 3/7 | 2/5 | 15/38 | 15/38 | 7/18 | 7/18 | 19/49 | 13/34 | 34/89 |

This agrees with THM-4486 and THM-4508. It also agrees with §5.4, where the Markov engine evaluates the k = 12..20 optimal rules to exactly these values.

## 9. Failures and caught mistakes

- **Caps.** The first from-scratch search used caps `N·F.denominator`. Divergence detection then cost a factor 4 per level. It was replaced by small caps; this is sound, because the value needs both certificates.
- **Inverted adversary.** The first d-scaled adversary run had the lift inverted: it followed positive rationals `c/d` and got values of about 0. This was corrected before any use.
- **Stuck wait loop.** A shell wait loop `until ! pgrep -f name` matched its own command line and hung until its timeout. Only scratch was affected.
- **Leaf counting.** Itinerary trees are counted in normalized leaves; each is two residue classes. The 783-leaf tree has 1566 classes.
- **Not in the runner (EMPIRICAL, scratch only):**
  - the CEGAR prototype, abandoned after 400 rounds (§5.6);
  - the 350-move continuation and the re-solve on the tree's own partition (§5.3);
  - the random automata samples (§5.2);
  - the 3000-trial defect-state exploration (§4).

## 10. Reproduction

- **Command.** `python3 04-computation/experiments/procgen_seven3_20260926_run.py > 05-knowledge/results/procgen_seven3_20260926.out`.
- **What it re-checks.** Every FINITE-EXACT claim above (sections C, M, F, L, G, O, I, R, B, P/U, D), ending with `ALL CHECKS PASSED`.
- **Resources.** 94 checks, 380 s wall time on the shared 8 GB machine; the runner process peaked at 152 MB RSS, and every C engine run stays far below that (the largest array is `2^20` int32 at k = 20). There is one process at a time, under `nice`.
- **Build.** Engines are compiled into `scratch/procgen_seven3/bin` by the lib, with source hashes in the file names.
- **Data.** `procgen_seven3_20260926_rules.json` holds the itinerary trees, the greedy trees and the itinerary-prefix adversary. They are re-evaluated, not trusted.
- **Output hash.** sha256 of the .out (raw bytes, includes timing fields): `5ab858af8f684469c133d7555a62f336134911df91add8646bc376b1de482b61`.

