# Reverse-tree pieces: the three-row automaton, the fibre-uniformity lemma, and the exchange of directions

**Status:** PROVED (scoped, elementary): Theorem 1.1 (row of the `j`-th child is `rho(2^(h0+1+2j) u mod 9)`, base row from `u mod 18`, rotation `1 -> 5 -> 3` of exact period `ord_9(4) = 3`), Theorem 1.2 (rows of the first `D` descendants along a fixed index path are a function of `u mod 3^(D+1)`), Lemma 2.1 (fibre uniformity: the infinite fibre of every internal node is uniform mod 9), Proposition 3.1 (branching factor `4/3` of the pruned `T`-tree and the exact `(4/3)^k` expectation), Proposition 3.2 (a split-independent modulus-9 difference inequality with rigorous exponent `0.2462`), Proposition 4.1 (Collatz `<=>` surjectivity of the guarded `D/E` closure of `{1}` onto the odd non-multiples of 3), Proposition 4.2 (the exact exchange of directions between Collatz and the E-graph question Q2), Proposition 5.1 (the minimal-child map `m` is injective with image the odd `n != 5 mod 8`, fixed point `1` only), Proposition 5.2 (survival law: exactly `(2/3)^L` of the classes `u mod 3^(L+1)` survive `L` minimal-child steps). FINITE-EXACT: tree to depth 6 with `j <= 8` (69,984 nodes), sharpness of Theorem 1.2 for `D <= 4`, depth census of all odd `u <= 10^6`, `(4/3)^k` for `k <= 7`, depth counts of the pruned tree to depth 32, `m`-orbits of all `333,333` odd non-multiples of 3 below `10^6`, survival classes for `L <= 8`. CITED: Krasikov-Lagarias (Acta Arith. 109, 2003, exponent 0.84), Applegate-Lagarias (Math. Comp. 64, 1995, exponents 0.643 and 0.81), Lagarias-Weiss (Ann. Appl. Probab. 2, 1992, branching model); UNCITED-RECOLLECTION: Krasikov (1989) exponent 3/7 and all page numbers. HEURISTIC: the adversarial constant-split modulus-9 system, whose vertex minimum `0.4366` is a ceiling on that scheme, not a theorem. OPEN: Collatz; Q2; whether every `m`-orbit ends at a leaf. No novelty claim for the inverse fibre, the guarded word model, the `4/3` branching heuristic, or the density-bound literature. SCOPE: no map found to THM-4139/THM-4146/THM-3341. Results note of session collatz-mod6-20260917 (machine mac-mini), lane `reverse_tree_pieces`, written 2026-09-21; not a reserved canon ID.

## Inheritance and concept board

Object and coordinates come from [the Collatz braid note](arithmetic_braids_20260917_collatz.md): rows `1, 3, 5 mod 6`, the inverse fibre `n_j = (2^(h0+1) 4^j u - 1)/3` (its (B1)), the braid `R(n) = 4n+1` (its (B2)), and the entry rules (`u = 1 mod 6` gives sources `1 mod 8`, `u = 5 mod 6` gives sources `3 mod 4`). The guarded word model `D(x) = 2x`, `E(x) = (2x-1)/3` with guard `x = 2 mod 3` is [the blueprint audit's affine note](collatz_blueprint_20260921_affine.md) (section 1, "guarded words `E o D^(k-1)`"). The wave-one `row_braid` lane proved `row(u,h) = rho_F(2^h u mod 9)` in the `F = (3n+1)/2` convention (its `.out`, section "row law", and the scratchpad typing JSON `w1_row_braid_typing.json`); this note restates it in the `3n+1 = 2^(h+1) u` convention of the task and builds the automaton on it. The wave-one `E-graph` lane ([note](collatz_mod6_20260917_extended_collatz_scc.md), `.out` S2-S4) supplies the greedy map `G`, its mod-9 Markov chain, the stationary law `2/9` on `{1,4,7}` and `1/9` on `{2,5,8}`, the sharp `m mod 3^(J+1)` dependence of the first `J` greedy letters, and the outsider lists; none of that is re-proved here. **Closest proved mechanism:** the row law `row = rho(3n+1 mod 9)` is a one-line residue identity (`3(6i+r)+1 = 18i + 3r + 1`), and everything in section 1 is bookkeeping on top of `ord_9(4) = 3` and `ord_27(4) = 9`. **Canonical hostile:** the loop `child(1,0) = 1`: a tree generator that re-expands it duplicates `5` at level 2 (found and fixed while writing the script; the loop edge is dropped from the frontier). **Corrected near miss:** the draft claim "`G` and the minimal odd child coincide iff `u = 1, 2 mod 9`" was refuted by `u = 17` (`8 mod 9`, both give `k = 1`); the true rule is `u = 1, 2, 8 mod 9`. A second near miss: the `m`-orbit length law was first compared with `(1/3)(2/3)^L`; the truth is `(1/3)(2/3)^(L-1)` because the leaf step is counted. **Least-used sidecar:** the phase of the 9-periodic residue cycle inside a fibre (the residue of `child(u,0) mod 9`, a function of `u mod 27`), which is what a truncated tree sees and the infinite tree does not. **Typed analogy (minimal-child orbit -> greedy `G`-orbit):** source = `u, m(u), m^2(u), ...` in the odd tree; target = `m, G(m), G^2(m), ...` in `E`; map = both choose the minimal exponent `k` (`m`: minimal odd-admissible `k`; `G`: minimal `k` with an internal image, any parity); preserved = the one-step residue law driven by `u mod 27 -> image mod 9`, and coincidence of the two maps exactly on `u = 1, 2, 8 mod 9`; lost = `G` never dies and empirically reaches `1`, while `m` dies at a leaf with density one and never returns to `1`; sidecar = the step-count and peak laws; test = survival `(2/3)^L` (this note) against the stopping-time bound `(7/9)^(J-1)` of the E-graph lane. No map found (SCOPE) to [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md), [THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md), or [THM-3341](../../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md): only the mod-6/mod-9 row typing is shared.

Notation. `T(n)` = odd core of `3n+1` for odd `n`. For odd `u` with `3 !| u`: `h0 = 1` if `u = 1 mod 3`, `h0 = 0` if `u = 2 mod 3`; `child(u,j) = n_j = (2^(h0+1+2j) u - 1)/3`, `j >= 0`, the complete list of odd `n` with `T(n) = u` (inherited (B1)). `row(n) = n mod 6 in {1,3,5}`. `rho(4) = 1`, `rho(1) = 3`, `rho(7) = 5`. The tree is rooted at `1`, with the loop edge `child(1,0) = 1` dropped.

## 1. The three-row automaton (PROVED, FINITE-EXACT)

**Theorem 1.1.** For odd `n`, `row(n) = rho(3n+1 mod 9)`. Consequently, for internal `u` and `j >= 0`,

    row(child(u,j)) = rho(2^(h0+1+2j) u mod 9) = base(u) + 4j mod 6,

where `base(u)` is the base row of the table below; the rows of successive children rotate `1 -> 5 -> 3 -> 1` with exact period `3 = ord_9(4)`, so exactly every third child is a row-3 leaf.

*Proof.* Write `n = 6i + r`, `r in {1,3,5}`; then `3n + 1 = 18i + (3r+1)` and `3r + 1 = 4, 10 = 1, 16 = 7 mod 9` for `r = 1, 3, 5`, which is `rho`. For the child, `3 n_j + 1 = 2^(h0+1+2j) u`, and passing from `j` to `j+1` multiplies by `4`, which permutes `4 -> 16 = 7 -> 28 = 1 -> 4` mod 9, i.e. rows `1 -> 5 -> 3 -> 1`. In the `row_braid` convention `F(n) = (3n+1)/2 = 2^h u` the same law reads `rho_F(2^h u mod 9)` with `rho_F(2) = 1, rho_F(5) = 3, rho_F(8) = 5` (multiply by `2^(-1) = 5 mod 9`). QED

**Base row from `u mod 18`** (equivalently `u mod 9` with `u` odd):

| `u mod 18` | `u mod 9` | `h0` | `2^(h0+1) u mod 9` | base row | rows of `j = 0..5` |
|---|---|---|---|---|---|
| 1 | 1 | 1 | 4 | 1 | 1, 5, 3, 1, 5, 3 |
| 5 | 5 | 0 | 1 | 3 | 3, 1, 5, 3, 1, 5 |
| 7 | 7 | 1 | 1 | 3 | 3, 1, 5, 3, 1, 5 |
| 11 | 2 | 0 | 4 | 1 | 1, 5, 3, 1, 5, 3 |
| 13 | 4 | 1 | 7 | 5 | 5, 3, 1, 5, 3, 1 |
| 17 | 8 | 0 | 7 | 5 | 5, 3, 1, 5, 3, 1 |

So `u mod 9 in {1,2}` gives base row 1, `{4,8}` gives base row 5, `{5,7}` gives base row 3 (the minimal child is a leaf). This agrees entry by entry with the `row_braid` table (its row column reads `1,3,3,1,5,5` for `u mod 18 = 1,5,7,11,13,17`).

**FINITE-EXACT.** The tree generated from `1` to depth 6 with `j <= 8` has `1, 8, 45, 270, 1620, 9720, 58320` nodes per level (69,984 in all); every node satisfies `T(n) = parent`, `row(n) = rho(2^(h0+1+2j) u mod 9) = base(u) + 4j mod 6`; the only repeated node is the loop `child(1,0) = 1`, all other nodes are distinct (`T` is a function). Row counts per level are `{1: 2, 3: 3, 5: 3}` at level 1 (nine children of `1`, one of them `1` itself) and then exactly one third per row: `{15,15,15}, {90,90,90}, {540,540,540}, {3240,3240,3240}, {19440,19440,19440}`, because each internal node contributes nine children with rows `base + 4j`, three per row.

**Theorem 1.2 (finite-state only to bounded depth).** Along a fixed index path `(j_1, ..., j_D)`, the rows of the first `D` descendants of `u` are a function of `u mod 3^(D+1)`, and this modulus is sharp.

*Proof.* `n_j = (2^a u - 1)/3` sends `u mod 3^(s+1)` to `n_j mod 3^s` (one exact division by 3), and the row of a node needs the node mod 9; induct on `D`. Sharpness (FINITE-EXACT, `D <= 4`, path `(0,...,0)`): a function of `u mod 9, 27, 81, 243` but not of `u mod 3, 9, 27, 81` (witness classes `1 mod 3^D` in each case; all odd `u < 2 * 3^(D+2)` checked). QED

This is the tree-side twin of the E-graph lane's S4.1 (the first `J` greedy letters are a function of `m mod 3^(J+1)`, sharp). The "three-row automaton" is therefore a finite automaton for one generation only; the row pattern of the `D`-th generation needs one more ternary digit per generation, exactly as the blueprint audit says ("iteration requires further digits").

## 2. Depth census of odd `u <= 10^6` and the fibre-uniformity lemma

All `500,000` odd `u <= 10^6` reach `1` (forward `T`-iteration; each `u` descends to a smaller odd already resolved). Depth = number of `T`-steps to `1`; the maximum is `195` at `u = 837799`.

| depth | count | frac `0 mod 3` | frac `1 mod 3` | frac `2 mod 3` | counts mod 9 (`r = 0..8`) |
|---|---:|---|---|---|---|
| 0 | 1 | 0.0000 | 1.0000 | 0.0000 | 0,1,0,0,0,0,0,0,0 |
| 1 | 9 | 0.3333 | 0.3333 | 0.3333 | 1,1,1,1,1,1,1,1,1 |
| 2 | 34 | 0.3529 | 0.3235 | 0.3235 | 4,4,3,6,5,5,2,2,3 |
| 3 | 78 | 0.3333 | 0.3333 | 0.3333 | 11,9,11,6,5,8,9,12,7 |
| 4 | 176 | 0.3239 | 0.3523 | 0.3239 | 20,23,15,20,22,20,17,17,22 |
| 5 | 282 | 0.3262 | 0.3156 | 0.3582 | 30,29,39,24,27,34,38,33,28 |
| 6 | 508 | 0.3445 | 0.3425 | 0.3130 | 58,49,55,58,63,54,59,62,50 |
| 7 | 871 | 0.3387 | 0.3261 | 0.3352 | 93,94,101,96,92,94,106,98,97 |
| 8 | 1056 | 0.3362 | 0.3400 | 0.3239 | 114,112,123,121,125,110,120,122,109 |
| 9 | 1657 | 0.3355 | 0.3313 | 0.3331 | 197,185,187,170,183,192,189,181,173 |
| 10 | 1860 | 0.3355 | 0.3360 | 0.3285 | 214,201,204,197,205,191,213,219,216 |
| 11 | 2593 | 0.3320 | 0.3301 | 0.3378 | 295,275,278,283,293,283,283,288,315 |
| 12 | 3415 | 0.3327 | 0.3271 | 0.3403 | 381,364,376,374,377,403,381,376,383 |
| 13 | 3514 | 0.3318 | 0.3349 | 0.3332 | 375,387,388,383,387,394,408,403,389 |
| 14 | 4882 | 0.3351 | 0.3359 | 0.3290 | 539,557,535,555,531,557,542,552,514 |
| 15 | 4295 | 0.3376 | 0.3374 | 0.3250 | 502,483,471,480,488,457,468,478,468 |

Bucketed (each entry of the last two columns is a fraction of the bucket):

| depths | count | `0 mod 3` | `1 mod 3` | `2 mod 3` | frac in 1, 4, 7 | frac in 2, 5, 8 |
|---|---:|---|---|---|---|---|
| 0-19 | 52108 | 0.3347 | 0.3336 | 0.3316 | 0.1102, 0.1114, 0.1120 | 0.1116, 0.1103, 0.1098 |
| 20-39 | 168230 | 0.3333 | 0.3334 | 0.3333 | 0.1110, 0.1112, 0.1112 | 0.1113, 0.1107, 0.1114 |
| 40-59 | 144141 | 0.3319 | 0.3338 | 0.3343 | 0.1122, 0.1113, 0.1104 | 0.1107, 0.1121, 0.1115 |
| 60-79 | 98798 | 0.3346 | 0.3321 | 0.3333 | 0.1106, 0.1104, 0.1110 | 0.1110, 0.1113, 0.1110 |
| 80-99 | 28539 | 0.3327 | 0.3332 | 0.3341 | 0.1102, 0.1102, 0.1128 | 0.1131, 0.1095, 0.1114 |
| 100-119 | 6751 | 0.3354 | 0.3360 | 0.3287 | 0.1077, 0.1154, 0.1129 | 0.1092, 0.1096, 0.1099 |
| 120-139 | 1234 | 0.3363 | 0.3485 | 0.3152 | 0.1321, 0.1126, 0.1037 | 0.0956, 0.1175, 0.1021 |
| 140-159 | 168 | 0.3810 | 0.3512 | 0.2679 | 0.1012, 0.1190, 0.1310 | 0.0833, 0.1071, 0.0774 |
| 160-179 | 27 | 0.4815 | 0.2593 | 0.2593 | 0.0741, 0.0741, 0.1111 | 0.1481, 0.0741, 0.0370 |
| 180-199 | 4 | 0.2500 | 0.5000 | 0.2500 | 0.0000, 0.2500, 0.2500 | 0.0000, 0.2500, 0.0000 |

Pooled over all depths the counts mod 9 are `55556, 55556, 55555, 55556, 55555, 55556, 55555, 55556, 55555` (trivially uniform: every odd `u <= X` is a node). Among internal nodes at depth `>= 1` the fraction `1 mod 3` is `0.5000`. The largest deviation of any mod-9 fraction from `1/9` over depths with at least 1000 nodes is `0.0194` (depth 91). The depth-by-depth law is uniform mod 9 to within truncation noise, with the leaves (`0 mod 3`) at one third; it is **not** `G`'s stationary law (`2/9` on each of `1,4,7`, `1/9` on each of `2,5,8`, `0` on multiples of 3, i.e. `2/3` on `1 mod 3`).

**Lemma 2.1 (fibre uniformity, PROVED).** For every internal `u`, the residues `child(u,j) mod 9`, `j = 0..8`, are the nine residues mod 9 in some order; the infinite fibre is uniform mod 9 with period 9 in `j`.

*Proof.* `child(u,j) = (c 4^j - 1)/3` with `c = 2^(h0+1) u = 1 mod 3`. Since `ord_27(4) = 9`, `<4> mod 27 = {1,4,7,10,13,16,19,22,25}` is the full subgroup of residues `1 mod 3`, so `c 4^j mod 27` runs over all residues `1 mod 3` and `(c 4^j - 1)/3 mod 9` over all nine residues, once each per nine consecutive `j`. Checked for `u in {1,5,7,11,13,17,19,23,25,29}`. QED

Hence a third of every fibre are row-3 leaves and one sixth lies in each of the six classes `1,2,4,5,7,8 mod 9`.

**Exact explanation of the discrepancy with `G`.** (i) `G` never visits multiples of 3 (its minimal `k` is chosen so that `2^k m = 4 or 7 mod 9`), whereas the tree counts each leaf once, so one third of the nodes are invisible to `G`. (ii) Among internal children the tree is uniform on the six classes (Lemma 2.1), while `G`'s greedy choice lands in `1 mod 3` from residues `{1,2,4,5}` and in `2 mod 3` only from `{7,8}` (E-graph lane S2.1), which under its stationary law gives `2/3` on `1 mod 3`; the tree gives `1/2`, and the census confirms `0.5000`. (iii) Truncation: the number of children `<= X` of internal `u` is `J_u(X) + 1` with `J_u(X) = floor(log_4((3X+1)/(2^(h0+1) u)))`; the sum over the `333,333` internal `u <= 10^6` is `416,666 = #{odd n <= X : T(n) <= X}` (cross-checked directly), mean truncated fibre `1.2500`. Small `u` contribute the phase-dependent prefix of the 9-periodic residue cycle (the phase is `child(u,0) mod 9`, a function of `u mod 27`), which is the whole source of the finite-depth deviations in the tables; the very deep buckets (`>= 140`) are small-sample.

## 3. Growth: branching factor, exact `(4/3)^k`, and a modulus-9 difference inequality

Here the tree is the pruned `T`-tree of the literature: nodes are integers `n >= 1` with `3 !| n` reaching `1` under `T(x) = x/2` (even) or `(3x+1)/2` (odd); children of `y` are `2y` (always) and `(2y-1)/3` iff `y = 2 or 8 mod 9` (`y = 5 mod 9` gives a multiple of 3, pruned).

**Proposition 3.1 (PROVED; no novelty claim).** The number of children by `y mod 9` is `1, 2, 1, 1, 1, 2` for `y = 1, 2, 4, 5, 7, 8`, mean `4/3` (mean `3/2` if multiples of 3 are kept). The child of a node uniform on the six classes mod `3^(s+1)` is uniform on the six classes mod `3^s` (doubling permutes classes; `(2y-1)/3` sends `9m+2 -> 6m+1` and `9m+8 -> 6m+5`, each uniform in its mod-3 class). Hence the expected number of depth-`k` descendants of a root uniform on the classes mod `3^(k+1)` is exactly `(4/3)^k`.

FINITE-EXACT: summing `#depth-k nodes` over the `2 * 3^k` root classes mod `3^(k+1)` gives `8, 32, 128, 512, 2048, 8192, 32768` for `k = 1..7`, i.e. averages `4/3, 16/9, 64/27, 256/81, 1024/243, 4096/729, 16384/2187`. In the concrete tree from root `8` (the `T`-loop `1 <-> 2` is cut at `4`, whose only preimage is `8`), the node counts at depth `0..32` are

    1, 2, 2, 2, 3, 4, 6, 10, 13, 14, 18, 25, 33, 46, 61, 77, 107, 144, 189, 253, 331, 441, 591, 802, 1066, 1412, 1876, 2492, 3345, 4453, 5936, 7925, 10563

with successive ratios `1.329, 1.325, 1.329, 1.328, 1.342, 1.331, 1.333, 1.335, 1.333` for `k = 24..32` and `10563^(1/32) = 1.3358` against `4/3 = 1.3333`. This is the Applegate-Lagarias / Lagarias-Weiss branching heuristic (CITED below) and nothing more: in size terms each doubling child is twice as large and each odd child two thirds as large, so the tree gains `4/3` nodes per level on average; the classical density exponent is not derived from this expectation alone.

**Proposition 3.2 (modulus-9 difference inequality, PROVED).** Let `N` be the pruned tree, `f_r(x) = #{n in N : n <= x, n = r mod 9}`, `A = f_1 + f_4 + f_7`, `B = f_2 + f_5 + f_8`. Then for all `x >= 1`

    A(x) >= A(x/4) + A(3x/128)/3 + A(3x/64)/3,

and consequently `A(x) >= c x^gamma` with `gamma = 0.246227` the root of `1 = 4^(-g) + (1/3)(3/128)^g + (1/3)(3/64)^g`.

*Proof.* Children are distinct (`T` is a function). Doubling swaps the mod-3 class, and the odd children of `n = 2 mod 9` are `1 mod 3`, those of `n = 8 mod 9` are `2 mod 3`, and are `<= x` iff `n <= (3x+1)/2`; so `A(x) >= B(x/2) + f_2(3x/2)` and `B(x) >= A(x/2) + f_8(3x/2)`. The doubling chains `1 -> 2` (1 step), `7 -> 5 -> 1 -> 2` (3), `4 -> 8 -> 7 -> 5 -> 1 -> 2` (5) give `f_2(y) >= max(f_1(y/2), f_7(y/8), f_4(y/32)) >= A(y/32)/3`, and `4 -> 8` (1), `1 -> 2 -> 4 -> 8` (3), `7 -> 5 -> 1 -> 2 -> 4 -> 8` (5) give `f_8(y) >= A(y/32)/3` (all `f_r` are nondecreasing). Substituting `B(x/2) >= A(x/4) + A(3x/128)/3` yields the displayed inequality. Since `1 in N`, `A(x) >= 1` for `x >= 1`; choosing `c <= 256^(-gamma)` makes `A(y) >= c y^gamma` for `1 <= y < 256`, and for `x >= 256` all three arguments are `>= 1`, so induction on `x` closes with `A(x) >= c x^gamma (4^(-gamma) + (1/3)(3/128)^gamma + (1/3)(3/64)^gamma) = c x^gamma`. QED

The bound `#{n <= x : n reaches 1} >= c x^0.2462` is rigorous and weak: it is the simplest split-independent scheme at modulus 9, obtained by throwing away which of the three lifts mod 27 an odd child lands in.

**Adversarial constant-split system (HEURISTIC).** Keeping the six functions and writing `f_s(x) = sum_{2r = s} f_r(x/2) + sum_{r in {2,8}} alpha_{r->s} f_r(3x/2)` with unknown splits `alpha_{2->.}` on `{1,4,7}` and `alpha_{8->.}` on `{2,5,8}`, the power ansatz `f_r ~ c_r x^g` needs Perron root `1`. The uniform split (the truth, by the uniformity of Proposition 3.1) gives `g = 1.000000`. The nine vertex splits give

| `2 ->` | `8 -> 2` | `8 -> 5` | `8 -> 8` |
|---|---|---|---|
| 1 | 1.244017 | 0.774576 | 1.500000 |
| 4 | 0.867659 | 0.576032 | 1.500000 |
| 7 | 0.602605 | 0.436588 | 1.500000 |

minimum `0.436588` at `2 -> 7, 8 -> 5`. Because a time-varying adversary is governed by a lower joint spectral radius, this vertex minimum is only a ceiling on what a constant-split modulus-9 Krasikov scheme could prove; it is close to but not equal to `3/7 = 0.428571`, and no identification with any published exponent is claimed.

**Literature (CITED).** Krasikov-Lagarias, *Bounds for the 3x+1 problem using difference inequalities*, Acta Arith. 109 (2003): `#{n <= x : n reaches 1} >= x^0.84` for large `x`. Applegate-Lagarias, *Density bounds for the 3x+1 problem*, Math. Comp. 64 (1995), part I (tree-search method, `0.643`) and part II (Krasikov inequalities, `0.81`). Lagarias-Weiss, *The 3x+1 problem: two stochastic models*, Ann. Appl. Probab. 2 (1992), for the branching random walk. UNCITED-RECOLLECTION: Krasikov (1989), exponent `3/7`; the page numbers of all of the above. Those results use moduli `3^k` with large `k` and a nonlinear program; they are not reproduced here and nothing beyond their statements is claimed.

## 4. How the three pieces fit, and the exchange of directions

**Generation.** With `D(x) = 2x` and `E(x) = (2x-1)/3` guarded by `x = 2 mod 3` (blueprint affine note), the odd `T`-preimages of odd `u` are `E(D^h(u))` with `h = h0 + 2j`, since `D^h(u) = 2 mod 3` exactly when `h = h0 mod 2`; `child(u,j) = E(D^(h0+2j)(u))` (checked for `u in {1,5,7,11}`, `j <= 2`), and `j -> j+1` is the braid `R(n) = 4n+1`. The tree is the closure of `{1}` under these guarded words.

**The pieces.** Row 1 (`u = 1 mod 6`, `h0 = 1`) and row 5 (`u = 5 mod 6`, `h0 = 0`) are internal, each with an infinite fibre whose rows rotate `1 -> 5 -> 3` from the base row of section 1. Row 3 (odd multiples of 3) are leaves: `3 !| 3n+1`, so no odd `T`-preimage. A leaf `n = 3 mod 6` belongs to the tree iff its parent `T(n)` does, and `T(n)` is internal.

**Proposition 4.1 (PROVED).** Collatz (every odd `n` reaches `1`) `<=>` every `u = 1 or 5 mod 6` is a node `<=>` every odd `u` with `3 !| u` equals `child(u', j)` for some node `u'` and some `j`, i.e. the guarded `D/E` closure of `{1}` is surjective onto the odd non-multiples of 3. The leaves then follow automatically. So the user's programme "show how the three pieces fit and connect" reduces to: rows 1 and 5 alone must be swept out by the two-row internal tree; row 3 is not an obstruction and never needs to be reached directly.

**Proposition 4.2 (exchange of directions, PROVED).** Let `C` be the deterministic Collatz graph (`n -> n/2` even, `n -> 3n+1` odd) and `E = C` plus the arrows `n -> 3n+1` for even `n` (the E-graph lane's relaxation). Then

- Q1 (Collatz): every `n` reaches `1` in `C` `<=>` the closure of `{1}` under the guarded inverse moves `x -> 2x`, `x -> (x-1)/3` (only for `x = 4 mod 6`) is everything. (Reverse each arrow of `C`; this closure is the Krasikov-Lagarias tree.)
- Q2 (E-graph lane): `1` reaches every `m` with `3 !| m` in `E` `<=>` every such `m` reaches `1` by the unguarded inverse moves `x -> 2x`, `x -> (x-1)/3` (any `x = 1 mod 3`) through non-multiples of 3. The greedy map `G(m) = (2^k m - 1)/3`, `k` minimal with `2^k m in {4,7} mod 9`, selects one such move, and `G^s(m) = 1` is a certificate of `1 ->* m` in `E`.

Thus Collatz is a "TO 1" statement in `C`, equivalently a "FROM 1" statement for the guarded inverse tree, whereas Q2 is a "FROM 1" statement in `E`, equivalently a "TO 1" statement for the unguarded inverse graph; reversing arrows exchanges the two, and the relaxation `C -> E` adds exactly the even `-> 3n+1` arrows, which the E-graph lane identified (its C1) as the fibre index `j = -1`. Q1 and Q2 together are the statement that the non-multiples of 3 form one strongly connected component of `E` (its conjecture `C_E`). Example `m = 7`: the Collatz path `7, 22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1` certifies Q1; the `G`-orbit `7, 2, 1` certifies Q2 through the `E`-path `1 -> 4 -> 2 -> 7`, which uses the even arrow `2 -> 7`.

`G` against the minimal odd child of section 5: on odd `u` they coincide iff `u = 1, 2, 8 mod 9` (both give `k = h0 + 1`; checked for all odd `u < 400`). For `u = 4, 7 mod 9`, `G` takes `k = 0` and its image `(u-1)/3` is even, not a tree node; for `u = 5 mod 9` the minimal odd child is the leaf `(2u-1)/3` and `G` skips to `k = 3`, i.e. `G(u) = child(u,1)`.

**FINITE-EXACT illustration (`X = 10^5`).** The guarded inverse closure of `1` inside `[1, X]` has `39,706` nodes (`26,478` not divisible by 3); the unguarded (reversed-`E`) closure of `1` inside `[1, X]`, multiples of 3 excluded, has `27,472` of the `66,667` non-multiples of 3, and contains the guarded closure minus multiples of 3 (checked). Neither is complete inside `[1, X]` because the forward peak or the inverse peak leaves the box: the smallest non-multiples of 3 missing from the guarded closure are `703, 871, 937, 1055, 1249, 1307, 1406, 1471` (the Collatz peak of `703` is `250,504 > X`), and from the unguarded closure `1535, 2047, 2207, 2287, 2303, 2527, 2687, 2815`, matching the E-graph lane's outsiders at `N = 10^5`. Both questions are FINITE-EXACT to `10^6` in that lane.

## 5. Wildcard: the minimal-child map `m(u) = child(u,0)`

**Proposition 5.1 (PROVED).** (a) `m` is injective: `2^(h0+1) u = 2^(h0'+1) u'` with `u, u'` odd forces `h0 = h0'` and `u = u'`. (b) `image(m) = {odd n : n = 1 mod 8} U {odd n : n = 3 mod 4}`, i.e. the odd `n != 5 mod 8`; the complement `5 mod 8` is exactly `R(odd) = {4n+1}`, the children with `j >= 1`. (For `u = 1 mod 3`, `4u - 1 = 3 mod 8` gives `n = 1 mod 8`, and conversely `n = 1 mod 8` gives `u = (3n+1)/4` odd and `1 mod 3`; for `u = 2 mod 3`, `2u - 1 = 1 mod 4` gives `n = 3 mod 4`, and conversely `u = (3n+1)/2` is odd and `2 mod 3`. Checked for odd `n < 2000`.) These are the inherited entry rules read as an image statement. (c) `T(m(u)) = u`, so `m` moves one level away from the root; `m(u) = (4u-1)/3 > u` for `u = 1 mod 3`, `u > 1`, and `m(u) = (2u-1)/3 < u` for `u = 2 mod 3`. (d) The only fixed point is `u = 1` (`3u = 2^(h0+1) u - 1` forces `h0 = 1` and `u = 1`), and an `m`-cycle would be a `T`-cycle. (e) `m(u) mod 3` is a function of `u mod 9`: `{1,2} -> 1` (row 1), `{4,8} -> 2` (row 5), `{5,7} -> 0` (leaf; the orbit stops).

**FINITE-EXACT.** For all `333,333` odd `u <= 10^6` with `3 !| u`, the `m`-orbit was followed until its first multiple of 3: no orbit reached the step cap `10,000`, no nontrivial cycle exists (impossible below `10^6` anyway, since every such `u` reaches `1` and the only `T`-cycle through `1` is `{1}`), the longest orbit has `L = 29` steps (`u = 774883`), and the largest peak-to-start ratio is `33.26` (`u = 86023`). Orbit lengths `L` (the leaf step counted; `L = 0` only for `u = 1`):

| `L` | count | fraction | `(1/3)(2/3)^(L-1)` |
|---|---:|---|---|
| 1 | 111112 | 0.33334 | 0.33333 |
| 2 | 74074 | 0.22222 | 0.22222 |
| 3 | 49383 | 0.14815 | 0.14815 |
| 4 | 32924 | 0.09877 | 0.09877 |
| 5 | 21949 | 0.06585 | 0.06584 |
| 6 | 14624 | 0.04387 | 0.04390 |
| 7 | 9754 | 0.02926 | 0.02926 |
| 8 | 6498 | 0.01949 | 0.01951 |
| 9 | 4336 | 0.01301 | 0.01301 |
| 10 | 2889 | 0.00867 | 0.00867 |
| 11 | 1931 | 0.00579 | 0.00578 |
| 12 | 1286 | 0.00386 | 0.00385 |
| `>= 13` | 2572 | 0.00772 | `(2/3)^12 = 0.00771` |

The law is geometric with ratio `2/3` to within `2e-4` at every row (checked by the script).

**Proposition 5.2 (survival law, PROVED + FINITE-EXACT `L <= 8`).** The first `L` minimal-child steps stay internal for exactly `(2/3)^L` of the classes `u mod 3^(L+1)` coprime to 3: `4/6, 8/18, 16/54, 32/162, 64/486, 128/1458, 256/4374, 512/13122` for `L = 1..8`. *Proof.* `m` needs only `x mod 3` to fix `h0` and divides by 3 once, so `m(x) mod 3^s` is a function of `x mod 3^(s+1)`, and as `x` runs over the three lifts of a class mod `3^s` to mod `3^(s+1)`, `m(x) mod 3^s` runs over the three lifts of its own mod-3 class (Lemma 2.1's mechanism one digit at a time); hence each step kills exactly one third of the surviving classes. QED

**OPEN.** Whether every `m`-orbit with `u != 1` ends at a leaf. The 3-adic set of `u` with an infinite internal orbit has Haar measure zero (survival `(2/3)^L`); an integer in it would be an odd `u` whose entire minimal-child line avoids `3Z`. This is a Collatz-type 3-adic question; FINITE-EXACT: none below `10^6`.

## 6. Reproduction

    cd /tmp/math-wt-collatz-mod6-b && python3 04-computation/experiments/collatz_mod6_20260917_reverse_tree_pieces.py > 05-knowledge/results/collatz_mod6_20260917_reverse_tree_pieces.out

Runtime about 8 s, memory well under 1 GB; every check is an explicit `raise` and the run ends with `ALL CHECKS PASSED`; identical output under `python3 -O`. Output: [collatz_mod6_20260917_reverse_tree_pieces.out](collatz_mod6_20260917_reverse_tree_pieces.out).

## 7. Stopping boundary / next question

Everything proved here is a residue identity, a counting argument, or an elementary graph equivalence; the two conjectures (Collatz, Q2) and the leaf-termination question for `m` remain OPEN, and the rigorous density exponent from modulus 9 (`0.2462`) is far below the literature's `0.84`, which needs large moduli `3^k`. The decisive structural fact for the user's programme is Proposition 4.1 with Lemma 2.1: the three pieces fit by a one-generation automaton (`base(u mod 18) + 4j mod 6`) whose leaves are a passive third of every fibre, so the whole difficulty is the two-row internal tree, and the residues of its `D`-th generation need `D+1` ternary digits (Theorem 1.2), exactly the digit-per-generation cost of the E-graph lane and of the blueprint's carry. Next question: the E-graph lane bounds the greedy stopping time by `(7/9)^(J-1)`; the minimal-child survival is exactly `(2/3)^L`. Is there a single 3-adic potential on `Z_3^*` that is monotone along both `G` (in `E`) and `m` (in the tree), so that the "FROM 1" and "TO 1" directions of Proposition 4.2 are controlled by one function? A first test is whether the coincidence locus `u = 1, 2, 8 mod 9` of `G` and `m` is closed under either map (it is not closed under `m`: `m(11) = 7`), which already rules out the naive candidate.
