# The merge depth of the Fermat-line pairs for L ≡ 1 (mod 4): a 2-adic tree, not a formula

2026-10-26, klein. **PROVED** (reduction to THM-4581's pair chain, forced prefix, 2-adic locality,
and the exact finite enumeration that decides every class); **FINITE-EXACT** receipts
`tree-normal.json`, `depth-normal.json`, `chain-normal.json` (byte-identical under `-O`). Extends
THM-4793's Theorem M. Setting: `x = 2^{2L+1}+1`, `y = 2^{2L+2}+1`, `T` the Terras map (`x/2` or
`(3x+1)/2`), `S` the Syracuse map; "merge at equal odd index" means `S^k(x) = S^k(y)` with the value
above the trivial cycle.

## 1. The pair is a THM-4581 pair from the start

`y = 2x - 1` is odd, so `T(y) = (3y+1)/2 = 3x - 1`. Hence `(u, v) = (T(y), x)` satisfies
`u = 3^1 v + (-1)`: it is an affinely related pair in the sense of THM-4581 with initial chain state
`(k_0, e_0) = (1, -1)`, and THM-4581's chain `(k_n, e_n)` tracks `T^{n+1}(y) = 3^{k_n} T^n(x) + e_n`
exactly. Absorption `(k, e) = (0, 0)` at Terras time `n` means `T^{n+1}(y) = T^n(x)` with equal
odd-step counts of `x` and `y` (the one extra `T`-step of `y` is the odd step `y -> T(y)`, and `k_0 = 1`
says `u` makes one fewer odd step than `v` from then on). So:

**Lemma A.** `x` and `y` merge at equal odd index iff the chain absorbs before `x` reaches the
trivial cycle, and the merge index is the number of odd steps of `x` up to absorption.

## 2. The forced prefix and the consecutive pair `(3^{L+1}+1, 3^{L+1}+2)`

By THM-4793 (Theorem F's run law), `x` makes `L` valuation-two steps and `y` makes `L` as well,
landing on `a = 2·3^L+1` and `b = 4·3^L+1` after `2L` Terras steps each. One more Terras step of `x`
gives `p := T^{2L+1}(x) = (3a+1)/2 = 3^{L+1}+2`; two more of `y` give `T^{2L+1}(y) = 2(3^{L+1}+1)` and
`T^{2L+2}(y) = 3^{L+1}+1 = p - 1`.

**Lemma B.** For every `L >= 2`, the chain state at Terras time `2L+1` is `(0, -1)`: the two orbits are
the consecutive integers `(p-1, p)` with `p = 3^{L+1} + 2`, and the merge question for the Fermat-line
pair is exactly THM-4581's "`z` and `z+1`" question for `z = 3^{L+1}+1`. (Receipt `chain-normal.json`:
state `(0,-1)` at `2L+1` and `(-1,-2/3)` at `2L+2` for every `L <= 33`.)

## 3. Why even `L` and `L ≡ 3 (mod 4)` are forced and `L ≡ 1 (mod 4)` is not

`3^{L+1} ≡ 3, 9, 11, 1 (mod 16)` for `L ≡ 0, 1, 2, 3 (mod 4)`, so `p ≡ 5, 11, 13, 3 (mod 16)`. The chain
state at `2L+3` is `(-1, -1/3)` for even `L` and `(-1, -2/3)` for odd `L`; it absorbs at Terras time
`2L+4` for even `L` (one odd step after `p`: index `L+2`) and at `2L+5` for `L ≡ 3 (mod 4)` (two odd
steps: index `L+3`), which is Theorem M read off the chain (receipt controls). For `L ≡ 1 (mod 4)`,
`p ≡ 11 (mod 16)` gives `T(p) ≡ 1 (mod 8)` and the chain leaves the forced region.

**Lemma D (three odd steps never suffice when `L ≡ 1 (mod 4)`).** Write `q' = (3^{L+1}+1)/2`
(`y`'s value at odd index `L+1`), so `p = 2q'+1`, and `w = v_2(3q'+1) = v_2(3^{L+2}+5) - 1`. For
`L ≡ 1 (mod 4)`, `3^{L+2} ≡ 11 (mod 16)` gives `w >= 3`. Then:
- `o = 1`: `S(p) = 3q'+2` is odd while `S(q') = (3q'+1)/2^w < 3q'+2`; no merge.
- `o = 2`: with `Q = S(q')`, `S(p) = 2^w Q + 1` and `S(2^w Q+1) = 3·2^{w-2} Q + 1 > S(Q)` for `w >= 3`;
  the `4x+1` rule needs `w = 2`; no merge.
- `o = 3`: put `P_2 = 3·2^{w-2}Q + 1`, `Q_2 = S(Q) <= (3Q+1)/2`. For `w = 3`, `S(P_2) = 9Q+2 > Q_2`'s
  successor bound `(9Q+5)/4`; for `w >= 5`, `S(P_2) = 9·2^{w-4}Q + 1 > (9Q+5)/4`; for `w = 4`,
  `S(P_2) = (9Q+1)/2^s` and equality with `S(Q_2) = (9Q+3+2^{w'})/2^{w'+w''}` would force either
  `-2 = 2^{w'}` or `2^{w'} >= 9Q-1 > 3Q+1 >= 2^{w'}`; no merge.
Hence any equal-index merge has at least four odd steps after `p`: `d(L) >= L + 5`.

## 4. Two-adic locality

**Lemma C.** The first `N` parity bits of the orbit of `p = 3^{L+1}+2` depend only on `p mod 2^N`,
hence (as `3` has order `2^{N-2}` modulo `2^N`, `N >= 3`) only on `(L+1) mod 2^{N-2}`. The chain's
state after `n'` post-prefix steps is a function of those bits. Therefore: *whether the pair merges
within `n'` post-prefix Terras steps, and at which step and after how many odd steps, is a function of
`L mod 2^{n'-2}`.* The merge depth is a 2-adically locally constant function on the set where it is
finite.

## 5. Theorem M* (the exact tree for `L ≡ 1 (mod 4)`)

Enumerate the `2^{N-2}` classes `r = (L+1) mod 2^{N-2}` and run the chain from `(0,-1)` on the bits of
`3^r + 2 mod 2^N` (`check_merge_tree.py`, `N = 20`, `65,536` classes with `r ≡ 2 (mod 4)`, plus all other
classes as controls). Since absorption at post-prefix time `n' <= N` is decided by these bits, the
enumeration is a proof for every `L` in each class (with `L` large enough that the orbit has not yet
reached `1`, which holds for `L >= N`).

**Theorem M*.** Let `L ≡ 1 (mod 4)` and let `d(L)` be the merge index (equal odd index), `∞` if the
orbits meet only in the trivial cycle.
1. (PROVED, Lemma D) `d(L) >= L + 5` for every `L ≡ 1 (mod 4)`.
2. (PROVED by the enumeration, "if"; FINITE-EXACT, "only if") `d(L) = L + 5` **if `L ≡ 81 (mod 256)`**:
   the chain absorbs after exactly 10 post-prefix Terras steps, 4 of them odd. Conversely, among all
   merges occurring within 20 post-prefix Terras steps (all `2^{18}` classes), only this class gives
   offset 5; directly, every `L <= 1500` with `d(L) = L+5` is `≡ 81 (mod 256)`.
3. (PROVED per class) The first rows of the tree, as `d(L) - L` on residue classes of `L`:
   `6` on `145, 445 (mod 512)`, `305 (mod 1024)` and `1793 (mod 2048)`; `7` on `529, 829 (mod 1024)`,
   `1073, 1277, 1485 (mod 2048)` and `453, 733, 1969, 3329 (mod 4096)`; `8` on `785, 1085 (mod 2048)` and
   `2349, 3633, 3837, 4045 (mod 4096)`; `9` on `273, 573 (mod 4096)`. (These are all absorptions with
   post-prefix Terras time `<= 14`; `tree-normal.json` lists the classes through time 20.)
4. (FINITE-EXACT, with THM-4581 for the law) The number of classes modulo `2^{18}` absorbed at
   post-prefix times `10, 11, …, 20` is `1024, 1024, 768, 768, 640, 544, 496, 448, 420, 388, 362`; `89.5%`
   are unabsorbed at time 20. Consequently no bound `d(L) <= L + c` holds on the progression
   `L ≡ 1 (mod 4)`: at every depth `<= 20` a positive-density union of classes is unabsorbed, and
   THM-4581's survival law `c·n'^{-1/2}` (sketch level there) says the same at every depth.
5. (Controls) Every class with even `L` absorbs at `n' = 3` with one odd step (`d = L+2`), every class
   with `L ≡ 3 (mod 4)` at `n' = 4` with two (`d = L+3`): Theorem M recovered from the chain.
6. (Direct check) For all `L <= 1500` the class prediction agrees with the big-integer orbits (`0`
   mismatches on the 375 values `L ≡ 1 (mod 4)`); `286` of them merge before reaching `1`, the fraction
   rising from `0.57` on `L < 250` to `0.87` on `1250 <= L <= 1500`.

*Proof.* Lemma A (THM-4581's exact chain identity), Lemma B (explicit Terras computation), Lemma C
(order of `3` modulo `2^N`), Lemma D for the lower bound, and the finite enumeration, which is exact
because absorption by time `n'` reads only the first `n'` bits: the `L+5` class is the unique absorbing
class modulo `2^8` at `n' = 10`, and the later rows are read from the same table modulo `2^{n'-2}`.
Statement 4's "no uniform bound" is exact through depth 20 and inherits THM-4581's rate beyond. QED.

## 6. What this is and is not

It is the exact answer: the merge depth for `L ≡ 1 (mod 4)` is the absorption time of a specific
2-adic automaton (THM-4581's chain started at `(0,-1)` on the bits of `3^{L+1}+2`), locally constant
in `L`, with explicit first levels. It is not a formula, and it cannot be one: the chain is a
(deterministic sample of a) random walk with return drift, and THM-4581 shows the unabsorbed mass
decays only like `n'^{-1/2}`. The Mersenne twin (THM-4556 (ii)) left its odd-exponent case at exactly
this stage ("lag 1 in stopping time", converse FINITE-EXACT); the Fermat line gets the stronger
statement because the landing pair is literally consecutive integers, which is THM-4581's own
`z, z+1` setting. Nothing here bears on Fermat-prime completeness; for how the landing points
`L_m = (3^{2^{m-1}}+1)/2` do carry Pépin (`L_{m+1} = L_m^2 + (L_m-1)^2`, `F_n` prime iff `F_n | L_{2^n}`),
see kind-pasteur's addendum `fermat_landing_points_pepin_20261026.md`.
