---
id: HYP-9469
title: "After THM-4876: (1) the minimal deletion increment k_min(h) = min {(H(T) - H(T - v))/2 > 0 : H(T - v) = h} grows linearly, k_min(h) >= h/4 for h >= 5 (census: k_min(h)/h in [0.27, 0.40] for 5 <= h <= 39, orders of T up to 9); equivalently a vertex lying on an odd cycle of T raises the count of T - v by at least a fixed fraction; (2) for every k the set of h with an increment of exactly 2k is finite: increment 4 forces h in {1, 5}, 6 forces {3, 5, 9}, 8 forces {1, 5, 9, 11}, and in general h <= 3k - 1 (?); (3) the one-centre value set {1 + 2N(A)} has density exactly 1/2 among the odd numbers and its complement 7, 13, 21, 29, 35, 41, 43, 49, 55, 57, 59, 67, ... (not in OEIS) has a recognisable structure; (4) the end-count bound: a non-sink vertex ends at most (h - 1)/2 of the h Hamiltonian paths, with equality in the reversed-extreme-arc family 2^k + 1 (census: max e_a/h over non-sinks is 2/5, 4/9, 8/17, 16/33, 32/65 at orders 4-8); (5) what remains for Chan-Pak: completeness needs unions over realisers (the digital banks) or products; the two prime lanes p and 7p of THM-4094 stand"
status: >
  (1) NUMERICAL (orders of T up to 9, h <= 39, plus the order-10 run for h <= 60 when complete); the only proved
  cases are k_min(h) >= 2 for h >= 5 (THM-4876 (2)) and k_min(h) >= 3 for h notin {1, 5} (THM-4876 (3), census).
  (2) PROVED for 2k = 2 (THM-4876 (2)), NUMERICAL for 2k = 4, 6, 8. (3) NUMERICAL (200,000). (4) NUMERICAL.
  (5) a reading of THM-4876 (5) and the inheritance, no claim. Pre-registered tests below.
session: mac-mini-2026-11-04-chanpak
source: 05-knowledge/results/chanpak_increments_thm013_thm017_merge_20261104.md
scripts:
  - 04-computation/experiments/chanpak_spine_20261104/increment_ladder.py, increment2_order9.py, ext_windows.py, spine_numeration.py
related:
  - 01-canon/theorems/THM-4876-*; THM-4873 (end-count lemma: no vertex ends exactly H - 1 paths, exactly H - 2 iff TT_j => C3), THM-1370, THM-4094, HYP-9453, HYP-9420, HYP-9029
---

# HYP-9469 — the local structure of the Hamiltonian-path spectrum

## Statements

1. **Linear minimal increment.** Let `k_min(h)` be the least positive value of `(H(T) - H(T - v))/2` over all tournaments
   `T` and vertices `v` with `H(T - v) = h`. Then `k_min(h) >= h/4` for `h >= 5`. Census (THM-4876 (5)): `k_min = 2, 3, 4, 5,
   5, 5, 6, 7, 9, 9, 8, 10, 9, 11, 12, 11` at `h = 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39`. By Claim
   A the increment is `sum_C H(T - V(C))` over the odd cycles through `v`; the conjecture says that a vertex on an odd
   cycle cannot have all those complements small while `T - v` has many paths. The sole-exit cut gives the increment
   `(h - e_a)/2`, so (1) would follow from (4) for those cuts, but the minimum is not always a sole-exit cut.
2. **Finiteness of each increment level.** For each `k >= 1` the set `{H(T - v) : H(T) - H(T - v) = 2k}` is finite,
   with largest element at most `4k` (the restatement of (1)). Proved for `k = 1` (`{1, 3}`); census (orders of `T`
   up to `9`, `endcount_ratio.py`): `k = 2`: `{1, 5}`; `k = 3`: `{3, 5, 9}`; `k = 4`: `{1, 5, 9, 11}`; `k = 5`:
   `{1, 3, 5, 9, 13, 15, 17}`; `k = 6`: `{3, 5, 11, 13, 15, 17, 19}` (stable from order `6` on). A proof would classify
   the configurations of `k` odd cycles through `v` whose complement counts sum to `k`, as THM-4876 (2) does for one
   triangle.
3. **One-centre density.** The set `{1 + 2 N(A)}` of THM-4876 (6) contains exactly half of the odd numbers in the
   limit (`50,006` of `100,000` below `200,000`), and the missing set `7, 13, 21, 29, 35, 41, 43, 49, 55, 57, 59, 67,
   71, 75, 77, 83, 89, 93, 101, 105, 109, 111, 113, 119, 125, ...` has a description in terms of the binary word `A`.
   Neither the attained nor the missing set is in OEIS (2026-11-04).
4. **End-count bound.** If `a` is not a sink of `X` and `H(X) = h`, then `e_a(X) <= (h - 1)/2`, where `e_a` counts the
   Hamiltonian paths ending at `a`; equality holds in the reversed-extreme-arc family (`h = 2^k + 1`, `e_a = 2^(k-1)`
   at the vertex `k + 1` of `TT_(k+1)` with the arc `(1, k+1)` reversed). Census (`endcount_ratio.py`): the maximum of
   `e_a/h` over non-sinks is `2/5, 4/9, 8/17, 16/33, 32/65` at orders `4, 5, 6, 7, 8`, always `(h - 1)/(2h)`.
   THM-4873 (i) gives the first two steps (`e_a != h - 1`; `e_a = h - 2` only at `h = 3`). For the sole-exit cut the
   increment is `(h - e_a)/2`, so this bound gives `k_min(h) >= (h + 1)/4` on sole-exit cuts.
5. **What remains for Chan-Pak.** THM-4876 (5) rules out every induction that extends one realiser by one vertex. The
   open statement is unchanged: strong realisation of the two lanes `p` (odd prime `p != 7`) and `7p` (THM-4094), or
   the tiling law HYP-9029, or the density of `B`-sets HYP-9453 (2); the digital bank proves it to `10^8`.

## Pre-registered decisive tests

- (1): the order-10 census (`increment2_order9.py`, all tournaments `X` of order `9` with all cuts) extends `k_min(h)` to
  `h <= 60`; the prediction is `k_min(h) >= h/4` throughout, with `k_min(h)/h` staying in `[0.25, 0.45]`.
- (2): the increment-`14` and `-16` levels at orders `<= 9`: predicted finite, with largest `h <= 28` and `<= 32`
  (the computed levels `k <= 6` have `h_max = 3, 5, 9, 11, 17, 19`).
- (4): the end-count bound over all tournaments of order `9` (`max e_a/h` over non-sinks, predicted `(h - 1)/(2h)` with
  the maximum at the Fermat-type count `129`).
