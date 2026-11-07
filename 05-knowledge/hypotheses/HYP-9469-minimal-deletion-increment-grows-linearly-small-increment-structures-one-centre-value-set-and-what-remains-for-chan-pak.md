---
id: HYP-9469
title: "After THM-4876: (1) the minimal deletion increment k_min(h) = min {(H(T) - H(T - v))/2 > 0 : H(T - v) = h} grows linearly, k_min(h) >= h/4 for h >= 5 (census: k_min(h)/h in [0.26, 0.40] for 5 <= h <= 59, orders of T up to 10; audit C: holds for every h attained at order 10, least ratio 33/129); equivalently a vertex lying on an odd cycle of T raises the count of T - v by at least a fixed fraction; (2) for every k the set of h with an increment of exactly 2k is finite: increment 4 forces h in {1, 5}, 6 forces {3, 5, 9}, 8 forces {1, 5, 9, 11}, and in general h <= 4k (levels k = 7..12 have h_max = 23, 29, 33, 33, 39, 37 at orders <= 10, audit C); (3) the one-centre value set {1 + 2N(A)} has density exactly 1/2 among the odd numbers and its complement 7, 13, 21, 29, 35, 41, 43, 49, 55, 57, 59, 67, ... (not in OEIS) has a recognisable structure; (4) the end-count bound e_a <= (h - 1)/2 for non-sinks, now PROVED as THM-4876 (7) (equality on 2^k + 1; census maxima 2/5, 4/9, 8/17, 16/33, 32/65, 64/129 at orders 4-9), which proves the sole-exit half of (1); (5) what remains for Chan-Pak: completeness needs unions over realisers (the digital banks) or products; the two prime lanes p and 7p of THM-4094 stand; (6) the ordered-path lemma: for every cut (A, B) with an arc from A to B, the number q of Hamiltonian paths of X listing all of B before all of A is at most (h - 1)/2 (every cut of every tournament of orders <= 8, 1.7 million cuts at order 8, maximum at the sole-exit cuts of R_n) - since H(X + v) >= 2h - q is PROVED (slot counting), the lemma gives (1) sharply and (2) with h <= 4k - 1, and THM-4876 (7) is its case |A| = 1"
status: >
  CURRENT as of Session214 (2026-11-07): (1) and (6) PROVED in all orders by the independently audited forced-shuffle
  injection. For a positive one-vertex increment 2k from parent count h, k >= ceil((h+1)/4), hence h <= 4k-1.
  (2) PROVED: every fixed positive increment has finitely many possible parent counts; the exact all-order parent
  sets for increments 2,4,6,8 are {1,3}, {1,5}, {3,5,9}, {1,5,9,11}. Larger explicit census lists below retain their
  original finite scope. (3) remains NUMERICAL. (4) remains PROVED by THM-4876(7). (5) is unchanged: unrestricted
  Chan-Pak completeness remains OPEN. See the Current resolution section and linked proof/audit.
historical_status_before_session214: >
  (1) NUMERICAL (orders of T up to 10: h <= 59 by the session's order-9 run, every h attained at order 10 by audit C's
  census, least ratio 33/129); the only proved cases are k_min(h) >= 2 for h >= 5 (THM-4876 (2)) and the sole-exit
  cuts, k >= (h + 1)/4 (THM-4876 (7)); k_min(h) >= 3 for h notin {1, 5} is FINITE-EXACT (orders of T up to 10), not
  proved. Audit C (2026-11-04, chanpak_increments_audit_C_20261104.md) PASS WITH CORRECTIONS, applied (C8-C13,
  MISTAKE-715).
  (2) PROVED for 2k = 2 (THM-4876 (2)), NUMERICAL for 2k = 4, 6, 8. (3) NUMERICAL (200,000). (4) PROVED (THM-4876 (7),
  written after this hypothesis was filed; audits C and D checked it on every class of orders <= 9, resp. 10, maximum
  64/129 and 128/257; audit D (chanpak_endcount_audit_D_20261104.md) PASS WITH CORRECTIONS, minor, applied). (6) NUMERICAL (orders of X up to 8, every cut; the inequality H(X + v) >= 2h - q is PROVED).
  (5) a reading of THM-4876 (5) and the inheritance, no claim. Pre-registered tests below.
session: mac-mini-2026-11-04-chanpak
source: 05-knowledge/results/chanpak_increments_thm013_thm017_merge_20261104.md
scripts:
  - 04-computation/experiments/chanpak_spine_20261104/increment_ladder.py, increment2_order9.py, ext_windows.py, spine_numeration.py
related:
  - 01-canon/theorems/THM-4876-*; THM-4873 (end-count lemma: no vertex ends exactly H - 1 paths, exactly H - 2 iff TT_j => C3), THM-1370, THM-4094, HYP-9453, HYP-9420, HYP-9029
---

# HYP-9469 — the local structure of the Hamiltonian-path spectrum

## Current resolution: statements (1), (2), and (6)

**PROVED, independently audited, 2026-11-07.**
[Session214 forced-shuffle proof](../results/session214_chan_pak_combined_proof_20261107/forced-shuffle/PROOF.md)
and [independent audit](../results/session214_chan_pak_combined_proof_20261107/forced-shuffle-audit/FORCED_SHUFFLE_AUDIT.md).

1. **Linear positive increment, all orders.** Every positive one-vertex extension from `h` paths has increment
   `2k` with `k >= ceil((h+1)/4)`. Equivalently, its count is at least `(3h+1)/2`. The multiplicative constant
   `3/2` is asymptotically optimal. Every induced extension, of any size, either preserves `h` or satisfies
   that same lower bound.
2. **Finite increment levels, all orders.** A positive increment `2k` forces `h <= 4k-1` and child count
   at most `6k-1`. SCC factorization and the cited strong-order floor reduce increments `2,4,6,8` to strong
   marked cores of order at most six. Their exact parent sets are respectively `{1,3}`, `{1,5}`,
   `{3,5,9}`, and `{1,5,9,11}`. The first level was already proved; the other three are now exhausted
   in all orders. The larger numerical lists below are not promoted to complete classifications.
6. **Ordered-path half-bound, all orders.** For a nonempty cut `(A,B)` containing an `A->B` arc, put
   `p=H(X[A])H(X[B])` and let `q` count Hamilton words listing all `B` before all `A`. A forced shuffle
   injects all `p` monochromatic path pairs into nonordered words, proving `h >= p+q` and
   `q <= (h-1)/2`. Combining this with the inherited slot count proves (1). Without an `A->B` arc,
   the extension preserves `h`; this hypothesis is essential.

Statements (3), (4), and (5) keep their original scope. In particular, this does not settle the one-centre
value-set density conjecture or unrestricted Chan–Pak completeness. Generic directed-word merging and
restriction recovery are inherited from Busch via THM-4745; the new step forces a nonordered image.
No external literature-priority claim is made.

## Historical statements and original census

The following records preserve the original November 4 filing, audits, computations, and proposed tests.
Their numerical-only descriptions and open proof language for (1), (2), and (6) are superseded by the
current resolution above. Other statements retain their stated status.


1. **Linear minimal increment.** Let `k_min(h)` be the least positive value of `(H(T) - H(T - v))/2` over all tournaments
   `T` and vertices `v` with `H(T - v) = h`. Then `k_min(h) >= h/4` for `h >= 5`. Census (THM-4876 (5)): `k_min = 2, 3, 4, 5,
   5, 5, 6, 7, 9, 9, 8, 10, 9, 11, 12, 11, 14, 15, 14, 13, 13, 15, 14, 18, 18, 17` at `h = 5, 9, 11, ..., 39, 41, ..., 59`
   (orders of `T` up to `10`; least ratio `13/49 = 0.265`). The sole-exit cuts satisfy it by THM-4876 (7). By Claim
   A the increment is `sum_C H(T - V(C))` over the odd cycles through `v`; the conjecture says that a vertex on an odd
   cycle cannot have all those complements small while `T - v` has many paths. The sole-exit cut gives the increment
   `(h + sigma - e_a)/2` with `sigma >= 0` (THM-4862 (a); `sigma = 0` in only 53 of the 420 `(T, a)` of orders `2-6`,
   audit C), so (1) follows from (4) for those cuts (THM-4876 (7)), but the minimum is not always a sole-exit cut.
2. **Finiteness of each increment level.** For each `k >= 1` the set `{H(T - v) : H(T) - H(T - v) = 2k}` is finite,
   with largest element at most `4k` (the restatement of (1)). Proved for `k = 1` (`{1, 3}`); census (orders of `T`
   up to `9`, `endcount_ratio.py`): `k = 2`: `{1, 5}`; `k = 3`: `{3, 5, 9}`; `k = 4`: `{1, 5, 9, 11}`; `k = 5`:
   `{1, 3, 5, 9, 13, 15, 17}`; `k = 6`: `{3, 5, 11, 13, 15, 17, 19}` (complete from order `7` of `T`, i.e. order `6` of `X`, on; unchanged at order `10`, audit C; `endcount_ratio_n8.out`). A proof would classify
   the configurations of `k` odd cycles through `v` whose complement counts sum to `k`, as THM-4876 (2) does for one
   triangle.
3. **One-centre density.** The set `{1 + 2 N(A)}` of THM-4876 (6) contains exactly half of the odd numbers in the
   limit (`50,006` of `100,000` below `200,000`), and the missing set `7, 13, 21, 29, 35, 41, 43, 49, 55, 57, 59, 67,
   71, 75, 77, 83, 89, 93, 101, 105, 109, 111, 113, 119, 125, ...` has a description in terms of the binary word `A`.
   Neither the attained nor the missing set is in OEIS (2026-11-04).
4. **End-count bound (now THM-4876 (7), PROVED).** If `a` is not a sink of `X` and `H(X) = h`, then `e_a(X) <= (h - 1)/2`, where `e_a` counts the
   Hamiltonian paths ending at `a`; equality holds in the reversed-extreme-arc family (`h = 2^k + 1`, `e_a = 2^(k-1)`
   at the penultimate vertex `k + 1` of `TT_(k+2)` with the arc `(1, k + 2)` reversed). Census (`endcount_ratio.py`): the maximum of
   `e_a/h` over non-sinks is `2/5, 4/9, 8/17, 16/33, 32/65, 64/129, 128/257` at orders `4, ..., 10` (`9` by the session,
   `10` by audit C), always `(h - 1)/(2h)` at `h = 2^(n-2) + 1`.
   THM-4873 (i) gives the first two steps (`e_a != h - 1`; `e_a = h - 2` only at `h = 3`). For the sole-exit cut the
   increment is `(h + sigma - e_a)/2 >= (h - e_a)/2` (THM-4862 (a)), so this bound gives `k >= (h + 1)/4` on sole-exit
   cuts (THM-4876 (7)).
5. **What remains for Chan-Pak.** THM-4876 (2) rules out the induction through `[h, 2h - 1]` from one realiser at every order, and THM-4876 (5)
   leaves no window induction through order `8` (FINITE-EXACT; the realisers of `h = 5` are the exception). The
   open statement is unchanged: strong realisation of the two lanes `p` (odd prime `p != 7`) and `7p` (THM-4094), or
   the tiling law HYP-9029, or the density of `B`-sets HYP-9453 (2); the digital bank proves it to `10^8`.
6. **Ordered-path lemma (NUMERICAL, orders of `X` up to `8`; the reduction is PROVED).** For a cut `(A, B)` of `X`
   with at least one arc from `A` to `B` (equivalently `v` lies on an odd cycle of `X + v`), let `q` be the number of
   Hamiltonian paths of `X` that list all of `B` before all of `A`. PROVED: `H(X + v) >= 2h - q` - every Hamiltonian
   path of `X` has a slot for `v` (between a consecutive `B`, `A` pair, in front if it starts in `A`, at the back if it
   ends in `B`), at least two slots unless it is ordered, and the slot map is injective because deleting `v` returns the
   path. Conjectured: `q <= (h - 1)/2`. Checked on every cut of every tournament of orders `4-8`
   (`ordered_paths_lemma.py`, subset tables; at order `8`: 1,746,432 cuts with an arc from `A` to `B`, no violation,
   maximum `q/h = 32/65` at the sole-exit cut of `R_8`, maximum `26/53` over cuts with `|A|, |B| >= 2`, and the bound
   `2h - q` is attained at `1921` cuts). The lemma would give (1) in the sharp form `k >= (h + 1)/4` for every cut
   and (2) with `h <= 4k - 1`; THM-4876 (7) is its case `|A| = 1` (`q = e_a`). Partial proof: moving the last vertex
   `x_r` of an ordered path into the `B`-block before its first out-neighbour there is an injection on the ordered
   paths whose last vertex has an out-neighbour in `B` (checked to order `7`), and the dual map handles those whose
   first vertex has an in-neighbour in `A`; ordered paths with neither property exist (up to `4/9` of `h` at order
   `7`) and the two images can collide, so no injection is in hand.

## Pre-registered decisive tests

- (1): the order-10 census (`increment2_order9.py`, all tournaments `X` of order `9` with all cuts) extends `k_min(h)` to
  `h <= 60`; the prediction is `k_min(h) >= h/4` throughout, with `k_min(h)/h` staying in `[0.25, 0.45]`.
- (2): the increment-`14` and `-16` levels at orders `<= 9`: predicted finite, with largest `h <= 28` and `<= 32`
  (the computed levels `k <= 6` have `h_max = 3, 5, 9, 11, 17, 19`).
- (4): the end-count bound over all tournaments of order `9` (`max e_a/h` over non-sinks, predicted `(h - 1)/(2h)` with
  the maximum at the Fermat-type count `129`).
