# Prior sessions on the +2 and ×2 recursions (survey, with the rows/digits reading)

Session klein-2026-11-16-splice; compiled from the repository by a read-only survey (58 items; the four groups below
keep the files' own status words). The unifying reading is COLLATZ_ROWS.md §3: both recursions are directions of the
one array `h + 2^k B` of the root response (THM-4865).

## +2 (additive) moves — the `B`-direction and the dials

- THM-4876: `H(T) − H(T − v) = 2 Σ_(odd cycles C ∋ v) H(T − C)`; a deletion lowers the count by exactly 2 only from
  `H ∈ {3, 5}` (PROVED). HYP-9469: increment floors `k ≥ ⌈(h + 1)/4⌉` (partly PROVED). THM-5124: two raising vertices
  with zero interaction only as `c, 3c, 3c, 5c` (PROVED).
- Ports and banks: THM-4872 / HYP-9465 (the unit-start ±2 port pair of order 7, wrappers `h + 2^m F` differing by 2 at
  every clock; FINITE-EXACT), THM-5272 (bank rigidity, no complete 3-bank through order 9), THM-5282 (C) (`K + 90, K + 92,
  K + 94` for all `A, B`; PROVED), THM-5302 (line banks never longer than `h`; PROVED), HYP-9822 / HYP-9842 (a)
  (no three arcs give eight consecutive odd counts, exhaustive to order 8), HYP-9852 (additivity limit: local dials stop
  at two bits; OPEN).
- Exact +2 seams and shields: chanpak_signed_seam_bank_20261116 (`+2` with three flips, the complete ordered +2 graph of
  the mirrored transitive core; proved all orders), wheels_spine14_shield_20261116 / THM-5255 (F) (factor-1 shields keep
  the +2/+4 seed swaps; 2.2 % of addresses have an exact +2 partner), physical_radix_carry_20261116, radix_three_shields_20261116
  (a width-four window with all its +2 moves), chanpak_mixed_defects_20261115 (721, 723, 725, 727 on one tournament),
  chanpak_bott_cayley_periods_20261115 and chanpak_qrh_residual_20261115 (three-count windows with both +2 steps within
  `O(N^(3/8))`, then `O(N^(1/3))`), chanpak_exact_seams_douglas_20261115 (no strong tournament of order ≥ 8 has
  `H = f(n) + 2`), chanpak_waves_20261113 wave 37 (fibre swaps in steps of 24).
- Two-cycle gluings: THM-4877 (the two-cycle converse) and HYP-9472 (every odd `N ≥ 23` but 67 with at most two
  Hamiltonian cycles).

## ×2 and 2d − 1 (doubling) lifts — the `k`-direction, the rows

- The signed descent: six_tutte_collatz_20261116/COUNTS.md and THM-5262 (inverse branches `2d − 1`, `(4d − 3)/3`,
  `(4d + 1)/3`; `(h + 1)/2` is a parent for 38 of 54 counts; PROVED/FINITE-EXACT); chanpak_signed_descent_lifts_20261115
  (the halving branch is not always reversible by adding vertices; HYP-9812 `|H(T) − 2H(T − x)| ≤ C`, OPEN);
  inverse_fences_landing_20261116 (one-apex words realise each branch on unbounded families, e.g. `63·2^k + 17 →
  63·2^(k+1) + 33`; indefinite `2h − 1` on one rooted interface, GUARDED_LIFTS.md).
- THM-5272 (C): the dyadic step `D: h ↦ 2h − 1` and the triadic step `T` satisfy `D(T(h)) − T(D(h)) = 2`, rays meeting
  only at 9 (PROVED) — the shear of COLLATZ_ROWS.md §3 in one line.
- Rows: THM-5126 (every odd `N > 1` lies on one line `L_c = {1 + c 2^s}`; the 65 lines with `c ≤ 129` are in the
  spectrum by explicit rays; PROVED), THM-5145 (`c 2^e + r` and `1 + 2^a + 2^b`), HYP-9658 (every odd `c` has a
  unit-offset seed; OPEN), HYP-9662 (seed and anti-seed `1 ± 2M`), THM-5130 (sign-flip towers `±1 + 2^t M_k`),
  THM-5109 (`B = 2^j` makes the near-join a pure cone `h + 2^K`, never reaching a binary expansion beginning 11),
  THM-5129 / HYP-9661 (one-vertex sources `D + 2^t U`, slopes `Slopes_0(D) 2^j`, at least `N/64` sources on a line).
- Near-join: spectrum_analogies_nearjoin_20261102 and HYP-9453 (`h + 2^k B(T, t_0)`; every odd number except
  1, 7, 11, 13, 21, 47, 67 to 22,227 from bases ≤ 8; OPEN). session222_tiling_recursions_20261107 (four genuine
  doubling mechanisms; the odd-chain grid `n = 2^k (2j + 1)` realised as `F(j, k)`), session265_converse_dyadic_20261111,
  coalescent_27_63_20261116 (no 27-path tournament embeds in a 63-path witness).
- HYP-9842 (c): a count-doubling generator of Price's Fibonacci-box tree (OPEN QUESTION).

## 3H + 2σ and the ternary children — the shear between rows

- THM-4862 (one-exit law `2h + σ − e_a`; every odd number in `[307, 26551]` is `3H + 2σ`), HYP-9454 (sweep; OPEN),
  THM-4865 (tribonacci tower, `3h + 2σ = h + 2B_a`), THM-4873 / HYP-9466 (digit heredity: no C2 grandchild has a digit-1
  mark), THM-4875 / HYP-9468 (exact sink digit; one-apex joins give every odd number to `2^32` except 7 and 21),
  clone_clocks_negative_collatz_20261113 (`3H + 2 = T(2H + 1) = 2T(H) + 1`; digits `3r − 1, 3r, 3r + 1` = the `3x − 1`
  step, pure tripling, the `3x + 1` step), prime_ternary_lift_20261113 (`H = 15·2^m − 1` with strong children `3H + 2`,
  `3H + 4`), boundary_carry_returns_20261103 (`H_m = 8·3^m + 6·2^m − 5` realises all four two-step `3H + 2/+4` branches),
  golden_inverse_defects_20261103 / target_directed_cores_20261112 (guarded singletons, a selective bridge to ready `3H + 2`),
  THM-4881 / HYP-9772 (exact tripling `H(u B^k) = 3^(k−3) H(u B^3)`; every odd `N ≥ 36835` is `3H(A) + P_2(A)`).
- The 2026-11-16 digit programme: fresh_unit_prime_gluing_20261116 (triple splice `3P + (d, 0, 0, 0)`, `d ∈ {0, 2}`;
  quadruple splice `d ∈ {0, 4}`; ternary deposit `d ∈ {0, 2, 4}` on one family; the obstruction `H = 2E − g` of the
  source family; CARRY.md's literal second-step boundary; DECODER.md: no three consecutive odd counts share `(E, O, 3)`
  with `E ≠ O`). CONTINUATION.md of this packet continues exactly there.

## Chains / rows of odd numbers

- THM-5126 (B) (the lines `L_c`), THM-5272 (D) / THM-5262 (D) (dyadic rays `55, 109, 217, …` and `47, 93, …, 2945`; the
  backward descent tree with Haar-average level size `(5/3)^k`), clone_clocks (dyadic rays shadowing Collatz cycles,
  `462·8^t − 5 → 462·9^t − 5`), session222 §3 (the odd-chain grid). COLLATZ_ROWS.md §2 adds: the Collatz map in count
  coordinates is the Syracuse map on row indices, and the seven exceptions are the first two rungs of `L_3`, `L_5` plus
  the heads of `L_23`, `L_33`.

## Verdict

Every "+2" mechanism in the repository is a change of the root response `B` (or of a port/dial) at fixed height, and
every "×2" mechanism is a change of height `k` at fixed `(h, B)`; THM-5272 (C)'s commutator `D T − T D = 2` is the
statement that the two directions do not commute by exactly one unit step. The missing piece in both directions is the
same (HYP-9453's density of `B`-sets, HYP-9658's seeds): which pairs `(h, B)` occur — the continuation problem.
