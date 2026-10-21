# Session mac-mini-2026-10-20-debtcones: debt cones covering 63, the delta = k zero-heavy layer, the order-ten and (n-5,4) endpoint layers, the 3-4-3 / 4-2-4 reflections, the annulus increment at finite horizons, and the idoneal dictionary of 105

Status: pre-audit. Canon written: THM-4733, THM-4734, THM-4735, THM-4736, HYP-9348, HYP-9349.
Experiments: `04-computation/experiments/{endpoint_layers_20261020, debt_cones_20261020,
zero_heavy_20261020, annulus_horizon_20261020}/`. All inherited packet engines were imported
unchanged (endpoint_three_repairs_20261020h/probe*.py, native_cancellation_105_20261018g.py).

## 1. The owner's directive and the repo state it referred to

The prompt named objects of the 2026-10-20 codex-lane packets (charge_palettes_boundary_g,
zero_heavy_ports_g, endpoint_two_repairs_g, endpoint_three_repairs_h, kernels_games_detours_f,
trijection_reset) and of the 2026-10-18 Collatz ROOT lane (guarded_children_105_g,
native_cancellation_105_g) plus THM-4566 (idoneal planes). Dictionary recovered: debt cone =
`0 <= D <= c(lambda - 1)` in the affine first-exit ledger; "covering 63" = source 63's prefix
fails the `c = 1` cone; zero-heavy boundary = the mass layers of pair-free weighted cores at
feedback budget `k` (exact through deficit `k-1`, `42 != 45` at deficit `k`); `(n-4,3)` detours
= three-black K moves for closed blue lines with endpoint scores `(n-4,3)` (proved `n >= 11`,
order 10 unclassified); `(n-5,4)` the next layer; 3-4-3 / 4-2-4 = rank and sign levels of the
ten proper divisors of `p^2qr`; movable cut = the packet's section-9 word; black-core
connectivity = all non-pure-blue classes in one black component; 153N = exact new-annulus source
norm `S(64N) - S(N)`; guards mod 15 / mod 35 = the compatible registers encoding `s mod 105`.

## 2. Results

**THM-4735 (debt cones).** The payment theorem holds for cone constants up to
`c_max(a) = (2^a + 4)/(2^a - 3)` at final valuation `a >= 3` (`< 8` at `a = 2`). Source 63
(34 odd steps, 5 resets, final valuation 4, prefix ratio `1.0573`) is covered by the
valuation-4 cone (`20/13`). Counts (exact rationals): `c = 1` cone 47,424 / 50,000 (the
packet's value) and 473,943 / 500,000; valuation cones 48,741 and 487,405; optimal terminal
splits 49,997 and 499,958. Obstruction: every split value is `>= D_w`, so `D_w >= 2` blocks
every threshold-two certificate: 3 sources below `10^5` (4575, 42703, 65567), 42 below `10^6`
(max debt 3.125 at 437,247). All 500,000 first exits pay directly (`D_w < n - h`). Universal
coverage by threshold-two cones is false; the first-exit payment conjecture is HYP-9349 (a).

**THM-4736 (zero-heavy deficit k).** `Z_{k,k} = (2k-1)!! (C(3k,k) - 1)`,
`Z^s_{k,k} = I_k C(3k-2,k) - (2k-3) I_{k-1}`: the deficit identity forces a perfect matching
on `2k` unit endpoints and two families (all gaps present with deficits summing to `k`; one
gap missing), and family (B) loses exactly the records whose deleted gap separated a matched
pair (homogeneous pair): `(2k-1)!!` collapses, `(2k-3) I_{k-1}` in the strong case. Census:
42 (k = 2, the packet's compiler value), 1245 (k = 3); strong 11, 344, 15490 (k = 2, 3, 4),
with exactly the predicted 3, 15 / 1, 6, 50 collapses and zero wrong-feedback records.

**THM-4734 (endpoint layers).** Order ten, `(6,3)`: all 86,016 presentations have three-black
detours (86,013 by the inherited bank, 3 by search with independent replay: closed middles,
`I = {3,6} = J`, `{5,6}`, `{7,8}`) -- the `(n-4,3)` layer is now settled for every `n >= 6`.
Layer `(n-5,4)`: all 10,240 (n = 9) and 143,360 (n = 10) presentations have detours (325 and
2,310 by search); the misses of the inherited bank sit on the consecutive-port patterns
`{2,3,4}`, `{2,n-3,n-2}`, `{n-4,n-3,n-2}`, and no single template covers such a pattern
(best: exactly half).

**THM-4733 (3-4-3 / 4-2-4).** On tilings the Klein group `<kappa, sigma>` is {complement
line, converse, path reversal `R_P = sigma kappa`}. The divisor Klein set `(C, Q, CQ)` of
`p^2qr` is the fold of the four mirror-moved tiles of the order-five cube along one mirror
orbit (`a` = backward count in the folded orbit; rank = Hamming weight gives 3-4-3; the sign
of folded minus unfolded weight gives 4-2-4) and also embeds as an invariant ten-pattern
subset; it cannot embed in the full cube of any order in either assignment, because `kappa`
is fixed-point free and path reversal has no fixed presentation (the mirror-fixed tile) while
`Q` fixes four divisors and `CQ` fixes `pq, pr`.

**Annulus increments (HYP-9349 (b)).** With the checker's own rows: `S(64N) - S(N) = 153N`
reproduced at `N = 4096, 5997, 11002, 16384`; the finite-horizon increment
`E(64N,D) - E(N,D)` is at most `0.68, 0.65, 0.62, 0.59` times `153N` over ALL horizons `D`,
and the cross term `sum_u cross(W_u, H_u)` is negative at every horizon on every scale.

**Idoneal dictionary (HYP-9349 (c)).** `105` is idoneal; discriminant `-420` has eight
one-class genera separated by `(n/3), (n/5), (n/7)` and the mod-4 character; each form
represents exactly six `(mod 15, mod 35)` register pairs. The lane's guards are the genus
characters of `-420`; no energy transport is claimed.

## 3. Open (HYP-9348, HYP-9349)

Structural `(n-5,4)` theorem (two ledgers per consecutive-port pattern; orders 11-14);
black-core connectivity via a `k`-independent movable cut in trijection coordinates; the
deficit-`(k+1)` zero-heavy layer; `p^k q r` versus the order-six cube; first-exit payment for
all sources; the annulus increment constant `c < 1` for all `N` and a sign proof for the
cross term; whether the ROOT terminal sees the class of discriminant `-420`.

## 4. Files

endpoint_layers_20261020: order_ten_bank.py/.out/.json, order_ten_residuals.py/.out/.json,
endpoint_four_bank.py/.out, endpoint_four_n9.json, endpoint_four_n10.json,
uniform_ledger_test.py/.out, divisor_plane_check.py/.out. debt_cones_20261020:
cone_splits.py/.out/_1e6.out/.json, cone_obstruction_check.py/.out/_1e6.out, genus_420.py.
zero_heavy_20261020: zero_heavy_delta_k.py/.out, zero_heavy_delta_k_strong.py/.out/_k4.out.
annulus_horizon_20261020: annulus_horizon_increment.py/.out/_b.out, annulus_horizon_N*.json.
