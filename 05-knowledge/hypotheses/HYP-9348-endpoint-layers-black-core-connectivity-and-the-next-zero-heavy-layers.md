---
id: HYP-9348
title: "After THM-4733/4734/4736: the structural (n-5,4) endpoint theorem and its two-ledger patterns, black-core connectivity through the trijection coordinates, the deficit-(k+1) zero-heavy layer, and the fixed-tile-free comparison of arithmetic Klein sets"
status: OPEN (pre-audit); every item carries its finite evidence and the first untested case
source: mac-mini-2026-10-20-debtcones
depends_on:
  - 01-canon/theorems/THM-4734-endpoint-three-layer-closed-at-order-ten-and-the-endpoint-four-layer-at-orders-nine-and-ten.md
  - 01-canon/theorems/THM-4736-the-delta-equals-k-zero-heavy-layer-exact-counts-missing-gaps-and-the-adjacent-matched-pair-collapse.md
  - 01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md
related:
  - 05-knowledge/results/endpoint_three_repairs_20261020h/README.md
  - 05-knowledge/results/charge_palettes_boundary_20261020g/zero_heavy.md
  - 05-knowledge/hypotheses/HYP-9345-merged-metagraph-open-laws-after-thm-4729-4730.md
---

# HYP-9348 -- the next layers

## (a) The structural (n-5,4) theorem

**Claim.** Every closed blue marked line with endpoint scores `(n-5, 4)` has a three-black
detour to its literal target for every `n >= 7`. Evidence: exact at `n = 9, 10` (THM-4734);
the packet's general-`k` lemmas (cyclic cut and common final word black for `n > 14`,
movable cut) cover the non-adjacent closed rows from `n = 15`. What is missing: the rows
where the three ports are consecutive at the start (`I = {2,3,4}`), at the end
(`{n-4,n-3,n-2}`, closed and open) or split (`{2,n-3,n-2}`); in each of these the best single
word template covers exactly half of the pattern (uniform_ledger_test.out), so the proof
needs two complementary ledgers per pattern, separated by one middle arc, as in the `k = 3`
`epsilon` branches; and the orders `11..14`, where the cyclic-cut bound fails. First test:
find the splitting arc for closed `{2,3,4}` at `n = 10` from the 64 distinct certificate
templates of the search.

## (b) Black-core connectivity

**Claim** (the lane's conjecture; THM-4731 verifies it for `n <= 8`): all non-pure-blue
classes lie in one black component of the line metagraph. The endpoint layers remove one
family of potential bridges at a time (`k = 2` all `n`, `k = 3` all `n`, `k = 4` at
`n <= 10`). In the trijection coordinates (graphs on ordered `[n-1]`; THM-4733 Proposition A)
a K move along `Q` is "change the marked path, then complement"; a three-black detour is a
closed walk of three coordinate changes whose complement bits cancel to `E(P)`. Question:
is there a `k`-independent movable-cut word, written in these coordinates, that repairs every
closed line whose port set is not an interval? The packet's section-9 words are the case
`k = 3`; the half-covering template of THM-4734 section 2 is their `k = 4` shadow.

## (c) The order-ten remainder is empty; what about order 17

The packet's second fixed-word obstruction (`n = 17`, reflection-product orbits `3,3,5,3,3`)
is bypassed by the movable cut at `k = 3`. For `k = 4` the analogous orbit obstructions of
the common final words should be listed (the orders where `r_Q r_P` has an orbit that a port
set of size three can balance); THM-4734's data suggest `n = 11..14` only.

## (d) The deficit-(k+1) zero-heavy layer

THM-4736 gives `Z_{k,k}` and `Z^s_{k,k}`. At deficit `k + 1` the identity admits a shared
endpoint (`t = 2k - 1`, defect graph a path `P_3` plus a matching), two missing gaps, a missing
gap plus reduced gaps, and three-endpoint unit intervals. **Claim.** `Z_{k,k+1}` is still a
finite sum of matching-type counts with explicit collapse corrections; the first unknown value
is `Z_{2,3}` (weight 11 at `k = 2`), computable by the packet's complete `k = 2` compiler, and
`Z_{3,4}` (weight 23) by an extension of zero_heavy_delta_k.py with shared endpoints. The
strong sequence `11, 344, 15490, 1520100` (`k = 2..5`, the last predicted) is not in OEIS
(not checked yet).

## (e) Arithmetic Klein sets and the fixed-tile-free cube

THM-4733 shows `D(p^2 q r)` is a fold and a sub-object of the four mirror-moved tiles of the
order-five cube and cannot embed in the full cube. **Question.** For `N = p^k q r` (the
packet's `(2a - k(b+c), b-c)` chart, `4k + 2` points) is the corresponding object the fold of
the `2k` mirror-moved tiles of a suitable order along `k` orbits, and which orders `n` have
`floor((n-1)/2) = f` fixed tiles compatible with the sign split `2k, 2, 2k` (even `k`) versus
`2k+1, 0, 2k+1` (odd `k`)? The odd-`k` case has no `CQ`-fixed divisor and therefore has no
fixed-tile obstruction: it may embed in genuine presentations. First test: `p^3 q r` against
the order-six cube.
