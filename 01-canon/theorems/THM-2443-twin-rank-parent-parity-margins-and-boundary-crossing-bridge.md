---
id: THM-2443
title: "Twin-rank parent parity, dyadic margins, and the boundary-crossing bridge"
status: >
  PROVED (parity lemma; bridge implication) + FINITE-EXACT
  (census corollaries at centers <= 10^8) + VERIFIED-EXACT
  (independent-path reproduction of THM-2422 eq 36). The ordered
  parent count R(k)=#{(a,b) in K^2 : a+b=k} of a twin rank k is odd
  exactly when k is even and k/2 is a twin rank. Under THM-2422's
  open statement (37), the OEIS A002822 boundary-crossing conjecture
  (J. Seymour, 2026-05-20) is equivalent to the Twin Prime
  Conjecture; unconditionally it is verified here for every
  1 <= V <= 16,666,597. An FFT convolution independently reproduces
  THM-2422's census with zero failures and exhibits growing dyadic
  representation margins (minimum ordered parent count 5,176 in the
  top window). The theorem does not prove statement (37), the Twin
  Prime Conjecture, or any all-n additive law.
source: kind-pasteur-2026-07-26-S131
depends_on:
  - THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry
  - THM-2413-prime-index-affine-drift-and-twin-center-weld
related:
  - THM-2433-operation-fibre-deletion-incidence-and-startup-scar
  - HYP-1994-twin-goldbach-necklace
  - HYP-9025-twin-gap-singular-series-partner-law
external: "OEIS A002822 (comment: Jonathan Seymour, May 20 2026), https://oeis.org/A002822"
script: 04-computation/twin_rank_parity_margins_bridge_thm2443.py
output: 05-knowledge/results/twin_rank_parity_margins_bridge_thm2443.out
script_sha256: 9ec7a9fdfb61f0ec90a188dba8670a744273a6068ba3b32a386d236443fcad69
output_sha256: 3f89d1f2cc3252bb0c5440a1199f6dea77fc4373867590aa5c63cd000b6baa51
hash_basis: working-tree bytes (LF)
---

# THM-2443 -- parent parity, margins, and the boundary-crossing bridge

**PROVED + FINITE-EXACT + VERIFIED-EXACT** as itemized below.

Typing per MISTAKE-268: `C_all = A014574`, `C_6 = C_all \ {4}`,
`K = C_6/6 = A002822` (twin ranks: `6k-1` and `6k+1` both prime),
ordered increasingly. THM-2422 SS4 defines the distinct parent fibre
`P(k)` and leaves the all-`n` statement

```text
K \ {1,2}  subset  K +_distinct K                            (37)
```

**OPEN**, with the census (36): every `k in K`, `k >= 3`, `6k <= 10^8`
has distinct twin-rank parents. This theorem adds the parity law of
the full ordered fibre, a quantitative margin profile, and the exact
logical relation between (37) and a very recent external conjecture.

## 1. Parent parity lemma (PROVED)

For `k in K` let `R(k) = #{(a,b) in K^2 : a+b = k}` count ordered
parent pairs. Then

```text
R(k) is odd  <=>  k is even and k/2 in K.                     (1)
```

*Proof.* The swap `(a,b) -> (b,a)` is an involution on the ordered
representation set; its fixed points are the diagonal pairs
`(k/2, k/2)`, of which there is exactly one when `k` is even with
`k/2 in K` and none otherwise. An involution preserves cardinality
parity off its fixed points. QED.

This is the thin-set survival of THM-2422's dense-fibre alternation
law `r(z) = floor((z-1)/2)` (whose period-two boundary behaviour is
the owner's "alternating object"): after thinning from `N` to `K`,
the alternation survives exactly as the parity bit (1) -- the
in-degree parity of the twin-rank summand graph is the indicator of
half-membership. Verified with zero violations on all `440,309`
census targets (output line 6).

Corollary (raw units, THM-2433 companion): center `12` is the unique
center whose only representation is the doubled `6+6` (`R(2)=1`,
checked exactly), and center `4` is summand-inert: `4 + 6K` meets no
center `>= 12` because `4 !== 0 (mod 6)`.

## 2. The boundary-crossing bridge (PROVED)

OEIS A002822 carries the comment (Jonathan Seymour, May 20 2026):

```text
Conjecture: for every positive integer V there exist u, v, w in
A002822 with u <= v <= V < w and u + v = w.  This implies the Twin
Prime Conjecture.                                             (S)
```

**Bridge.** Statement (37) implies: for every `V >= 1` such that some
twin rank exceeds `V`, the least such rank `w` is a Seymour witness.
Hence, given (37),

```text
(S)  <=>  K is infinite  <=>  Twin Prime Conjecture.          (2)
```

*Proof.* Fix `V >= 2` and suppose some rank exceeds `V`; let
`w = min{k in K : k > V}`. Then `w >= 3` (since `1,2 in K`), so by
(37) there are `a < b` in `K` with `a + b = w`. Both are `< w`, and
every rank `< w` exceeding `V` would contradict minimality of `w`;
hence `b <= V` and `(u,v,w) = (a,b,w)` witnesses (S) at `V`. For
`V = 1` take `(1,1,2)`, allowed since (S) permits `u = v`. This
proves (S) for every `V` below the supremum of `K`; if `K` is
infinite, (S) holds for all `V`. Conversely (S) trivially forces
ranks above every `V`, i.e. infinitude, i.e. TPC. QED.

The direction "(S) => TPC" is Seymour's own remark; the content of
the bridge is that **modulo the repo's open (37), Seymour's
conjecture is not merely sufficient for TPC but equivalent to it** --
(37) is the local mechanism that converts bare infinitude into
boundary-crossing sums. Chain: `(37) + TPC => (S) => TPC`.

**FINITE-EXACT corollary.** From the census universe (all ranks with
center `<= 10^8`, max rank `16,666,598`): statement (S) holds for
every `1 <= V <= 16,666,597` (output line 32).

## 3. Independent-path census and dyadic margins (VERIFIED-EXACT)

The companion recomputes every `R(k)` by a real-FFT convolution of
the `K` indicator (exact after rounding; counts `<< 2^52`),
cross-checked against a direct double loop on all `810` ranks
`<= 10^4` with zero mismatches. This is an independent algorithmic
path from THM-2422's forward/reverse scans and reproduces (36) with
zero failures, both in weak (parents may coincide) and distinct form.

The dyadic margin table (output lines 7-30) shows the minimum ordered
parent count per window `[2^j, 2^{j+1})` rising monotonically from
single digits (min `2` at `k = 18, 38, 268`) to `5,176` at
`k = 8,405,513` in the top window, with window means growing like the
Hardy-Littlewood prediction `k/log^4 k` up to bounded singular-series
factors. The self-additive law is not merely unbroken at `10^8`; its
failure margin is growing without interruption at every scale.

Sidecar: the doubling subsequence `{k in K : 2k in K}` has `2,547`
members below `16,666,598/2`, first terms
`1, 5, 110, 135, 355, 425, 555, 565, 975, 1045, 1755, ...`
(output line 31); by (1) these are exactly the ranks whose double has
odd parent count.

## 4. Scope

Nothing here proves (37), TPC, or any all-`n` law. The bridge is a
conditional equivalence; the census statements are exact over the
stated universe only. The gap-partner distribution observed in the
companion's normalized table is recorded separately as HYP-9025.

## 5. Reproduction

```bash
python 04-computation/twin_rank_parity_margins_bridge_thm2443.py
```

Deterministic; runs in a few minutes; every check is a hard
`SystemExit` on failure; final line `ALL CHECKS PASSED`.
