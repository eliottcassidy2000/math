---
id: THM-3478
title: "Dual seven-exit block closure through r=2603"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Two
  separately pinned coefficient/Newton implementations give identical full
  degree traces and close all 38 seven-exit residuals in 2501<=d<=2600.
  Two formula-only implementations of THM-3475 then agree that divisor-place
  pair ledgers close d=2601,...,2605 and leave d=2606 with exactly the four
  necessary degrees 521,1042,1563,2084.  Consequently every exact-support
  quadratic three-moment window beginning at 1<=r<=2603 contains a nonzero
  moment.  This is a bounded exact-support theorem, not SFC(3) or FC(3).
audit: >
  The fixed block reconstructs the inherited seven-exit universe, all three
  resonant rows, full zero-root-aware barcodes at the ordered 15-prime bank,
  and every progressive survivor set.  Primary and independent routes agree
  on four chunk digests and one global digest, retain planted linear factors,
  and have no final survivor.  On the next six rows, primary and independent
  digit compilers agree on every divisor polygon and local degree set.  The
  stored transcript records the completed route outputs; the expensive final
  wrapper was syntax/AST/render audited but was not rerun end-to-end after
  packaging the already completed routes.
source: root/factorial-jacobian-alternation/2026-08-15
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
  - THM-3201-seven-exit-factorial-newton-euclidean-closure-through-r2403
  - THM-3475-factorial-all-divisor-digit-polygon-and-pair-ledger-compiler
related:
  - THM-3467-factorial-seven-exit-newton-barcode-extension-after-r2498
  - THM-3474-factorial-binary-submask-polygon-and-prime-power-reset-families
script: 04-computation/factorial_dual_seven_exit_block_closure_thm3478.py
output: 05-knowledge/results/factorial_dual_seven_exit_block_closure_thm3478.out
script_sha256: 65c651ce611317bc44de472c79d246c0a1bcb90ddfec8145cdd0ef33304073e7
output_sha256: 0bda5199917f16e59f0c6b3319cabe6184e7d93dc2556752d15166b71f653935
hash_basis: raw bytes
---

# THM-3478 -- dual seven-exit block closure through r=2603

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and scope

Let

```text
L(t^m)=m!,                    q(t)=a+bt+ct^2,               (1)
```

with `abc!=0`.  For every integer

```text
1<=r<=2603,                                                (2)
```

the three moments

```text
L(q^r),                 L(q^(r+1)),                 L(q^(r+2)) (3)
```

cannot all vanish.

This is an exact-support `{0,1,2}` result.  It does not cover a missing
coefficient, translated or arbitrary support, all of `SFC(1)`, `SFC(3)`, or
the three-variable Factorial Conjecture `FC(3)`.

## 2. Resonance pair, Euclidean row, and inherited exits

Put

```text
d=r+2,                    N=d-1.                           (4)
```

THM-3124 reduces a hypothetical bad window to the normalized resonance and
a common root of

```text
P=A_(d-2)^(d),            Q=A_(d-1)^(d),
A_n^(d)(v)=L((d-t+v t^2)^n).                              (5)
```

With `n=d-2`, THM-3152's first full integral Euclidean row is

```text
R_1=(2n-1)(Q-2(n+1)(2n+1)vP)+2d(n+1)P,                   (6)
```

and

```text
gcd_Q(P,Q)=gcd_Q(P,R_1).                                  (7)
```

THM-3201 supplies the conclusion through `d=2500` and the seven uniform
exits

```text
d prime; d-1 a prime power; d-j prime for 2<=j<=6.        (8)
```

A **seven-exit residual** is a row on which none of (8) applies.  The task is
therefore to exclude a rational common factor at every residual row from
`d=2501` onward.  A common complex root would give a nonconstant rational
gcd, so rational coprimality is the required target.

## 3. The exact 38-row fixed universe

On the complete interval `2501<=d<=2600`, applying the filters (8) in their
displayed order gives the progressive counts

```text
(89,78,69,60,52,44,38).                                   (9)
```

The last `38` rows, and no others, are

```text
2501,2502,2510,2511,2512,2513,2514,2515,2516,2517,
2518,2519,2520,2528,2529,2530,2538,2564,2565,2566,
2567,2568,2569,2570,2571,2572,2573,2574,2575,2576,
2577,2578,2586,2587,2588,2589,2590,2600.                 (10)
```

For each row, both companions construct `P,Q,R_1` and intersect THM-3152's
zero-root-aware necessary common-factor degree barcode over the ordered bank

```text
S=(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47).            (11)
```

The resulting positive-degree intersection is empty on every row:

```text
intersection_(p in S) D_p(P,Q,R_1) intersection Z_(>0)
 =empty.                                                   (12)
```

By the degree-barcode theorem, (12) excludes a nonconstant rational factor
common to the three rows.  Equation (7) then gives `gcd_Q(P,Q)=1`.

The first prime that empties each progressive trace is grouped below.  This
is the complete killer map, not a sample:

```text
p=2:  2565,2569,2577;
p=3:  2513,2515,2517,2519,2529,2567,2571,2573,2575,2587,2589;
p=5:  2501,2502;
p=7:  2511,2512,2514,2516,2518,2520,2528,2530,2600;
p=11: 2510;
p=13: 2538,2564,2566,2570,2574,2576,2578,2586,2588,2590;
p=17: 2568,2572.                                          (13)
```

Thus the killer histogram is

```text
((2,3),(3,11),(5,2),(7,9),(11,1),(13,10),(17,2)).        (14)
```

The nontrivial intermediate packets matter.  In particular, THM-3475's
divisor-only pair compiler leaves exact addresses at
`d=2516,2564,2571,2576,2586`; the present three-row fixed traces eliminate
all five.  Hence (12) is not the result of an observer that was already
vacuous before the fixed-bank pass.

## 4. Independent coefficient replay

The primary route uses the canonical FLINT coefficient engine and Fraction
hulls.  The independent route uses the separately pinned determinant/hull
engine.  Each rebuilds the exact universe (10), the three rows, every local
barcode, and the full progressive degree trace.  Both retain degree one for
planted factors `v` and `v+1`, and both return no survivor.

The four deterministic chunks and their common semantic digests are

```text
[2501,2516] d8ddf55e638c503d653927422899c4bddc381ac2bdcc0c53c9248bf161679da6
[2517,2565] f3575666ce9261b9186bbfe6bcc961dcd1600c484fe36e35e310084c45e99b32
[2566,2575] e58d91f4c2637b0555e079f5abc47703049f29e2c7e9258f792157d4df77e1cb
[2576,2600] 8ec99dd199f66f764a163e8fb52d562b151a939a124e5a1efc8c4c6404fe95f2.
                                                                    (15)
```

Under the frozen schema

```text
tuple((d,tuple((prime,sorted_post_intersection),...)),...), (16)
```

the two global records are byte-for-byte equal and have digest

```text
28bdc29a6dadfc941ccab7b8eddafd77bd6bbcec9a2574455d4ba3a6dd439b9f. (17)
```

This is a dual finite-exact coefficient result.  The agreement does not turn
the bounded scan into an all-height structural theorem.

## 5. Digit-only extension through d=2605

Both inherited implementations independently verify that each of

```text
d=2601,2602,2603,2604,2605,2606                         (18)
```

is again a seven-exit residual.  For these six rows, no large coefficient
construction is needed.  THM-3475 proves that, at every prime `p|N=d-1`, the
complete necessary pair barcode of `P,Q` is the intersection of two explicit
base-`p` digit polygons.  Two formula-only implementations independently
compute all divisor polygons, capacities, denominators, and progressive
degree sets.

Their complete boundary summary is

```text
d      factorization of N=d-1       first empty divisor / final address
2601   2^3*5^2*13                    p=13: empty
2602   3^2*17^2                      p=17: empty
2603   2*1301                        p=1301: empty
2604   19*137                        p=137: empty
2605   2^2*3*7*31                    p=7: empty
2606   5*521                         {521,1042,1563,2084}.  (19)
```

The two independently shaped records agree exactly.  With the frozen schema

```text
tuple((d,tuple((p,blocks,sorted_post),...),sorted_final),...), (20)
```

their common semantic digest is

```text
602f0fac54c487114457683a3264d2a095a7f048f9a0b3769332d3ead0e61289. (21)
```

At `d=2606`, the `521`-adic block has reduced denominator `521` and capacity
`2084`; after the `5`-adic intersection it leaves exactly the four multiples
shown in (19).  They are necessary degree addresses, not constructed
factors.  This is the first row not closed by this divisor block.

## 6. Consequence and exact boundary

THM-3201 closes through `d=2500`.  Sections 3--4 close every residual through
`d=2600`, while (8) closes all nonresidual rows.  Section 5 closes the next
five residuals.  Therefore every `3<=d<=2605`, equivalently every
`1<=r<=2603`, is closed, proving (2)--(3).

The first unclosed row for the present block is

```text
d=2606,             r=2604,
necessary degrees={521,1042,1563,2084}.                  (22)
```

No assertion is made that one of these addresses is realized by a factor or
that `r=2604` supports a bad moment window.

## 7. Reproduction and audit boundary

The orchestration companion pins all four dependency hashes, has no Python
`assert` nodes, and emits the complete universes, schemas, digests, killers,
controls, and boundary packet.  Reproduce both expensive routes plus the
cheap dual formula extension with

```bash
python3 04-computation/factorial_dual_seven_exit_block_closure_thm3478.py
```

or rerun one coefficient route at a time with

```bash
python3 04-computation/factorial_dual_seven_exit_block_closure_thm3478.py \
  --fixed-engine primary
python3 04-computation/factorial_dual_seven_exit_block_closure_thm3478.py \
  --fixed-engine independent
```

Repeating a command with `python3 -O` exercises the same explicit
`require` gates.  Reference runtimes on the session host were about `6,300`
seconds for the primary route and `4,020` seconds for the independent route,
each with four workers; the dual six-row formula extension took under one
second.

For audit clarity: those two costly routes were completed separately and
gave the exact common records (15)--(17).  After the final wrapper was
packaged, it was syntax-checked, AST-checked, and its rendered transcript was
matched to the stored output; the wrapper itself was not then rerun for
another end-to-end multi-hour pass.  The completed route records and pinned
digests, rather than that packaging step, are the evidence for Section 4.

## 8. Failure boundaries and information contract

- Empty Newton degree intersections are sufficient to exclude a common
  factor; nonempty intersections, including (22), are only necessary
  addresses.
- The fixed-bank closure is bounded and coefficient-exact.  It does not imply
  that the same bank closes later blocks.
- The divisor compiler in Section 5 is structural, but its six-row
  application and the first-survivor assertion are finite-exact.
- Dropping exact support invalidates the resonance pair (5); no
  arbitrary-support or multivariate FC consequence follows.
- The full progressive degree trace, prime labels, capacities, and reduced
  denominators are load-bearing.  A final survivor count alone would not
  distinguish a real obstruction from an inactive observer.

The information contract is

```text
source:    a hypothetical exact-support quadratic three-moment zero window
target:    a rational common factor of P,Q (and the Euclidean row R_1)
map:       resonance reduction -> local Newton degree barcodes
preserved: d, prime, full degree set, zero-root capacity, slope denominator
lost:      factor existence, root residues, arbitrary support, all heights
sidecar:   seven-exit invoice; complete progressive trace; independent route
hostile:   d=2606 with four nonempty necessary degree addresses.             (23)
```

**QED.**
