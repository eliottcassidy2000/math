---
id: HYP-1860
status: OPEN
source: codex-2026-05-31-S392
related:
  - HYP-1841
  - HYP-1842
  - HYP-1854
  - HYP-1857
  - HYP-1858
  - HYP-1859
---

# HYP-1860: n=16 LRC counterexamples must satisfy an antipodal fan certificate

## Statement

At denominator `n=16`, a primitive open-cover counterexample must satisfy two
simultaneous certificates:

1. the antipodal pair quotient

```text
even-forbidden OR (odd-forbidden AND shifted-odd-forbidden)
```

must have no open or boundary gap on `[0,1/2]`;

2. any maximal dyadic endpoint branch closed by lower speeds must either use
the nine-speed scaled odd fan

```text
(v/2) plus (v/32)*{1,3,5,7,9,11,13,15}
```

or pay at least as much endpoint debt elsewhere.

The conjectural proof is that these two certificates are incompatible once
primitivity is imposed: the fan is imprimitive at scale `v/32`, while the
antipodal quotient needs odd double-coverage on every pair not already covered
by even speeds.

## Proved / Exact Branch

For the antipodal quotient, even speeds are invariant under `t -> t+1/2`.
An odd speed is one-sided on every antipodal pair, because shifting by `1/2`
moves `vt` by `1/2 mod 1`.  Therefore a full open cover implies that every
half-turn pair is covered by either an even speed or by odd coverage on both
sides of the pair.

For the local dyadic fan, S392 verifies exactly for `v=16,32,64,128,256,512`
that the nine-speed fan covers all `2v` endpoints of speed `v`.  For
`v>=32`, the normalized fan is stable:

```text
(16,1,3,5,7,9,11,13,15).
```

Its gcd is `v/32`, so a large maximal fan is not primitive by itself.

## Evidence

`lrc_n16_antipodal_fan_proof_search_s392.py` adds three tests beyond S390.

Antipodal-pair audits:

```text
initial 1..15:
  pair_measure=1/2, pair_gap=0, pair_boundary_witnesses=4

drop 15 add 16:
  pair_gap/th=0.038462, boundary=16

best 8-ladder:
  pair_gap/th=0.007576, boundary=80

fan64 plus odd breakers:
  pair_gap/th=0.202083, boundary=36

fan128 plus small breakers:
  pair_gap/th=0.101042, boundary=86
```

One-gate antipodal scan:

```text
{1,...,15} - {r} + {16q}, 1<=r<=15, 1<=q<=32
```

All `480` rows had positive pair gaps.  The closest rows had
`pair_gap/th=0.027344`.

Maximal fan audit:

```text
v=16:  fan covers 32/32 endpoints, gcd=1
v=32:  fan covers 64/64 endpoints, gcd=1
v=64:  fan covers 128/128 endpoints, gcd=2
v=128: fan covers 256/256 endpoints, gcd=4
v=256: fan covers 512/512 endpoints, gcd=8
v=512: fan covers 1024/1024 endpoints, gcd=16
```

The exact candidate audits still have ordinary positive gap, unprotected
endpoints, positive pair gap, and terminal `coreE=0` for every non-initial row.

## Predictions

1. The antipodal quotient should be a useful first-pass branch-and-bound
   certificate for `n=16`, filtering gated candidates before full endpoint
   peeling.
2. Any near-counterexample with a large maximal dyadic speed should contain an
   approximate scaled fan plus a small number of gcd-breaker speeds.
3. The proof target is an incompatibility lemma: gcd-breakers needed for
   primitivity cannot also supply enough odd double-coverage in the antipodal
   quotient without creating private endpoint debt.
4. A cycle-first search for `n=16` should track two labels on each arrow:
   dyadic fan residue and antipodal side.

## Sources

- `04-computation/lrc_n16_antipodal_fan_proof_search_s392.py`
- `05-knowledge/results/lrc_n16_antipodal_fan_proof_search_s392.out`
- `07-reflections/lrc-n16-antipodal-fan-proof-search-s392.md`
- HYP-1857
- HYP-1858
