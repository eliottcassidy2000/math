---
id: THM-3613
title: "Exponent-two three-by-four size-seven ray parity gate"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.  The seven
  sumset-seven words in the exponent-two 3 x 4 support atlas are primitive
  one-dimensional rays.  At scales n>=4 their scalar-arm placements and
  singleton sign graphs stabilize, and the UFD exponent of the negative arm
  coefficient depends only on n modulo two.  The exact gate rejects 31 of 64
  tail placements at even scales and 18 of 64 at odd scales; scales 1,2,3 are
  exhausted separately.  This is a necessary nonentry gate only and does not
  assert that any surviving placement is a Darboux pair.
source: kps-s188 / THM-3609 continuation, 2026-08-21
audit: >
  Author proof and exact companion verification.  The companion hash-pins
  THM-3606, proves rank four for each collision-equality system, identifies its
  primitive integral kernel generator, enumerates all exceptional scales,
  certifies affine sign stabilization for every tail candidate, bridges the
  symbolic tail ledger record-for-record to the original integer evaluator in
  both parity classes, and agrees byte-for-byte under ordinary and optimized
  Python.  Independent hostile audit remains pending.
depends_on:
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
related:
  - THM-3609-three-by-four-size-nine-euler-factor-nonentry
script: 04-computation/jc2_three_by_four_size_seven_ray_parity_gate_thm3613.py
output: 05-knowledge/results/jc2_three_by_four_size_seven_ray_parity_gate_thm3613.out
script_sha256: 2b4423460e89696f95b5a046affeeaf36920a59f529ee08762440cbd0260daed
output_sha256: 0b6600e995c310c5279242e2d3b69c793989b488ed697eed17117936fc7de9b8
semantic_sha256: 3e30e66899f6c1f49dd82ad41fa1a19a73f4fe247fa8f0b30fbcfd4ab7fe0ac6
hash_basis: LF-normalized bytes
---

# THM-3613 -- exponent-two three-by-four size-seven ray parity gate

**PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.**  This theorem
is a global arithmetic refinement of the scalar-arm/singleton gate on the
smallest unresolved additive layer.  It is not a bounded search: after three
explicit exceptional scales, the proof treats every positive scale at once.

All rings are over `C`.  Let `Sigma in C[b]` be squarefree with
`deg Sigma>=2`, and consider a normalized exponent-two Darboux pair with three
retained homogeneous pieces against four as in THM-3606.  Write its positive
support gaps as

```text
(X,Y,U,V,W),
A=(0,X,X+Y),                 B=(0,U,U+V,U+V+W).          (1)
```

We retain exactly the same necessary gates as THM-3606: a collision fibre is
chosen as the scalar row, at least one address there has weights `(-2,1)` or
`(1,-2)`, every other singleton row passes the same-sign/zero test, and the
scalar fibre is not rectangle-exposed.

## 1. The seven primitive rays

The sumset-seven words are precisely the following seven oriented words.  The
table gives the primitive gap generator `rho`.

```text
word    rho=(X,Y,U,V,W)
W002    (1,1,1,1,2)
W003    (1,2,1,1,1)
W004    (1,2,1,1,2)
W005    (2,1,2,1,1)
W006    (1,1,1,2,1)
W007    (2,1,1,1,1)
W008    (1,1,2,1,1).                                    (2)
```

For each word, the equalities imposed by its collision fibres have rank four
in the five gap variables, and `rho` spans their kernel.  Since every `rho`
is primitive and has a coordinate equal to one, every positive integral point
in that equality cone is uniquely

```text
(X,Y,U,V,W)=n rho,                         n in Z_{>0}.  (3)
```

Positive scaling preserves all strict fibre-order inequalities, so `(3)` is
the full integral word cone, not merely a sampled subfamily.

## 2. Why only parity survives

Fix a word and a scalar-arm placement.  Anchoring an address at `(-2,1)` or
`(1,-2)` makes every displayed weight an affine function

```text
r(n)=alpha n+c,                              alpha,c in Z.          (4)
```

The companion enumerates the complete finite list of anchor choices.  For
every singleton weight in every candidate that survives at scale four, a
positive-slope form is already positive at `n=4`, a negative-slope form is
already negative there, and a zero-slope form is constant.  Therefore no sign
can change for any `n>=4`.  Candidate survival, the singleton graph, and the
unique eligible scalar-arm address are consequently constant throughout the
tail.  Direct evaluation at `n=4` and `n=5` agrees record-for-record with this
symbolic affine ledger, covering both parity classes.

Let `v` be the negative coefficient at the unique eligible tail address.  Its
weight is `-2`.  In the connected same-sign singleton component containing
`v`, the zero Wronskian equations and unique factorization express all
coefficients as powers of one base polynomial.  If the component weights are
`r_1,...,r_s`, the exponent of the coefficient at `v` is

```text
R_v = 2/d,               d=gcd(2,|r_1|,...,|r_s|).       (5)
```

Thus `R_v` is one exactly when every component weight is even, and is two
otherwise.  By `(4)`, this condition depends only on `n mod 2`.  If `R_v=2`,
the negative coefficient is a square up to a nonzero scalar, so all of its
zeros have even order.  It cannot provide the simple zero required at any
root of the squarefree arm polynomial `Sigma`.  Such a placement is therefore
impossible.  At an exceptional scale with several eligible addresses, the
same conclusion is used only when every eligible negative coefficient has
exponent two.

## 3. Exact classification

Each entry below is `total/rejected`, where `total` is the number of placements
surviving the inherited THM-3606 gates at that fixed ray point and `rejected`
is the number eliminated by `(5)`.

```text
word      n=1     n=2      n=3      even n>=4    odd n>=5
W002      0/0     3/1      4/1         6/4          6/2
W003      1/1     6/4      6/0         8/4          8/1
W004      3/1    12/5     10/1        12/5         12/2
W005      4/3     7/5     10/3        12/5         12/3
W006      2/1     7/5      8/0        12/5         12/2
W007      2/2     4/4      6/3         8/4          8/4
W008      0/0     2/1      4/2         6/4          6/4
-----------------------------------------------------------------
tail total                              64/31         64/18.       (6)
```

In particular, the inherited gates plus the parity obstruction leave no
placement at `W002,n=1`, `W003,n=1`, `W007,n=1,2`, or `W008,n=1`.  On every
even tail point the new gate removes 31 of 64 placements; on every odd tail
point it removes 18 of 64.

The exact companion derives the seven rays from the full atlas, verifies the
rank-four kernels, enumerates all scales `1,2,3`, performs the universal affine
tail analysis, records every singleton component and parity exponent, and
checks the symbolic/integer bridge in both residue classes.  Reproduce with

```bash
python3 04-computation/jc2_three_by_four_size_seven_ray_parity_gate_thm3613.py
python3 -O 04-computation/jc2_three_by_four_size_seven_ray_parity_gate_thm3613.py
```

The boundary is important.  Table `(6)` classifies a necessary
scalar-arm/singleton obstruction; it does not solve the remaining collision
differential equations.  A count such as `6/4` leaves two candidates, not two
Darboux pairs.  The theorem proves no planar Jacobian counterexample and no
case of the planar Jacobian conjecture.

**QED.**
