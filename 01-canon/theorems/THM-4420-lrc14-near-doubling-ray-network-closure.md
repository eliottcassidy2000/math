---
id: THM-4420
title: "LRC14 near-doubling ray network closure"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4414 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. Every ternary-unit
  near-doubling triple (1,m,2m+sigma), sigma=+/-1, has a one-ray live-carrier
  dictionary with every third natural address deleted. Each of its three
  degree-zero component-network capacities is at most 6/77, with equality in
  this family only at (1,5,11). Arbitrary triples, entry, synchronization, and
  LRC(14) remain open.
source: root + network_universal / LRC14 continuation session, 2026-09-05
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
related:
  - THM-4409-lrc14-third-sheet-component-network-certificate
  - THM-4413-lrc14-owner-transversality-gap-and-complete-norm-eighteen-empty-atlas
script: 04-computation/lrc14_near_doubling_ray_network_closure_thm4420.py
output: 05-knowledge/results/lrc14_near_doubling_ray_network_closure_thm4420.out
script_sha256: b2535a9b685ca1b068b2027d8d3825a20b17f2530c697e0ccfa5f1f6026f718c
output_sha256: 5e8fc3926fb8d93cacfdc42cd78f1c13feb5136263b343e93be32925b9f2b013
semantic_sha256: d8ff287083ba8fcab1fa5c1dd4b787bea70b4779940cfa41c49a758eed4eeccb
independent_audit_script: 04-computation/lrc14_near_doubling_ray_network_closure_thm4420_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_near_doubling_ray_network_closure_thm4420_independent_audit.out
independent_audit_script_sha256: 337ee03884fe22777a8993e879cc39b1b43c3ffcb7800eb3d0f9763a49ee887f
independent_audit_output_sha256: 50bfd22c0620da0d19717b938e0c285405d9b7c0e3d8eca918d61c4c6ec426b7
independent_semantic_sha256: cefbe0c18d557cc6ec5786e3be317d8f8c2af9fe20fa7143bcf706b89217d798
hash_basis: raw LF bytes
audit: >
  PASS. The exact verifier checks 832 family rows through largest speed 4999,
  every strict endpoint and residue deletion, the six small threshold rows,
  and equality uniqueness. On fifteen hostile controls an independent path
  enumerates the full relation box, constructs all six literal sheet graphs,
  runs exact rational max flow, and compares literal physical mass, for
  1,080,839 explicit gates. Normal and optimized outputs byte-match.
  A clean-room script checks 333,333 family rows through largest speed
  2,000,003 and sixteen definition-level lattice boxes with 1,690,591 gates;
  its normal and optimized outputs also byte-match.
---

# THM-4420 -- LRC14 near-doubling ray network closure

**PROVED ELEMENTARY RELATIVE TO THM-4414 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** This closes an infinite pair of degree-zero local-certificate rays.
It proves no arbitrary three-speed entry statement, fourteen-runner
synchronization, or `LRC(14)`; the conjecture remains **OPEN**.

## 1. Statement

Let `sigma` be `+1` or `-1`, and suppose

```text
w=(1,m,2m+sigma)                                        (1)
```

is sorted, distinct, positive, odd, and ternary-unit. Equivalently,

```text
sigma=+1:  m=6t+5, t>=0;
sigma=-1:  m=6t+1, t>=1.                                (2)
```

Let `E_i(w)` be the three exact degree-zero raw projection capacities from
THM-4414. Then, for every `i`,

```text
E_i(w)<=6/77.                                           (3)
```

Equality in any coordinate occurs in this family only at

```text
w=(1,5,11),                                             (4)
```

where all three capacities and the physical comb mass equal `6/77`.

More precisely, the complete live raw-carrier set is

```text
C=k(-sigma,-2,1),   0<|k|<=K,   3 does not divide k,    (5)

K=strict_floor(3(m+1)/14),   sigma=+1,
K=strict_floor(3m/14),       sigma=-1,                  (6)
```

where `strict_floor(x)` is the greatest integer strictly below `x`. Hence its
cardinality and a simultaneous upper bound for all three projections are

```text
N=2(K-floor(K/3)),
E_i(w)<=3N/[7(2m+sigma)]<=6/77.                         (7)
```

The strict floor in `(6)` is essential: a carrier on the endpoint has a zero
roof facet and is not live.

## 2. Integer-width collapse to one relation ray

Let `C=(C_1,C_2,C_3)` be live. The kernel equation is

```text
m(C_2+2C_3)=-(C_1+sigma C_3).                           (8)
```

Positivity of the first and third THM-4413 roof margins gives

```text
|C_1|<(3/14)(3m+sigma),
|C_3|<(3/14)(m+1).                                      (9)
```

Consequently

```text
|C_2+2C_3|
 <(3/14)(4+(sigma+1)/m).                               (10)
```

For `sigma=-1`, the right side is `6/7`. For `sigma=+1`, condition `(2)`
gives `m>=5`, so it is `6/7+3/(7m)<1`. The left side is integral. Therefore

```text
C_2=-2C_3,       C_1=-sigma C_3,                        (11)
```

and every live carrier lies on `(5)`.

Conversely, every `k(-sigma,-2,1)` lies in `ker(w)`. Its three strict support
conditions are

```text
|k|<(3/14)(3m+sigma),
2|k|<(3/14)(2m+sigma+1),
|k|<(3/14)(m+1).                                       (12)
```

For the plus sign, the last two coincide and are strongest. For the minus
sign, the middle condition is strongest. This gives exactly `(6)`. All three
carrier coordinates are nonzero modulo three exactly when `3` does not divide
`k`, proving the full dictionary `(5)` in both directions.

This is also a lossless natural address. List positive integers not divisible
by three in increasing order, truncate at `K`, and retain one sign bit. The
planar-looking carrier atlas has become a one-dimensional ordinal list with
every third address deleted; no geometric datum needed by `(7)` is lost.

## 3. Counting closes every projection

Among `1,...,K`, exactly `K-floor(K/3)` integers survive the residue deletion.
Reflection `k -> -k` proves the formula for `N` in `(7)`.

For sorted `w`, every summand in every THM-4414 projection is bounded by the
shortest one-sheet cap

```text
q=3/[7(2m+sigma)].                                      (13)
```

There is one raw summand per live carrier, so

```text
E_i(w)<=Nq.                                             (14)
```

It remains to prove

```text
11N<=2(2m+sigma).                                       (15)
```

Writing `A=K-floor(K/3)`, blocks of three give

```text
A<=(2K+2)/3,       N<=(4K+4)/3.                         (16)
```

For the plus family, `(6)` and `(16)` imply

```text
N<2(m+1)/7+4/3<=2(2m+1)/11                             (17)
```

for `m>=19`; the first eligible value there is `m=23`. For the minus family,

```text
N<2m/7+4/3<=2(2m-1)/11                                 (18)
```

for `m>=20`; the first eligible value is `m=25`. The complete finite heads
are

```text
plus:   m=5,11,17,   (K,N)=(1,2),(2,4),(3,4),
minus:  m=7,13,19,   (K,N)=(1,2),(2,4),(4,6).           (19)
```

Direct substitution proves `(15)` in all six cases. Equality occurs only at
the plus row `m=5`. Every later row is already strict in `(17)` or `(18)`, or
in the finite check. THM-4414 evaluates `(1,5,11)` exactly, proving `(3)--(4)`.

## 4. Boundaries and audit

The proof uses all of the following hypotheses.

- The inequalities defining the roof are strict; ordinary floor is wrong at
  such plus rows as `m=41` and `m=83`.
- `k=0` and every multiple of three are excluded by the distinct-owner gate.
- Both orientations `k` and `-k` are retained, producing the factor two in
  `N`.
- The claim bounds each projection, not merely their minimum. This is stronger
  than needed for the universal THM-4414 target but only on the declared rays.

The exact verifier checks the formula through largest speed `4999`. On fifteen
controls it ignores `(5)`, enumerates the complete relation-lattice box,
constructs all six physical sheet graphs, computes exact rational maximum
flows, and compares the raw projection and literal mass calculations. These
paths agree on all controls, including the six finite heads, both endpoint
hostiles, changing active roof facets, and both terminal rows.

An independent implementation imports no repository or primary code. It
checks `333,333` family rows through largest speed `2,000,003` and independently
enumerates sixteen complete relation-lattice boxes. This also confirms that
the count-to-capacity step is specifically THM-4414 equations `(10)--(11)`;
the raw-carrier identification is a load-bearing dependency, not an optional
presentation.

Run

```powershell
python -B 04-computation/lrc14_near_doubling_ray_network_closure_thm4420.py --height 4999
python -B -O 04-computation/lrc14_near_doubling_ray_network_closure_thm4420.py --height 4999
python -B 04-computation/lrc14_near_doubling_ray_network_closure_thm4420_independent_audit.py
python -B -O 04-computation/lrc14_near_doubling_ray_network_closure_thm4420_independent_audit.py
```

The reusable mechanism is

```text
short relation + transverse integer width below one
  -> every live carrier is a multiple of that relation
  -> residue deletion becomes exact one-dimensional counting.             (20)
```

The next extension is `w=(a,m,2m+sigma a)`. There, primitivity only forces
`a | C_2+2C_3`; proving a sharp all-ratio count and the `6/77` ceiling requires
additional work. It is not part of this theorem.
