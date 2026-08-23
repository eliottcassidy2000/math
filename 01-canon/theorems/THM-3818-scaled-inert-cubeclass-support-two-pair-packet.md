---
id: THM-3818
title: "Scaled inert cube classes recover the support-two pair packet"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the 5,855
  admissible support-two ratios from THM-3793, the scaled cube value recovers
  the common scale and primitive ratio for every positive common scale: the
  5,855 primitive cube sums have distinct rational-cube classes.  The placed
  primitive covector recovers
  the oriented labelled runner pair and opposite exposed facets, and the
  thirteen ambient residues modulo the pair sum recover exactly the 1/14-
  safety mask on the pair-sum time lattice.  There are 456,690 unoriented
  supports and 913,380 oriented assignments.  This is an all-scale address
  theorem, not an LRC exclusion: the grid can be empty while an off-grid
  loneliness time exists.
source: root + lrc_reversible_address / incoming-signal extension, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (thm3791-hostile-audit, 2026-08-23).  The
  all-scale decoder, 94/5,855 census, support/orientation quotient, literal
  signed covector, eleven-dimensional opposite facets, inclusive residue
  boundary, same-packet off-grid information loss, AP13 hostile, and split/
  exponent-three collision controls were independently rederived.  The audit
  caught the provisional factor-two half-atlas defect recorded as MISTAKE-457;
  the repaired primary and independent assertion-free companions enumerate
  all 913,380 oriented assignments in 5,504,236 and 2,752,607 active gates.
  A separate factor/scan audit proves all 5,855 rational-cube classes distinct,
  checks 17,137,585 pairs, and retains taxicab and arrival hostiles.  A
  17,596-gate sidecar finds one square, eight triangular and no
  square-triangular primitive addresses.  Normal/-O/frozen streams agree.
depends_on:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
related:
  - THM-778-centered-christoffel-endpoint-skew-product
script: 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.py
output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.out
script_sha256: 282e79f865d943d6b7435d3375d21e28de1f7fe2ed438979cc913c1c1ce338ca
output_sha256: 9ca40f84ce8b9ce02cb22197832700699c5674a39578acc50898e276eaa29cdd
semantic_sha256: 03be043c6db1e7679e0bb858ee4e419f482133b0d459c2ee0bbfbf7237b0237f
independent_script: 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_independent_audit_thm3818.py
independent_output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_pair_packet_independent_audit_thm3818.out
independent_script_sha256: aaab79b1a08cad7eed01d23412c1c286662b63679d14b7e4a22dbc96592cc2c2
independent_output_sha256: a859f2863660389caf706c7892cf8f564712519539dea744d11437cb4391aeb5
independent_semantic_sha256: 998891530735604ea2ede50d0296cf003dabdfb777d585431192759c48f6989c
cubeclass_script: 04-computation/lrc14_scaled_inert_cubeclass_rational_class_extension_thm3818.py
cubeclass_output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_rational_class_extension_thm3818.out
cubeclass_script_sha256: e37ae5e8adf3414081e76d56eb01607172a6364703dca76b41ac096c9f6d77c3
cubeclass_output_sha256: 6199147668a19c0a2ed403b14d889caf65078c675068f000856202cb68a861a8
cubeclass_independent_script: 04-computation/lrc14_scaled_inert_cubeclass_rational_class_extension_independent_audit_thm3818.py
cubeclass_independent_output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_rational_class_extension_independent_audit_thm3818.out
cubeclass_independent_script_sha256: 0bb6b60238335695addce0010182f5222da4c4a86c6f59912bd0f36123d063ca
cubeclass_independent_output_sha256: 37aa867d653d7123435f9841cdc772c3cd71c289162006f80aeb6079a17ad7bc
cubeclass_independent_semantic_sha256: 86c481b03158ba3cb7024ef8739fde640be331d9856ce5d4fc0c1b3b4fcc06cb
square_triangular_script: 04-computation/lrc14_cubeclass_square_triangular_sidecar_thm3818.py
square_triangular_output: 05-knowledge/results/lrc14_cubeclass_square_triangular_sidecar_thm3818.out
square_triangular_script_sha256: 0050033a2782e0adf4ab86aaca74db43b7e879caf48efdc0f2fcb7da31a9e963
square_triangular_output_sha256: 673a764dd14e1b06c549a80cf72e747376d0fba6d28f17da3ac778335f3f9b15
square_triangular_semantic_sha256: e80c76bcbbf9a83c6cb87983ababe96906c1160581aaee8079b6325051eb1e53
hash_basis: raw LF bytes
---

# THM-3818 -- the scaled cube address retains one exact finite time grid

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Packet and statement

Let `n=(n_1,...,n_13)` be a row of positive integer speeds.  Choose distinct
labels `i,j` and suppose, after orienting the pair, that

```text
n_i=gp,                 n_j=gq,                 1<=p<q,
gcd(p,q)=1,             p+q<=356.                         (1)
```

Put

```text
D=g(p+q)=n_i+n_j,       M=n_i^3+n_j^3,
a=q e_i-p e_j,          rho_l=n_l mod D.                  (2)
```

Assume every prime divisor of the primitive sum `p+q` is `2 mod 3` and has
exponent at most two.  The common scale `g` is arbitrary.  Then the packet

```text
(M,a,D,rho_1,...,rho_13)                                  (3)
```

recovers exactly:

1. the common scale `g` and primitive ratio `(p,q)`;
2. the oriented labelled pair `(i,j)`, its primitive relation `a`, and the
   opposite exposed facets selected by the signs of `a`; and
3. the complete pair-sum-lattice safety set

```text
K_D(n)={k in Z/DZ : ||k n_l/D||>=1/14 for every l}.        (4)
```

There are exactly `5,855` primitive ratios satisfying these hypotheses.  They
give `5,855*C(13,2)=456,690` unoriented label supports and
`5,855*13*12=913,380` oriented assignments at every allowed scale.  This is an
exact decoder for the displayed packet, not for the other eleven speeds,
their integer quotients by `D`, or any time outside `(1/D)Z/Z`.

## 2. The cube value recovers scale and ratio

First quotient by rational cubes.  On the finite set of all `5,855`
admissible primitive ratios, prime-exponent reduction modulo three and an
independent maximal-cube-divisor scan prove

```text
(p^3+q^3)/(r^3+s^3) in (Q_(>0)^x)^3
  iff (p,q)=(r,s).                                         (5a)
```

Since `M=g^3(p^3+q^3)`, its rational-cube class recovers `(p,q)` and then
`g=(M/(p^3+q^3))^(1/3)`.  Thus finite-atlas decoding works for arbitrary
positive `g`.

There is also a table-free divisor decoder on the stronger physical
all-inert subbranch.  If every prime divisor of `D=g(p+q)` is `2 mod 3`,
apply THM-3793 to

```text
x=gp,                  y=gq,
gcd(x,y)=g,            (x+y)/g=p+q.                       (5)
```

The hypotheses above are precisely its all-inert pair-sum condition and its
cube-free primitive-quotient condition.  Hence `M` has exactly one positive
distinct two-cube representation:

```text
M=(gp)^3+(gq)^3.                                           (6)
```

The representation is effectively recoverable.  For every possible sum
`d=x+y`, necessarily `d|M`, and

```text
xy=(d^3-M)/(3d).                                           (7)
```

Enumerating the positive divisors of `M`, retaining integral positive (7)
with square discriminant, therefore finds the unique pair.  Taking its gcd
recovers `g`; division by `g` recovers `(p,q)`.  This proves the all-scale
table-free decoder without bounding the inert prime powers in `g`.

## 3. The covector recovers placement and the facet pair

The vector `a` in (2) is primitive and satisfies

```text
a dot n=q(gp)-p(gq)=0,                ||a||_1=p+q<=356.    (8)
```

It has one positive coordinate `q` at `i` and one negative coordinate `-p`
at `j`, so its signed support recovers the oriented placement.  On the
projected lonely-runner zonotope of THM-3743, pairing with `a` is the same as
pairing before projection.  Thus its maximum fixes the `i` cube coordinate
at the upper endpoint and the `j` coordinate at the lower endpoint; its
minimum fixes the reverse pair.  The other eleven coordinates are free.
Because the positive speed vector cannot lie in that free-coordinate
subspace, projection preserves dimension eleven, so these are opposite
exposed facets.  No owner or facet information is inferred from `M` alone;
it is retained by `a`.

## 4. Ambient residues recover the pair-sum grid

For every label `l` and residue class `k mod D`, let

```text
r_(l,k)=k rho_l mod D,          0<=r_(l,k)<D.              (9)
```

Then

```text
||k n_l/D||=min(r_(l,k),D-r_(l,k))/D.                    (10)
```

Consequently (4) is given by the integer-only rule

```text
k in K_D(n)
iff 14 min(r_(l,k),D-r_(l,k))>=D for every l.            (11)
```

This proves that the thirteen residues, together with `D`, recover the whole
finite safety schedule.  Adding arbitrary multiples of `D` to any ambient
speed leaves (11) unchanged, which is exactly the information destroyed by
the residue quotient.

## 5. Failure boundary

The schedule is not a loneliness certificate.  In the safe arithmetic
progression row

```text
n=(1,2,...,13),              (p,q,g,D)=(1,4,1,5),        (12)
```

runner `5` lies at the origin at every time `k/5`, so `K_5(n)` is empty.
Nevertheless `t=1/14` has

```text
||l/14||>=1/14                 for 1<=l<=13.              (13)
```

Thus even the complete pair-sum grid can miss a valid off-grid time.  The
cube packet also destroys the other eleven exact speeds, other relations,
global wall chronology, arrival, owner, and every off-grid value.

The arithmetic hypotheses are sufficient, not necessary.  The controls

```text
1729=1^3+12^3=9^3+10^3                       (split prime),
515375=15^3+80^3=54^3+71^3                   (5^3 pair sum) (14)
```

show the two first failed mechanisms when the inert and exponent conditions
are dropped.  No converse or classification outside the stated packet is
claimed.

## 6. Exact verification and scope

The companion independently enumerates every coprime `p<q`, `p+q<=356`,
checks complete coordinate fibres and the divisor decoder, and reconstructs
all 156 oriented covectors for every admissible ratio.  It also tests arbitrary
inert powers in 32 scale controls and compares (11) with direct rational
evaluation.  Normal and `python -O` streams match the frozen output.

The rational-class extension independently factors and scans every primitive
cube sum, then exhausts all `17,137,585` unordered ratio pairs and scaled
taxicab/arrival hostiles.  Its separate audit reaches the same zero-collision
answer.  The square/triangular sidecar exhausts its stated Pell-capped shells;
these classical labels add no LRC predicate.

The precise connection ledger is

```text
source:      support-two Graver branch of THM-3743
target:      scaled singleton two-cube fibre plus pair-sum time grid
map:         (row,pair) -> (M,a,D,ambient residues mod D)
preserved:   scale, primitive ratio, placement, facets, exact grid safety
destroyed:   ambient quotients, other relations, off-grid times, loneliness
sidecar:     full speed row and owner/phase chronology for any LRC use
next test:   intersect the 5,855 ratios with rank-eleven star spaces.       (15)
```

In particular, (3) supplies no exclusion of a hypothetical LRC(14)
counterexample and proves no case of the Lonely Runner Conjecture.  **QED.**
