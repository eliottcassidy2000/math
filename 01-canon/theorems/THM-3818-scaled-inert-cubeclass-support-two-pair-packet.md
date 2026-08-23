---
id: THM-3818
title: "Scaled inert cube classes recover the support-two pair packet"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
  INDEPENDENT HOSTILE AUDIT.  On the 5,855 admissible support-two ratios from
  THM-3793, the scaled cube value recovers the common scale and primitive
  ratio, the placed primitive covector recovers the labelled runner pair and
  opposite exposed facets, and the thirteen ambient residues modulo the pair
  sum recover exactly the 1/14-safety mask on the pair-sum time lattice.
  This is an all-scale address theorem, not an LRC exclusion: the grid can be
  empty while an off-grid loneliness time exists.  The candidate is not a
  proved dependency until independent audit.
source: root + lrc_reversible_address / incoming-signal extension, 2026-08-23
audit: >
  PROVISIONAL EXACT CANDIDATE.  The assertion-free companion exhausts all
  5,855 ratios and 456,690 labelled placements, replays 32 nonprimitive scale
  controls by the sum-divisor decoder, checks direct versus residue schedules,
  and retains split-prime, exponent-three, and off-grid hostiles in 2,764,096
  active gates.  Normal and optimized streams match the frozen output.
  Independent hostile audit is pending.
depends_on:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
related:
  - THM-778-centered-christoffel-endpoint-skew-product
script: 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.py
output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.out
script_sha256: 0bb6b60238335695addce0010182f5222da4c4a86c6f59912bd0f36123d063ca
output_sha256: 37aa867d653d7123435f9841cdc772c3cd71c289162006f80aeb6079a17ad7bc
semantic_sha256: 86c481b03158ba3cb7024ef8739fde640be331d9856ce5d4fc0c1b3b4fcc06cb
hash_basis: raw LF bytes
---

# THM-3818 -- the scaled cube address retains one exact finite time grid

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
INDEPENDENT HOSTILE AUDIT.**

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

Assume every prime divisor of `D` is `2 mod 3` and every prime exponent in
the primitive sum `p+q` is at most two.  Then the packet

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

There are exactly `5,855` primitive ratios satisfying these hypotheses and
`5,855*C(13,2)=456,690` labelled placements at every allowed scale.  This is
an exact decoder for the displayed packet, not for the other eleven speeds,
their integer quotients by `D`, or any time outside `(1/D)Z/Z`.

## 2. The cube value recovers scale and ratio

Apply THM-3793 to

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
decoder without bounding the inert prime powers in `g`.

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
all 78 placed covectors for every admissible ratio.  It also tests arbitrary
inert powers in 32 scale controls and compares (11) with direct rational
evaluation.  Normal and `python -O` streams match the frozen output.

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
counterexample and proves no case of the Lonely Runner Conjecture.  **QED
candidate, pending independent audit.**
