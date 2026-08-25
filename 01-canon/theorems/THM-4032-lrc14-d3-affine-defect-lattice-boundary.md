---
id: THM-4032
title: "LRC(14) d=3 affine defect-lattice boundary"
status: >
  PROVED RELATIVE TO THM-4024/4004 + VERIFIED-EXACT + INDEPENDENTLY
  REFEREED. At the remaining d=3,c_3=8 equality boundary, full spoilage of
  the three labelled lifts is equivalent to a bounded support-three defect
  circuit with exact gcd, mod-3, and strict-window sidecars.
  Certificate-negative rows close. Certificate-positive rows remain
  unresolved because the test forgets the divided pack's safe-phase set.
  This does not prove LRC(14).
source: root + tournament_sequence_bridge / sequence continuation, 2026-08-24
audit: >
  PASS. Direct rational wall-cell subdivision, bounded affine-centre
  enumeration, and bounded residue-defect enumeration agree on all 560
  distinct triples in {1,...,23} coprime to three. The audit checks 3,627
  semantic gates, the first false pairwise-gcd converse, the minimal distinct
  geometric hostile, and a typed ten-pack hostile. An independently written
  literal-open-interval checker agrees on 1,540 triples through 32, including
  1,065 positives and 14,334 gates, and isolates separate circuit and mod-3
  hostiles. All normal and optimized outputs are byte-identical.
depends_on:
  - THM-4024-lrc14-complete-divisor-incidence-envelope
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
related:
  - THM-4030-lrc14-d4-affine-defect-lattice-boundary
  - THM-2059-crt-fiber-product-phase-packet
script: 04-computation/lrc14_d3_affine_defect_lattice_boundary_thm4032.py
output: 05-knowledge/results/lrc14_d3_affine_defect_lattice_boundary_thm4032.out
script_sha256: be09b2ff98503ce470a4cf233098c39ec72fd74c193d2f57e3a488f5ad052a43
output_sha256: 60909f655d86293216b001377b222431a36bb8a66b5cd1187b014a96c21a91ca
independent_script: 04-computation/lrc14_d3_affine_defect_independent_referee_thm4032.py
independent_output: 05-knowledge/results/lrc14_d3_affine_defect_independent_referee_thm4032.out
independent_script_sha256: 2935a6976d19b5ef99bdaae45be7b77e905e4f9edb8946af06af7441b843ffb2
independent_output_sha256: e7b72cded0da828beaf8c237341bb13dca8eba22741fa3b027ac8a7085ec4cc5
hash_basis: raw LF bytes
---

# THM-4032 -- the d=3 affine defect-lattice boundary

**PROVED RELATIVE TO THM-4024/4004 + VERIFIED-EXACT + INDEPENDENTLY
REFEREED.** The theorem converts the order-three selector obstruction into a
finite exact Diophantine test for each exception triple. It is necessary for
a counterexample, not sufficient. LRC(14) remains open.

## 1. Typed boundary and inheritance

Retain THM-4024's `d=3,c_3=8` boundary in the rank-eleven `11+2` branch.
The eight body speeds divisible by three and the divided pair give a cited
ten-speed pack. The three remaining body speeds are distinct positive integers

```text
a,b,c,                 3 does not divide abc.                 (1)
```

All thirteen LRC speeds are pairwise distinct. For a pack phase
`y in R/Z`, its three lifts are

```text
x_j=(y+j)/3,                    j=0,1,2.                       (2)
```

Call `y` fully spoiled by the exceptions when every label `j` satisfies at
least one of

```text
||a x_j||<1/14,  ||b x_j||<1/14,  ||c x_j||<1/14.             (3)
```

The closest proved mechanism is THM-4024's exact bad-count equality
`1+1+1=3`. Its canonical hostile is THM-4004's typed order-three selector.
The corrected near miss is to retain only those three counts. The least-used
sidecar is the affine centre lattice together with its nontrivial residue
classes modulo three.

## 2. Exact equivalence

Fix the name `a` for any one exception. The following are equivalent.

1. Some phase `y` is fully spoiled.
2. For one of the two orientations `(a,b,c)` or `(a,c,b)`, there are integers
   `A,B,C` such that, with

   ```text
   N_ab=3bA-3aB+ab,
   N_ac=3cA-3aC+2ac,
   N_bc=3cB-3bC+bc,                                  (4)
   ```

   one has

   ```text
   14|N_ab|<3(a+b),
   14|N_ac|<3(a+c),
   14|N_bc|<3(b+c).                                  (5)
   ```

3. For one orientation there are integers `N_ab,N_ac,N_bc` satisfying `(5)`,

   ```text
   gcd(a,b)|N_ab,  gcd(a,c)|N_ac,  gcd(b,c)|N_bc,     (6)

   N_ab == ab  (mod 3),
   N_ac == 2ac (mod 3),
   N_bc == bc  (mod 3),                              (7)

   cN_ab-bN_ac+aN_bc=0.                              (8)
   ```

Because `a,b,c` are units modulo three, `(7)` makes all three defects
nonzero. It also makes equality in any line of `(5)` arithmetically
impossible: the left side would be nonzero modulo three while the right side
is zero. Thus this branch has an automatic no-wall gap despite originating
from open danger intervals. Quantitatively, every passing pair has overlap
slack at least `1/(42xy)` for its two speeds `x,y`.

The cited ten-speed result supplies at least one safe phase of the divided
pack. If the full row were non-lonely, the exceptions would have to spoil all
three lifts of that phase. Consequently every hypothetical counterexample at
this boundary forces a genuine support-three relation `(8)` with

```text
|N_ab|+|N_ac|+|N_bc| < 3(a+b+c)/7.                    (9)
```

Absence of a defect triple satisfying `(5)--(8)` closes the row. Existence
says only that some abstract selector phase can be spoiled; it does not say
that a safe phase of the divided pack is spoiled, and it does not imply
non-loneliness.

## 3. From three labels to three intervals

Every exception in `(1)` has order three on the lifts and can spoil at most
one label. Full spoilage therefore forces the three bad-label sets to be
singletons partitioning `{0,1,2}`.

Rotate the labels so that the one killed by `a` is zero. Swapping `b,c` if
necessary leaves one real lift `z` with

```text
a dangerous at z,
b dangerous at z+1/3,
c dangerous at z+2/3.                                 (10)
```

For `D_w={z:||wz||<1/14}`, this is

```text
z in D_a intersect (D_b-1/3) intersect (D_c-2/3).     (11)
```

Conversely, `y={3z}` reconstructs a fully spoiled three-lift phase; the
integer part of `3z` only rotates the labels. This proves that the two
orientations exhaust all six label assignments after the designated `a` is
used as origin.

## 4. Interval-to-lattice proof

Lift the selected open components in `(11)` to the real line. Their centres
and radii are

```text
c_a=A/a,        R_a=1/(14a),
c_b=B/b-1/3,    R_b=1/(14b),
c_c=C/c-2/3,    R_c=1/(14c).                          (12)
```

Three real intervals intersect iff they intersect pairwise. Direct
subtraction gives

```text
c_a-c_b=N_ab/(3ab),
c_a-c_c=N_ac/(3ac),
c_b-c_c=N_bc/(3bc),                                  (13)
```

so the three pairwise radius tests are exactly `(5)`. This proves `1 iff 2`.

Expanding `(4)` gives `(7)` and `(8)`, while every common divisor in `(6)`
divides the corresponding formula. Hence `2 implies 3`.

For the converse, put

```text
K_ab=(N_ab-ab)/3,
K_ac=(N_ac-2ac)/3,
K_bc=(N_bc-bc)/3.                                    (14)
```

The residue conditions make these integers. Since every relevant gcd is
coprime to three, `(6)` passes through division by three and gives

```text
gcd(a,b)|K_ab,  gcd(a,c)|K_ac,  gcd(b,c)|K_bc.        (15)
```

Equation `(8)` becomes

```text
cK_ab-bK_ac+aK_bc=0.                                  (16)
```

It remains to solve

```text
bA-aB=K_ab,
cA-aC=K_ac,
cB-bC=K_bc.                                           (17)
```

For a prime power `p^alpha||a`, write the full valuations
`beta=v_p(b)` and `gamma=v_p(c)`. Each of the first two congruences for
`A` is separately soluble by `(15)`. Suppose `beta<=gamma`; the `b`
congruence is the stronger one. If it holds, use the exact identity

```text
b(cA-K_ac)=c(bA-K_ab)-aK_bc.                           (18)
```

The first term on the right has valuation at least `alpha+gamma`, while
`gcd(b,c)|K_bc` makes the second have valuation at least
`alpha+beta`. Division by `b`, with its full valuation `beta`, proves
the `c` congruence modulo `p^alpha`. The case `gamma<beta` is symmetric.
CRT supplies `A`; then

```text
B=(bA-K_ab)/a,         C=(cA-K_ac)/a
```

are integers, and `(16)` gives the third equation. Thus `3 implies 2`.

This prime-local step is the exact local/global lesson: the single rational
circuit is not a substitute for the three gcd and residue sidecars, while the
sidecars alone do not enforce a common affine lift.

## 5. Finite decision and cheap gates

For a fixed exception triple it is enough to enumerate

```text
0<|N_ab|<=(3(a+b)-1)//14,
0<|N_ac|<=(3(a+c)-1)//14,
0<|N_bc|<=(3(b+c)-1)//14,                              (19)
```

retaining `(6)--(8)`. Equivalently, simultaneous translation

```text
(A,B,C)->(A+a,B+b,C+c)
```

is a complete gauge, and pairwise overlap allows the finite box

```text
0<=A<a,  0<=B<=2b,  0<=C<=2c.                         (20)
```

Indeed `0<=A<a` gives `c_a in [0,1)`; overlap gives
`|c_b-c_a|,|c_c-c_a|<1/7`, whence `0<B/b<31/21<2` and
`0<C/c<38/21<2`.

Every defect is nonzero and gcd-divisible, so a fully spoiled phase requires

```text
3(a+b)>14gcd(a,b),
3(a+c)>14gcd(a,c),
3(b+c)>14gcd(b,c).                                    (21)
```

These gates are not sufficient. The first distinct positive triple passing
all three but failing full spoilage is `(1,4,7)`; no defect satisfies every
sidecar. More sharply:

```text
omit mod 3: (a,b,c)=(2,7,8), N=(-1,-2,-3), circuit residual 0;
omit circuit: (a,b,c)=(1,4,7), N=(1,-1,-2), residual 9. (22)
```

Both weakened predicates pass their other stated conditions, while neither
profile is semantically spoiled. Thus the incoming local residue and global
circuit coordinates are separately load-bearing. The decision is finite for
each triple, not uniformly finite over the unbounded `d=3` branch. A height,
ratio, owner, or pack-safe-phase bound is still needed.

## 6. Sharp controls

In the pairwise-distinct exception setting, the smallest possible maximum
speed of a fully spoiled triple is five:

```text
(a,b,c)=(1,4,5),  y=23/112,
bad labels=({0},{2},{1}).                              (23)
```

An affine certificate for the same profile is

```text
(A,B,C)=(0,1,3),
(N_ab,N_ac,N_bc)=(1,1,-1),
5-4-1=0.                                               (24)
```

The typed divided-pack hostile is sharper about scope. Take

```text
body=(6,9,12,15,18,21,24,27,1,5,11),
pair=(3,30).                                           (25)
```

After division by three, the ten-speed pack is `H={1,...,10}`. Its safe
phase `y=2/11` has exception bad-label partition

```text
(1,5,11): ({0},{1},{2}).                               (26)
```

Yet the full row has clearance exactly `1/14` at `x=1/14`. Thus `(26)` is a
hostile to arbitrary selector choice, not an LRC counterexample.

## 7. Preservation, weighted tournament, and scope

The map from three labelled lifts to the defect triple preserves cyclic label
assignment, strict endpoints, pairwise interval margins, and exact
common-phase compatibility. It destroys the safe-phase set `G(H)` of the
divided pack, the other body geometry, height filters, and which cited pack
witness was selected.

There is an intrinsic weighted tournament on the three selected interval
centres: orient each edge from the lower centre to the higher and decorate it
by

```text
|c_u-c_v|/(R_u+R_v).
```

A common-centre realization forces a transitive orientation, and full
spoilage is exactly the condition that all three weights are below one. The
ordinary tournament forgets those thresholds and therefore is not an
equivalent certificate. Even the edge-weight order is insufficient. For
`(a,b,c)=(1,4,5)`, both centre triples below have order `b<c<a` and the same
ranking of normalized edge weights:

```text
(A,B,C)=(0,1,3): weights=(14/15,7/9,14/27), fully spoiled;
(A,B,C)=(0,0,2): weights=4(14/15,7/9,14/27), not spoiled. (27)
```

The absolute threshold decoration is load-bearing. This weighted tournament
is a controller for the affine selector only; it still loses `G(H)`.

The needed next sidecar is the intersection of the certificate-positive
spoiled-phase set with `G(H)`. Certificate absence is a loneliness
certificate; certificate presence has no non-loneliness consequence.

This theorem is relative to THM-4004's cited ten-speed pack and THM-4024's
`d=3,c_3=8` reduction. It does not close every triple uniformly, the remaining
`d=2` equality boundary, or LRC(14).

## 8. Replay

```text
python3 -B 04-computation/lrc14_d3_affine_defect_lattice_boundary_thm4032.py
python3 -B -O 04-computation/lrc14_d3_affine_defect_lattice_boundary_thm4032.py
python3 -B 04-computation/lrc14_d3_affine_defect_independent_referee_thm4032.py
python3 -B -O 04-computation/lrc14_d3_affine_defect_independent_referee_thm4032.py
```

Both runs reproduce the frozen raw-LF output. **QED.**
