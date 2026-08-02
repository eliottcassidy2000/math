---
id: THM-3090
title: "Affine/projective prime-power handshake and septimal counterfeit"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Equal natural permutation degrees between AGL_2(F_q) on F_q^2 and
  PGL_2(F_r) on P1(F_r), for prime powers q,r, occur only at (q,r)=(2,3)
  and (3,8).  Equal group order, permutation-group equivalence, or
  three-transitivity of the affine action selects only (2,3), yielding the
  exceptional S4 clutch.  The sharp (3,8) counterfeit has degree nine and
  two-transitivity on both sides but orders 432 and 504, two affine triple
  orbits, and a residual C7 scalar fibre of internal contrast rank 54.
  No canonical modular generator, tree, quartic owner, Keller, or LRC map is
  asserted.
source: root-affine-projective-handshake-2026-08-01
depends_on:
  - THM-3088-punctured-projective-direction-algebra-and-exceptional-parity-saturation
related:
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
  - THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness
  - THM-3083-exceptional-binary-point-ternary-direction-s4-tomography-clutch
script: 04-computation/affine_projective_prime_power_handshake_thm3090.py
output: 05-knowledge/results/affine_projective_prime_power_handshake_thm3090.out
script_sha256: e4eead4d76cee45d9895ff88b6e73fc5757269cab23d6ada6198350f849331d0
output_sha256: c8660421ce2bf79be662a62a06658364f4a1f3c672e95dfbafb1e1aa78ed0fe5
hash_basis: LF-normalized bytes
---

# THM-3090 -- affine/projective prime-power handshake and septimal counterfeit

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Exact classification

[THM-3083, the exceptional binary/ternary clutch](THM-3083-exceptional-binary-point-ternary-direction-s4-tomography-clutch.md)
uses the chosen identification

```text
AGL_2(F_2)=S4=PGL_2(F_3)                                  (1)
```

between the four binary affine points and the four ternary projective
directions.  [THM-2996](THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness.md)
proves a different uniqueness statement: only the binary two-dimensional
module admits the displayed faithful `S3` action transitive on its nonzero
vectors.  The present theorem classifies the cross-field permutation-degree
handshake itself.

Let `q,r>=2` be prime powers.  Compare the natural actions

```text
A(q)=AGL_2(F_q) acting on F_q^2,
P(r)=PGL_2(F_r) acting on P^1(F_r).                       (2)
```

Their degrees agree if and only if

```text
q^2=r+1.                                                  (3)
```

Among prime powers, `(3)` has exactly two solutions:

```text
(q,r)=(2,3), (3,8).                                       (4)
```

On this two-element atlas, each of the following additional conditions is
equivalent to `(q,r)=(2,3)`:

1. `|A(q)|=|P(r)|`;
2. the two natural permutation groups are isomorphic;
3. the affine action is three-transitive; or
4. the projective punctured-direction algebra already exhausts its parity
   algebra.

Thus `(1)` is the unique full affine/projective handshake.  Equal degree and
two-transitivity alone leave the sharp counterfeit `(3,8)`.

## 2. Prime-power factorization

Equation `(3)` gives

```text
r=(q-1)(q+1).                                             (5)
```

If `q` is even, the two odd factors in `(5)` are coprime.  Since their
product is a prime power, one factor is one.  Hence `q-1=1`, giving
`(q,r)=(2,3)`.

If `q` is odd, then `r` is even and therefore a power of two.  Both `q-1`
and `q+1` are powers of two.  Write

```text
q-1=2^u,             q+1=2^v,             1<=u<v.       (6)
```

Subtracting gives

```text
2^u(2^(v-u)-1)=2,                                         (7)
```

so `u=1`, `v=2`, `q=3`, and `r=8`.  This proves `(4)`
without a finite search.

## 3. Group order and transitivity gates

The exact orders are

```text
|A(q)|=q^2(q^2-1)(q^2-q),
|P(r)|=r(r^2-1).                                          (8)
```

Substituting `r=q^2-1` gives

```text
|P(q^2-1)|-|A(q)|=q^2(q^2-1)(q-2).                       (9)
```

Therefore equal order holds exactly at `q=2,r=3`.  The boundary rows are

```text
(2,3): degree 4,       orders 24 and 24;
(3,8): degree 9,       orders 432 and 504.                (10)
```

Both families in `(2)` are two-transitive: translate the first affine point
to zero and use `GL_2(F_q)` on the nonzero difference; on the projective
side use a fractional-linear map.  The projective action is in fact sharply
three-transitive for every `r`, since one fractional-linear map is uniquely
determined by the images of three distinct points.

The affine action is three-transitive only at `q=2`.  For `q>2`, ordered
triples split into collinear and noncollinear orbits.  At the counterfeit
`q=3`, their exact sizes are

```text
72 and 432,                                                (11)
```

whereas `PGL_2(F_8)` has one ordered-triple orbit of size `504`.  At
`q=2,r=3`, both actions have one orbit of size `24`, and both groups are the
full symmetric group `S4` in degree four.  This proves the remaining
equivalences in Section 1 except the scalar clause.

## 4. The septimal counterfeit and the missing coordinate

At `(q,r)=(3,8)`, equal degree and two-transitivity can therefore imitate the
outer combinatorics of the binary/ternary clutch.  They do not identify the
groups or the direction atoms.

Indeed, by [THM-3088](THM-3088-punctured-projective-direction-algebra-and-exceptional-parity-saturation.md),
one punctured `F_8` direction is the scalar orbit

```text
F_8^* ~= C7.                                               (12)
```

Characteristic two makes parity trivial, so projectivization erases this
entire seven-point fibre.  On all nine directions,

```text
parity-ideal rank = 8^2-1=63,
direction rank    = 8+1=9,
internal contrast rank = 54.                              (13)
```

This is the sharp first failed implication:

```text
equal degree + two-transitivity
  does not imply equal group, equal triple geometry,
  or a saturated parity/direction clutch.                  (14)
```

The surviving `C7` is an internal scalar coordinate, not automatically the
septimal carry, owner, or root action in LRC.  A physical identification
would have to transport that coordinate rather than merely count nine outer
states.

## 5. Relation to the modular two/three frame

[THM-2768](THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks.md)
owns the quotient maps from `C2*C3` to `A4,S4` and their Bass--Serre cycle
ranks.  Equations `(3)--(13)` prove no new quotient of the modular group.
They instead classify when two natural finite permutation realizations can
even have the same outer set and show which extra predicates single out the
four-sheet `S4` realization.

The transfer contract is therefore:

| source | target/map | preserved | destroyed | needed sidecar |
|---|---|---|---|---|
| `F_q^2` | `P^1(F_r)`, degree match | number of outer states | group order and affine lines | order or triple geometry |
| `AGL_2(F_q)` | `PGL_2(F_r)` | two-transitivity | collinearity at `(3,8)` | three-point frame |
| punctured `F_r` lines | projective directions | direction label | scalar fibre `F_r^*` | retained within-line coordinate |
| exceptional `(2,3)` | quartic `S4` frame | permutation module after a chosen gauge | owner and phase | physical cover/order data |

Nothing here chooses modular generators, identifies binary and ternary trees,
or supplies a quartic root owner, Keller graph order, LRC carrier, or current.

## 6. Exact evidence

Run

```text
python 04-computation/affine_projective_prime_power_handshake_thm3090.py
python -O 04-computation/affine_projective_prime_power_handshake_thm3090.py
```

Both executions byte-match the stored transcript.  The companion uses only
explicit `require` gates.  It checks the identities `(3),(8),(9)`, enumerates
prime powers through `4096` as a finite hostile control, constructs exact
`F_2,F_3,F_8` arithmetic, enumerates all `24,432,24,504` elements in the four
natural actions, and computes every ordered-pair and ordered-triple orbit.

```text
PROVED IN THE CANDIDATE:
  the complete prime-power degree atlas (2,3),(3,8);
  exact group-order separation and exceptional S4 equality;
  two-/three-transitivity hierarchy and triple-orbit sizes;
  the C7 scalar fibre and rank-54 parity/direction defect at (3,8).

NOT PROVED:
  independent hostile audit or promotion;
  a canonical generator or cross-field gauge;
  a common PSL2(Z), Farey, partial-cube, or tree carrier;
  a quartic-owner, Keller, GMC, or LRC consequence.                       (15)
```

QED (candidate).
