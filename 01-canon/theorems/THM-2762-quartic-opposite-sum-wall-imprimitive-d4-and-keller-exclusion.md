---
id: THM-2762
title: "Quartic opposite-sum wall, imprimitive D4 descent, and the alternating Keller sextic"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For a separable quartic in characteristic different from two,
  T=e1^3-4e1e2+8e3 vanishes iff its root set is centrally symmetric about
  e1/4.  The resulting antipodal pairing is Galois-invariant, so an
  irreducible T=0 quartic has transitive monodromy C4, V4, or D4.  THM-2633
  excludes all three for a polynomial Keller map.  Hence every primitive
  quartic presentation of a generic degree-four Keller extension has
  T nonzero.  Its six pair sums are therefore distinct; the live A4/S4
  monodromy is transitive on them, and their irreducible sextic has square
  discriminant.  This is a field-level auxiliary resolvent, not an affine
  cover, inverse, JC(2), or DC(2) proof.
source: root/quartic-opposite-sum-wall-2026-07-28
depends_on:
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2758-quartic-pair-sum-sextic-resolvent-pullback-and-discriminant-square
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
script: 04-computation/quartic_opposite_sum_imprimitive_keller_thm2762.py
output: 05-knowledge/results/quartic_opposite_sum_imprimitive_keller_thm2762.out
script_sha256: f529cf4586e0373cf03381c9d209f892e3876877ebb702650da04cf43c103f07
output_sha256: baef2224994002c00e6491392b620431ceb0777768fdb69d451b5d9957a67dfa
hash_basis: LF-normalized bytes
---

# THM-2762 -- the opposite-sum wall is the imprimitive quartic wall

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2758 identifies one additional factor in the discriminant of the six
pair sums of a quartic.  That factor is not another mysterious resolvent
invariant.  It is precisely the obstruction to central symmetry of the four
roots.  This makes its Keller consequence much cleaner than a discriminant
comparison: the wall forces an imprimitive block system, and the existing
monodromy theorem already excludes every transitive group on that wall.

## 1. The wall is central symmetry

Let `K` be a field with `char(K)!=2`, and let

```text
f(X)=product_(i=1)^4 (X-r_i)
    =X^4-e1 X^3+e2 X^2-e3 X+e4                         (1)
```

be separable.  Put

```text
c=e1/4,                    y_i=r_i-c,
T=e1^3-4e1e2+8e3.                                         (2)
```

The centered polynomial has the form

```text
f(Y+c)=Y^4+pY^2+qY+s,             q=-T/8.                (3)
```

Indeed its cubic coefficient is zero, and the elementary translation gives
`q=-(e1^3-4e1e2+8e3)/8`.  Equivalently, Newton's identity gives

```text
sum_i y_i^3=3e3(y_1,...,y_4)=3T/8.                       (4)
```

Consequently the following are equivalent:

```text
T=0;
f(Y+c) is even;
{r_1,r_2,r_3,r_4} is stable under r |-> e1/2-r;
one opposite pair of pair sums is equal.                 (5)
```

The last equivalence can also be read directly from THM-2758:

```text
T=(r1+r2-r3-r4)(r1+r3-r2-r4)(r1+r4-r2-r3).              (6)
```

Separability makes the involution in `(5)` fixed-point-free.  If a centered
root were zero, evenness would make zero a root of multiplicity at least two.
Thus the four roots split into two honest antipodal pairs.

## 2. Central symmetry forces imprimitive monodromy

Let `G` be the Galois group of the splitting field of `f` over `K`, acting on
the roots.  The involution

```text
iota(r)=e1/2-r                                             (7)
```

is defined over `K`; hence every `g in G` satisfies

```text
g iota = iota g.                                          (8)
```

Therefore `G` preserves the two-block partition into the orbits of `iota`.
The stabilizer in `S4` of such a perfect matching is

```text
C2 wreath S2 = D4,                 order 8.               (9)
```

If `f` is irreducible, then `G` is transitive.  The exact transitive
subgroups of this matching stabilizer are

```text
C4,                         V4,                         D4. (10)
```

This is sharp algebraically.  The irreducible polynomial

```text
X^4-2                                                     (11)
```

lies on `T=0`; its splitting field is `Q(2^(1/4),i)` and its Galois group is
`D4`.  Thus no group-theoretic contradiction exists without the Keller
input.

## 3. Degree-four Keller consequence

Let

```text
F:A_C^n -> A_C^n                                         (12)
```

be a polynomial Keller map of finite generic degree four.  Write

```text
K=C(F_1,...,F_n) subset L=C(x_1,...,x_n)                 (13)
```

and choose any primitive element `theta` for the separable quartic extension
`L/K`.  Its minimal polynomial is an irreducible separable quartic over `K`.
If its invariant `T_theta` vanished, Sections 1--2 would put its transitive
geometric monodromy in `C4`, `V4`, or `D4`.  THM-2633 excludes all three as
Keller monodromy groups.  Hence

```text
T_theta != 0                                              (14)
```

for every primitive quartic presentation of the extension.

THM-2633 leaves only `A4` and `S4`.  Both act transitively on the six
unordered pairs of four sheets.  By `(14)` and THM-2758, the six pair sums of
the conjugates of `theta` are distinct.  Their polynomial

```text
G_theta(U)=product_(i<j)(U-theta_i-theta_j)               (15)
```

is therefore irreducible and separable of degree six.  Moreover

```text
disc(G_theta)=disc(f_theta)^2 T_theta^2                   (16)
```

is a square in `K*`.  Thus every degree-four Keller survivor has a canonical
primitive alternating degree-six **field-level pair-sum sidecar**.

This conclusion does not construct a finite polynomial cover of affine
space: the pair sums are elements of the Galois closure, and their integral
or boundary behavior is not controlled here.  It also does not replace the
`V4` Kummer unit/class-group gate of THM-2655.

## 4. Exact verification and boundary

Run

```bash
python 04-computation/quartic_opposite_sum_imprimitive_keller_thm2762.py
python -O 04-computation/quartic_opposite_sum_imprimitive_keller_thm2762.py
```

The exact companion uses explicit exceptions and no truth-bearing Python
assertions.  It verifies the coefficient translation on a `7^4` coefficient
grid, checks `(5)--(6)` on all `210` distinct integer root quadruples in
`[-4,5]`, enumerates the ten subgroups of one matching stabilizer, and finds
exactly one transitive `C4`, one transitive `V4`, and the full `D4`.  It also
checks that `A4` and `S4` are transitive on the six edges and that their edge
actions are even.

```text
PROVED HERE (candidate):  T=0 iff centered quartic is even/centrally symmetric;
                          invariant antipodal block system;
                          transitive wall groups C4,V4,D4;
                          exclusion of T=0 for every primitive quartic
                          presentation of a degree-four Keller extension;
                          irreducible separable alternating pair-sum sextic.

NOT PROVED:               an affine realization of the pair-sum sextic;
                          extension across the Jelonek divisor;
                          exclusion of A4 or S4;
                          degree-one conclusion;
                          JC(2), general JC, or DC(2).                    (17)
```

QED (candidate).
