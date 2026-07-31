---
id: THM-2821
title: "Sextic e=3 power-wreath and Chebyshev-dihedral block lattice"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  THM-2817's
  power carrier has monodromy
  (C3 x C3) semidirect C2 of order 18 and exactly one nontrivial block
  system, of block size three.  Its Chebyshev carrier has dihedral
  monodromy of order 12 and exactly two nontrivial block systems, of
  sizes three and two.  These give all rational functional decompositions
  up to intermediate-coordinate Mobius equivalence.  This is a scalar
  cover theorem, not a Keller-map decomposition or a result on JC(2).
source: root/sextic-e3-monodromy-block-lattice-2026-07-28
depends_on:
  - THM-2817-sextic-e3-maximal-pole-power-chebyshev-accessory-classification
related:
  - THM-2643-degree-five-six-keller-stabilizer-and-regular-block-quotient-census
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
script: 04-computation/jc_sextic_e3_monodromy_block_lattice_thm2821.py
output: 05-knowledge/results/jc_sextic_e3_monodromy_block_lattice_thm2821.out
script_sha256: 571d6293559be66faf36b364d64b13f962de55018bef41bf100baf560c8a56ec
output_sha256: 4d2456d8a9396536a8d34e7bcd16a1e3426b031b325d5749643a4682cd37c9c4
hash_basis: LF-normalized bytes
---

# THM-2821 -- the two sextic carriers have different 2/3 block lattices

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
phi(z)=z^2/(z^2-1),

F_pow(x)=phi(2x^3-1),
F_cheb(y)=phi(T_3(y)),             T_3(y)=4y^3-3y.       (1)
```

These are the two unmarked response carriers proved complete in THM-2817.
Their common outer formula hides different intermediate-field lattices.
The power carrier has only the inner-degree-three/outer-degree-two
decomposition.  The Chebyshev carrier also has the reverse
inner-degree-two/outer-degree-three decomposition.  The distinction is
exactly the distinction between the subgroup intervals

```text
C3 < C3 x C3 < (C3 x C3) semidirect C2

and

C2 < C2 x C2, Dih(C3) < Dih(C6).                    (2)
```

Here `Dih(C6)` has order `12`; the two middle subgroups in its interval
are incomparable.

## 1. Branch fibres and explicit branch permutations

All permutations below act on `{1,...,6}`, products are composed
right-to-left, and the branch values are ordered `(0,infinity,1)`.

For the power carrier put `P=2x^3-1`.  Exact factorization gives

```text
P^2-1=4x^3(x^3-1),             F_pow-1=1/(P^2-1).    (3)
```

Thus the fibre over zero consists of three points of index two; the fibre
over infinity consists of one point of index three and three simple
points; and the fibre over one is the single point at source infinity of
index six.  Compatible sheet labels give

```text
sigma_0        =(1 4)(2 5)(3 6),
sigma_infinity =(1 2 3),
sigma_1        =(1 4 3 6 2 5),                       (4)

sigma_0 sigma_infinity sigma_1=1.
```

For the Chebyshev carrier,

```text
T_3(y)^2-1
 =(y-1)(y+1)(2y-1)^2(2y+1)^2.                       (5)
```

Its fibre partitions over `(0,infinity,1)` are respectively
`(2,2,2)`, `(2,2,1,1)`, and `(6)`.  One compatible branch tuple is

```text
tau_0        =(1 6)(2 5)(3 4),
tau_infinity =(2 6)(3 5),
tau_1        =(1 2 3 4 5 6),                         (6)

tau_0 tau_infinity tau_1=1.
```

This convention agrees exactly with THM-2817's Nielsen convention.  There
one writes `rho=sigma_1^-1` and `tau=sigma_0`, so the pole permutation is

```text
tau rho=sigma_0 sigma_1^-1=sigma_infinity,            (6a)
```

where the last equality is the branch-product relation.  Thus the
`tau rho` used in that theorem is not a left/right product mismatch.

There are no omitted branch values.  Direct differentiation shows that
the only finite critical points away from numerator zeros are `x=0` in
the power case and `y=+-1/2` in the Chebyshev case; all of these already
occur in `(3)` or `(5)`.  Equivalently, the defects in each displayed
passport sum to `10=2(6)-2`.

The tuple `(4)` can also be obtained without drawing a dessin.  Apply the
target transformation

```text
M(t)=1/(1-t).
```

Then

```text
M(F_pow)=1-(2x^3-1)^2=4x^3(1-x^3).                  (7)
```

A loop around the quadratic outer branch interchanges the two triples of
cube roots, giving `sigma_0`; a loop around the image of `x^3=0` rotates
only one triple, giving `sigma_infinity`.  This lifting constructs `(4)`
for the actual cover rather than merely matching its cycle types.

For `(6)`, use

```text
M(F_cheb)=1-T_3(y)^2=(1-T_6(y))/2.                  (8)
```

The inverse branches `y=cos((+-arccos z+2 pi k)/6)` show directly that
loops about the two finite Chebyshev branch values are the two reflection
classes of a regular hexagon.  They are the fixed-point-free reflection
`tau_0` and the two-fixed-point reflection `tau_infinity` in `(6)`.

## 2. The power group and its complete block lattice

Set

```text
a=sigma_infinity=(1 2 3),       s=sigma_0,
b=sas=(4 5 6).                                      (9)
```

Then `a` and `b` commute, both have order three, and

```text
N=<a,b> = C3 x C3,                 |N|=9,
sas=b,   sbs=a,                    s^2=1.            (10)
```

Consequently

```text
G_pow=<a,s>=N semidirect <s>
     =(C3 x C3) semidirect C2,     |G_pow|=18,        (11)
```

where the involution exchanges the two `C3` factors.

The stabilizer of sheet `1` is

```text
H_pow=<b>,                         |H_pow|=3.          (12)
```

Let `K` be any subgroup with `H_pow <= K <= G_pow`.  If `K` lies in `N`,
then viewing `N` as the two-dimensional vector space over `F_3` shows
that `K` is either the line `H_pow` or all of `N`.  If `K` does not lie
in `N`, it contains an element `ns` with `n in N`.  Since `N` is
abelian,

```text
(ns)b(ns)^-1=a.
```

Hence `K` contains `N`, and then it contains `s`; so `K=G_pow`.  Therefore
the entire interval is

```text
H_pow < N < G_pow,                 orders 3,9,18.     (13)
```

By the subgroup/block correspondence, the only nontrivial block system is

```text
{{1,2,3},{4,5,6}}.                                  (14)
```

In particular there is no invariant partition into three pairs.

## 3. The Chebyshev group and its complete block lattice

Put

```text
r=tau_1=(1 2 3 4 5 6),          s=tau_infinity.      (15)
```

The branch tuple gives

```text
r^6=s^2=1,                 srs=r^-1.                  (16)
```

The twelve permutations `r^i` and `r^i s`, `0<=i<6`, are distinct, so

```text
G_cheb=<r,s>=Dih(C6),              |G_cheb|=12.       (17)
```

Its sheet-`1` stabilizer is `H_cheb=<s>`.  Any subgroup strictly larger
than `H_cheb` contains a rotation: if its new element is `r^i s`, multiply
it by `s`.  Since every subgroup of `<r>=C6` is unique at its order, every
intermediate subgroup is

```text
<s,r^d>,                   d in {6,3,2,1}.            (18)
```

The complete interval and the orbits of sheet `1` are therefore

```text
subgroup                 order       orbit

<s>                        2         {1}
<s,r^3> = C2 x C2          4         {1,4}
<s,r^2> = Dih(C3)          6         {1,3,5}
<s,r>   = Dih(C6)         12         {1,...,6}.       (19)
```

Thus there are exactly two nontrivial block systems:

```text
size 2: {{1,4},{2,5},{3,6}},
size 3: {{1,3,5},{2,4,6}}.                           (20)
```

This proves completeness, not merely the existence of two visible
factorizations.

## 4. All rational functional decompositions

For a separable rational cover `F in C(x)`, a nontrivial decomposition

```text
F=A after B,                  deg A,deg B >1           (21)
```

gives the intermediate field

```text
C(F) proper-subset C(B) proper-subset C(x).
```

Conversely, the intermediate-field/block correspondence followed by
Luroth's theorem produces `(21)`.  Replacing

```text
B by nu after B,             A by A after nu^-1,
```

for `nu in PGL_2(C)` is exactly the ambiguity in choosing a generator of
the same intermediate field.  Hence the complete block lattices
`(13)` and `(19)` classify all rational decompositions up to Mobius
equivalence.

For the power carrier, the unique nontrivial class is

```text
F_pow=phi after P,        P(x)=2x^3-1,   deg(P,phi)=(3,2).
                                                               (22)
```

Equivalently, using `u=x^3`,

```text
F_pow=A_pow(u),
A_pow(u)=(2u-1)^2/[4u(u-1)].                         (23)
```

The affine change `u -> 2u-1` identifies `(22)` and `(23)`.

For the Chebyshev carrier the size-three system gives

```text
F_cheb=phi after T_3,                    degrees (3,2). (24)
```

The size-two system gives the reverse factorization

```text
F_cheb=Psi after (y^2),                  degrees (2,3),

Psi(u)=u(4u-3)^2/[u(4u-3)^2-1]
      =u(4u-3)^2/[(u-1)(4u-1)^2].                    (25)
```

The structural reason for `(25)` is the commuting Chebyshev identity.
With `T_2(z)=2z^2-1` and `chi(w)=(w+1)/(w-1)`,

```text
T_6=T_2 after T_3=T_3 after T_2,
F_cheb=chi after T_6
      =(chi after T_2) after T_3
      =(chi after T_3) after T_2.                    (26)
```

Since `T_2(y)` and `y^2` differ by an affine change, the second expression
in `(26)` is exactly `(25)` up to the allowed Mobius equivalence.  There
are no further nontrivial decompositions of either sextic map.

## 5. Exact controls

The companion script independently rebuilds the two carriers from
THM-2817's `(D,E,A)` formulas and checks:

1. the fibre factorizations `(3)` and `(5)`;
2. the branch products, cycle types, transitivity, and group orders;
3. the semidirect-product and dihedral relations;
4. every uniform set partition of a six-element fibre;
5. every subgroup of each finite permutation group, then the complete
   point-stabilizer interval;
6. the subgroup-orbit/block correspondence; and
7. all five noncrossing matchings of the hexagon, their two rotation
   orbits, and their two group/block signatures; and
8. `(7)`, `(8)`, and every decomposition `(22)--(26)`.

The group routines use raw permutation tuples and exhaustively generate
subgroups by adjoining every ambient element.  They use no group-theory
CAS and contain no Python `assert` node.  Run

```text
python 04-computation/jc_sextic_e3_monodromy_block_lattice_thm2821.py
python -O 04-computation/jc_sextic_e3_monodromy_block_lattice_thm2821.py
```

The normal, optimized, and stored transcripts agree byte for byte.
An independent hostile audit separately rechecked the branch products and
cycle types, both point-stabilizer intervals, the block orbits, the
Luroth correspondence, and the explicit `Psi` factorization; it found no
omitted block system or inequivalent decomposition.

## 6. Keller and Jacobian boundary

This theorem classifies intermediate fields of the scalar extensions

```text
C(F_pow) subset C(x),          C(F_cheb) subset C(y). (27)
```

It does not construct a second coordinate, preserve a constant Jacobian,
verify Faber flux, enter a Keller chart, or decompose a two-variable
polynomial endomorphism.  Equations `(7)` and `(8)` do give honest
one-variable polynomial decompositions: respectively through `x^3`, and
through both `T_3` and `T_2`.  But those target-polynomialized covers are
still ramified one-variable maps of degree six.  A scalar polynomial
decomposition is not a decomposition of a two-variable polynomial Keller
map, and target Mobius polynomialization does not supply the missing
Keller coordinate.

Accordingly, neither the extra Chebyshev block system nor the missing
power block system proves anything by itself about `JC(2)` or `DC(2)`.
Transport to those problems requires a sidecar retaining the second
coordinate and the Keller constraints.
