---
id: THM-3632
title: "Russell-cylinder formal-pair algebraization triple-fibre obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  THM-3630's two side germs are restrictions of one rational fold on x!=0;
  their opposite normalized primitives prove that its nonzero formal F is
  not target-local algebraic. Separately, for every fixed retained-collision
  polynomial Q, the pullback ring meets C[w]+Cx exactly in C[w], so its three
  analogous nonzero primitives admit no target-local or global algebraic F.
  Arbitrary mixed target pairs and Keller claims remain open.
source: root / audit_thm3627 THM-3630 globalization boundary, 2026-08-21
audit: >
  PASS -- an independent hostile audit rederived the common rational-side
  continuation and the separate fixed-polynomial triple-fibre theorem,
  checked the local-to-completion and constant-field steps, and verified
  byte-identical normal, optimized, and stored 96-gate transcripts.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3630-russell-cylinder-noneven-formal-survivor-no-finite-jet-bound
related:
  - THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary
  - THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary
script: 04-computation/jc2_russell_cylinder_formal_pair_algebraization_triple_fibre_obstruction_thm3632.py
output: 05-knowledge/results/jc2_russell_cylinder_formal_pair_algebraization_triple_fibre_obstruction_thm3632.out
script_sha256: fee4dd9e9fc7b2b70b8e49e010a9baafb29b0525968a5c79986d84041738a08b
output_sha256: a54e5be841e342d42267d3d87079ee521bf23bca85097fcdb667f6c839225675
hash_basis: raw LF bytes
---

# THM-3632 -- formal-pair algebraization triple-fibre obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem has two distinct algebraization statements.  First, it applies
directly to THM-3630's abstract formal Newton interpolant: its two side
`Q`-germs are restrictions of the same rational function on `x!=0`, and their
opposite normalizations prevent the nonzero formal `F` from being a
target-local algebraic germ.  Second, for every one fixed retained-collision
polynomial `Q`, its analogous three normalized branch primitives admit no
one target-local or global algebraic `F` with second coordinate `w`.

Both obstructions occur before any late Newton denominator.  Algebraic
continuation forces an affine separator in `x`, while regularity over the
retained fibre forces equal values.  The formal completed-ring construction
itself remains valid.

All rings, local rings, completions, and differentials are over `C`.

## 0. Setup and strict scope

Use the THM-3561 exponent-two Danielewski compiler

```text
D=1+x^2q,
b=(D-1)(D+2)^2,
c=xD(D+2),
e=q(D+3),                    c^2e=b(b+4).                (1)
```

Let `Q in C[x]` satisfy only the retained-collision values

```text
Q(-1)=Q(1)=-3,                  Q(0)=-3/4,               (2)
```

and put

```text
q=Q(x)+w^2.                                                (3)
```

Write

```text
M_Q:A2_(x,w) -> Y_2 x A1_w                               (4)
```

for the resulting fixed polynomial fold, and denote its pullback target
ring by

```text
B_Q=C[b_Q,c_Q,e_Q,w] subset C[x,w].                       (5)
```

At

```text
p_-=(-1,0),             p_0=(0,0),             p_+=(1,0),
```

the three compiler values are exactly

```text
M_Q(p_-)=M_Q(p_0)=M_Q(p_+)=(b,c,e,w)=(0,0,-3,0)=y_0.    (6)
```

The principal non-even `rho=0` stratum of THM-3624 additionally has

```text
Q'(-1)=-9/2,       Q'(0)=u,       Q'(1)=9/2,
u notin {0,9/4,-9/4}.                                    (7)
```

The obstruction below needs only `(2)`, so in particular applies to every
fixed polynomial on `(7)`.  It concerns pairs whose second target coordinate
is exactly the stable coordinate `w`.  It says nothing about a pair `(F,G)`
in which both functions genuinely mix the Danielewski and stable variables.

## 1. Direct rational-side obstruction for the THM-3630 packet

The two side germs in THM-3630 are not unrelated formal series.  Both are
restrictions of the one rational function

```text
Q_infinity(x)=-3/4-9/(4x^2)                             (A1)
```

on the connected open set `x!=0`.  Put

```text
g=(1-4w^2/3)^(-1/2).                                    (A2)
```

The rational fold `q=Q_infinity(x)+w^2` is regular near `x=-1` and `x=1`.
Its two moving `c=0` sections are `x=-g` and `x=g`.  The normalized side
primitives for a common scalar `kappa` pull back exactly to

```text
A_-=kappa(x+g),                    A_+=kappa(x-g).       (A3)
```

Suppose that the formal Newton interpolant `F` of THM-3630 were the
completion of an algebraic target-local germ at `y_0`.  Pull it back along
this one rational fold.  The result is a single rational function

```text
r in C(x,w)                                               (A4)
```

regular near both side preimages.  Its completed germ at `p_-` equals
`A_-`, so `r_x=kappa` in that completion.  Injectivity of the local ring into
its completion makes this an identity in `C(x,w)`.  As the constants of
`partial_x` are `C(w)`,

```text
r=kappa x+h(w),                  h in C(w).              (A5)
```

Regularity at the two side points makes `h` regular at `w=0`.  Both
normalized primitives vanish at their base points.  Equations `(A3)--(A5)`
therefore give

```text
x=-1: h(0)= kappa,                  x=1: h(0)=-kappa.    (A6)
```

Thus `kappa=0`.  In particular the `kappa=12` formal `F` constructed in
THM-3630 is **not** the completion of a target-local algebraic germ.  This
argument uses only the common rational continuation of the two side germs;
it does not pretend that the independent central germ and the side germs
come from one fixed polynomial `Q`.

## 2. Fixed-polynomial global pullback-ring intersection

The affine-in-`x` slice of the pullback ring is

```text
B_Q intersect (C[w]+C x)=C[w].                           (8)
```

Indeed, every element of `B_Q` has equal values at the three points in
`(6)`.  If

```text
p(x,w)=kappa x+h(w),                                     (9)
```

then its values there are

```text
h(0)-kappa,                 h(0),                 h(0)+kappa.
                                                                    (10)
```

Equality of `(10)` forces `kappa=0`.  Conversely every `h(w)` belongs to
`B_Q` because `w` is already a target coordinate.  This proves `(8)` with
both inclusions.

Now let `F` be a global regular target function and put

```text
p=M_Q^*F in C[x,w].                                      (11)
```

If, in the completed source local ring at even one of the points `p_i`,

```text
d p wedge d w=kappa d x wedge d w,                       (12)
```

then `p_x=kappa` in that completion.  The polynomial `p_x-kappa` has zero
Taylor series at a closed point and is therefore the zero polynomial.  Thus

```text
p=kappa x+h(w).                                          (13)
```

Equations `(8)` and `(13)` force `kappa=0`.  In particular no global regular
target function can realize a nonzero constant primitive on all three
branches with `w` as the second coordinate.

The same conclusion holds after the Russell-cylinder isomorphism of
THM-3605/THM-3630: there the stable coordinate is the global polynomial

```text
w=Y(B+2)/8-CS.                                           (14)
```

## 3. Fixed-polynomial target-local germ obstruction

The failure is not merely global.  Let

```text
F in O_(Y_2 x A1,y_0)                                    (15)
```

be an algebraic target-local germ.  Its pullback is one rational function

```text
p in C(x,w)                                               (16)
```

regular at all three points `p_i`: a denominator representing `(15)` is a
unit at `y_0`, hence remains a unit at every point over `y_0`.

Suppose `(12)` holds in one completed source branch.  The source local ring
injects into its completion, so `p_x=kappa` already in `C(x,w)`.  Viewing
`C(x,w)` as `C(w)(x)`, the constants of `partial_x` are exactly `C(w)` in
characteristic zero.  Hence again

```text
p=kappa x+h(w),                 h in C(w).               (17)
```

Regularity of `p` at `p_0` makes `h` regular at `w=0`.  Since `(15)` has one
residue at `y_0`, its three pullback residues are equal.  Evaluating `(17)`
at `(6)` gives `(10)` and therefore `kappa=0`.

Equivalently, any rational target expression that pulled back to a nonzero
affine separator `kappa x+h(w)` would need its denominator to vanish at
`y_0` (or fail to represent the separator there).  It cannot be an element
of the target local ring.  This is the first genuine algebraization
obstruction: an order-zero fibre/conductor condition, not a late pole in
`D_21` or `D_321` from THM-3630.

## 4. Fixed-polynomial normalized branch primitives

On the `rho=0` stratum `(7)`, one has `c_x=3` at all three collision points.
Thus each fold branch has a unique formal `c=0` section

```text
x=chi_i(w),                    chi_i(0)=x_i.             (18)
```

For a desired common scalar `kappa`, write the branch area form in `(c,w)`
coordinates and take the normalized primitive

```text
J_i=kappa/c_(x,i),
A_i(c,w)=integral_0^c J_i(r,w)dr.                        (19)
```

Along the source branch, the change of variables `dc=c_(x,i)dx` gives the
exact identity

```text
M_Q^*A_i=kappa(x-chi_i(w)).                              (20)
```

For reference, if

```text
r_-=Q''(-1),                         r_+=Q''(1),          (21)
```

the two side sections begin

```text
chi_-(w)=-1-(2/3)w^2-[4(r_-+18)/27]w^4+O(w^6),
chi_0(w)=0,
chi_+(w)= 1+(2/3)w^2+[4(r_++18)/27]w^4+O(w^6).          (22)
```

For the THM-3630 comparison germs, `r_-=r_+=-27/2`, and `(22)` is the
beginning of `(-g,0,g)` with `g=(1-4w^2/3)^(-1/2)`.

If one algebraic target germ restricted to all three normalized primitives
for this same fixed polynomial `Q`, its single continuation `(17)` and `(20)`
would require

```text
h(w)=-kappa chi_i(w)                    for every i.     (23)
```

Already at `w=0`, the three required values are

```text
kappa,                         0,                       -kappa.
                                                                    (24)
```

Thus `(23)` is compatible if and only if `kappa=0`.  Allowing additive
branch constants does not change the nonzero boundary: any actual target
function has the same constant value at the common target `y_0`, after which
the three affine continuations still differ by `kappa(x_i-x_j)`.

The formal construction of THM-3630 is not contradicted.  An element of the
completed target ring need not pull back to one rational function on the
connected rational side fold.  Section 1 proves precisely that its nonzero
formal `F` is not algebraic target-local.  Sections 2--4 give the separate
statement obtained after replacing all three germs by one fixed
retained-collision polynomial `Q`.  Equality through a prescribed finite jet
does not turn `p_x-kappa` into a rational or polynomial identity.

## 5. Sharp zero boundary and exclusions

For `kappa=0`, equation `(8)` gives the exact positive boundary: every
stable-only `h(w)` is a global target function, and the normalized primitives
in `(19)` are all zero (so `F=0` is a representative).  Functions vanishing
on the fixed fold image may of course be added without changing its
restriction.  This boundary has zero source Jacobian and produces no Keller
map.

Therefore the two proved conclusions are precisely

```text
(a) THM-3630's specific nonzero formal Newton F is not target-local
    algebraic, by continuation across its common rational side fold.

(b) For every one fixed retained-collision polynomial Q, its analogous
    three normalized nonzero branch primitives with second coordinate w
    cannot be realized by one target-local or global algebraic F.       (25)
```

This theorem does **not** prove any of the following:

- that the formal Newton series of THM-3630 is divergent;
- that every formal multi-germ interpolant with different side data is
  nonalgebraic;
- that arbitrary decomposable pairs `(F,G)` with both coordinates mixed fail;
- that every fixed non-even polynomial fails the arbitrary-two-form system;
- that one global Keller pair exists or does not exist; or
- `JC(2)` or any general polynomial all-order no-go.

## 6. Exact reproduction

The deterministic companion verifies

- the compiler relation and exact three-point fibre;
- the common rational THM-3630 side fold and its opposite primitive constants;
- a genuine non-even `rho=0` polynomial control;
- the affine-slice equalizer and its `kappa=0` ideal;
- the polynomial derivative kernel and target-local unit-denominator gate;
- the normalized side-section expansions and primitive constants;
- the THM-3630 `(-g,0,g)` comparison; and
- an AST gate excluding inactive Python `assert` statements.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_formal_pair_algebraization_triple_fibre_obstruction_thm3632.py
python3 -O 04-computation/jc2_russell_cylinder_formal_pair_algebraization_triple_fibre_obstruction_thm3632.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_formal_pair_algebraization_triple_fibre_obstruction_thm3632.out`.
