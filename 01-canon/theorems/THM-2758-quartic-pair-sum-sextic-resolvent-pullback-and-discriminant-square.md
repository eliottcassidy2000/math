---
id: THM-2758
title: "Quartic pair-sum sextic, resolvent pullback, and discriminant square"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a monic
  quartic with elementary root functions e1,e2,e3, the
  sextic whose roots are its six pair sums becomes, after translation by
  e1/2, the quadratic pullback of the standard cubic resolvent translated by
  e1^2/4-e2.  With lexicographic edge ordering its Vandermonde is exactly the
  quartic discriminant times T=e1^3-4e1e2+8e3; hence its discriminant is
  disc(f)^2*T^2.  The cubic-resolvent discriminant remains disc(f).  Thus a
  separable quartic has separable pair-sum sextic iff no opposite pair sums
  agree.  The separable hostile x^4-10x^2+9 has T=0 and pair sums
  -4,-2,0,0,2,4.  This is a universal quartic identity, not a Keller
  exclusion.
source: a4-resolvent-next-gate-scout-2026-07-28
audit: >
  root/2026-07-28 (independent pullback, resolvent convention,
  Vandermonde sign, characteristic scope, separability boundary, hostile,
  and normal/-O/LF-hash replay: ACCEPT)
depends_on: []
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2743-c2-c3-off-diagonal-projector-rank-and-s3-s4-compatibility-defect
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch
script: 04-computation/quartic_pair_sum_resolvent_discriminant_thm2758.py
output: 05-knowledge/results/quartic_pair_sum_resolvent_discriminant_thm2758.out
script_sha256: 5675989aeb6b02475e2699fdf63ad6315e243ea0123caf2e99e0ca83899bbc16
output_sha256: a7ecc7ed38b40a3ecba3bae0fe108c72b834ebe7370a78807c06634cdffd7b86
hash_basis: LF-normalized bytes
---

# THM-2758 -- the quartic pair-sum sextic is a resolvent pullback

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2753 and THM-2756 describe the six edges of a four-sheet set at the level
of finite representations.  Edge complementation produces three opposite
pairs, the matching quotient is the cubic-resolvent `S3` set, and the ambient
six-edge sign is trivial.  The present theorem is the exact polynomial avatar
of that structure.

The six pair sums of four roots come in three pairs symmetric about half the
total root sum.  Squaring the centered pair sums produces a translate of the
three classical resolvent roots.  At the discriminant level the twelve
adjacent-edge differences give the quartic discriminant, while the three
opposite-edge differences give one extra cubic factor.  Nothing asymptotic or
Galois-generic is used.

## 1. Conventions

Let `K` be a field of characteristic different from `2`, and let

```text
f(X)=product_(i=1)^4 (X-r_i)
    =X^4-e1 X^3+e2 X^2-e3 X+e4.                          (1)
```

The identities are symmetric and therefore descend from a splitting field to
`K`.  Order the six edges lexicographically as

```text
12,13,14,23,24,34,                                         (2)
```

put `s_ij=r_i+r_j`, and define the pair-sum sextic

```text
G(U)=product_(i<j)(U-s_ij).                                (3)
```

The three oriented opposite-edge differences are

```text
d_1=s_12-s_34,       d_2=s_13-s_24,       d_3=s_14-s_23, (4)
```

and

```text
T=d_1 d_2 d_3=e1^3-4e1 e2+8e3.                            (5)
```

The sign in `(5)` is fixed by `(2)--(4)`.  For the coefficient convention

```text
f(X)=X^4+aX^3+bX^2+cX+d,
```

one has `T=-(a^3-4ab+8c)`, so the wall `T=0` is

```text
a^3-4ab+8c=0.                                              (6)
```

The cubic factor is also the first odd centered root invariant.  Direct
Newton expansion gives the integral identity

```text
sum_i (4r_i-e1)^3=24T.                                    (6a)
```

Over characteristic zero, equivalently
`sum_i(r_i-e1/4)^3=(3/8)T`.  Thus the opposite-sum collision wall is exactly
the zero centered-cubic-skewness wall; for a depressed quartic it is simply
the vanishing of its linear coefficient.

## 2. Exact centered pullback of the cubic resolvent

Use the standard three resolvent roots

```text
z_1=r_1r_2+r_3r_4,
z_2=r_1r_3+r_2r_4,
z_3=r_1r_4+r_2r_3,                                        (7)
```

and let

```text
C(Z)=product_(m=1)^3 (Z-z_m)
    =Z^3-e2 Z^2+(e1e3-4e4)Z
       -(e1^2e4-4e2e4+e3^2).                              (8)
```

For each matching,

```text
d_m^2/4=z_m+e1^2/4-e2.                                   (9)
```

Indeed, for example,

```text
(r_1+r_2-r_3-r_4)^2
 =e1^2-4e2+4(r_1r_2+r_3r_4),                              (10)
```

and the other cases follow by relabelling.  Since the two pair sums in one
matching are `e1/2 +/- d_m/2`, equations `(3)` and `(9)` give the polynomial
identity

```text
G(e1/2+V)
 =product_m (V^2-d_m^2/4)
 =C(V^2+e2-e1^2/4).                                      (11)
```

Thus `G` is not merely associated with the cubic resolvent: after its canonical
translation, it is the degree-two pullback of a translated resolvent cubic.
Equivalently, if

```text
R(Y)=C(Y-(e1^2/4-e2)),                                    (12)
```

then `G(e1/2+V)=R(V^2)` and the roots of `R` are `d_m^2/4`.

## 3. The two discriminant identities

The resolvent differences factor as

```text
z_1-z_2=(r_1-r_4)(r_2-r_3),
z_1-z_3=(r_1-r_3)(r_2-r_4),
z_2-z_3=(r_1-r_2)(r_3-r_4).                              (13)
```

Squaring their product gives every quartic root difference exactly once:

```text
disc(C)=disc(f).                                          (14)
```

For the sextic, let

```text
V_E=product_(a<b)(s_a-s_b),                               (15)
```

where `a,b` refer to the lexicographic edge order `(2)`.  Of its fifteen
factors, twelve come from adjacent edge pairs and three from opposite pairs.
Every sheet difference `r_i-r_j` occurs twice in the adjacent product, with
the resulting sign positive in convention `(2)`.  The opposite product is
exactly `(5)`.  Hence

```text
V_E=disc(f) T.                                            (16)
```

Squaring `(16)` proves

```text
disc(G)=disc(f)^2 (e1^3-4e1e2+8e3)^2.                    (17)
```

The pullback form gives an independent derivation.  If `t_m=d_m^2/4`, then

```text
disc(product_m(V^2-t_m))
 =4^3 product_m(t_m) disc(R)^2.                           (18)
```

Here `product_m(t_m)=T^2/4^3`, translation preserves discriminant, and
`disc(R)=disc(C)=disc(f)`, again yielding `(17)`.

## 4. Sharp separability boundary and hostiles

Equation `(17)` gives the exact iff:

```text
G is separable  iff  f is separable and T != 0.           (19)
```

The two failure mechanisms are geometrically different.

- If `f` has a repeated root, two adjacent pair sums collide; this is the
  `disc(f)` factor.
- If `f` is separable but `T=0`, one opposite pair has equal sums; this is the
  extra cubic factor invisible to `disc(f)`.

The second boundary is load-bearing.  Take

```text
f(X)=X^4-10X^2+9=(X+3)(X+1)(X-1)(X-3).                  (20)
```

Then `disc(f)=589824` is nonzero but `T=0`, and the six pair sums are

```text
-4,-2,0,0,2,4.                                           (21)
```

Thus

```text
G(U)=U^2(U^2-4)(U^2-16)
```

has a double root.  This refutes any claim that quartic separability alone
makes the pair-sum cover etale.  Conversely, a repeated quartic root such as
`0,0,1,2` gives adjacent collisions and tests the other factor in `(17)`.

## 5. Representation meaning and Keller boundary

The degree-six pair-sum cover has monodromy inside the even image
`S4 -> A6` of THM-2753; `(17)` supplies its explicit discriminant square.  The
three quadratic fibres in `(11)` are the opposite-edge pairs, and their base
is exactly the three-matching resolvent set.  The extra factor `T` records
whether one quadratic fibre ramifies over its center.  This is the polynomial
form of THM-2756's ternary matching block plus binary opposition involution.

For a hypothetical degree-four Keller fibre, these identities apply to its
quartic polynomial exactly as they apply to every quartic.  They do **not**
exclude such a fibre: one must first identify a separately constrained graph
sextic or boundary cover with `(3)`, then transport its ramification data.
Neither a graph-quartic realization nor that transport is proved here.  In
particular, the fact that `disc(G)` is a square is universal and cannot by
itself contradict the Keller condition.

## 6. Exact verification

Run

```bash
python 04-computation/quartic_pair_sum_resolvent_discriminant_thm2758.py
python -O 04-computation/quartic_pair_sum_resolvent_discriminant_thm2758.py
```

Both executions byte-match the stored `17`-line transcript
`05-knowledge/results/quartic_pair_sum_resolvent_discriminant_thm2758.out`.
The companion uses explicit exceptions and no truth-bearing Python assertions.
It works in a custom exact polynomial ring to prove `(5)`, `(6a)`, all three
instances of `(9)`, and the factorizations `(13)`; it classifies all fifteen edge pairs,
including the exact Vandermonde sign in `(16)`.  It then checks `(9)`, `(14)`,
and `(16)--(17)` on all `210` distinct quadruples from `[-4,5]`, including
`50` opposite-sum walls, and on `65` repeated-root controls from `[-2,2]`.
The generic row `(0,1,3,7)` and hostile `(20)` are printed explicitly.

## 7. Boundary ledger

```text
PROVED HERE:              exact centered pair-sum/resolvent pullback;
                          explicit integral standard resolvent convention;
                          disc(resolvent)=disc(quartic);
                          signed edge Vandermonde=disc(quartic)*T;
                          disc(pair-sum sextic)=disc(quartic)^2*T^2;
                          T as the normalized centered cubic root invariant;
                          separability iff and two collision mechanisms;
                          separable-quartic/multiple-sextic hostile.

NOT PROVED:               a new restriction on arbitrary quartics;
                          graph-quartic or Keller-fibre realization;
                          transport to a Jelonek/boundary sextic;
                          exclusion of A4 or S4 Keller monodromy;
                          a physical LRC six-packet polynomial;
                          JC(2), DC(2), or LRC(14).                       (22)
```

QED.
