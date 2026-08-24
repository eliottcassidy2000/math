---
id: THM-3939
title: "Two-boundary elliptic resolvents have at most one new three-character"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED. Let a normal affine
  surface over A1 have integral closed fibres and generic fibre an elliptic
  curve with exactly two rational points deleted. If the difference of those
  boundary points has infinite order, then the surface has scalar units and
  its class group is the Mordell--Weil group modulo that one difference. Its
  three-torsion rank is at most the rational Mordell--Weil three-torsion rank
  plus one. In particular, without a rational three-torsion section there is
  at most one smooth-locus C3 character, independent of the free Mordell--Weil
  rank. A nonzero Cardano class is then the unique character up to inversion.
source: root / post-THM-3937 global character-rank reframe, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root independent audit, 2026-08-24). The
  audit reconstructed the principal vertical localization, elliptic divisor
  quotient and unit calculation, the one-relation three-torsion exact
  sequence without a finite-generation assumption, and the Hartogs/Kummer
  step. Torsion boundary, rational E(K)[3], extra boundary points, and split
  vertical fibres were checked as genuine out-of-scope hostiles.
related:
  - THM-3912-even-degree-split-boundary-cusp-three-torsion-sieve
  - THM-3922-affine-plane-completion-free-boundary-class-group-obstruction
  - THM-3935-linear-conic-resolvent-class-group-unique-cubic-character
  - THM-3937-linear-conic-fold-three-family-uniform-resolvent-rigidity
  - THM-3940-i7-rank-two-linear-cross-term-resolvent-unique-character
---

# THM-3939 -- free Mordell--Weil rank cannot supply the second character

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.** Work over an algebraically closed
field `k` of characteristic zero and put `K=k(t)`. Let
`B` be a normal affine `k[t]`-domain of dimension two. Assume:

1. every closed fibre `B/(t-lambda)` is integral;
2. the generic fibre is a smooth affine curve

   ```text
   C=E minus {P_0,P_1},                                  (1)
   ```

   where `E/K` is an elliptic curve and `P_0,P_1 in E(K)` are distinct;
3. the boundary difference

   ```text
   D=P_1-P_0 in E(K)                                     (2)
   ```

   has infinite order.

Write `S=Spec(B)` and `U=Sreg`. Then

```text
B^*=k*,
Cl(B)=Pic(C)=E(K)/<D>,                                   (3)
```

and there is a natural exact sequence of `F_3`-vector spaces

```text
0 -> E(K)[3] -> Cl(B)[3] -> I_D -> 0,                   (4)
```

where `I_D` is either zero or a subspace of `F_3`. Consequently

```text
dim_F3 Cl(B)[3] <= dim_F3 E(K)[3]+1.                    (5)
```

If `E(K)[3]=0`, the conclusion sharpens to

```text
Cl(B)[3] = 0                    if D notin 3E(K),
Cl(B)[3] = Z/3                  if D in 3E(K).            (6)
```

Normal Hartogs and Kummer theory then give

```text
H^1_et(U,mu_3)=Cl(B)[3].                                 (7)
```

Thus an exact resolvent in this grammar with no rational Mordell--Weil
three-torsion has at most one nontrivial `C3` cover up to inversion. If its
natural Cardano divisor is nonzero, it spans `(7)` and no alternative normal
`S3` cubic field can be obtained from a second codimension-one-unramified
character. This last statement concerns the exact branch and quadratic
resolvent; descent, normality, and monogenicity of the natural cubic remain
separate gates.

## 1. Integral vertical fibres remove the localization sidecar

Let `V_lambda` be a vertical height-one prime. Hypothesis 1 says there is
exactly one such prime over each closed base point and that its fibre is
reduced. Therefore

```text
V_lambda=(t-lambda),                  div(t-lambda)=V_lambda. (8)
```

The Weil localization sequence is

```text
B^* -> Gamma(C,O_C)^* -> direct_sum_lambda Z[V_lambda]
    -> Cl(B) -> Pic(C) -> 0.                              (9)
```

Every generator in the middle group is the divisor of `t-lambda`, so the
middle map has zero cokernel. Hence

```text
Cl(B)=Pic(C).                                             (10)
```

This hypothesis is load-bearing. A reducible closed fibre can contribute
component classes not visible on the generic elliptic curve.

## 2. Two deleted points impose exactly one Mordell--Weil relation

Choose any rational origin on `E`. The usual degree decomposition is

```text
Pic(E)=Z + E(K).                                         (11)
```

The two deleted rational points represent `(1,P_0)` and `(1,P_1)`. The
divisor localization sequence for `(1)` therefore gives

```text
Pic(C)=Pic(E)/<[P_0],[P_1]>
      =E(K)/<P_1-P_0>
      =E(K)/<D>.                                         (12)
```

There are no nonconstant generic-fibre units. Indeed, a unit on `C` has a
degree-zero divisor supported on `P_0,P_1`, necessarily

```text
n(P_1-P_0)=nD.                                           (13)
```

If it were principal for `n!=0`, then `nD=0`, contradicting Hypothesis 3.
Thus

```text
Gamma(C,O_C)^*=K^*.                                      (14)
```

A member of `K^*` that is a unit of `B` has valuation zero at every prime
`(8)`. Since `k` is algebraically closed, a rational function on `A1_t`
without a finite zero or pole is constant. This proves `B^*=k^*` and all of
`(3)`.

## 3. The one-relation three-torsion exact sequence

The group-theoretic mechanism does not require finite generation. Put

```text
M=E(K),                         Q=M/<D>.                  (15)
```

Because `D` has infinite order, `<D>` meets `M[3]` trivially, so the natural
map

```text
M[3] -> Q[3]                                              (16)
```

is injective. If `m+<D> in Q[3]`, there is a unique integer `n` such that

```text
3m=nD.                                                    (17)
```

Send this class to `n mod 3`. Changing `m` by a multiple of `D` changes `n`
by a multiple of three, so this is well-defined. Its kernel consists exactly
of the image of `M[3]`: if `n=3j`, then

```text
3(m-jD)=0.                                                (18)
```

This proves `(4)` with `I_D` the image in `Z/3`, and hence `(5)`.

Suppose now that `M[3]=0`. If `D=3R`, then `R+<D>` is a nonzero element of
order three in `Q`. Conversely, if `Q[3]` is nonzero, its image in `I_D` is
nonzero. Equation `(17)` has `n=1` or `2` modulo three. In the second case,

```text
3(2m-((2n-1)/3)D)=D;
```

in the first case,

```text
3(m-((n-1)/3)D)=D.
```

Hence `D in 3M`. This proves `(6)`.

Finally, normality identifies `Pic(U)=Cl(B)` and gives
`Gamma(U,O_U)=B`. Since `k^*` is cube-divisible, the etale Kummer sequence
proves `(7)`.

## 4. Counterexample-design consequence

THM-3935 and THM-3937 sit in the sharp second row of `(6)`: their boundary
difference is `3Q`, their Mordell--Weil group has no rational three-torsion,
and their unique nonzero class is the Cardano class. The present theorem
shows that increasing their **free** Mordell--Weil rank would not change that
conclusion. Quotienting by one boundary difference creates at most one new
three-primary invariant factor.

A search for a genuinely different cubic character must therefore change at
least one of the following coordinates:

1. produce a rational `E(K)[3]` section that survives the deck-action and
   descent gates;
2. introduce additional boundary or nonprincipal vertical relations whose
   Smith form has more than one factor divisible by three;
3. move to a higher-genus generic fibre whose Jacobian has a larger rational
   three-torsion packet;
4. leave the two-boundary elliptic localization grammar entirely.

These are necessary design routes, not sufficient counterexample criteria.
In particular, rational three-torsion, extra boundary components, or a larger
Jacobian still has to coexist with scalar units, a one-place branch,
anti-invariant quadratic descent, normal nonmonogenic cubic completion, and
the affine-plane boundary constraints.
