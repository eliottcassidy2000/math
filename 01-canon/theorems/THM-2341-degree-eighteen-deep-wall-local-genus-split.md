---
id: THM-2341
title: "Degree-eighteen deep-wall local genus split"
status: >
  PROVED + VERIFIED-EXACT. Of the three H_4 square-class orbits on
  THM-2338's deep common-root wall, the rational ratio
  t_0=-361/30618 has a genuine totally ramified three-cycle and
  normalization genus one. It cannot carry a rational Keller
  trajectory. At the two quadratic ratios t_+,t_-, the double
  polynomial-discriminant root is instead a split unramified node; each
  normalization has genus zero. Thus the exact deep-wall residual bank
  falls from three square-class candidates to two rational spectral
  curves. Neither remaining curve is proved to satisfy the Keller
  sidecars, and JC(2) remains open.
source: codex-2026-07-25-deep-wall-local-genus
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2338-degree-eighteen-deep-common-root-wall-hurwitz-quartet
related:
  - THM-2335-degree-eighteen-cyclic-square-class-stratum-empty
script: 04-computation/jc2_degree18_deep_wall_local_genus_thm2341.py
output: 05-knowledge/results/jc2_degree18_deep_wall_local_genus_thm2341.out
script_sha256: e2e2fc35515a06e50631ac84bc2cdd932601a0e8fb009443b8dfa395bb003e24
output_sha256: afda8759089ad093250ea4387c9aaad2a496a66875d61f8ccf42286938a00730
hash_basis: working-tree bytes (LF)
---

# THM-2341 -- one deep-wall ratio is elliptic and two are rational

**PROVED + VERIFIED-EXACT.**

THM-2338 isolates three four-support square-class ratios on

```text
126D=25B^2,                 21W=-20BC.              (1)
```

After weighted normalization `B=1`, they are

```text
t=C^2,

t_0=-361/30618,

t_+=-250/13041 + (500/117369)sqrt(3),

t_-=-250/13041 - (500/117369)sqrt(3).              (2)
```

The square-class calculation alone cannot tell whether its even
discriminant root is true three-cycle ramification or merely the square
of a nonmaximal-order index. Local normalization separates the two:

```text
t_0:   one e=3 point + four e=2 values,    genus 1;

t_+:   one split unramified node
       + four e=2 values,                  genus 0;

t_-:   one split unramified node
       + four e=2 values,                  genus 0.  (3)
```

Thus `t_0` is excluded from any rational Keller trajectory. The two
quadratic ratios remain genuine rational spectral-curve candidates.

## 1. Blow up the false sixfold discriminant at y=0

In THM-2332's depressed coordinate, the spectral equation is

```text
f(v,y)=v^3+p(y)v+q(y)=0,

p=(16/964467)P,
q=(64/703096443)Q.                                 (4)
```

On (1),

```text
P=35y^2(54B+7y^2),

Q=7y^3(1620By+26244C+77y^3).                       (5)
```

The polynomial discriminant contains `y^6`, but this is not
ramification. Put

```text
v=yz
```

and divide (4) by `y^3`. The normalized affine chart is

```text
Phi(z,y)
 =z^3+(16/964467)(P/y^2)z
     +(64/703096443)(Q/y^3)=0.                     (6)
```

Its discriminant is exactly

```text
Disc_z(Phi)
 =-(4096/96873331012983)G(y),                      (7)
```

where THM-2338's residual sextic is

```text
G
 =7889y^6+211680By^4+1047816Cy^3
  +1814400B^2y^2+22044960BCy
  +2916000B^3+178564176C^2.                        (8)
```

At `B=1`, its value at zero is

```text
G(0)=2916000+178564176t.                           (9)
```

For `t=t_0`,

```text
G(0)=810648 !=0.                                   (10)
```

For the two roots of

```text
K(t)=62500+7654500t+199644669t^2,
```

the exact resultant is

```text
Res_t(G(0),K)
 =-295233008801472000000 !=0.                      (11)
```

Hence (6) has three distinct roots over `y=0` at every ratio in (2).
Hensel's lemma gives three smooth branches, each with `y` as a local
parameter. The entire `y^6` factor in the original polynomial
discriminant is an order-index square created by the collapsed
coordinate `v=yz`; it contributes no field ramification.

## 2. The rational ratio is a true three-cycle

Remove the displayed powers of `y` from (5):

```text
P_2=54B+7y^2,

Q_3=1620By+26244C+77y^3.                           (12)
```

Their exact resultant is

```text
Res_y(P_2,Q_3)
 =7715736(361B^3+30618C^2).                        (13)
```

Thus a nonzero common root of `P` and `Q` occurs precisely at the first
factor of THM-2338's sextic discriminant, namely `t=t_0`.

With `B=1`, that common root is

```text
r_0=-486C/19.                                      (14)
```

Direct reduction modulo `C^2=t_0` gives

```text
P(r_0)=Q(r_0)=0,

P'(r_0)= (1837080/19)C !=0,

Q'(r_0)=-(4251528/19)C !=0.                        (15)
```

Put `s=y-r_0`. At `s=0`, equation (4) is `v^3=0`.
All nonleading coefficients of the monic cubic in `v` are divisible by
`s`, while its constant coefficient has valuation exactly one by
`Q'(r_0)!=0`. Therefore (4) is Eisenstein over `C[[s]]`.

There is one point of the normalization over `r_0`, with

```text
e=3,                   ramification contribution 2. (16)
```

This is a genuine three-cycle, not an index square.

## 3. The two quadratic ratios have split nodes

Equation (13) and

```text
361+30618t_+ !=0,
361+30618t_- !=0                                  (17)
```

show that `P,Q` do not have a common root at either quadratic ratio.
At each ratio, THM-2338's exact Euclidean algorithm gives

```text
deg gcd(G,G')=1.                                   (18)
```

Thus `G` has exactly one double root `r` and four simple roots. The
corresponding depressed cubic at `y=r` has one double root `v_0` and a
distinct simple root; in particular `v_0!=0`.

Work over the complete local ring `C[[s]]`, `s=y-r`. Hensel's lemma
first splits off the simple root, leaving a monic quadratic factor. The
discriminant of that quadratic differs from the cubic discriminant by
the square of its resultant with the simple factor, which is a unit.
By (18), its valuation is exactly two:

```text
quadratic discriminant=s^2 times a unit.            (19)
```

Every unit of `C[[s]]` has a square root. The quadratic therefore splits
into two distinct linear branches in `C[[s]][v]`; their difference has
valuation one. Both branches use `s` as local parameter. Together with
the already separated simple branch, the normalization is unramified
over `r`.

This is the missing-coordinate distinction:

```text
same polynomial discriminant exponent 2,

t_0:       triple-root Eisenstein ramification;

t_+,t_-:   split-node index square.                 (20)
```

The square class sees neither difference.

## 4. Riemann--Hurwitz gives the exact genera

All three parameter points have full support, so THM-2332 proves that
the spectral cubic is absolutely irreducible. Its normalization is a
connected degree-three cover of the `y`-line.

THM-2338 gives four simple roots after removal of the one double root.
At each simple discriminant root the inertia is a transposition and the
ramification contribution is one. Equation (7) and (10)--(11) make
`y=0` unramified, and THM-2332's separable infinity cubic

```text
Disc(1127-138915z+1607445z^2-26040609z^3)
 =-153384762202971019112448 !=0                    (21)
```

makes infinity unramified.

The total ramification is therefore

```text
R(t_0)=4+2=6,

R(t_+)=R(t_-)=4.                                   (22)
```

For a connected degree-three cover,

```text
2g-2=3(-2)+R.
```

Consequently,

```text
g(t_0)=1,

g(t_+)=g(t_-)=0.                                   (23)
```

A nonconstant rational Keller trajectory would give a nonconstant map
from `P^1` to the normalization, as in THM-2332. Such a map cannot land
on the genus-one curve. Hence `t_0` is eliminated.

## 5. Frontier effect

The deepest common-root wall now has the exact ledger

```text
three H_4 square-class points                    THM-2338;

one genuine three-cycle / genus-one point        ELIMINATED here;

two split-node / genus-zero points               SURVIVE locally.     (24)
```

The two survivors are Galois conjugate over `Q(sqrt(3))` at the
square-class-ratio level, but they are distinct complex
weighted-projective orbits. Their normalization is rational; no claim
is made that the retained parameterization satisfies the Faber
equations, Keller one-form, first flux, nonsplit deck, or polynomial
origin conditions. Those sidecars are the next exact tests.

The mixed `H_2 S_5^2` and four-simple `H_4 S_4^2` strata away from the
deep wall, other Newton edges, `JC(2)`, and `DC(2)` remain open.

## 6. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_deep_wall_local_genus_thm2341.py
python3 -O 04-computation/jc2_degree18_deep_wall_local_genus_thm2341.py
```

Both transcripts are byte-identical to the stored output. The companion
checks (7), the complete `y=0` etale data, the common-root resultant
(13), the point and derivatives (14)--(15), all three exact
multiple-root gcds, separable infinity, and the Riemann--Hurwitz ledger.
The Eisenstein, complete-local-ring splitting, and
Riemann--Hurwitz arguments are the mathematical proof above, not
computer assumptions. No executable check uses Python `assert`.

Independent audit is pending. QED.
