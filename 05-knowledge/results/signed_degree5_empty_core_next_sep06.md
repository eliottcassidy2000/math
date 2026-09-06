# The sharp degree-five signed duplication margin

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED (root).**
Bounded wildcard continuation at inherited `6dd59c9c41`, 2026-09-06.
No theorem ID, Git mutation, or external-priority claim.

## 1. Result and recovered boundary

Let `G(s)=a_0 prod_(i=1)^5(1+r_i s)` be a real polynomial of exact degree
five, with `a_0` real and nonzero,
and only real roots. Thus every reciprocal-root factor `r_i` is real and
nonzero. Suppose `[s^2]G=0`, equivalently `e_2(r)=0`. Put

```text
E=a_0^2 e_2(r_1^2,...,r_5^2)>0,
D=-[s^4]G^2,
c_5=(81-8sqrt(6))/87=0.7057940466406273... .
```

Then the sharp inequality is

```text
D >= c_5 E >0.                                             (1)
```

Equality holds exactly when, up to permutation and a common nonzero real
factor, the reciprocal roots are

```text
(1,1,-t,-t,-t),    t=1-sqrt(6)/3,    1/t=3+sqrt(6).          (2)
```

Thus the extremal polynomial has one repeated root of multiplicity two
and an oppositely signed root of multiplicity three. It can be represented
as `G(s)=(1+s)^2(1-ts)^3`. In particular the equality is attained at exact
degree five, away from every zero-factor boundary.

The closest mechanism is
[THM-4440 — signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md),
whose terminal square gives `D/E>=1/3`. The immediately preceding
[sharp quartic margin and pair stability](signed_duplication_stability_empty_core_sep06.md)
gives `7/9` in degree four and already records the degree-five hostile
`(1,1,-1/6,-1/6,-13/60)`, whose margin is `18501/26101<7/9`.
The new constant lies strictly between `1/3` and `7/9`, and below that
hostile's margin.

Targeted incoming/repository searches recovered no degree-five margin
theorem. Incoming
[endpoint-27 quartic-channel return](overnight7_20260906_laurent_quartic_carry.md)
was inspected: its degree refers to a channel polynomial, with lower carry,
and does not assert real-rootedness of the actual ordinary carrier used here.
No conclusion from that different object is transported.

The live board is: subset-product energy; two moment constraints; stationary
root values; multiplicities and tangent curvature; zero-factor boundaries;
real versus complex gauge. The least-used sidecar is the multiplicity of
each stationary root value. The successful map is

```text
ordinary square coefficient -> normalized power sums
 -> compact two-moment sphere -> stationary cubic
 -> multiplicity-sensitive tangent directions -> four exact value cases.
```

It preserves the global minimum of a polynomial objective, and retains the
root multiplicities. Merely listing the cubic's three possible root values
would lose the negative tangent direction that excludes three-value minima.

## 2. Normalize the cancellation before taking a minimum

Divide out `a_0`. Write `p_j=sum_i r_i^j`. Since
`2e_2=p_1^2-p_2=0`, a nonzero list has `p_1!=0`. Replacing every `r_i`
by `r_i/p_1` changes `D` and `E` by the same positive fourth-power factor.
We may therefore assume

```text
sum_i r_i=1,              sum_i r_i^2=1.                    (3)
```

Direct elementary-symmetric expansion, or Newton's identities through
degree four, gives

```text
E=(1-p_4)/2,
D=(5-8p_3+3p_4)/6,
D/E=(5-8p_3+3p_4)/(3(1-p_4))       when E>0.                (4)
```

The identities retain the actual ordinary-square coefficient. In
coefficient coordinates they are also
`E=-2a_1a_3+2a_0a_4`, `D=-2a_1a_3-2a_0a_4` before dividing out `a_0`.
They hold in higher degree too; it is the optimization universe below
that depends on degree five.

Rather than minimize a quotient with a vanishing boundary denominator,
for a proposed constant `c>-1` minimize the polynomial

```text
F_c(r)=3(1+c)p_4-8p_3+5-3c=6(D-cE)                         (5)
```

on the full compact sphere cut out by (3), **including zero entries**.
The only zero-energy points there are permutations of `(1,0,0,0,0)`;
indeed `p_4<=p_2^2=1`, with equality precisely when at most one entry
is nonzero. At those points `F_c=0` for every `c`. Discarding this boundary
before optimizing a quotient would be unjustified.

## 3. An elementary two-value minimum lemma

For every `n>=4` and every `A>0`, the function

```text
Phi(r)=A sum_i r_i^4-8 sum_i r_i^3
```

attains its minimum on `sum r_i=sum r_i^2=1` at a list with exactly two
distinct values. In fact **every local minimum** has at most two distinct
values. This is the sole extremal reduction needed here.

Proof. The constraint surface is compact and smooth: the gradients of
`sum r_i` and `sum r_i^2` are independent, since a constant list cannot
satisfy both constraints for `n>1`. At a local minimum, Lagrange multipliers
give one cubic equation for every occupied value:

```text
P(u)=4A u^3-24u^2-2mu*u-lambda=0.                          (6)
```

The constrained Hessian is diagonal with entries `P'(r_i)`. Thus there
are at most three distinct values. If three occur, write them `x<y<z`.
Then `P(u)=4A(u-x)(u-y)(u-z)`, so `P'(y)<0`. If the middle value occurs
twice, the vector supported as `+1,-1` on two such coordinates is tangent
to both constraints and has strictly negative second variation. This is
impossible at a minimum.

It remains to exclude a singleton middle value. Let the multiplicities
of `x,z` be `a,b>=1`; then `a+b+1=n>=4`. Put `u=y-x>0`, `v=z-y>0`.
Give every `x` coordinate velocity `b v`, the `y` coordinate velocity
`-ab(u+v)`, and every `z` coordinate velocity `a u`. Both the velocity sum
and its scalar product with the root list vanish, so it is a nonzero tangent
vector. Its constrained second variation is exactly

```text
-4Aab*u*v*(u+v) [a(b-1)u+b(a-1)v] <0.                     (7)
```

The strict inequality uses `a+b>=3`. This excludes the final three-value
case. One occupied value is incompatible with the constraints, so a
minimum has exactly two. The usual nonnegative constrained-Hessian
condition applies because the constraint surface was checked smooth;
zero coordinates do not create extra faces or exceptions. This proves
the lemma.

The existence of a two-value minimum also follows from the directly
eligible [Riener Lemma 4.2](https://arxiv.org/pdf/1001.4464), with `e_1,e_2`
fixed and the objective linear in `e_3,e_4`. The
[all-degree addendum](signed_degree5_empty_core_next_sep06_uniform.md)
records this degree-principle antecedent and the precise parameter map.

The lemma is an analytic finite reduction, not a numerical inference from
root sampling. Its dimension-five use below involves only four exact
orbits. No general degree-dependent best-constant formula is asserted here.

## 4. Complete degree-five orbit evaluation and equality

If one of the two values has multiplicity `m`, the other has multiplicity
`5-m`. Interchanging their names permits `m=1,2`. Solving (3) gives

```text
x=(1 +/- sqrt(4(5-m)/m))/5,
y=(1 -/+ sqrt(4m/(5-m)))/5.                              (8)
```

The complete list is:

| Multiplicity of `x` | `x` | `y` | `D/E`, or boundary |
|---:|---|---|---|
| 1 | `1` | `0` | `E=D=0` |
| 1 | `-3/5` | `2/5` | `7/3` |
| 2 | `(1+sqrt(6))/5` | `(3-2sqrt(6))/15` | `(81-8sqrt(6))/87` |
| 2 | `(1-sqrt(6))/5` | `(3+2sqrt(6))/15` | `(81+8sqrt(6))/87` |

Every nonzero energy in the table is positive, and the third row has
strictly the smallest ratio. Apply the minimum lemma to (5), with
`A=3(1+c_5)>0`. Its value is nonnegative in every table row, including
the zero-energy boundary; therefore it is nonnegative everywhere on
the constraint sphere. This proves (1).

If equality holds with `E>0`, its point is a global minimum of `F_(c_5)`.
The lemma therefore forces two values, and the table forces the third
row. Their opposite-sign magnitude ratio is `3+sqrt(6)`. Undoing the
normalization yields exactly (2), including all permutations and common
negative scalings. Conversely direct substitution in (2) attains equality.
This classifies all equality cases rather than only producing a candidate.

## 5. Boundaries, gauges, and the exact consumer

Zero entries in the reciprocal-root list reduce the degree of `G`; they
are the compactification boundary corresponding to ordinary roots escaping
to infinity. They are retained in the proof. If energy is positive there,
the inequality is strict, since every positive-energy equality list (2)
has five nonzero entries. With four nonzero entries the inherited sharper
quartic constant `7/9` is consistent with this boundary. The one-nonzero
case has `D=E=0` and no defined normalized ratio.

Zero **ordinary** roots are a different operation: write `H=s^ell G` and
shift the coefficient index to `k=ell+2`. The effective core must have
degree five for the sharp theorem above. As in the quartic note, define

```text
J(H,k)=|[s^(2k)]H^2| /
       (|G(0)|^2 e_2(|r_1|^2,...,|r_5|^2)).                (9)
```

If the core admits a real gauge, (1) applies and `J>=c_5` at its vanishing
indicated coefficient. Under `H(s)->lambda*s^b H(sigma*s)`, with nonzero
complex `lambda,sigma`, integer `b`, and `ell+b>=0`, shift `k->k+b`.
Both numerator and denominator in (9) scale by
`|lambda|^2 |sigma|^(2k)`. Thus the magnitude statement and equality shape
are gauge invariant. The real negative sign belongs to the declared real
gauge; arbitrary non-real-rooted cores remain outside this theorem.

The actual Laurent example

```text
f(u)=u^(-2)(1+u)^2(1-tu)^3,       t=1-sqrt(6)/3,
```

has `CT(f)=0`, `CT(f^2)<0`, and normalized ratio exactly `c_5` by literal
coefficient extraction. This uses the actual ordinary core, with complete
endpoint support. No general channel-polynomial or LRC positivity transfer
is inferred.

The inherited Hermite `k=2` limit has normalized ratio `5/3` at a first
coefficient zero, so it does not supply the new sharp finite-degree margin.
The canonical degree-five hostile is retained unchanged and passes (1).
The analytic minimum argument, rather than extrapolation from that single
hostile or from the Hermite limit, provides the improvement.

## 6. Exact controls and audit manifest

[Source](../../04-computation/signed_degree5_empty_core_next_sep06.py) and
[output](signed_degree5_empty_core_next_sep06.out).

```bash
python3 -B 04-computation/signed_degree5_empty_core_next_sep06.py
python3 -B -O 04-computation/signed_degree5_empty_core_next_sep06.py
```

The standard-library producer uses exact rationals and `Q(sqrt(6))`.
It checks all four stationary value orbits, the literal sharp factor list,
and 81 exact negative-curvature controls for the three possible singleton-
middle multiplicity pairs `(1,3),(2,2),(3,1)`. These computations verify the
written identities; the full Hessian argument proves their unbounded use.

The rational hostile universe is all 126 four-entry multisets from
`{-3,-2,-1,1,2,3}`. Solve for a fifth entry using `e_2=0`, reject 16
undefined or zero-fifth cases, and keep all 98 distinct sorted five-entry
lists. A separate boundary universe uses all 35 four-entry multisets from
`{-2,-1,0,1,2}` containing zero: four are infeasible, and 23 distinct
boundary lists remain, including the all-zero list. Each list is tested
by literal factor multiplication and ordinary squaring, root-product
energy, moment normalization when eligible, scalar/variable gauges, and
monomial index shifts. The inherited `18501/26101` hostile is additional.

All **1,341 explicit gates** pass, and normal/optimized outputs agree.
These small exact controls are not the proof of the sharp constant.
Frozen raw-byte hashes:

```text
source SHA256 6cd33af5ddc4ea90976a37009e498d17f6eddc7c13a83e0c29cb64fd238cd608
output SHA256 88056219b75e78a3609f0de878c166219184ae49121e6917bc035106f5114b46
trace  SHA256 2cd366878ba3d4aab4977b4e3ad79f7a2db8ccd6457f768e188412de01705137
```

The root independently audited the complete normalized Newton identities,
compact polynomial objective, smooth constraint surface, both negative
tangent directions, four value orbits, equality classification, zero-energy
boundary, and gauge normalization: **PASS**. No mathematical correction was
required, and the source/output remained unchanged.

The independent `certificate_audit` referee separately passed the full
degree-five proof and replayed its exact source, checking the frozen hashes.
Their sole clarification was to state the real polynomial/scalar gauge
explicitly in the opening hypotheses; it is incorporated above.
