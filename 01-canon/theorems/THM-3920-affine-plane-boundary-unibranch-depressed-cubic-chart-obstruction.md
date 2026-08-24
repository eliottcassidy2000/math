---
id: THM-3920
title: "Affine-plane boundary unibranchness closes every irreducible depressed-cubic radial chart"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Every
  irreducible boundary curve in a normal surface completion of A2 is rational
  and unibranch. Consequently a boundary divisor in a finite-flat degree-d
  completion has at most d normalization addresses over any target point.
  In degree three this closes the affine four-address THM-3906, seven-address
  THM-3913, and eight-address THM-3915/THM-3918 branch packets; the three-place
  projective-infinity packet of THM-3907 is the sharp scope boundary, not an
  affine application of the cap.
  More strongly, for the rational radial chart A=F/4, C=sF/4 with
  F=z^3-3p(s)z+2h(s) irreducible, no same-field affine-plane Keller atlas
  exists. If p is nonsquare, one-address pressure makes the normalized cubic
  cyclic and a character argument forces p=0. If p=q^2 with q nonconstant,
  the two derivative components force incompatible pure-power attachments.
  If q is a nonzero constant, the unique Chebyshev survivor has third-Veronese
  normalization and a deleted total-ramification ray makes a target coordinate
  a forbidden unit. The p=0 pure-Kummer field is cyclic and is excluded by
  THM-3801. Thus this whole chart grammar is closed, not JC(2).
source: jc_zero_debt_lift / post-THM-3917 boundary-address and cyclic-character lane, 2026-08-23
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASS (root and jc_degree6_one_place,
  2026-08-23). The audits independently reconstructed the boundary-tree
  multigraph argument and finite-flat address injection; the nonsquare
  one-address/C3-character contradiction; the square-profile Mason split;
  the exact Chebyshev third-Veronese normalization, ramification rays, and
  forbidden-unit conclusion; and the THM-3801 pure-Kummer seam. The
  assertion-free exact companion verifies the universal
  discriminant/resultant/Jacobian identities, the Chebyshev cover, low-degree
  Kummer controls, and the complete 97-gate symmetric-family stratification.
  Normal and optimized runs byte-match the frozen output; raw hashes and
  documentation checks pass. No repair was required.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3917-quintic-parameter-rational-collapsed-cubic
related:
  - THM-3906-degree-six-common-zero-normal-cubic-two-place-boundary
  - THM-3907-unit-ideal-nonmonogenic-cubic-six-place-boundary
  - THM-3913-moving-triple-root-one-place-decic-normal-nonmonogenic-cubic
  - THM-3915-rational-decic-cube-resolvent-index-debt-euler-tariff
  - THM-3918-rational-moving-triple-root-one-place-decic-and-anchored-color-transition
  - THM-3916-positive-genus-collapsed-valuation-keller-obstruction
  - THM-3887-binary-cubic-common-zero-quintic-one-tangent-obstruction
  - THM-3921-quintic-genus-collapse-decic-degeneration-packet
script: 04-computation/jc2_affine_plane_boundary_unibranch_depressed_cubic_chart_thm3920.py
output: 05-knowledge/results/jc2_affine_plane_boundary_unibranch_depressed_cubic_chart_thm3920.out
script_sha256: aeab7d44be76d9ee50b11ab352701251ad8adc0aee7a85b10efa88b8c281b286
output_sha256: 11e01119ba16805131c3912bd89a96b733cf962d9f511198271d7d561dadd00b
semantic_sha256: d7669931acc942440998ad33d9a378ab5a658e63a34ff603c96fddc5dd653495
hash_basis: raw LF bytes
---

# THM-3920 -- the radial depressed-cubic chart is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero.

This theorem separates three layers that were previously conflated:

```text
boundary topology:  an A2 boundary curve cannot acquire two branches;
finite degree:       a rank-d fibre can carry at most d boundary addresses;
radial cubic chart:  an irreducible collapsed fibre forces all derivative
                     addresses through one point, so the true cap is one.
                                                                    (1)
```

The cap one, rather than genus alone, closes the rational chart that produced
the genus-zero near miss of THM-3917.

## 1. Boundary curves of an affine plane are rational and unibranch

Let `X_0` be a normal surface containing a dense open

```text
U isomorphic to A2.                                                (2)
```

Then every irreducible curve `D_0 subset X_0 minus U` has rational
normalization and is unibranch at every point of `X_0`.

Choose a normal projective completion `Xbar` of `X_0` and a smooth projective
common resolution

```text
                 Y
              /     \
          Xbar       P2,                                           (3)
```

whose two maps restrict to the chosen identity on `U`. Resolve the map to the
standard pair `(P2,L_infinity)` by point blowups supported on the boundary.
The reduced SNC divisor `Y minus U` has a tree as its **incidence
multigraph**: a blowup at a smooth boundary point adds a leaf, while a blowup
at a crossing subdivides one edge. In particular every component of
`Y minus U` is rational. The strict transform of `D_0` is one such component,
so the normalization of `D_0` is rational.

Suppose two distinct normalization branches of `D_0` lay over a point `q`.
On a sufficiently high relative log resolution, the strict transform
`Dtilde_0` meets the fibre over `q` in two distinct branch attachments. The
proper birational fibre is connected because `Xbar` is normal. Hence a path
of contracted boundary components joins the two attachments. Adding the two
edges from that path to the vertex `Dtilde_0` creates a cycle in the boundary
multigraph. If both attachments meet the same exceptional component, the two
parallel edges already form a multigraph two-cycle. Both alternatives
contradict the tree property. Thus `D_0` is unibranch.

Affineness causes no gap here: the point `q` remains in the normal open
`X_0 subset Xbar`, and every component of its resolution fibre is disjoint
from `U`, hence belongs to the same boundary tree.

### 1.1 The finite-flat address cap

Suppose in addition that

```text
pi:X_0 -> A2_y                                               (4)
```

is finite flat of degree `d`. For an irreducible boundary curve `D_0`, let
`D_0^nu` be its normalization. For every target point `y`, unibranchness gives
an injection of sets

```text
(D_0^nu)_y  ->  (X_0)_y.                                    (5)
```

The number of closed points in a finite fibre is at most its scheme length,
and flatness makes that length `d`. Therefore

```text
#(D_0^nu)_y <= #(X_0)_y <= length((X_0)_y)=d.                (6)
```

This is the universal degree-`d` normalization-address cap. It counts
branches, not merely distinct points on a singular plane model.

### 1.2 The cubic branch-address corollary

For a hypothetical degree-three Keller open, THM-3801 gives decomposition
type `(2,1)` over every irreducible affine branch component `Gamma`: the
ramified boundary divisor `E_Gamma` has residue degree one over `Gamma`.
Consequently their normalizations have the same function field. If `Gamma`
has `r` distinct normalization addresses above an affine target point `y`,
then `E_Gamma^nu` has at least the same `r` addresses above `y`. Equation
`(6)` gives the sharp necessary condition

```text
r <= 3.                                                       (6a)
```

No extra `+1` for the unramified companion is available at a special closed
point: the companion can meet the ramified sheet in the finite completion or
lose that point on the etale open. The safe obstruction is exactly `(6a)`.

This yields immediate alternative no-atlas proofs for several existing
finite cubic completions:

| theorem | affine branch packet | address-cap consequence |
|---|---:|---|
| THM-3906 | ordinary `4`-address origin | excluded (`4>3`) |
| THM-3913 | `7`-address origin | excluded (`7>3`) |
| THM-3915 | ordinary `8`-address origin | excluded (`8>3`) |
| THM-3918 | `8`-address origin | excluded (`8>3`) |
| THM-3921 | provisional `6`-address degeneration packet | excluded once that packet's pending audit promotes |

THM-3907 has three normalization places over each of two **projective
infinity** points. It numerically meets the cubic cap, but those are not
points of the affine target `A2`; `(6)` does not close that candidate. This is
the exact equality/scope boundary of the present corollary.

## 2. The depressed-cubic radial chart and its finite completion

Take `p,h in k[s]` and put

```text
F(s,z)=z^3-3p(s)z+2h(s),
r(s)=p(s)^3-h(s)^2.                                       (7)
```

Assume throughout that `F` is irreducible in `k[s,z]`. Define

```text
A=F/4,                       C=sF/4=sA,
K=k(s,z),                    R=k[A,C].                     (8)
```

Since `s=C/A` in the fraction field and `A` has a unique pole of order three
at `z=infinity` over `k(s)`,

```text
[K:Frac(R)]=3.                                              (9)
```

Let `B` be the integral closure of `R` in `K` and `X=Spec(B)`. The ring `B`
is finite over `R`. Normality in dimension two makes `B` Cohen--Macaulay, so
miracle flatness makes

```text
X -> A2_(A,C) finite flat of degree three.                  (10)
```

There is also an everywhere-defined birational morphism

```text
phi:A2_(s,z) -> X.                                          (11)
```

Indeed, an element of `B` is integral over `R`, hence integral over the
normal ring `k[s,z]`; therefore `B subset k[s,z]`.

The exact universal identities are

```text
disc_z(F)=108r,
Res_z(F,z^2-p)=-4r,
F mod (z^2-p)=2(h-pz),
Jac_(s,z)(A,C)=-F F_z/16.                                  (12)
```

Suppose, toward a contradiction, that the same field admitted plane
coordinates `x,y` with

```text
K=k(x,y),                A,C in k[x,y],
Jac_(x,y)(A,C) in k*.                                      (13)
```

The normalization ring `B` lies in the normal ring `k[x,y]`. The induced
birational map `A2_(x,y)->X` is quasi-finite because its composite with the
finite map `(10)` is etale and quasi-finite. Normalization-form Zariski Main
makes it an open immersion. We may therefore apply Section 1 to this
hypothetical plane open. Every intrinsic ramification divisor of `(10)` is
part of its boundary.

## 3. An irreducible derivative curve gives a one-address invoice

First suppose

```text
p is not a square in k[s].                                   (14)
```

Then

```text
D_1:z^2=p(s)                                                 (15)
```

is irreducible. Its image `D subset X` is an intrinsic ramification divisor:
at its generic point the minimal cubic has derivative `F_z=0`, while `F` is
not identically zero there. Moreover `D_1->D` is generically birational. On
`D_1`,

```text
A=(h-pz)/2,
s=C/A,                         z=(h-2A)/p                   (16)
```

away from the indicated proper exceptional set, so the two curves have the
same function field.

The irreducible curve `F=0` maps under `(11)` into the finite fibre of `X`
over `(A,C)=(0,0)`. Its image is therefore one point `o`.

Every finite support point `alpha` of `r` supplies at least one normalization
address of `D_1` above `o`. If `p(alpha)!=0`, the address is uniquely

```text
z=h(alpha)/p(alpha),                                       (17)
```

and `(12)` verifies both `z^2=p` and `F=0`. If `p(alpha)=0`, then also
`h(alpha)=0`; every normalization branch of `(15)` above `(alpha,0)` lies on
`F=0`. Singular residual roots therefore do not evade the invoice--they can
only add addresses.

Distinct `s`-values give distinct normalization addresses. Since `D` must be
unibranch at `o`, Section 1 forces

```text
# Supp_fin(r) <= 1.                                         (18)
```

Irreducibility of `F` also guarantees `r!=0`: a cubic with zero discriminant
over `k(s)` has a repeated factor.

## 4. Two branch values force a pure C3 character

Let `C_F` be the normalization of `F=0`. The projection

```text
C_F -> P1_s                                                 (19)
```

is a connected separable cover of degree three. Its finite branch locus is
contained in `Supp(r)`, so `(18)` says that it is branched over at most one
finite point, with infinity the only other possible branch value.

A connected tame cover of `P1` of degree three cannot be branched over fewer
than two points. With exactly two possible values, the two inertia
permutations are inverse and generate the transitive monodromy group. A
transitive cyclic subgroup of `S3` is `C3`, generated by a three-cycle.
Consequently `(19)` is cyclic, both branch values are totally ramified, and
Riemann--Hurwitz gives `C_F isomorphic to P1`.

After translating and scaling the affine coordinate and choosing a parameter
`t`, write

```text
s=t^3.                                                       (20)
```

Monicity of `F` makes `z` integral over `k[s]`. The normalization of `k[s]`
in `k(t)` is `k[t]`, hence

```text
z=z(t) in k[t].                                             (21)
```

Let `omega` be a primitive cube root of unity. The cubic is depressed, so its
trace equation is

```text
z(t)+z(omega t)+z(omega^2 t)=0.                            (22)
```

Thus `z(t)` contains only powers whose exponents are `1` or `2 modulo 3`.
By `(18),(20)`, the power discriminant is nonzero on `G_m`. Therefore

```text
z(t)-z(omega t)                                             (23)
```

has no nonzero root. It is a unit in `k[t,t^(-1)]`, hence `c t^m`. In `(23)`
every allowed coefficient of `z(t)` is multiplied by the nonzero scalar
`1-omega^j`; consequently `(23)` can be a monomial only if

```text
z(t)=c_0 t^m,                 3 does not divide m.          (24)
```

The minimal polynomial of `(24)` over `k(t^3)` is

```text
Z^3-c_0^3 s^m.                                             (25)
```

Comparing `(25)` with `(7)` gives `p=0`, contradicting `(14)`. Hence no
Keller plane atlas `(13)` exists when the derivative curve is irreducible.

## 5. A split derivative with nonconstant charge also fails

Now suppose

```text
p=q^2,                         q in k[s] nonconstant.       (26)
```

The reduced derivative divisor has two components

```text
D_+:z=q,                         D_-:z=-q.                  (27)
```

They are generically distinct intrinsic ramification divisors. Their
restrictions to the collapsed fibre are

```text
F|_(D_+)=2(h-q^3),               F|_(D_-)=2(h+q^3).        (28)
```

Each image is a boundary curve and all of its intersections with `F=0` map
to `o`. Unibranchness therefore says that each nonzero polynomial in `(28)`
has at most one root. Neither is identically zero, because then `z=q` or
`z=-q` would be a polynomial root of `F`. Thus

```text
h-q^3=c(s-a)^m,                 h+q^3=d(s-b)^n             (29)
```

with `c,d!=0` and `m,n>=0`; exponent zero denotes a nonzero constant.

If the two nonconstant supports have different centres--and a constant may
be assigned any centre--the three terms in

```text
c(s-a)^m+2q^3=d(s-b)^n                                  (30)
```

are pairwise coprime. Indeed, a common factor of `q` and `s-a` would force
the right side to vanish at `a`, hence `a=b`, and similarly at `b`.
Mason--Stothers gives

```text
3 deg(q) <= deg(rad((s-a)q(s-b)))-1 <= deg(q)+1,           (31)
```

impossible for `deg(q)>=1`.

It remains to consider a common centre with `m,n>0`. If, say, `m<n`, then

```text
2q^3=(s-a)^m(d(s-a)^(n-m)-c).                              (32)
```

The second factor has simple nonzero roots in characteristic zero, so `(32)`
cannot be a cube. Hence `m=n`. Equation `(30)` then makes `m` a multiple of
three and

```text
q=q_0(s-a)^(m/3),                    h=lambda q^3.          (33)
```

After `z=qw`, the polynomial `F/q^3=w^3-3w+2lambda` has a constant root over
the algebraically closed field, so `F` is reducible. This final contradiction
closes every nonconstant split derivative.

## 6. The constant split seam is the third Veronese, not an escape

Let `q` be a nonzero constant. Scale it to `q=1`. Equations `(28),(29)` now
say that `h-1` and `h+1` each have at most one root. If `h` were constant,
`F` would split over `k`. Otherwise write the two nonzero forms as in `(29)`
and differentiate their difference `2`. If one exponent were zero, both
derivatives would vanish and `h` would be constant. Hence both are positive,
and

```text
cm(s-a)^(m-1)=dn(s-b)^(n-1).                               (34)
```

Degree comparison gives `m=n`. If this common exponent exceeded one, the
unique derivative roots and leading coefficients would give `a=b,c=d`,
making the original difference zero. Thus `m=n=1`. An affine change of `s`
reduces the only surviving profile to

```text
F=z^3-3z+2s.                                               (35)
```

Its two derivative lines do indeed pass the one-address test:

```text
F(s,1)=2(s-1),                    F(s,-1)=2(s+1).           (36)
```

The finite completion nevertheless exposes a different obstruction. Use the
source coordinates

```text
a=A=F/4,                         z,
```

which are polynomial coordinates because
`s=2a-(z^3-3z)/2`, and make the triangular target change

```text
b=4A^2-2C=a(z^3-3z).                                     (37)
```

Adjoin `t,w` with

```text
a=t^3,                             z=w/t.                  (38)
```

Then

```text
b=w^3-3wt^2,                                                 (39)
```

and `mu_3` acts diagonally by `(t,w)|->(omega t,omega w)`. The integral
closure of `k[a,b]` in `k(a,z)` is exactly the third-Veronese invariant ring

```text
B=k[t^3,t^2w,tw^2,w^3].                                   (40)
```

Indeed `k[t,w]` is integral over `k[a,b]`, `(40)` is normal as a finite-group
invariant ring, and its fraction field is `k(a,z)`.

Upstairs the Jacobian is

```text
Jac_(t,w)(t^3,w^3-3wt^2)=9t^2(w^2-t^2).                   (41)
```

The scalar `mu_3` quotient is etale away from the origin, so the images of

```text
t=0,                         w=t,                         w=-t    (42)
```

are the three intrinsic ramification divisors. A Keller plane open must
delete all three, in particular the divisor `t=0`. But the zero locus in
`X` of the nonconstant target coordinate

```text
a=t^3                                                       (43)
```

is precisely that divisor. Thus `a` restricts to a nonconstant unit on the
hypothetical `A2` open, contradicting `k[x,y]^*=k*`. The minimal Chebyshev
seam is therefore closed.

## 7. The nonreduced p=0 seam is cyclic, not S3

Finally let `p=0`. Over `Frac(R)=k(s,A)`, the cubic relation is

```text
z^3=4A-2h(s).                                               (44)
```

If `F` is irreducible, the right side is not a cube. Since `mu_3 subset k`,
the degree-three field in `(44)` is already a cyclic `C3` Kummer extension.
THM-3801 proves that the finite normalization of any constant-unit
degree-three Keller open has nontrivial quadratic resolvent, `S3` normal
closure, and only `(2,1)` codimension-one inertia. A cyclic cubic with total
ramification cannot occur. Thus `(13)` is impossible here as well.

As a hostile concrete control, for `h=s` the same triangular target change
as above gives `b=az^3`; after `a=t^3,z=w/t`, its normalization is again the
third Veronese and the two coordinate rays are totally ramified. The cited
S3 gate, not this one example, is the uniform proof of Section 7.

Combining Sections 4--7 proves:

```text
for every p,h in k[s] with F=z^3-3pz+2h irreducible,
the radial chart A=F/4, C=sF/4 has no same-field A2 Keller atlas.    (45)
```

## 8. Exact symmetric-family classification

The family that produced THM-3917 is now a corollary of `(45)`, because its
`p` is nonsquare. Its full parameter stratification is nevertheless useful:
it shows exactly how address count changes before the conceptual closure.

Put

```text
p=(s^2-1)(s^2+b)^2,

h=s^9+(3b-3/2)s^7+(3b^2-9b/2+3/8)s^5
 +(b^3-9b^2/2+9b/8+1/16)s^3
 +(-3b^3/2+9b^2/8+3b/16+3/128)s.                         (46)
```

Thus `h` is the polynomial part at infinity of `p^(3/2)`. For `x=s^2`, set

```text
S(x)=16384(p^3-h^2)(sqrt(x)).                              (47)
```

Introduce

```text
f_1=4b+1,                         f_2=8b^2+2b+1,
E_2=48b^2+16b+3,
E_3=64b^3+48b^2+24b+5,
K_5=2304b^5+10176b^4+4064b^3+996b^2+84b+5.                (48)
```

The six polynomials `b,f_1,f_2,E_2,E_3,K_5` are squarefree and pairwise
coprime. Exact elimination gives

```text
S(0)=-16384b^6,                S(-b)=bE_2^2,
S(1)=-E_3^2,

disc_x(S)=-27*2^14 E_2^6 E_3^3 K_5,

disc_s(S(s^2))=-2^57*3^7 b^6 f_1 f_2 E_2^12 E_3^6 K_5^2. (49)
```

Call an address *good* when it is a distinct finite normalization address
with `p!=0`. The exact classification is

| parameter stratum | residual shape | good addresses |
|---|---:|---:|
| none of the six factors in `(48)` vanishes | four simple nonzero `x` roots | `8` |
| `b=0` | one bad root `x=0` plus a squarefree cubic | `6` |
| `f_2=0` | squarefree cubic; the fourth root moved to infinity | `6` |
| `E_2=0` | `(x+b)^2` times a squarefree quadratic | `4` |
| `E_3=0` | `(x-1)^2` times a squarefree quadratic | `4` |
| `K_5=0` | one double good root plus two simple good roots | `6` |
| `f_1=0`, equivalently `b=-1/4` | collapsed `F` reducible | outside `(45)` |

For example, on the two common-zero strata the exact quotient identities in
the corresponding residue fields are

```text
E_2=0:
 S=(16/27)(x+b)^2
   [144(10b+3)x^2-(1392b+531)x+84-128b],

E_3=0:
 S=-16(x-1)^2
   [(192b^2+144b+36)x^2+(48b^2-60b-31)x
     -(60b^2+39b+5)].                                    (50)
```

Their quadratic discriminants are respectively

```text
10125(32b-3),                    2(288b^2+6b-37),          (51)
```

and the companion verifies that neither quotient acquires a zero at
`0,1,-b`. On `K_5=0`, it freezes the exact double root, disjoint quadratic,
and every required noncollision resultant.

Finally, the polynomial-root coefficient chain gives

```text
F is irreducible iff b!=-1/4.                               (52)
```

At the excluded value, with `q=s^3-3s/4`,

```text
F=(z-q)(z^2+qz+q^2-3p).                                   (53)
```

Thus every irreducible member of the symmetric family has at least four
good addresses--and the general one-address theorem already rules it out.

## 9. Scope and next design rule

The theorem closes a representation, not the planar Jacobian conjecture.
Its hypotheses use both radial target coordinates

```text
A=F/4,                         C=sF/4.                     (54)
```

A prospective counterexample must break at least one of the mechanisms that
make `(54)` rigid: the irreducible collapsed fibre, the common radial factor,
the one-variable depressed cubic, or the derivative-divisor attachment
grammar. Nonmonogenic binary-cubic orders remain the natural next place to
look, because they do not supply a global `z` whose conjugate differences
collapse to a single character.
