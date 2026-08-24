---
id: THM-3971
title: "Canonical-debt determinantal completions fail the exact-volume gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every
  m>=1, the algebra generated inside k[x,t] by x, z=1+xt, and zt,...,zt^m
  is a smooth normal affine surface X_m containing A2_(x,t) with one smooth
  A1 boundary D, scalar units, Cl(X_m)=Z[D], and
  div(dx wedge dt)=(m-1)D. Thus m=2 and m=3 pay exactly the simple and total
  cubic canonical-different debts. The m=2 ring is the 2-by-3 determinantal
  surface with minors of [[z,p,x],[p,q,z-1]]. Nevertheless
  H^2_dR(X_m)=k, and dx wedge dt is its nonzero generator: the primitive
  -t dx=(1-z)dx/x has residue one on D. Hence no X_m admits any global
  Darboux pair J(A,C)=1, in any degree. This is a sharp completion near miss,
  not a planar Jacobian counterexample.
source: jc-zero-debt-lift / post-THM-3968 positive completion design, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-degree6-one-place, 2026-08-24). The
  audit rederived the two-chart affine-plane gluing, boundary presentation
  and arc, smoothness, units and class group, canonical order, Gysin rows
  and residue-one volume class, determinantal kernel, bracket table,
  t-not-in-B endpoint, and bounded search. It also independently verified
  the exponent-two THM-3572 exact-volume escape. Normal and optimized
  115-gate runs byte-match the frozen output; hashes and docs pass.
depends_on:
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
related:
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
  - THM-3789-higher-pole-hermite-spectral-completion
  - THM-3968-canonical-vector-different-affine-plane-boundary-obstruction
script: 04-computation/jc2_canonical_debt_determinantal_completion_thm3971.py
output: 05-knowledge/results/jc2_canonical_debt_determinantal_completion_thm3971.out
script_sha256: 8d956b4b3753ca82a177caf957b8cfb59d871ed0e035528501467294a46245d2
output_sha256: 10a52f241fdfecb1bd76f942040e428f750262fd46fe53dd4d920ccca46d5afa
semantic_sha256: 3ef5c52adacd658670473c447c7e5dcad9f435bdb3942269af4631b86074f40d
hash_basis: raw LF bytes
---

# THM-3971 -- canonical debt is paid, but the volume is not exact

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. For an integer
`m>=1`, put

```text
z=1+xt,                    y_j=z t^j       (1<=j<=m),       (1)
B_m=k[x,z,y_1,...,y_m] subset k[x,t],
X_m=Spec(B_m).
```

Then `X_m` is a smooth normal affine surface with an open affine-plane chart
and one smooth boundary curve:

```text
U=D(x) union D(z) isomorphic to A2_(x,t),
D=X_m minus U=V(x,z) isomorphic to A1_(y_m).               (2)
```

Its unit, class-group, and canonical passports are

```text
B_m^*=k^*,              Cl(X_m)=Z[D],
div_Xm(dx wedge dt)=(m-1)D.                               (3)
```

In particular, the canonical vector for `m=2` is `[D]`, exactly the
different vector required by THM-3968 for one simple `(2,1)` cubic
ramification prime. For `m=3` it is `2[D]`, exactly the total-cubic vector.
These are compatibility statements only: no finite map to a target plane is
asserted.

The decisive second invoice is de Rham exactness. One has

```text
H^1_dR(X_m)=0,              H^2_dR(X_m)=k[dx wedge dt],    (4)
```

and the displayed two-form has nonzero class. Consequently there are **no**
elements `A,C in B_m` such that

```text
J_(x,t)(A,C)=1.                                             (5)
```

Thus the entire family is an all-degree no-Darboux near miss. It proves no
case of `JC(2)`; rather, it separates a new necessary gate from the
previously isolated boundary-tree, class-group, and canonical-vector gates.

## 1. The affine-plane chart and the sole boundary

The fraction field of `B_m` is `k(x,t)`, because

```text
t=(z-1)/x.                                                  (6)
```

The ring is therefore a finite-type integral surface. Localizing its two
distinguished generators gives

```text
B_m[x^-1]=k[x,x^-1,t],
B_m[z^-1]=k[x,t,z^-1].                                     (7)
```

The two identifications agree on their overlap. Moreover `x` and
`z=1+xt` have no common zero in `A2_(x,t)`. They therefore glue to the open
immersion and affine-plane identification in `(2)`.

Near the complement, `z-1` is a unit. The elementary identities

```text
x^(m-j)y_m=(z-1)^(m-j)y_j,                   1<=j<=m,
x^m y_m=z(z-1)^m                                             (8)
```

give an exact local presentation

```text
B_m[(z-1)^-1]
 = k[x,z,y_m,(z-1)^-1]/(z(z-1)^m-x^m y_m).                 (9)
```

Indeed `(8)` expresses every lower `y_j` through `x,z,y_m` after inverting
`z-1`, so the map from the right side of `(9)` is onto. Its defining
polynomial is primitive and irreducible: it is linear in `y_m`, and
`x^m` is coprime to `z(z-1)^m`. Both sides are two-dimensional domains with
fraction field `k(x,t)`, so the surjection is an isomorphism.

Since `z-1` is already invertible modulo `(x,z)`, equation `(9)` yields

```text
B_m/(x,z)=k[y_m].                                          (10)
```

Thus the complement in `(2)` is exactly the nonempty prime curve
`D isomorphic to A1`; no hidden codimension-one component is present. An
explicit closure control is the Laurent arc

```text
x=s,      t=-s^-1+(-1)^m q s^(m-1).                       (11)
```

Along `(11)`, `z` and every `y_j` with `j<m` tend to zero, while `y_m`
tends to the arbitrary value `q`. This checks every boundary address in
`(10)` rather than only its generic point.

## 2. Smoothness, units, and class group

The open `U` is smooth. On the boundary chart `(9)`, write

```text
F=z(z-1)^m-x^m y_m.
```

Then

```text
F_z=(z-1)^(m-1)((m+1)z-1),
F_z|D=(-1)^m !=0.                                         (12)
```

Hence `X_m` is smooth along `D` as well, and therefore smooth and normal
everywhere.

Now apply THM-3922 to the open immersion `(2)`. Because the open is the
ordinary affine plane and `D` is its only prime divisorial complement, its
localization sequence gives exactly

```text
B_m^*=k^*,                 Cl(X_m)=Z[D].                   (13)
```

In particular `D` is primitive. The family passes both the unit and free
boundary-basis invoices; the later obstruction does not come from a hidden
class relation or a punctured boundary component.

## 3. The exact canonical vector

On `U`, let

```text
eta=dx wedge dt.                                           (14)
```

Differentiating the top generator at fixed `x` gives

```text
partial y_m/partial t
 =t^(m-1)((m+1)z-1).                                      (15)
```

Using `(6)`, equation `(15)` rewrites near `D` as

```text
eta=
 x^(m-1) / ((z-1)^(m-1)((m+1)z-1))  dx wedge dy_m.        (16)
```

The denominator in `(16)` is a unit along `D`, and `x` is a local equation
for `D`. On `U`, `eta` has neither zero nor pole. It follows that `eta`
extends to a global regular two-form and has the exact divisor

```text
div_Xm(eta)=(m-1)D.                                       (17)
```

This proves `(3)`. Notice the equality, not only the class: no unlisted zero
or pole can occur because `U` and `D` exhaust `X_m`.

The values `m=2,3` therefore reproduce the two positive tame cubic entries
of THM-3968. The endpoint `m=1` has trivial canonical module, so THM-3968
already prevents it from being a nontrivial finite Keller completion. The
next section is stronger for this family: it excludes even a Darboux pair
without assuming finiteness.

## 4. The residue-one de Rham obstruction

Use the Gysin localization sequence for the smooth divisor `D` in the
smooth surface `X_m`. Since

```text
U isomorphic to A2,       D isomorphic to A1,
H^i_dR(A2)=H^i_dR(A1)=0             for i>0,              (18)
```

the degree-one and degree-two rows reduce to

```text
H^1_dR(X_m)=0,
0=H^1_dR(U) -> H^0_dR(D) -> H^2_dR(X_m) -> H^2_dR(U)=0.   (19)
```

Thus the Gysin arrow is an isomorphism

```text
k=H^0_dR(D) isomorphic to H^2_dR(X_m).                    (20)
```

The class of `eta` can be identified without a dimension-only argument.
On `U`, take the primitive

```text
alpha=-t dx=(1-z) dx/x.                                   (21)
```

Although the second expression appears to have poles on `V(x)`, its factor
`1-z` cancels the interior component `x=0,z=1`; only the boundary `D`
remains logarithmic. Since `x` is a local equation for `D`, its residue is

```text
Res_D(alpha)=(1-z)|D=1.                                   (22)
```

Direct differentiation fixes the sign:

```text
d alpha=dx wedge dz/x=dx wedge dt=eta.                    (23)
```

In the logarithmic de Rham realization of `(19)`, the connecting/Gysin map
sends the residue in `(22)` to the class of `(23)`. Therefore

```text
[eta]=Gysin(1) !=0.                                       (24)
```

This also supplies a direct hostile check on `(20)`: the canonical form is
the generator, not the zero class.

If `(5)` held for global regular functions on `X_m`, then on the dense open
`U`

```text
dA wedge dC=eta=d(A dC).                                  (25)
```

Both sides of `(25)` are global regular two-forms, so equality on `U`
implies equality on all of `X_m`. Equation `(25)` would make `[eta]=0`, in
contradiction with `(24)`. This proves the all-degree no-Darboux assertion.

The proof exposes a reusable necessary condition. If a smooth affine
surface `X` contains `A2_(x,t)` and its source volume `dx wedge dt` extends
regularly to `X`, then any global pair with `J_(x,t)(A,C)=1` forces

```text
[dx wedge dt]=0 in H^2_dR(X).                             (26)
```

Canonical effectivity or the correct canonical divisor is not enough; the
actual volume class must be exact.

## 5. The `m=2` determinantal passport

For `m=2`, rename

```text
p=zt,                       q=zt^2.                        (27)
```

Then

```text
B_2 = k[x,z,p,q]/I_2([[z,p,x],[p,q,z-1]]),                (28)
```

or explicitly the kernel is generated by

```text
zq-p^2,              z(z-1)-xp,              p(z-1)-xq.  (29)
```

To see that `(29)` is the full kernel, localize first at `z`. There
`t=p/z`, and the quotient is

```text
k[x,t,(1+xt)^-1].                                         (30)
```

Localize instead at `z-1`. Eliminating `p` leaves

```text
k[x,z,q,(z-1)^-1]/(x^2q-z(z-1)^2),                        (31)
```

which is exactly the `m=2` case of `(9)`. The opens `D(z)` and `D(z-1)`
cover because `z-(z-1)=1`, so the proposed quotient and `B_2` are globally
isomorphic. This proves the determinantal presentation without a
set-theoretic or saturation gap.

The generator brackets are

```text
J(x,z)=x,                 J(x,p)=1+2xt,
J(x,q)=t(2+3xt),          J(z,p)=p,
J(z,q)=2q,                J(p,q)=zt^3.                    (32)
```

They show why this surface looked counterexample-positive: constants occur
locally in the bracket table, the boundary is a single primitive `A1`, and
`K_X=[D]` has the desired simple-cubic color. Equation `(24)` is the first
global invoice it fails.

## 6. The exact finite-completion search and bounded hostile control

The strongest search question on `B_2` would be to find `A,C in B_2` with

```text
J_(x,t)(A,C)=1,
B_2 finite over k[A,C].                                   (33)
```

Because `B_2` is normal, `(33)` would make `X_2` the actual finite normal
completion of the resulting plane map. Even without the finiteness row, a
pair satisfying the first equation would already be a planar Jacobian
counterexample. Indeed, at the generic point of `D`,

```text
ord_D(x)=1,           ord_D(z-1)=0,        ord_D(t)=-1.    (34)
```

Thus `t` does not lie in `B_2`. If `k[A,C]=k[x,t]`, the inclusions
`k[A,C] subset B_2 subset k[x,t]` would force `t in B_2`, contrary to
`(34)`. The de Rham argument proves that the first equation in `(33)` is
already empty.

The exact companion retains a bounded search as a hostile control, not as
the proof. Give `x,z,p,q` filtration one and let `F_N` be the span of their
monomials of total generator degree at most `N`. Exact row reduction gives

```text
dim(F_0),...,dim(F_4)=1,5,12,22,35.                        (35)
```

It then solves the full linear mate problem `J(A,C)=1`, with `C in F_4`,
for:

1. all forty nonzero coefficient rows in `{-1,0,1}^4`, modulo overall
   sign, for `A` in the span of `x,z,p,q`; and
2. the eleven nonconstant members of the deterministic monomial basis

   ```text
   q,p,z,x,q^2,pq,p^2,zp,z^2,xz,x^2                       (36)
   ```

   of `F_2`.

There are no survivors in either census. These 51 null rows independently
match `(24)`, but no finite cutoff is used to infer the all-degree theorem.

Reproduce with

```bash
python3 04-computation/jc2_canonical_debt_determinantal_completion_thm3971.py
python3 -O 04-computation/jc2_canonical_debt_determinantal_completion_thm3971.py
```

Both runs must print `CHECKS=115` and byte-match the frozen output.

## 7. Scope and design consequence

This theorem does not say that every smooth affine-plane completion with
nonzero `H^2_dR` is impossible as a finite cover; the conclusion concerns a
source volume whose restriction is the Keller volume and whose global class
is nonzero. It also does not say that the boundary class and canonical
vector are irrelevant. THM-3922 and THM-3968 remain necessary and eliminate
many candidates before de Rham cohomology is computed.

The exact lesson is a conjunction of invoices:

```text
boundary tree/A1 + free primitive boundary basis + correct canonical vector
             does not imply a Keller affine-plane atlas;
the extended source volume must additionally be exact.                  (37)
```

For this one-boundary family the independence can be recorded as two scalar
coordinates attached to the same form:

```text
kappa=ord_D(eta)=m-1,             rho=[eta]/Gysin(1)=1.    (38)
```

The first is the canonical/different coordinate; the second is the
logarithmic-residue coordinate. A finite cubic design asks for
`kappa=delta`, while a Keller pair asks for `rho=0`. Thus `X_2` has
`(kappa,rho)=(1,1)`: it pays the simple cubic canonical color and fails the
independent exactness color. This is more informative than recording only
`dim H^2_dR(X_m)=1`.

The residue-one failure is not automatic for an `A1` boundary basis.
THM-3572 supplies the sharp route-around. On

```text
Y_Sigma={c^2e=Sigma(b)},
```

retain one line over a root `beta_1` in the affine-plane chart and delete
the other root lines. With `b=beta_1+c^2t`, the physical volume becomes

```text
db wedge dc/c^2=dt wedge dc=d(t dc).                      (39)
```

Its boundary primitives have no logarithmic residue even though the
boundary components are affine lines and their classes give the expected
free basis. In other words, changing the affine modification from a simple
logarithmic pole to an exponent-two principal part can kill the volume class
without destroying the boundary passport. This is an actual design route,
not a claim that every such modification supports a single Darboux pair.

The `m=2` surface is therefore a deliberately sharp hostile: it is smooth,
rational, determinantal rather than globally complete-intersection, has the
right single simple-cubic canonical debt, and still fails because the
logarithmic residue is one. A positive completion must change the gluing so
that its canonical volume loses this residue while preserving the boundary
and different passports. **QED.**
