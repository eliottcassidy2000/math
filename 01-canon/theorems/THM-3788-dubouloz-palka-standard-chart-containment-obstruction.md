---
id: THM-3788
title: "Dubouloz--Palka factorizations cross every standard plane chart"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every
  alpha-zero Dubouloz--Palka etale endomorphism of S(3,3,1), including every
  allowed degree 6N+3 and every source/target automorphic translate, factors
  through the universal cover with image crossing at least two components of
  its special fibre.  The factor image is therefore contained in no standard
  affine-plane chart.  This blocks the standard pseudo-cover chart from
  turning the known generalized-Jacobian counterexamples into planar Keller
  maps; nonstandard plane opens and non-equivariant maps remain open.
source: jc_sparse_direct_search / pseudo-plane chart-entry lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  Equations (8)--(10)
  were checked against the primary text of Dubouloz--Palka, and the standard
  chart inverse, actual t=1 and R2-root fibres, triple valuations, N=0 seam,
  and full target-automorphism argument were rederived.  The exact
  consequence of essential uniqueness of the ruling is now written out in
  Section 4.  Equations (8)--(10) are imported exactly
  from Dubouloz--Palka, arXiv:1701.01425, Theorem D and Lemma 5.1.  The proof
  checks every standard component, the N=0 boundary, all simple R2 roots,
  contact multiplicities, and source/target automorphism scope.  The exact
  companion verifies the degree packet through N=512, the genuine N=0
  cyclotomic control, pairwise separability resultants, the component-choice
  census, and hostile endpoint/simple-contact deletions.  Normal and optimized
  runs byte-match the frozen transcript.
depends_on: []
related:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
script: 04-computation/jc2_dubouloz_palka_standard_chart_obstruction_thm3788.py
output: 05-knowledge/results/jc2_dubouloz_palka_standard_chart_obstruction_thm3788.out
script_sha256: cc77519095810b56ebbc1367b53b6f14994b98c24032385b47ffdd5ce6a48466
output_sha256: 2c8c95174c91fa5f0474c9a192aa0500d73a05c8e8a41e4cf490901237999d75
semantic_sha256: c27b57e1ead3a8277566076a1ea7e152811dd41f037a9f370801f1b532a347bf
hash_basis: raw LF bytes
---

# THM-3788 -- the known pseudo-plane counterexamples cannot enter one standard plane chart

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This
theorem isolates a precise failure of the most direct attempt to turn the
known generalized-Jacobian counterexamples on the cubic pseudo-plane into a
counterexample on the affine plane.  Their lift to the universal cover is
etale, but it necessarily crosses several components of the split special
fibre.  No single standard `A2` chart contains it.

Work over `C`.  Put

```text
T=C[X,Y,Z]/(X^3Y-Z^3+1),                              (1)
B=C[u,v,w]/(u(1+uv)-w^3),                             (2)
```

and let

```text
pi:T -> B,             (X,Y,Z) |-> (X^3,Y,XZ)         (3)
```

be the universal cyclic cover of `B=S(3,3,1)`.  Its deck action is

```text
epsilon.(X,Y,Z)=(epsilon X,Y,epsilon^(-1)Z).          (4)
```

For `omega in mu_3`, write

```text
L_omega=V(X,Z-omega) subset T.                        (5)
```

The special fibre `X=0` is the disjoint union of the three affine lines
`L_omega`.  The standard affine-plane chart retaining `L_omega` is

```text
U_omega=T minus union_(gamma in mu_3, gamma!=omega) L_gamma ~= A2.  (6)
```

These are the charts obtained by deleting all but one component of the split
fibre in the standard affine pseudo-cover construction.  Explicitly,

```text
q_omega=(Z-omega)/X^3,
Z=omega+X^3q_omega,
Y=3omega^2q_omega+3omega X^3q_omega^2+X^6q_omega^3,  (6a)
```

so `(X,q_omega)` identifies `U_omega` with `A2`.  This also proves directly
that the complement in `(6)` is exactly the union of the other two lines.

Let `eta:B -> B` be a `C*`-equivariant nonproper etale endomorphism respecting
the quotient `A1_*`-fibration and lying in the `alpha=0` branch of
Dubouloz--Palka's Theorem D.  Then every factorization

```text
eta=pi o j,                     j:B -> T,              (7)
```

given by their Lemma 5.1 satisfies

```text
j(B) is not a subset of U_omega             for every omega in mu_3. (8)
```

The same conclusion holds after arbitrary precomposition by an automorphism
of `B` and arbitrary postcomposition by an automorphism of `T`.  In
particular it holds for every deck translate and for the deformations
`Theta_P o j` used in Dubouloz--Palka's Theorem C.

## 1. The exact imported factorization packet

The cited theorem applies to

```text
S(k,r,a)=S(3,3,1),                    alpha=0.          (9)
```

Its degree congruence and degree formulas say, for a unique `N>=0`,

```text
d=deg(eta)=6N+3,
deg R2=N,             deg R1=2N+1,       deg R0=3N+1. (10)
```

The polynomials obey

```text
1-R1(t)^3=t(1-t)R0(t)R2(t)^3,                       (11)
R1(0)=R2(0)=1,       R0(0)!=0,                       (12)
(1-t)R0(t)R1(t)R2(t) is separable.                   (13)
```

Writing

```text
t=-uv,                                                (14)
```

and absorbing the harmless nonzero scale into `lambda`, Lemma 5.1 gives

```text
j(u,v,w)=
 (lambda wR2(t), lambda^(-3)vR0(t), R1(t)).           (15)
```

For reference, `(10)` follows directly from the cited general formulas:

```text
d2=(d-r)/(k(r-1))=N,
d1=d2(r-1)+r/k=2N+1,
d0=(kd2+1)(r-1-r/k)=3N+1.                            (16)
```

No existence or classification assertion beyond the hypotheses of Theorem D
and Lemma 5.1 is imported.

## 2. Every possible chart forces the same component address

Equation `(2)` and `(14)` give

```text
w^3=u(1-t).                                           (17)
```

At `u=w=0` one has `t=0`, and `(12),(15)` give

```text
j(0,v,0)=(0,lambda^(-3)R0(0)v,1).                    (18)
```

Since `R0(0)!=0`, the arm `V(u,w)` maps onto the whole line `L_1`.  Therefore
containment in `U_omega` is already impossible for `omega!=1`.  It remains
only to exclude `U_1`.

At `t=1`, equation `(17)` gives the nonempty curve

```text
(u,v,w)=(u,-u^(-1),0),                 u in C*.       (19)
```

Its image has first and third coordinates

```text
(X,Z)=(0,R1(1)).                                      (20)
```

Equation `(11)` gives `R1(1)^3=1`.  If `j(B)` were contained in `U_1`, the
only allowed component above `X=0` would be `L_1`, so necessarily

```text
R1(1)=1.                                             (21)
```

Now let `beta` be a root of `R2`.  By `(12),(13)`, the `N` such roots are
simple, mutually distinct, and disjoint from `0,1` and from the roots of
`R0R1`.  The fibre `t=beta` is nonempty: choose `u in C*`, put
`v=-beta/u`, and choose `w` with `w^3=u(1-beta)`.  Equation `(15)` again has
first coordinate zero there, while `(11)` gives `R1(beta)^3=1`.  Containment
in `U_1` would force

```text
R1(beta)=1               for every beta in V(R2).     (22)
```

This is the complete component-address ledger.  It uses actual points of
`B`, not merely roots of the base polynomial.

## 3. Triple contact exceeds the available degree

Factor

```text
1-R1^3=(1-R1)(1+R1+R1^2).                            (23)
```

At a root `beta` of `R2`, the right side of `(11)` has order exactly three by
separability.  Under `(22)`, the second factor of `(23)` specializes to `3`,
so

```text
ord_beta(R1-1)=3.                                    (24)
```

Similarly, `(11)--(13),(21)` show that `R1-1` has a simple zero at each of
`t=0` and `t=1`.  All these points are distinct.  Consequently

```text
deg(R1-1)>=1+1+3 deg(R2)=3N+2.                       (25)
```

But `(10)` gives

```text
deg(R1-1)=deg R1=2N+1,                               (26)
```

because `deg R1>0`.  The strict inequality `3N+2>2N+1` is the desired
contradiction.

The smallest degree is included rather than hidden in the general notation.
For `N=0`, `R2` is constant, while `R1-1` would have the two distinct roots
`0,1` despite `deg R1=1`.  Indeed, with a primitive cube root `zeta!=1`, the
degree-three positive control is

```text
R1=1+(zeta-1)t,
R2=1,
R0=3(2 zeta t+t-zeta+1),                             (27)
```

which satisfies `(11)--(13)` and has `R1(1)=zeta`, visibly landing on a
second component rather than `L_1`.

There is also a useful planar description of this smallest failure.  Put

```text
s=x^3y,
F=x(1+s),
G=y(s^2+3s+3)/(3(1+s)^3).                            (27a)
```

Then

```text
G=1/(3x^3)-1/(3F^3),                 J(F,G)=1.        (27b)
```

The target shear by the function of `F` cancels the pole along `x=0`, where
`G` specializes to `y`, but transfers a third-order pole to the other factor
`1+s=0` of `F=0`.  This is the rational shadow of the two component addresses
`R1(0)=1` and `R1(1)=zeta`: the debt is moved to the wrong special arm rather
than removed.

## 4. Automorphisms do not restore containment

Precomposing `j` by any automorphism of `B` does not change its image.  For
the target, the affine ruling `X:T -> A1` is essentially unique because the
Danielewski exponent is `3>1`; this is the uniqueness statement immediately
following Flenner--Zaidenberg's Theorem 1.1.  Hence every automorphism of `T`
gives a ruling `X o Phi` with the same general fibres, so `X o Phi=h(X)`.
A general fibre of `X o Phi` is connected, forcing `deg h=1`; and its
reducible fibre is the inverse image of `X=0`, forcing the translation term
of `h` to vanish because `X` has only one reducible fibre.  Thus every
automorphism preserves `X=0` and permutes its three irreducible components
`L_omega`.  From definition `(6)`, it therefore permutes the three standard
charts.

If `Phi in Aut(T)` and `sigma in Aut(B)` satisfied

```text
Phi o j o sigma(B) subset U_omega,                    (28)
```

then `j(B)` would be contained in the standard chart
`Phi^(-1)(U_omega)`, contradicting `(8)`.  This includes deck transformations
and the explicit triangular automorphisms `Theta_P` of the universal cover.

There is an independent conceptual cross-check from the arm geometry.  On
`U_1`, a standard second plane coordinate is

```text
q=(Z-1)/X^3=Y/(Z^2+Z+1),              q|L_1=Y/3.      (29)
```

Equation `(18)` would therefore give

```text
(q o j)(0,v,0)=lambda^(-3)R0(0)v/3,                  (30)
```

an injective linear restriction on the arm.  THM-3790 proves that every
Darboux morphism from `B` to `A2` restricts noninjectively on that arm.  The
etale factor has constant nonzero symplectic multiplier because `B` has only
constant units, so a harmless rescaling of `q` makes the induced plane map
Darboux without changing injectivity.  A hypothetical containment in `U_1`
would therefore also fail its arm gate.  This is a second proof route and is
not used in the Shabat degree argument above.
Moreover `Theta_P` sends `q` to `q+P(X)`, so it preserves the linear arm
restriction.  THM-3790 is independently hostile-audited, but the present
theorem does not depend on it.

## 5. Consequence and exact boundary

Dubouloz--Palka construct genuine nonproper etale self-maps of `B`, but this
theorem proves that their `alpha=0` universal-cover factors cannot be captured
inside one of the canonical copies of `A2` obtained by retaining a single
special-fibre component.  The obstruction is a component-address deficit:
the cube on `R2` demands triple return to one component, while `R1` has only
two-thirds of the required degree.

This does **not** prove `JC(2)`.  It does not exclude:

- an `A2` open subset of `T` not automorphic to a standard `U_omega`;
- an etale map of `B` outside the equivariant/fibration-preserving class;
- the `alpha=1` branch, which does not have this factorization; or
- a new Darboux pair on `B` unrelated to the cited self-maps.

## 6. Primary sources and exact verification

The imported formulas are from A. Dubouloz and K. Palka,
[*The Jacobian Conjecture fails for pseudo-planes*](https://arxiv.org/abs/1701.01425),
Theorem D, equations `(4.3)--(4.7)`, and Lemma 5.1, equations `(5.5)--(5.6)`.
The target-automorphism scope uses H. Flenner and M. Zaidenberg,
[*On a result of Miyanishi--Masuda*](https://arxiv.org/abs/math/0506235),
Theorem 1.1 and the following uniqueness remark.

The exact companion is

```text
04-computation/jc2_dubouloz_palka_standard_chart_obstruction_thm3788.py
```

Run

```bash
python3 04-computation/jc2_dubouloz_palka_standard_chart_obstruction_thm3788.py
python3 -O 04-computation/jc2_dubouloz_palka_standard_chart_obstruction_thm3788.py
```

and compare both byte for byte with

```text
05-knowledge/results/jc2_dubouloz_palka_standard_chart_obstruction_thm3788.out
```

The transcript records `CHECKS=4654`.  It checks the exact degree packet for
`0<=N<=512`, a cyclotomic `N=0` member, its separability resultants, and the
component-choice census.  It also verifies the chart inverse `(6a)`, `(27b)`,
and the exact pole transfer.
Its hostile controls show that the `t=1` address is essential at `N=0` and
that weakening the `R2` contacts from triple to simple destroys the all-degree
contradiction.  The computation verifies the arithmetic certificate; the
proof above is uniform in `N`.
