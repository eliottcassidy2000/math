---
id: THM-2102
title: "Power-free weighted faces are tame and proper-power faces obey a first-defect descent"
status: >
  PROVED. If a planar Keller component has a top form that is not a proper
  power for any positive integer weight, a weighted centralizer descent
  forces its full polynomial to be triangular and hence a coordinate.
  Therefore every hard component has a proper-power top face for every
  positive weight. For a common top pair A H^m,B H^n, the first lower layer
  obeys an exact Hamiltonian first-defect equation. A terminal defect makes H
  weighted-linear; an earlier nonresonant defect either lifts to a common
  approximate root or leaves an explicit divisibility quotient obstruction,
  while resonant layers have an exact target-shear criterion. This corrects
  HYP-8955/S110's linear-face and DvdK mechanism claims; it does not prove
  that every quotient obstruction vanishes or settle JC(2).
source: codex-2026-07-22-JC2-power-free-face-descent
depends_on: []
related:
  - THM-2045
  - THM-2063
  - HYP-8950
  - HYP-8955
  - THM-2113
script: 04-computation/jc2_power_free_face_descent_codex_20260722.py
output: 05-knowledge/results/jc2_power_free_face_descent_codex_20260722.out
script_sha256: e836e3c589ef388b5ab8d01e27ba18cb39337972aa3ebc473214e1d337066163
output_sha256: 26d29a591c1454b0e222d78420fc046649507aea538d29161938db518a8e4c87
hash_basis: repository blobs with LF line endings
---

# THM-2102 -- power-free faces and the first proper-power defect

Work over `C[x,y]`. For positive integer weights `w=(a,b)`, write

```text
deg_w(x^i y^j)=ai+bj,          W=a+b.                  (1)
```

The Jacobian bracket is `J(f,g)=f_x g_y-f_y g_x`.

## 1. Weighted common-power lemma

Let `F,G` be nonconstant weighted-homogeneous polynomials of degrees `D,E`.
If `J(F,G)=0`, then

```text
F=A H^m,        G=B H^n                               (2)
```

for a weighted-homogeneous `H` that is not a proper power and nonzero
constants `A,B`.

Indeed, contract `dF wedge dG=0` with the weighted Euler vector field
`a x partial_x+b y partial_y`. Euler's identity gives

```text
D F dG=E G dF.                                        (3)
```

Thus the rational function `G^D/F^E` has both partial derivatives zero and
is constant. Unique factorization applied to

```text
G^D=c F^E                                               (4)
```

gives (2), after taking the gcd of all irreducible-factor multiplicities.
If `F` itself is not a proper power, then `m=1`, `D|E`, and

```text
G=c' F^(E/D).                                          (5)
```

This proof uses only weighted Euler and UFD arithmetic.

## 2. Power-free top-face theorem

> **Theorem.** Suppose `J(f,g)=c in C*`. If, for some positive weight, the
> top weighted form `F=in_w(f)` is not a proper power up to scalar, then `f`
> is a coordinate. More precisely, after possibly swapping `x,y`,
>
> ```text
> f=A(x)+lambda y,          lambda!=0.                 (6)
> ```

Let `D=deg_w f`. Subtract the constant term of `g`. If its top form `G` has
degree `E` with `D+E>W`, then the top bracket must vanish. Equations (2)--(5)
give `G=c'F^k`; replace `g` by `g-c'f^k`. This leaves `J(f,g)=c` unchanged
and strictly lowers `E`. Repeat.

The process cannot end below `D+E=W`: if `D+E<W`, every possible bracket
monomial has negative weight, so the bracket is zero. It also cannot retain a
commuting top form at equality, because another subtraction would put it
below the equality line. Hence it terminates with

```text
D+E=W,          J(F,G)=c,          E>0.                (7)
```

In particular `D<W`. Every mixed monomial `x^i y^j`, `i,j>=1`, has weight
at least `W`, so

```text
f=A(x)+B(y).                                           (8)
```

If `p=deg A` and `q=deg B` are both at least two, then `pa<=D<W` and
`qb<=D<W` imply

```text
(p-1)a<b,             (q-1)b<a,                       (9)
```

whose product is impossible. Thus one side is linear. If one side of (8) is
constant and the other has degree at least two, its top monomial is a proper
power, contrary to the hypothesis. This proves (6), whose triangular inverse
is explicit. QED.

Consequently, a hypothetical noncoordinate Keller component has a proper-
power top form for **every** positive weight. Square-free, irreducible, and
more generally power-free multi-root faces are all excluded.

## 3. The exact first-defect equation

The residual begins with a common power. Suppose

```text
F=f_D=A H^m,          G=g_E=B H^n,
deg_w H=kappa,
Delta=D+E-W=(m+n)kappa-W>0,                            (10)
```

where `H` is power-free. Let `r>0` be the first weight drop at which either
polynomial has a nonzero component, and write

```text
P=f_(D-r),        Q=g_(E-r),                           (11)
```

using zero when one component is absent. Then `r<=Delta`. The bracket part of
weight `Delta-r` is exactly

```text
J(H,L),
L=A m H^(m-1)Q-B n H^(n-1)P.                          (12)
```

To see this, the only contributions at the first drop are

```text
J(F,Q)+J(P,G)
=J(H,A m H^(m-1)Q-B n H^(n-1)P).                     (13)
```

If `r>Delta`, no lower-layer pair can supply weight zero, contradicting the
constant bracket. Therefore

```text
r<Delta  implies J(H,L)=0,
r=Delta  implies J(H,L)=c.                            (14)
```

The degree of `L` is `(m+n-1)kappa-r`. In the first line of (14), the
common-power lemma gives

```text
kappa not|r  implies L=0,
r=t kappa    implies L=lambda H^(m+n-1-t).             (15)
```

In the second line, `deg_w L=W-kappa` and `L` is a mate for `H`. Applying
Section 2 to the quasi-homogeneous polynomial `H` shows that `H` is a
weighted-linear coordinate. Thus the terminal first defect exposes the tame
primitive root directly.

## 4. Approximate-root and target-shear sidecars

Assume `r<Delta` and `kappa not|r`, so `L=0`. A simultaneous first-order
approximate root

```text
K=H+R,           deg_w R=kappa-r                      (16)
```

absorbs both layers in (11) if and only if

```text
H^(m-1)|P       and       H^(n-1)|Q.                  (17)
```

When it exists it is unique and

```text
R=P/(A m H^(m-1))=Q/(B n H^(n-1)).                   (18)
```

This follows by expanding `A(H+R)^m` and `B(H+R)^n` to the first dropped
weight. Equation `L=0` makes the two quotients agree. The divisibility is
automatic when `min(m,n)=1`; for `m,n>=2` its failure is a genuine quotient
class, such as `[P] mod H^(m-1)`.

At a resonant drop `r=t kappa<Delta`, (15) leaves a central scalar class. A
target shear `g -> g-lambda f^k` can change that layer exactly when

```text
n-t=k m             for some integer k>=0.             (19)
```

Indeed, the shear's top weight is `kD=km kappa`, whereas the layer to be
changed has weight `E-r=(n-t)kappa`. If (19) fails, no polynomial in `f`
reaches the layer.

The smallest synchronization-without-lift control is

```text
w=(1,2),       H=y+x^2,       (m,n)=(2,3),       r=1,
P=x^3,         Q=(3/2)H x^3.                                (20)
```

Here `L=0`, but `H` does not divide `P`, so no common approximate root exists;
no target shear has the required weight. In fact

```text
J(H^2+P,H^3+Q)=(9/2)x^5.                               (21)
```

Thus first-layer synchronization is necessary but not sufficient. The next
honest JC(2) target is to show that later constant-bracket equations kill
these quotient classes, or to classify a surviving class.

## 5. Repair of the single-face program

For a quasi-homogeneous `f` of degree `delta`, decompose a mate by weight.
The constant bracket directly forces

```text
delta<W.                                               (22)
```

Then no mixed monomial occurs, and the axis-degree argument in (9) forces a
linear axis term. This proves HYP-8955/S110's base theorem, but not by the
mechanism recorded there.

Three corrections are essential:

1. `xy` is primitive; a monomial factor is not polynomial compositeness.
2. The face polynomial for `x^3+y^2` at weights `(2,3)` is linear, but
   `delta=6>5` and no mate exists. “Linear face” is not equivalent to a
   weighted-linear coordinate; the genus/degree deficit (22) is independent.
3. On the torus, the face bracket is a Laurent Wronskian/Bézout equation, not
   a DvdK constant-term moment. If `D+E=W`, all its monomials have weight zero
   and

   ```text
   J(F,G)=[x]F [y]G-[y]F [x]G.                        (23)
   ```

   Hence the terminal face map has rank one exactly when `F` has a linear
   `x` or `y` term. No DvdK transfer is needed or presently justified.

The lawful program is therefore

```text
power-free face -> triangular coordinate;
proper-power face -> first-defect Hamiltonian equation;
earlier defect -> approximate-root divisibility or resonant shear class.
```

Repeated proper-power faces really occur for coordinates—for example

```text
f=x+(y+x^2)^m,       g=y+x^2,       J(f,g)=1.           (24)
```

so they cannot be discarded by a multi-root slogan. This theorem sharpens
the residual but does not make the quotient obstruction vanish and does not
settle JC(2).

## Exact controls

The companion verifies (12), (18), (21), (23), and the coordinate controls
with exact symbolic arithmetic. These checks illustrate the identities; the
proof above is algebraic and does not depend on computation.
