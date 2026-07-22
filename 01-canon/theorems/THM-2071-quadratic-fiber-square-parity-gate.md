---
id: THM-2071
title: "Quadratic-fiber square/parity gate for planar Keller pencils"
status: >
  PROVED. If one member P of a complex planar Keller pair has degree two
  along a linear source fiber, then its leading fiber coefficient is a square
  polynomial. After minimizing the complementary component modulo target
  shears Q -> Q-H(P), its fiber degree is odd. If the quadratic leading
  coefficient is constant, the reduced degree is one and the map has an
  explicit tame normal form. Consequently a hypothetical JC(2)
  counterexample with a quadratic pencil member must have a nonconstant
  square leading coefficient and odd reduced complementary fiber degree at
  least three. This is a pencil-direction gate, not a reduction of JC(2), a
  statement about generic cover degree, or a Jelonek classification.
source: codex-2026-07-21-JC2-quadratic-fiber
depends_on: []
related:
  - THM-2063
  - THM-1345
  - THM-1330
  - THM-2045
  - MISTAKE-229
script: 04-computation/jc2_quadratic_fiber_square_gate_codex_20260721.py
output: 05-knowledge/results/jc2_quadratic_fiber_square_gate_codex_20260721.out
script_sha256: 147ade74f38ed8285cbfd5ab9960f6622626cb36db8df5eadf743beb716984fa
output_sha256: 1bf684ba0a0588702df2bae7c5709586a6c15da07c32e579833a3d5b84ff2d2d
hash_basis: repository blobs with LF line endings
---

# THM-2071 -- quadratic-fiber square/parity gate

Let

```text
P=A(x)y^2+B(x)y+C(x),       A!=0,
Q in C[x,y],                J(P,Q)=P_x Q_y-P_y Q_x=kappa in C*.
```

For this fixed source-fiber direction define the reduced complementary degree

```text
mu_y(P,Q) = min_H deg_y(Q-H(P)),       H in C[T].        (1)
```

Every polynomial in the minimum is nonzero, because `Q=H(P)` would have zero
Jacobian with `P`. Thus (1) is an ordinary minimum of a nonempty subset of
the natural numbers.

The conclusions are:

1. `mu_y(P,Q)` is odd.
2. `A=a U^2` for some `a in C*` and `U in C[x]`. Over `C` the scalar may be
   absorbed into `U`, so `A` itself is a square.
3. If `mu_y(P,Q)=1`, the map is tame and has the explicit normal form in
   Section 3 below.
4. More strongly, if `A` is constant, then `mu_y(P,Q)=1`; hence every
   constant-leading quadratic-fiber Keller pair is tame.

In particular, a nonsquare `A` admits no polynomial Jacobian mate at all. A
simple root, or any root of odd multiplicity, in the top fiber coefficient is
already an obstruction.

## 1. Top-coefficient law

Choose `H` attaining (1), put `Q_0=Q-H(P)`, and write

```text
Q_0=sum_(j=0)^n q_j(x)y^j,       q_n!=0,       n=mu_y(P,Q).
```

Target shear does not change the Jacobian, so `J(P,Q_0)=kappa`. For `n>=1`,
the coefficient of `y^(n+1)` is

```text
n A'(x)q_n(x)-2A(x)q_n'(x)=0.                         (2)
```

Equivalently,

```text
(q_n^2/A^n)'=0,
```

and characteristic zero gives

```text
q_n^2=c A^n,                     c in C*.             (3)
```

The case `n=0` is impossible: then `Q_0=S(x)` and

```text
J(P,S)=-(2Ay+B)S'.                                    (4)
```

The coefficient of `y` in (4) forces `S'=0`, contradicting `kappa!=0`.

## 2. Even descent, odd stopping degree, and the square gate

Suppose `n=2k` is positive and even. From (3),

```text
q_n=d A^k                                  (d in C*),  (5)
```

because `C(x)` is a field and `C` contains a square root of `c`. But `P^k`
has top fiber term `A^k y^(2k)`. Therefore

```text
Q_0-dP^k
```

has smaller `y`-degree and the same Jacobian with `P`, contradicting the
minimality of `n`. Hence `n` is odd.

Now factor (3) in the UFD `C[x]`. For every irreducible `pi`,

```text
2 ord_pi(q_n)=n ord_pi(A).
```

Since `n` is odd, every `ord_pi(A)` is even. Thus `A=aU^2`. This proves the
square gate.

The same argument is an explicit terminating descent, not just a minimality
argument: whenever the current complementary degree is even, (5) identifies
the unique top target shear `dP^k`; subtract it and continue. The descent
cannot stop at degree zero, so it stops at an odd degree. This is one rigorous
quadratic-face instance of the leading-form propagation sought in THM-1345;
it says nothing about arbitrary Newton faces.

## 3. Reduced degree one: exact tame normal form

Suppose `n=1`, so after a target shear

```text
Q_0=L(x)y+M(x).
```

The coefficients of `y^2`, `y`, and `1` in `J(P,Q_0)=kappa` are respectively

```text
A'L-2AL'                         =0,                   (6)
B'L-2AM'-BL'                    =0,                   (7)
C'L-BM'                         =kappa.               (8)
```

Here `L!=0`, since `L=0` makes (7)--(8) incompatible with `kappa!=0`. Equation
(6) gives `L^2/A` constant, so write

```text
A=aL^2,                         a in C*.              (9)
```

After division by `L^2`, equation (7) becomes

```text
(B/L)'=2aM'.                                           (10)
```

Thus `B/L-2aM` is constant in `C(x)`. In particular

```text
V=B/L in C[x],                 M=V/(2a)+m_0.          (11)
```

Substituting (9)--(11) into (8) yields

```text
kappa=L (C-V^2/(4a))'.                                (12)
```

Both factors on the right are polynomials and their product is a nonzero
constant. Therefore

```text
L=ell in C*,       A=a ell^2 in C*,
C-V^2/(4a)=(kappa/ell)x+d_0.                          (13)
```

Set

```text
Y=y+B(x)/(2A).
```

Equations (11)--(13) give

```text
P=A Y^2+(kappa/ell)x+d_0,
Q-H(P)=ell Y+m_0.                                     (14)
```

For target coordinates `(u,v)=(P,Q)`, the inverse is

```text
Y=(v-H(u)-m_0)/ell,
x=(ell/kappa)(u-A Y^2-d_0),
y=Y-B(x)/(2A).                                        (15)
```

It is polynomial. The construction consists of a triangular target shear, a
triangular source shear, and affine/elementary triangular maps, so the map is
tame. This degree-one endpoint can also be fed to THM-2063 after the nonlinear
target shear; (6)--(15) record the stronger normal form native to the
quadratic descent.

## 4. Constant quadratic coefficient forces the tame endpoint

Assume now only that `A in C*`; the reduced odd degree `n` may initially be
arbitrary. Completing the square by the polynomial source shear

```text
Y=y+B/(2A),             P=A Y^2+D(x),
D=C-B^2/(4A),
```

preserves the reduced fiber degree. Write the reduced complement as

```text
Q_0=sum_j q_j(x)Y^j.
```

The Jacobian equation is

```text
D' (Q_0)_Y-2AY (Q_0)_x=kappa.                         (16)
```

The top coefficient in (16) makes `q_n` constant. For the odd coefficient
chain, comparison of the coefficient of `Y^(j-1)` gives

```text
q_(j-2)'=(j/(2A))D' q_j,          j=n,n-2,...,3.      (17)
```

Starting from nonzero constant `q_n`, recurrence (17) shows inductively that
`q_1=f(D)` for a polynomial `f` of degree `(n-1)/2` with nonzero leading
coefficient. The constant term of (16) is

```text
D' q_1=kappa.                                         (18)
```

Let `G'=f`. Then (18) gives

```text
G(D(x))=kappa x+constant.                             (19)
```

If `n>=3`, then `deg G=(n+1)/2>=2`; equation (19) is impossible by degree
multiplicativity. (`D` cannot be constant, since then (18) reads `0=kappa`.)
Hence `n=1`, and Section 3 proves tameness.

## 5. Pencil and coordinate invariance

The gate is intrinsic to a pair consisting of an output direction and a
linear source foliation.

Fix `P` and replace a linear target complement by

```text
Q'=alpha Q+beta P,                  alpha!=0.
```

As `H` ranges over `C[T]`, so does

```text
K(T)=(H(T)-beta T)/alpha,
```

and

```text
Q'-H(P)=alpha (Q-K(P)).
```

Thus `mu_y(P,Q')=mu_y(P,Q)`. Scaling `P` also changes neither the algebra
`C[P]` nor the square class of `A`, because every nonzero complex scalar is a
square. Affine changes preserving the chosen source foliation send

```text
x -> ax+b,                 y -> cy+linear(x),
```

so they preserve fiber degree, reduced degree, and the square property of the
top coefficient. For a different source direction, one simply applies the
theorem in a new affine coordinate system.

Consequently, for an arbitrary planar Keller map and any nonzero output-pencil
member `R`, any linear complement defines the same reduced degree. Combining
this theorem with THM-2063 gives the all-directions gate for a hypothetical
counterexample:

```text
deg_fiber R <= 1    is impossible by THM-2063;

deg_fiber R = 2     forces
                    (top quadratic coefficient) = nonconstant square,
                    mu_fiber(R, complement) in {3,5,7,...}.             (20)
```

## 6. Newton/Jelonek frontier and the exact remaining lemma

In the Newton picture, `A(x)` is the coefficient polynomial on the top
`y=2` row. The theorem says that its zero divisor must be even. Hence the
generic squarefree top row, and in particular every simple top-row zero, is
excluded. The constant top row is not a residual: Section 4 closes it. The
first unresolved cell is therefore

```text
A=U(x)^2 with U nonconstant,       mu_y=3.             (21)
```

The next coefficient equation already centers this cell. After absorbing a
scalar so `A=U^2`, writing the reduced cubic top as `q_3=cU^3`, and using a
remaining linear target shear to normalize the quadratic coefficient, one
gets

```text
q_2=(3c/2)B U.                                         (22)
```

Indeed, the coefficient of `y^3` is

```text
2A'q_2-2Aq_2'+3B'q_3-Bq_3'=0.
```

After substituting `A=U^2` and `q_3=cU^3`, this becomes

```text
(q_2/U^2)'=(3c/2)(B/U)',
```

so `q_2=eU^2+(3c/2)BU`; subtracting `eP` gives (22).

These are exactly the first two coefficients of

```text
c (U y+B/(2U))^3,
```

but the centered expression need not be polynomial when `U` does not divide
`B`. A precise missing lemma is therefore:

> decide whether the lower two coefficient equations in the cubic residual
> (21)--(22) force a polynomial centering/divisibility, force a contradiction,
> or admit a genuinely new Keller-pencil normal form.

This is a rank-one-at-infinity coefficient problem, not VC(4). The integer
`mu_y=3` is also not the generic function-field cover degree `3`: no map
between those two degree strata is proved. THM-1330's Jelonek/monodromy
constraints remain separate necessary data, and MISTAKE-229 forbids calling
this descent equivalent to the Jelonek, Vanishing-Conjecture, or
Lame-for-polygons programs.

## 7. Exact symbolic referee

The companion script verifies:

- (2) for complementary degrees `1` through `7`;
- the even target-shear cancellation for powers `P^1` through `P^4`;
- the three equations (6)--(8);
- the odd-chain primitive degrees in Section 4 for `n=3,5,7`;
- a nontrivial rational instance of (14), its constant Jacobian, and both
  compositions of (15);
- bounded hostile linear-system searches against a nonsquare top coefficient
  and a square-but-nonconstant control.

The bounded searches are controls only; the proof above is degree-uniform.
The stored output ends in `RESULT: PASS`. QED.
