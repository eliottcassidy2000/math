---
id: THM-3790
title: "Cubic pseudo-plane arm nodal-immersion gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every
  Darboux morphism on the THM-3785 surface restricts on its triple arm to
  a noninjective polynomial immersion.  Its first normal jets satisfy an
  exact Bezout determinant law.  Consequently neither boundary coordinate
  is constant or linear, the two degrees cannot both be quadratic, and the
  first possible bidegree is (2,3).  The nodal cubic profile realizes that
  boundary gate, but its canonical first-order lift has seven critical
  points and is not a Darboux pair.  Arbitrary Darboux pairs remain open.
source: root / multiple-fibre arm and injective-line lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS by jc_zero_debt_lift, 2026-08-23.  The audit
  rederived the conormal quotient and Poisson reduction, checked
  Gwozdziewicz's theorem against the primary source with all hypotheses
  typed, reconstructed the quadratic-centre obstruction, and proved the
  seven-point critical family exhaustive by equality with a triangular
  Groebner ideal.  It also verified the universal nodal normal coefficient
  and residual Picard-class argument.  The companion was hardened from
  Python assertions to 32 active gates plus an AST assertion audit.  Normal
  and optimized runs byte-match the frozen transcript.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
related:
  - THM-3788-dubouloz-palka-standard-chart-containment-obstruction
external:
  - "Gwozdziewicz, Injectivity on one line, arXiv:alg-geom/9305008, Theorem 1.1."
script: 04-computation/jc2_cubic_pseudoplane_arm_nodal_gate_thm3790.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_arm_nodal_gate_thm3790.out
script_sha256: d6cbd4b83fbc46a8914bbc63cadf9f679b87b74566ea185909df900ed5e2573d
output_sha256: abe3b659d56097d402c27ff169f84d53afa0f2c57efc465617243a08570e0434
hash_basis: raw LF bytes
---

# THM-3790 -- every live arm map is a nodal-type immersion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over
`C`, fix `c!=0`, and retain the THM-3785 smooth symplectic surface

```text
Y=Spec B,                  B=C[r,z,e]/(r^2e-z^3+c^3r),                (1)
```

with

```text
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3c^3+6re.                  (2)
```

Its triple arm is

```text
L=V(r,z)=Spec C[e].                                                   (3)
```

Suppose `A,C in B` satisfy `{A,C}=1`.  Write their restrictions as

```text
gamma(e)=(a_0(e),b_0(e))=(A|L,C|L): A1 -> A2.                        (4)
```

Then `gamma` is an immersion but is not injective.  In particular it must
parametrize a singular rational plane curve by identifying at least two
distinct arm points; a cusp is forbidden because the derivative never
vanishes.

## 1. The exact first-normal-jet law

Let `I=(r,z)` be the ideal of `L`.  The surface relation gives

```text
c^3r=z^3-r^2e in I^2,                                                   (5)
```

so

```text
B/I^2 = C[e,z]/(z^2).                                                    (6)
```

There are unique polynomials `a_1,b_1` with

```text
A=a_0(e)+z a_1(e) mod I^2,
C=b_0(e)+z b_1(e) mod I^2.                                              (7)
```

Modulo `I`, only `{z,e}=3c^3` survives from `(2)`.  Taking the bracket of
`(7)` and using `{A,C}=1` gives the polynomial identity

```text
3c^3[a_1 b_0'-a_0' b_1]=1.                                             (8)
```

Thus `a_0'` and `b_0'` cannot vanish at the same point.  This proves that
`gamma` is a polynomial immersion and records the exact normal sidecar that
an arm parametrization must carry.  Merely guessing a singular rational
curve is insufficient; its tangent row must admit the unimodular completion
`(8)`.

## 2. Injective arm images would prove too much

Let

```text
phi:A2_(x,y) -> Y
```

be the cubic affine-plane atlas of THM-3785.  On the source line `x=0`,

```text
phi(0,y)=(0,0,3c^2y),                                                    (9)
```

so `phi` restricts to an isomorphism from that line onto `L`.  The pullback

```text
(A o phi,C o phi):A2 -> A2                                               (10)
```

is a polynomial Keller map.  If `gamma` were injective, `(10)` would be
injective on the line `x=0`.  Gwozdziewicz's injective-line theorem then
forces `(10)` to be a polynomial automorphism.  This is impossible: its
function-field degree is the cubic atlas degree times the positive degree of
`Y -> A2_(A,C)`, hence is divisible by three.  Therefore

```text
gamma is noninjective.                                                    (11)
```

This argument uses the actual consequence of line injectivity, not a
heuristic about the nonproperness curve.

## 3. The first possible boundary bidegree

Neither `a_0` nor `b_0` can be constant.  If, say, `a_0` were constant,
immersion would make `b_0'` nowhere zero; over `C` this forces `b_0` to be
linear, making `gamma` injective.  A linear nonconstant coordinate also
makes `gamma` injective.  Hence

```text
deg a_0>=2,                         deg b_0>=2.                          (12)
```

They cannot both be quadratic.  Complete their squares.  A collision
`a_0(s)=a_0(t)` with `s!=t` forces `s+t` to be twice the centre of `a_0`;
the same collision for `b_0` forces the two quadratic centres to coincide.
Both derivatives then vanish at that centre, contradicting immersion.
Consequently

```text
max(deg a_0,deg b_0)>=3,             deg a_0+deg b_0>=5.                (13)
```

The first unordered bidegree is `(2,3)`.  It is genuinely present at the
arm-jet level:

```text
gamma_*(e)=(e^2,e^3-e).                                                   (14)
```

Its derivative `(2e,3e^2-1)` never vanishes simultaneously, while
`gamma_*(1)=gamma_*(-1)=(1,0)`.  Thus the first live boundary geometry is a
node, not another embedded line or a cusp.

## 4. The canonical nodal lift fails globally

For `(14)`, one Bezout completion in `(8)` is

```text
a_1=-1/(3c^3),                    b_1=-e/(2c^3).                         (15)
```

The resulting smallest global lifts are

```text
A_*=e^2-z/(3c^3),
C_*=e^3-e-ez/(2c^3).                                                     (16)
```

Their bracket factors exactly as

```text
{A_*,C_*}=(c^3+2er)(2c^3+z)/(2c^6),                                    (17)
```

so the correct arm jet does not propagate to a Darboux pair.  More
decisively, `A_*` has exactly seven displayed critical points.  For every
root

```text
8zeta^7+9c^15=0,                                                         (18)
```

put

```text
r=2zeta^3/c^3,                  e=-c^6/(4zeta^3).                        (19)
```

These points lie on `(1)`, satisfy `c^3+2re=0`, and make all three brackets
of `A_*` with `r,z,e` vanish.  The seven roots in `(18)` are distinct and
nonzero.  Conversely, the critical equations and the surface equation have
the exact triangular Groebner basis, over `C(c)`,

```text
9c^9e-2z^4,                 c^3r-2z^3,
8z^7+9c^15.                                                        (19a)
```

Thus there are no other critical points: `(18)--(19)` are exhaustive, not
merely a displayed family.  Hence no polynomial correction of only the
second output can repair `(16)`; a live construction must change the carrier
itself and must retain the nodal arm data plus the normal Bezout sidecar.

## 5. Every node forces a residual divisor

There is a further global sidecar even before choosing higher jets.  For a
hypothetical Darboux pair whose arm restriction is the normalized nodal
profile `(14)`, let

```text
Delta(A,C)=C^2-A(A-1)^2.                                                  (20)
```

This is the irreducible equation of the nodal cubic.  Substituting `(7)`
and using `(8)` gives, independently of the chosen Bezout completion,

```text
Delta(A,C)=-z(e^2-1)/(3c^3) mod I^2.                                    (21)
```

Indeed, the coefficient of `z` is

```text
(e^2-1)[2e b_1-(3e^2-1)a_1]=-(e^2-1)/(3c^3).                            (22)
```

Thus the pullback divisor of the nodal cubic contains `L` reduced at its
generic point, but at each of the two preimages `e=+/-1` of the node it has
another local branch.  This is also forced directly by etaleness: locally
the Darboux morphism is an isomorphism and the inverse image of a node has
two transverse branches.  Globally write

```text
div(Delta(A,C))=L+D.                                                       (23)
```

The residual effective divisor `D` is nonempty, meets `L` at both nodal
preimages (possibly through different global components), and satisfies

```text
[D]=-[L]=2[L] in Pic(Y)=Z/3.                                              (24)
```

So a construction cannot specify only a noninjective arm parametrization:
it must simultaneously build the residual class-two curve that realizes
the other local branch at every collision.

The exact companion named in the metadata verifies `(6),(8),(14)--(22)`, the
complete critical ideal over `Q(c)`, and the assertion-free optimized replay.
**QED.**
