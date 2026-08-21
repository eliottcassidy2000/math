---
id: THM-3604
title: "Compiler multiarm Hermite--Pade pencil-degree floor"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On a nontrivial
  critical-value compiler with h arms, every nonzero member of a hypothetical
  global Darboux pencil has central-chart a-degree at least 2h, and every
  target basis has a member of a-degree at least 3h.  If
  T=h+n(h-1), the corresponding central-chart total degrees are composite and
  at least 2T, while the pencil height is at least 3T.  In particular these
  imply the weaker floors max(2h,6) and max(3h,8).  The result is only a
  necessary degree invoice; no Darboux pair and no counterexample to JC(2) is
  constructed.
source: root / planar-Jacobian arm-chart degree-invoice session, 2026-08-21
audit: >
  PASS.  An independent current-byte hostile audit rederived the target-basis
  completion quantifier, both central-arm degree invoices, the full
  Hermite--Pade leading-coefficient tariff, its ordinary-total-degree
  conversion, the independent composite-degree input, the A13 specialization,
  and every equality hostile.  It also ran an independent 1,200-row global-ring
  cancellation probe.  Normal and optimized runs are byte-identical to the
  stored 1,849-gate transcript; the AST has no assertion gates, and the
  documentation and committed diff checks pass.
depends_on:
  - THM-3581-critical-value-multiarm-keller-compiler-and-A13-carrier
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
  - THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor
related:
  - THM-3544-planar-keller-target-pencil-total-degree-six-floor
script: 04-computation/jc2_compiler_multiarm_pade_pencil_floor_thm3604.py
output: 05-knowledge/results/jc2_compiler_multiarm_pade_pencil_floor_thm3604.out
script_sha256: eaf4e83f12bf7918ac5338ede266e92b24d2e27227205b46fbed9afb931dcb3d
output_sha256: c5e65493ec62eb16ec24544fbac529b038aa2976c6cbe423f36862e3af798d06
hash_basis: raw LF bytes
---

# THM-3604 -- Compiler multiarm Hermite--Pade pencil-degree floor

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
combines three already proved mechanisms.  It constructs no polynomial
Darboux pair and proves no case of `JC(2)` false.

All rings, maps, degrees, and brackets below are over `C`.

## 1. Setup and statement

Take a nontrivial critical-value compiler from THM-3581.  Thus `n>=2`, the
compiler polynomial `S` is nonconstant, and the smooth target is

```text
Y=Spec A,
A=C[b,c,e]/(c^n e-Sigma(b)),
Sigma(b)=bE(b).                                         (1)
```

The roots of `Sigma` are distinct.  Put

```text
h=deg Sigma=1+|V|>=2,                                  (2)
```

Here `V` is THM-3581's deduplicated set of nonzero critical values, so `h` is
exactly the number of target arms; it is not one plus the number of critical
points when values repeat.  Retain the central arm `b=0`.  Its plane chart
from THM-3600 is

```text
U_0~=A2_(a,c),
a=b/c^n,
b=c^n a,
e=aE(c^n a).                                           (3)
```

Suppose, hypothetically, that

```text
F,G in A,                    {F,G}=kappa in C*.         (4)
```

For `(s,t)!=(0,0)`, write

```text
H_(s,t)=sF+tG,
d_a(H)=deg_a H(a,c),
d_tot(H)=deg_(a,c) H(a,c).                              (5)
```

Also put

```text
T_(n,h)=h+n(h-1).                                      (5a)
```

The degrees in `(5)` are degrees of the polynomial restriction to the
central chart.  They are not degrees of a chosen representative in the
quotient presentation `(1)`, and they are not the generic cover degree.

Then every nonzero target-pencil member satisfies

```text
d_a(H_(s,t))>=2h.                                      (6)
```

For every target basis `(H,K)`, meaning an invertible `GL_2(C)` combination
of `(F,G)`, one has

```text
max(d_a(H),d_a(K))>=3h.                                (7)
```

In ordinary central-chart total degree, every nonzero member satisfies

```text
d_tot(H_(s,t)) is composite,
d_tot(H_(s,t))>=2T_(n,h)>=max(2h,6).                   (8)
```

The total-degree pencil height

```text
height_(a,c)(F,G)
 =max_((s,t)!=(0,0)) d_tot(H_(s,t))
 =max(d_tot(F),d_tot(G))                               (9)
```

satisfies

```text
height_(a,c)(F,G)>=3T_(n,h)>=max(3h,8).                (10)
```

Equations `(6)--(10)` are necessary conditions on a hypothetical pair.
They are not an existence assertion at the boundary values.

## 2. Every target line pays two arm degrees

Fix a nonzero pencil member `H=sF+tG`.  Complete `(s,t)` to an invertible
target matrix and let `K` be the other output.  After rescaling `K`, the pair
`(H,K)` again satisfies a constant-one Darboux equation.  Target `GL_2`
changes preserve the nontrivial compiler and its two central source lines.

THM-3589 therefore applies to `(H,K)`.  On the central target arm

```text
L_e={b=c=0}~=A1_e,                                    (11)
```

each coordinate has degree at least two:

```text
deg_e(H|L_e)>=2,               deg_e(K|L_e)>=2.        (12)
```

Moreover, for every such target basis,

```text
max(deg_e(H|L_e),deg_e(K|L_e))>=3.                     (13)
```

The quantifier over all nonzero pencil directions in `(12)` is supplied by
target-basis completion.  Applying THM-3589 only to the originally displayed
coordinates would not prove `(6)` for the whole pencil.

At `c=0`, formula `(3)` gives

```text
e=E(0)a,                         E(0)!=0.               (14)
```

Consequently restriction to `L_e` changes its parameter only by a nonzero
scalar, and

```text
deg_a H(0,a)=deg_e(H|L_e).                             (15)
```

## 3. The multiarm Hermite--Pade multiplier

THM-3600 gives, for every nonzero global function `H in A`,

```text
d_a(H)>=h deg_a H(0,a).                               (16)
```

Combining `(12),(15),(16)` gives `(6)`.  Combining `(13),(15),(16)` gives
`(7)` for every target basis.  The factor `h` is the degree of the full
squarefree arm divisor; it is lost if one tests only the central plane.

The same theorem gives the sharper leading-coefficient invoice

```text
ord_c(lc_a H)>=n(d_a(H)-floor(d_a(H)/h)).              (17)
```

Hence every nonzero pencil member obeys

```text
ord_c(lc_a H)>=2n(h-1),                               (18)
```

and every target basis contains a member `K` with

```text
d_a(K)>=3h,
ord_c(lc_a K)>=3n(h-1).                               (19)
```

The exact tariff `(17)`, rather than only the coarse values `(18),(19)`,
should be used when an actual degree exceeds its floor.

## 4. Leading order becomes total degree; compositeness is independent

Let `d=d_a(H)` and write

```text
lc_a H=c^q u(c),                 u(0)!=0.              (20)
```

The monomial `c^q a^d` occurs in `H` with nonzero coefficient.  Therefore
`(17)` gives the stronger ordinary-degree inequality

```text
d_tot(H)>=d+q
 >=d+n(d-floor(d/h)).                                  (21)
```

The last expression is nondecreasing in `d`.  At `d=2h` it equals
`2T_(n,h)`, and at `d=3h` it equals `3T_(n,h)`.  Equations `(6),(7)` and
`(21)` therefore give the total-degree parts of `(8),(10)`.

This is the useful connection between the Hermite--Pade tariff and ordinary
degree: recording only `d_tot>=d_a` would discard the boundary vanishing that
all other arms force into the leading coefficient.

By THM-3600, the restriction

```text
(F,G):A2_(a,c)->A2                                    (22)
```

is a polynomial Keller map which is not a polynomial automorphism: if it
were an automorphism, its inverse would put the non-global chart coordinate
`a` in `A`.  The hypothesis `h>=2` in `(2)` is load-bearing here.

THM-3550 now applies to `(22)`.  It says that every nonzero target-pencil
member has composite total degree at least six, and that the total-degree
pencil height is at least eight.  For `n,h>=2`, the stronger numerical floors
`2T_(n,h),3T_(n,h)` already dominate `max(2h,6),max(3h,8)`.  The independent
content imported from THM-3550 is that **every** pencil degree is composite.
That arithmetic conclusion does not follow from the Hermite--Pade multiplier.

## 5. The degree-thirteen A13 carrier

For the A13 row of THM-3581,

```text
n=3,                 h=3,
kappa_13=(72/91)^3,
E(b)=b^2+kappa_13^2,
e=a(c^6a^2+kappa_13^2).                               (22a)
```

Thus every nonzero hypothetical pencil member satisfies

```text
d_a(H)>=6,
ord_c(lc_a H)>=12,
d_tot(H)>=18 and is composite.                        (23)
```

Every target basis contains a member `K` satisfying

```text
d_a(K)>=9,
ord_c(lc_a K)>=18,                                    (24)
```

and the total-degree pencil height is at least `27`.  The numbers `6,9` are
central-chart `a`-degrees; `12,18` are lower bounds on the `c`-orders of the
corresponding leading `a`-coefficients; `18,27` are the resulting ordinary
chart total-degree floors.  The weaker requested shadows are total-member
degree at least six and height at least nine.

## 6. Equality controls and hostile boundaries

The Hermite--Pade invoices are sharp as ring-theoretic statements.  From
`(3)`, for every `r>=1`,

```text
d_a(e^r)=rh,
deg_a(e^r|_(c=0))=r,
ord_c(lc_a(e^r))=nr(h-1),
d_tot(e^r)=rT_(n,h).                                  (25)
```

Thus `e^2` and `e^3` attain exactly the two floors and tariffs in
`(6)--(7)` and `(18)--(19)`.

The every-line arm boundary is simultaneously sharp.  The global functions

```text
N_2=e^2-1,                 N_3=e^3-e                   (26)
```

restrict, after the nonzero rescaling `s=e`, to the immersed noninjective
nodal parameterization

```text
s |-> (s^2-1,s(s^2-1)),                               (27)
```

of degrees `(2,3)`.  Every nonzero member of the pencil generated by
`N_2,N_3` has arm degree two or three, and every target basis has a
degree-three member.  In the central chart the corresponding `a`-degrees are
exactly `2h` and `3h`, with leading-`c` orders `2n(h-1)` and `3n(h-1)`.
Their ordinary total degrees are exactly `2T_(n,h)` and `3T_(n,h)`.

This sharp control is deliberately hostile: both functions lie in `C[e]`, so

```text
{N_2,N_3}=0.                                           (28)
```

It is not a Darboux pair.  It shows that neither the nodal arm lemma nor the
multiarm degree invoice alone can be strengthened numerically; the constant
Jacobian equation must supply any further exclusion.

Other load-bearing boundaries are:

1. **Nontrivial compiler.**  If `S` is constant, the noninjective every-line
   input from THM-3589 disappears.
2. **At least two arms.**  For `h=1`, the target is a plane and the
   chart-automorphism contradiction in Section 4 fails.
3. **Squarefree distinct arms.**  The Hermite--Pade multiplier uses the full
   squarefree divisor of degree `h`; repeated roots change the chart atlas.
4. **Global pair.**  A local Keller pair on one arm plane need not survive the
   other singular shears and does not inherit `(16)`.
5. **Degree types.**  The `a`-degree, ordinary chart total degree, leading
   coefficient `c`-order, and generic cover degree are different invariants.
6. **No Jacobian-conjecture conclusion.**  The theorem assumes a global
   Darboux pair in order to invoice it.  It neither constructs that pair nor
   proves that the remaining degree region is nonempty.

## 7. Preservation and loss ledger

```text
source       hypothetical global Darboux pencil on the compiler target
target       its restriction to the central arm plane A2_(a,c)
map          target GL2 completion, arm restriction, then all-chart descent
preserved    constant Jacobian, each target-pencil direction, arm degree
destroyed    other-arm regularity under central-chart restriction alone
sidecar      the h-arm Hermite--Pade intersection and nonautomorphism gate
cheap test   e^r tariff equality and the nodal (2,3) hostile pair
```

## 8. Exact verification contract

The companion uses exact symbolic arithmetic and explicit failure gates to
check:

- the integer degree and leading-coefficient invoices for finite ranges of
  `n,h,d`;
- equality in `(25)` for several split squarefree targets and powers `e^r`;
- the immersed noninjective nodal control, finite target-`GL_2` hostiles, and
  the zero-Jacobian warning in `(28)`;
- the A13 chart degrees `6,9`, leading-`c` orders `12,18`, and resulting total
  degrees `18,27`; and
- the arithmetic separation between the multiarm floors and the independent
  composite/six/eight planar floors.

The universal every-line theorem, Hermite--Pade multiplier, nonautomorphism
gate, and prime/composite theorem are inherited proof inputs.  The finite
rows are controls, not extrapolations.

Reproduce with

```bash
python3 04-computation/jc2_compiler_multiarm_pade_pencil_floor_thm3604.py
python3 -O 04-computation/jc2_compiler_multiarm_pade_pencil_floor_thm3604.py
```
