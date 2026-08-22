# The coupled critical divisor has an exact inverse graph, but not a Keller cofactor

**Status:** NON-CANONICAL RESEARCH REFLECTION.  The universal identities below
are proved symbolically, and the two named accessory controls are
**FINITE-EXACT**.  While this probe ran,
[THM-3289 -- affine transverse `C_0,E_0` coupled clutch](../01-canon/theorems/THM-3289-affine-transverse-c0-e0-coupled-clutch-critical-no-go.md)
was independently promoted to PROVED.  This note reconstructs its two named
controls but does not claim an independent audit of its universal parameter
argument or reserve a theorem ID.

The matching exact artifact is
[`jc_critical_inverse_cover_cofactor_jacobian_probe_agent.py`](../04-computation/jc_critical_inverse_cover_cofactor_jacobian_probe_agent.py),
with frozen
[`output`](../05-knowledge/results/jc_critical_inverse_cover_cofactor_jacobian_probe_agent.out).

## Inheritance pass and portfolio

- **Anchor:** the missing branchwise cofactor/inverse-cover information after
  [THM-3279 -- affine transverse clutch critical no-go](../01-canon/theorems/THM-3279-affine-transverse-clutch-critical-no-go.md).
- **Closest proved mechanism:** THM-3279 retains only the saturated scalar
  critical resultant after localizing by `V`; affine `C_0` leaves at least
  `50` units off the owner divisor.
- **Infinity-staircase inheritance:**
  [THM-3237](../01-canon/theorems/THM-3237-degree-nine-jacobian-infinity-wall-and-square-root-escape.md),
  [THM-3257](../01-canon/theorems/THM-3257-degree-eight-tuned-cubic-infinity-wall-and-three-root-critical-escape.md),
  [THM-3263](../01-canon/theorems/THM-3263-degree-seven-retuned-quartic-infinity-wall-and-odd-critical-resultant-inertia.md),
  [THM-3265](../01-canon/theorems/THM-3265-degree-six-retuned-quintic-infinity-wall-and-five-root-critical-escape.md),
  and
  [THM-3276](../01-canon/theorems/THM-3276-degree-at-most-eight-polynomial-clutch-critical-no-go.md)
  divide boundary multiplicity and follow the first nonzero reciprocal row;
  they explicitly discard root labels and Keller cofactors.
- **Canonical hostile:**
  [THM-3066](../01-canon/theorems/THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor.md)
  proves that a symmetric product of branch cofactors is blind to the
  sheetwise Keller equation.
- **Least-used sidecar:**
  [THM-3064](../01-canon/theorems/THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary.md)
  retains a fixed inverse sheet and its primitive-element cofactor.  That is
  the true Keller sidecar; it must not be confused with the elimination
  cofactor constructed below.
- **Niche:** recover the critical root before asking for a polynomial mate.
- **Wildcard:** separate three conditions that the scalar resultant blends:
  reducedness, existence of a selected-root chart, and a unit cofactor pair.

The live concept board was:

| object | predicate | operation | lost coordinate / cheapest test |
|---|---|---|---|
| saturated scalar `H(x)` | critical divisor nonempty/reduced | boundary quotient | loses `y`; test `gcd(H,H')` |
| quadratic row `h` | cubic cancellation | `4R_2+V'R_1` | leading pivot `a`; test `gcd(H,a)` |
| linear row `ell` | selected critical root | fraction-free pseudo-division | denominator `ell_1`; test `gcd(H,ell_1)` |
| pair `(U,W)` | elimination cofactors generate the unit ideal | retain the Bezout row | not a Keller cofactor; test `U-(V'/4)W=a^2` |
| primitive Keller cofactor | sheetwise constant Jacobian | inverse normalization | absent; impose `Jac(P,Q)=1` before scalarization |

The method cards used were **refine and saturate before transporting a
resultant shadow**, **controlled forgetting requires a sidecar**, and **divide
exceptional multiplicity before judging a wall**.  No new meta-pattern is
proposed from one lane.

## Universal cubic-to-linear sidecar

Retain the coupled affine gradient pair from THM-3289,
but treat it here as an algebraic object rather than a proved universal
no-go.  On `V!=0`, with `y=Vz`, write

```text
L=y^2+y+CV,

R_1=2L(2y+1)+VA,
R_2=V^3 k+V^2y+L(-V'y+2V^2d),                         (1)
```

where `d=C'` and `k=E'` are constants.  Expand

```text
R_1=4y^3+6y^2+q y+r,
q=2+4CV,                 r=V(2C+A).                   (2)
```

The constant leading coefficient `4` permits a denominator-free cubic
cancellation:

```text
h=4R_2+V'R_1=a y^2+b y+c,                              (3)

a=2V'+8V^2d,
b=2V'+4V^2(1+2d),
c=V(4V^2k+8CV^2d+V'(2C+A)).                            (4)
```

Put

```text
D=6a-4b,
ell_1=a(aq-4c)-Db,
ell_0=a^2r-Dc,
Q=4ay+D.                                                (5)
```

One fraction-free pseudo-division gives the exact linear row

```text
ell=ell_1 y+ell_0
   =a^2R_1-Qh
   =U R_1+W R_2,                                       (6)

U=a^2-V'Q,                    W=-4Q.                    (7)
```

The cofactor pair has the unexpectedly simple unit identity

```text
U-(V'/4)W=a^2.                                         (8)
```

Thus, wherever `a` is a unit, `(U,W)` is a unimodular **elimination**
cofactor pair.  This is useful branchwise information, but it is not the
primitive-element cofactor `q_X/f_X'` of THM-3064 and does not encode a
second polynomial coordinate.

There is also an independent norm identity.  A literal `6 x 6` Sylvester
determinant gives the standard-oriented resultant, while direct substitution
of the root `-ell_0/ell_1` into `(3)` gives

```text
a ell_0^2-b ell_0 ell_1+c ell_1^2
  =16a^2 Res_y(R_1,R_2).                                (9)
```

Equation `(9)` was expanded and checked identically.  It is an independent
reconstruction of the scalar resultant from the quadratic row and its linear
pseudo-remainder; no black-box resultant participates in `(9)`.

## Exact inverse-graph lemma

Let `H` be a saturated off-owner resultant factor over a field, and suppose
in the quotient `K[x]/(H)` that

```text
Res_y(R_1,R_2)=0,          a in units,          ell_1 in units.           (10)
```

Then set

```text
Y=-ell_0/ell_1 in K[x]/(H).                              (11)
```

Equation `(9)` gives `h(x,Y)=0`.  Equations `(6)` and `(10)` then give
`R_1(x,Y)=0`, and `(3)` gives `R_2(x,Y)=0`.  Conversely, every common zero
of `R_1,R_2` satisfies `ell_1y+ell_0=0`, so `(11)` is its unique `y`
coordinate.  Therefore the critical scheme over `H` is the graph

```text
Spec(K[x]/(H)) -> A^2,             x |-> (x,Y(x)).       (12)
```

If `H` is squarefree, this is a reduced inverse critical cover.  Notice what
the lemma does and does not say:

- it restores the root coordinate erased by scalar elimination;
- `(8)` restores a unit **pair** of gradient-elimination cofactors;
- it does not produce a polynomial inverse of a generically finite map;
- it does not produce the sheetwise primitive Keller cofactor; and
- it does not solve `P_xQ_z-P_zQ_x=1` for a polynomial `Q`.

## Two exact degree-53 controls

The companion independently reconstructs both cubic accessory fields from
the literal THM-3212 formulas and specializes

```text
C=1+x,                       d=1,                       k=1.              (13)
```

It uses a literal Sylvester determinant for the `40`-term universal factor,
divides the exact degree-`44` passport boundary, and obtains the same monic
degree-`53` `H` digests printed by the THM-3289 companion.  It then
computes the new sidecars directly in characteristic zero:

| passport | `deg H` | `deg a` | `deg ell_1` | `deg ell_0` | exact unit tests |
|---|---:|---:|---:|---:|---|
| `(4,1,1,1)` | 53 | 32 | 80 | 88 | `gcd(H,H')=gcd(H,ST)=gcd(H,a)=gcd(H,ell_1)=1` |
| `(3,2,1,1)` | 53 | 32 | 80 | 88 | `gcd(H,H')=gcd(H,ST)=gcd(H,a)=gcd(H,ell_1)=1` |

The common-scale-normalized projective pair digests `(ell_1,ell_0)` are

```text
(4,1,1,1):
5ea209d5455f7fb13488cb11fd1e82aded97ded1f6c3681b27b3b5ebc6904a95,

(3,2,1,1):
7436a32bbdd517c27a9ed7a977784ff64706f7ab5f7e108891d13e68995cbae2.
```

The pair is normalized by one common scalar (the leading coefficient of
`ell_1`), not by making both entries monic separately; the digest therefore
retains the relative scale and freezes the rational section `-ell_0/ell_1`.

Hence each fixed control has exactly the graph description `(12)` and the
unimodular elimination-cofactor pair `(7)--(8)`.  This is stronger data than
the scalar statement “`H` has 53 squarefree roots,” but only for the two
named controls.  It neither proves nor independently audits THM-3289's universal affine
`C_0,E_0` PRS argument.

## Hostiles and failure anatomy

Two small exact examples separate the three gates.

### Nonreduced resultant destroys the linear unit

Take

```text
f=4y^3,                       h=y^2+x.
```

Then

```text
Res_y(f,h)=16x^3,             ell_1=-4x.                (14)
```

The reduced support has one base point, but three coincident branches are
stored in the resultant multiplicity and `ell_1` is not a unit.  Replacing
the resultant by its radical before retaining the linear row loses exactly
the information needed for the inverse chart.

### Squarefree resultant does not force the quadratic pivot to be a unit

Take

```text
f=4(y^3-1),                   h=xy^2+y-1.
```

Then

```text
Res_y(f,h)=16x(x^2+3),        a=x.                      (15)
```

The resultant is squarefree, and root recovery at `x=0` survives through a
linear chart, but the quadratic pivot `a` vanishes there.  Thus squarefree
scalar data alone do not justify the unimodular shortcut `(8)`; `a` is a
separate sidecar.  The actual accessory controls pass both tests.

The first failed implication is therefore

```text
squarefree scalar resultant
  -> one fixed inverse/cofactor chart without further unit tests          FALSE.
```

The repaired statement is the exact two-pivot lemma `(10)--(12)`.

## SymPy sign audit

The installed SymPy `1.14.0` gives

```text
resultant(x-a,x^3+c,x)=-(a^3+c),
resultant(x^3+c,x-a,x)=-(a^3+c).
```

Under the standard convention the first sign is positive and the second is
negative.  The literal Sylvester implementation reproduces those two
standard signs and confirms broken swap antisymmetry in the installed
routine.  No load-bearing calculation in the companion calls
`sympy.resultant`; the call occurs only in this explicit hazard test.

## A canonical divergence class for the missing Keller sidecar

The distinction between an elimination cofactor and a Keller mate can be
made exact before choosing any ansatz.  Let `R=K[x,z]` in characteristic
zero, let

```text
D_P=P_x partial_z-P_z partial_x,
```

and suppose first that the gradient ideal is the unit ideal.  Choose any
Bezout row

```text
A P_x+B P_z=1.                                           (16a)
```

Define

```text
mu(P)=[A_x+B_z] in coker(D_P)=R/D_P(R).                 (16b)
```

This refines the repo's HYP-8950 Hamiltonian-cokernel formulation
`1 in im(D_P)`: after gradient unimodularity supplies any transverse Bezout
row, the remaining obstruction has the canonical lower-degree representative
`mu(P)`.

This class is independent of the chosen row.  Indeed, coprimality of
`P_x,P_z` says every other row is uniquely of the form

```text
(A+hP_z, B-hP_x),
```

and its divergence is `(A_x+B_z)-D_P(h)`.  Moreover,

```text
P has a polynomial mate Q with Jac(P,Q)=1
  iff (P_x,P_z)=R and mu(P)=0.                           (16c)
```

For the forward implication take `(A,B)=(Q_z,-Q_x)`, whose divergence is
zero.  Conversely, if `D_P(h)=A_x+B_z`, the adjusted Bezout row has zero
divergence.  The polynomial one-form

```text
-(B-hP_x) dx+(A+hP_z) dz
```

is therefore closed, hence exact on affine two-space.  Its primitive `Q`
satisfies `Q_x=-(B-hP_x)`, `Q_z=A+hP_z`, and `(16a)` becomes
`Jac(P,Q)=1`.

Thus Keller entry splits into two typed gates:

```text
no critical point / gradient-unimodularity,
then vanishing of the canonical divergence class mu(P). (16d)
```

The second gate is not cosmetic.  For the standard punctured-fibre hostile

```text
P=x+x^2 z,
```

the gradient is unimodular because

```text
(1-2xz)P_x+4z^2P_z=1,
```

but this row has divergence `6z`.  Here

```text
D_P=(1+2xz)partial_z-x^2partial_x,
```

and `6z` is not in its image.  To see this, write
`h=sum_j a_j(x)z^j`.  If its top `z`-degree is `m>1`, the top equation forces
`a_m=c x^(2m)`; subtracting `cP^m` lowers the degree without changing
`D_P(h)`.  Reduction therefore reaches degree at most one.  At degree one,
the constant and linear equations would give

```text
a_1=x^2a_0',               -x^4a_0''=6,
```

which has no polynomial solution.  Hence `mu(P)=[6z]` is nonzero.  This
recovers the global obstruction missed by critical-point tests on a minimal
example and turns the vague request for an "integrable cofactor" into one
explicit Hamiltonian-cokernel class.

The class does not rescue the two degree-53 controls above: those controls
already have critical points, so their gradient ideals are not units.  Its
role is upstream.  For any future deformation that kills the critical
divisor, compute one global gradient Bezout row and then test `(16b)`, rather
than searching blindly for both coefficients of a mate.

## Connection contract and integration recommendation

The new typed connection is

```text
source:      saturated critical resultant H;
target:      reduced critical graph plus elimination-cofactor pair;
map:         cubic cancellation followed by one fraction-free pseudo-division;
preserved:   every off-owner common gradient zero and its x-coordinate;
restored:    the unique y-coordinate and a unimodular gradient cofactor pair;
destroyed:   primitive inverse-sheet label, Keller cofactor, mate Q, global chart;
sidecars:    a and ell_1;
test:        gcd(H,a)=gcd(H,ell_1)=1.                   (17)
```

This changes the search priority.  On a squarefree boundary-disjoint control
that passes `(17)`, selecting and gluing critical roots is no longer the
hard part: the scalar divisor already has a canonical rational graph after
one subresultant row.  Those roots are obstructions, not candidate inverse
sheets.  The unresolved planar-JC information lies earlier and in a different
module: a polynomial mate `Q` and its primitive branchwise cofactors.

Recommended integration:

1. Carry the compact triple `(H,a,ell_1)` in future `B/C/E` scouts.  The new
   walls `gcd(H,a)!=1` and `gcd(H,ell_1)!=1` distinguish chart degeneration
   from multiple critical branches without computing an entire inverse.
2. Do not treat “missing inverse cover” as one undifferentiated debt.
   Critical-root inversion is cheap on the squarefree locus; Keller-sheet
   inversion remains open.
3. Put the next serious calculation before scalar elimination.  If the
   gradient ideal is a unit, compute one Bezout row and test the canonical
   class `mu(P)` in `(16b)` by bounded Hamiltonian reduction.  This replaces
   a two-component mate ansatz by the single equation
   `D_P(h)=A_x+B_z`; only after that gate survives should one reconstruct
   `Q`.
4. If a future broader deformation removes all critical points, retain the
   primitive-element cofactor of THM-3064 immediately.  The pair `(U,W)` here
   is a useful audit control but cannot substitute for it.

## Reproduction

From the repository root run

```text
python3 04-computation/jc_critical_inverse_cover_cofactor_jacobian_probe_agent.py
python3 -O 04-computation/jc_critical_inverse_cover_cofactor_jacobian_probe_agent.py
```

and compare LF-normalized output with
`05-knowledge/results/jc_critical_inverse_cover_cofactor_jacobian_probe_agent.out`.
The computation is exact, deterministic, has no floating literals or Python
assertions, and reconstructs both accessory fields internally.

```text
script SHA-256:    a719b2582b93a0a6d110b1f13b65e9d54800e8669914da9f21a9371545bbae31
LF output SHA-256: 67067d9448caa6a809520b190208a561ab4cc14517455d6da0eef9210ccce1ff
```

## Honest frontier

- **Direct progress:** a universal linear subresultant/Bézout identity and an
  exact inverse-critical-graph criterion.
- **Finite-exact progress:** both named degree-53 coupled controls pass the
  two pivot-unit gates and have explicit selected-section digests.
- **Niche progress:** scalar-resultant loss is repaired on these controls.
- **Refuted/narrowed:** squarefreeness alone does not provide a uniform
  cofactor chart, and an inverse critical cover is not a Keller inverse cover.
- **Still open:** broader simultaneous `B/C/E` deformations beyond THM-3289,
  a polynomial mate, primitive branchwise Keller units, `JC(2)`, and `DC(2)`.
