# A classical affine surface realizes the surface constraints but collapses the carrier

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.
Classical surface identification is CITED. No Keller realization
or identification of an actual finite envelope is asserted.**
September 6, 2026. Construction supplied by root and independently checked
by `orthogonal_returns`.

## 1. Inheritance and precise object

Work over `k=C`.

The [alternating-envelope note](planar_jc48_sep06_alternating.md) forces
an actual six-sheet candidate's normal affine envelope to contain `A2`
with one boundary prime, `Cl=Z[D]`, canonical Weil class `2[D]`, and
compactly supported Euler characteristic two. The class-group assertion
inherits [THM-3922, affine-plane open boundary basis](../../01-canon/theorems/THM-3922-affine-plane-open-boundary-basis-class-group-obstruction.md).
These are necessary surface conditions; they do not specify a finite map.
The source carrier comes from
[the Hamiltonian theorem](planar_jc48_sep06_hamiltonian.md) and its
[discrete-descent obstruction](planar_jc48_sep06_discrete_carrier.md).

Here is an explicit smooth affine surface realizing all the listed
**bare-surface** conditions, even with a globally exact canonical form:

```text
Y=P1_s x P1_z,          Gamma={z=s^2},
W=Y\Gamma.                                                (1)
```

The graph includes the point `(s,z)=(infinity,infinity)`. Relative to
the two product hyperplane classes it has class `(2,1)` and square four.
It is an ample section of the ruled surface `Y=Sigma_0`.

This is a classical Danilov–Gizatullin surface `V4`. The primary source
Flenner, Kaliman and Zaidenberg,
[*On the Danilov–Gizatullin isomorphism theorem*](https://archive.mpim-bonn.mpg.de/2481/1/preprint_2008_83.pdf),
§1, defines these surfaces as complements of ample sections of Hirzebruch
surfaces. Its Theorem1.1 identifies their isomorphism type by the section's
self-intersection. The actual primary PDF was read. We use that result
only to identify this familiar surface, not to classify possible Keller
envelopes. The construction and the calculations below are direct.

The concept board is: a specified affine completion; class and canonical
divisors; exact differentials; extension of the source carrier; and
boundary separation. The inherited near miss is to treat surface
invariants as sufficient to identify or exclude an envelope. The new
hostile is (1). The previously missing sidecar is a global function
that separates boundary points while the chosen carrier does not.

## 2. Complete affine charts and the boundary

Over finite `s`, use

```text
U0=A2_(x,t),        x=s,        t=1/(z-s^2).
```

This includes `z=infinity` as `t=0`; its projective second coordinate
is `[1+x^2t:t]`. Over finite `r=1/s`, use

```text
Uinfty=A2_(r,b),    b=z/(1-r^2 z),
z=[b:1+r^2b].
```

Both projective coordinate pairs have no common zero. They cover all
of `W`, since the two base charts cover `P1_s` and each removes exactly
the graph point from its fibre. On their intersection,

```text
x=1/r,          t=-r^2-r^4 b,
r=1/x,          b=-x^2-x^4t.                              (2)
```

The inverses are exact. The remaining fibre over infinity is

```text
D={r=0} isomorphic to A1_b,        W\D=U0=A2.             (3)
```

The surface is smooth. It is affine without any classification theorem:
the line bundle `O(2,1)` is very ample, via the degree-two Veronese
embedding of the first factor followed by the Segre embedding. Its
section defining `Gamma` is a hyperplane coordinate. Thus `W` is the
closed image of `Y` intersected with that affine projective chart.

## 3. Class group, canonical class, Euler characteristic, and exact form

Let `H_s=(1,0)` be a fibre of the projection to `s` and `H_z=(0,1)`
the other product hyperplane class. Localization of divisor classes
and smoothness give

```text
Pic(W)=Cl(W)=Z^2/Z(2,1)=Z[D].                             (4)
```

The primitive generator is `H_s|W=[D]`: the quotient map is
`(a,b)->a-2b`, and the columns `(1,0),(2,1)` form an integral basis.
In particular the generator is not a nontrivial multiple of a hidden
class. Restricting `K_Y=(-2,-2)` gives

```text
K_W=2[D].                                                (5)
```

The boundary-basis theorem independently gives (4); the canonical value
(5) uses the displayed restriction of the canonical class.
Every unit of `W` restricts to a unit of its dense open `A2`, so is
constant. Euler additivity gives independently

```text
chi_c(W)=chi(P1 x P1)-chi(Gamma)=4-2=2
        =chi_c(A2)+chi_c(A1)=1+1.
```

There is a stronger differential control. The ordinary source volume
form on `U0` extends to a global regular two-form with

```text
omega=dx wedge dt = r^2 dr wedge db,
div(omega)=2D.                                           (6)
```

It is globally exact. Its primitive on `U0` is `alpha=x dt`, and (2)
gives the regular expression on `Uinfty`

```text
alpha=(-2-4r^2b)dr-r^3 db,        d alpha=r^2 dr wedge db.
```

These forms agree on the overlap. Exactness does not make the canonical
line bundle trivial: `omega` has its specified zeros on `D` and is not
a nowhere-vanishing global trivialization. Thus even adding this exact
differential to (3)–(5) does not by itself contradict existence of a
smooth affine surface.

## 4. Exact extension and loss of boundary information

On the original source chart use precisely the inherited functions

```text
u=x^2t,       p=t(1+u),       y=xtp.
```

Their expressions in the other chart are

```text
u=-1-r^2b,
p=r^4b(1+r^2b),
y=-r^5b(1+r^2b)^2.                                      (7)
```

They therefore belong to `O(W)` and generate a subalgebra
`B=k[u,p,y]`. On the entire boundary,

```text
(u,p,y)|D=(-1,0,0),
G_H|D = (-u/2+H(p,y))|D = 1/2+H(0,0)                    (8)
```

for every polynomial `H`. Consequently **no morphism `W->A2` whose
two coordinate functions both belong to `B` is finite**. Indeed every
polynomial in those generators is constant on `D`, so that morphism
has a fibre containing the whole affine line. A finite morphism has
finite, zero-dimensional fibres.

This loss comes from the chosen subalgebra, not the surface. The global
regular function

```text
v=x^2(1+u)=-b                                            (9)
```

is nonconstant on `D`; in particular `v` lies outside `B`. It separates
the boundary points that (7)–(8) identify. The first failed implication
would be to assume that extension of the source form ensures that its
carrier contains enough functions to define a finite map.

There is a second essential distinction. Original candidate coordinate
polynomials `A,C` are not known to belong to `B`. Knowing only

```text
P(A,C)=C^2-A^3+(3/4)A+1/4=G_H
```

does not allow the preceding finite-map obstruction to be applied to
them. Restriction to `D` gives one equation for their boundary values,
not individual constancy. For the genuine normalization `H(0,0)=0`,
the boundary level is a nodal cubic and has the nonconstant parametrization

```text
A_D=q^2-1,       C_D=q(q^2-3/2),       P(A_D,C_D)=1/2.
```

This is an exact hostile to the boundary-level inference that a constant
nodal-defect value forces both coordinates constant. It is not a solution
of the original global equation on `W`, not a bracket-one pair, and not
a proposed Keller map.

## 5. A complete polynomial shear family on the surface

There is a genuine all-order source operation on this completion which
**changes** the old carrier. For a polynomial `f(t)`, let

```text
Phi_f(x,t)=(x+f(t),t).                                    (10)
```

**Shear theorem.** This source automorphism extends to an automorphism
of `W` if and only if `f(0)=0`. In that case it fixes `D` pointwise,
preserves `omega`, and changes `alpha` by a globally exact one-form.
For every nonzero `f`, it fails descent into the old `A=k[p,y]`.

To prove extension, put `a=1+r^2b` and `delta=f(-r^2a)`. The exact new
coordinates on the second chart are

```text
r'=r/(1+r delta),
b'=b+2(delta/r)(1+2r^2b)+delta^2(5+6r^2b)
     +4ra delta^3+r^2a delta^4.                          (11)
```

These formulas follow directly from `r'=1/(x+f(t))` and
`b'=-(x+f(t))^2-(x+f(t))^4t`. If `f(0)=0`, then `delta` is divisible
by `r^2`. Thus `b'` is a polynomial in `r,b` and restricts to `b` on
`D`; `r'` is regular on `V={1+r delta!=0}`, which contains all of `D`.
The two domains `V` and `U0` cover `W`. Their maps agree on the overlap
and glue. Replacing `f` by `-f` supplies the inverse, since `t` is
unchanged; the inverse identities on the dense source chart hold on
all of `W`. This is a global automorphism, not merely a formal chart map.

Conversely, if `f(0)=c!=0`, then `delta=c+O(r^2)` and (11) has the
nonzero pole `2c/r` in `b'`. The function `b=-x^2-x^4t` is globally
regular on `W`, so its pullback could not have that pole under a regular
extension. This already excludes even a morphism extension of (10).
In particular it excludes an automorphism extension; a global source
automorphism extending in both directions would also have to preserve
the complement `D` of `U0`.

The source Jacobian of (10) is one. Equality of the regular forms on
the dense source chart therefore gives `Phi_f^*omega=omega` globally.
Moreover `t=-r^2a` is a global regular function, and

```text
Phi_f^*alpha-alpha=f(t)dt=dF(t),         F'(t)=f(t).       (12)
```

For each fixed admissible `f`, the maps `Phi_(c f)` form a polynomial
additive-group action for all `c in C`, on both the source and `W`.
Joint regularity follows from the same cover `U0` and
`{1+c r delta!=0}` in `C x W`, which contains all of `C x D`.
This genuinely integrates the locally nilpotent derivation
`f(t) partial_x`, whose source Hamiltonian is `F(t)`. It does not retry
the excluded class of Hamiltonians in `A`: a nonconstant polynomial in
`t` cannot lie in `A`, as the next restriction test also shows.

Here is the exact failure of old-carrier descent. On the component
`M={1+x^2t=0}` of the old carrier's collapsed fibre, put `x=q`,
`t=-q^-2`. The transformed `p` is

```text
Phi_f^*(p)|M=2q^-3 f(-q^-2)+q^-4 f(-q^-2)^2.             (13)
```

If `f` has degree `d>=0` and nonzero leading coefficient `c_d`, its
second term has the unique highest pole, order `4d+4`, with coefficient
`c_d^2`; the first term has pole order at most `2d+3`. Thus (13) is
nonconstant. Every element of `A` is constant on this collapsed fibre,
so `Phi_f^*(p)` is not in `A`. A nonconstant `F(t)` likewise varies on
`M`, proving that its Hamiltonian lies outside the previously excluded
carrier. This is compatible with the one-way rigidity theorem in the
[discrete carrier note](planar_jc48_sep06_discrete_carrier.md).

There is a sharp remaining obstruction even after this carrier change.
Let `B_orbit` be the subalgebra of `O(W)` generated by all
`Phi_f^*(B)` with `f(0)=0`. Since every such automorphism fixes `D`
pointwise, every translated generator still has its constant value
from (8) on `D`. Sums and products preserve this property. Consequently

```text
B_orbit subset {h in O(W): h|D is constant},
```

and no map to `A2` with both coordinates in `B_orbit` is finite.
In particular none of these operations, even used to generate a larger
algebra, recovers the separator (9). The shear family preserves the
polynomial source chart, the surface and its canonical form; it does
not preserve the original affine source family `-u/2+H(p,y)` for all
`H`, and does not by itself pay any genuine fixed-source or finite-cover
obligation. This is a positive changed-carrier operation with a precise
global stopping obstruction.

## 6. Connection contract and bounded conclusion

Source: the explicit classical surface (1) with its specified `A2` chart.
Target: the bare-surface invariants from the actual-envelope reduction,
and the extension of the inherited carrier. Map: the literal chart
transition (2) and restriction of global functions and forms. Preserved:
the open affine plane, the boundary divisor, canonical order, and exact
source formulas. Lost if only the carrier is retained: the separator
(9) and hence boundary point information. Lost if only the surface is
retained: any finite map, its mapping degree, inertia, or monodromy.

The construction shows that the listed surface/class/differential
constraints are simultaneously consistent. It does **not** identify
an actual finite envelope with `W`, and no ramification index or `A6`
cover has been constructed on it. The strongest proved obstruction is
the conditional finite-map statement for **both coordinates in `B_orbit`**,
which includes `B` and all admissible shear-translated carriers.
The next decisive information is the actual global coordinate algebra
and boundary map, not another comparison of these surface invariants.

The [standalone source](../../04-computation/planar_jc48_sep06_dg_surface.py)
checks both chart inverses, projective coordinate pairs, form and primitive
gluing, all carrier extensions, primitive Picard lattice classes, the
Euler counts, a boundary separator, the nodal-level hostile, and the
general shear formula with its necessary pole and inverse. Named shear
controls are `f=t,t^2,t+t^3`. These
are symbolic identities for the one declared surface, not a census.
The universal restriction assertion follows from the generator values
(8), rather than from finitely many sampled polynomials.

```sh
python3 -B 04-computation/planar_jc48_sep06_dg_surface.py
python3 -B -O 04-computation/planar_jc48_sep06_dg_surface.py
```

Both modes pass **62 always-active gates** with byte-identical output,
filed as [the frozen transcript](planar_jc48_sep06_dg_surface.out).

```text
source SHA-256:
f41184dbb2e0f6cc4159747e3149e2ccdb8d7af0c011aaff6ddb1024ff6fcb57
output SHA-256:
d28d37bf908fbd196b679117eb4a43f2ec7b37b76f375ca540d822bb09150641
semantic SHA-256:
79ed292dbd9f4fabb7a73b968e54ed1fc8cec97a24c0be7994e667a9ccbd0b95
```

Source and output are frozen. The [independent root audit](planar_jc48_sep06_dg_surface_audit.md)
passes the geometry, all-degree shear extension, and entire orbit-algebra
obstruction. No theorem ID is reserved by this result.
