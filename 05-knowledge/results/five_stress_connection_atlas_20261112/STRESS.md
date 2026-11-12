# Why an isothermal stress net has five parameters

**PROVED IN THIS NOTE:** the local classification below, its period tests,
sensor failures, and translation algebra. **FINITE-EXACT:** symbolic controls
in `stress_receipt.json`. These are self-contained derivations, checked by the
same agent using different representations, not an independent-agent audit.
The classical local family is attributed to Neményi; no historical priority
is claimed for the derivation or its geometric interpretation.

## Scope and source boundary

Assume homogeneous isotropic linear plane elasticity, no body force, and a
connected open plane domain. Local compatibility is distinguished from a
globally single-valued displacement. No boundary tractions are prescribed in
advance. Fix an isothermal orthogonal net locally written as the level curves
of the real and imaginary parts of a holomorphic function `w=Z(z)`, with
`Z'(z) != 0`. Put `h=(Z')^2`. The axes are unordered: exchanging the larger
and smaller principal stresses preserves the net. At isotropic points the
stress itself does not specify unique principal axes; the prescribed net is
their continuation. If unique axes are required everywhere, restrict to
parameters/domains on which the stress difference never vanishes.

Primary bibliographic records identify P. Neményi, *Stromlinien und
Hauptspannungstrajektorien*, ZAMM 13(5), 364–366 (1933),
[DOI](https://doi.org/10.1002/zamm.19330130506), and his
[1951 review](https://www.sciencedirect.com/science/article/pii/S0065215608703004).
The publisher's [Sáenz paper](https://link.springer.com/article/10.1007/s00707-008-0118-3)
states a broader 2–5 parameter theorem on a simply connected cross section.
Its [2010 erratum](https://doi.org/10.1007/s00707-010-0370-1) is confirmed by
Crossref metadata; its text was not accessible in this session. We therefore
do not import the broader theorem's proof or assert what the erratum changes.
The restricted result used here follows from the displayed equations.

## The local classification

Write the real stress tensor using

    p=(sigma_xx+sigma_yy)/2,
    q=(sigma_xx-sigma_yy)/2 - i*tau_xy.

The two equilibrium equations are equivalent to `q_barz=-p_z`. Elastic
compatibility is `Delta p=0`, so `q_barzbarz=0`. A tangent to a level curve of
Z has direction `theta=-arg Z'` or its perpendicular. Therefore the prescribed
principal-axis condition is exactly

    q=h*lambda,     lambda real.

Since h is holomorphic and nonzero, `lambda_barzbarz=0`. Reality gives
`lambda_zz=0` too. In real coordinates these say `lambda_xy=0` and
`lambda_xx=lambda_yy`. Integrating locally, and then gluing the constant
coefficients over the connected domain, gives precisely

    lambda=A|z|^2+Bz+conj(B)*conj(z)+C,
    A,C real, B complex.                                  (1)

Conversely every (1), with

    F=h*(Az+conj(B)),
    p=p0-2 Re integral F dz,      p0 real,                 (2)

satisfies equilibrium, compatibility and the axis condition wherever the
integral is defined. Thus the local solution space has **four real amplitude
parameters plus one hydrostatic parameter**. The parameterization is injective
on every nonempty open domain, so its dimension is exactly five, not merely
at most five.

For a simply connected domain the primitive in (2) exists globally. Local
strain compatibility then also gives a global displacement, up to rigid
motion. This is the setting in which the informal statement needs no extra
period data.

## Two different global period tests

On a multiply connected domain, first ask only for a single-valued stress
that is locally compatible. The necessary and sufficient test is

    Re integral_gamma F dz = 0 for every closed loop gamma.   (3)

These are real linear constraints on A, Re B, Im B; C and p0 survive. If their
rank is r, this stress space has dimension `5-r`, where `0<=r<=3`.

Single-valued displacement is a stronger condition. For the standard complex
elastic potentials choose

    phi''=-F,       psi'=-G,       G=h*(Bz+C),
    2 mu (u+i v)=kappa*phi-z*conj(phi')-conj(psi),           (4)

where for plane stress `kappa=(3-nu)/(1+nu)` and `mu>0`. Direct differentiation
of (4) recovers the stress and isotropic strain relations. For a physical
Poisson ratio, `kappa+1 != 0`. Global displacement exists exactly when

    integral_gamma F dz = 0,
    kappa*integral_gamma zF dz
        + conjugate(integral_gamma G dz) = 0              (5)

for every loop. These are real linear constraints on A, Re B, Im B, C; p0
remains free.

To prove necessity, (3) first makes the continuation increment of phi' purely
imaginary. In (4) it contributes a rigid-rotation increment proportional to
`(kappa+1)z`, which cannot vanish as a displacement function unless that
increment is zero. This is the first condition in (5). Now phi' is globally
single-valued and integration of `d(z phi')` gives
`integral phi' dz = integral zF dz`. The remaining translation increment of
(4) is the second condition. Conversely these conditions remove both
increments and make (4) single-valued. No displacement conclusion is inferred
from (3) alone.

For rational h, write at each excluded pole

    R=Res h, S=Res zh, V=Res z^2 h.

Then (3) is `Im(A*S+conj(B)*R)=0`; (5) is

    A*S+conj(B)*R=0,
    kappa*(A*V+conj(B)*S)-conj(B*S+C*R)=0.                (6)

The exact examples below use `kappa=3` (`nu=0`). Take a connected bounded
domain enclosing the indicated poles as holes, also excluding zeros of Z'.
Those additional holes have zero residue and add no constraints.

|Conformal coordinate Z|Single-valued, locally compatible stresses|Stresses with global displacement|
|---|---:|---:|
|z|5|5|
|z+1/z|5|2|
|z+i/z|4|2|
|z+1/z+i/(z-1)|3|1|
|z+1/z+i/(z-1)+(1+i)/(z-i)|2|1|

For `Z=z+i/z`, `h=1-2i/z^2-1/z^4`. The loop integral of F around zero is
`4*pi*A`, so stress single-valuedness already forces `A=0`. For either
one-pole example, global displacement additionally forces `B=0`, leaving C
and p0. For the displayed two-pole example the complex F periods force
`A=B=0`, and the G period forces `C=0`. The remaining one-dimensional family
is hydrostatic and does not itself determine a unique trajectory net.

This is a global qualification, not a refutation of a simply connected
theorem. Boundary conditions can reduce the family further.

## Circle geometry and what measurements retain

Equation (1) is the Hermitian form

    lambda=[conj(z),1] [[A,conj(B)],[B,C]] [z,1]^T.       (7)

The four-dimensional real space of these matrices has determinant
`AC-|B|^2`, a quadratic form of Lorentz signature after a linear coordinate
change. The zero set of lambda is a circle, line, point, empty set or the
whole plane, according to the coefficients. The fifth, hydrostatic coordinate
does not affect this locus or the principal axes.

A fractional linear substitution transforms (7) by Hermitian congruence,
with the usual positive denominator `|c*z+d|^2`. This is an exact statement
about the amplitude forms. It is **not** ordinary conformal covariance of the
entire physical stress tensor: pressure must be reconstructed from (2), and
periods must be checked again.

Four signed lambda evaluations recover its four coefficients iff the
evaluation matrix with rows `(x^2+y^2,2x,-2y,1)` has rank four. For four distinct
sites, failure means that a nonzero generalized circle contains them. One
additional pressure measurement then fixes p0. Directions alone fix none of
these five parameters.

Even full stress measurements at three noncollinear sites need not determine
the field. For the Cartesian net h=1,

    sigma_xx=2A y^2-4 Im(B)y+C+p0,
    sigma_yy=-2A x^2-4 Re(B)x-C+p0, tau_xy=0.             (8)

The nonzero field `(sigma_xx,sigma_yy)=(2y(y-1),-2x(x-1))` vanishes at **all
four** corners of the unit square. Three generic sites such as
`(0,0),(1,1),(2,3)` do give rank five. The exact sensor ranks are 4,5,4,4 for
two generic sites, those three sites, three square corners, and four square
corners respectively.

The connection to tournament endpoint work is methodological and precisely
typed: a scalar observation is a linear map with a kernel; continuation may
depend on that kernel. It does not identify stress parameters with path counts.

## A natural five-dimensional algebra that is not the compiler algebra

Cartesian translations act on (8). Their infinitesimal generators D_x,D_y
satisfy

    D_x D_y=0, D_x^3=D_y^3=0,
    span(1,D_x,D_y,D_x^2,D_y^2) has dimension five.

Consequently their operator algebra is exactly

    R[X,Y]/(XY,X^3,Y^3).                                (9)

The quadratic field `(2y^2,-2x^2)` is cyclic: its successive x and y derivatives
span the whole stress space. Thus the stress space is the regular module of
(9). This is commutative and local, with radical layers of dimensions 2 and2
after its one-dimensional scalar quotient.

The recovered tournament compiler algebra
`[[d,o,j],[0,l,i],[0,0,d]]` is also five-dimensional but is noncommutative, has
a two-dimensional semisimple quotient, and radical layers 2 and1. These are
isomorphism obstructions. Merely matching the number five would erase the
operation each space is designed to describe.

For a further exact scale family, replace the biharmonic Airy equation by
`Delta^m Phi=0`, retain Cartesian principal axes, and let `m>=2`. Equilibrium
gives `sigma_xx=f(y),sigma_yy=g(x)`, while the new compatibility says
`f^(2m-2)+g^(2m-2)=0`. Hence f,g are polynomials of degree at most `2m-2` with
opposite top derivatives. The dimension is `4m-3`, and the translation algebra
is `R[X,Y]/(XY,X^(2m-1),Y^(2m-1))`. This explains the sequence 5,9,13,17,...
inside a specified mathematical model. For m>2 this is not ordinary linear
elastic compatibility.

Translation along a fixed direction has minimal polynomial dividing
`(t-1)^3` on the ordinary stress space. A tournament clone response containing
`2^a` or `3^a` cannot be an intertwining linear image of that translation
representation. Both systems admit finite descriptions, but their recurrence
spectra differ. The useful bridge retains the operators, not just dimensions.

## Reproduction and board update

Run `python3 -X utf8 .../verify_stress.py` and its `-O` version. The symbolic
engine checks equilibrium and compatibility for three nets, exact residue
matrices for all five examples, two sensor methods, and the derivative algebra.
The proofs above establish the quantifiers; finite checks are controls.

Board consequences: stress fibres supply an explicit missing-coordinate model;
global periods show why local existence does not establish a global object;
the Hermitian carrier gives actual circle tests; the compiler comparison
rejects an algebra identification; and the history search should classify
fives by their operations. The tournament map in `SAMPLING.md` is a separate,
literal construction rather than this methodological analogy.
