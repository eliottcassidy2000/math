# Rerouting from a closed carrier to the exact clock's divisor

2026-10-11. **PROVED, scoped and self-audited:** the divisor and area-change
calculation below. **INHERITED:** the exclusion of all polynomial mates for
the displayed carrier. **PROPOSED:** the clock/divisor-first search outside
that carrier. **OPEN:** the planar Jacobian conjecture. The exact new clock
is not a Jacobian mate for the original area form.

## Recover the stronger closures first

The preceding [transfer-word packet](genus_two_transfer_words_20261011.md)
already recovered
[THM-3772, variable-flank three-charge classification](../../01-canon/theorems/THM-3772-variable-flank-three-charge-rational-exact-classification.md)
and the holomorphic time-form obstruction that survives every birational
symplectic transport. Two other stronger closures apply to the control here:
[THM-3418, one-monomial nonlinear fibre](../../01-canon/theorems/THM-3418-one-monomial-nonlinear-fiber-keller-classification.md)
and [THM-3716, monomial Broughton family](../../01-canon/theorems/THM-3716-monomial-broughton-hamiltonian-obstruction-family.md).
We do not present a finite ansatz failure as a new exclusion of this family.

The least-used sidecar is the divisor of a prospective **exact meromorphic
clock**, including which zeros are inside the actual affine source. The
hostile is an exact differential obtained only after changing the area form.
The corrected near miss is thinking that unchanged genus means unchanged
exactness, or that exactness alone has kept the polynomial Jacobian problem.

## A genus-two clock with one pole and seven zeros

Consider the inherited smooth polynomial

    H=x+x^4 t^5.

On its generic fibre H=E over C(E),

    t^5=(E-x)/x^4.                                       (1)

The degree-five map x to P1 is totally ramified at x=0,E,infinity.
Riemann-Hurwitz gives 2g-2=-10+3*4=2, hence g=2. Denote those points
P0,PE,Pinf. Their function orders are

| point | order of x (or x-E at PE) | order of t |
|---|---:|---:|
| P0 | 5 | -4 |
| PE | 5 for x-E | 1 |
| Pinf | -5 | 3 |

The standard-area time differential, up to orientation, is

    omega=dx/H_t=dx/(5x^4 t^4).

Its orders at those points are0,0,2 respectively. It is holomorphic and
nonzero; its divisor is2Pinf. Thus it cannot be the differential of a
nonconstant rational function on the projective curve.

But the rational function t itself is a perfectly good meromorphic clock:

    div(t)=PE+3Pinf-4P0.

Its exact differential has divisor

    div(dt)=Z1+...+Z5+2Pinf-5P0,                         (2)

where the five distinct affine points Zi satisfy

    x=4E/3,          256E^3 t^5+27=0.                   (3)

To verify the complete divisor, differentiate (1), or use

    dt=-(4E/x-3) omega.                                 (4)

The multiplier has a pole of order5 at P0, a simple zero at each point in
(3), and no other zero or pole. In particular the zero count is7 and the
pole count is5, giving canonical degree2. The five affine zeros are simple
over C(E); this is a generic-fibre statement, not the E=0 degeneration.

The full x-resultant is -t^15(256E^3t^5+27). Its t^15 factor comes from a
degree drop, not additional clock zeros: at t=0, H_x=1. The verifier retains
that factor and separately checks its exclusion instead of equating the
raw resultant with the affine zero polynomial.

This shows exactly how exactness is purchased on the *same* genus-two curve:
the clock acquires a pole, and its derivative acquires affine zeros. It does
not require or demonstrate a collapse to genus zero.

## The global bill is the area multiplier

On the actual affine surface,

    dH wedge dt = (1+4x^3 t^5) dx wedge dt.               (5)

Thus t becomes a unit-Jacobian companion only if one changes the source
area form to the right side of (5). That form degenerates along
1+4x^3 t^5=0, precisely the five affine clock zeros on each generic fibre.
Deleting this divisor gives a legitimate different open-surface problem;
it does not give a polynomial mate on A2.

The general accounting rule is elementary. If a nonconstant rational clock
J on a smooth projective genus-g curve has poles of orders m1,...,ms,
then dJ has poles of orders m1+1,...,ms+1, so

    degree(zero divisor of dJ)=2g-2+sum_i(mi+1).          (6)

In characteristic zero the leading coefficient at a pole cannot vanish.
Equation (6) does not locate those zeros: all may lie in a chosen boundary
on some different affine curve. The affine source's actual deleted points
are therefore indispensable. For (1) they do not contain the Zi in (3).

This sharpens the prospective search object from a bare fibre genus to

    (curve, affine boundary, rational clock, zero divisor, pole divisor,
     source-area multiplier).

It also explains why old transfer words cannot repair this example. A
birational symplectic change pulls back the original nonexact omega. Using
dt instead deliberately changes that invariant and incurs (5). Relabeling
the curve or composing the closed transfer alphabet cannot make (5) constant.

## An explicit carried defect rather than another failed ansatz

THM-3418 already proves the all-degree obstruction. Its character-one chain
specializes here to

    J_K=sum_(k=0)^K c_k x^(3k)t^(5k+1),
    c0=1,      c_(k+1)=-(5k+4)/(5k+6) c_k.

Direct differentiation gives the exact receipt

    Jac(H,J_K)=1+c_K(5K+4)x^(3K+3)t^(5K+5).             (7)

Every c_K is nonzero. Increasing K moves the defect outward and never pays
it. This is inherited, not a newly discovered proof of nonexistence. The
companion verifier checks seven truncations only to protect the specialization
and signs in (7).

The next proposed search should change the polynomial fibre support beyond
the excluded three-charge and top-plus-constant strata, then compute the
clock divisor and all affine components before enlarging coefficient ansatzes.
A useful positive candidate must both admit an exact clock and locate its
ramification entirely outside the retained affine source, with a constant
standard-area Jacobian. No such candidate is produced here.

The cross-thread connection is an operation with a retained bill: a denominator
change in [unit-path geometry](cyclotomic_unit_component_20261011.md) has prime
support, a [ROOT terminal completion](collatz_terminal_completion_20261011.md)
has a source guard and exponent phase, and a clock change has a divisor.
Each representation allows new operations; each keeps the obstruction that
those operations would otherwise hide. These are analogous proof interfaces,
not an identification of their invariants.

[Exact verifier](../../04-computation/experiments/source_context_localization_20261011.py)
and [output](source_context_localization_20261011.json) check the Bezout
identity, fibre multiplier, resultant (3), divisor degree, genus count and
the inherited truncation defects. No independent peer audit is claimed.
