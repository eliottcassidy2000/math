# Session 63: integral observers, root constructions and effectivity

## The proposed unification

I propose organizing these repository strands as **integral observers of
decorated constructions**. The common object is not the number seven, a
scalar average, or a rank. It is a construction equipped with an integral
carrier, the operations we need to preserve, and a positive/geometric
realization condition.

This is a proposed research framework, not a theorem that all repository
problems are equivalent. Its support is a collection of exact compatibility
results and hostile examples. The longer theory note states its universal
continuation principle, precise limits and falsifiable next steps.

The most important distinction is that at least four different things can
change under specialization: invariant lifting, multiplication, polynomial
jets, and effectivity. The examples below separate them rather than fold
them into a single vague notion of degeneration.

## 1. The 50-vertex example is recovered, and its full spectrum is explained

The supplied example is now located in the concurrent
`unprotected_ports_cancellation_20261025` packet. Our independent Cayley
construction is explicitly isomorphic to it by negating the F5² coordinate.
The full automorphism group is C5², and transitive-pair inflation gives
exactly 25 canonical gaps in one group orbit. The invariant constant-XOR
rank is therefore one.

For a regular C_p^k action, with p an odd prime, k>=1, n=p^k
and h=p^(k-1),

    P(t)=((1+t)^n+(n-1)(1+t^p)^h)/n.

Integer division must happen BEFORE reducing coefficients. The recovered,
independently checked spectrum is h over Q, n-p+1 in characteristic p,
and n in a prime characteristic dividing n-1; other primes retain h.
For the 25-gap example these are 5 over Q,21 modulo 5, and 25 modulo 2 OR 3.
At k=1, the characteristic-p depth remains one; not every matching prime
causes a jump.

There is a wider exact law: for EVERY odd-order action on n gaps,
P(t)≡(1+t)^n modulo 2, because each subset orbit has odd cardinality.
The maximal characteristic-two depth does not produce n legal toggles.
Diagonal copies and independently acting copies have different spectra;
the operation that joins the examples remains part of the theorem.

The integral jet makes the jump auditable. For P(-1+x)=sum b_jx^j, the
gcds of successive initial coefficient segments have nonzero levels
3000,600,120,24,1 in the 25-gap example. They also explain the depths
21,17,13,5 modulo 5,25,125,625. A field-rank statistic alone cannot retain
that filtration. The incoming cyclic C25 action is a different carrier:
it has rational depth 1 and modulo 5 depth 5, despite matching an ordinary
derivative of the elementary-abelian positive polynomial modulo 5.

## 2. New modular cancellation need not mean new invariants

For a signed permutation lattice, each basis-line orbit contributes either
R or R[2] to invariants over a coefficient ring R, depending on whether its
stabilizer preserves or reverses the chosen line orientation. The exact
integral H1 defect is a direct sum of Z/2, one per reversing orbit.

For an odd group acting on exterior powers of a permutation lattice,
all these stabilizer signs are trivial. Invariants commute with EVERY
coefficient ring. Thus the 50-vertex depth jumps occur without any new
invariant classes. They are cancellations among coefficients of an
already-defined integral polynomial.

Multiplication has a separate failure mode. In the smallest cyclic
three-coordinate exterior example, the invariant algebra has degreewise
bases 1,x,y,z but x y=3z. All four graded ranks survive reduction; the
product vanishes modulo 3. The invariant complex with differential
x-wedge has top cohomology Z/3. This is torsion created by an additional
operation, not failure of degreewise invariant base change.

This makes an actual bridge to the root-cover algebra: its character
pieces multiply with branch factors. Forgetting those multiplication
constants loses precisely the additive equations the user emphasized.

## 3. Three geometric carriers give three different twist averages

For the hourglass lines E,E±x,E±z,E±(x+z), let U be their complement,
V its free rooted cover, X the singular projective cover and Y its node
resolution. In odd characteristic other than 3, fresh 64-twist means are

    V: q²−6q+12,
    X: q²+q+1,
    Y: q²+4q+1.

The first is the count of U. The second is the count of P². The third
is the count of the blowup of P² at the three collinear triple points.
An explicit local A1 quotient calculation proves the resolved statement;
it is not guessed from a point-count correction.

In characteristic 3, the seven lines remain distinct but acquire two more
triple points. The open and resolved means become q²−6q+10 and q²+6q+1.
Characteristic 2 destroys the quadratic étale interpretation and is not
silently included.

The fixed-family extension law is more revealing. For a constant finite
abelian group G and a FIXED collection of twists over Fq, extension of
degree m has average count #(W/G^m)(Fq^m). For exponent 2, odd extensions
retain the quotient; even extensions collapse all twists to the original
cover. Over F25, for example, the fixed family from F5 has open mean 832,
not the fresh mean 487.

Consequently an average of the ordinary motives of the forms is not the
quotient motive. For the hourglass's order-64 exponent-two group, the exact zeta repair is

    product_g Z(W_g/Fq,T)
      = Z(B/Fq,T)^64
        * [ Z(W/Fq²,T²) / Z(B/Fq²,T²) ]^32,
    B=W/G.

The missing factor is exactly the even-extension character block. This
is a precise failure of commuting trivial-character projection with the
Adams character-power operation; Frobenius itself still commutes with the
constant group action.

## 4. Multiplicative roots need the additive branch map

In characteristic different from 2, over the branch complement U, the
open construction is genuinely

    V = U x_(Gm^6,[2]) Gm^6.

The power map is multiplicative. The map U->Gm^6 is defined by the six
ratios of the ADDITIVE branch forms. Dropping that map replaces the
problem by a torus and loses its geometry.

Its character algebra makes the retained equations visible:

    r_S r_T = r_(S symmetric_difference T)
              * product_(i in S intersection T) f_i.

At zero coordinates a unit-locus Laurent decoder is insufficient. The
classical iterated root stack of the seven branch divisors is naturally
[X/G], preserving local stabilizer/inertia data. The coarse space of [X/G] is P². Resolving the source first gives
[Y/G], whose coarse space is the blowup of P² at the branch triple points.
These are different functors,
not alternative names for the same object.

The incoming marked-root coordinate theorem supports this viewpoint:
invertibility modulo the root order controls marked root algebras on the
unit locus, whereas the specified monomial map of ambient tori is an isomorphism
only with determinant±1.
For quadratic roots a unimodular lift can repair a chosen Boolean exponent
presentation, but it cannot erase an intrinsic Smith torsion factor of a
different physical map or extend automatically across branch zeros.

## 5. The real connection to seven is character-resolved

The 63 nontrivial global even-subset characters split into blocks of sizes
21,35,7, of weights 2,4,6. The seven weight-six characters are precisely
the canonical/K3 double-plane sectors. In the affine six-bit coordinates,
they are the six weight-five modes plus the weight-six mode.

The mean retains NONE of these 63 nontrivial modes. Thus geometric genus 7
and the average polynomial's formal value 7 come from different observers.
For seven distinct complex lines with only double/triple intersections, the complement's
Euler characteristic is 10 minus the number of triple points, while the
stated root-resolution geometric genus remains 7. Generic seven lines
already give Euler 10, not 7. The hourglass equality is incidence-dependent,
not a universal identity between the invariants.

The full character transform restores the marked twist profile. Its
64-dimensional integer Walsh matrix has Smith factors 2^j with multiplicity
binom(6,j), total index 2^192, and rank 1 modulo 2. An abstract character vector
A is the Fourier transform of a nonnegative integral histogram exactly
when its inverse transform
is nonnegative integral, equivalently when it has the required lattice
membership and its associated character-correlation matrix is positive
semidefinite. These conditions do not by themselves prove geometric
realizability.

Exact Legendre controls show successive losses: over F7 the same branch
incidence and mean give twist vectors(4,4) and(0,8); over F5, vectors(4,0)
and(0,4) have the same mean and variance but different marked phases.
Full covariance is still not a complete marked boundary response.

## 6. The negative coefficient is a sharp effectivity test, not a ban on F1 ideas

For the open base, q²−6q+12 becomes T²−4T+7 at T=q−1. Split-torus
strata would give nonnegative integer coefficients in T, so this cannot
be such a decomposition. The general coefficient is 3−n for a complement
of n distinct complex projective lines, independently of concurrency. Complex line
complements are split-torifiable exactly for n≤2 or three nonconcurrent
lines; the three-line pencil fails by its negative constant term.

This excludes that particular positive model, not every framework called
geometry over F1. The formal q=1 value is an Euler specialization, not a
literal field-point count. Moreover positive coefficients do not solve
all lifting problems: over C the resolved quotient blowup torifies, while its
root-cover source has nonzero geometric genus and does not.

The analogous warning in tournament work is that extra polynomial factors
need not be extra legal toggles or positive exterior factors. The common
principle is a lifting obligation from formal algebra to the actual
positive/geometric category, not an identification of the categories.

## What is proved, recovered and proposed

The exact modular, signed-lattice, twist-power, root-stack and Fourier
statements have independent proof and finite-check audits. The incoming
spectrum, explicit 50-vertex witness, marked-root coordinate repair and
integral-wrap results are credited as recovered work, not rediscovered
under new names. The comprehensive organizing framework is PROPOSED:
it is supported by these maps and counterexamples but is not a theorem
that it classifies every repository object.

No global Hadwiger–Nelson, LRC, black-core connectivity, integer-square
hourglass or universal F1-descent theorem is claimed. The useful next
research is to make the observer operation-complete and then solve the
remaining realizability or boundary-lifting problem, rather than infer
it from a repeated number or scalar count. That is a future direction,
not automatic research after this session ends.

## Packaged navigation (presentation-only)

The audited body above is preserved byte for byte. This suffix only supplies
portable links and is not an additional mathematical claim.

- [Session overview and status boundaries](README.md) · [synthesis](SYNTHESIS.md)
  · [proposed framework, not a universal classification](PROPOSED_UNIFYING_THEORY.md).
- [Cancellation proof](cancellation/RESULT.md) · [cancellation audit](audits/cancellation/REPORT.md).
- [Twist carriers](twists/RESULT.md) · [root stack](twists/ROOT_STACK.md)
  · [zeta/effectivity](twists/ZETA_AND_EFFECTIVITY.md) · [twist audit](audits/twists/AUDIT.md).
- [Integral specialization](specialization/RESULT.md) · [Fourier criterion](fourier-effectivity/FOURIER_EFFECTIVITY.md)
  · [integral/Fourier audit](audits/specialization/REPORT.md).
- [Current source reconciliation](CURRENT_SOURCE_UPDATE.md) · [historical source map front door](recovery/README.md).
- [Final integral readback](audits/specialization/FINAL_ROOT_READBACK.md)
  · [final geometric/visual readback](audits/twists/READBACK_SYNTHESIS_THEORY_FIGURE.md).

![Integral carrier, multiplication, scalar jets and effectivity remain distinct](figures/integral-observer-map.svg)

[Deterministic figure generator](figures/make_observer_map.py) · [verification scope](PUBLICATION_CHECKS.md).
