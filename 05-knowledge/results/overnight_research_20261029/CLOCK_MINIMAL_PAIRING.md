# Minimum pairing completion with the actual seven-cycle clock

2026-10-29. **PROVED finite-module completion theorem / FINITE-EXACT
native source modules.** This extends [THREE_RUN_GLUE R6](THREE_RUN_GLUE.md),
where an abstract group repair did not yet retain its clock.

**Continuation:** [PRIMARY_CLOCK_COMPLETION](PRIMARY_CLOCK_COMPLETION.md)
now closes all primary exponents for the specified C7 action, proves an
all-height minimum recurrence, and distinguishes injective completion
from ordinary-quotient-compatible completion. Its P6 also verifies full
vertex/run invariant forms at5 and17. The following continuation addresses
the higher-primary pieces; the six elementary primes are complete.

The next [full vertex/run continuation](FULL_NATIVE_COMPLETION.md) now
closes the dyadic minimum at index2^40. At3 it proves the stronger lower
bound3^8; joint-clock attainment and ordinary-quotient-compatible minima
remain open. The original C7 formula retains its stated scope.

## C4. The clocked minimum is controlled by reciprocal constituents

Let p be odd, p not7, and V an elementary abelian p-group with a specified
order-dividing7 clock C. Require an injective clock-preserving map of V
into a finite abelian p-group W with a perfect alternating pairing and
an isometric clock extending C. No original source element may be
identified or discarded. Write the Fp[C7] decomposition as

    V = trivial^(m0) + sum_f V_f^(m_f),

where f ranges over irreducible nontrivial factors of Phi7 and
f*(X)=X^(deg f) f(1/X)/f(0) is its reciprocal polynomial. For elementary
abelian completions the exact minimum additional dimension is

    delta_C7(V)=(m0 mod2)
                 + sum_{unordered pairs f!=f*} deg(f)|m_f-m_f*|.  (C6)

Self-reciprocal nontrivial factors add nothing. Semisimplicity follows
because p does not divide7. An invariant pairing can pair a non-self-dual
irreducible only with its dual, so their multiplicities must match.
The trivial part must have even dimension. Adding precisely the missing
dual copies and, when needed, one trivial line attains the lower bound.

For completeness, a self-reciprocal nontrivial constituent has even
degree2r. Identify it with F_(p^(2r)), with clock multiplication by a
seventh root zeta. The involution bar(x)=x^(p^r) sends zeta to zeta^-1.
Choose nonzero c with bar(c)=-c. Then Tr(c*x*bar(y)) is a nondegenerate
alternating form preserved by the clock. This supplies the sufficiency
claimed for that case without assuming a chosen matrix basis is orthogonal.

The same formula and proof hold for EVERY odd cyclic clock order n with
p not dividing n: use all nontrivial irreducible factors of X^n-1 in
place of those of Phi7. Thus(C6) is an all-order classification of
elementary abelian pairing completions with a supplied odd cyclic clock.
The native calculation below specializes it to the actual order7 action.

An especially sharp case is V with no constituent isomorphic to any
constituent of V*. Every invariant bilinear form on V then vanishes.
In ANY paired finite p-group W containing V equivariantly, V is isotropic,
so |W|>=|V|^2. The canonical hyperbolic V+V* attains equality. This
bound permits higher-exponent ambient groups and needs no direct-summand
assumption. It also proves optimality for any larger clock group acting
on V: the lower bound already follows from C7, and hyperbolic doubling
extends every source automorphism contragrediently on the dual.

This theorem concerns a finite source module and a specified pairing
completion. It neither realizes the result as Sha nor bounds an arithmetic
inverse limit. The reciprocity mechanism is classical representation
theory; the calculations below supply the actual native source and maps.

## C5. The six elementary odd-primary pieces of Paley7

Use the actual210-by210 three-run matrix M from THREE_RUN_GLUE. The
translation tau:v->v+1 modulo7 permutes its named port triples and
commutes with M. The full Smith certificate shows that its primary
groups at p=5,11,17,79,5861,392809 all have exponent p.

The checker constructs a literal quotient projection B: Fp^210 -> V
from the left nullspace of M. Its selected free columns form the identity,
and it verifies BM=0 and B R_tau=C B on EVERY named port column.
The receipt retains B, its free columns, and C. Hence the following are
induced clock modules, not arbitrary vector spaces of the same dimension.

| p | dim V | Clock constituents | Abstract extra dimension | C7-equivariant extra dimension |
| ---: | ---: | --- | ---: | ---: |
| 5 | 12 | two copies of the self-dual degree6 constituent | 0 | 0 |
| 11 | 6 | two copies of one cubic, zero of its dual | 0 | 6 |
| 17 | 4 | four trivial constituents | 0 | 0 |
| 79 | 3 | one cubic, zero of its dual | 1 | 3 |
| 5861 | 6 | two copies of one cubic, zero of its dual | 0 | 6 |
| 392809 | 3 | one cubic, zero of its dual | 1 | 3 |

At the four non-self-dual rows every clock-invariant bilinear form is
zero. A separate linear-system calculation verifies that its entire
solution space has dimension0, without relying on a factor-count formula.
The explicit doubled clock diag(C,C^(-T)) preserves the standard
hyperbolic pairing. Thus doubling is the minimum pairing completion
for these four primary pieces, even while preserving ALL original
source clocks and names. For the5 and17 rows, the continuation above
now supplies literal perfect forms for the full vertex/run clock group.
The2- and3-primary groups have higher exponents; their C7 completion is
also resolved there, with the extra clocks still separate.

The11 row is the strongest warning: its underlying group(F11)^6 already
has even dimension and square order. Neither test reveals that the
supplied clock forbids every nonzero invariant bilinear pairing.

The same reciprocal-orbit mechanism appears in Session97's
[cyclic Frobenius Lagrangian count](../session97_overnight_arithmetic_carriers_20261029/supersingular/CYCLIC_KERNEL_RECURRENCE.md).
That problem starts with a symplectic binary module and chooses one
side of each reciprocal pair for an invariant Lagrangian. Here the
source may occupy only one side and the task is to add the missing dual.
The maps are different. In particular(characteristic,clock length)=(11,7)
has Frobenius degree3, whereas(2,11) has degree10. Keeping those two
parameters prevents turning the repeated11 into a false identification.

## C6. Frobenius phase explains the5/11 distinction

The irreducible degrees of nontrivial seventh roots over Fp are ord_7(p).
For p=5 this degree is6 and the Frobenius orbit contains inverse roots,
so the constituent is self-dual. For p=11, p=4 modulo7 has order3;
the two inverse orbits{1,2,4} and{3,5,6} remain distinct. The actual
source occupies just one of them twice. Clock-preserving pairing requires
adding the missing reciprocal orbit. The same mechanism applies to
79,5861 and392809, whose residues are2,2 and4 modulo7.

The factors are certified explicitly, without an opaque large-field
factorization step. For the cubic case choose a square root r of-7 in Fp
and set eta=(-1+r)/2, bar(eta)=(-1-r)/2. Then

    Phi7(X)=(X^3-eta X^2+bar(eta)X-1)
             (X^3-bar(eta)X^2+eta X-1).                 (C7)

The product is checked literally. Since each factor has degree ord_7(p),
each is irreducible. Evaluating these polynomials on C determines every
constituent multiplicity. At p=79 the occupied cubic is

    X^3+13X^2+12X-1,                                   (C8)

and its missing dual is X^3+67X^2+66X-1. The numbers12 and13 arise
because12^2+12+2=158=2*79. This is a specific quadratic-period equation
and an actual native clock module. It is distinct from the sampled
Ramanujan unit's11/12/13 clock in
[CLOCKS](../weighted_crt_clocks_20261028/CLOCKS.md).

The older [Ramanujan comparison, R5](../cursor_cm_phase_20261028/RAMANUJAN.md)
also contains indices79 and158, but uses CM discriminants-79 and-632.
Here the period field is Q(sqrt(-7)). Equal integers do not identify
these fields or their marked algebraic elements; no such transport has
been supplied. This records the tempting connection and its stopping
reason rather than promoting it from a numerical coincidence.

## Audit and next test

The independently computed full Smith form establishes the six exponent-
one source groups. Row reduction, the literal intertwiner BR_tau=CB,
C^7=I, cyclotomic factor products, complete nullities and invariant-form
linear systems then certify the clock obstruction. The source projections
make the calculation reproducible without losing names under a basis change.

The continuation now computes every2/3 C7 layer with reversible source
charts. Next extend the minimum or pure completion to the remaining
affine/run clocks at those primes. The canonical hyperbolic double is
already available. An abstract summand repair, a count recurrence and
a spectral factorization do not justify a smaller carrier without the
actual action and the required transition maps.
