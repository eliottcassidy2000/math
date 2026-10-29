# Overnight arithmetic carriers: verified checkpoints, research continuing

2026-10-29. Initial source pin `5d6eb91dd76eee86b91e04d856daa32c32b3d61d`.
**Status: ACTIVE overnight session.** This index records published checkpoints;
continued investigations use separate files until independently reviewed.
No general Collatz, rank-boundedness, Sha-finiteness, or moonshine theorem is
claimed. Classical inputs and existing curve discoveries are credited.

## Checkpoint 1: an arithmetic exceptional set is actually excluded

[Integer termination](controllers/INTEGER_TERMINATION.md) proves that the
cheaper ternary controller terminates for **every positive odd integer**.
It follows the native accelerated Collatz orbit, adding two queue units
when3 divides the valuation and subtracting one otherwise. The potential
(n+1)(3/2)^Q weakly decreases, and every branch contracts by at most27/28.
The exact positive-integer guard n>=13 is load-bearing; the rational fixed
source1/5 is a hostile to dropping it.

[Five-cylinder classification](controllers/FIVE_EXCEPTIONAL_CYLINDERS.md)
then determines precisely which first stops increase the source:

    n = 3 mod4;
    n = 125 mod128;
    n = 1357,2269,2877 mod4096.

Equality holds only at1; every other positive odd integer decreases.
This is a complete arithmetic comparison classifier, not a proof that
repeated failures eventually descend below the original source. Full
source carries and division/odd-step clocks remain necessary for retries.
[Independent proof audit](audit/CONTROLLER_AUDIT.md) and two independent
checkers accompany the results.

[Hasse quotient certificates](frobenius/HASSE_QUOTIENT_CERTIFICATES.md)
keep j-Frobenius distinct from Frobenius on a named Legendre parameter.
At5, supersingular j=0 is rational but its Legendre parameters are not.
At71, all seven supersingular j-invariants are rational, while21 rational
Legendre parameters sit over only four of them. This is a genuine typed
7/21 occurrence, not the tournament Hamilton-count exclusion. The finite
prime scan through199 does not replace the cited all-prime Ogg theorem.

## Checkpoint 2: exact rank witnesses and genuine fifteen-object geometry

[Rank and family investigation](ranks/README.md) reads the supplied ICARM
curve273 and current adjacent records. The number21 in the PPVW model is a
heuristic **all-but-finitely-many threshold**, not a universal upper bound;
finite examples of higher rank are not counterexamples to that prediction.
Its source is a height-count/tail-exponent balance, not the tournament
Hamilton-path exclusion. The selected leaderboard is not an unbiased sample
from which rank densities can be inferred.

For curves273,302,4239, exact local quadratic characters independently
certify **30,31,32 independent rational points**, respectively. No BSD,
GRH or floating regulator is used for these lower bounds. The certificate
retains each good prime, named cubic root, curve model and point. Exact
coordinate transport preserves those characters; forgetting the root
translation corrupts them. [Independent certificate audit](audit/INDEPENDENCE_CERTIFICATE_AUDIT.md)
reconstructs the matrices by different finite-field arithmetic. Conditional
upper bounds and generic family ranks remain separate statements.

We acknowledge the ICARM Elliptic Curve Rank Leaderboard, its contributors,
and NSF Grant DMS2425401 for the public curve records. Only the three
necessary curve snapshots and compact mathematical certificates are included,
not a copy of the full selected database or a research paper.

[Richelot, Frobenius and Fano](supersingular/RICHELOT_FROBENIUS_AND_FANO.md)
gives the real15-object bridge: six branch points,15 perfect matchings,
15 maximal isotropic two-torsion kernels, and Richelot isogenies. Frobenius
stability is the rationality guard. The fact that Ogg's list also has15
primes gives no canonical identification. At the Ogg prime11, two supersingular
elliptic squares have respectively15 and3 rational Lagrangian kernels.

The next genus reveals missing structure: eight branch points give135
Lagrangians, only105 from matchings. The remaining30 are extended-Hamming
structures; puncturing at a named point gives labelled Fano planes.
Frobenius cycle type(1,7) fixes two of these nonmatching kernels, with
seven-cycle normalizers of order21. The explicit finite-field example and
independent binary-subspace audit retain the arithmetic and combinatorial
maps. This does not identify their21 with Hamilton-path counts.

## Checkpoint 2: the fruit equation needs a cubic splitting guard

[Fraction-cubic investigation](fractions/README.md) recovers

    A/(B+C)+B/(A+C)+C/(A+B)=4

with exact excluded denominators, primitive scaling and positivity rules.
Its cyclic-symmetry quotient is an explicit3-isogeny. On the quotient curve,
a named infinite-order point H satisfies3H=phi(G), while H is outside the
rational image. Consequently nH lifts rationally exactly when3 divides n.

A rational orientation and three positive real roots are not sufficient:
35H lies in that positive-real chamber but its root cubic has no Q2 root.
There are infinitely many such nonlifting real-positive points. This is a
local obstruction, not a claimed Sha class. The [second-checkpoint audit](audit/SECOND_CHECKPOINT_AUDIT.md)
independently checks the isogeny arithmetic and discharges infinite order by
a nonintegral point on an integral short model and Nagell–Lutz.

## Checkpoint 3: retain the source through failures

[Source-preserving retries](controllers/SOURCE_PRESERVING_RETRIES.md) repair
four half-rays and a separate3 mod16 class below their **original** inputs.
The exact compressed valuation-one run keeps its unbounded integer clock.
For sources2^t-1 its first compressed return descends only at t=2,4,8;
this is an all-t theorem, not a finite cutoff inference.

The [25-chart failure transducer](controllers/FAILURE_TRANSDUCER.md) keeps
congruence pullbacks, full carries and both clocks. Its rank-one Haar
quotient describes random sources, not resampled deterministic integers.
An [all-period certificate](controllers/PERIODIC_EXCEPTION_CERTIFICATE.md)
shows that-1 is its only integer periodic failure cycle, and-1,-3 are its
only eventually periodic integer failure sources. This does not classify
all native Collatz cycles.

The residual [exceptional Cantor set](controllers/EXCEPTIONAL_CANTOR_RECURRENCE.md)
has exact dyadic-cell recurrence

    N_L=N_(L-1)+N_(L-6)+3N_(L-11).

Its dimension s satisfies2^(-s)+2^(-6s)+3*2^(-11s)=1, with s about0.441171.
The number of positive integers at most X whose cheap returns expand
forever is at most(17/8)X^s. This is neither an emptiness theorem nor a
bound for all nonconvergent Collatz trajectories. The exceptional set is
specified by this particular failure language.

## Checkpoint 4: changing moduli and retaining generic sections

[Rank phase2](ranks/phase2/README.md) supplies compact source records and
certificates for every one of the66 rank-at-least30 point lists in the
selected snapshot:64 work with binary characters; two require a ternary
finite-group certificate. Altogether2005 named independent point directions
are verified. Failure of a sufficient binary certificate was not mistaken
for dependence. This is a certificate for the listed curves, not a census
of all high-rank curves or a proof of an exact upper rank.

The classical Elkies948 K3 model has17 independently reconstructed sections
with exact positive height Gram determinant948. Fifteen secant choices
can be flipped by unimodular changes of the generated subgroup basis;
these bits can be forgotten for independence but not for named-point
reconstruction. Saturation in the full Mordell–Weil group is not reproved.

The three supplied record fibrations have distinct good-reduction point
counts at233:58509,58404,58374. They are therefore pairwise non-isomorphic
**over Q**, not necessarily over its algebraic closure. Nevertheless all
three retain the same verified [A3(2) three-section sublattice](ranks/phase2/COMMON_A3_CARRIER.md),
with determinant32. That small carrier does not determine the ambient K3
surface. [Separate audit](audit-k3/REPORT.md) reconstructs the arithmetic
from independently transcribed mathematical inputs, without copying a paper.

## Checkpoint 5: cyclic kernels, four ports and arithmetic image guards

The [cyclic-kernel recurrence](supersingular/CYCLIC_KERNEL_RECURRENCE.md)
classifies Frobenius cycle type(1,n), n odd: invariant Lagrangians exist
exactly when every prime p dividing n has odd ord_p(2). Reciprocal factors
are paired and one is chosen from each pair. Reflection exchanges these
choices rather than preserving them. The7-case gives the Fano/Hamming
carriers, the23-case the Golay carrier. Among the15 Ogg primes, the exact
intersection is7,23,31,47,71; [explicit finite-field models](supersingular/OGG_INTERSECTION_MODELS.md)
retain branch count and characteristic as separate quantities.

At n=31, eight kernels share their complete weight enumerator and all named
one-, two- and three-port minimum-word counts. [Four-port data](supersingular/REVERSAL_AND_FOUR_PORT_HOSTILE.md)
distinguishes every marked kernel and separates two coordinate-isomorphism
classes. In six cases it [reconstructs the entire affine F2^5 geometry](supersingular/AFFINE_RECONSTRUCTION_AND_ISOSPECTRAL_LATTICES.md).
Construction A then gives two nonisometric rank32 even unimodular lattices
with identical full theta series. This lattice rank32 is not identified
with elliptic Mordell–Weil rank32. Plain Golay Construction A gives
Niemeier A1^24, not the Leech lattice without further gluing data.

For the fruit curve, [full rational image](fractions/FULL_IMAGE.md), using
the explicitly cited source-generator theorem, yields a constant-size
arithmetic guard. For a **rational point on the quotient cover**, q=a/d
in lowest terms lifts rationally iff3 divides d or a+d. The rational-cover
condition and forbidden fruit denominators remain essential. The
[positive-real lift density](fractions/DENSITY.md) is exactly one third
under source-index or canonical-height ordering; its proof uses irrational
circle rotation, not a finite sample or an effective return-time bound.

## Checkpoint 6: an all-size reconstruction theorem and the Leech guard

[Affine port reconstruction](supersingular/ALL_M_AFFINE_PORT_RECONSTRUCTION.md)
works uniformly on F2^m for m>=4. For affine codimension-two flats, the
number containing a named set S depends exactly on its affine-span dimension:
R(S)=[m-d choose2]_2. All one-, two- and three-port responses are constant,
while four ports reconstruct the entire affine geometry after a named origin
is chosen. The number of labelled geometries forces26 and90 bits at m=4,5.
At m=3 the blocks are merely all pairs, so even the complete incidence bank
cannot distinguish the30 affine structures. The self-dual code/lattice
application remains special to m=5, not all dimensions.

The classical [Golay-to-Leech construction](leech/GOLAY_LEECH_SUM_GUARD.md)
then retains a genuinely additional arithmetic condition: a modulo-eight
sum guard plus the odd glue coset. Keeping only the code and common parity
reintroduces short vectors. With all guards present, exact enumeration gives
196560 minimal vectors in three families1104+97152+98304, and the lattice
has minimum norm4 and is even unimodular. This recovers the actual Leech
lattice rather than misidentifying the A1^24 lattice from plain Construction A.
The route from the Ogg prime23 to its cyclic binary code and then the
classical guarded lattice is explicit; it does not prove Ogg's classification
or construct the Monster from an arbitrary fifteen-element set.

## Ongoing work after these checkpoints

The current bounded extensions investigate the two binary-certificate gaps
through exact rational halving, and an older primary source on tournament
solutions. They are not automatically included in the completed claims above.
The aperiodic positive-integer failure boundary and the missing generic
section banks for newer rank-record fibrations remain open here. New results
are published only after their own proof and computational audits.

## Reproduction

Each published script states its own finite universe and uses exact arithmetic.
Checks remain active under `python -O`. Most use the standard library;
six rank-phase2 symbolic checks and one separate K3 audit declare SymPy1.14.0. Mathematical JSON
inputs are adjacent to their consumers and use stable file-relative paths.
The first controller/Frobenius checkpoint was replayed from its committed
checkout in both modes with identical output. Further checkpoint replay
receipts are recorded as publication proceeds. A script pass is not a proof
of an unstated open conjecture or an untested larger universe.
