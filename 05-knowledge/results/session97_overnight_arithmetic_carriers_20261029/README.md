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

## Ongoing work after these checkpoints

The active extensions concern source-preserving failure retries and their
periodic versus aperiodic boundary; cyclic Frobenius kernels and higher-port
code invariants; and exact high-rank K3 section/family certificates. These
are not automatically included in the completed claims above. New results
are published only after their own proof and computational audits.

## Reproduction

Each published script states its own finite universe and uses exact arithmetic.
Standard-library checkers remain active under `python -O`. Mathematical JSON
inputs are adjacent to their consumers and use stable file-relative paths.
The first controller/Frobenius checkpoint was replayed from its committed
checkout in both modes with identical output. Further checkpoint replay
receipts are recorded as publication proceeds. A script pass is not a proof
of an unstated open conjecture or an untested larger universe.
