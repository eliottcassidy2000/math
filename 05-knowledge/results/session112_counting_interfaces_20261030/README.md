# Session 112: exact counting interfaces and constructive bank ranking

30 October 2026. **Final research packet: proved statements, explicit dependencies and finite checks.**

The inspiration was Rudich–Rousseau's [*Deterministic Factoring in Time
N^{o(1)}*](https://hexagonmath.org/pdf/2610.00600v1). The useful transferable
pipeline is precise: construct a count with honest witness multiplicity,
retain the information needed by scaled modular observations, then turn exact
counts into search. In symmetry problems, stabilizers, shared residues and
canonical ordering must survive every step.

[Paper scope](paper/PAPER_SCOPE.md) distinguishes its stated theorem from a
stronger conjectured rate and reports inspected formal source and external
registry checks. We did not independently build its Lean proof or reproduce
its general counting engine. None of the following results imports a new
factoring complexity bound.

## 1. A constructive rank/unrank algorithm for the shared-residue bank

[Bank ranking](counting/BANK_RANKING.md) turns the previous total-count
recurrence into an implemented canonical selection algorithm. At completed
cycle blocks only k+2 stabilizer states are needed. A finite threshold
automaton counts primitive necklaces below each value; exact continuation
weights then rank and unrank whole bank states in polynomial time in the
number of bank coordinates D and log Q.

The implementation retains the named cycle order, cycle bases and fixed
orientation origin. It recovers 692,224 orbits at 24 binary coordinates and
the earlier 124-coordinate count, now with canonical selection. Its input
rank operation requires a canonical tuple and checks that promise.
The bound is in D=p^k-1, not in k. This is a specified bank action, not a full
tournament-isomorphism algorithm.

[The general reduction](counting/README.md) explains why ordinary prefix
orbit counts are not additive: an orbit can meet both children. Saturated
cumulative cuts repair rank semantics, but their overlap arrangement can be
huge. The specialized bank algorithm avoids that unrestricted arrangement
through its actual stabilizer structure; it does not make every such
arrangement compact.

## 2. Prime-bank counts retain one essential inner twist

[Prime-bank convolution](symmetry/PRIMENESS_INTERFACE.md) applies the inherited
canonical ideal-module decomposition before averaging over the acting group.
The same element acts on quotient and inner bank, with an exceptional twist
at the first ideal layer. Dropping it produces a negative false count.
Keeping it gives exact all-depth prime-bank orbit counts using k+2 action
types and polynomial bit cost in the explicit bank dimension.

At p=5, depth 2 has 665,856 prime orbits for the specified order-25 action.
At depth 3 some nonidentity actions fix prime banks, so filtering only free
orbits is invalid. A short rational recurrence explains these scalar-fixed
prime banks at arbitrary depth and yields prime banks with unbounded extra
scalar symmetries. The complete universal group is the intersection over all
banks; individual prime banks can exceed it by an unbounded factor. The
displayed subgroups do not determine their full automorphism groups. Prime counts are not yet a prime rank/unrank algorithm,
and these are not full-normalizer or graph-isomorphism counts.

## 3. Modular observations and reversible arithmetic coordinates

[Scaled modular observers](arithmetic/MODULAR_OBSERVERS.md) identify exactly
what readings of sD retain: D modulo lcm(m_i)/gcd(lcm(m_i),s). Avoiding the
scale's prime support gives a small-modulus zero test whose cutoff depends
on the bound for D, not the exponent of the scale. An exact deletion criterion
handles corrupted readings. Burnside's numerator supplies a genuine
stabilizer-sensitive application. Computing the observed residues still has
its own cost.

[Paired aliquot residues](arithmetic/ALIQUOT_MODULAR.md) recover n modulo p^k
from s(n) and s(p^t n) modulo p^k, even when a large valuation cannot be read
fully. The decoder is a Latin-square coordinate map. Retaining the scaled
divisor sum of the p-free core conjugates every prime-ray update to simple
scaling of the source coordinate. Shifted probes incur a sharp precision
debt that arbitrarily many later readings cannot repair.

[Weighted observers](arithmetic/POWER_OBSERVERS.md) recover n^w instead.
Their exact root fibres expose genuine ambiguities, including four prime
sources with identical square-weighted dyadic readings along the whole ray.
Two coprime weights recover n on units. These are exact information and
dynamics results, not cheap methods for evaluating divisor sums.

## 4. Signed norms meet the old divisor-sum code

[Norm signs](norms/NORM_SIGNS.md) express Euclidean norm composition in the
normalized signed-XOR model as binary rectangle equations. Four probes
recover a promised-valid eight-square table. Five explicit equations give
a contradiction in dimension sixteen, within this model.

An exact coordinate map identifies the 16 valid tables with one affine
syndrome class in the old divisor-sum Hamming bank: squareclass 30 selects
the composing orientations. This is an explicit finite correspondence, not
an identification of integers with octonions. [Norm defects](norms/NORM_DEFECT.md)
show that every invalid promised line table fails exactly 24 of 112 rectangle
tests and has mean-square defect 384 for independent centered variance-one
inputs. Unsigned support loses the relevant information; coherent sign-gauge
changes preserve it.

## 5. One pin move repairs the remaining reflection phase

[Random selectors](symmetry/PROOFS.md), conditional on the inherited universal
colored-group classification, separate residue coverage from the remaining
reflection phase. At 625 vertices exactly 1/125 of random minimal transversals
retain that phase. Moving any one pin by delta repairs every bad transversal
without adding pins. Exact probabilities and stopping-time formulas preserve
the separate depth-one boundary. These are colored-carrier base results,
not automorphism classifications for arbitrary individual tournaments.

## Incoming synthesis

The adjacent [pinned-attachment packet](../pinned_ramified_counts_20261030/README.md)
combines the random-axis theorem with attachment rigidity. At 625 vertices,
a random minimal transversal followed by independent fair binary colors
retains an unsigned reflection with probability exactly 2^(-300)/125:
the exceptional transversal has probability 1/125, and its unique reflection
exchanges 300 vertex pairs. This credited elementary deduction was reviewed
against its source; it is separate from this packet's nine checker receipts.

## Reproduce and continue

Run `python3 verify_session.py` for all nine mathematical checkers, normally
and with optimization. It compares all 18 parsed outputs with the finite
receipts. [Verification scope](CHECKS.md) records their universes, proof
dependencies and external-source limits. No CI success or historical priority
is inferred from these local checks.

Remaining targets are prime-prefix counting/ranking, full graph transporters
beyond the specified bank actions, and further lawful observation interfaces.
General factoring speedups, unrestricted tournament classification and Collatz
convergence are not established here.
