# Algorithm sources and a concrete precision test

**CITED source claims / PROVED local precision example / OPEN transfer to
general strong-atom composition.** Sources inspected 2026-10-29. This note
does not certify the full correctness of a new complexity breakthrough.

## What the published algorithm claims actually measure

[Big O Speedrun's matrix problem](https://bigospeedrun.com/matrix-multiplication)
uses division-free arithmetic circuits over C, with constants free. The
listed 9/4+epsilon exponent and verification badge concern that model;
they are not automatically bit-operation bounds for integer matrices.

[Its 3-SAT problem](https://bigospeedrun.com/3sat) uses a worst-case bit-time
target with a specified success guarantee. A positive subexponential saving
of the stated kind would contradict ETH. A label on the site is evidence
about the submitted claim, not our own proof audit.

The pinned primary manuscript
[The Exponential Time Hypothesis is False, commit 7ef5c0cc](https://github.com/openai/math/tree/7ef5c0cc5ca4b0e072b7a050fca09ab24bc6fa3a/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026)
claims deterministic exact counting for fixed-width CNF in subexponential
bit time. The proof route combines sparse parity constraints, interpolation,
and tensor contraction. Its claimed exact counting also supplies conditional
counts for self-reduction and exact sampling. The source expressly keeps
size blowup in reductions separate from the original variable count.

## The useful implementation mechanism

The inspected [contraction appendix](https://github.com/openai/math/blob/7ef5c0cc5ca4b0e072b7a050fca09ab24bc6fa3a/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex)
retains a tensor's creation precision when lifting entries. Its Cauchy
deformation groups ports, then uses local interpolation and a valuation
budget to truncate crossing modes. Tapered Hermite multiplicities control
precision loss at different dyadic depths. This is the portion inspected
in detail, through the crossing-mode expansion; the complete recursion and
complexity theorem are not audited here.

These mechanisms suggest a useful question for our carriers: can the target
be evaluated from a bounded precision budget without reconstructing every
hidden state? A large minimal faithful state space does not by itself prove
a matching time lower bound for one scalar target. Conversely a fast scalar
contraction need not preserve enough data for the next attachment operation.

The exact [precision certificate](PRECISION.md) gives a small hostile to
naive extrapolation. Our quadratic incidence matrix has binary rank five,
but its Smith factors are 1^5,2^4,4. The missing directions become observable
at higher powers of two. Thus an implementation that lifts only a five-row
parity summary cannot silently treat it as the full integer response.
Canonical residue lifts are not ring homomorphisms: already lifting
1+1=0 from F2 gives 0, whereas adding the two lifted ones gives 2.

## A composition principle with a testable boundary

The same object can support three distinct reductions:

1. A realization equation recovers an omitted coordinate for each genuine
   object, as t^2=ad+bc does for our symplectic carrier.
2. A linear relation removes a redundant observation at a fixed coefficient
   field, as the third split-prime elliptic character follows from two.
3. A precision budget postpones directions whose integer invariant factors
   have high valuation, as the Smith filtration records exactly.

Each reduction needs a different closure proof under the next operation.
The integer return controller adds another kind: a source-sensitive height
and queue potential certifies termination, while endpoint descent still
needs a separate guard. These are compatible parts of a composition theory,
provided their source, operation and target remain explicit.

The immediate algorithmic experiment should use the actual native partner
path tensor from [ROOKS R1--R3](../charged_recurrence_memory_20261028/ROOKS.md),
keeping endpoint injectivity and cycle exclusion. Compare exact boundary
response matrices over F2, modulo 4 and over Q; then test whether a
precision-aware contraction evaluates the same target faster. A collapsed
cycle count or a separate-port histogram is not an adequate control.

## Next source checks

Continue with the remaining contraction recursion and its input tensor
access model, then the matrix-multiplication manuscript's variable-separation
conditions. Only after those checks should a complexity consequence be
imported into the tournament enumeration lane. Exact output size and the
need to name all isomorphism classes remain separate costs from counting
A000568 by a recurrence.
