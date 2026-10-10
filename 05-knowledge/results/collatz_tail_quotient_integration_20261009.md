# What parity conjugacy preserves, and what the Collatz target adds

2026-10-09. **CITED + PROVED HERE:** the classical parity conjugacy, its
eventual-tail equivalence interpretation, and the elementary measurable
selector obstruction below. No priority is claimed for these classical
facts. **PROVED:** the marked-cycle image `1/(4-a)` and the resulting
failure to preserve the positive-integer target. **FINITE-EXACT:** the
declared controls. **OPEN:** universal Collatz coverage.

[Program](../../04-computation/experiments/collatz_tail_quotient_integration_20261009.py)
and [output](collatz_tail_quotient_integration_20261009.out).

## Inheritance and the object being represented

The closest mechanism is exact parity-cylinder coding, also used in
[THM-4476, thin-divergent-orbits-reciprocal-sums-finite](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md).
The classical primary source is Bernstein and Lagarias,
[The 3x+1 Conjugacy Map](https://websites.umich.edu/~lagarias/doc/bernstein.pdf),
Introduction and Appendix B. Its conjugacy is named Phi; our Q is its inverse.
We use only this conjugacy, not its permutation-cycle estimates.

The corrected near miss is the
[October5 Vitali-selector note](vitali_selector_certificate_20261005.md):
a finite observation fibre or a discrete integer equivalence relation is
not a non-Lebesgue-measurable set. The hostile here is a genuine change of
ambient space: positive integers have Haar measure zero in Z_2. The
least-used sidecar is the embedding of the original integer domain and its
distinguished terminal inside the ambient dynamical system.

Anchor: source-specific ROOT certificates. Niche: actual versus symbolic
coalescence. Wildcard: the precise infinite-measure-space version of the
earlier Vitali analogy. Board: **source embedding / finite word / clock /
equivalence class / measure / distinguished terminal**.

## 1. All odd multipliers have the same unmarked parity model

For any odd integer a, put

    T_a(x) = x/2 (x even), (a*x+1)/2 (x odd), on Z_2,
    Q_a(x)_j = T_a^j(x) mod2.

For a binary word beta of length m, let r count its ones. Iterating the
specified affine branches gives

    T_a^m(x) = (a^r*x+b)/2^m.

Each word is legal on exactly one residue modulo2^m. One proof reverses
the branches `y -> 2y` and `y -> (2y-1)/a`; division by odd a is a unit
operation in Z_2, and the prescribed branch fixes the parity. Thus each
length-m parity cylinder has Haar mass2^-m. Compatible finite inverses
give a measure-preserving homeomorphism Q_a to binary sequence space, and

    Q_a T_a = shift Q_a.                                      (1)

This identifies a dynamical system with a shift, but does not identify its
ordinary positive integers or its specified ROOT.

Indeed, `Q_3(1)=101010...`. The unique T_a point with that itinerary is

    x_a = 1/(4-a),           x_a -> 2*x_a -> x_a.               (2)

The denominator is odd, so these are legitimate Z_2 points. Substitution
proves both arrows and uniqueness follows from the cylinder bijection.
Consequently the conjugacy `Q_a^-1 Q_3` sends1 to `1/(4-a)`.
For a=5 it sends the positive cycle `(1,2)` to the negative cycle `(-1,-2)`.
The positive T_5 cycle through1 instead is `(1,3,8,4,2)`.

**Mechanism:** an unmarked conjugacy preserves orbit equations while losing
the embedding of N_+ and the terminal name. A transfer of the integer
Collatz assertion needs those extra predicates preserved. The exact
cycle calculation is the cheapest hostile to transferring integer
convergence from the shared shift model alone.

## 2. Synchronous coalescence really is eventual-tail equivalence

For T=T_3 define

    x E_sync y iff some n>=0 has T^n(x)=T^n(y).

By (1), this is precisely

    Q_3(x)_j=Q_3(y)_j for all sufficiently large j.              (3)

The binary relation in (3) is usually called E_0. Its classes are the
orbits of the countable group of finite bit flips. The action is free and
preserves the fair product measure. Each class is countable and null.

There is no Haar-measurable set containing exactly one representative
from each E_sync class. If V were such a set, its translates by all
finite flips, transported through Q_3, would be disjoint, measurable,
equal in measure, and cover Z_2. Measure(V)=0 contradicts the cover;
positive measure contradicts finite total mass. This is a genuine
Vitali-style argument in an explicitly specified infinite probability
space, unlike the earlier finite-selector analogy.

There is also no Haar-measurable real-valued complete invariant for this
relation. To prove this without a classification theorem, let A be a
measurable invariant set. Invariance under flips of the first m bits
makes its intersection with every length-m cylinder have measure
`2^-m*measure(A)`. Hence A is independent of every cylinder, then of the
generated completed sigma-algebra, including itself. Thus its measure
is0 or1. Apply this to the rational sublevel sets of an invariant real
function: the function is constant almost everywhere. A complete
invariant would then put a conull set in one countable, null class,
which is impossible.

For asynchronous common future, allow `T^i(x)=T^j(y)` with i and j
different. This is a different equivalence relation. Its classes are
still countable and contain E_sync classes; the invariant-function
argument still excludes a measurable real complete invariant. We do not
identify asynchronous coalescence with E_0 in (3).

## 3. What this does and does not say about integer certificates

The T_3 basin of the cycle `(1,2)` is countable: it is a countable union
of finite preimage sets. It has Haar measure zero and comprises exactly
two synchronous classes, the two eventual alternating phases. Its
membership predicate is measurable. A finite path to1 is still a fully
valid witness; the obstruction concerns a *complete classifier of all
ambient classes*, not a classifier of this one chosen basin.

The integer conjecture becomes

    Q_3(N_+) is contained in the eventually alternating sequences.

Neither the shift model nor a Haar-almost-everywhere conclusion gives
that containment. Conversely, the measurable-selector obstruction gives
no integer undecidability claim: on N_+ the least representative exists
by well-ordering and every map is measurable in the discrete sigma-algebra.
If Collatz is true, the entire positive domain has just the two synchronous
classes represented by1 and2. No contradiction with the ambient no-go arises.

The [atomic source ledger](atomic_prefix_certificate_20261005.md) supplies
a different, target-sensitive measure: put positive mass at every integer.
An exactly zero unresolved mass then rules out every integer exception.
There is no conversion from a Haar-null exception set to zero atomic
mass without additional source-specific information.

## 4. Reuse without erasing the question

The useful representation is a pointed system `(space,T,source domain,ROOT)`
together with finite witnessed arrows. A quotient can compress a verified
finite certificate graph, but its class label alone does not give a ROOT
arrow. The next primitive is therefore an actual common-endpoint receipt
carrying its original sources, its two clocks, and its terminal witness.
This agrees with the independent
[fibre/clock packet](collatz_fibre_integration_20261009.md),
[specialized-join packet](collatz_debt_observer_integration_20261009.md), and
[cycle-observer packet](collatz_basin_observer_20261009.md).

The connection ledger is explicit:

| Source -> target | Map / preserved predicate | Lost data / required sidecar | Decisive test |
|---|---|---|---|
| Z_2 dynamics -> binary shift | Q_a; time and orbit equations | N_+ and chosen ROOT; marked embedding | Equation(2), a=5 |
| Synchronous trajectories -> tail classes | Q_3; equal-time meeting | finite meeting witness; retain words | Decode both prefixes |
| All ambient classes -> real label | proposed complete invariant | impossible with Haar measurability | finite-flip proof |
| Ambient measure -> integer coverage | proposed implication | positive mass on each original source | Haar-null N_+ hostile |

## 5. Exact controls and scope

Run from the repository root:

    python -B 04-computation/experiments/collatz_tail_quotient_integration_20261009.py
    python -B -O 04-computation/experiments/collatz_tail_quotient_integration_20261009.py

Universe: all40,950 residues for a in `{1,3,5,7,11}` and1<=m<=12.
Independent inverse paths use an affine numerator and reversed branches.
The finite-flip control has64 classes of size16 at width10; finite
selectors exist there, so it is not presented as a finite proof of the
infinite obstruction. Rational and modular checks agree on (2), and the
positive T_5 cycle is evaluated directly. All166,956 checks pass.
The infinite assertions above have proofs; the finite controls do not
extend their universe by extrapolation. No new ROOT coverage is claimed.

Independent peer proof/API review: PASS. Normal and optimized executions
match the saved LF output; SHA256
`29e4352a6a78e0a163428610adb8edc14ebd7afc049fded114179611e9ccef9e`.
