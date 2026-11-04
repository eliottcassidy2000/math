# Independent audit of the symmetric seven prime classification

**PASS after one small empty-exterior correction.** This review was performed
by the sibling counting agent, independently of the bridge producer. The
full prime/7p existence problem remains OPEN.

## Proof review

The following implications and boundaries were checked directly.

* The Hamilton-frame fibre has size H/Aut, and its class-count identity term
  is B/H. Session142's inherited Aut restrictions justify both the prime
  division formula and the C7 correction. The unrestricted strong-to-padded
  series uses one nontrivial SCC; it is valid because a second factor would
  force the forbidden count7.
* A unique nontrivial7-orbit is a uniform module, characteristic under
  isomorphisms. The disjoint path-product bound excludes two such orbits
  below19467. Every inner automorphism extends over a uniform module, so
  the189-path module with Aut21 is excluded. The six175-path circulants
  form one class; the producer now checks the unit multipliers explicitly.
* Maximal-run contraction gives the exact ordered-cover coefficient law.
  Inserting a new exterior vertex preserves the visit number even when the
  named vertex occurs repeatedly. Deleting that new vertex recovers the
  source word, proving injectivity coefficient by coefficient.
* The maximal proper strong-set proof retains an arbitrary named root.
  When one exterior vertex remains the desired deletion is already found;
  otherwise uniformity and a two-vertex exterior bridge give the claimed
  strong deletion. Repetition reaches quotient order7. The longest-cycle
  proof also supplies a Hamilton path starting at the named root, so the
  fixed-spine universe covers every rooted strong quotient class.
* The complete order7 Pareto rows imply the weighted minimum. For every
  module of order at least2, C2>=2H, making27H+9C2 the smaller row. Thus
  every larger strong quotient has H(Q[M])>=9(3H(M)+C2(M)). The strict
  cutoff19467 is used correctly; equality is not accidentally excluded.
* Rooted quotient automorphisms would extend independently of the C7 inner
  action, so the filtered rooted quotients are asymmetric. This justifies
  both reciprocal root-start weights and division of labelled counts by
  (q-1)!. The unique7-orbit makes the quotient-class to full-class map
  injective and surjective in the classified family.
* The proposed all-order two-visit inequality is explicitly OPEN. Its
  conditional weighted consequence is algebraically correct and is not
  used in the proved five-fibre inventory.

The subset kernel recursion was reviewed: a return to root preserves the
exterior mask but changes the last vertex to root; the next transition must
remove an exterior vertex. Thus the recursion is acyclic and counts exactly
the permitted repeated-root words.

## Boundary correction

The note defines A0=H(Q-v), with H(empty)=1. The initial helpers returned
A0=0 at the one-vertex quotient, because their word loops omitted the empty
word. Both helpers now return the conventional empty-exterior value, with
explicit one- and two-vertex controls. Every production quotient has order
at least3, so the classification and all numerical inventory rows are
unchanged. The correction affects the helper's boundary contract, not the
all-order proof.

## Independent computation

[counting/bridge_audit.py](../counting/bridge_audit.py) imports no mathematical
function from the producer. For each permutation of exterior vertices, it
classifies each possible root gap as forced, optional, or forbidden. The
contribution is x^forced*(1+x)^optional. This polynomial method independently
computes the whole visit kernel.

The[receipt](../counting/bridge_audit.json) covers:

* All959 strong rooted Hamilton frames through quotient order6. Literal
  root-preserving canonicalization reproduces the five full counts4151,
  5789,14791,16079,17717 with multiplicities2,2,2,4,2 respectively.
* All six displayed Pareto source kernels through quotient order7, including
  both incomparable order7 rows. This is a source check, not a second full
  order7 dominance census.
* All five actual full tournament witnesses, independently counted by a
  predecessor-mask Hamilton DP, together with their literal7-cycle
  automorphism.

The producer supplies the exhaustive32768-frame order7 dominance scan; this
review checks its complete-universe argument, recursion and two extremal
sources. It does not mislabel the targeted audit as a second exhaustive
order7 run. The producer's normal and optimized receipts, and the independent
normal and optimized audit receipts, are compared separately.
