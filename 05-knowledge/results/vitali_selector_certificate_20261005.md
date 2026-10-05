# Quotient selectors, finite certificates, and three different stopping questions

2026-10-05. **PROVED:** the finite factorization criterion; the computability
classification of least representatives of a computably enumerable
equivalence relation; the halting-pair obstruction; and the distinction
between shortest explicit certificates and shortest universal programs.
**FINITE-EXACT:** the tournament witnesses and checked component banks below.
**CITED:** the stated connections to Reimann's survey. **OPEN:** universal
Collatz completion and the decidability of its full common-future relation.
No undecidability assertion about Collatz is made.

[Program](../../04-computation/experiments/vitali_selector_certificate_20261005.py)
and [output](vitali_selector_certificate_20261005.out). This package also
documents the scoped repair to the historical
[Vitali sample script](../../04-computation/vitali_nonmeasurable.py).

## 1. Recovery: what the repository's Vitali language actually records

The closest older routes are
[THM-168, lambda-completeness](../../01-canon/theorems/THM-168-lambda-completeness.md),
[THM-169, vitali-atom-characterization](../../01-canon/theorems/THM-169-vitali-atom-characterization.md),
and [THM-172, c5-lambda-determined](../../01-canon/theorems/THM-172-c5-lambda-determined.md).
They concern a finite tournament observation: the labelled pair-coverage
array lambda_ij counts directed-triangle vertex sets containing i and j.
The older documents mix finite verification, sampled characterization,
and suggested analogies; this note does not promote the sampled all-case
claims. In particular it does not use THM-169's proposed full
characterization or any all-order sufficiency of pairwise overlap data.

The exact mechanism is the following elementary theorem. If X is a finite
set, pi:X->Y an observation, and f:X->Z an observable, then

    f=g∘pi for some g on pi(X)
    iff f is constant on every pi-fibre.                         (1)

The reverse implication defines g(y) using the common value on the fibre.
Equivalently, for discrete Z, f is measurable for the coarse sigma-algebra
{pi^(-1)(A): A subset pi(X)} iff (1) holds. Every f remains measurable for
the full power-set sigma-algebra on X. No choice axiom or non-Lebesgue-
measurable set is involved. A least binary tournament mask is a computable
fibre representative, but its H value cannot stand in for every member
when H varies inside that fibre.

**Independent exact hostile.** Encode labelled tournaments by the bits on
lexicographically ordered pairs i<j, with bit1 meaning i->j. At order7,
the masks

    1354195 and 1351956

differ exactly by reversing the six internal arcs on {0,1,2,3}. Their
labelled lambda arrays are identical, while their Hamiltonian-path counts
are109 and111. Their simple directed7-cycle counts are8 and9. Both path
counts are checked by subset dynamic programming and direct permutation
enumeration. This proves the needed non-factorization without using any
statistical threshold or unproved characterization.

Our board is **observed fibre; global equivalence; canonical label;
finite witness; stopping oracle; source identity**. Anchor: retain the
certificate when selecting a Collatz component. Niche: recover the exact
finite meaning of the earlier Vitali analogy. Wildcard: identify precisely
which shortest-description problem is uncomputable.

The recent proof-carrying mechanisms are
[finite seed receipts](finite_seed_receipt_compiler_20261005.md) and
[fair frontier extension](fair_frontier_extension_20261004.md): observed
joins are independently checked, and a component is grounded only when
a retained certificate connects it to ROOT. These supply the finite
positive counterpart of the computability statements below.

## 2. Classical Vitali, discrete measurability, and computability

The classical construction chooses one representative from each equivalence
class x~y iff x-y is rational, among classes meeting [0,1]. With a choice
principle supplying such a transversal V, its translates V+q for rational
q in[-1,1] are disjoint, cover[0,1], and lie in[-1,2]. If V were Lebesgue
measurable, measure zero would contradict the covered interval, while
positive measure would give infinite measure inside a bounded interval.
This is an obstruction to Lebesgue measurability on an uncountable space.

An equivalence relation E on positive integers has a different situation:

    mu_E(n)=min{m>=1: m E n}                                    (2)

exists by well-ordering, without selecting from uncountably many unordered
sets. In the discrete topology every subset of the domain is Borel, so
every such map is measurable. This says nothing about an algorithm for (2).
Changing the topology or sigma-algebra would require a separately stated
problem; neither classical Vitali nor discrete measurability resolves
Collatz stopping.

**The second repository Vitali lane: covering versus boundary witnesses.**
[THM-398, lrc-reduction-to-Cprime-and-dominance-dodge](../../01-canon/theorems/THM-398-lrc-reduction-to-Cprime-and-dominance-dodge.md)
and [HYP-2104, lrc-vitali-handoff](../hypotheses/HYP-2104-lrc-vitali-handoff.md)
used a different analogy, concerning periodic danger arcs. Its sound core
is elementary: a connected interval covered by separated danger arcs must
fit in one arc. A component longer than an arc cannot be covered. The
fixed positive radius does not give a fine Vitali cover at arbitrary small
scales, so a Vitali covering theorem cannot be invoked from bounded
eccentricity alone. The corrected notes now retain this distinction.

The same audit exposed a separate scaling obstruction to the old
unrestricted conjecture C′, "a speed divisible by n implies M(S)>1/n."
For n=3, S={3,6}, one has M(S)=1/3, since M(dS)=M(S) under the surjective
circle map t->dt and M({1,2})=1/3. Already n=2, S={2} gives M=1/2.
The repaired *primitive* candidate gcd(S)=1 remains OPEN. Its sufficient
reduction to LRC first divides all speeds by their gcd; the elementary
no-multiple clock witness and the n>=3 dominance criterion survive.
The endpoint arc test must retain one common lifted integer centre and
ordinary absolute distances, not separate circle distances which erase
that centre. Those clauses are repaired with explicit correction lineage.

For the strict safe set G={t:min_v||vt||>1/n}, positive Lebesgue measure
is equivalent to nonemptiness. At the weak threshold >=1/n, a null safe
set can still be nonempty: S={1,2}, n=3 has exactly the two witnesses
1/3 and2/3. Null measure alone cannot distinguish this tight set from an
empty one. The current [atomic input ledger](atomic_prefix_certificate_20261005.md)
instead assigns every integer a positive atom, so an exact global zero
residual excludes each integer exception. This is an explicit change of
space and measure, not a Vitali theorem or a conclusion from a finite
sample. A small residual requires the stated minimum-atom bound to certify
a specified finite interval.

## 3. Certified minima stabilize, but need not announce stabilization

Call E *computably enumerable* if an algorithm enumerates exactly its true
pairs. The enumeration supplies positive witnesses; it need not decide
that a pair is absent.

**Theorem.** For every computably enumerable equivalence relation E on
the positive integers:

1. mu_E is the pointwise limit of a computable nonincreasing sequence of
   positive integers mu_s(n), with at most n-1 strict changes on input n.
2. mu_E is computable **iff** E is decidable.
3. More generally, a computable class-constant selector r satisfying
   r(n) E n exists **iff** E is decidable.
4. A halting oracle computes mu_E uniformly from an enumerator for E.

**Proof.** At stage s take the finite pairs enumerated so far, add their
undirected edges, and form finite connected components. Unseen points are
singletons. Let mu_s(n) be the smallest point in n's current component.
Every edge is E-valid, so its transitive closure remains a subset of E;
hence mu_E(n)<=mu_s(n)<=n. More edges can only lower the minimum. The
pair(n,mu_E(n)) eventually appears, so the limit is exactly mu_E(n).
The descending positive-integer sequence has at most n-1 strict changes.

If E is decidable, query the finitely many m<=n and take the least true
one. Conversely n E m iff mu_E(n)=mu_E(m). The same converse holds for
any class-constant representative selector, proving the third statement.
Finally, for each pair(m,n) a program can halt exactly when its E-witness
is enumerated. A halting oracle decides these finitely many queries.

The first assertion is unconditional for c.e. E. It must not be confused
with the stronger **computable iff decidable** assertion in the second.
A computable number of possible changes is not a computable upper bound
on when the last change occurs.

**One-change hostile.** Fix an effective enumeration of programs. Partition
the positive integers into pairs {2e+1,2e+2}. Merge a pair iff program e
halts on its designated input; otherwise leave its two points separate.
This defines one c.e. equivalence relation E_H, all of whose classes have
size at most2. Yet

    mu_EH(2e+2)=2e+1 iff program e halts.

A computable selector would decide the halting problem. For completeness,
if a total halting decider existed, the program which on input e loops
exactly when the decider predicts that e halts on input e would contradict
the prediction on its own index. Thus no selector algorithm follows from
the mere existence of a least representative, measurable selector,
eventually stable approximation, or even a one-change bound. This is an
abstract hostile relation; it is not a reduction from halting to Collatz.

## 4. What these statements say for actual Collatz receipts

On positive odd integers put

    n E_C m iff U^i(n)=U^j(m) for some finite i,j>=0.

This is an equivalence relation because U is deterministic; two shared
futures can be advanced until their times on the middle trajectory agree.
It is c.e.: dovetail the finite trajectories and enumerate witnessed
intersections. It has the least-representative approximation of section3.
No assumption of convergence is used to enumerate valid pairs.

Since U(1)=1,

    mu_EC(n)=1 iff n has a finite actual route to1.              (3)

ROOT is therefore one distinguished class, with its unique least label1.
Finite checked component graphs can safely recognize any witnessed part
of that class. Their label1 must come with an actual route or checked join
chain; the code reconstructs and independently replays one. A nonroot
label in a finite bank is only a current minimum, not a certificate that
its full class excludes1. A label greater than1 does not prove divergence.

The full quotient cannot be used to bypass the open step. A computable
selector for E_C would be equivalent to a decision algorithm for E_C,
not to an already proved assertion that every value of that selector is1.
If Collatz holds then E_C is the single class and the constant selector1
is computable; there is no contradiction with the general hostile above.

An oracle-selected true relation may still be used as a suggestion:
dovetail a search for its finite witness and accept only the witness when
it appears. The witness is independently checkable without the oracle.
There is no promise of termination for a false suggestion. This is a
precise useful interface for a heuristic or oracle, preserving source
identity and the proof obligation rather than trusting a class label.

## 5. Shortest explicit proof versus shortest universal program

Let V(n,c) be any decidable verifier of finite certificate strings, with
only finitely many strings of each bounded length. Searching strings in
shortlex order and applying V computes the shortest valid certificate
whenever one exists. The search is partial: if no certificate exists it
may run forever. For explicit Collatz first-hit valuation words, each
check is finite integer arithmetic. Such a search is total on all
positive odd inputs exactly when every input has a valid ROOT certificate.
Literal deterministic trajectories have a unique first-hit valuation
word, but richer common-future proof syntax can have many certificates;
the same search statement applies to either syntax.

The test script uses the self-delimiting syntax

    code(a1,...,ar)=1^a1 0 ... 1^ar 0 0,   code(empty)=0.

It is prefix-free: each valuation block begins1, while the final marker0
appears where a next block would begin. The code retains both the complete
finite word and its stopping marker. Prefix-freeness makes concatenations
decodable; it does not prove that every source has an accepted certificate.

For an optimal universal prefix-free machine M, a computable selector of
a shortest program for every output x would instead compute K_M(x).
That is impossible: given k, computably choose the first string with
K_M(x)>k. One exists by counting the at most2^(k+1)-1 programs of length
at most k. The chosen string would then have a description consisting
of this fixed algorithm and a self-delimiting encoding of k, of length
O(log k), contradicting K_M(x)>k for sufficiently large k.

The failed analogy is precise. Checking whether an arbitrary program
eventually outputs x is not decidable finite-certificate verification.
Adding a terminating execution trace makes verification decidable, but
changes the object and potentially its length. Retain that trace or a
checked runtime bound when transporting program compression into proof
compression. Neither theorem gives a useful general efficiency bound.

Reimann's [Information vs Dimension: An Algorithmic Perspective,
arXiv2408.05121v1](https://arxiv.org/pdf/2408.05121v1) treats prefix-free
coding in section2 and oracle-relative effective dimension in section4;
Theorem4.15 states the point-to-set principle using a minimizing oracle.
That oracle is part of the information supplied to the description
machine, not an automatically computable certificate selector. Applying
the principle to an infinite Collatz coding space requires specifying
that space, metric and measure; a finite graph quotient alone provides
none of those dimension hypotheses.

## 6. Historical correction and exact scope of the experiments

The historical `vitali_nonmeasurable.py` computed trace(A^7)/7 under the
name c7. This statistic is not the number of simple directed7-cycles.
The four-vertex tournament with arcs

    0->1, 0->2, 1->2, 1->3, 2->3, 3->0

has trace(A^7)=14 but no seven-vertex simple cycle. Padding by a transitive
three-vertex block preserves this witness at the script's order7. At the
prime lengths3 and5 used there, every closed walk in a loopless tournament
is simple: a nonsimple one would split into cycles of length at least3,
requiring total length at least6. At length7 a3+4 concatenation is possible.

The scoped repair renames the historical quantity `walk7`, restricts the
helper to its prime lengths3,5,7, and explains the finite coarse
measurability convention. It also corrects the conditional-variance
normalizer: singleton fibres contribute zero variance but their samples
must remain in the total denominator. The toy groups {0,2} and {100}
give within-fibre variance2/3, not1. The empirical entropy calculations
are distinct from the displayed count of observed fibres; that count's
logarithm is not itself an entropy estimate.

The historical20,000-draw seed42 sample was rerun after repair. Its
19,242 observed fibres and five walk7-varying fibres among624 repeated
fibres are retained; the within-fibre variance now prints0.0001 instead
of0.0018. These rounded floating sample summaries are **not** an exhaustive
theorem, a simple-cycle census, or a measure-theoretic Vitali result. The
old saved historical output is provenance; the corrected reproduction
summary appears in this package's output.

Run from the repository root:

    python -B 04-computation/experiments/vitali_selector_certificate_20261005.py
    python -B -O 04-computation/experiments/vitali_selector_certificate_20261005.py

The new exact controls use the standard library. Reproducing the historical
sample additionally uses its existing NumPy dependency. Universes:

* every labelled tournament of orders3,4,5:8+64+1024 objects, with lambda
  fibres2,11,208 and agreement of two Hamiltonian-path counters;
* the fixed order7 reversal pair, independently counted simple cycles,
  and the order4/order7 trace-versus-simple-cycle hostiles;
* the nonprimitive C′ counterexamples, and82 exact piecewise-linear
  scale-invariance controls from subsets of{1,...,6} of sizes1,2,3 scaled
  by2 and3; the maximizer search includes all affine-piece intersections;
* all64 positive odd starts below128, at depths0,1,2,4,8,16,32,64,
  with finite graph closure and512 independently replayed representative
  receipts. All64 are rooted by depth32 of this shared bank; this is a
  finite result, not a global stabilization bound;
* a finite delayed-pair illustration with its final change at stage20;
  the general no-selector conclusion uses the halting reduction, not this
  illustration;
* the explicit local sibling BREAK at x=7, with state (17,1), sign+ and
  exponent4. The literal ROOT words Q21=(6), Q3=(1,4), and
  Q7=(1,1,2,3,4) independently give
  boundary(Q21-Q3-Q7)=[21]+[1]-[3]-[7]. Thus failure of that local
  automaton is not a global obstruction to actual-edge repair; this
  corrects the unscoped converse in the incoming
  `collatz_two_sheet_receipts_20261005.md` section5 (commit9eb19b215),
  without disputing the sufficient local-success certificate;
* all2584 certificate codes of length at most18, exact prefix-free parsing,
  and bounded shortlex searches for odd inputs below32. Missing certificates
  within that budget are reported as absent from that finite search only.

## 7. Transfer ledger

| Source -> target | Preserved predicate and map | Lost information / required sidecar |
|---|---|---|
| Tournament -> labelled lambda fibre | Equality of directed-triangle pair coverage | H and simple c7 can vary; retain a fibre member or enough additional observables |
| Finite checked graph -> component minimum | A witnessed common-future relation | Retain the edge receipts and original source; the label is not a global negative certificate |
| c.e. E -> limiting least representative | True equivalence and minimum | Limit time is not supplied; a halting oracle adds information rather than proving convergence |
| Finite certificate -> prefix-free code | Lossless syntax and decidable verification | Coverage/termination of certificate search remains separate |
| Program -> output | Partial computable execution | Shortest-program selection cannot discard the unknown halting behavior |

The constructive survivor is an anytime proof-carrying component selector:
every accepted improvement has a finite witness, and ROOT labels are
grounded immediately when their witnesses arrive. Its safe finite steps
do not require a computable promise that no later improvement exists.
