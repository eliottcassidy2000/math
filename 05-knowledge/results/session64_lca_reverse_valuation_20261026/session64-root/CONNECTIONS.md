# Three different ancestor structures, and the summaries that survive them

Session64, 2026-10-26. This is a scoped synthesis, not a new universal
classification or a Collatz/nonhalting theorem. The operative proposal is
**continuation-exact summaries with explicit guards**. Prior Session56 already
developed this viewpoint; the present work supplies native applications.
Session63's proposed "integral observers of decorated constructions" likewise
separated an integral carrier, its operation bank and its realizable states.
Here that proposal becomes testable: a coding may preserve a 2-adic carrier
but lose integer realization, while a Newton observation kernel is useful
precisely because the actual integral transfer preserves it.

## 1. Recovering, rather than conflating, the LCA threads

The prior `session56_rmq_lca_connections_20261025/CONNECTIONS.md` supplies
three distinct constructions:

* Card6: least-significant-digit trie LCA has depth `v_p(x-y)`. Sorting by
  those digit addresses reduces pairwise depth queries to adjacent RMQ, with
  labels and ordering retained.
* Card8: the lowest strong tournament module locates two vertices' separating
  children. Their child ports and the quotient tournament determine the arc;
  the module-tree shape alone does not.
* Card12: an ordered update word has a total increment and a minimum prefix.
  This is an associative range summary, but overlapping blocks are unsafe.

The older `collatz_bugs_20260925.md` uses a fourth, dynamical ancestry notion:
actual forward coalescence in the known Collatz root basin, after removal of
the root cycle. Universal membership in that basin is not assumed.

These objects are not interchangeable merely because each admits tree
language. A shared low-bit prefix is not an eventual common future, and a
module-tree LCA is not a legal path-move certificate.

## 2. Exact transport: valuation coding preserves the p-adic trie

On odd 2-adics outside the respective countable numerator-zero preimage
sets, ordinary Syracuse and twisted Syracuse have the same positive-integer
valuation-word shift. If two words first differ at p versus q after common
cumulative valuation A, the corresponding points have difference valuation
`A+min(p,q)` in either coding. Thus the same-code conjugacy preserves the
actual 2-adic metric and its prefix hierarchy, not merely an abstract tree.

Lost coordinates are real order, positive integrality, the integer inverse
branch structure, and Antihydra's resource. The all-2 word gives ordinary
Collatz point1 but twisted point3/5. The all-1 word gives ordinary point-1
but twisted point1. Positive Collatz integers3 and13 both map to5, whereas
a positive twisted odd integer has at most one positive odd predecessor.
There is therefore no induced isomorphism of their positive-integer ancestor
trees. The general shift conjugacy does preserve coalescence on its stated
2-adic domains; restricting to positive integers is the invalid step.

## 3. Exact transport: Antihydra's guard is the prior debt monoid

For a finite valuation word w=(k_1,...,k_m), use updates `d_i=2k_i-3`.
Let Delta be their sum and let mu be the minimum of all prefix sums,
including the empty and complete prefix. The exact positive resource
requirement is `a>=1-mu`. Concatenation obeys

    (Delta,mu)*(E,nu)=(Delta+E,min(mu,Delta+nu)).

This is the additive specialization of Session56's guarded affine-debt
summary, not a newly discovered monoid. Its reverse form is

    R(empty)=1,   R(k w)=max(1,R(w)+3-2k).

It answers finite guarded-word queries and supports disjoint range
composition. It cannot establish that a chosen integer has that valuation
word, or that every prefix of its infinite word is safe. The endpoint-only
hostile k=(1,2) has total increment0 but requires a>=2; the reversed word
(2,1) requires only a>=1. Hence net drift alone discards a decisive guard.

The prefix-inclusive convention means the resource at the final odd section
also remains positive, ready for another odd step. If only execution of a
finite number of steps is requested, omit the final continuation guard
explicitly rather than silently changing the threshold.

## 4. Exact transport: the reverse valuation is a finite-DAG gauge

The ordinary A000568 graph is finite and stage decreasing. An edge v->w
carries an odd integer factor T and a nonnegative binary exponent J. The
minimum suffix cost V gives reduced cost `J'=J+V(w)-V(v)>=0`. With
`B=V(initial)`, the reachable coordinate change
`Z(v)=2^(V(v)-B)Y(v)` is integral path by path. Because V belongs to the
destination state, actual sums into the same state commute with this gauge.

This is a real connection to guarded continuation bookkeeping: both methods
retain enough suffix/prefix data to justify a transformation. It is NOT a
transfer of an infinite-orbit termination argument. The resource summary
composes a specified word using prefix extrema; the counting potential
minimizes across alternative suffixes in a finite graph.

The forward minimum reduced prefix cost P adds an exact arithmetic test:
an edge lies on a complete path of reduced cost<K iff `P(v)+J'<K`.
This reduces modular arithmetic updates after constructing the full graph.
It does not merge structural states or remove preprocessing cost.

Equal potential values alone never justify a merge. Two states can both
have suffix valuation0 while their only terminal edges have odd factors1
and3. Their minima coincide but their weighted continuations differ.
The tropical shadow preserves a minimum, not multiplicity or the final sum.

## 5. Exact ports precede any clone quotient

The repaired endpoint-four ledgers work on literal labelled tournaments.
Each step retains its actual Hamiltonian word, flipped pairs, blackness
witness, endpoints and clone-block intersections. The cyclic-shift bridge
gives legal moves and an exact XOR target identity, but blackness must be
checked in every successive source. A quotient with the right endpoint
degrees or tree shape cannot supply those tests.

Session56 card4's correlated gluing explains the same danger abstractly:
reusing one shared choice is not two independent choices. At a shared port,
retain the literal owner and correlations. This is methodological transport,
not a claim that a join-tree theorem proves the native black-core problem.

The twenty order12/14 exceptional paths are clone-free at all eighty states.
This is an observed consequence, not a restriction on the original source
families. The subsequent order11 closure has an actual repair that forms
the ordered clone pair (5,4) after its second move and splits it at the
target. Every arc and block intersection is retained throughout. Thus clone
transport is tested in the native repair, not only on an external positive
control. No contraction or quotient lifting is used in the repair proof.

There is also an all-size sufficient symmetry gate behind most certificates.
For a closed blue P, let k=d_T(z), B=K_P(T), and choose a directed Hamiltonian
word W of B minus vertex0, beginning at1 and ending atz. The three bridge
words are (1,...,z,0), (W,0), (0,W). If d_T(1) differs from k, the first
two words are black by their endpoint-degree sums. If the third were blue,
its reversal composed with P's reversal would be an automorphism of B,
fixing0,z and preserving the labelled port set I union {1}. Thus failure
to preserve that set, or even order of the composite permutation, proves
the third word black. A tournament automorphism cannot have even order:
its involution power would swap a pair while preserving its directed arc.

This transports a symmetry constraint into a native move certificate, not
merely a numerical analogy. It is sufficient, not necessary: the order11
clone-merger source has d_T(1)=k=4 and still admits valid repairs. No
existence theorem for the required endpoint-fixed W follows from the gate.

The counting-side clone hostile is equally concrete. If A is the ordinary
unlabelled tournament generating series, naive scalar substitution
A(z/(1+z)) predicts0 pair-free classes at order4, whereas the complete class
census gives2; at order6 it predicts20 rather than26. Exact unlabelled clone
substitution needs its automorphism cycle inventory. A module-tree shape or
scalar class series is not enough. The independent Newton audit supplies
the literal orbit census, keeping ordinary and pair-free counts distinct.

## 6. A new direct LCA connection: precision loss in the Newton bank

At the evaluation nodes 2,4,8,... the low-bit trie gives
`v_2(2^r-2^i)=min(r,i)`. The j-th Newton product
`q_j(X)=product_(i=1)^j(X-2^i)` therefore loses exactly
`tau_j=1+...+j` bits at every later evaluation node. This loss is the sum
of pairwise LCA depths with preceding nodes. Earlier nodes are exact zeros.

The normalized triangular evaluation matrix has units on its diagonal. Its exact
image modulo2^K has coordinates with precisions K,K-1,K-3,K-6,... until
nonpositive. At K=64, a sufficiently long full polynomial bank thus has
eleven mixed-precision coordinates, totalling484 residue bits, rather than
the previous64 full-precision coordinates totalling4096 bits. This is an
observation-image statement, not a Python-memory or runtime measurement.

The signed clone-excess transfer descends because each destination query
reads a source query in its retained positive-node horizon and multiplies
by integers and nonnegative powers of two. Division by two or by X would
break the claim and is explicitly excluded. Sparse exponent lattices can
have smaller images; the full-lattice invariant factors are not a claim of
minimality for every reachable bank. Ordinary A000568 was already scalar
at this point, so this bank reduction belongs to the signed core lane.

The separate `newton-lca/REPORT.md` supplies the complete prime-power-node
proof and independent finite controls, with explicit loss of unit factors
and no unsupported arbitrary-tree generalization.

## 7. An actual safe structural merger: forced unit skip contraction

When a counting stage ell exceeds remaining degree S, only the zero-cycle
skip is possible. Its factor is1 and exponent0, and the affordable defect
coordinates stay unchanged. Such a chain can be contracted to the next
affordable stage, with early S=0 sent to the common terminal. Deleting and
uniquely reinserting those forced edges gives a weight-preserving bijection
of complete paths. This supplies a valid continuation-exact structural
reduction, unlike equating arbitrary equal-potential states.

Stage still belongs to every surviving branching vertex. A k=0 branch at
an affordable stage is not necessarily forced or zero-cost: at N=8 the
actual state (3,3,(0,2)) has two branches and its zero-cycle edge costs3.
It cannot be contracted by this rule. Suffix potential is constant along
forced chains, but prefix minima and accumulated weights must be compared
at surviving aggregates, not separately at each removed incoming history.

This is elementary administrative compression implicit in the existing
recurrence, not a newly discovered equivalence of active defect states.

## 8. Proposed organizing principle and falsification rule

A proposed summary must specify its queries AND its allowed continuations.
To promote a quotient, exhibit a commuting operation rule and a kernel that
is stable under those continuations; preserve required weights and guards.
An exact query observer can still be invalid for new operations.

The cheap tests in this session are different by design:

1. Equal valuation code, different integer locus: 1 versus3/5.
2. Equal net resource, different prefix feasibility: (1,2) versus(2,1).
3. Equal suffix minimum, different odd weights:1 versus3.
4. Equal intended K target, failed intermediate blackness: the old six
   ledgers, now replaced on their complete offending-blue loci.

These support a usable research discipline, not an assertion that all repo
results are instances of one proved universal theorem. Whole black-core
connectivity, universal Collatz convergence, Antihydra's infinite survival,
and further compression of active ordinary defect states remain separate
obligations unless their native proofs are supplied. The proved unit-skip
contraction is a smaller implementation graph, not their resolution.
