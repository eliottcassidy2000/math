# A proposed framework: recognizable frames and guarded composition

**PROVED elementary criteria; independently reviewed synthesis.**
[Independent review](audit_integration/COMPOSITION_AUDIT.md). The elementary
criteria are proved here for their stated abstract models; no priority is claimed.
Their application is limited to the audited lane theorems linked by the
final packet. This does not assert that all repository results, or famous
open problems, are equivalent.

## 1. Three obligations, not three interchangeable invariants

Specify a source class, its boundary, the permitted operations and the
target observation. A compositional representation must answer three
different questions:

| Obligation | Exact question | A failure mechanism |
|---|---|---|
|Joint incidence|Which concrete boundary occurrences belong to the same object or witness?|Separate port or endpoint totals can mix incompatible realizations|
|Transport|Which identified boundary occurrences correspond under the next operation?|A permutation, character, phase or root-moving sign can be lost|
|Realization guard|Does the proposed glued state belong to the actual source image and satisfy acceptance?|A formal diagram can decode locally yet create forbidden new paths|

These are data and predicates of a declared experiment. Their names alone
do not supply maps between tournaments, integers, machines or dissections.

## 2. Elementary descent criterion for partial composition

Let X be a set of fully specified states, q:X→Y an observation, and let
each named operation have a specified arity. Include all boundary transport
choices as arguments or explicit operation labels. Write G_f for the
operation's legality predicate, and f for its result when legal.

**C1.** The operations descend to well-defined partial operations on q(X)
if and only if, for every operation and any two argument tuples with the
same coordinatewise q-values:

1. their legality predicates agree;
2. when legal, their results have the same q-value.

Necessity follows by applying the descended operation. Conversely define
its guard and value using any representatives; conditions1–2 make both
independent of that choice. This proves the criterion for all finite
arities, including nullary constants. Equality of a scalar statistic on
the source objects checks neither condition by itself.

For a weighted transition protocol, equality of the pushed-forward rows
is the corresponding condition. It is weaker than retaining each named
transition. The fresh-reset theorem in the incoming unnamed-carrier packet
is an example: its pair of observed source/target types closes that reset
alphabet but loses the alignment tested by coloured assembly. Its exact
success and exact failure are compatible applications of C1.

## 3. Recognizable encodings transport whole operation theories

Suppose E:X/~→Z is injective on the intended isomorphism classes, has a
decoder d on its image, and an exact image predicate R. Assume each operation
and its legality predicate already descend to X/~ by C1, and that d composed
with E is the identity there. Define the transported operation by

    f_E(z1,...,zr)=E(f(d(z1),...,d(zr)))

when all R(zi) and the original operation guard hold. Then E intertwines
every finite composition of the named operations, by induction on its
expression tree. Every identity between the original partial operations,
including its domain, transfers to this image. This is an elementary
conjugation theorem, not a claim that a native operation on Z has changed.

The four-marker encoding realizes this scheme for all distinct two-port
tournaments, with an explicit quadratic image guard. Its product on
unrooted encodings is decode–hinge–encode. Native Hamilton dynamics on the
larger encoded tournament need not agree with that transported product.
Hamilton support realizes the same scheme with a more expensive image
guard. Its deletion is decode–delete–resupport; induced deletion fails.

Isomorphism-class injection does not automatically preserve automorphisms.
The support encoding preserves them exactly. Four-marker encoding adds the
explicit kernel C3^epsilon at sink/source boundaries. Thus an automorphism-
weighted enumeration must remove or account for that kernel. Equal
unlabelled counts alone do not specify an attachment-compatible encoding;
in particular they do not supply a labelled bijection, as labelled counts
already differ at two vertices.

## 4. Compose compatible witnesses, then take counts or characters

For a finite witness system, form the set of tuples consisting of the
input witnesses and their actual boundary transports. Restrict to tuples
whose joint incidence and realization guards hold, then map each tuple to
its output. Pushforward of weights is the sum over that exact output
fibre. This is the precise meaning of a compositional response coefficient.

Associativity requires a hypothesis: the two parenthesizations must give
weight-preserving bijections with the same set of admissible multi-input
witnesses and the same output. Under that hypothesis, finite distributivity
proves associativity of the weighted convolution over any commutative
semiring. An arbitrary pairwise guard need not satisfy it. Counting recipes
also need not count distinct output objects; multiplicities or symmetry
normalizations must match the chosen target.

The path-cover expansion gives an exact witness-to-output bijection: split a Hamilton word at
the missing arcs. Vertex disjointness, exhaustive coverage and actual
connector incidences are retained. It yields a polynomial with nonnegative
coefficients, so evaluation at1 equals evaluation at0 exactly when no new
path is realized. The four-cycle control has defect z_02*z_13: neither
connector alone produces a new path, but the pair does.

Signed root-card transport has another evaluation. Marked deletion-card
isomorphisms extend to genuine automorphisms; their signs multiply under
composition. The character sum then cancels odd children and retains the
correct even-child orbit mass. Zero is intentional cancellation here,
not absence of witnesses. Keeping the realization guard separate from
the signed weight prevents precisely that confusion. The Legendre determinant
character on PGL2(F7) in the eight-point observer has the same multiplicative transport role, with a different
source, target and acceptance predicate.

## 5. Finite lifting and one common realization

A compatible residue tower or an arbitrarily long legal clock path may
have a realization at each finite stage without one realization in the
required positive-integer domain. The two quantifier patterns are

    for every L, some x_L realizes the first L requirements;
    some x realizes every requirement.

They are not equivalent. The growing-wheel class−1 and the all-odd clock
branch are exact boundary examples in the arithmetic lane. A finite clock
cycle is therefore not automatically a periodic positive integer orbit,
and a finite-state machine observation is not automatically a recurring
full machine configuration. The source and its guards must lift together.

There is a precise bounded repair. If all realization sets A_L are nested
subsets of one fixed finite source set X_H and every A_L is nonempty, then
their intersection is nonempty. Otherwise every x∈X_H has a stage at which
it disappears; the maximum of those finitely many stages empties a later
set, contradiction. There are at most |X_H| strict deletions, but this proof
gives no bound on the waiting times between them. A growing source-height
bound is not the same hypothesis.

The sieve supplies an effective cutoff on each fixed height interval.
Clock lifting supplies exact congruence cylinders and their positive
integer counts in such an interval. These are constructive finite
realization tests. Neither alone establishes unbounded twin primes,
Collatz termination, or a Busy Beaver classification.

## 6. What a larger theorem would still need

The supported proposal is a **theory of operation-relative recoverability**:
first seek an intrinsic decoder; otherwise retain a composable frame,
transport coordinate or increasing certified horizon. Then audit the image
and the actual next operation. This extends the earlier factorization
method to partial multi-input composition without declaring every boundary
representation minimal.

Three remaining obligations are concrete:

- Compress the joint path-cover/image guard while preserving its operations,
  rather than merely giving another exact exhaustive encoding.
- Find an attachment-compatible tournament/even-graph construction beyond
  the proven signed counting operator; scalar equinumerosity does not do it.
- Show that the realized black-line transports connect all terminal basins.
  Faithful representations and local common hubs do not establish that
  global reachability statement.

These are candidate directions for stronger unifying theorems. The proved
criteria explain why the present constructions succeed or fail, while
leaving their existence, efficiency and connectivity obligations visible.
