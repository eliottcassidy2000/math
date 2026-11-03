# A common seventeen-coordinate two-port observer

**PROVED by the endpoint-cover bijection and exact finite matrix identities:
the old static two-feedback language and the new two-live-port adaptive
language sit inside one twelve-command operation with the same minimal
rational observation dimension17.** Seventeen static continuation queries
already determine every future Hamiltonian-count response, including
commands that replace a retained port. This identifies a common response
quotient, not the graph languages or their numerical images. No spectrum
coverage claim is made for the twelve-command language.

This result follows the earlier comparison in
[live_port_audit section4](live_port_audit.md). There the proper common
two-letter sublanguage had finite-field rank13, which did not explain the
two full ranks17. The coordinator's larger twelve-command operation supplies
the missing boundary map. The independent audit below verifies that map on
the entire25-coordinate carrier, rather than inferring it from equal ranks.

## 1. The operation and its universal endpoint carrier

Keep two named vertices in ordered slots0,1. A command `(m,d)` consists of
an incoming mask m in `{0,1,2,3}` and a discarded slot d in `{0,1,2}`. Add a
new vertex x that loses to retained vertex i exactly when bit i of m is set,
and dominates every other old vertex. Then discard slot d, where slot2 means
the newborn. If d<2, x replaces that old slot. Discarded vertices remain in
the graph but can never be named again. There are twelve commands, encoded
`3m+d`.

As before, retain the unordered multiset of paired `(start,end)` labels in
a spanning path cover, replacing forgotten endpoints by `*`. At most one
forgotten end can survive: no forgotten vertex has an edge to a future
newborn. A component entering x may end at any retained port in m; one
leaving x may start at a retained port outside m or at `*`. At most one of
each is selected, and selecting the same component twice is forbidden.
Forgetting is then applied. The path-cover trajectory bijection proves the
exact nonnegative integer transfer for all command words.

There are exactly25 possible nondegenerate signatures, not merely25 found
in one experimental source. The two-port abstract count is28, including
the empty signature. An active port occurs at most once among all endpoints,
except that `(i,i)` denotes its singleton component. The empty cover is
impossible on a nonempty source; a whole-cover signature `((0,0))` or
`((1,1))` is impossible when there are at least two vertices. Removing these
three cases leaves25. The independent audit enumerates this entire abstract
set and checks that every member is already reachable in the static
four-command sublanguage from TT2.

Consequently the25-coordinate carrier and all its matrix identities apply
to **any actual tournament of order at least two with two named vertices**,
using its actual initial path-cover vector. No assertion that every such
tournament was generated from TT2 is needed. The finite carrier describes
the future operation, not a classification of its sources.

## 2. The seventeen static queries close under all twelve commands

Let `M_(m,d)` be the25-by-25 transfer matrices, f the one-component observer,
and u the actual TT2 initial vector with its arc from slot0 to slot1. A
static command keeps both old ports, so d=2. Write `S_m=M_(m,2)`.

The exact certificate uses only the following seventeen static query words,
written in chronological mask order:

```
empty, 1, 2, 3,
01, 11, 21, 31, 02, 12, 22, 32,
201, 301, 021, 221, 321.
```

Let C be the25-by-17 matrix with columns `S_q f` for these words q. Its
columns are independent, its first column is f, and the frozen certificate
contains explicit17-by-17 rational matrices A_(m,d) such that

`M_(m,d) C = C A_(m,d)` for every m,d.                 (C1)

These are coefficient identities on **all25 coordinates**. They are not
fits to the TT2 reachable orbit, a finite set of Hamiltonian values, or a
test modulo one prime. The independent audit reconstructs every positive
matrix from component-chain joins, reconstructs C from the static words,
and verifies every entry of all twelve identities over Q.

For any actual initial cover row v, define

`X(v)=vC=(H after query q_1,...,H after query q_17)`.

Then (C1) gives the exact update `X(vM_(m,d))=X(v)A_(m,d)`; H is its first
coordinate. Thus the span of observations under the four static commands
equals the span under all twelve commands. In particular changing a port
requires no extra rational observation coordinate in this operation.

This is the coarsest deterministic future-H observer in the following
precise sense: two actual two-port sources have equal X if and only if all
twelve-command continuations have equal H. Equality of X implies equality
of all responses by (C1); the reverse implication holds because the
seventeen listed queries are themselves permitted continuations. This does
not identify the original graphs or claim that every positive integral
X is an actual source.

From TT2 the static four-command reachable row space has dimension25 and
its scalar Hankel rank is17. The twelve-command reachable space also has
dimension25 and its scalar rank is17. The upper bound is (C1); a17-rank
actual prefix/query matrix gives the lower bound. These statements concern
minimal **linear** dimension, not the number of numerical states or the
size of a nonlinear encoding.

## 3. The old static language is mapped literally, with reversal retained

The old language in
[October30, two feedback vertices and minimal memory](../square_sum_spectrum_overnight_20261030/README.md)
fixes vertices x->y and a chain `v_1->...->v_n`. A letter a has bit i set
exactly when marked vertex i points to that chain vertex. Its positive
gap-assignment model has36 coordinates and rational Hankel rank17.

Starting the new operation at TT2 and using only static commands produces
exactly that graph family after reversing the chronology of the unmarked
vertices. Every newborn dominates all older unmarked vertices, so newest
first is their transitive-chain order. Its incoming mask is **the same** old
incidence code, with no bit complement. Hence

`new chronological masks m_1...m_n`
` -> old chain word m_n...m_1`.                       (C2)

The map fixes the two marked vertices and reverses the other labels. It
preserves every arc and hence the full Hamiltonian-path multiplicity. The
independent audit checks all85 four-letter words of lengths0 through3 as
literal adjacency matrices.

Reversal is essential for interpreting the observer. Appending a new static
query in chronological order **prepends its reversal** to the old chain
word. The seventeen queries above therefore describe old left-continuation
responses. They are not asserted to be the identical coordinate list of
right-suffix probes printed in October30. Reversing a scalar series
transposes its Hankel matrix after exchanging prefix and suffix indices,
so minimal rank is preserved. The new25-coordinate model is a distinct
positive realization of the reversed old series, with an explicit
extension to eight additional port-replacement commands.

This supplies the structural explanation for the shared17: both responses
factor through the same two-port continuation quotient after the explicit
reversal, and that quotient is stable under port replacement. Equality of
two dimension numbers alone was not the proof.

## 4. The adaptive six commands embed, while graph languages stay distinct

Command8 is `(m=2,d=2)`. From TT2 it creates the literal C3 rows `(2,4,1)`
and keeps ports0,1. Thereafter the adaptive two-port operation is precisely
the six commands with m=1 or m=2. A root-slot command `(r,d)` becomes

`(m=2^r,d)`, encoded `3*2^r+d`.                       (C3)

Thus the earlier24-coordinate C3 carrier embeds in the universal25-coordinate
space, and (C1) supplies a17-coordinate observer for it. Its independently
proved Hankel rank17 shows that this restriction loses no observation
dimension. It does lose physical source states and command choices.

The static and adaptive graph languages still differ. Static words admit
a fixed pair whose deletion leaves a transitive tournament. The adaptive
newborn-only sublanguage has arbitrarily many vertex-disjoint directed
triangles and unbounded feedback vertex number, as proved in
[ACTIVE_PORTS section6](../adaptive_metagraph_reformulation_20261103/ACTIVE_PORTS.md).
The common observer does not turn that family into a fixed-feedback one.
The full twelve-command language also includes nonstrong sources and
zero-growth commands. No unrestricted graph-class equivalence is asserted.

## 5. Zero-growth commands destroy the bounded-H finiteness argument

Repeat `(m=0,d=2)` from TT2. Each new vertex is a source and is immediately
forgotten. The graph remains transitive and H stays1 for arbitrarily many
births. The positive cover vector is not bounded: after n>=2 such births,
the signature

`{(*,*),(*,0),(1,1)}`

has coefficient `2^n-2`. Partition the n forgotten chain vertices into two
nonempty subsets. One forms the forgotten-end component and the other
forms the path ending at port0; port1 is a singleton. Their orders are
forced by transitivity. Every nonempty proper subset gives exactly one
such cover, proving the coefficient formula.

Therefore even the set of reachable positive carrier vectors with H=1 is
infinite. The adaptive language's strict-growth termination and its finite
target-saturation receipts cannot be transferred to the twelve-command
envelope. A terminating inverse would need a new normalization or descent
that preserves all needed source/continuation data. No such inverse or
coverage theorem is supplied here.

## 6. The observer still needs an actual-source guard

The TT2 seed cover vector u is nonnegative integral and has H=1. The formal
vector7u, and its seventeen queries7uC, are also nonnegative integral and
satisfy every homogeneous transfer identity. Every formal continuation
count is odd, since it is seven times the odd count of an actual TT2-word
tournament. Nevertheless7u cannot be an actual tournament's cover vector:
its H coordinate is7, excluded by
[THM-1370, h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).

This specializes the inherited October30 hostile `7X(empty)` to the common
observer. It is not a new exclusion proof. Exact linear closure, integrality,
nonnegative observation values and parity still do not characterize the
actual graph image. A target constructor must retain an actual word or
another proved realization certificate.

## 7. Verification and source boundary

The producer is `two_port_envelope.py/json`. The independent
`two_port_envelope_independent_audit.py/json` imports no repository module
and verifies:

* all25 abstract signatures and all twelve matrices;
* all1885 programs through three births by literal permutation H counts;
* all157 programs through two births by full literal path-cover vectors;
* all425 proposed full-prefix/static-query Hankel entries by an independent
  backward path recursion, on321 distinct graphs of order at most8;
* static reachable rank25, the17 independent static query columns, and
  all twelve exact rational coefficient identities on25 coordinates;
* the85 static graph maps, literal C3 entrance, and zero-growth hostile.

No target-support scan is run. Normal and optimized audit receipts agree:
SHA256 `54b7bb046da4286fc76d0f274ea05fa03e19a2fe255e3cdf9b386bd91c3d9094`.
The receipt pins the producer JSON hash and retains the actual static basis
words. Independent finite controls audit the implementation; the all-word
conclusions follow from the transfer bijection and exact invariant-subspace
identities.
