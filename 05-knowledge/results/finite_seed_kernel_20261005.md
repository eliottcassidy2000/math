# Finite root assumptions, recursive families, and an exact certificate kernel

2026-10-05. **PROVED:** the finite-assumption kernel and the scoped logical
assembly rules below. **PROVED in the linked packages:** the recursive rooted
tree, guarded symbolic family compiler, and bootstrap boundaries.
**FINITE-EXACT:** the declared point experiments. **OPEN:** covering every
positive integer by these constructions. No general termination oracle or
literature-priority claim.

The user's proposed move is sound: first prove

    [Root(s1) and ... and Root(sk)] implies [every n in D reaches1],

then supply and verify k actual finite routes. The work lies in constructing
the universal implication without leaving an unrecorded family of additional
assumptions. We now have an explicit recursive D with k=1, plus a compiler
that tracks and discharges finite assumption manifests.

## 1. Inheritance and the two kinds of leaf

Closest proved mechanism: the bidirectional guarded home equivalence and
rooted least closure in [fair frontier extension](fair_frontier_extension_20261004.md),
section2; and the [sibling common-future grammar](creative_sibling_20260925.md).
Canonical hostile: an ungrounded component cannot prove itself by cycling
through valid equivalences. Corrected near miss: a descending family may
leave the family, so its finite boundary certificates alone do not close it.
Least-used sidecar: a finite explicit assumption manifest, rather than a
Boolean label saying that the destination is presumed rooted.

Anchor: discharge finite root obligations. Niche: minimal kernels in a frozen
observation graph. Wildcard: recursive refuelling that repairs a depleted
ternary congruence. Board: **point assumption / universal family law / native
guard / common future / decreasing rank / retained prime precision**.

Two leaves must remain different:

* A point leaf Root(27) requests one actual finite certificate. Its verification
  terminates for any supplied finite word; finding that word remains a separate
  search, with an explicit cap in the experiment.
* A schema leaf, for example 'all odd children e+Pk reach1', has an unbounded
  parameter. Writing it on one line does not make it one point obligation.

The small hostile is D={positive odd n: n=1 mod4}. Every n>1 in D satisfies
U(n)<n, and the boundary1 is certified. But9->7 leaves D. Thus these facts
alone do not give a standalone proof of Root(D). Strong induction over all
positive integers would need rules for the outside children too.

## 2. Minimal finite assumption sets for a fixed recorded graph (PROVED)

Let U(n)=oddpart(3n+1), with Root(n) meaning a finite route to1. A checked
common-future record consists of two actual finite words

    u --a--> z <--b-- v.

It proves Root(u) iff Root(v). Each word retains its exact source, valuations
and endpoint; arbitrary symmetric pairs are not accepted. The implementation
requires first-hit-safe words: no padding after reaching1.

Let G be a finite graph of such records, X its finite requested vertex set,
and A a set of vertices carrying already verified root certificates; include1
in A whenever1 is a vertex. For the following fixed-graph count, a newly
supplied point certificate grounds its attached vertex only. Importing its
internal edges or suffix vertices changes the deduction graph and belongs to
the adaptive version below.

**Finite-kernel theorem.** One new point certificate in each component that
meets X but not A is necessary and sufficient to ground all of X. In
particular the smallest number of additional point assumptions, relative to
this fixed record set, is exactly the number of those components.

Sufficiency: choose a representative and transport its certificate along a
spanning tree. If the known route from v begins with b, remove that prefix
and prepend a to obtain the route from u. Determinism ensures prefix agreement;
literal validation enforces it. Necessity: none of the recorded implications
crosses a component, so without a supplied certificate an ungrounded component
has no derivation from this rule set. This is a proof-system lower bound, not
a claim that its integers really fail to reach1.

Using a finite graph never requires a cyclic proof to be accepted. The two
sound implications Root(3)<->Root(5), taken alone, still leave one pending
component. Supplying the checked word(4) for5 or(1,4) for3 grounds it.

If G contains one actual outgoing observation for every odd n<=N, every
least representative of an ungrounded component meeting that interval is
3 mod4. Its minimum is <=N. If it were1 mod4 and greater than1, its observed
successor would be smaller and in the same component, a contradiction.
This identifies a precise reduced growth-state frontier; it does not bound
how many such components persist as N increases.

## 3. A finite observation experiment and actual discharge

Freeze four odd steps from each requested odd source, stopping at1. No seed
other than1 is initially rooted. The declared universe is separate for each row.

| requested odd sources | distinct observed edges | unresolved components | one minimal manifest |
|---|---:|---:|---|
|1..31|21|1|27|
|1..127|86|5|27,63,111,123,127|
|1..255|171|7|27,111,127,159,223,231,255|

Independent literal discovery supplies certificates for every named seed;
the last row's odd-step lengths are41,24,15,18,24,46,15. Spanning-tree transport
then exports exact first-hit words for all128 requested odd integers in1..255.
Every exported word agrees with an independent forward replay.

Certificates can themselves add observations. In a separate adaptive run on
the same255 snapshot, querying seeds27,127,159,223,231,255 inserts26 previously
unrecorded edges, for197 total. The separate111 obligation disappears because
the new27 route connects its component to the root. Six queries here do not
contradict the fixed-graph minimum of seven: the record set changed. No minimum
search cost or speed improvement over ordinary finite verification is claimed.

A cap-one attempt at27 stays unresolved. An empty alleged certificate for3,
a false common endpoint, a Boolean source, and a padded root word are rejected.
Neither a pending assumption nor a search timeout is promoted to a theorem.

## 4. One checked seed supports an infinite recursively closed family

The [refuel-tree package](finite_seed_refuel_tree_20261005.md) supplies an
all-height example with the single seed3 and its checked route3->5->1.
Let

    S(u)=4u+1,
    H(n)=(729n+669)/1024,   native source n=155 mod2048.

For a parent u=3 mod8, choose the unique kappa(u) in{1,...,729} satisfying

    4^kappa(u) (3u+1)=334 mod2187.

It exists because4 has order729 in the subgroup1+3Z modulo2187. For any branch
t>=0, set k=kappa(u)+729t and

    child(u,t) = [1024*S^k(u)-669]/729
               = [1024*4^k*(3u+1)-3031]/2187.          (1)

This is a positive odd integer in the exact H guard, and
`u < S^k(u) < child(u,t)`. The sibling identity U(S^k(u))=U(u), together with
H's checked common-future rule, transports a supplied proof of Root(u) to
Root(child(u,t)). It is not necessary to assume every intermediate integer
rooted: the two actual paths to their common future provide the transport.

Starting from3, repeat (1) along any finite tuple of nonnegative branch indices.
This gives a countably branching rooted tree at every depth. The whole family
is proved from one finite base certificate and one universal constructor law.
No new root assumption is added at a generation. The odd first-hit rank is
exactly6d+2 at generation d, independently of the branch sizes.

The decoder recovers the parent from a supplied literal child: compute c=H(n),
then k=(v2(3c+1)-1)/2 and u=((3c+1)/4^k-1)/3. Require the native guard,
positive integral k, the parent phase, the kappa congruence and u<n; repeat
until reaching3 or a failed guard. Hence membership is a total decision,
with actual route extraction on acceptance. Native guard membership alone
is insufficient:155 has H(155)=111 and k=0, so it is rejected.

The family has an additional exact address property. For fixed u and distinct
branches t,s,

    v3(child(u,t)-child(u,s))=v3(t-s).

Thus t mod3^a bijects onto every ternary source residue at precision a.
This constructs rooted members in every ternary address; it does not identify
a prescribed integer with one of those members. Their binary source guard
remains155 mod2048. Six parent ternary digits per generation must be retained
when reading an unexpanded source: a one-residue observer is not a complete
recursive state.

This example establishes infinite root closure, not new generic first-descent
coverage: its nonseed members first descend at the seventh accelerated odd step.
Its new use here is finite-base proof reuse, reversible addresses, and a
precise account of the extra congruence information consumed by recursion.

## 5. Conditional symbolic programs and finite fuel are complementary

The [finite-seed receipt compiler](finite_seed_receipt_compiler_20261005.md)
retains explicit assumption leaves, verifies acyclic proof dependencies, and
substitutes checked seed certificates. Its symbolic family compiler closes
selected controller programs against an inverse ray to a supplied seed,
retaining the whole binary/ternary guard and a first-hit-safe route.

For example LG with seed7 gives the family

    a=46+162t,   n=(896*2^a-287)/243,   t>=0.

The actual word `(1,2,1,1,a+2)` reaches7; its checked suffix `(1,1,2,3,4)`
then reaches1. This is one finite root obligation reused over infinitely many
symbolically represented inputs. The program language and its selected ray
parameters are explicit; arbitrary-source coverage is not inferred.

The [bootstrap-boundary package](finite_seed_bootstrap_boundary_20261005.md)
explains why refuelling matters. For the paid L map, z=7n+3 transforms by
z->9z/16; its native guard is v2(z)>=8. Every inverse L step consumes two
factors of3 from7u+3. A finite endpoint-seed set therefore has only finitely
many L-only ancestors. An arbitrary L guard cannot be closed by finitely
many L-only seed chains. H plus sibling refuelling changes the operation,
rather than assuming that an exhausted inverse chain can continue.

There is also a strong coverage test: the boundary package constructively
proves that rooting every member of any complete odd dyadic progression tail
would already root every positive integer. It preserves a common future,
not a claim that each forward orbit visits that progression. This prevents
mistaking a seemingly narrow whole-cylinder home claim for an easy partial
result. The rooted tree and displayed LG family are selected subsets of
their guards, not whole rooted dyadic tails.

## 6. How the proof pieces can be assembled (PROVED conditional rule)

For a domain D, retain a finite certified base set S, a well-founded rank rho
on D, and a rule that for every n in D outside S provides a sound common-future
dependency to m in D with rho(m)<rho(n). Then all of D is rooted by
well-founded induction. The finite family of seed certificates discharges
the point leaves; the universal rule and closure in D discharge the schema
leaves. A parameter rank is allowed even if an intermediate orbit value rises.
Each implication must still retain its actual guarded receipt.

Several domains can share certified seeds and share joins. Apply the finite
kernel theorem to their finite point interfaces; merge components when a new
checked receipt connects them. Their union is proved rooted once each domain
has its own universal closure argument. To conclude Collatz one must also
prove that this union includes every positive odd integer. Finitely many
schema names or successful residue probes do not supply that final cover.

The next productive target is explicit entry from an uncovered source into
one of these finitely grounded recursive domains, or a second closed domain
with a different decoder. Each new rule should report its new point seeds,
its universally quantified leaves, and its exact change to coverage.

Incoming `2498aa9be` repairs the Fourier sampling interpretation and an index
in the published injectivity argument. This session uses no statistical drift
claim; the prior [constructive phase decoder](translation_phase_decoder_20261005.md)
and its independently proved Haar second moment remain compatible with the
repaired statement. The changed board emphasizes source-level closure over
phase abundance or a numerical tail statistic.

## Reproduction

Run `python -X utf8 -B 04-computation/experiments/finite_seed_kernel_20261005.py`;
add `--write` to refresh the [saved output](finite_seed_kernel_20261005.out).
The companion uses only the standard library. Normal and optimized Python
keep explicit checks. Finite universes, source manifests, caps and hostiles
are specified above; the linked packages independently test the infinite
construction laws through exact symbolic readers and literal controls.
