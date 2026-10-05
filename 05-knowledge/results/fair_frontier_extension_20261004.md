# Fair extension of supplied Collatz states with retained proof work

2026-10-04. **PROVED:** guarded-word home equivalence, rooted least closure,
and conditional completeness of the fair scheduler. **FINITE-EXACT:** the
frozen 239-source and separate 16-source controls below. **OPEN:** termination
for every positive odd input, and any generally superior symbolic policy.
No historical-priority claim.

Artifacts: [program](../../04-computation/experiments/fair_frontier_extension_20261004.py)
and [saved output](fair_frontier_extension_20261004.out).

## 1. Inheritance and the change of obligation

The closest mechanism is the [grounded observation union](adaptive_observation_union_20261004.md):
unfinished checked edges survive, and certificates are propagated from ROOT1.
The [frontier family compiler](frontier_family_compiler_20261004.md) turns a
checked word into an all-height source guard, without making every member
an unproved root seed. The incoming
[join shields](collatz_join_shields_and_lifts_20261004.md) explain why merely
deepening an inverse branch at a fixed forward window can waste work: 222
of its 239 frozen obligations have no smaller ancestor at any inverse depth
through the first six forward positions. Its mixed binary/ternary family
lifts and the [sixteen-row guard bank](collatz_binary_ternary_guard_fusion_20261004.md)
remain usable as candidate rules.

Canonical hostile: 27's first smaller-source coalescence is at odd step37,
so repeatedly increasing inverse depth at its earlier joins cannot replace
forward progress. The corrected near miss is treating completion of one
decoder as necessary for extension of the supplied integer. The concurrent
[half-child obstruction](half_child_extension_obstruction_20261004.md) makes
that failure exact at7: the fixed affine half-child representation can fail
at every length even though7 has a finite route. The least-used sidecar here
is a permanently retained cursor on each *original* request.

The concept board is **original request / exact word / rooted component /
resumable job / cost counters / first-hit certificate**. Anchor: fair supplied
state extension. Niche: bidirectional certificate transport through compressed
actual words. Wildcard: interrupt and resume with exactly the same proof state.
No tournament or phase label replaces a source identity.

## 2. Guarded actual words give bidirectional home equivalences

Write $U(n)=\operatorname{oddpart}(3n+1)$ on positive odd integers.
A nonempty positive valuation word $w$, of length $j$, cost $A$, and
carry $B$, has $P=3^j,Q=2^A$ and formal endpoint

\[
 z=(Pn+B)/Q.
\]

**PROVED exact guard.** For positive odd $n,z$, equality $Pn+B=Qz$
is sufficient as well as necessary for $w$ to be the exact actual word.
Equivalently the exact source cylinder is $Pn+B=Q\pmod{2Q}$.
This uses the endpoint's oddness, not merely divisibility by $Q$.
To see sufficiency, reverse the divisions from $z$. Integrality of the
complete inverse and reduction modulo3 force the last inverse to be an
integer; divide out that factor and repeat. Each integer inverse is odd
and positive. Thus every prescribed two-adic valuation is exact.

Store an arc $n\xrightarrow{w}z$, its source, target, word, and origin.
Its insertion checks the complete affine identity and exact odd endpoint.
A literal observation also receives an independent while-even division check.
Guarded macro insertion uses the proved word guard; its expanded literal
validation is separately counted after search.

For every actual finite word,

\[
 n\text{ reaches }1\quad\Longleftrightarrow\quad z\text{ reaches }1.
\]

The reverse implication prepends the word. For the forward implication,
cut the supplied route after $j$ steps. If its first root occurs earlier,
all later actual steps stay at1, so $z=1$. Consequently the least home
closure of a finite guarded-word graph is the **undirected component of
ROOT1**. This is still a least fixed point, using the two valid implications
of each checked word. An unrooted cycle cannot certify itself.

This differs from merely labelling an arbitrary symmetric relation a proof.
Every edge here retains a directed, exactly guarded actual Collatz word.
The common-future diagram $n\to z\leftarrow m$ joins proof components,
and either supplied home certificate transports across it. No smaller-child
premise is silently assumed. The graph initially contains no certified
nonroot vertex.

**First-hit export.** Starting at ROOT, traverse a rooted spanning tree.
Backward traversal prepends the checked inverse word. Forward traversal cuts
the existing certificate, checking each actual exponent. Formal trailing
root steps have exponent2 and are discarded. This constructs canonical
first-hit inverse ASTs, not merely a Boolean home label. Source1 self-arcs
are omitted. The implementation supplies the codec's conservative full-suffix
bit bound for each expansion; its post-search comparison with an independent literal encoder
is validation, never an input seed or a search oracle.

## 3. A fair scheduler with an explicit budget boundary

The API `FairSearch(requested, rules=True)` takes a finite, nonempty list of
distinct positive odd integers. Every request has an immutable original
label and a current cursor connected to it by retained actual arcs. A
request's cursor is never replaced by an unproved smaller child.

One scheduler round performs:

1. One unresolved original request, chosen round-robin, advances by a known
   actual word if available, or queries and validates exactly one new odd
   edge. At least one actual odd step is represented by an unresolved visit.
2. At most one resumable auxiliary action performs one arithmetic family
   guard, advances one candidate-child cursor, or advances one unfinished
   decoder job. Alternating auxiliary turns give priority to already-guarded
   child walks; the other turns retain FIFO progress for the rest.

The production family actions use the incoming sixteen debt rows, the
reset-at-least-three switch, four stored coarse words (including the learned
2287 tail and the two 223-derived first-descent words), and the incoming
ternary sibling guards. A source-specific sibling index is tested in one
action; the proved positivity bound supplies its finite index cap. Every
returned actual arc is checked before it enters closure. Applicable children
are additional work, never replacement root seeds.

**PROVED conditional completeness.** Let the fixed request list have size
$N$. If one supplied source has actual first-hit odd rank $T<\infty$, it
is grounded after at most $NT$ rounds, unless already grounded sooner.
This remains true if other requests never terminate or an auxiliary job
never finishes.

Indeed an unresolved request receives a cursor turn at least once every
$N$ rounds. Each turn advances along its actual trajectory by at least one
odd step. A guarded jump that crosses the first root ends at1 and grounds
the request. Otherwise at most $T$ such visits reach ROOT. Extra correct
arcs can only enlarge the grounded component. Thus every finite batch whose
members have finite routes eventually completes; a mixed batch eventually
certifies every member that has a finite route while retaining the rest.

The bound is in logical rounds, not machine instructions or elapsed time.
Each individual auxiliary action must terminate; the scheduler does not
preempt an arbitrary nonreturning Python function. Inverse or decoder
searches must expose finite resumable actions. Arbitrary-precision arithmetic,
graph scans, guard costs, and word lengths are not constant-time operations.

`run(rounds)` returns `COMPLETE` only when all original requests belong to
the rooted component. Otherwise it returns `PENDING`, preserving exact arcs,
cursors, jobs, and counters. JSON `snapshot`/`restore` retains the schedule;
restoration rechecks arcs and each cursor's retained directed path from its
original request. A zero budget adds no evidence. `PENDING` is not evidence
of divergence and is not a failed certificate supplied under another label.

There is no claim that this theorem makes every $T$ finite. Universal
termination of this sound, conditionally complete procedure is equivalent
to the positive odd Collatz convergence assertion. It is not a new proof of
that assertion. The useful new implementation guarantee is that an optional
symbolic branch cannot suppress the complete literal route search.

## 4. Frozen 239-source benchmark, with honest work accounting

The exact ordered requests are the 239 inherited seeds from
`checked_switch_phase19_20261004.json`, already frozen by the incoming
shield experiment. Their compact JSON SHA256 is

    c1a01b7d890f9bcf6205572a65f112b4038c290b67c10ef86e5d1eba40df31b7

This experiment starts again from ROOT1 and **no supplied child routes**.
It does not rerun the old online learner or reinterpret smaller-child joins
as completed sources. The declared cap is10000 rounds.

| Search | All239 complete at round | Original literal queries | Auxiliary literal queries | Guard tests | Guarded arcs / stored letters |
|---|---:|---:|---:|---:|---:|
| Literal only |4101|3482|0|0|0 / 0|
| Guarded extension |3708|3332|135|2588|28 / 170|

At round2000 the literal policy has grounded0 requests and the guarded policy
22; the latter has used2071 literal queries against1980. At round4000 the
literal policy has grounded223, while the guarded run has already completed.
These round comparisons do not equate unequal amounts of work.

The completed guarded run uses3467 fresh literal queries, fifteen fewer
than3482, plus2588 exact arithmetic guard tests. It stores3495 proof arcs
and3478 vertex labels; the literal run stores3482 arcs and3483 vertices.
After search, independent expanded checks replay3637 arc letters for the
guarded graph, against3482 for the literal graph. Certificate construction
uses3619 versus3482 splice letters. Both export all239 canonical ASTs and
agree with independently replayed first-hit routes totalling9300 odd edges.

This is a small reduction in fresh literal observations traded for additional
symbolic and validation work. It is not a runtime speedup, a minimum-query
claim, or a larger convergence theorem. The old834-edge minimum concerned
a different source universe and a literal-edge-only representation; it
does not constrain these guarded word arcs.

Here a literal query specifically means an uncached one-edge action by an
original or child cursor. Guard tests themselves also calculate exact
valuations and endpoints and can certify an entire word. The difference of
fifteen is therefore not a count of fifteen fewer arithmetic evaluations of
$3n+1$.

## 5. Separate positive lifts, interruption, and adversarial controls

The second prespecified universe has16 sources. For each frozen debt row,
take its recorded first tested source $n_0$ and actual parent word $w$,
then request $n_0+2^{\sum w+1}$, the next exact-cylinder lift. These sources
are disjoint from the239 benchmark and their list is printed in stdout.
They were selected by source guards, not by checking for fast home routes.
Their prior sample's certificate is not supplied. The cap is again10000.

| Search | Completion round | Original / auxiliary literal queries | Guard tests | Stored macro letters |
|---|---:|---:|---:|---:|
| Literal only |571|551 / 0|0|0|
| Guarded extension |565|449 / 102|208|298|

Both use551 fresh literal queries. Expanded post-search validation uses551
and849 letters respectively; independent requested routes total729 odd edges.
This control prevents presenting every successful family match as a saving.

Additional exact controls are:

- The full guarded239 search is identical after a zero budget,137 rounds,
  JSON round-trip,863 rounds, a second round-trip, and9000 further rounds,
  compared with one10000-round call. Equality includes arcs, cursors, queue,
  counters, and query trace, not just the final certificate count.
- An intentionally unfinished auxiliary decoder yields one empty finite
  action forever. Requests27 and703 still finish in102 rounds and102 fresh
  literal queries, while that decoder advances102 times and remains pending.
  This is a scheduling hostile, not a Collatz oracle or an actual nonterminating
  arithmetic decoder claim.
- Wrong source/endpoint data, a root self-arc, and Boolean, floating, zero,
  negative, or even supplied sources are rejected. A checked disconnected
  edge is not grounded. A separate generic reachability control leaves an
  abstract unrooted cycle unsupported; it does not assert a positive Collatz
  cycle exists. Independent graph traversal agrees with component closure.
- A singleton27 export regression retains the codec's conservative suffix
  bound. An initial local implementation used endpoint bit lengths instead,
  which understated that API bound and rejected a completed certificate.
  The repaired bound changes validation bookkeeping, not the mathematical
  word guard or any reported search count.

Reproduction:

    python 04-computation/experiments/fair_frontier_extension_20261004.py
    python -O 04-computation/experiments/fair_frontier_extension_20261004.py

The optional export writes exactly the239 guarded-run requested certificates,
as first-hit valuation words with source identities, ranks, provenance,
status and the frozen-universe hash:

    python 04-computation/experiments/fair_frontier_extension_20261004.py --export-certificates 05-knowledge/results/fair_frontier_extension_20261004.json

The serialized words are independently replayed with strict first-root
stopping, then the written JSON is read back and compared. Export happens
only after completed search; these words are never search inputs. Without
the option, the program creates no output file.

All mathematical checks use explicit exceptions. The current source domain
is positive odd integers; even-source normalization and signed basins are
separate typed interfaces, not implicit extensions of this ROOT1 search.

## 6. What to improve next

The guarantee is now independent of the quality of the symbolic policy.
The present fixed bank saves little literal work on these universes and adds
substantial guard work. The next controlled experiment should schedule a
guard according to the unresolved requests its endpoint component could
connect, while preserving the same reserved original-request turns. Compare
fresh queries, guarded word lengths, and validation costs separately.
The incoming shield and the finite supplied-child decoder bound can eliminate
specific unproductive jobs; neither may delete the original literal cursor
or rule out different child maps and unequal-length joins.
