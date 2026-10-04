# Keep unfinished Collatz observations: monotone closure and an exact edge frontier

2026-10-04. **PROVED** finite-graph soundness, monotonicity, and the scoped
literal-edge minimum. **FINITE-EXACT** for the stated512-input experiment.
Universal Collatz coverage and any runtime or symbolic-proof optimality claim
remain **OPEN / NOT CLAIMED**. No priority claim is made for memoization,
rooted reachability, or least-fixed-point evaluation.

Keeping every unfinished checked prefix changes the earlier finite experiment:
the union already certifies509 of the512 odd inputs1..1023, starting from1
alone. Exactly32 additional checked odd edges then complete all512. The new
object is a shared graph of verified observations, including unfinished work;
completion does not depend on which source supplied an edge first.

## 1. Inheritance and the changed operation

The [adaptive selector](adaptive_boundary_selector_20261004.md) retains exact
macro prefixes, smaller common-future obligations, and completed suffixes. Its
increasing-source pass at budget8 certified350 inputs with the basic adaptive
policy,354 with completed-suffix memory alone, and508 after adding the learned
27-cylinder plus completed-suffix memory. Four learned cylinders eventually
gave512 certificates and a shared835-node proof graph.

The incoming [checked-switch note](checked_switch_phase19_20261004.md) adds
the guarded reset-word witness

    (1^r,a) at n  <->  (1^(r-1),2,a-2) at (n-1)/2,
    r>=1, a>=3,

with equal actual endpoints. Its reset2 counterexample is the immediate
hostile: exponent0 is not an odd-map edge. The earlier source-debt hostile
27->41->31 also remains relevant; neither a smaller moving checkpoint nor an
unfinished prefix certifies its source.

The present operation keeps incomplete checked prefixes even when their
owner receives `PENDING`, unions all exact primitive edges, and propagates
certification backward from1. It also retains every advertised one-step
sibling edge. This is different from storing only already completed routes
or choosing one available child in an increasing pass.

The concept board is: exact source; unfinished prefix; checked edge; actual
join; grounded rank; missing-edge demand. The anchor is certificate reuse,
the niche is order-independent closure, and the wildcard is scheduling work
at a shared missing edge rather than restarting each unresolved input.

Source object: a finite collection of guarded selected prefixes and joins.
Target: a finite partial functional graph `n -> (U(n),v2(3n+1))`.
The map preserves every realized transition, actual labels, ordinary/odd
clocks and finite first-hit routes. Flattening loses symbolic family coverage
and the compact macro representation; original source/word-position provenance
is retained, while the family guard remains in its defining compiler. This
literal graph is not proposed as a replacement for all-height symbolic rules.

## 2. The exact graph and least grounded closure

An observation is a triple `(n,m,a)` of positive odd integersn,m and positive
integera, withn>1 and `3n+1=2^a m`. Every inserted triple is checked by an
independent while-even calculation. Conflicting successors are rejected.
The root self-edge1->1 is excluded.

For a finite graphE, define

    C_0={1},
    C_(k+1)=C_k union {n : an edge n->m is stored and m is in C_k},
    C(E)=union_k C_k.

**PROVED: soundness and monotonicity.** Every member newly added at stagek+1
has one exact odd edge into an already certified member, hence a finite route
to1. Conversely, a stored finite path of lengthk to1 puts its source inC_k.
ThusC(E) is exactly rooted reachability in the stored graph. It is its least
grounded fixed point. IfE is a subset ofE', thenC(E) is a subset ofC(E').
No ungrounded cycle can prove itself. Because each nonroot source has at most
one successor, the first-hit rank obeys `rank(n)=1+rank(U(n))` whenever grounded.

The implementation computes this closure by reverse breadth-first propagation.
An independent forward walk from every observed vertex gives the same initial
closed set and ranks. After closure, rank-ordered inverse-ray construction
turns stored edges into canonical home ASTs. Independent literal encoding is
used only afterward for audit; it supplies no seeds or observations.

For reset witnesses that reach1 early, the remaining formal word is checked
usingU(1)=1 but no root self-edge is inserted. The resulting stored path is
truncated at the first1. This reconciles the incoming switch formula with the
first-hit codec rather than weakening its root convention.

## 3. Initial universe and redundant alternative observations

The requested universe is the512 positive odd integers1..1023. For each input,
call `select(n,8,'adaptive')` with macro cap128 and retain its entire checked
actual word plus each advertised sibling-to-join edge. All resulting aggregate
words in this experiment also meet the128-letter observation cap. No supplied
home certificates or learned orbit seeds enter this graph.

There are1920 submitted primitive edge claims,802 distinct nonroot edges and
805 vertices. The root closure contains761 vertices, of which509 are requested.
The three missing requested inputs have the following exact frontier:

| Unobserved outgoing source | Requested inputs whose stored path ends there |
|---|---|
|2125|871|
|2287|703,937|

The saved output contains a complete sorted `INITIAL_GRAPH_JSON` record:
requested universe, budgets, all802 edge triples and the frontier. Its hash
allows the initial checked graph to be reused without rerunning selection.

**Same-observation control.** The previous suffix-memory algorithm skips
inputs already certified, whereas the census above calls the selector on
every requested input. Merely giving both procedures the same universe and
macro budget would not control their observations or cost. To resolve this,
the script wraps the frozen original `seed_closure` without altering any
choice and records exactly its returned words and advertised sibling edges.
It makes353 selector calls (352 nonroot calls plus root validation), submitting
1388 edge claims, and certifies354 inputs. Unioning precisely that trace gives
the identical802-edge graph and certifies509 inputs. Thus this finite gain
uses the same checked transition information, not additional distinct edges.
Raw selector-call counts and duplicate replay work still differ from the
full-census construction; no runtime comparison follows.

The following additive controls contribute **zero new distinct edges**:

| Added observations on the same requested universe | Submitted edge claims after union | New distinct edges |
|---|---:|---:|
| Eight-macro inherited-policy words |3616|0|
| Eight literal odd steps, stopping at1 |5004|0|
| Both words of every eligible reset switch |5770|0|
| The learned27-cylinder word at its matching requested source |5807|0|

There are128 eligible reset witnesses; both their actual and smaller-source
paths are checked. Source7's reset exponent2 is rejected. The learned27 word
has37 odd steps but every one was already observed somewhere in the union.
This is a precise stopping reason for adding these particular controls in this
finite universe, not a general redundancy theorem for reset switches or
learned cylinders. The final expansion below starts from the original802-edge
graph, not from the larger duplicate-submission count.

## 4. Query the shared missing edge

For each unresolved requested input, follow its stored path to either a
missing outgoing edge or an ungrounded cycle. Count how many requested inputs
depend on each missing source. Choose the source with largest demand, breaking
ties by smallest integer label. Query exactly one new odd transition, verify
its valuation independently, insert it, and recompute grounded closure.
Stop when all requested inputs are certified, there is no missing-edge
frontier left, or128 new queries have been made. A cycle remains unresolved.

**FINITE-EXACT:** the demand rule makes32 new queries and closes all512 inputs.
The first31 follow the common2287 frontier for703 and937. The last edge of
that chain is `3349 --a=6-->157`, which joins the existing root closure and
certifies both requested inputs. The remaining query is
`2125 --a=3-->797`, which certifies871. The full per-query list, demand and
newly certified inputs are saved in the output.

The final graph has834 nonroot edges and835 vertices; every vertex is now
grounded. It matches the prior fully completed proof graph's size, but required
no learned cylinder to expose its missing edges. All512 requested ASTs equal
the canonical first-hit encoder's output. Their separate odd-route lengths
sum to11373, while the shared literal graph stores834 distinct edges.

For a clean sharing comparison, give each of the three missing inputs a fresh
copy of the same initial802-edge cache, and do not share its new queries with
the other two. The required query counts are31 for703,1 for871, and31 for937:
63 queries for32 distinct edges. The common frontier avoids repeating the31
queries shared by703 and937. This comparison concerns exact edge observations,
not total CPU operations, macro verification, or the cost of initially building
the802-edge cache.

## 5. A representation-specific minimum

**PROVED, with a FINITE-EXACT cardinality certificate:** for this supplied
initial graph and these512 requested inputs,32 is the minimum number of new
literal odd edges required for completion.

Each positive source has a unique actual successor. Thus any complete
literal first-hit route certificate must contain every edge on its canonical
route to1. After completion the script takes the union of all512 such paths
and checks that it is exactly the834-edge final graph. The initial802 edges
are a subset of this required union: they came from actual prefixes of
requested inputs or from advertised smaller odd children, which also lie in
the requested universe. Every added query lies on a still-unresolved requested
path. Exactly32 required edges were missing initially. Any literal-edge
completion must supply them; the scheduler supplies each once.

The minimum is independent of this scheduler's priority rule. It is not a
claim that largest demand minimizes waiting time, arithmetic work, memory,
or total proof length. A symbolic macro can certify many edges at once or
give a different proof without listing them. No lower bound against such
representations, or against all Collatz methods, follows.

## 6. Reproduction, hostiles and next boundary

    python 04-computation/experiments/adaptive_observation_union_20261004.py
    python -O 04-computation/experiments/adaptive_observation_union_20261004.py

The [script](../../04-computation/experiments/adaptive_observation_union_20261004.py)
and [saved output](adaptive_observation_union_20261004.out) retain the initial
graph and each new query. All checks remain active under-O. Hostiles cover an
ungrounded formal cycle (explicitly not claimed to be Collatz data), a root
self-edge, a wrong successor, a boolean source, zero query budget, reset2,
and reverse insertion order. Insertion order changes neither the closed set
nor its exact first-hit ranks.

The successful move is to store checked unfinished work and schedule its
shared missing edge. It is finite and auditable. A new input can still produce
an indefinitely growing missing frontier; a hypothetical nonroot cycle would
remain ungrounded. A global termination or well-founded selection theorem is
still needed for universal convergence. The802 initial edges and32 added edges
are observations with explicit provenance, not free consequences of a phase
label or a finite experiment.
