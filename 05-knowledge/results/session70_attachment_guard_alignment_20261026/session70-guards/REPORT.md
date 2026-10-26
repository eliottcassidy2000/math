# Legal re-marking guards reconstruct every finite tournament

**Status: PROVED relative to the explicitly CITED nontraceable-arc structure
 theorem of Busch–Jacobson–Reid.** The new deduction and decoder below are
 elementary. Finite tests are corroboration, not the all-order proof. No
 literature-priority or formal-assistant claim is made.

## Result first

For every finite labelled tournament `T`, its set `HP(T)` of directed
Hamilton vertex words determines **every labelled arc**. Consequently the
complete named-permutation legality guards separate all fixed-path states
`X_n`, at every size. The exact smallest deterministic continuation carrier
for this operation protocol therefore has

    |X_n| = 2^M states,    M = (n-1)(n-2)/2,    n >= 1.

The fixed-length binary state cost is exactly `M` (zero at orders one and
two). This is a lower bound for preserving the declared legality observations,
not merely for reconstructing raw data by stipulation.

There is a real distinction at the word “guards”:

* **All named permutation guards, observed once:** reconstruct at all sizes.
* **All elementary single-vertex relocation guards, observed only now:**
  do not reconstruct; the first failure is order four.
* **All accessible finite traces of those named elementary relocations:**
  reconstruct at all sizes, using the inherited relocation-connectivity
  theorem. Complement and colour observations are not needed for this lower
  bound.

This closes the all-order question left open in the pinned strong-atom
packet C5. It does **not** refute smaller aggregate, unlabelled, unnamed,
freshly resampled, or distribution-valued response interfaces.

## 1. Inheritance and fixed scope

Research pin: `fb872272125510c317b8114dd3f697b433c595f9`.
Exact repository source identities are in `SOURCE_PINS.json`.

The closest proved mechanism is strong-atom C1/C5 in
`05-knowledge/results/strong_atom_continuation_vitali_20261026/README.md`:
admissible source-target alignment reconstructs the agreement path, and
continuation equivalence is the coarsest operation-stable refinement of the
chosen observer. The canonical hostile is its C3 pair at order seven,
which agrees on both jointly rooted endpoint types but disagrees on a named
re-marking guard. C4 amplified this to an unbounded alignment obstruction.
The corrected near miss is “every arc appears in a spanning Hamilton path”:
it is false even in a strong tournament. The least-used useful sidecar is
**one-vertex-deleted reachability of the traceable-arc core**.

Session69 P3, in
`05-knowledge/results/session69_prime_ports_power_aliquot_20261026/session69-tournament/REPORT.md`,
aggregates full multi-rooted ordinary uncoloured responses recursively, with
ordered quotient-walk tuples and multiplicities. The present theorem concerns
individual states and named failure/success tests, not that aggregate theorem.

The concept board is: full named guard language; traceable-arc core;
cut-vertex order; aligned cover states; aggregate response counts. The map
`HP(T) -> core -> T` preserves the complete labelled legality predicate;
passing only to counts, or forgetting operation labels, destroys that
predicate. No LRC, prime-classification, or efficient atlas consequence follows.

## 2. The core and the primary structural input

For any directed graph `D`, define its traceable-arc core `U(D)` to have the
same labelled vertices and exactly the arcs used by at least one directed
spanning Hamilton word. Then

    HP(U(D)) = HP(D).                                      (G1)

Indeed `U(D)` is a subgraph of `D`, and by construction contains every arc of
every member of `HP(D)`. Thus the complete guard language and its directed
edge union contain **exactly the same information**. The core is not a lossy
summary of the language. What can fail is naïvely identifying the core with
the original tournament.

**CITED input (BJR).** If a strong tournament `T` has a nontraceable arc
`x->y`, there is a vertex `z` such that the ordered strong components of
`T-z` are `C1,...,Ck`, where `k>=4`, with `x` in `C1` and `y` in `Ck`.
Only this portion of the published theorem is needed below.

Primary source: Arthur H. Busch, Michael S. Jacobson, K. B. Reid,
[*On arc-traceable tournaments*, Journal of Graph Theory 53 (2006),
157–166](https://doi.org/10.1002/jgt.20173).
The statement was independently visually checked in the authors’ primary
restatement: Busch–Jacobson, *Extremal results on arc-traceable tournaments*,
JCMCC 63 (2007), 3–15,
[Theorem 1.1, printed pages 4–5](https://combinatorialpress.com/article/jcmcc/Volume%20063/vol-063-paper%201.pdf).
The cited theorem is not represented as new here.

We also use the elementary Hamilton-cycle fact: a strong tournament of size
at least three has a directed spanning cycle. For completeness, choose a
longest cycle. A vertex seeing it in both directions inserts between two
successive cycle vertices. Otherwise every outside vertex is uniform to the
cycle; strong connectivity gives an outside path from the vertices dominated
by the cycle to the vertices dominating it, and the first cross arc inserts
two outside vertices. Both contradict maximality unless the cycle spans.
A strong tournament has a cycle initially; every shortest cycle is a triangle.

## 3. The essential lemma: two traceable detours

**Lemma G2.** In a strong tournament, every nontraceable arc `x->y` has two
internally vertex-disjoint directed paths of length two from `x` to `y`
using only traceable arcs.

*Proof.* Use BJR and select `a` in `C2`, `b` in `C(k-1)`; these are distinct
because `k>=4`. Every vertex in an interior component of `T-z` is nonseparating
in `T`. To see this, delete such a vertex `v`. The unchanged first and last
components, together with `z`, remain strongly connected: strongness of `T`
forces an arc from `z` into `C1` and an arc from `Ck` into `z`, while all arcs
from `C1` to `Ck` point forward. Every other surviving vertex is reached from
`C1` and reaches `Ck`. This argument does not assume that the damaged interior
component itself stays strong.

If `v` is nonseparating, every arc incident with `v` is traceable. For
`v->w`, rotate a spanning cycle of `T-v` to a Hamilton path beginning at `w`
and prepend `v`. For `w->v`, use a path of `T-v` ending at `w` and append
`v`. Here `T-v` has at least three vertices in the application.
Thus the four arcs in `x->a->y` and `x->b->y` are traceable. ∎

**Corollary G3.** If `T` is strong, `U(T)-v` and `T-v` have exactly the same
labelled reachability relation, for every vertex `v`. Also `U(T)` and `T`
have the same reachability relation.

*Proof.* One inclusion is automatic from subgraph containment. Replace each
nontraceable arc of a path by one of G2’s two detours. At most one internal
detour vertex is deleted, so a detour avoiding `v` exists. This produces a
walk in the core; deleting repetitions gives reachability. Endpoints of a
path in `T-v` are never the deleted vertex. ∎

The strongness hypothesis matters. For the transitive path `0->1->2`, the
core omits `0->2`; after deleting `1`, core reachability is empty although
`T-1` retains `0->2`.

## 4. Reconstruction and explicit decoder

**Theorem G4 (strong case).** A strong tournament is uniquely determined by
its labelled traceable-arc core.

*Proof and decoder.* Keep every arc already in the core. For each missing
unordered pair `{x,y}`, inspect `U-v` for all `v` other than the endpoints.
Choose a `v` for which exactly one of `x` reaching `y` and `y` reaching `x`
holds. Orient the missing pair in that one-way reachable direction.

Such a `v` exists: apply BJR to the actual missing arc and take its `z`.
Its endpoints lie in the first and last components of `T-z`, so G3 gives
exactly the required one-way reachability in `U-z`. The direction is sound
for **every** eligible `v`: the true tournament arc survives deleting `v`,
so its direction must remain reachable in the core by G3. Hence two
eligible separators cannot prescribe opposite directions.
This recovers the exact labelled tournament without knowing it in advance. ∎

For a uniqueness proof phrased with another candidate `T'`, note that a
strong core forces `T'` to be strong. G3 applies to both candidates, so the
same deleted-vertex reachability forces every missing orientation in both.
It is not enough merely to assume that an arbitrary completion of the core
has the same deleted reachability.

**Theorem G5 (all tournaments).** The same reconstruction holds for every
finite tournament, including the singleton and two-vertex boundaries.

*Proof.* Let `T=A1->...->Ar` be its canonical ordered strong-component
factorization. Every Hamilton word traverses the components in this order
and contiguously. Every local Hamilton word extends to a global one, so the
core restricted to `Ai` is `U(Ai)`. A nontrivial strong component has a
spanning cycle; every one of its vertices can begin or end a local Hamilton
word. Thus every cross arc between consecutive components appears in some
global Hamilton word, and no cross arc between nonconsecutive components
does. Consequently the core has exactly the same SCCs, and its condensation
is the directed path `A1->...->Ar`. Recover that order, restore all forward
cross arcs, and apply G4 inside each nonsingleton strong component.
Singleton blocks need no internal reconstruction. ∎

A straightforward implementation computes reachability after each vertex
deletion in each strong block. With ordinary transitive closure this costs
`O(n^4)` time and `O(n^3)` space if all closures are retained, reducible in
space by processing deletions successively. This is the cost **given the
core**. It is not the cost of constructing the core from an implicit
membership oracle or enumerating `n!` guard names.

## 5. Exact continuation minimality, not an overstatement about aggregation

Fix the positions `0,...,n-1` and the identity directed path. `X_n` consists
of the tournaments containing that path. For every permutation `p`, the
named operation `R_p` succeeds exactly when `p` is a directed Hamilton word;
a failure has a distinct observable symbol. On success it reads the same
tournament in order `p`.

By G1/G5, distinct `x,y` have distinct success/failure vectors. There is
therefore a **one-operation** distinguishing test. Any deterministic
quotient preserving these guards must be injective. Conversely storing the
`M` free arc bits gives the state and all successors. Hence the smallest
closed quotient for this protocol has exactly `2^M` states, not merely at
most that many. At `n=1,2`, `M=0` and the single state is correct.

This conclusion is equivariant: relabelling all vertex names in the
language relabels the reconstructed tournament. If language is given only
up to one common relabelling, it determines the tournament only up to that
relabelling. Erasing permutation names altogether is a different observer.
The result asserts neither a small fixed test bank nor an optimal oracle
query algorithm; a binary query decision tree needs at least `M` queries
in the worst case, but attaining that bound is not proved here.

In particular an aligned source/target cover carrier may remove redundant
stored path words, but at arity one it cannot identify distinct literal
presentations while preserving the full named guards. A canonical encoding
of all `2^M` states is possible; a genuinely smaller exact quotient is not.
This is compatible with Session69 P3’s response coefficients and the incoming
fresh-re-marking multiset protocol: those observers retain different data.

**Corollary G7 (pointwise observation algebra).** Over any field `F`, the
unital pointwise algebra generated by the Boolean named-permutation guards
on `X_n` is the entire function algebra `F^(X_n)`, of dimension `2^M`.
For a chosen state `x`, take the product over all permutation names of the
guard itself when it is true at `x`, and of `1-guard` otherwise. Separation
makes this product exactly the point indicator of `x`; these indicators
span all functions. This concerns closure under pointwise products, not the
linear span of the raw guards. At order six the latter, even with a constant,
has at most `6!+1=721` generators for a `1024`-dimensional space and cannot
be all functions. Compact circuits, implicit descriptions, and different
aggregate response algebras are not ruled out.

## 6. Sharp one-step elementary hostile and the repaired dynamic statement

An elementary re-marking deletes one vertex from its present position and
reinserts it elsewhere, preserving the order of every other vertex. Adjacent
swaps can never be legal in a tournament, and doing nothing supplies only a
constant guard.

At order four, list free pairs `(0,2),(0,3),(1,3)` with bit one pointing
forward. Masks `0` and `5` have rows

    x = (2,4,9,3),      y = (6,12,8,1).

Every named elementary nonadjacent relocation has the same legality at
these two states. However `p=(1,3,0,2)` is illegal at `x` and legal at `y`.
The complete six elementary guards and exhaustive partition counts through
order four are in `normal.json`. Exhaustion of the one/two/three-vertex
universes proves four is the smallest possible failure.

Two elementary steps already separate them:

    p1 = (1,2,3,0),    p2 = (0,2,3,1).

The first succeeds at both. Their new rows are respectively
`(2,12,9,1)` and `(6,4,8,3)`; the second fails/succeeds. This two-step witness
is shortest because the initial one-step signatures agree.

**Theorem G6 (accessible elementary guards).** If the observer retains
success/failure for every finite word of named elementary relocation
operations, its smallest quotient is again the identity on `X_n`.

*Proof.* The inherited theorem
`session67_metagraph_enumeration_information_20261026/session67-root/inputs/THEOREM.txt`
proves that the labelled Hamilton words of every tournament are connected
by legal single-vertex relocations. Suppose `p` is legal in `x` but illegal
in `y`. Connect the identity word to `p` using these legal relocations in
`x`. Express each relocation in the positions of its current word. If the
same sequence succeeded throughout in `y`, the composite positional
permutation would be the same `p`, contradicting its illegality in `y`.
Thus some finite prefix has different success/failure. G5 supplies a
separating `p` for every distinct pair, choosing whichever side admits it. ∎

The theorem concerns **guard-labelled accessible dynamics**, not an
unlabelled abstract transition graph. A simple path in the labelled
Hamilton-word graph gives a finite length bound `H(T)-1` on the chosen
side; this is not a polynomial bound. Neither this statement nor G5 needs
literal complement lines, reflection colour, or a chosen tournament-class
output.

## 7. Exact checks, negative control, and remaining boundary

`check.py` is a portable standard-library verifier importing no producer
code. It reconstructs every fixed-path presentation through order seven,
every labelled tournament through order five, and 40 seeded samples at each
order eight through eleven. DFS builds the core from complete Hamilton words;
an independent subset-endpoint dynamic program cross-checks the core in the
explicit receipt’s smaller and seeded universes. It also checks every
one-vertex-deleted reachability relation in each tested strong tournament,
core idempotence in the independent controls, odd Hamilton count, SCC order,
and the order-four hostile.

The five-vertex strong ordinary-mask-10 control has rows
`(20,1,2,7,14)` and a nontraceable arc `4->1`. The decoder restores that arc;
the witness prevents silently replacing reconstruction by a false
arc-traceability assertion.

From this directory:

    python3 -I -B verify.py

The verifier checks the exact file/hash universe, replays `check.py` in
normal and optimized modes, compares both with the frozen receipt, and
requires an intentional decoded-arc corruption to fail in both modes.
All mathematical checks are explicit exceptions, not removable `assert`s.

Closed: the all-order labelled guard question and exact state-cardinality
minimality for its named-operation protocol. Still open here: optimal
query complexity; economical encodings for weaker aggregate/freshly
resampled observers; unrestricted prime-atom classification; efficient
full metagraph enumeration. The cited nontraceable-arc theorem remains an
explicit proof dependency; finite replay does not independently prove it.
