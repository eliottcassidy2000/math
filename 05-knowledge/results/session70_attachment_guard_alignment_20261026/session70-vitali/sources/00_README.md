# Strong-atom responses: relative alignment is unbounded continuation memory

2026-10-26. **PROVED** elementary statements C1-C5 below; **FINITE-EXACT**
complete joint-port responses through order seven and all run arities through
order five. **INHERITED PROVED** agreement-path reconstruction and the
operation-kernel criterion. This is an exact finite atlas and an all-order
closed carrier, not an efficient classification of all prime tournaments.
Local theorem labels belong to this file. No literature-priority claim.

The next compression barrier is stronger than a missing endpoint: for every
`k>=1`, there are `2^k` strong presentations at order `7k+2` with the same
two rooted endpoint tournaments and colour, distinguished by legal re-marking
guards. A fixed number of extra alignment bits cannot repair that quotient.
An admissible relative alignment does repair it. The full path can then be
reconstructed, rather than separately stored.

## 1. Inheritance, terminology and operations

Closest mechanism: [joint-port R1-R4](../rooted_lines_frobenius_sheets_20261026/README.md)
and the older [Session25 agreement and run theorems 1-3](../session25_spectrum_wiring_20261021/session25-wiring/RESULT.md).
Canonical hostile: the same labelled cyclic triples with automorphism orders
81 and 27. Corrected near miss: independently canonicalizing two endpoints
forgets their common vertex alignment. Least-used sidecars: the relative-frame
double coset and the legal-operation guard, recovered in [HISTORY.md](HISTORY.md).

The final integration also inherits [Session69 P2-P3](../session69_prime_ports_power_aliquot_20261026/session69-tournament/REPORT.md):
all-arity multi-rooted response multiplicities recursively compose ordinary,
uncoloured target incidence along ordered quotient-walk tuples. This is a
genuine aggregation theorem. Our continuation question additionally retains
individual choices, colour and named re-marking guards. C1 classifies the
individual alignment orbits inside those aggregate response coefficients;
C3-C5 measure what is lost when only their endpoint types are retained.

An **ordered-join atom** means a strongly connected tournament (with the
singleton as a formal unit). A **modular prime** has no proper non-singleton
module. These are different: the unique strong four-vertex tournament is
not modular prime. The atlas records both predicates.

Let `X_n` be tournaments on positions `0,...,n-1` containing the path in that
order. A free-bit mask lists pairs `i<j-1` lexicographically, least significant
bit first, with bit one meaning `i->j`. There are `2^M` states,
`M=(n-1)(n-2)/2`. The actual complement `tau` reverses the free arcs. A
re-marking `R_p` sends `x` to the same tournament read in permutation order
`p`, **only if p is a directed Hamilton path**; otherwise its result is the
distinguished failure state. Reflection `sigma` is converse followed by
position reversal. For `n>=3`, a literal line is blue iff `sigma(x)=x`.
Orders one and two are formal composition inputs, not literal paired lines.

The continuation statements retain the names of permutation operations.
Uniformly resampling an unspecified path is a different, weaker protocol.
Adding attachments to this operation menu cannot remove its obstructions.

## 2. C1: the all-arity response is an admissible double-coset count

Fix a source tournament A. An ordered spanning path cover
`F=(P_1,...,P_k)` has nonempty, disjoint directed paths covering every
vertex; singleton paths are allowed. Put

    B_F = A^op with every arc of F restored,
    r_F = (start(P_1),end(P_1),...,start(P_k),end(P_k)).

Retain the **joint** rooted isomorphism type `h=[B_F;r_F]`, including repeated
root names when a component is a singleton. Define `R_A^(k)(h)` to be the
number of `Aut(A)`-orbits of these ordered covers. The action is free: an
automorphism fixing an ordered list of directed paths fixes every vertex.
For k=1 this is exactly the previous F response. The blue subcount is its B
response. For k>1 colour is tested after gluing the full global path.

**C1 (PROVED).** Fix a representative `h=[B;r]`. For each bijection
`f:V(A)->V(B)`, form the agreement digraph

    D_f = {u->v in A : f(u)->f(v) in B}.

Call f admissible when D_f is the ordered spanning path forest whose
ordered endpoint pairs are `f^-1(r)`. Let `N_Ah` count admissible bijections.
Then

    R_A^(k)(h) = N_Ah / (|Aut(A)| |Aut(B;r)|).                 (C1)

More precisely, the individual cover orbits are the admissible relative
double cosets `Aut(B;r) \ {f} / Aut(A)`.

**Proof.** An admissible f recovers its cover by following the agreement
arcs from each specified start to its specified end. The paths must cover
all vertices and have no detached cycles. Conversely, for every cover
landing at h, its identifications with (B;r) form a torsor of size
`|Aut(B;r)|`. Quotienting the free source action proves the formula. The
two-sided action on f is free as well: a stabilizing pair induces an
automorphism of A preserving the ordered forest and hence every vertex.
Thus each admissible double coset has exactly the displayed denominator.
This is Session25's k=1 agreement argument extended to ordered forests.

**Closed carrier.** Store A, B on one common vertex set and the ordered
forest endpoint list, modulo **simultaneous** isomorphism; equivalently use
the two models and their admissible relative alignment orbit. Agreement
reconstructs all internal paths. Forest complement swaps A and B, retaining
the forest. Re-marking selects another legal cover of the current A and
constructs its new B. Substitution assembles these recovered paths in the
quotient walk and restores exactly the inter-run hinges. Thus the carrier
is closed under repeated complement, re-marking and the specified assembly.

This removes redundant storage of internal path words **only because their
source-target alignment is retained**. Separately storing A's and B's
isomorphism types does not permit this reconstruction. At k=1 the agreement
path itself orders every vertex; the aligned-pair carrier is equivalent to
the original literal presentation, not a new smaller family of tournaments.

In particular no canonical physical path has been selected equivariantly
from an unmarked symmetric tournament. Canonical encodings of whole orbits
and selecting one vertex from a nontrivial orbit are different operations.

## 3. Complete finite strong-atom responses

[strong-atlas.json](strong-atlas.json) contains every source representative,
H, automorphism order, modular-prime flag, and every joint-rooted target
with its total and blue multiplicities. [run-atlas.json](run-atlas.json)
contains every ordered-cover response at every arity for sources of size
at most five. There are no sampling or symmetry filters beyond exact
isomorphism and the explicitly stated strong-source restriction.

| n | Strong classes | Modular-prime classes | Strong presentations | Largest k=1 target fibre |
|---:|---:|---:|---:|---:|
|3|1|1|1|1|
|4|1|0|5|1|
|5|6|3|50|2|
|6|35|15|903|3|
|7|353|197|30773|8|

The singleton adds one formal source; order two adds no strong source.
The all-arity atlas has 3,967 raw ordered covers and 2,973 response states.
The independent audit checks every coefficient in that atlas by (C1),
using alignment enumeration rather than the producer's cut enumeration.

For a specified A and h, the missing exact discriminator has
`R_A^(k)(h)` possible values. An alignment-orbit index is a sufficient and
cardinality-minimal discriminator for **faithful cover classification**;
its fixed-length binary cost is `ceil(log2 R_A^(k)(h))`. At k=1 the largest
cost in the above finite universes is respectively 0,0,1,2,3 bits. This is
conditional on A and h being given; it is not the total encoding cost or a
claim about every weaker numerical observer.

Already at n=5, masks 25 and 42 have the same unrooted source, rooted target
and colour but different rooted sources. Thus the one-sided response state
fails even a second complement update. Recording both rooted endpoint
types makes complement a swap, but does not close re-marking.

## 4. C2: two outer ports fail inside a revisited child at order six

Take `T=C3[A,1,1]` with the strong child A having outgoing bit rows
`(2,4,9,3)`. Label A by 0,1,2,3 and the singleton children by 4,5.
The source rows are `(18,20,25,19,32,15)`. Compare

    P=(0,1,2 | 4 | 5 | 3),
    Q=(1,2 | 4 | 5 | 3,0).

Both use the quotient walk `(A,4,5,A)` and two child runs. Their child
complements are, respectively,

    (14,12,0,4), with outer roots (0,3),
    (4,13,0,5),  with outer roots (1,0).

These are the same jointly rooted tournament: both are transitive, with
outer-root outdegrees 3 and 1. But the full ordered port lists are
`(0,2,3,3)` and `(1,2,3,0)`. The singleton coincidence and the inner port
positions were lost. Their full partners have H=15 and H=23, respectively;
both presentations are black. These values are also checked by all 6!
vertex orders. Thus source type, quotient run word, run count, complemented
child relation and its two outer roots do **not** suffice.

This is a six-vertex hostile for the specified compression, not a claim
that six is universally minimal under every possible attachment interface.

## 5. C3: even both joint-rooted endpoints fail, then amplify

**C3 (PROVED by explicit witness).** At n=7 let x,y have free masks
1056,4160 and rows

    x=(2,12,41,17,39,75,31),
    y=(2,20,9,51,37,71,31).

Both sources are strong. Their source joint-rooted types `[T;0,6]` agree,
as do their target joint-rooted types `[tau T;0,6]`. Both are black.
An isomorphism from x to y fixing both roots is
`(0,1,4,2,3,5,6)`; one from tau(x) to tau(y) is
`(0,2,4,1,5,3,6)`. They are different relative alignments.

The re-marking word `p=(0,1,2,3,5,6,4)` is illegal in x and legal in y.
Indeed x lacks `3->5`, while all six required arcs are present in y.
This proves failure of the paired-rooted quotient for the declared guarded
operation. Independent exhaustive root-fixing permutations check both
isomorphisms and the colour/guard assertions.

**C4 (PROVED, unbounded alignment repair).** For each k>=1, put k copies
of the unmarked source of C3 in an ordered sum S, and form `C3[S,1,1]`.
In each copy independently choose presentation x or y; traverse S in
component order and then the two singleton children. There are 2^k
resulting presentations at n=7k+2, all of one strong source class, with
the same joint-rooted source and target types and the same black colour.
Each of k re-marking guards reads one independently chosen bit.

**Proof.** Source rooted isomorphisms from C3 extend blockwise and fix the
hinges. Target rooted isomorphisms do the same after complement, since
every inter-copy hinge uses the fixed endpoints. The quotient C3 makes
the whole source strong. Its first and last vertices have outdegrees
`7k-5` and `7k`; their sum differs from `n-1=7k+1`, whereas a reversing
anti-automorphism would force equality. Thus every presentation is black.
Extend p from C3 inside any chosen copy
and keep all other positions in order. All boundary transitions stay
legal by uniform substitution, and the one internal guard is legal
exactly when that copy used y. Thus every bit is separately observable.
The source and target isomorphisms need not be the same, which is exactly
the missing coordinate. This proves the count and the lower bound of k
additional binary bits for any exact deterministic repair on this fibre.
Controls check all 2^k variants for k=1,...,6; the proof is all-order.

These examples are strong ordered-join atoms, not asserted modular primes.
They disprove a constant-size repair of the paired-rooted observer, not
the existence of an integer encoding or a more efficient variable-length
representation. A single integer can of course store arbitrarily many bits.

## 6. C5: the precise smallest continuation quotient

Fix an observation O and a named family of guarded operations, with a
distinct failure symbol. Define

    x ~_O y iff O(w x)=O(w y) for every finite operation word w,
    including the empty word and failed operations.                 (C5)

**C5 (PROVED; inherited criterion applied).** This is the coarsest
operation-stable equivalence refining equality of O. It is obtained on
a finite universe by repeatedly refining each cell by its labelled
successor cells, including guard failure.

Proof: append any one operation before a continuation word to see stability.
Every stable quotient refining O must preserve all words, by induction,
and therefore refine (C5). Finite refinement stops because cells only split.
This is the operation-kernel mechanism of
[THM-840](../../../01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md)
and the older [three-sort continuation minimization](../three_sorted_metagraph_continuation_minimization_codex_S11.out),
with a different operation alphabet; its earlier face-deletion counts are
not reused as results for re-marking.

**FINITE-EXACT.** Through n=7, the complete legal re-marking guard signature
alone distinguishes every x in X_n. Thus for this named-operation protocol,
even before later output observations, the minimum state counts are
`1,1,2,8,64,1024,32768` for n=1,...,7. At seven, the paired-rooted observation
has only 32,272 cells. Its 496 lost states are separated by guards.
The exact fixed-length state cost in these universes is M free bits.

The all-order claim that guards alone distinguish every presentation is
**OPEN here**. Neither the finite census nor the exponential family proves
it. If the observation keeps only distributions under fresh uniform
re-marking, a smaller stochastic quotient may exist; it answers a different
question and must not stand in for the actual named path protocol.

## 7. The 81/27 control, phase and historical synthesis

[FROBENIUS-ABC.md](FROBENIUS-ABC.md) carries the retained ancestry phase,
the exact six-frame S3 bridge, the 4p distinction, and the three arithmetic
maps. [HISTORY.md](HISTORY.md) connects the March-to-October Vitali, moment,
selector, seam, continuation, carry and normalization mechanisms. The common
claim is about a **specified quotient and specified continuations**, not a
universal equivalence of all objects called atoms.

The 12-vertex control is rechecked exactly: identical 67 labelled cyclic
triples, H=462915 at both sources, automorphism orders 81 and 27, and
presentation totals 5715 and 17145. The endpoint-flux histograms differ.
The source mass and relative phase cannot be replaced by the C3 field.

## 8. Reproduction and remaining boundary

Run `python3 -B 05-knowledge/results/strong_atom_continuation_vitali_20261026/verify.py`.
It checks frozen hashes, replays all scripts normally and under -O, and
reconstructs both complete atlas files in memory for byte comparison.
`audit.py` imports no producer code: it tests all vertex orders and all
root-fixing permutations, checks all 29,820 strong response states, and
verifies all 2,973 all-arity coefficients using the agreement criterion.
Its legal-guard audit uses permutation constraints and Boolean subcubes,
independently of the producer's path DFS. There is no separate-agent audit.

The remaining problem is **useful compression of the admissible alignment
and efficient generation of new strong/prime atoms**. The sufficient closed
carrier and its exact finite coefficients are now explicit. No universal
constant-size root/interface, all-order guard-injectivity theorem, or
cross-domain Collatz/ABC consequence is claimed.
