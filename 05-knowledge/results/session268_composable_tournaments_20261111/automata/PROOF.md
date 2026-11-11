# From invisible tournament computers to visible interacting cellular automata

**Status: PROVED elementary family and compiler theorems; FINITE-EXACT seed data.**
This is an infinite, intrinsically recognizable family of tournament isomorphism
classes, not a census of all tournaments. Chan–Pak remains open. No computational
universality is inferred from a finite Rule 110 trajectory. The CA rule is
supplied; an undirected arc-flip graph does not spontaneously select a rule.

## Inheritance and the missing coordinate

The closest recovered mechanism is the literal five-vertex invisible bit of
`05-knowledge/results/session168_computational_structures_20261103/turing/NOTE.md`,
§1: masks 8 and 10 differ in one arc and share the complete unordered path-cover
profile `(9,33,31,10,1)`. That source compiles actual Turing/FRACTRAN steps but
keeps the full profile invariant. Its companion `cellular/NOTE.md` lifts
ternary reversible gates to actual tournament arcs. We inherit those constructions,
not a claim that their unmarked graphs are autonomous computers.

The closest hostile example is precisely an invisible bit: arbitrary complicated
internal evolution still has constant Hamilton count. The corrected near miss is
calling an ordinal sum of invisible cells an H-interacting machine. The least-used
sidecar is a *rooted two-path cover*, retained below. The sibling
`../modular/PROOF.md`, §1, formula(2) and its Contextual congruence corollary,
proves the full ordered-cover substitution interface. Its immediate consequence
is a no-go: replacing a leaf by any other tournament with the same complete
path-cover profile preserves the profile in every iterated tournament-substitution
context. Pure substitution cannot reveal these bits. This follows by induction
on the context tree, and requires no independence assumption about the other
children. Ordered counts are `q_j=j!p_j`; our seed tuple uses unordered `p_j`.

Live concepts: invisible computation; rooted path forests; deliberate module
breaking; pair-word automata; exact class-versus-observable compression. The move
below changes one external arc, so the two five-vertex blocks are deliberately
no longer modules. This is what lets information leave the old quotient kernel.

## 1. Exact seam law for arbitrary modules

For a nonempty tournament A and vertex a write:

* `H(A)` for its number of directed Hamilton paths;
* `e_A(a)` for paths ending at a and `s_A(a)` for paths starting at a;
* `L_A(a)` for ordered pairs `(P,Q)` of vertex-disjoint directed paths covering
  A, with Q nonempty and starting at a and P possibly empty;
* `R_A(a)` for ordered pairs `(P,Q)` covering A, with P nonempty and ending at a
  and Q possibly empty.

The empty path has one realization. Let `A ⋈_(a,b) B` be the ordinal sum A⊕B
with its single cross arc a→b reversed. Then

**Theorem 1 (one-seam interaction).**

`H(A ⋈_(a,b) B) = H(A)H(B) - e_A(a)s_B(b) + L_A(a)R_B(b).`

**Proof.** Paths not using the new b→a arc are exactly the old ordinal-sum paths
except those whose unique cross transition was a→b. Their count is the first two
terms. A path using b→a uses it exactly once and has a unique run pattern
`P R Q S`, where P,Q are in A, R,S in B, Q starts at a, R ends at b, and P,S
may be empty. No other B→A arc exists, so there can be no other run pattern.
Conversely every such pair of covers concatenates legally: P→R and Q→S,
when present, cannot be the missing a→b because a belongs to Q and b to R.
This is a bijection with the Cartesian product counted by L_A(a)R_B(b). ∎

For two binary families A_x,B_y with constant H individually, the mixed
Hamilton interaction is

`H11-H10-H01+H00 = (L_A1-L_A0)(R_B1-R_B0)
                 - (e_A1-e_A0)(s_B1-s_B0).`

Thus full unrooted cover equality does not force zero coupling. The sidecar is
explicitly rooted and therefore is not covered by the substitution congruence.
This theorem also holds for nonempty directed graphs with all cross arcs
initially forward; strong connectivity is not needed for the identity.

## 2. Four strong, intrinsically distinguishable interacting cell states

Use vertices 0,…,4 and lexicographic unordered pairs 01,02,03,04,12,13,14,23,24,34.
A 1 bit means the smaller vertex points to the larger. Let B_0 have mask8 and
B_1 mask10. They differ only on arc {0,2}. Define

`G_xy = B_x ⋈_(0,0) B_y`,  x,y∈{0,1}.

Direct subset DP and the independently derived seam formula give:

| xy | H(G_xy) |
| --- | ---: |
| 00 | 215 |
| 01 | 179 |
| 10 | 245 |
| 11 | 201 |

The receipt includes the complete rooted interface values. The additive mixed
interaction is −8. The multiplicative interaction is
`215·201−179·245 = −640`, so these weights do not factor as an independent
left weight times an independent right weight either. This is a genuinely
interacting two-state response, not just different unary labels.

**Lemma 2.** All four G states are strong, pairwise nonisomorphic, and changing
one bit changes exactly one literal arc within its corresponding five-vertex
half.

**Proof.** Each B state is strong: both contain
`0→4→3→2→1→0`. Each half can therefore reach all its own vertices. The seam
provides a B→A transition, and many untouched cross arcs provide A→B
transitions, so the union is strong. The four distinct Hamilton counts force
four distinct isomorphism classes. The seam is fixed, and the single internal
bit edge is the only changing arc. ∎

## 3. An intrinsic infinite family with exponentially many classes

For N≥1 and a binary word x=(x_0,…,x_(N−1)), indexed periodically, put

`T_N(x) = G_(x_0,x_1) ⊕ G_(x_1,x_2) ⊕ … ⊕ G_(x_(N−1),x_0).`

This has 10N vertices. A tournament is in the family iff its *ordered strong
components* are N members of the four G classes, whose decoded pairs have
second coordinate equal to the next pair's first coordinate, including the
last-to-first check. This is an intrinsic definition with no marked vertices.
The cyclic adjacency is a relation on the uniquely ordered component positions;
it is not a quotient of that linear order by rotation.

**Theorem 3 (exact class census).** The map x↦[T_N(x)] is injective. Consequently
the family has exactly 2^N tournament isomorphism classes. Each has

`H(T_N(x)) = ∏_(i=0)^(N−1) w_(x_i,x_(i+1)),`

where `w00=215, w01=179, w10=245, w11=201`.

**Proof.** The G blocks are strong, while every arc between blocks points
forward. They are therefore precisely the SCCs. Every tournament isomorphism
preserves the unique linear SCC order and each SCC's type. Lemma2 decodes each
pair and hence x. Conversely each word constructs a member. Hamilton paths in
an ordinal sum traverse its components in order, giving the product formula. ∎

Rotating x generally gives a different isomorphism class: it rotates the
*sequence* of component types, which the linear SCC order does not identify.
A periodic word may equal its own rotation; it was already the same word.
We count indexed words, not binary necklaces.

### Automorphisms and labeled mass

All four G states are rigid. The receipt checks every permutation preserving
outdegree, at most8,24,48,144 candidates respectively; all other permutations
are impossible automorphisms. Only the identity survives in each case. This
is an exhaustive finite certificate, not a heuristic graph-isomorphism guess.
Every T_N(x) is consequently rigid as well: an automorphism must preserve each
position in the SCC order and is identity inside each rigid component. Thus
the family contains exactly `(10N)! 2^N` labeled tournaments on a fixed10N-set.
An H-bin of ordinary multiplicity M has labeled multiplicity `(10N)! M`.
The normalized reciprocal-automorphism mass equals the ordinary class count.
Also the labeled sum of Hamilton counts is `(10N)! trace(W^N)` with the
matrix W specified below. These labeled statements concern all labelings of
this intrinsic family, not only the pinned10-vertex-block representatives.

## 4. A quadratic-size exact Hamilton spectrum and multiplicities

Write n0,n1 for the counts of 0 and 1. A nonconstant circular word has k≥1
0-runs and k 1-runs. It has n0−k pairs00, n1−k pairs11 and k pairs of each
of 01 and10. Thus

`h(n0,k)=215^(n0−k) 201^(N−n0−k) (179·245)^k`,

for `1≤n0≤N−1`, `1≤k≤min(n0,N−n0)`. The two constant words have H=215^N
and H=201^N.

**Theorem 4 (complete observable census).** These values are all distinct,
and their exact class multiplicities are

`M(N;n0,k) = (N/k) binom(n0−1,k−1) binom(N−n0−1,k−1).`

The spectrum contains exactly `floor(N²/4)+2` values, including N=1 and N=2.

**Proof.** Factor `215=5·43`, `201=3·67`, `245=5·7²`; 179 is prime.
In h the valuations at 43,67,179 are respectively n0−k,n1−k,k. Hence h
recovers all three coordinates. Constant words have k=0 and only one of the
43/67 exponents nonzero, so cannot coincide with any nonconstant case or each
other. To count words, distinguish a start of a zero-run. Choose its indexed
position in N ways and choose ordered positive compositions of n0 and n1 into
k parts. Each word is counted k times, once per zero-run start, giving M.
Finally
`Σ_(n0=1)^(N−1) min(n0,N−n0)=floor(N²/4)`.
For N=1 the sum is empty, leaving the two constant words; for N=2 there is
one mixed case with multiplicity2 and two constant cases. ∎

Precompute binomial coefficients and the three needed power arrays. The full
weighted spectrum then takes O(N²) arithmetic operations and O(N²) records,
versus enumerating 2^N distinct classes. Integers have O(N) bits, so this is an
arithmetic-operation claim, not unit-cost bit arithmetic. With a multiplication
cost M(N), a conservative output-sensitive bound is O(N² M(N)) bit operations
up to logarithmic factors for divisions. A single specified word's H takes
O(N) pair-counting operations followed by three exponentiations; a complete
explicit list of all classes still necessarily has exponential output size.

* N=100: 1000 vertices, `2^100=1267650600228229401496703205376` classes,
  exactly 2502 H values. The receipt computes all multiplicities and checks their
  sum against2^100 without enumerating those classes.
* N=1000: 10000 vertices, 2^1000 classes, exactly 250002 H values. This is a
  theorem consequence, not a claim that the receipt enumerated that output.

The 2×2 weighted transfer matrix `W=(w_ab)` gives the further exact identity
`Σ_x H(T_N(x)) = trace(W^N)`. Replacing each w_ab by a formal marker records
pair statistics. The closed multiplicity formula is the compressed coefficient
extraction of this literal word automaton. Hamilton count forgets arrangement
but the ordered SCC word retains it; no two arrangements are silently merged
in the ordinary class census.

## 5. Literal cellular dynamics and an exact local response law

Let f:{0,1}^(2r+1)→{0,1} be any radius-r binary cellular rule on the indexed
N-cycle. Apply it to x, then use the same construction T_N. Theorem3 gives a
well-defined map on these *ordinary tournament isomorphism classes*, conjugate
to the binary CA on this finite ring. Decode pair states from the ordered SCCs,
apply f, and rebuild the SCC types. The alphabet is four actual strong tournament
classes, with the overlap constraint above.

For radius1, the new pair at position i is
`(f(x_(i−1),x_i,x_(i+1)), f(x_i,x_(i+1),x_(i+2)))`.
The three neighboring pair states determine this, so it is literally a radius1
CA on the four-state SCC alphabet restricted to its consistent subshift.
For general r the analogous pair update has radius r in that alphabet.

**Theorem 5 (arc cost, light cone, exact response).**

1. If x and y differ at d positions, the pinned tournaments T_N(x),T_N(y)
   differ in exactly 2d arcs. A differing bit changes its two stored copies,
   one in each neighboring G block (for N=1 these are the two halves of the
   same block). Thus one CA step has a literal path of2d single-arc reversals.
   Intermediate pinned tournaments are allowed to violate the pair-overlap
   constraint; a synchronous step is not a single arc-flip edge.
2. If initial words differ only in a set S, after t steps their different bit
   positions lie in the cyclic distance-rt neighborhood of S. Different G
   components lie in that set and its predecessor. In particular, a single
   initial defect affects at most min(N,2rt+2) G positions. This is a true
   finite propagation law in SCC-position geometry, not sparse tournament
   graph distance (the underlying tournament is complete).
3. For N≥2, flipping bit x_j to b=1−x_j changes H by the exact factor
   `w_(x_(j−1),b) w_(b,x_(j+1)) /
    (w_(x_(j−1),x_j) w_(x_j,x_(j+1)))`.
   N=1 instead has ratio w_bb/w_(x_j,x_j). Hence a defect's response depends
   on its neighbors and can be maintained with two factor updates. For N≥2,
   d changed bits affect at most min(N,2d) pair factors. A whole CA round
   costs O(N) local decisions; updating only an already located defect cone
   costs O(its affected factor count) arithmetic operations.

**Proof.** (1) follows from Lemma2 and disjoint stored copies. (2) is induction
on t using the local dependency neighborhoods, then one extra predecessor for
the pair code. (3) cancels all unchanged factors in Theorem3. ∎

The integer H is a global readout: it can change after one local move. The
finite-speed theorem concerns where the encoded state changes, not delayed
propagation of a global numeric observable. A global simultaneous update need
not preserve H or its ordering. Finite rings are finite dynamical systems;
we make no claim of undecidability for a fixed ring. Any finite spacetime
window on an infinite binary CA is reproduced if its dependency cone does not
meet the chosen finite ring's wrap boundary. This is the precise bridge to
longer tape simulations, not a universality proof.

## 6. The scalar spectrum is not a closed dynamical state

The quadratic spectrum compression is exact for the *static* census, but does
not make H a sufficient state for subsequent CA evolution. Under Rule 110,
width6 words `101000` and `100100` both have
`H=88902740880625`, but after one step their H values are respectively
`76567725723825` and `77701668671025`. They have the same pair statistics;
the missing coordinate is longer-word arrangement. The receipt reproduces
these values. This is a new obstruction to replacing the literal class state
by its scalar observation, not a failure of the static census theorem.

**Proposition 6 (finite-horizon interaction hierarchy).** For a radius-r rule,
the length-(2rt+2) cyclic block inventory of x determines H after t steps.

**Proof.** The pair of bits at positions i,i+1 after t steps is a deterministic
function of the initial interval from i−rt to i+1+rt, of length2rt+2.
Tabulate this local function and sum initial block multiplicities to recover
the four final pair counts, then use Theorem3. This remains valid when the
block length exceeds N: read the periodically repeated word. ∎

Thus the state interfaces grow with prediction horizon. At t=0 the length2
inventory suffices; at t=1 for radius1 length4 suffices. This is a concrete
connection between finite propagation, interacting word automata, and the
failure of Hamilton-count-only lumping. It does not claim that the listed
block inventory is minimal or that no other uniform compression can exist.

## 7. Receipts, scope, and next question

Run from this directory:

    python3 probe.py
    python3 interacting_ca.py

`probe.py` computes all 25 one-seam choices on the two invisible bits.
`interacting_ca.py` checks both seed profiles, strong connectivity and seam
formula for the chosen four states, all 2^N words for1≤N≤12 against the closed
census, a24-cell Rule 110 trace, finite cones for five rules/all 128 inputs on a
7-cycle, and the polynomial100-cell inventory. These are targeted checks of
new objects, not a wholesale audit of inherited research.

This family makes progress on exact structural enumeration and composable
Hamilton evaluation, not a proof of the global Chan–Pak spectrum conjecture.
It supplies explicit arbitrarily large interacting CA families, an exact
finite-speed law, ordinary-isomorphism classification, and an exponential to
quadratic compression of their Hamilton-value distribution.

**OPEN:** Can rooted interfaces produce a similarly compact complete spectrum
for *strong* interaction chains without ordinal SCC factorization? A single
seam already makes a strong10-vertex cell, but this family is not globally
strong when N>1. A fixed strong wrapper preserves much computational structure
but need not preserve the simple product spectrum. Its rooted-cover response
must be derived, not inferred from the present formula.
