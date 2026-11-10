# Rerooting restores the missing response: a ternary sidecar and a decreasing potential

**PROVED scoped constructions and inverse bounds; FINITE-EXACT seed grammars.**
This second round preserves the first-round packet and pursues its positive
survivor. It does not prove Chan–Pak. The main gain is not a larger finite
cutoff: it is an explicit rerooted family, a proof of its hidden ternary
response, and a target-dependent finite inverse allowing free root changes.

## 1. What multiplication repairs, and what it cannot repair

Ordinal joins satisfy H(A ▷ B)=H(A)H(B), by concatenating Hamilton paths in
the two components. This is inherited from THM-1862, cited in
`01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`.
The THM-4094 reduction is inherited: given the attained carry atoms 49,63,343,
the remaining universal construction obligations may be restricted to p
for odd primes p!=7, and 7p for odd primes p!=3. Endpoint universality does
not independently discharge either obligation, since E alone does not
control the full-count response D+U-E.

For the first-round bank of all seeds of order <=5 followed by retained-root
O/Z tails, exact product closure through 2001 leaves only

    7, 21, 173, 631, 1897.

In particular 105 is repaired by 3*35. Here 173 and 631 are primes and
1897=7*271. This finite list is not a global classification of the grammar.
It illustrates exactly why the two primitive lanes survive multiplicative
repairs. In a multiplicative grammar of attained positive odd counts, a
prime absent from its generator image cannot first appear as a nontrivial
product. Likewise 7p cannot first appear as a product when neither factor7
nor the full generator7p is available. Therefore 173 is a decisive hostile
for multiplication alone, not merely an inconvenient composite.

The missing facts in the gap argument are still physical. A positive
endpoint selector E=N says nothing about whether its source also meets
`D+U-E=h` at a specified primitive target h.

## 2. Exact target-bounded selection also permits J

Recall J sends `(D,E,U,d)` to `(D+U-E,D-E,D+U,D)`.
All three moves O,Z,J strictly increase U on the reachable-root sources in
this packet. A last O can be stripped without affecting the full count;
a last Z or J producing count h has `h=D+U-E>=U` at its predecessor.
Consequently the first-round finite inverse with cutoff U<=h extends
unchanged to arbitrary O/Z/J words. There are finitely many profiles with
this cutoff because every coordinate is a nonnegative integer at most U.
Equivalent profiles may safely be merged for these three operations.

This permits a sharp statement about *all word lengths*, not a depth cap.
The verifier exhaustively constructs all rooted seed profiles from labelled
tournaments of orders2 through b, with the root reaching every vertex:

- b=5: 42 profiles. Arbitrary O/Z/J tails miss exactly 7,21,173 through2001.
- b=6: 264 profiles. Arbitrary O/Z/J tails miss exactly 7,21 through2001.

Thus the minimum initial seed order for obtaining 173 in this particular
O/Z/J grammar is exactly6. This does not assert minimum order8 for a
173-count tournament, nor minimality among unrestricted constructions.
The witness below explains what the sixth seed vertex contributes.

## 3. A six-vertex seed with count15 and a new response59

Start with edge0->1. Letters create vertices in order, beginning at2.
Apply ZOOZ while retaining root0. The resulting graph has adjacency masks

    [26,60,57,48,32,1].

Its root0 profile is `(15,1,23,1)`. Change only the distinguished root to
vertex2, the vertex created by the first Z. Its profile becomes

    (15,1,59,9).

Now apply Z twice at this root:

    (15,1,59,9) -> (73,9,109,9) -> (173,9,209,9).       (1)

The final eight-vertex graph is strong and has adjacency masks

    [218,252,57,240,224,193,132,4].

The standard-library verifier counts its Hamilton paths and all listed
profiles directly from the graph. Rerooting changed no edge and preserved
H=15; it exposed the response59 and deletion count9 needed for the target.
Thus the missing resource was not a new scalar seed count, but a different
physical profile of an already constructible graph.

The analogous smallest repair also handles the first-round one-apex hostile:
edge.ZZ, reroot at vertex2, then Z produces a five-vertex graph of count13.
Before the final move its rerooted profile is `(5,1,9,3)`. The output masks
are `[18,28,9,17,4]`. The output has left the d=1 class, as it must.

## 4. An infinite rerooted family: dyadic count, ternary response

For k>=1, let A_k be edge.Z O^k Z, rerooted at the first Z newborn.
Then its exact profile is

    D = 3*2^k+3,
    E = 1,
    U = 4*3^k+6*2^k-1,
    d = 2^(k+1)+1.                                   (2)

The k=0 profile is `(5,1,9,3)`; the sole exception in the displayed formula
is D, which is5 rather than6. All these sources are strong.
The family demonstrates a specific sidecar inflation: H is of dyadic
scale, whereas the rerooted response contains a ternary term. A bounded
scalar or a fixed-root binary description can hide this larger response.

### Core-slot proof, not a fitted recurrence

Remove the k transitive chain vertices inserted by O. The four remaining
core vertices are a=0,b=1,r=2,z=3, with directed triangle
`a->b->r->a`, and z beating only a. Every inserted chain vertex is beaten
by a,b,r and beats z. Chain vertices are mutually ordered transitively.

In any fixed order of the core vertices along one or two paths, chain
vertices may occupy these slots:

- between a core vertex other than z and z;
- at the end of a path ending in a core vertex other than z;
- at the beginning of an ordinary path beginning with z;
- in an entirely core-empty ordinary path.

There is no prefix slot in the distinguished path, which must start at r.
For endpoint E, the path must end at r, so its last suffix slot is forbidden.
If a core adjacency is illegal and there is no admissible intervening slot,
discard that core order. If such a slot exists, require it to be nonempty.

With s available slots and q required-nonempty slots, the exact number of
assignments of k ordered chain vertices is

    sum_{j=0}^q (-1)^j binom(q,j) (s-j)^k.            (3)

Each chain vertex chooses a slot; within each slot its order is forced.
This is a bijection, and inclusion-exclusion enforces the required slots.
At k=0 use 0^0=1, as required by empty assignments.

For D and E enumerate the24 core permutations, and for d the6 permutations
after deleting r. For U enumerate the24 ordered core-cover configurations:
choose which core vertices accompany r in its distinguished path, order
them after r, and order the complementary core path. Grouping their slot
contributions by the number of core vertices in the distinguished path:

    1 core vertex:   2^k + 2*3^k
    2 core vertices: 2*2^k
    3 core vertices: 2*2^k - 1^k + 3^k
    4 core vertices: 2^k + 3^k.

The resulting exact exponential polynomials are

    D(k) = -0^k + 3*1^k + 3*2^k,
    E(k) = 1^k,
    U(k) = -1^k + 6*2^k + 4*3^k,
    d(k) = 1^k + 2*2^k.

These prove (2), including the k=0 boundary. The verifier independently
recomputes the finite core-slot coefficients, then compares the formulas
against literal graph profiles for k=0,...,7.

### Exact continuing family

For t>=1 copies of Z at the new root, the inherited recurrence gives

    D_t = D+(2^t-1)U-E-(2^t-2)d.

Substitution in (2), for k>=1, gives the strong count family

    H(A_k.Z^t)
      = 4(2^t-1)3^k + (4*2^t+1)2^k + 5-2^(t+1).     (4)

At t=2 this becomes `12*3^k+17*2^k-3`, giving67,173,457 at k=1,2,3.
No primality infinitude or arbitrary-target coverage is claimed.
The novelty here is the explicit recovered sidecar and its counting proof,
not that these individual full counts were previously unknown.

## 5. A decreasing potential for arbitrary root choices

An arbitrary reroot can increase or decrease the designated U. Therefore
one cannot use a single marked response as a monotone pruning coordinate
when allowing free root changes. The correct finite state is the actual
graph, together with its complete menu of reachable-root profiles.

For a nontrivial tournament T having a vertex reaching every vertex, define

    mu(T) = min { U(T,r) : r reaches every vertex of T }.

At each extension choose any such r and append either a sink O or a vertex
beating exactly r (Z). Rerooting itself changes neither T nor mu(T).

**Theorem (PROVED).** Along this freely rerooted grammar:

1. An O move leaves H unchanged and doubles mu.
2. A Z move strictly increases H by a positive even integer and never
   decreases mu.

**Proof.** For O, the eligible old roots remain exactly the same; the new
sink is ineligible. For every old root, removing the new sink from a rooted
two-path cover has exactly two inverse placements, at the end of either
path, including an empty ordinary path. Thus all eligible U values double.

For Z at r, write D=H(T). Always U(T,r)>=D, by splitting each Hamilton path
at r into its distinguished suffix and complementary prefix. Also E(T,r)<D:
since r has an outgoing neighbor, insert r before the first vertex it beats
in any Hamilton path of T-r. The predecessor, if present, beats r by
minimality. This gives a Hamilton path not ending at r. Hence
`H(new)=D+U-E>D`, and odd parity makes the increase at least2.

Moreover `H(new)=D+U-E>=U>=mu(T)`. Every root response in the new graph is
at least H(new), again by splitting paths. The new graph is strong: every
old nonroot reaches the new vertex, which reaches r and hence all old
vertices; r reaches the new vertex via any old outgoing neighbor. Therefore
all new roots are eligible and `mu(new)>=H(new)>=mu(T)`. QED.

**Corollary (PROVED complete target cutoff).** Fix an initial seed T0 with
D0=H(T0) and mu0=mu(T0). Consider any nonempty construction yielding target
h, and strip all trailing O moves. If any moves remain, its last move is Z.
The numbers z and o of Z and O moves then obey

    z <= floor((h-D0)/2),
    o <= floor(log_2(h/mu0)).                       (5)

If h<mu0, no such nonempty reduced construction exists. The empty-word
seed count is handled separately. Indeed each Z raises H by at least2.
Before the last Z, mu is at least `2^o mu0`, while the final count is at
least that preceding mu, proving the second bound.

Thus a finite target-dependent enumeration of graphs and all eligible root
choices is a complete inverse for this grammar. It is not necessary to
assume an arbitrary order cutoff or to mistake a truncated search for an
all-length result. A lexicographic potential `(H,mu)` strictly increases
at every graph extension: Z increases its first coordinate, O increases
its second. In an inverse certificate, deletion of the last Z decreases H,
and deletion of a final sink halves mu.

This is a genuine sufficient *termination* invariant for recursively
producing seed banks, not a proof that every target passes the inverse.
Merging arbitrary graphs solely by H or `(H,mu)` would lose the root menu
and is not justified. Full four-coordinate states suffice only while
future choices are restricted to O/Z/J at the marked root.

## 6. Exact remaining coverage obligation

Combining this grammar with ordinal products gives a concrete, well-founded
route to a universal theorem. For each primitive lane target h from THM-4094,
find a smaller constructed graph T and an eligible root r such that

    U(T,r)-E(T,r) = h-H(T),                          (6)

or provide a base witness. Equation (6) then gives a physical extension
with full count h. Its predecessor has strictly smaller H, so recursion
cannot be justified by circularly assuming the same target. Product nodes
also have strictly smaller positive factors. Carry atoms49,63,343 and the
finite bases must be separately supplied as in THM-4094.

What remains OPEN is the existence of these target-dependent root-menu
intersections for all primitive lane targets. Universal even differences
forget H(T), and endpoint universality forgets U(T,r), so neither theorem
proves (6). The explicit 13 and173 constructions show that free rerooting
can restore precisely this missing sidecar in cases where earlier complete
subfamily inverses failed. They motivate the invariant without claiming
that every target is covered.

## 7. Reproduction, controls, and research-round ledger

Run from this directory:

    python3 -I -B verify_round2.py > ROUND2.replay.json
    cmp ROUND2.replay.json ROUND2.json
    python3 -I -B -O verify_round2.py > ROUND2.optimized.json
    cmp ROUND2.optimized.json ROUND2.json

Only standard-library modules are used, plus the first-round sibling loaded
by an explicit local path. Checks survive optimization. The verifier:

- exhaustively computes all rooted seed profiles of orders2 through6;
- compares its optimized seed implementation with the first-round order<=5
  implementation;
- completely explores the specified O/Z/J grammar up to response2001,
  with the proof in §2 supplying all-word-length completeness for target
  values at most2001;
- computes product closure and the 105 factorization;
- verifies all core-slot coefficients exactly,8 literal profiles, and15
  literal continuing-family counts;
- recounts both rerooted hostile repairs and strongness of the173 output;
- checks the mu transformation on O and Z at42 chosen representative seed
  graphs/roots,84 checks total. These are controls for the proof, not an
  exhaustive census of every graph sharing one of those profiles.

Round sequence: multiplication repairs105 but preserves prime hostiles;
J repairs additional residuals but leaves173; six-vertex seed inspection
finds a count15 source; deleting its degree-one vertex reveals it was
already edge-generated; rerooting exposes response59; varying the O block
reveals the ternary coefficient; core-slot counting proves it; replacing
marked U by the graph minimum mu restores a rigorous inverse potential.
No step promotes the finite prefix to universal odd coverage.
