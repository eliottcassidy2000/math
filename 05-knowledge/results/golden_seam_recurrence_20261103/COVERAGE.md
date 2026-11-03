# Realizable value search and the prime lane obstruction

**PROVED search reduction; FINITE-EXACT closed-state certificates.** All
source words use actual labelled strong atoms and the same path-preserving
complement. No choice of arbitrary positive coordinates is accepted as a
tournament source.

## V1 Why a finite search covers every size

Every certified matrix P has nonnegative entries, P_11>=1, and nonnegative
other entries in its first column. Therefore h(vP)>=h(v) for every positive
state. The cone

    C_K={ (h,x,y): h>=0, 0<=x<=Kh, 0<=y<=Kh }

is invariant for K=4 in the golden alphabet, K=8 after adjoining the strong
five-vertex charts, and K=16 for the six-vertex bank. Its four extreme rays
are(1,0,0),(1,K,0),(1,0,K),(1,K,K). Applying every matrix and checking the
linear inequalities on these rays is an exact finite certificate of cone
invariance for every state. The verifier does this independently of the
reachability search.

For a value bound B, every state on a word ending at h<=B also has h<=B,
because h never decreases. Its other integer coordinates lie in0..KB.
Thus it belongs to a finite state set. Starting at(1,0,0), saturate all
matrix transitions whose new first coordinate is at most B. Saturation
contains exactly the states of every source word with h<=B, at **every
vertex order**. The proof is induction on the word length in both
directions: each discovered transition constructs a literal source word,
and every word ending below B stays within the retained region throughout.

For a minimal source order, assign each transition the vertex size of its
source atom and run Dijkstra. The independent checker instead saturates by
a LIFO stack, ignoring all costs, parent pointers and discovery order.
Both implementations give the same sorted-state SHA-256 digest. Every
stored value witness retains its full source word and three coordinates.
It can be reconstructed into a labelled tournament without solving an
isomorphism problem.

The minimum source cost is only within the declared alphabet. It is not a
new global minimum order theorem. Empty initialization contributes h=1;
the singleton a is a nonempty witness for that value.

## V2 The golden family is insufficient

At B=100000 the golden family has119590 reachable triples and30194 odd
outputs. Its prime lane has4907 of9590 eligible primes, and its7p lane has
851 of1674 eligible targets. The first missing prime is37; the first
missing7p target is133=7*19. These remain absent at every size by V1.
`exclusions.json` includes the complete small closed state sets at B=37
and133, with every outgoing transition checked to stay in the set or
exceed the bound. No argument from a finite order cutoff is used.

The first missing odd values are

    7,21,37,39,43,47,63,67,71,79,83,103,113,115,117,... .

Thus neither a finite matrix recurrence nor Fibonacci provenance implies
the global two-gap conjecture. For example(37,1,1) satisfies positivity,
oddness, primitivity and C_4, but is not reachable. The missing condition is
membership in the actual matrix-generated orbit of(1,0,0).

This does not refute the unrestricted tournament spectrum claim. The
restriction is the source SCC alphabet, followed by a single actual
complement. The target can contain more complicated strong substructures,
but it must still have such a source certificate.

## V3 The first repair needs a five-vertex source atom

Let r13 be the strong five-vertex chart with mask13 on the six nonpath
pairs in lexicographic order. Its out-neighbor rows are

    (22,12,8,17,6).

Its own Hamilton count is11; after actual complement it has count13.
The compiled matrix is

    P_r13 = [[13,6,9], [28,18,24], [24,16,22]].

The exact words

    r13 a -> (37,22,31),       q0 r13 -> (133,178,193)

give strong targets at orders6 and9. Both are independently checked by
subset Hamilton DP on the full labelled relations. They contain one
five-vertex source atom. Since every source with strong atoms of size at
most4 belongs to the golden family, size5 is the smallest source-atom
extension that can repair either obstruction under this operation.

This is an operation-scoped minimality claim. The five-vertex chart is not
a new tournament isomorphism class or a previously unknown realization of
37 or133. It identifies the missing legal source coordinate in this family.

At B=100000, the eight-symbol family attains35661 odd outputs,6000 eligible
primes and1042 eligible7p targets. It still misses47, among many others.

## V4 Complete five-vertex bank and eight six-vertex patches

Adjoin all50 strong five-vertex path charts, keeping their full labelled
relations even when their43 counting matrices agree. This57-symbol alphabet
attains4990 of the5000 positive odd values at most10000. Its complete
missing list is

    7,21,47,83,103,491,709,833,3227,4619.              (V1)

All these absences hold at every order within this alphabet. In particular,
source atoms of size at most5 cannot realize47 by this construction.

Every strong six-vertex chart also satisfies the same integer intertwining
identities. There are903 such charts and544 distinct matrices. The following
eight actual six-vertex charts suffice to patch(V1):

| Value | Source word | Source order |
|---:|---|---:|
|47|a s14|7|
|83|a s268|7|
|103|a s330|7|
|491|a a s248 a|9|
|709|a s261 a a a|10|
|833|a s85 b|10|
|3227|a r35 s262|12|
|4619|a s453 r29|12|

Here sj means mask j on the ten nonpath pairs of a six-vertex atom, with the
same orientation convention as qj; rj is its five-vertex counterpart.
`atom_atlas.json` stores every row list, so the masks have no unstated gauge.
The listed words are witnesses, not claims of globally minimum order.

The resulting fixed65-symbol alphabet attains every positive odd number
through10000 except7 and21. The57-symbol bank already handles4990 values,
so the eight explicit witnesses give a separate sufficient coverage
argument. No further charts were selected after this range.

The eight charts were selected **after** inspecting the holes of the
five-vertex bank; they are not an independently prescribed universal set.
No minimality in their number is asserted. The holdout test, with exactly
the same65 charts, attains14998 of15000 positive odd values through30000;
the only holes are7 and21. Its complete saturation has644767 states and
is independently replayed by Dijkstra and LIFO reachability, with the same
sorted-state digest. Every positive hit retains a source word. This
extension beyond the selection range is evidence for the fixed alphabet,
not a proof of global completeness.

For every reported prime or7p hit, the retained witness is also checked to
be strong. Whenever a word contains at least two components and at least
three vertices, its source is nonstrong; complement reverses the endpoint
arc and closes its retained Hamilton path into a cycle, proving target
strongness. The two-singleton word is excluded from this argument: its
endpoint arc is the path arc and is retained. Single-atom boundary cases
are checked directly.

## V5 The remaining proof problem

The inherited
[THM-4094, two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
reduces the unrestricted claim to strong realizations of primes p!=7 and
7p with prime p!=3, with its stated carry witnesses. The old full spectrum
banks already extend much farther numerically than10000; no new overall
coverage record is claimed. The gain here is a **fixed finite alphabet,
three exact integer coordinates and explicit realization guards**.

A possible route to a proof is a guarded inverse for these matrices. For
example the singleton rule sends

    (h,x,y) -> (h+2x, y, h+x).

Its inverse at(H,X,Y) is

    (2Y-H, H-Y, X),

requiring H/2<Y<=H and nonnegative coordinates, together with actual source
reachability. This inverse decreases the vertex cost on a known word; the
cone inequalities alone do not certify that a desired prime H has a legal
preimage. The golden37 obstruction demonstrates the missing implication.

Useful precise next targets are: a complete semigroup membership invariant
for the65-symbol rule; an inverse covering theorem for all prime target
fibres; and a proof or counterexample to the same three-coordinate
intertwining for strong atoms at arbitrary orders. None is asserted here.
