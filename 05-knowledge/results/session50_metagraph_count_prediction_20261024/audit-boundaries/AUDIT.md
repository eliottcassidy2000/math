# Session 50 independent audit: information boundaries

Status: **PASS** for the stated proof claims and independently checked finite
scope. No native next-order prediction theorem is claimed or obtained here.

The audited producer texts are exact snapshots in
`sources/identifiability/THEOREM.md.txt` and `sources/inertia/THEOREM.md.txt`.
`bytepinned.txt` and `SOURCE_PINS.json` identify every snapshot by SHA-256 and
byte length. Producer source code is archived inertly and is not imported or
executed by the independent checker. No producer files or Git state were edited.

## 1. Relaxed complement/reflection realization

**PROVED, proof review passed.** The realization iff is correct for n >= 3 and
the declared partition category: nonempty odd ordinary fibres, with R acting
on those fibres, but no requirement that C act on the ordinary-fibre set.
This last distinction matters: C supplies graph edges, not an additional
ordinary-class permutation that the construction must respect.

The actual cube has B/2 two-word Klein orbits and (N-B)/4 four-word orbits.
Both C and CR are fixed-point-free, whereas R is not. Thus the *whole* Klein
action must not be described as free. Assign a two-word orbit to each blue
edge and a four-word orbit to each K edge. A positive-blue vertex becomes one
R-stable fibre of size b + 2 deg_K, and a zero-blue vertex receives two
R-exchanged sheets of size deg_K each. Since the latter is positive odd, no
empty or even ordinary fibre is produced. The chosen sheet assignments are
unconstrained by C and can be made independently for each R pair.

Loops have the required incidence factors:

- A blue loop adds two R-fixed words and degree two at its vertex.
- A K loop contributes two literal black loops, degree four, and four words.
- At an NS endpoint of a K loop, two words belong to each ordinary sheet.

**FINITE-EXACT independent exhaustive control.** At n=4 the checker enumerates
all 4,140 set partitions of the eight actual cube words, retains all 61 odd
R-stable ordinary partitions, and forms their colored quotient graphs. It
separately generates every possible assignment of the six endpoint slots of
two blue edges and one K edge, filters by the theorem's parity conditions,
and compares graph isomorphism types. Both routes give exactly 11 types.
Forty eligible endpoint assignments have loops. This is a genuine exhaustive
small-instance iff check, not just replay of the producer's examples.

The order restriction is essential: at n=1 or n=2 there is no free coordinate
and C is not fixed-point-free, so the same orbit decomposition does not apply.

## 2. Hostile controls and exact ambiguity scope

**FINITE-EXACT independent checks passed.** The auditor classifies every native
fixed-path word through n=6 using all actual Hamiltonian permutations of an
unclassified representative. It does not use the producer's score-block
canonicalizer, a precomputed class assignment, or its graph-generation code.
It recovers (A,V,S) = (2,2,2), (4,3,2), (12,10,8), (56,34,12) at n=3..6,
and reproduces the full native n=6 joint degree/blue-degree histogram and
blue/black loop totals.

Every one of the 1,024 words in each of the four saved hostile partitions is
checked against actual C and coordinate R. The check includes exact ordinary
sheet sets, R fix/exchange behavior, odd masses, literal degrees with loops
counted twice, edge multiplicities, and a concrete nonnative score mismatch.
Both bounded models have N=1024, B=64, V=34, M2=48724, two singleton blue
vertices, two blue loops, 24 black loops, and connected nonempty support in
both layers. Every signature belongs to the independently recovered native
palette. Their S values are 12 and 20, giving abstract A values 56 and 48.

The exact identification set {48,56} has the advertised restricted scope.
Exactly two blue-degree-one vertices and the remaining allowed positive
blue degrees 3..9 imply 3S-4 <= 64 <= 9S-16. Together with even S this gives
10 <= S <= 22. M2 = 4V-3S modulo 8 gives S = 4 modulo 8, leaving S=12 or 20.
The two saved witnesses realize both possibilities. Loop totals and layer
connectivity introduce no further exclusion because both witnesses have them.

The conclusion does **not** preserve the full palette as a set, its
multiplicities, the literal degree histogram, colored incidence, or native
tournament realization. In fact the score witnesses refute native realization
for both arbitrary lifts, including the lift with the native histogram.
Only the actual native partition is a tournament census. The separate
same-uncolored-support pair has its deliberately nonnative degree-1003 center;
it must not be promoted to a native mass-bounded example.

## 3. Finite moments and degree normalization

**PROVED, proof review passed.** Rational interpolation gives an integer signed
measure D delta_2 - sum_j D L_j(2) delta_(2j+1) annihilating moments 0..r.
The zeroth and first moment equations imply D is even. Its positive and
negative parts have equal vertex counts and moments, but odd-vertex counts
separated by D. A common degree-one leaf, when needed, makes their common
leaf mass N0 odd. Adding the common center of degree N-N0 preserves all
moment equalities and gives an odd center degree.

The previously easy-to-miss integrality point is sound: N0-S_leaf is even
because each odd leaf contributes d-1 and the only even leaf degree is two.
For a sufficiently large n, N-B is divisible by four. Therefore
(N-B-2N0+2S_leaf)/4 is an integer; the stated strict padding inequalities make
both this K-loop count and (B-2S_leaf)/2 blue-loop count positive. Every leaf
is connected to the center and both centers have loops, so the uncolored
support is the same star-plus-loop graph. The actual-cube realization theorem
then applies to both, at the same order n.

The independent checker recomputes interpolation and all parity/integrality
conditions for r=1..12, comparing the saved relations and chosen orders for
r=1..8. It does not expand exponentially large multisets for this check.
The proof, rather than that finite run, supplies the all-r statement.

The merged-degree congruences are correct: d^2 is 1 modulo 8 on odd SC
degrees and 4 modulo 8 on NS degrees 2 modulo 4. Thus
S = 4V-3M2 and A = 3M2-2V modulo 8. For j>=2, odd units satisfy
d^(2^j)=1 modulo 2^(j+2), while the NS contribution vanishes. These are
merged-degree moments, not power sums of ordinary odd fibre masses. The
counterexamples need orders and degree ranges growing with r, so they do not
contradict bounded-support recovery from sufficiently many moments.

## 4. Owner-labelled holonomy and the count sandwich

**PROVED, proof review passed.** An authenticated arrow is an isomorphism from
the target word's entire source to the current word's source. Tree transports
r_x consequently give actual root-source Hamiltonian orders. Each
r_x p inverse(r_y) is a root automorphism. No score-only or path-fragment
comparison is sufficient to authenticate this arrow.

For b distinct words, subgroup order g, total source automorphism order a,
and full ordinary word-fibre size h=H/a, both b<=h and g<=a hold. Hence bg<=H,
and equality holds iff both inequalities are equalities. The bg root-source
orders gamma composed with r_x are distinct. With reverse-closed edges they
are exactly the connected marked-path lift: any walk differs from the tree
path by a loop product, and every loop product followed by the tree path is
a lift. Thus the refined covering argument is valid.

Changing only the spanning tree at fixed root and gauge preserves the
holonomy subgroup; changing root/gauge conjugates it. The producer repaired
the original imprecise conjugacy wording before the frozen snapshot.

Summing bg/H over components covering X_n gives a lower bound because each
word's certified subgroup has size at most its true automorphism group.
Component count gives an upper bound because no authenticated component
crosses a native isomorphism class. The deficit is exactly
sum_x (a(T_x)-g_B(x))/H(T_x). In particular the lower bound can already be
exact on a disconnected rigid fibre: isotropy completeness is not connectivity.
The adaptive rescue strictly increases the positive integer bg bounded by H,
provided the new global path lies outside the existing lifted component.
Its termination is correct and makes no efficiency or bounded-locality claim.

The identity-loop mass, full inertia mass, cycle-type correction, and raw-loop
count distinctions are correct. Aut(T) acts freely on actual ordered paths,
so the marked observer stabilizer is one even when closed normalized-word
reobservations retain nontrivial source automorphisms. In particular the
triangle has (H,h,a)=(3,1,3). Source inertia *cardinality* must not be confused
with its number of isomorphism classes. None of these identities is itself a
new cross-order recurrence.

## 5. Independent native reobservation checks

**FINITE-EXACT, independent scope n<=6.** For each native source representative,
the auditor constructs the graph on literal Hamiltonian permutations, joining
two paths by one valid single-vertex relocation. Every such path graph is
connected through n=6. Direct full-permutation automorphism checks match the
producer's subgroup order; the saved generators authenticate and close to
the complete independently computed automorphism group.

For span-three moves the auditor constructs connected components in that
marked-path graph, projects each to its normalized word set, and recovers
g from marked-component size divided by projected-word count b. It does not
reuse the producer's tree transports or fundamental-loop extraction. It
independently reproduces:

| n | lower bound | A(n) | normalized components | unsaturated components |
|---|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 0 |
| 4 | 4 | 4 | 6 | 3 |
| 5 | 176/15 | 12 | 35 | 31 |
| 6 | 80758/1485 | 56 | 420 | 414 |

The unrestricted directed-arrow totals independently agree with
2^m(n^2+n-4)/8 through n=6. The all-size proof correctly treats the two
endpoint relocations of each interval of length at least three; each requires
one interior free arc and one additional distinct free arc for each existing
outside boundary. Adjacent swaps are never valid. This edge total does not
imply connectivity or determine the class count.

The producer's full n=7 result is retained as its **FINITE-EXACT** result,
not falsely advertised as a separate auditor rerun. No all-size relocation
connectivity theorem has been established by either finite run.

## 6. All-size window hostiles

**PROVED, proof review passed.** In the backward-nonpath source T*_n any
proper rerouted interval has either a left outside neighbor forcing its
first owner, or a right outside neighbor forcing its last. Induction then
forces the full original order. Global cyclic rotations remain valid, and
for n>=4 the score classes have sizes at most two; internal arcs prevent
swapping a two-owner score class. Hence it is rigid but has at least n
distinct normalized observers. The auditor exhaustively checks all proper
window permutations and all source automorphisms through n=8.

For C3[TT_k], every window of at most k+1 owners meets at most two consecutive
chains and is transitive with a unique Hamiltonian order. Yet a global block
rotation is a genuine order-three automorphism. Therefore the concatenated
word is locally isolated with g=1<a and has missing inertia weight at least
2/H. The whole fixed-order lower bound is strictly deficient. Independent
window permutations and rotation checks pass for k=1,2,3,4.

These are contiguous-span obstructions. A single removed-and-reinserted
vertex can span the whole path, so they do not refute unrestricted-span
single-vertex relocation connectivity. Nor are reobservation edges the
original complement edges: reobservation keeps the source tournament fixed
and adds authenticated owner permutations that an unmarked complement
metagraph does not contain for free.

## 7. Replay and final boundary

`check.py` is standard-library-only and uses exception-based checks active
under `python -O`. Normal and optimized stdout are byte-identical.
`verify_manifest.py` validates all packet hashes and source byte pins;
`verify_portable.py` copies only manifest-listed files to a new temporary
directory, reruns both modes there, and checks byte-identical output.

No blocking mathematical defect remains within the stated scopes. The
remaining predictive gap is native cross-order closure: producing enough
next-order incidence, moments, or authenticated observer information without
already reconstructing that whole next-order native universe. Relaxed
countermodels delimit one information level, not the possibility of every
native recurrence.
