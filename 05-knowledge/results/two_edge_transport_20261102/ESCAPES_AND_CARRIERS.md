# Exact escapes and the cost of retaining every named operation

**PROVED:** finite-group selector lemma and consequences stated below.
**FINITE-EXACT:** complete six-vertex escape classification and fixed-path
carrier minima through seven. Same-agent audit; no global priority claim.

## The complete target17 trap classification

Encode the fifteen edges of a six-vertex tournament in lexicographic
pair order, bit one meaning the smaller vertex dominates the larger.
Let H count directed Hamiltonian vertex words and let E=(H-17)^2.

The inherited768 isolated vertices of the E<=4 single-edge graph split
into precisely four ordinary isomorphism classes. Every one has H15 and
is obtained by adjoining a source or sink to a five-vertex H15 tournament.

| Five-vertex core scores | Core automorphisms | Source versions | Sink versions | Targets at distance2 | Least distance |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1,2,2,2,3 | 3 | 240 | 240 | 3 per source | 2 |
| 2,2,2,2,2 | 5 | 144 | 144 | 0 | 3 |

Thus480 traps have exactly three two-edge escapes;288 have none.
Every successful two-edge pair consists of disjoint edges. For the
order3 class the two intermediate counts are9 and27, depending on
which edge is flipped first. The best two-step maximum error is64.

All six-vertex H17 targets have sorted scores(1,1,2,3,4,4). Half the
minimum l1 score distance from each trap is2 or3 respectively. One edge
flip transfers only one score unit, so this explains the sharp lower
bounds independently of a search radius. The incoming concurrent
[prime-inverse move packet](../prime_inverse_ports_20261102/MOVES.md)
records the same complete escape classification and additionally gives
common-path-preserving packets and native guarded-complement distances.
The distance/barrier conclusions are independently reproduced here;
they are not claimed as discoveries exclusive to this session.

All768 traps nevertheless have a three-step route with maximum error36.
No route of any length can keep the maximum below36. Exact flood fill
from all H17 vertices at cap16 reaches none of the traps; odd path counts
make16 the largest possible squared error below36. At cap36, breadth-first
search reaches every trap in exactly three steps. Examples are:

    source72: 72 -> 200 -> 456 -> 968
    counts:   15 -> 13  -> 23  -> 17

    source140: 140 -> 12 -> 28 -> 30
    counts:    15 -> 11 -> 23 -> 17

Every arrow flips one actual edge. The first trap has a shorter route,
72 -> 104 -> 16488, with counts15 -> 9 -> 17, but its maximum error64
is larger. Shortest length and least barrier are different objectives.

The complete15-bit cube has distance distribution to H17:

| Distance | Number of labelled orientations |
| ---: | ---: |
| 0 | 1440 |
| 1 | 13520 |
| 2 | 16080 |
| 3 | 1728 |

This gives a finite Boolean inverse: decode17.py accepts any six-vertex
orientation and returns a shortest sequence of at most three edge flips
ending at H17. With --max-error36, it returns a sharp three-step escape
from every trap. The radius statement is an exhaustive finite theorem,
not an all-order inverse or an efficiency result for arbitrary targets.

### Independent certificates

scout.py enumerates every distance-two shell of all768 traps.
escape_depth.py computes full-cube distances, minimax barriers, and
independent constrained breadth-first searches. verify.py independently
enumerates permutations along every representative trajectory and target
shell. It reconstructs the four complete S6 orbits and verifies their
union is exactly the trap set. The finite verification is an exhaustive
proof in this universe; no extrapolation in n is used.

## Why a deterministic inverse needs a symmetry-breaking port

**Selector lemma, PROVED.** Let a finite group G act on source and target
sets X,Y, and fix x in X. A G-equivariant map from the orbit Gx into Y
exists with image of x equal to y iff
\[
 \operatorname{Stab}_G(x)\subseteq\operatorname{Stab}_G(y).
\]
Necessity follows from f(gx)=gf(x). For sufficiency define f(gx)=gy:
if gx=g'x, then g'^{-1}g fixes x and hence y, so the rule is well-defined.

For tournaments on the same named vertex set, an automorphism fixing a
Hamiltonian word fixes each vertex, so the action on Hamiltonian words is
free and |Aut(U)| divides H(U). Hence a relabelling-equivariant constructor
D with H(D(T))=h requires
\[
 |\operatorname{Aut}(T)|\mid h.                 (1)
\]
This is necessary, not generally sufficient. With added vertices one
must specify the extended permutation action before invoking this test.

For every present trap, |Aut(T)| is3 or5, neither dividing17. Thus **no
deterministic relabelling-equivariant inverse exists on any of these
unmarked source orbits**, regardless of the allowed edit distance.
The labelled decoder deliberately uses external vertex order to choose
a branch. Marking the already distinguished universal source/sink does
not help, since every automorphism fixes it.

A full Hamiltonian-path marking already fixes every vertex and is itself
symmetry-breaking input. This obstruction concerns an unmarked source
orbit, not a constructor that has retained such a marking all along.

One additional vertex in a moved core orbit destroys the entire
stabilizer: there are three eligible ports in the order3 case and five
in the order5 case. The selector lemma then gives an equivariant choice
on each such marked orbit. This is an existence statement for a marked
selector, not a claim that every arbitrary vertex mark is sufficient.

The free-path action is inherited; the incoming
[Session142 prime-lane analysis](../session142_tournament_odd_spectrum_20261102/minimal_gap/PROOF.md)
also proves the stronger all-order asymmetry of prime-count targets>3.
Only the elementary divisibility statement is needed here.

### Positive laws retain symmetry at a measurable price

**Support lemma, PROVED.** Put K=Aut(T). For a K-invariant probability law
supported on actual tournaments U with H(U)=h, each support orbit has size
\[
 [K:K_U]\ge |K|/\gcd(|K|,h),
 \qquad K_U=K\cap\operatorname{Aut}(U).
\]
Indeed |K_U| divides both |K| and h. Every positive support atom brings its
entire K-orbit with it. This bounds support cardinality, not merely a mean
or a formal signed sum.

For the order3 traps, the three distance-two targets form one free orbit.
For the order5 traps, the twenty distance-three targets form four free
orbits of size five. Thus invariant target-supported laws need at least
three or five atoms respectively, and these bounds are sharp. They have
zero realization energy but still require a symmetry-breaking choice to
return a deterministic orientation.

### The81/27 control gives sharp support sizes3 and1

The old twelve-vertex pair still has67 identical labelled cyclic triples
and automorphism orders81 and27. On the same named vertices, form the
ordered SCC chain

    (0), (1,2,3), (4), (5,6,7), (8), (9,10,11),

orienting each triple as its directed3-cycle. Its H is3^3=27.
Every one of the27 source automorphisms preserves this target; thus an
equivariant selector exists on that source orbit by the selector lemma.
No such deterministic selector can exist on the81 source orbit because
81 does not divide27.

The target has an orbit of exactly three objects under the81 source
group. Its uniform law attains the support lower bound81/gcd(81,27)=3.
Thus the hostile pair has **minimum equivariant target-support sizes3
and1**, respectively. Equal cyclic triples do not preserve even this
inverse-selection cost. The actual groups and target arcs are replayed.

## Closed carriers for declared command languages

Use the fixed-path chart
\[
 X_n=\{0,1\}^{\binom{n-1}{2}},
 \qquad 0\to1\to\cdots\to n-1.
\]
Let C flip every free edge. Let rho be all-arcs converse followed by
reversing the vertex order; on chart bits it permutes
(i,j) to (n-1-j,n-1-i). C and rho are everywhere defined.

At n7 let R re-mark along the vertex word
(0,5,6,3,4,1,2), provided that word is a directed Hamiltonian path.
This permutation is an involution. R maps tile24511 to31; both have H17,
whereas their C-images have H71 andH33. The guard and its path port are
part of the operation.

An exact deterministic observer must preserve the current output H, the
legality of each named command, and its destination observer class.
Iterative partition refinement produces the coarsest such congruence:

| Commands | Minimal observer states at n7 |
| --- | ---: |
| C | 3008 |
| C and guarded R | 10455 |
| C, guarded R, and rho | 16149 |
| C and every named single-vertex relocation and its inverse | 32768 |

The concurrent [selected-marking proof](../two_selector_boolean_20261102/CARRIER.md)
now gives an all-size formula for the second row's language. Its states
retain (H,H after C), a phase recording membership in D, C(D), or neither,
and one conditional response after R and C. Here D is the legality domain
of the specified nonidentity involutive marking. D and C(D) are disjoint,
which closes the transition formula. The same packet independently
obtains10455 at n7. Its second marking S gives15874 states; S is a
different operation from the converse-and-reversal rho used here, so that
number is not in conflict with16149. These all-size restricted-carrier
and second-selector conclusions are CITED, not new results of this packet.

A right relocation moves the vertex at position i to position j, j>=i+2,
retaining the relative order of all other vertices; it is legal iff the
resulting word is a Hamiltonian path. Each inverse is also a separately
named partial command. No adjacent swap is silently allowed. Adding rho
to the last language leaves the full32768-state result unchanged.

More strongly, the **legality output alone**, with no H observation,
already requires all orientations through n7:

| n | Cube size | Classes from current legality | After one refinement |
| ---: | ---: | ---: | ---: |
| 3 | 2 | 2 | 2 |
| 4 | 8 | 7 | 8 |
| 5 | 64 | 49 | 64 |
| 6 | 1024 | 563 | 1024 |
| 7 | 32768 | 10532 | 32768 |

Thus a command followed by a legality query separates any remaining
pair at these orders. Every compiled re-marking is independently checked
against literal arcs:1,004,340 transition cases. The all-order injectivity
claim remains OPEN.

**Why these minima are certified.** Initial blocks share the requested
outputs and guard vector. Each refinement splits by destination block
for every legal named command. A stable partition supports every future
word by induction. Every split has a distinguishing future word, obtained
by prefixing the command causing the split. Therefore no smaller exact
observer for that language can merge its final classes.

These are minima for fixed command names and fixed source coordinates.
They do not exclude equivariant compression that also transports command
names, or smaller carriers for restricted attachment families. A full
unmarked metagraph is still the inherited span of classes, marked tiles
and complement lines; it is not a scalar deterministic transition system.

