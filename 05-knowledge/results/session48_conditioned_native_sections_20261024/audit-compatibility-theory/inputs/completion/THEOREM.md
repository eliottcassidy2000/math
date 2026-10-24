# Exact global-observer completion from coherent local records

Status: **PROVED all-size labelled completion theorem and generating function**;
**FINITE-EXACT** small exhaustive controls. This is a conditional realization
result, not all-size native metagraph wiring or harmonic-bank transport.

## 1. Inheritance and the change of question

Incoming `compatible_charged_helpers_20261023/controller/PARTNER_GLUING.md`
classifies when local K-partners glue after source coherence. Session47's
`global-observer/THEOREM.md` retains required observer arcs, contracts their
path components, and tests the remaining transition graph in a FIXED full
source host. The current source of the externally pinned controller,
`external_pin_sections_20261024/observer.py`, has conflicting required
successors4->6 and4->11. Adding private cells cannot remove that conflict.

Here we distinguish another legitimate problem: the local source cells and
observers are fixed, but uncovered source arcs have NOT yet been chosen.
Exactly when can we complete them to a global source/observer/partner record,
and how many labelled completions are there? The answer does not need a
search through all host tournaments.

## 2. Full data and coherence

Let V be a nonempty finite NAMED owner set. For each cell C_i⊆V, supply its
entire source tournament T_i and an actual Hamiltonian order P_i. Assume all
sources agree on every shared pair. Thus the local source arcs specify one
orientation on each covered pair; all other pairs are free. Let e be the
number of uncovered pairs.

Each local record prescribes its K-partner by keeping its consecutive path
arcs and reversing its other arcs. We seek a full source S on V and a global
Hamiltonian order P such that S[C_i]=T_i and K_P(S)[C_i]=K_(P_i)(T_i) for all i.
The second requirement is stronger than merely recovering P_i as a subsequence
of P: every local consecutive pair must also be GLOBAL consecutive.

All owner names and local paths are retained. This theorem counts labelled
marked records, not isomorphism classes or independent normalized capacities.
If these cells also constitute a native packet, its resource disjointness,
kernel-lift restrictions and electrical guards must be checked separately.

## 3. The necessary-and-sufficient forest gate

Let F be the union of all directed consecutive arcs of the P_i, on all V,
including any isolated owners.

**Completion theorem.** Such a global record exists if and only if F is a
directed linear forest: every owner has indegree and outdegree at most one,
and F has no directed cycle.

Necessity follows because every required local consecutive arc must belong
to the global path P. For sufficiency, let B_1,...,B_p be the directed path
components of F (singletons included). Every cell lies wholly in one block,
and its local path is a contiguous segment of that block. Indeed its required
arcs already form a directed path, and no extra arc can branch or insert an
owner inside one of those required edges. Thus every local negative adjacency
flag is automatically retained too.

Every pair joining two different blocks is uncovered: if a cell contained
both owners, its own consecutive path would connect them in F. Choose any
ordering of the blocks, orient the p-1 joining endpoint pairs forward along
that order, and orient every other uncovered pair arbitrarily. The expanded
block order is a Hamiltonian path in the completed source. Its global K
restricts to every local K because the consecutive-pair flags agree exactly.
This constructs the required global record. QED.

This iff is not the claim that arbitrary local records can be realized. A
coherent directed star can force two successors at an owner; a coherent
cycle of two-owner observers can force a directed cycle. Both fail the gate.
It is also not an iff for an already FIXED host: its cross-block endpoint
arcs may block every admissible block ordering. Session47's transition-graph
criterion remains necessary there.

## 4. Exact completion count

If the forest gate holds, let p be the number of blocks. Then e>=p-1, and

    number of labelled marked pairs (S,P) = p! * 2^(e-p+1).       (1)

For a chosen block order its p-1 connecting pairs are distinct uncovered
pairs whose source orientations are forced. The other e-p+1 pairs may be
chosen independently. Different block orders give different P, even if
they sometimes give the same source S. Hence summing counts MARKED pairs,
not necessarily distinct source hosts. Since the global partner is determined
by S,P, the same formula counts the triples (S,P,K_P(S)).

For example two disjoint two-owner path cells on four owners have p=2,e=4:
there are16 marked completions. It would be incorrect to infer16 distinct
source hosts without forgetting and quotienting P separately.

Absolute GLOBAL clock pins replace p! by the number of block orders satisfying
the specified cumulative block-length equations. An owner at offset j within
its block, pinned to global clock c, forces that block to start at c-j. The
factor2^(e-p+1) is unchanged for each allowed order. Local cell clocks are not
silently identified with these global clocks. If no order satisfies the global
pins, there is no completion even when the forest gate passes.

## 5. Observer completion polynomial and the Eulerian specialization

Fix a reference total ordering of the named owners. Weight a completion by
z raised to the number of forward-oriented UNCOVERED pairs of its source S.
For an allowed block order O, let a(O) count its p-1 connecting endpoint arcs
that point forward in this reference order. Then the exact polynomial is

    G(z) = (1+z)^(e-p+1) * sum_(allowed O) z^a(O).                (2)

The (1+z) factors represent genuinely unprescribed binary arc orientations,
not extra resources or fractional occupancy. Formula(1) is G(1). If all
covered source arcs are also to be counted, multiply by their fixed forward-
arc monomial. This distinction prevents a degree/count mismatch.

For the empty cell family on n>=1 owners, every block is a singleton. Let
A_n(z)=sum_(permutations π) z^(number of ascents of π), the Eulerian ascent
polynomial. Equation(2) becomes

    sum_(labelled S) H(S) z^(number of forward arcs of S)
       = (1+z)^(choose(n,2)-n+1) A_n(z).                         (3)

The left side counts every actual directed Hamiltonian path once. This is a
bijection: specify an observer order, keep its directed path arcs, and freely
choose every nonpath arc. In particular the total is

    sum_(labelled S) H(S) = n! * 2^(choose(n,2)-n+1).

The Eulerian name concerns ASCENTS of permutations, not Eulerian trails.
Neither formula is an enumeration of tournament isomorphism classes such as
A000568. The projection forgetting the marked observer is many-to-one.

For n=4, A_4(z)=1+11z+11z²+z³. Thus the two middle11s have a precise observer-
permutation meaning in this formula. This does NOT identify their sum22 with
the unrelated order-eight clock bottleneck: that would require another map.
The actual source-to-target correspondence here is completely explicit.

## 6. Converse and preserved information

Converse reverses every local source, every local observer order and each
completed global observer. Reflect any global clock c to |V|-1-c. The forest
condition and completion count are invariant. The uncovered orientations all
reverse, giving the exact reciprocity

    G_converse(z) = z^e G(1/z).

This is a comparison with the CONVERSE input data. It is not an assertion
that one conditioned polynomial is palindromic when its fixed records are
not converse invariant.

Connection contract:
- Source: coherent actual named source cells with fixed local path records.
- Target: directed path blocks, their allowed orderings, and free arc bits.
- Preserved: existence and exact count of global marked completions and
  all local source/partner restrictions, with declared global clocks.
- Lost if P is forgotten: the observer multiplicity, so no unmarked count
  or native normalized-port capacity can be inferred.
- Separate obligations: original-resource capacity, kernel records, current
  equations, and any requested harmonic transport to the larger host.

## 7. Exact controls

`check.py` independently enumerates coherent one- and two-cell families of
sizes2/3 on at most four owners, every orientation of their covered pairs,
every legal local observer choice, and every completion of uncovered arcs.
It compares full global Hamiltonian orders with the forest/count/polynomial
formulas, with additional star, cycle, clock-pin and converse controls.
The empty-cell identity is checked against every labelled tournament through
order5, independently counting its Hamiltonian orders. Normal and optimized
runs agree; no floating-point arithmetic or hypothetical resource is used.

Run `python3 check.py --output normal.json` and repeat with `python3 -O`.
