# Two re-marking selectors are necessary and sufficient for this quadratic code

**PROVED finite realization and an all-size parity obstruction.**
The encoding is chosen explicitly; it is not a canonical isomorphism of
the quadratic sign quotient and tournament class quotient.

## N1. A general obstruction from odd presentation fibres

Partition a finite state set X into blocks of odd size. A re-marking selector
R is any function with R(x) in the same block as x; it need not be bijective.
Suppose F:X->X has even cardinality in every point fibre, including empty ones.
Then F cannot be expressed as P2 R P1 for any fixed permutations P1,P2 of X.

For any block A, (P2 R P1)^(-1)(P2 A)=P1^(-1)(A), of odd size. But an
inverse image of any set under F is a sum of even point-fibre sizes.
Contradiction. Any fixed word of bijections before or after the one selector
is included by composing it into P1 or P2.

The hypothesis applies at every size to the fixed-Hamilton-path tournament
presentation blocks. The inherited fibre law is H(T)/|Aut(T)|
(LEM-003, `aut-acts-freely-on-ham-paths`; THM-1430,
`the-tiling-class-metagraph-dictionary-and-which-tricks-pay`). H(T) is odd
by Redei parity (THM-001, `redei`), and |Aut(T)| is odd: an even-order
permutation group would contain an involution, whose swapped vertex pair
cannot preserve its tournament arc. Hence every block has odd size.

The theorem is about deterministic selectors on the full presentation set.
Restricting the source set, using partial maps, adding auxiliary states,
or interleaving two selectors changes its hypotheses.

## N2. Actual four-vertex realization

Use the fixed path0->1->2->3. Bits b0,b1,b2 orient the three free pairs
(0,2),(0,3),(1,3) forward when1 and backward when0. The class blocks are

    strong: {0,1,2,4,5}; cycle with source:{3};
    cycle with sink:{6}; transitive:{7}.

An actual selected-path complement retains the fixed path and reverses
the other three arcs, so C(b)=b xor7. Re-marking within a class is realized
by a Hamilton path in its current tournament, followed by relabelling the
path positions0,1,2,3. The explicit legal path words are in verification.json.

Encode the quadratic numerator by this table:

|a|-7|-5|-3|-1|1|3|5|7|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|b(a)|3|0|6|1|2|7|4|5|

Two fixed selectors, listed on b=0,...,7, are

    R_before = [2,4,4,3,2,5,6,7],
    R_after  = [0,0,4,3,4,1,6,7].

They preserve every class block. Directly, at every a in the eight-state set,

    b((a^2-29)/4) = R_after(C(R_before(b(a)))).           (1)

Thus the source trajectory, its transient tail and its clock are retained
by a literal re-mark/complement/re-mark implementation. The cycle is
3->4->1->3 in this chart. Each of the four image numerators has exactly
two preimages, so N1 proves that two selectors are minimal if only fixed
bijections and class-preserving selectors are permitted. This does not
minimize the number of vertices or auxiliary memory in other encodings.

There is also a phase-placement guard. For any block A and fixed
permutation P, an exact r-cycle lying wholly in A under R_after P R_before
needs at least r states in A intersect P^(-1)(A). Those are the only
possible intermediate return states, and the final deterministic selector
must produce r different outputs. Here the strong block has only two
internal complement return states,2 and5. Placing all three cycle points
inside it therefore fails. The displayed encoding places two there and
one in the singleton source-plus-cycle block. The clock placement is
load-bearing even though the overall source has only eight points.

The construction needs named frame choices. Quotienting directly to the
four classes loses (1); their fibre sizes5,1,1,1 cannot be the quadratic
sign fibres2,2,2,2. Conversely, the four-state magnitude factor in H2 is
exact for its own future-value observer. These are two different carriers.

## N3. Literal machine realizations and their source guards

Assign the eight a values in increasing order the control primes
2,3,5,7,11,13,17,19. The ordered FRACTRAN list

    17/2, 7/3, 3/5, 2/7, 2/11, 3/13, 7/17, 17/19

realizes one F step per instruction, on inputs consisting of exactly one
control prime. Multiplying every legal state by a fixed integer coprime to
all eight control primes gives another faithful clocked copy. Products of
several control primes are outside this state encoding; priority then
chooses a rule and does not simulate parallel independent cells.

An eight-symbol tape with one active head has the equally literal TM rule
"replace a by F(a), stay". A strict radius-zero CA applies the same rule
at every site, so it commutes with translation and evolves each cell
independently. These simple finite machines are an exact representation
of the bounded graph, not universal machines or a conjugacy to PRIMEGAME.
The actual unbounded PRIMEGAME stream is used separately in ARITHMETIC.md.

This makes the information loss visible: a sign bit is erased by F after
one step, but the odd presentation-block partition prevents erasing that
bit by only one legal re-marking selector plus permutations. Two selectors
with an intervening actual complement suffice for this explicit source.
