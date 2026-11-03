# Recovering the recursion modes without dropping their carriers

**INHERITED PROVED mechanisms; elementary restatements below.**
This is a synthesis of current sources, not a claim that an analogy closes
the unrestricted H spectrum or LRC(14).

## 1. Seven terms: additive coordinates, multiplicative assignments

For finite coordinate sets A,B,C, inclusion-exclusion says

    |A union B union C|=|A|+|B|+|C|
       -|A intersect B|-|A intersect C|-|B intersect C|
       +|A intersect B intersect C|.                         (10)

If every coordinate is an independent bit, exponentiating gives an exact
multiplicative quotient for the assignment counts2^|set|. This is the
literal A+B+C-D-E-F+G law and its addition/multiplication duality.

In the fixed-Hamilton-spine model, the free coordinate count is
q_n=binom(n-1,2), so

    q_n=3q_(n-1)-3q_(n-2)+q_(n-3), n>=4;
    2^q_n=(2^q_(n-1))^3 (2^q_(n-3))/(2^q_(n-2))^3.

[THM-553, tile addresses](../../../01-canon/theorems/THM-553-tile-address-recurrence.md)
gives the precise birth/mirror clocks and the three-region cover.
These are coordinate/presentation counts, not isomorphism-class counts or
Hamilton-path counts. The half-tiling coordinates parametrize reflection
fixed points; [THM-549 with its THM-830 correction](../../../01-canon/theorems/THM-549-halftiling-complement-quotient.md)
does not identify them with all quotient classes.

The next operation matters. For a group element preserving a cover, the
same inclusion-exclusion can count its coordinate **orbits**, then
exponentiate its fixed assignments. Averaging these fixed counts by
Burnside does not commute with multiplying their separate averages.
The joint permutation and its intersection action are the missing
sidecar. This explains why an exponent recurrence need not descend to
the quotient census.

## 2. +2, times2, and affine port growth

Distinguish adding two vertices, adding two to H, doubling a class count,
doubling an orientation count, and doubling the charge(H-1)/2.
Every tournament has odd H, so literal H -> 2H is impossible. Charge
doubling would instead send H to2H-1.

The recent birth grammar adds a new vertex beating exactly one named
old port. Its child values depend on the source and the chosen port.
[Tree constructions](../primegame_port_spectrum_20261103/TREE_CONSTRUCTIONS.md)
give two H17 sources with respective child lists

    33,27,41,59,85,47
    49,29,31,35,77,71.

Thus one supplies2H-1 and the other2H+1; scalar H17 cannot choose a
common source with both branches. The stronger
[induced-source obstruction](../metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md)
excludes31->63 under *any* old-edge-preserving extension.
The earlier pattern of “all odd targets except7,21” stays OPEN;
a universal source-preserving affine induction is already disproved.

Equation(9) in FOUR_BITS.md relates +2 in **order** to doubling the
**class count modulo a stated power of two**. It supplies an exact bridge
at a different observable, not a way around the Boolean source obstruction.

## 3. The complete structural recurrence with an explicit unknown input

The incoming [joint spectrum packet](../tournament_joint_parity_20261102/spectrum.md)
recovers the all-order SCC and square-central laws. Here is an equivalent
coefficient form. Let A(n,h), S(n,h), C(n,h), U(n,h) count ordinary,
strong, self-converse, and strong self-converse classes with H=h.
Set A(0,1)=1 and all other A(0,h)=0. Then

    A(n,h)=sum_(k=1..n) sum_(d|h) S(k,d) A(n-k,h/d).          (11)

Every tournament is a unique ordered word of its strong components.
All cross arcs point forward; a Hamilton path traverses components in
that order and independently chooses one Hamilton path in each. This
proves both the class bijection and H multiplicativity.

Converse reverses the word and converses each component. A self-converse
word is uniquely L, optional central strong Q, L^op. Therefore

    C(n,h)=sum_(2m<=n) sum_(d^2|h) A(m,d) U*(n-2m,h/d^2),   (12)

where U*(0,1)=1 and U*=U at positive order. This is an ordinary class
recurrence with multiplicative H labels, not an exponential generating
function or independent random halves. Formula(12) exposes the exact
doubling/squaring relation: order of the outside half doubles, its H
squares. It has no need to make two independent choices for the halves.

Our full class census rechecks(11),(12) through6, including actual
converse classes. These recurrences determine every aggregate once the
strong inventories are supplied. They do not supply the missing
all-order strong inventories themselves. The finite membership method
and the prime/7-prime lane reduction are credited to the incoming packet,
not claimed solved.

## 4. Toothpicks and compatible axes

The incoming [toothpick audit](../tournament_joint_parity_20261102/toothpick.md)
recovers the Applegate--Pol--Sloane rules with their actual exposed-port
predicate. Its primary source is [The Toothpick Sequence and Other
Sequences from Cellular Automata](https://arxiv.org/abs/1004.3036).
At powers of two its cumulative count obeys U_(k+1)=4U_k-1, with
U_k=T(2^k)=(2*4^k+1)/3. The off-dyadic recurrence has its indispensable
-1 correction and the two clocks i,i+1.

[Session149's correction](../session149_tournament_spectrum_20261103/toothpick/CORRECTION.md)
proves that half-turn fixes just the seed, whereas converse fixes at
least two tournament classes at every order>=3. Thus there is no
equivariant injection from the full tournament class set to toothpicks
with that half-turn. Our compatible-action theorem provides the
general test: fixed-locus sizes constrain compatible symmetries.
It does not make a new toothpick/tournament bijection.

## 5. Even graphs and the limit of numerical transfer

The theorem of [Royle, Praeger, Glasby, Freedman and Devillers](https://arxiv.org/abs/2204.01947)
equates the number of tournament classes with the number of “even”
graph classes. “Even” here concerns parity of orientation reversals by
automorphisms, not Eulerian degree or edge-count parity.
Thus every arithmetic congruence for a_n transfers to that class count.
An H-preserving map, counterpart of converse, or port operation on those
graphs does not follow from the counting theorem.

Connection contract throughout: source and target are declared above;
the maps preserve counts only at their specified level. Quotienting
loses actual incidences, individual orientations and often symmetry
fixed loci. The cheapest decisive tests are order3 fixed points,
the order5 non-tiling blocks, and the actual12-vertex81/27 source.
They reject stronger identifications while leaving the exact recurrences intact.

