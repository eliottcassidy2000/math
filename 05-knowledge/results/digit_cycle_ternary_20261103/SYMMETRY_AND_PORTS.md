# The digit clock acts on the hostile relation; orientations retain a subgroup

**PROVED for the explicitly defined infinite family; FINITE-EXACT literal
group checks at twelve vertices and affine checks at larger levels.**
Local labels G1--G6. Standard modular substitution and Burnside mechanisms
are inherited; no claim of classifying all tournaments is made.

## G1 An infinite family containing the 81/27 pair

Let A_1 be a directed triangle with a new source dominating it. Let B_1
be its converse. Define

`A_k=C3[A_(k-1),A_(k-1),A_(k-1)]`,
`B_k=C3[B_(k-1),A_(k-1),A_(k-1)]`, k>=2.

The three blocks follow the indicated cyclic orientation. Keep all vertex
names. Equivalently, B_k reverses every internal edge of one bottom
four-vertex block of A_k. They have `n_k=4*3^(k-1)` vertices and exactly
the same **labelled** cyclic triples. A triple in one bottom block is
cyclic exactly when it is the three inner triangle vertices. A triple
using two vertices in one module and a third outside is never cyclic.
All other triples depend only on the unchanged quotient. This proves
equality of the whole labelled relation, not only its size.

Write a_k=|Aut(A_k)| and b_k=|Aut(B_k)|. Then

```text
a_1=b_1=3,
a_k=3 a_(k-1)^3,       b_k=b_(k-1) a_(k-1)^2,
a_k=3^((3^k-1)/2),    b_k=3^((3^k-1)/2-(k-1)).
```

Proof: the three top blocks are exactly the maximal proper modules.
A module meeting two of the blocks must contain all of the third,
because an omitted vertex there sees the two included blocks in opposite
directions. Applying the same argument forces all three blocks to be
entirely included. Thus every proper module lies in one top block.
Automorphisms consequently induce only cyclic permutations of the three
blocks. For A_k all three rotations are possible. In B_k the exceptional
block is not isomorphic to the other two, so only the identity rotation
is possible. A_j and B_j are nonisomorphic: A_j has degrees
`2*3^(j-1)+1` at its added sources and `2*3^(j-1)-1` at all inner vertices;
one flipped block changes its source degree by -3 and its three inner
degrees by +1. This proves the recurrences and their closed forms.

The common number of cyclic triples also has a closed recurrence:

`c_1=1`, `c_k=3 c_(k-1)+(4*3^(k-2))^3`.

With B=3^(k-1), it is `c_k=(8B^3-5B)/3`.

|k|Vertices|Cyclic triples|Aut(A_k)|Aut(B_k)|Ratio|
|---|---:|---:|---:|---:|---:|
|1|4|1|3|3|1|
|2|12|67|81|27|3|
|3|36|1929|1594323|177147|9|

The familiar 81/27 control is exactly k=2, up to relocating the flipped
bottom block by the outer cyclic symmetry. The higher rows follow from
the proof, not from a full automorphism census at those orders.

## G2 Translation and doubling on the actual labelled relation

Let N=3^k and B=N/3. Label the inner vertices by x in Z/N and the added
sources by s_b, b in Z/B. The bottom block of x is b=x mod B. Within its
three inner vertices use the cyclic order of the quotient digit x/B.
Between different bottom blocks, orient b toward c when the first nonzero
ternary digit of c-b is 1. The source s_b dominates its own three inner
vertices. This is A_k, with a fixed orientation gauge and no ties.

Every affine map `x -> a x+d` induces the simultaneous source map
`s_b -> s_(ab+d)` modulo B. It preserves all cyclic triples. If a=1 mod 3,
it preserves the quotient and the inner triangle orientations. If
a=-1 mod 3, it reverses both kinds of directed cycles, preserving their
cyclic-triple status; the source-to-triangle spokes are kept. Therefore
the whole group of size `2*3^(2k-1)` acts on the common triple relation.

Its intersection with Aut(A_k) has size `3^(2k-1)`: exactly a=1 mod 3.
With the flipped block at zero, its intersection with Aut(B_k) additionally
requires d=0 mod B, giving size `3^k`. The radial source/sink direction
forces the flipped block to stay fixed. At k=1 the same formulas hold
with B=1.

Thus at twelve vertices the cycle/doubling group of order 54 has oriented
intersections of sizes **27 and 9** in the two hostile tournaments.
[hostile_embedding.json](hostile_embedding.json) gives literal permutations
in the old twelve-vertex labels, not just an abstract group description.
The complete 54-element group is checked against all 67 shared triples.

There are two different groups of order 27 here. The affine odd subgroup
has element-order counts `{1:1,3:8,9:18}` and is nonabelian. The second
tournament's actual automorphism group is C3^3, with counts `{1:1,3:26}`.
The digit clock's order-nine translation lies in the first tournament's
automorphism group, but not in the second's. This gives a direct reason
why the shared triple representation cannot authorize the same oriented
source operation on both.

## G3 All-order shared-port recurrences

Adjoin p **named new ports**, fixing their mutual orientations. Each old
vertex has q=2^p possible incidence bit-vectors to the ports. Automorphisms
fix the ports pointwise, so their restriction to the old source is exactly
the stabilizer of this q-colouring. These formulas use an unmarked old
source; fixing additional old roots requires the corresponding subgroup.
No strong-output or Hamilton-count filter is imposed.

Let O_k(q), E_k(q) count these attachment orbits for A_k, B_k respectively.
Let R_k(q), F_k(q) count labelled attachments with trivial port-fixing
automorphism group. Then

```text
O_1=E_1=(q^4+2q^2)/3,
R_1=F_1=q^4-q^2,

O_k=(O_(k-1)^3+2 O_(k-1))/3,
E_k=E_(k-1) O_(k-1)^2,

R_k=R_(k-1)^3-a_(k-1)^2 R_(k-1),
F_k=F_(k-1) R_(k-1)^2.
```

For the orbit counts, apply the three-cycle to the three child colour
orbits; the one-flip case permits no top rotation. For regular attachments,
each child must first be regular. A top rotation then stabilizes the
colouring exactly when all three child colourings belong to the same
orbit. There are R/a such orbits, each with a elements, so this bad set
has `(R/a)*a^3=R a^2` labelled triples. This is the subtraction above.

At k=2:

|Named ports|A_2 attachment orbits|B_2 attachment orbits|A_2 regular attachments|B_2 regular attachments|
|---|---:|---:|---:|---:|
|1|176|512|1620|1728|
|2|294976|884736|13821840|13824000|

The two-port discrepancy from independent orbit counts has a general
explanation. For any finite permutation group G and c(g) its vertex-cycle
count, Burnside gives

`O_G(q^2)-O_G(q)^2 = Var_(g uniform in G)(q^c(g))`.

Both ports use the same automorphism, which is the joint datum. For q=2
the excess is 264000 on A_2 and 622592 on B_2. This is an explicit
application of the earlier shared-port variance mechanism, not an
identification of this finite group average with a branching process.

Independent controls use Burnside on the actual automorphism lists and
inclusion-exclusion on orbit partitions of order-three elements. Every
nontrivial stabilizer in a 3-group contains such an element. The latter
calculation confirms the regular-attachment polynomial without using the
recursive R formula. See [attachments.json](attachments.json).

## G4 Which groups and prime powers can these attachments retain?

For A_k, the coloured automorphism group has an exact recursive form.
If the three child coloured types are not all isomorphic, it is the direct
product of their three groups. If they are all isomorphic to one coloured
type with group H, it is `H wr C3`. At a bottom block the group is C3
when all three inner colours agree and is trivial otherwise. Actual
labelled transporters are needed to implement an abstract wreath-product
identification. In B_k the exceptional top branch cannot rotate, so its
group is always the corresponding product, recursively.

This classifies the stabilizers by the complete coloured-source grammar.
Only the prime 3 occurs. More sharply, **one named port already realizes
every order**

`3^j, 0<=j<=(3^k-1)/2` on A_k,
`3^j, 0<=j<=(3^k-1)/2-(k-1)` on B_k.

Proof by induction. Each attainable stabilizer order has at least two
distinct binary-coloured types: toggle the colour only on all added
source vertices. Every automorphism preserves that distinguished subset,
so the stabilizer stays unchanged. Its size B is odd, so the number of
ones on it changes and the new colouring is not isomorphic to the old.
At level k, any exponent from 0 to three times the previous maximum is
a sum of three attainable exponents. If all three are equal, choose two
different coloured types to suppress a top rotation. The final exponent,
three times the previous maximum plus one, is obtained by three identical
maximal types. The B recurrence is a direct product, so its exponent
interval is the sum of its three child intervals. At the bottom both
orders 1 and 3 occur. This proves the assertion. Extra named ports may
be given constant incidences and do not remove these realizations.

For B_k the distinguished added-vertex set remains recognizable from the
recursive block decomposition: ordinary source and exceptional sink
vertices are both the unique singleton components of their bottom blocks.
Its size is still odd, so the same toggle argument applies.

## G5 Scope of the classification

The recurrences classify this explicit cyclic composition family and its
named-port attachments. They neither enumerate every tournament class
nor provide a complete A000568 reconstruction. The index-two clock
operations preserve triples but not all orientations. The marked 2107
return uses actual supplied paths and full endpoint observers separately;
an affine relation clock cannot replace those guards.

## G6 A native normalization makes two structural halvings possible

Every order n_k=4*3^(k-1) has valuation two, but the original A_k and B_k
have **no two-vertex modules**. A bottom source-plus-triangle block has
none, and vertices in different blocks are distinguished by the regular,
pair-free quotient. Size divisibility alone therefore supplies no pair
contraction.

There is an explicit repair. In each of the B=3^(k-1) bottom triangles,
reverse one directed edge. Each becomes transitive, and each whole bottom
block becomes TT4, with its designated added vertex still at the source
or sink end. The result is `W_B[TT4]`, where W_B is the regular ternary
quotient on B vertices.

Every one of these repairs has a legal two-native-complement implementation.
Concatenate Hamilton paths in the blocks ordered 0,...,B-1; consecutive
quotient vertices differ by 1 and have the needed direction. In the block
being repaired, choose its internal path to contain the chosen triangle
arc with the second vertex not last. The helper construction in
[R2](../return701_algebra_20261103/RETURN701.md) then reverses exactly that
arc using two supplied-path complements and deletes its temporary helper.
Repeat for all blocks. All old names and the quotient frame survive.
The cost is `2B` native moves, B helpers, peak order `4B+1` and no old
deletions. B reversals are optimal for this prescribed bottom-block
transitivization; no global repair minimum for other interfaces is claimed.

The pair-module graph of the result is exactly B disjoint paths P4.
Inside TT4 only consecutive pairs are modules. Across blocks the regular
quotient has no pair module: a pair module would give its two endpoints
outdegrees differing by one, contrary to regularity. Its distinguishing
third vertex therefore also distinguishes any proposed cross-block pair.
Each even path has a unique perfect matching, forced from an endpoint.
Contract that matching to obtain `W_B[TT2]`, then contract the remaining
forced pairs to obtain W_B. This is precisely the prior intrinsic
halving mechanism, now reached by an explicit native normalization.

The receipts retain the original-name fibres after both contractions.
The halving clock is exactly two, and the final odd-core order is B.
This is halving of **vertex order in a framed encoding**, not halving
of Hamilton-path counts. In particular it does not replace the separate
2107 return or supply a Collatz convergence argument.
