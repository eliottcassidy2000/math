# Three concrete maps, and why they need not agree

## A. The Sylow-3 map on three vertices

There are two tournament types: the transitive tournament and the directed triangle. There are two sign-even graph types: the empty graph and the path P3. In those row/column orders the signed Walsh matrix is

    [ 1   1 ]
    [ 1  -3 ].

Modulo three its only full nonzero matching is

    transitive tournament -> P3;
    directed triangle    -> empty graph.

The tournament automorphism orders are 1 and 3; the graph orders are 2 and 6. Their 3-parts match, and the actual three-cycle action survives on the empty graph. The previous SCC/Lyndon rank codec made the opposite pairing. Both are bijections; they preserve different information. “The bijection” is therefore incomplete terminology unless the desired invariant is specified.

## B. The circulant map on seven vertices

Use vertices F_7 and primitive root g=3. A connection set S specifies arcs x→x+s. For the Paley selector S={1,2,4}, its exponent word along 1,3,2,6,4,5 is

    f = 101010.

The binary derivative D=I+shift gives 111111; restricting to the first three coordinates gives source necklace 111. One more derivative is 000. The graph output is the empty graph.

For S={1,2,3}, the exponent word is 111000. Its first derivative has first-half necklace 001, and the next derivative has necklace 011. The symmetric graph connection set is {2,3,4,5}, the complement of the seven-cycle with differences ±1. These are the two prime-circulant types and their two target types.

The phase of the actual derivative word is retained when returning a labelled connection set. Least rotation is used only to decide a cut or return a type-valued inverse. This is what makes the forward map affine-equivariant. It is not a labelled bijection: complementary tournament selectors are already identified by a multiplier and can have the same labelled output.

## C. The primitive map on four vertices

Let a denote the singleton tournament, b the directed triangle, and c the strong four-vertex tournament. The Eulerian primitive lift of b is

    e(b)=b-a^3.

With the declared canonical orientations,

    F(P3) = -e(b)/4;
    F(K1,3) = [a,e(b)]/8;
    F(C4) = e(c)/8.

Thus the connected star witnesses the Lie bracket of a singleton with a primitive triangle. The connected four-cycle witnesses the corrected strong tournament. These are equalities of rational linear combinations, not claims that a graph literally is a tournament or a Lie bracket as an unmarked set. Nonzero coefficients then certify the primitive-support matching used by the type codec.

## A useful choice rule for future work

Choose the map according to the operation it needs to preserve:

- subgroup action at one odd prime: use the modular Walsh matching;
- cyclic multiplier data and an atlas-free inverse: use necklace differentiation and routing;
- SCC word factorization versus graph components: use primitive-support atom matching;
- legal neighborhood extensions: retain the sign-aware stabilizer automaton, including odd parent states.

No one of these maps preserves all four structures. The five-vertex simultaneous-prime obstruction and ten-vertex mixed-prime count obstruction show why demanding all symmetry data is impossible.
