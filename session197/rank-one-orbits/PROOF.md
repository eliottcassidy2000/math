# A rooted-orbit correspondence inside the self-converse universe

Session 197. The rank-one square construction gives an exact bridge between ordinary tournament enumeration, self-converse path-orbit mass, and a physical class that excludes every 7p value. This is a structural bijection and exact generating-function result, not an inference of spectrum coverage from moments.

## 1. Physical objects and the bijection

For a nonempty rooted tournament (S,a) on t vertices, form F(S,a) on 2t+1 vertices:

- O is a copy of S, and I is a copy of its converse.
- I -> r -> O, where r is a new vertex.
- Every crossing arc goes I -> O except the single arc o_a -> i_a.

The involution tau fixes r and exchanges corresponding vertices of O and I. It reverses every arc. The crossing graph at r has matching rank one.

Conversely every self-converse tournament with an anti-involution fixing a root r of crossing matching rank one has this form. Its out-neighborhood is intrinsically S, and its unique out-to-in crossing arc identifies a. This is the symmetric rank-one classification in the companion binary-rigidity packet.

Therefore F is a bijection between rooted tournament isomorphism classes (S,a) on t vertices and rooted self-converse rank-one classes (T,r) on 2t+1 vertices. No particular anti-involution is retained as part of either class.

The same correspondence is injective after forgetting the output root. Indeed every anti-involution of an odd-order tournament has exactly one fixed vertex. The previously proved conjugacy theorem says all anti-involutions of a self-converse tournament are conjugate by automorphisms. Given an unmarked isomorphism between two outputs, conjugate the transported anti-involution to the chosen anti-involution of the target. The resulting isomorphism sends the fixed root to the fixed root. Thus no additional identifications occur upon forgetting r.

We obtain a bijection with ordinary, unmarked self-converse tournaments having rank one at the fixed vertex of an anti-involution. This property is independent of which anti-involution is chosen. The correspondence is not asserted for a root that is not fixed by such a symmetry.

## 2. Exact automorphisms, fixed words and mass

Let G=Aut(S,a) and g=|G|. Any automorphism of F(S,a) fixing r preserves O and I and fixes the two endpoints of the unique crossing arc. Conversely any independent pair of rooted automorphisms of the two sides extends. Hence

    Aut(F(S,a),r) = G x G.                              (1)

This is the root stabilizer, not necessarily the full unmarked automorphism group. For example F(K1,a) is C3, whose root stabilizer is trivial but whose full automorphism group has order three.

With the canonical tau, conjugation exchanges the two factors in (1). Its centralizer is the diagonal copy of G. Every automorphism commuting with tau fixes its unique fixed vertex r, so this is also the centralizer in the full unmarked automorphism group:

    |C_Aut(F)(tau)|=g.                                 (2)

There are exactly g anti-involutions fixing r: choose a rooted converse isomorphism O -> I, with the reverse map then forced by involutivity. They form one conjugacy orbit under the root stabilizer.

Put h=H(S). Exactly h Hamilton words of F(S,a) are fixed by reversal followed by tau. In such a word, r occupies the middle position, and the first half chooses exactly one vertex from each tau pair. Its last vertex is in I, because it must point to r. If the first half contained an O vertex, reaching its final I vertex would require the unique O-to-I arc o_a -> i_a. But those two vertices form a tau pair and cannot both occur in the first half. Thus the first half is entirely I. Conversely every Hamilton word of I gives precisely one symmetric Hamilton word, by appending r and the reversed tau-image. Their number is h.

The exact self-converse normalization weight from the inherited mass theorem is consequently

    c(F(S,a)) = h/g.                                   (3)

It is a positive odd integer, also directly because G acts freely on Hamilton words of S. Together with the companion count H(F)=h^2+2q^2, (1) yields the full root-stabilizer path-orbit count (h/g)^2+2(q/g)^2. Here q is the two-bin endpoint response defined in the companion packet; G acts freely on its full-vertex path covers, so g divides q too. Formula (3) does not substitute the root stabilizer for the full group in the ordinary unmarked Hamilton-path quotient.

## 3. Closed all-order restricted mass and exact sampling fraction

Let L_t be the sum of c(T) over the ordinary unmarked rank-one self-converse output classes of order 2t+1. By the bijection and (3),

    L_t = sum_[(S,a)], |S|=t H(S)/|Aut(S,a)|
        = t * 2^((t-1)(t-2)/2),     t>=1.              (4)

To prove the last equality directly, distinguish a Hamilton word of S and normalize its vertices into word order. Its t-1 consecutive arcs are forced; the other binom(t,2)-(t-1) arcs are independent. There are t choices for the marked position a. Free automorphism action makes each rooted class contribute exactly H(S)/|Aut(S,a)| normalized objects.

The full self-converse normalized universe at order 2t+1 has size 2^(t^2). Thus the exact fraction occupied by this entire rank-one physical family is

    t * 2^(-(t^2+3t-2)/2).                             (5)

Every member excludes all counts 7p by the companion rank-one square/norm theorem. Equation (5) therefore measures a completely described arithmetic-obstructed portion of the normalized universe. It is not the fraction among uniformly chosen unlabelled classes or uniformly chosen labelled tournaments.

At t=1 the family is C3 and has weight one out of the full order-three mass two. At t=4 its weight is 32 out of 65,536. The empty side t=0 is excluded: it gives a singleton of rank zero, not rank one.

## 4. Strong outputs and A000568 generating functions

The output F(S,a) is strong exactly when a lies in the terminal strong component of S. If some vertex of S cannot reach a, its O copy cannot reach the only exit o_a -> i_a. Conversely, when every S vertex can reach a, every O vertex reaches o_a, and i_a reaches every I vertex by converse symmetry. Together with I -> r -> O these paths make the whole output strong. In a tournament, the reachability condition is exactly membership of a in the terminal strong component.

Let

    A(x)=sum_(t>=0) a_t x^t, a_t=A000568(t), a_0=1,
    R(x)=sum_(t>=1) r_t x^t,

where r_t counts ordinary rooted tournament isomorphism classes, and let C(x) count the strong rank-one output classes by side order t (the actual output order is 2t+1).

The ordered strong-component decomposition is unique. Marking one vertex marks exactly one strong component. If R_str(x) counts rooted strong tournament classes, the blocks before and after that component are arbitrary tournaments. Hence

    R(x)=A(x)^2 R_str(x),
    C(x)=A(x) R_str(x)=R(x)/A(x).                       (6)

The second equation uses the terminal-component condition, so the post-marked block is empty. These are ordinary generating functions of isomorphism classes, not exponential generating functions. In particular R(x) is not x A'(x): at t=3 there are four rooted classes, not six, because the directed triangle has one vertex orbit rather than three. This supplies a concrete counterexample to replacing orbit-pointing by ordinary differentiation in A000568.

For the weighted analogue define

    M(x)=1+sum_(t>=1) 2^((t-1)(t-2)/2) x^t.

Its coefficient is the sum of H(S)/|Aut(S)| over unmarked tournament classes. Both H and automorphism order multiply across an ordered sequence of strong components. Orbit-stabilizers give the rooted weighted series x M'(x), in contrast to the unweighted R(x). The exact total c-mass of the strong rank-one output classes is therefore

    C_mass(x)=x M'(x)/M(x).                            (7)

The first coefficients, for t=1,2,3,4, are:

    r_t:                   1,2,4,12;
    strong rank-one count: 1,1,2,7;
    full rank-one c-mass:  1,2,6,32;
    strong rank-one mass:  1,1,4,25.

These initial checks are consequences of the all-order bijection and factorization, not an extrapolated recurrence.

## 5. Sources and verification

The companion `../binary-rigidity/PROOF.md`, frozen SHA256 `65102468958b198f8538022288238c40d4e3ca33684cf0d419bc86b0993681a9`, proves the rank-one source form, square count and 7p obstruction. The anti-involution conjugacy and exact c(T) interpretation are in `../CONVERSE_MASS.md`, SHA256 `07f4b20dde9f9d23b7711f65bd484d7468bdbef3cc6da0cc2888aa91d402a620`. The ordinary path-orbit normalization comes from the named S196 symmetry proof and is reproved for the relevant rooted objects above. The present bijection, root-stabilizer product, explicit fixed-word count, restricted mass and strong-class series are derived here. No literature-priority claim is made.

The standalone verifier enumerates all labelled sides through order four and canonicalizes their rooted classes. For every distinct rooted class it builds the physical output, independently enumerates its root stabilizer and root-fixing anti-involutions, counts symmetric Hamilton words by half-word permutations, and checks strongness. It verifies the weighted and unweighted series identities against those physical inventories and gives C3 as the explicit full-group-versus-stabilizer boundary. Normal and Python -O receipts must match. No code assertion is needed for validity.
