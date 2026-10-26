# The carousel palette meets every ordinary realization orbit

**Status:** PROVED using the explicitly CITED structural theorem of Boussaïri–Chergui–Ille–Zaidi, [2018, Theorem10, Proposition51 and Theorem24](https://arxiv.org/html/1805.04607v1). This addendum completes the ordinary-orbit scope deliberately left separate in the arithmetic producer. It does not say the palette contains every labelled realization.

**Concurrent credit.** The complete carousel-fibre classification and its hypergraph group/mass were also proved in the incoming [companion packet](https://github.com/eliottcassidy2000/math/blob/c9a8823436ed8c49a3effafd40353880a9909a76/05-knowledge/results/cyclictriple_composition_20261026/README.md), commit c9a8823436ed8c49a3effafd40353880a9909a76. This independently reconstructed addendum integrates that result; it is not a claim of a new carousel classification or priority. The distinct Session68 additions are the recursive ternary stabilizer spectrum, its comparison with the carousel action, and the normal-basis/CRT interfaces.

The cited inputs are: all realizations of one labelled cyclic-triple hypergraph share its strong modules; their quotient realization choices at those nodes are independent; a prime tournament quotient has just itself and its converse as labelled realizations. Our deductions below use these facts, not the erroneous historical assertion that automorphism order is constant on a fibre.

## 1. Orbit completeness

Let m>=3 be odd, R_m the carousel, A=source->C3, and

    H_m = C3(R_m[A,...,A]).

The m order-four blocks are canonical, since R_m is prime. Inside each block the common strong-module tree has one linear node with children of sizes1 and3, and one prime triangle node. A labelled realization therefore has:

- one of two outer quotients, R_m or R_m^op;
- either ordering of the two children in each linear node;
- either orientation of the triangle in each block.

There are exactly 2*4^m such labelled realizations. The 2^m standard binary decorations are fewer than these.

Nevertheless each realization is isomorphic to a standard binary decoration. Normalize each terminal triangle by an isomorphism; its linear-node order leaves precisely the two ordinary types A and A^op. If the outer quotient is reversed, a reflection of its m positions brings it back to R_m, carrying the decorated blocks to their new positions. This gives a word epsilon with T_epsilon in the palette. Between two normalized palette objects an ordinary tournament isomorphism must preserve the oriented R_m quotient, hence be a rotation, not a reflection. Thus the arithmetic lane's binary-necklace classification is exactly the ordinary isomorphism classification of ALL realizations of H_m.

## 2. Hypergraph symmetries and the complete mass check

The hypergraph C3(A) consists of a three-vertex edge and an isolated vertex. Its automorphism group is S3, fixing that singleton. At the prime outer quotient, any hypergraph automorphism sends R_m to itself or its converse by the cited prime realization theorem. Therefore the outer hypergraph automorphism group is exactly the dihedral group of order2m. All internal hypergraphs are identical, so

    Aut(H_m) = S3^m semidirect D_(2m),
    |Aut(H_m)| = 6^m * 2m.

Here D_(2m) means the group of order2m. This is a hypergraph group; it is allowed to have even order. Tournament automorphism groups remain odd.

For a necklace with least period d|m, the arithmetic theorem gives

    |Aut(T_epsilon)|=3^m*m/d.

Its orbit under Aut(H_m) inside the labelled realization fibre therefore has

    |Aut(H_m)|/|Aut(T_epsilon)| = 2^(m+1)*d

elements. If M_d(2) counts primitive binary necklaces of length d, then

    sum_(d|m) M_d(2) * 2^(m+1)*d = 2^(m+1)*2^m = 2*4^m.

This verifies the whole fibre mass, retaining each different stabilizer. Dividing every lift by a single assumed automorphism order would be wrong.

At m=3, the four ordinary types have automorphism orders81,27,27,81 and labelled hypergraph-orbit sizes16,48,48,16. The sum is128. The independent composition auditor enumerates all128 node choices and all1296 hypergraph symmetries, rather than only checking the total.

## 3. Same size and triangle count do not identify the compositional carrier

At36 vertices the carousel m=9 and the depth-two ternary family both have cyclic-triple COUNT1929. Their hypergraphs are nevertheless not isomorphic: their canonical maximal strong blocks have sizes4 (nine blocks) and12 (three blocks), respectively. This distinction follows from the cited common-strong-module theorem.

| Carrier | Ordinary realization types | Full labelled fibre | Tournament automorphism orders |
|---|---:|---:|---|
| Prime carousel R9 with A/A^op blocks | 60 | 2^19 | 3^9,3^10,3^11, with56,2,2 types |
| Two-level C3 tree with A/A^op leaves | 24 | 2^22 | 3^9 through3^13, with2,10,8,2,2 types |

The ternary proof and polynomial appear in the classification lane. Its512 standard binary palettes meet all24 ordinary orbits but do not exhaust its2^22 labelled realizations. The normal-basis field-degree dictionary applies to the cyclic action; the independently checked wreath/odometer hostile explains why it cannot silently replace the larger hierarchical action.

## 4. Scope of the metagraph claim

These are complete classification statements for named realization fibres, not a complete atlas of all prime tournaments or all metagraph structures. Literal complement lines can leave the family, as proved in METAGRAPH_BOUNDARY.md. The extra path marking and its exact ports must be retained for that operation. A good ordinary classifier does not automatically supply a closed line grammar.
