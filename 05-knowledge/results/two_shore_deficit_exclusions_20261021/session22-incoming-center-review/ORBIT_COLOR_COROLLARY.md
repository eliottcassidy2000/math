# Automorphism-orbit coloring: remaining orders 25 and 35 are vertex-transitive

2026-10-21 UTC. Verdict: PROVED, independently checked direct extension of the inherited score-run obstruction. This is a necessary condition, not an existence assertion, group classification, or census. The reconciliation and incoming certificates remain unchanged.

## Statement

Let T be strong and self-converse with black path-orbit deficit at most six. If a marked Hamilton cycle has exact cut index 5 <= q <= 8 and odd shore size m >= 5, then Aut(T) acts transitively on vertices.

Consequently every still-possible deficit-six center of order 25 or 35 is vertex-transitive. The reconciled parameter rows at those orders are (q,m)=(5,5),(5,7),(7,5), respectively.

## Intrinsic colors and the converse permutation

Color each vertex by its orbit under the full group A=Aut(T). Every automorphism preserves each color individually, so unequal color words certify inequivalent Hamilton words. Since R^q belongs to A, each rotation shore is contained in one A-orbit. Thus a color class comprises c entire shores and has size cm.

If tau is any anti-automorphism, tau A tau^{-1}=A: conjugating an orientation-preserving map by two orientation reversals preserves orientation. Therefore tau induces a permutation pi of the A-orbits. It preserves cardinalities and hence their shore multiplicities c. Also tau^2 belongs to A, so pi is an involution. On orbit-color words the converse/reversal action is reversal followed by pi. No assertion that tau preserves individual marked shores is needed or used.

## Long runs separate four additional black orbits

Assume some color class has 1 <= c <= 4 shores. It is a proper class because q>=5. Along every cut of the marked Hamilton cycle its color appears with cyclic period q, precisely once per shore of that color. Every monochromatic run of that color has length at most c. The same bound holds for its pi-image color, which also has c shores.

Choose a shore A of this color and a directed four-cycle S inside A, available because the odd regular shore has m>=5 vertices. Deleting S leaves a strong complement by the inherited proper-shore deletion lemma. The inherited half-supported prefix argument applies since n=qm>=25.

If c=1, choose any other shore B. Exact odd-shore counting gives a vertex v in B having at least three exits from S or at least three entries to S: otherwise the total number of arcs S->B would force 4k=2m, impossible for odd m. The three associated prefix cuts supply at least two black path orbits, pairwise not converse mates, and their two mates. All four are distinct. Original words have a four-vertex run of the chosen color; mates have a four-vertex run of its pi-image color. Neither is possible in a marked-cycle cut, whose corresponding run bound is c=1.

If 2<=c<=4, choose B to be another shore in the SAME A-orbit. The same counting gives at least three exits or entries. In the exit case, each four-cycle prefix is immediately followed by terminal v in B, so it contains a monochromatic run of length five. Its mate contains a run of length five of the pi-image color. Both exceed c, so the four distinct black orbits are again disjoint from every marked-cycle cut orbit.

In the entry case, apply this construction in the opposite tournament, using the reversed marked Hamilton cycle. Its automorphism group and orbit colors are identical, and its cycle cuts are reversals of the original cycle cuts. Reversing the resulting Hamilton words back to T preserves run lengths, blackness, orbit inequality, and disjointness from cycle cuts. Equivalently one can use suffixes. Thus no unjustified claim that a transformed shore is one of the original marked shores enters the argument.

The marked cuts supply at least q-gcd(2,q) black orbits, so

    Delta >= q-gcd(2,q)+4 >= 8,

contradicting Delta<=6. Therefore every A-orbit contains at least five shores. Since there are at most eight shores total, there is only one A-orbit. This proves vertex transitivity.

## Control scope

orbit_color_control.py enumerates every set partition of the q marked shores for 5<=q<=8, every class of size c<=4, and every potential image class of the same size. For m=5 and7, it checks the cyclic run bound in every marked cut, constructs the abstract four-run or five-run color prefix, and checks that neither its color nor any cardinality-preserving image can have the required run in any cut. It also checks the black-orbit arithmetic and the order25/35 parameter rows.

This is a finite color-word/permutation control, not a construction of tournament Hamilton words or a center census. Literal word existence and four distinct black orbit counting are the all-parameter inherited half-supported prefix theorem; the proof above supplies the new intrinsic coloring and normalizer argument. The finite tests do not replace either analytic ingredient.

Fresh normal and optimized runs agree: 170,181 always-active checks over 41,384 color-class cases. Receipts are orbit-color-normal.json and orbit-color-optimized.json.
