# A perfect isolator for a marked birth cannot be inherited from an unmarked child

The attached [FMCAD paper](https://www.cs.cmu.edu/~mheule/publications/FMCAD22.pdf) fixes a transitive subtournament by moving its vertices into prescribed positions, then repeats on the complement (Section IV). It enforces one representative per *unmarked* orbit by ExactOne (Eq. 8), and extends isolators on disjoint free blocks (Section VI-A). A reverse Collatz lift has stronger data: a fixed old tournament, named ports, and a newborn whose deletion recovers that old tournament. The following describes exactly which of those operations still apply.

## Attachment cube and rooted orbit count

Fix a labelled tournament T and an ordered tuple P of old ports. A new vertex is specified by a binary word x on the old vertex set, with bit 1 for vertices it beats. If the newborn and all old ports remain distinguished, two such extensions are isomorphic precisely when their words are in the same orbit of G=Aut(T;P), the pointwise stabilizer of the old ports. Necessity follows by deleting the marked newborn; sufficiency extends the automorphism by fixing the newborn. This also applies to a succession of births that must preserve named prior vertices: use the current stabilizer, not the stabilizer of the unmarked new tournament. It does not justify identifying words by count, degree, or position in the same unrestricted child class.

Thus the exact number of distinct marked births of a fixed T and P is

    (1/|G|) sum_{g in G} 2^(number of cycles of g on V(T)).

More generally suppose the allowed words prescribe specific bits on named old ports and otherwise are unrestricted. Then the same formula holds on the remaining vertices, whose permutations are induced by G. A fixed g admits exactly one bit per free cycle. If further invariant conditions constrain the word, replace the power of two by the number of constant-on-cycles assignments satisfying those conditions. This is Burnside counting; the input cut labels and deletion witness are essential to its application.

## No additional universal unit on the attachment bits

Every unrestricted coordinate can be zero in one admissible word and one in another: take the free bits all zero or all one, with any already prescribed ports unchanged. Each such assignment is fixed by the stabilizer. Therefore an isolator that covers every marked orbit cannot enforce even a single new attachment unit. This is the opposite of TT-fixing, which gains many units by moving the **underlying** old vertices freely before the attachment cube is fixed.

It is possible to recover the paper's constructive unit reduction when the old tournament is not fixed, but port incidences are. Partition all old nonports by their binary adjacency pattern to the ordered ports. Permutations inside a type preserve every incidence to the ports. Within each class separately, fix disjoint transitive subtournaments at the positions allowed in that class, as in Section IV; all cross-type arcs remain variables. Thus for type sizes n_α, at least (1/2) sum_α sum_{j=1}^{n_α} floor(log₂ j) old arc units can be fixed, without moving or permuting any port. For a fixed number p of ports the number of types is at most 2^p, so this is (1/2)n log₂ n−O_p(n) for n old nonports. The bound follows also from sum n_α log₂ n_α ≥ n log₂(n/2^p), plus the paper's elementary lower bound per class. These are units in the *underlying tournament search*, not in a fixed tournament's attachment word.

## Small exact hostile and the extension-graph consequence

With old core TT₂, all four attachment words are distinct as marked births since Aut(TT₂)=1. Three of them give the same unmarked transitive triangle (the newborn can be first, middle or last); only one gives the cycle. This is already a 3-to-1 failure if the isolator strips the deletion witness. Likewise a canonical Hamilton path on an unmarked tournament need not extend any old chosen order.

For the elementary illustration of stabilizer compression, old C₃ with no fixed port has four marked attachment orbits, with weights 0,1,2,3; its group has cycle counts 3,1,1, giving (8+2+2)/3=4. If one old vertex is a named port, the group becomes trivial and there are eight distinct attachment words (before fixing any port bit). Merging them by their old unmarked symmetry would silently drop legal continuation candidates.

For a proof search use two levels: isolate the *actual marked state* by its stabilizer, and only then apply a contextual response quotient proved to preserve the desired future class of operations. The paper's ExactlyOne orbit condition remains valid at the first level; its gain in an unmarked model count does not prove any safety at the second. For an exact reverse Collatz block, keep the newborns/ports needed at the endpoint as part of that state, or certify an explicit replacement transporting them.
