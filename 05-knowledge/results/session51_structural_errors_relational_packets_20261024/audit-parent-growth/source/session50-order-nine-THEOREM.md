# Exact order-nine prediction by structural refinement-error correction

24 October 2026. Status: PROVED structural classification and counting identities;
FINITE-EXACT input certification and parent-only computation. This predicts the
ordinary tournament count at order nine without an order-nine class atlas,
order-nine canonicalization, unrestricted order-nine re-rooting, or using the
known order-nine answer. It is not an all-order scalar recurrence or a proof
that the unmarked blue/black metagraph alone determines the answer.

## 1. Source boundary and inherited results

The input is the complete ordinary order-eight parent atlas with its actual
vertex automorphism permutation actions, reconstructed from the archived
Session48 merged-node table. Converse pairs are separated: ordinary oriented
isomorphism classes are counted, never classes modulo converse. The inherited
Session50 equitable-cell packet supplies the following proved statements:

1. Stable directed color-refinement cells induce odd regular tournaments.
2. Crossings between cells of coprime sizes are uniform.
3. For cells of sizes one and three, the full automorphism group is the
   signed-matching kernel: a triangle cell is transitive precisely when its
   signed-support connected component is balanced. Each triangle in an
   unbalanced component contributes six to refinement error.
4. Every nonregular tournament of order seven and every strong tournament
   of order eight has zero refinement error; the unique non-vertex-transitive
   regular-seven class has error fourteen.
5. Refinement error is additive under ordered strongly connected component
   decomposition, so E(z)=A(z)^2 E_s(z).
6. The three-triangle transitive-majority matching family has sixteen
   ordinary classes, twelve with cyclic automorphism group of order three
   and four rigid. This family count and the signed affine mechanism are
   inherited from Session39, not newly discovered here.

The archived theorem gives their proofs and the exhaustive regular-seven
control. This packet's new step is the exhaustive classification of all
nonregular strong order-nine errors, the almost-regular-eight regular-sector
reduction, and their combination into an exact next-order predictor.

## 2. The root-refinement counting identity

For an ordinary tournament T and a vertex orbit O under Aut(T), let C(O) be
its stable directed color cell. Define

    e(T) = sum over vertex orbits O of (|C(O)| - |O|).

Every summand is nonnegative, because automorphisms preserve stable colors.
A pointed extension is an order-eight parent together with an incoming
attachment subset, modulo that parent's actual automorphism action. Deleting
the distinguished root identifies these states bijectively with ordinary
vertex-pointed order-nine tournaments. Thus the states over a class T are
exactly its vertex orbits, not its individual vertices. Consequently

    U9 = sum over pointed extension states p of |C(root_p)|
       = 9 a9 + E9,                     E9 = sum_[T, |T|=9] e(T).

U9 needs no target isomorphism. The bulk predictor visits each parent
attachment orbit once and computes only the stable root cell size.

## 3. Exhaustive nonregular strong-nine classification

The odd partitions of nine are

    9; 7+1+1; 5+3+1; 5+1+1+1+1;
    3+3+3; 3+3+1+1+1; 3+1+1+1+1+1+1; 1+1+1+1+1+1+1+1+1.

These are possible stable-cell size patterns, not an assertion that every
pattern is realized. The first pattern is exactly the regular-nine sector:
all vertices have outdegree four, so refinement never splits them. Every
regular tournament is strong. Indeed different SCCs occupy disjoint
outdegree intervals, which is impossible when every degree is equal.

For a seven-cell and two singletons, all crossings are uniform. The quotient
is a three-vertex tournament with one distinguished size-seven vertex. The
whole tournament is strong if and only if this quotient is a directed cycle:
the seven-cell is a regular, hence strong, tournament, and a transitive
quotient would have a nontrivial SCC cut. There is exactly one directed-cycle
quotient up to swapping the two singleton vertices. Its stable cells remain
separated by degrees four, one and seven, and its full automorphism action
on the seven-cell is precisely the internal action. Thus there is one strong
nine-vertex substitution class for each regular-seven class, with the same
error. Their total error is fourteen.

A five-cell is the unique regular cyclic tournament C5. Its crossings to
all remaining cells are uniform by coprimality, so its transitive automorphism
action extends independently. Any remaining triple likewise has independent
cyclic rotations. These patterns contribute zero. In the patterns with at
most two triples, the signed nonuniform-crossing support is a forest, so all
its components are balanced and the signed-kernel theorem gives zero error.

It remains to consider three actual stable triple cells. Positive error
requires an unbalanced signed cycle. With only three cells this forces all
three crossings to be nonuniform and their matching signs to have negative
product. Each nonuniform crossing has row sum one or two. Orient a majority
edge from cell A to cell B when each A vertex beats two B vertices. These
three majority orientations form a three-vertex tournament.

If the majority tournament is cyclic, each vertex has internal outdegree one,
one outgoing majority crossing contributing two, and one incoming majority
crossing contributing one. Every vertex therefore has degree four. Such a
tournament belongs to the regular-nine sector, whose stable partition is one
nine-cell, not three cells. It must not be counted again here.

If the majority tournament is transitive, order the three cells by increasing
vertex degree. Those degrees are exactly three, four and five. For every
ordered pair of cells, arcs from the lower-degree cell to the higher-degree
cell form a perfect matching. Hence the tournament is exactly the inherited
three-triangle matching family. Conversely every member of that family has
those three actual stable cells and is strong. All isomorphisms preserve the
three distinct degrees and are rotations within their directed triangles.

For completeness, its class count is elementary. There are three matching
permutations, each affine x -> epsilon x+c over F3. The eight sign assignments
are rotation invariants. Four have positive triangle product and rotation
kernel dimension one; each has three phase orbits. Four have negative product
and zero-dimensional kernel; each has one phase orbit. Thus the sixteen
classes split as twelve transitive-on-cell and four rigid. Every rigid class
has nine singleton vertex orbits in cells of size three, so its error is
nine times (three minus one), namely eighteen. Their total error is seventy-two.

The cases exhaust every nonregular stable-cell pattern. Therefore

    E_s,9 = E_regular9 + 14 + 4*18 = E_regular9 + 86.

Using E_s,k=0 for k<=6 and k=8, E_s,7=14, and [z^2]A(z)^2=3, the ordered-SCC
identity yields

    E9 = E_regular9 + 86 + 3*14 = E_regular9 + 128.       (1)

## 4. The regular-nine sector is an order-eight pointed problem

Call an order-eight tournament almost regular when every degree is three or
four. Its total degree is twenty-eight, so it has four vertices of each degree.
It has a unique extension to a regular order-nine tournament: each degree-three
vertex points to the new root; the root points to each degree-four vertex.
The new root then also has degree four. This attachment subset is invariant
under every parent automorphism, hence gives exactly one pointed state for
each almost-regular parent class.

Conversely deleting any vertex of a regular-nine tournament leaves an almost-
regular-eight tournament. An isomorphism between such deleted parents extends
uniquely over the roots because root adjacency is forced by parent degree.
Thus the almost-regular parent classes are in bijection with regular-nine
pointed states. No additional attachment label or hidden equivalence is needed.

Let P be their number. For a pointed state p, form its unique regular extension,
delete each of its nine vertices, and identify the resulting parent within
the complete almost-regular-eight atlas. Let D[p,q] count deletions leading
to state q. Then D[p,p] is precisely the size of the distinguished root orbit.
Each regular-nine class contributes nine to tr(D), so

    a_regular9 = tr(D)/9.

A regular-nine class with q vertex orbits has one stable cell of size nine,
so its refinement error is 9q-9. Summing over classes gives

    E_regular9 = 9(P - a_regular9).                    (2)

All isomorphism tests in this reduction involve eight vertices, never nine.
The standalone checker also verifies that rows of D are identical exactly
within the recovered child components, each component has trace nine, and
its destination multiplicities sum to nine. This is an independent check
of the trace arithmetic, not a target canonicalization.

## 5. Finite-exact computation and result

Input authentication finds 6,880 pairwise nonisomorphic ordinary parents.
The supplied action group for each parent is complete because it is obtained
by exhausting all permutations within its invariant stable color cells.
The sum of their labelled class sizes 8!/|Aut(T)| is exactly 2^28, the entire
labelled tournament universe at order eight. Distinctness plus this equality
certifies completeness independently of trusting a list-length target.
The orbit-stabilizer mass is an INPUT audit, not a Burnside computation of a9.

Exactly 85 parents are almost regular. The regular-sector calculation makes
765 parent canonicalization calls and inspects only 2,130 candidate parent
labellings. It obtains trace 135, hence fifteen regular-nine classes and

    E_regular9 = 9(85-15) = 630.

Equation (1) gives E9=758, consisting of regular-sector error 630, seven-cell
substitution error 14, triple-holonomy error 72, and lower-SCC padding error 42.

The bulk pass visits 1,716,608 pointed extension states. Stable root cell sizes
and their multiplicities are

    size 1: 1,712,964; size 3: 3,491; size 5: 48;
    size 7: 20; size 9: 85.

Their weighted sum is U9=1,724,582. Therefore

    a9 = (1,724,582 - 758)/9 = 191,536.

A second structural check resolves the error by stable root cell size. A
seven-cell can occupy any of three ordered-SCC positions among two singleton
components, or the unique cyclic three-cell quotient. Each of the two
vertex-transitive regular-seven types contributes one pointed seven-cell state
per placement; the remaining type contributes three, with orbit sizes three,
three and one. Thus there are exactly 4*(2+3)=20 size-seven root states, whose
true orbit-size sum is 4*(7+7+3+3+1)=84 instead of the refinement sum 140. Their
error is 56=42+14. Of the 3,491 size-three states, precisely 4*9=36 come from
the rigid three-triangle classes and have true orbit size one; the others
have orbit size three. This contributes error 72. Size-five and singleton
cells are already exact; the size-nine stratum contributes error 630. Thus

    E9 = E(root size 7) + E(root size 3) + E(root size 9)
       = 56 + 72 + 630 = 758.

The true root-orbit sums by stable cell sizes 1,3,5,7,9 are respectively
1,712,964; 10,401; 240; 84; 135, which sum to 1,723,824=9*a9. This verifies the
correction against the entire computed root-cell histogram, not only its
aggregate weighted total.

The known order-nine sequence value is not an input, correction calibration,
class-list cardinality check, or code constant. It can be compared only after
this calculation. No full order-nine representative census is produced.

## 6. Information retained, costs, and limits

The actual parent vertex actions and attachment subsets are essential input;
class counts, automorphism orders alone, or bare unmarked line-graph adjacency
do not specify this computation. The method retains enough parent-level
incidence to construct pointed extensions, but it avoids full target-level
isomorphism by classifying the gap between refinement and actual root orbits.
The large sector uses no deletion and no isomorphism test. Only 85 regular
pointed states require rerooting, compared with 1,716,608 total pointed states.

The compiled full-refinement control visits all states in about six elapsed
seconds in this run, with 4,131,890 refinement rounds. This is an observed
machine-dependent runtime, not an asymptotic complexity claim. The Python
predictor uses independently structured block refinement and a safe early
exit when the root becomes a singleton. Normal and -O replays plus independent
compiled and separate audit controls are reported in VALIDATION.md.

Order nine closes the next coefficient by a special complete small-cell
classification and a tiny regular sector. At larger orders, larger equal-cell
crossings and more complicated internal regular-cell actions appear. Neither
this theorem nor the positive examples imply that all future corrections are
known, cheap, scalar, or obtainable from the unmarked metagraph alone.
