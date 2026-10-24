# A structural classification of nonrigid regular tournaments on eleven vertices

24 October 2026. Status: PROVED exhaustive structural reduction; FINITE-EXACT
full automorphism and ordinary isomorphism certification of nineteen small
candidate presentations, representing eighteen distinct classes. The count of
rigid regular-eleven classes and the almost-regular-ten parent count are not
computed by this packet. These are local input boundaries, not global research
status claims.

## 1. Result and the improved prediction interface

Ordinary orientation-preserving isomorphism is used throughout. Converse is
not an equivalence. A regular eleven-vertex tournament has outdegree five at
every vertex and one stable directed-refinement cell. If its full automorphism
group has q vertex orbits, its refinement error is e=11(q−1).

There are exactly eighteen nonrigid regular-eleven classes:

    4 vertex-transitive classes, error0, full Aut orders11,11,11,55;
    7 classes with vertex orbits5,5,1, error22, full Aut order5;
    5 classes with vertex orbits3,3,1,1,1,1,1, error66,
      four with full Aut order3 and one with full Aut order9;
    2 classes with vertex orbits3,1,1,1,1,1,1,1,1, error88,
      both with full Aut order3.

In particular error44 is impossible. If r denotes the number of rigid
regular-eleven classes, the entire regular-eleven positive-error polynomial is

    7 t^22 +5 t^66 +2 t^88 +r t^110.

The nonrigid classes contribute78 vertex-pointed classes and total error660.
Let P10 be the number of ordinary ten-vertex tournaments with degree multiset
4,4,4,4,4,5,5,5,5,5. The unique regular extension establishes a bijection from
these almost-regular-ten classes to vertex-pointed regular-eleven classes.
Therefore

    P10 =11r+78,
    a_regular11 =r+18 =(P10+120)/11,
    E_regular11 =110r+660 =10P10−120.                       (1)

The source-pinned nonregular audit proves E_s,11=E_regular11+9166 and
E_11=E_regular11+19480. Substituting(1) yields

    E_s,11 =10P10+9046,
    E_11   =10P10+19360,
    a_11   =(U11−10P10−19360)/11.                          (2)

Thus the regular correction needs only the count P10, not rerooting all of its
regular extensions. This packet does not supply P10 or U11 and makes no
numerical prediction for the unrestricted order-eleven class count. It proves
P10≡1 modulo11. Also E_regular11 is divisible by110 and E_11≡10 modulo110;
this does not imply the same modulo110 statement for U11, which includes11a_11.

## 2. Full-orbit restrictions and exhaustive partition list

The vertices here are full automorphism orbits, not stable-refinement cells
and not a chosen equitable presentation. Every automorphism orbit induces a
regular tournament: transitivity within the orbit makes its internal outdegree
constant. It therefore has odd size. Crossings between two full orbits have
constant row and column sums. For sizes a,b the usual double-counting identity
ar=bs applies, so coprime-sized orbit crossings are uniform.

Every tournament automorphism group has odd order. A nonidentity involution
would swap a pair and reverse their edge; Cauchy's theorem rules out an even
order group. Equivalently every orbit size is odd, as also proved above.
Since eleven is odd, the number q of full vertex orbits is odd. The complete
odd partitions of eleven, grouped by q, are

    q=1: 11;
    q=3: 9+1+1, 7+3+1, 5+5+1, 5+3+3;
    q=5: 7+1^4, 5+3+1^3, 3+3+3+1+1;
    q=7: 5+1^6, 3+3+1^5;
    q=9: 3+1^8;
    q=11: 1^11.

The q=11 case is exactly rigidity and remains the free count r. The following
eliminations are exhaustive among the remaining positive-error patterns.

- 9+1+1: a fixed vertex has degree0,1,9,or10, never5.
- 7+3+1: the fixed vertex has degree0,3,7,or10, never5.
- 5+3+3: a vertex in either three-orbit has internal degree1, either zero or
  five outgoing neighbors in the five-orbit, and between zero and three in
  the other three-orbit. Its degree is at most4 or at least6, never5.
- 7+1^4: a fixed vertex has at most3 outgoing neighbors if it loses to the
  seven-orbit and at least7 if it beats it. Neither allows degree5.
- 5+3+1^3: the degree d of a fixed vertex within the fixed triple satisfies
  5b+3c+d=5 with b,c in{0,1}, so d is0 or2. But the sum of the three internal
  degrees is3, not an even number. This is impossible.
- 3+3+3+1+1: a fixed vertex has degree0 or1 modulo3 according to its edge
  within the fixed pair, whereas5 is2 modulo3.
- 5+1^6: a fixed vertex must be a source or a sink of the induced six-vertex
  tournament, according to whether it loses to or beats the five-orbit.
  A tournament has at most one source and at most one sink, not six.

Consequently the only positive nonrigid possibilities are5+5+1,3+3+1^5,
and3+1^8. In particular q=5 and error44 cannot occur. The checker generates
all twelve odd partitions and independently tests the elementary degree
obstructions, including every six-vertex tournament in the last tiny control.

## 3. Three orbits: seven classes of error22

Both five-orbits induce the unique regular five-vertex tournament C5. For
completeness, fix a vertex v, order its outneighbors a→b and inneighbors c→d.
Regularity forces b→c,d, a→c, and d→a, producing the cyclic tournament.
Its full group consists of translations: a vertex's successor is the unique
source in its ordered two-vertex outneighborhood, so every automorphism
preserves successor. The fixed singleton s has degree5 and a
uniform relation to each five-orbit. Thus name them intrinsically by s→A and
B→s. If M records A→B, its row and column sum d must satisfy

    2+d=5, hence d=3.

Any full automorphism preserving these orbits restricts to C5^2. Since both
orbits are transitive, this action is nontrivial and projects nontrivially to
both factors. A degree-three matrix has no nonidentity coordinate-axis
stabilizer: such an element would force identical rows or columns and hence
uniformity. A nontrivial translation stabilizer is therefore exactly a
noncoordinate order-five line, not the whole two-dimensional group.

For a slope m in{1,2,3,4}, such invariant matrices are exactly

    M[x,y]=1 iff y−mx belongs to S,

where S is a three-element subset of F5. There are4*binom(5,3)=40 matrices.
Distinct slopes cannot give the same matrix, because two different stabilizer
lines generate the whole translation group and would force uniformity. Their
C5^2 orbits have size5, giving exactly eight presentation orbits.

Recognition of actual FULL vertex orbits must still be checked: these are
regular tournaments and refinement never exposes the proposed blocks. The
checker enumerates all forty matrices and eight translation orbits, then
computes the unrestricted eleven-vertex automorphism groups. Seven have
Aut order5 and full vertex orbits5,5,1. The eighth has Aut order55 and is
vertex-transitive. It is a required hidden-symmetry control, not an eighth
error22 class. The full ordinary isomorphism checks show all eight candidates
are distinct, and identify the order55 control with the vertex-transitive
class in Section6.

Completeness follows from the FULL-orbit argument at the start: every regular
class with q=3 supplies exactly this intrinsic5+5+1 decomposition and appears
in the forty-matrix universe. No regular class having some arbitrary C5
presentation is assumed to have only the displayed block stabilizer.

## 4. Seven orbits: five classes of error66

The surviving orbit pattern is two triangles A,B and five fixed vertices F.
For any v in F, the uniform relations to A and B contribute0,3,or6 to its
outdegree. Its degree inside F lies between0 and4. To have total degree5,
it must beat exactly one triangle and have degree2 within F. Thus F induces
the unique regular C5.

Let d be the constant row sum from A to B, and let S be the fixed vertices
beaten by A. B must have the opposite relation to each fixed vertex. A's
regularity gives

    1+d+|S|=5, or |S|=4−d.

Swapping A and B changes d to3−d and S to its complement. The unique canonical
normalization is therefore d=0 or1. This uses crossing degree, not a supposed
stable color distinguishing the two triangles.

If d=0, |S|=4. The C5 rotations have one orbit on these subsets. The crossing
is uniform and imposes no coupling on the triangle rotations. This gives
one candidate with subgroup C3^2.

If d=1, |S|=3. The ten subsets form two free C5 rotation orbits. The crossing
is a three-by-three permutation matrix. In oriented triangle coordinates its
matching is x→epsilon*x+c. Independent internal C3 rotations remove the phase
c but not the sign epsilon. Each of the two attachment orbits therefore gives
two candidates, four altogether, each with a displayed diagonal C3 subgroup.

These arguments exhaust every q=7 class: the full orbits force the internal
C5 and the complete d/attachment/matching data just described. The checker
computes the full ordinary automorphism group of all five candidates, finding
one group of order9 and four of order3, always with precisely the intended
seven vertex orbits. Pairwise unrestricted isomorphism checks separate all
five. Their errors are all66, totaling330.

## 5. Nine orbits: two classes of error88

Here the only nontrivial orbit is a triangle C; the other eight vertices F
are individually fixed. Each fixed vertex uniformly beats or loses to C.
Its degree within F is respectively2 or5. Since F has28 internal edges,
there must be exactly four of each degree. Call these sets L and H. All
vertices of L beat C, and C beats all vertices of H.

The four low vertices have total degree8 inside F, of which six edges lie
inside L. Hence there are exactly two edges L→H. If M records these edges,
then the internal degrees satisfy

    deg_L(v)=2−row_sum_M(v),
    deg_H(w)=1+column_sum_M(w).

Thus every internal L degree is at most2 and every internal H degree at least1.
There are only two possibilities for L: the unique strong four-tournament,
with scores1,1,2,2, or a directed triangle pointing to a sink, with scores0,2,2,2.
Similarly H is strong4 or a source pointing to a directed triangle, with
scores1,1,1,3. These elementary four-vertex classifications are also exhausted
by the checker over all64 labelled tournaments.

The nonstrong L case cannot occur with F fixed pointwise under full Aut. Its
three cycle vertices have internal degree2, so all have row sum zero in M.
They lose uniformly to H, beat C uniformly, and have the same relation to
their fourth L vertex. Rotating this triangle is consequently an automorphism
of the whole eleven-vertex tournament, contradicting that they are fixed.
The nonstrong H case is excluded dually: its three cycle vertices have column
sum zero and their rotations extend independently. This is an essential
full-orbit exclusion, not a claim that those raw regular constructions do
not exist. They belong to a smaller-q sector instead.

Therefore L and H must both be the unique rigid strong four-tournament. The
two degree-one L vertices each require one outgoing edge to H; the two
degree-two H vertices each require one incoming edge from L. All other entries
of M vanish. There are exactly two perfect matchings between these named pairs.
No internal quotient factor is missing because strong4 has trivial full Aut.

The checker builds both candidates, computes their unrestricted full groups,
and verifies that both have group C3 with full vertex orbit sizes3+1^8. The
two are ordinarily nonisomorphic. Both errors are88, totaling176.

## 6. The four zero-error nonrigid classes

A regular-eleven class has error zero iff its full group is vertex-transitive.
Such a group has order divisible by11 by orbit-stabilizer. Cauchy's theorem
therefore supplies an element of order11, necessarily an eleven-cycle. In
coordinates along this cycle, the tournament is a cyclic Cayley tournament:
choose one member of each pair{d,−d} for d=1,...,5. Exactly32 connection sets
occur, with no target census as input.

The nonzero multipliers in F11 give a permitted relabeling group. Their orbits
on those32 sets have sizes10,10,10,2. This produces four candidates, but a
multiplier quotient alone need not be trusted as the full ordinary isomorphism
quotient. The checker therefore independently tests all six pairs for full
ordinary isomorphism and finds none. Full automorphism orders are11,11,11,55;
all four are vertex-transitive. Consequently there are exactly four classes.
No circulant-isomorphism theorem or unverified CI assumption is used.

The Aut55 class is the hidden transitive case from Section3, and is counted
once, not twice. The combined list has exactly eighteen distinct nonrigid
classes; all153 pairs are independently tested for ordinary isomorphism.

## 7. Pointed mass, inheritance, and the remaining scalar

Deleting any vertex of a regular-eleven tournament yields an almost-regular
ten-vertex tournament with five degrees4 and five degrees5. Conversely the
unique regular extension is obtained by directing each degree-four vertex
to the new root and the root to each degree-five vertex. Parent isomorphisms
extend uniquely over the marked root. Thus P10 counts vertex-pointed regular
classes, equivalently the sum of the full vertex orbit counts q.

For the eighteen nonrigid classes this sum is

    4*1 +7*3 +5*7 +2*9 =78.

A rigid regular class supplies11 pointed classes. The ordinary class total
is r+18, proving the first two formulas(1). The error of the nonrigid part is

    7*22 +5*66 +2*88 =660,

and every rigid class contributes110. This proves the final formula(1).
Equivalently the nonrigid unpointing deficit is

    4*(11−1)+7*(11−3)+5*(11−7)+2*(11−9)=120.

The exact nonregular/SCC constants in(2) are inherited from the pinned audit,
not regenerated by this addendum. No known numerical value of P10, r, the
regular-eleven class count, or A000568(11) appears as input or calibration.

## 8. Reproduction and evidence limits

`check.py` builds forty special5-by5 matrices, five q7 candidates, two q9
candidates, and thirty-two circulant presentations. It computes full ordinary
automorphism groups on nineteen representatives; one appears in two candidate
families, leaving eighteen true classes. It never enumerates all regular-eleven
tournaments or any unrestricted order-eleven universe.

`support.py` is an exact source-pinned copy of the independent main audit's
standard tournament primitives. Only its graph, refinement, group-search,
relabeling, and tiny matrix helpers are used here; no inherited parent atlas
or main audit run is invoked. The mathematical completeness argument is the
full-orbit partition proof above, not the number of generated candidates.
All checks are explicit exceptions and remain live under Python optimization.
