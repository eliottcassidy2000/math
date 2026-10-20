# Fibonacci rectangles, parity-path matchings, and arithmetic fibres

**PROVED elementary:** the explicit matching bijection and weighted transfer
below. The underlying SCC language and its generating function are inherited
from [golden SCC words](../scc_basins_20261020d/arithmetic.md), not newly
discovered here. **FINITE-EXACT:** the adjacent checker verifies the maps and
their stated hostile boundaries. No unrestricted tournament-spectrum,
black-connectivity, or Collatz theorem follows from these arithmetic models.

## 1. The sequence has an actual tournament carrier

Let a,b,c stand for SCCs 1,C3,R4, of orders1,3,4. The unique ordered SCC
decomposition identifies words in these letters with distinct ordinary
tournament types. Write W_n for the words of total order n. Inherited:

    A(z,y)=sum_W z^|W| y^c3(W)
          =1/(1-z-y z^3-y^2 z^4),

where c3 counts directed triangles; b contributes1 and c contributes2.
This is a restricted class, not A000568.

**Parity-path bijection.** On vertices1,...,n put undirected edges
{i,i+2}. This is the disjoint union of the odd-position and even-position
paths. Encode a word, left to right, as follows on its consecutive blocks:

* a: one unmatched vertex;
* b: three vertices, matching the first to the third;
* c: four vertices, matching first to third and second to fourth.

This gives every matching of the two paths exactly once. To decode, inspect
the first remaining vertex i. If unmatched, remove a. Otherwise it is
matched to i+2. If i+1 is unmatched, remove b; otherwise its only available
partner is i+3, so remove c. In each case no edge leaves the removed block.
Induction proves both existence and uniqueness. Each matching edge is one
unit of c3. The carrier, positions, and matching are retained, not merely
the number of edges.

Let f_0(y)=0, f_1(y)=1 and f_(r+1)(y)=f_r(y)+y f_(r-1)(y). The matching
polynomial of a path on r vertices is f_(r+1)(y), by examining its first
vertex. Independence of the two parity paths proves

    sum_(W in W_(2m)) y^c3(W) = f_(m+1)(y)^2,
    sum_(W in W_(2m-1)) y^c3(W) = f_m(y) f_(m+1)(y).       (1)

At y=1 these give Fibonacci squares at even orders and

    1,2,6,15,40,104,273,... = F_m F_(m+1)

at odd orders. Also

    sum_(j=1)^m F_j^2=F_m F_(m+1),

since F_m F_(m+1)-F_(m-1)F_m=F_m^2. Thus the user's square sum is literally
the order-thirteen class count of this SCC language. The odd-order sequence
obeys s_m=2s_(m-1)+2s_(m-2)-s_(m-3), m>=4. Equation(1), rather than a fit,
determines the entire sequence.

### The split 63+210 retains a marked edge

At n=13 the odd path has seven vertices and21 matchings; the even path
has six vertices and13 matchings. Mark its second edge {4,6}. Matchings
containing this edge leave paths on one and three vertices, giving
1*3=3 choices on the even path. Therefore

    containing {4,6}: 21*3=63,
    omitting {4,6}:   21*(13-3)=210,
    total:           21*13=273.                            (2)

This is an explicit partition of tournament types through the bijection,
with a retained position. Other marked edges give other splits: the first
even-path edge gives105+168. An equality with the seven-port bank size105
does not identify that bank with these matchings or preserve its S7 action.

### Why the same rectangle does not settle Hamiltonian counts

For a word W with b triangle atoms and c four-vertex atoms, ordered SCC
multiplication gives H(W)=3^b 5^c. Consequently the sum of Hamiltonian
counts over this restricted class has generating function

    sum_W H(W) z^|W| = 1/(1-z-3z^3-5z^4).                 (3)

The independent parity-path factorization would require the two-edge c
block to have weight3^2=9; its actual weight is5. Thus Hamiltonian weight
introduces an interaction between the rows. The pair of order-six words
bb and caa has the same matching cardinality c3=2, but H=9 versus5.
Keeping the full matching and its unique block decoding preserves H;
keeping only its cardinality does not.

The known black moves also have precise images. The inherited interchange
ab<->ba moves one matching edge between the parity rows, with its stated
outside-context/activity hypotheses. The inherited R4 dissolution c<->ba
deletes one of two matching edges and preserves feedback b+c, while c3
falls by one. Those moves are paths in the ambient tournament metagraph,
not an assertion that every elementary matching toggle is a black edge.
See [fixed-inventory normalization](../scc_basins_20261020d/README.md).

## 2. The commutative half-plane has a fixed-point correction

The faithful operation carrier consists of triples (a,b,a+b) or (a,b,ab),
with positive integer parents. Swapping the parents gives an involution;
the fundamental region a<=b has boundary a=b. The boundary outputs are
doubles and squares. The exact fibre counts are

    unordered positive sum pairs for n: floor(n/2),
    unordered positive product pairs for n: (tau(n)+[n is square])/2.

The extra term is the fixed-point contribution, not a removable cosmetic
diagonal. The older [operation-fibre result](../arithmetic_braids_20260917_summand.md)
already proves that deleting equal-parent witnesses removes exactly the
doubling and squaring forests inside their respective ambient graphs.
We retain that corrected definition; ordinary reachability in a shadow
graph does not certify that both parents are available.

Keeping both outputs S=a+b and P=ab reconstructs the unordered pair:

    (a-b)^2=S^2-4P,
    {a,b}={(S-sqrt(S^2-4P))/2,(S+sqrt(S^2-4P))/2}.

The discriminant must be a nonnegative square of the correct parity and
the smaller root positive. This is a precise sum/product coordinate
system. Sum alone merges(1,4) with(2,3); product alone merges(1,6) with(2,3).
Nor do sum and product suffice for triples: (1,6,6) and(2,2,9) both have
sum13 and product36. More inputs need the missing symmetric coordinates.
This is a concrete analogue of retaining external defect ports, not a
proof that those ports can be replaced by arithmetic scalars.

## 3. Infinite half-planes, rational trees, and the missing braid height

The shears L(a,b)=(a+b,b), R(a,b)=(a,a+b) have disjoint images a>b and b>a.
They preserve gcd. Their inverse is subtract the smaller coordinate from
the larger, strictly decreasing a+b until the diagonal (g,g). Thus every
positive pair has a unique gcd root and a unique L/R address. The coprime
component is the usual rational-pair tree rooted at(1,1). Commutativity's
reflection exchanges the branches; throwing away its address loses the
order of additions. These statements are proved by the explicit inverse.

This supplies a genuine self-similar half-plane decomposition. It does not
supply a paradoxical finite dissection into congruent pieces: the maps
are shears, the region has infinite extent, and finite height cutoffs are
not preserved. Banach--Tarski's original theorem distinguishes dimensions
at least three from the plane for finite congruent decompositions;
see the [primary paper and abstract](https://pldml.icm.edu.pl/pldml/element/bwmeta1.element.bwnjournal-article-fmv6i1p27bwm).
No use of that theorem is made in the finite counting arguments above.

The braid/Farey connection also already has an exact correction in
[THM-2646, braid-three modular central pullback](../../../01-canon/theorems/THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre.md):
the modular braid image must be supplemented by its integer exponent
sum. The Farey/modular image alone forgets central full twists.
For the elementary matrices

    X=[[1,1],[0,1]], Y=[[1,0],[-1,1]],

one checks XYX=YXY and (XY)^6=I. The three-braid (sigma1 sigma2)^6
therefore has trivial matrix image but is nontrivial: its exponent sum
is12. Its closed braid has three components with pairwise linking number2,
whereas the identity braid closes to the unlink. Thus even the exact
matrix does not determine the closed link. This explicit small control
reuses the canon's missing-height mechanism, rather than claiming a new
braid classification. Distinguish this Y from the positive shear R=X's
transposed matrix: here Y=R^(-1).

The practical transfer to tournament research is to keep operation words,
ports, and the corresponding integer information when quotienting. The
Fibonacci matching bijection succeeds because it has an explicit inverse;
the Hamiltonian weight obstruction shows exactly where a smaller scalar
state stops preserving the target.

## Reproduction

Run `python 05-knowledge/results/rotation_owners_spectrum_20261020e/arithmetic_check.py`
and the same command with `python -O`. The checker enumerates both sides
of the matching bijection through order18, tests every preserved triangle
weight, checks literal Hamiltonian counts through order8, and reproduces
the marked-edge splits, operation fibres, Euclidean inverses and braid
matrix/linking hostile. These finite universes accompany the all-order
proofs and do not substitute for them.
