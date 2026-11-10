# The first strong tournament carrying simultaneous seven and eleven

Session240, 2026-11-10. **PROVED from the orbit argument and independently audited FINITE-EXACT inputs.**
Here H counts ordered directed Hamiltonian paths. An order-77 symmetry means
an actual automorphism whose order is 77, not merely 77 dividing H or the
order of an unrelated incidence group.

## Theorem: minimum order and its exact path floor

A strong tournament admitting an order-77 automorphism has at least 19
vertices. The bound is attained. Among all such tournaments on exactly 19
vertices, the minimum H is

    3,049,061,631 = 77 * 39,598,203.

For a specified permutation with cycle structure (7)(11)(1), exactly 120
labelled invariant tournaments attain this minimum. This is a count of
orientations on specified labels, not isomorphism classes. No all-order
minimum or full automorphism-group classification is asserted.

## Structural proof

A permutation of order 77 has cycle lengths in {1,7,11,77}, with least common
multiple 77. Below order 77 it must contain both a 7-cycle and an 11-cycle.
Thus below order 19 the only possible faithful profile is (7)(11), on 18
vertices. On any pair of the two orbits the diagonal generator acts
transitively: its pair orbit has length lcm(7,11)=77, the full Cartesian
product. Hence all cross arcs agree, giving an ordinal sum of the two
internal tournaments, which is not strong. Orders below 18 cannot support
the required permutation at all.

At order 19 the profile is necessarily (7)(11)(1). The same coprime argument
makes each moving orbit a module; the singleton is trivially a module.
The induced tournaments are circulants and the three-vertex quotient must
be a directed triangle. Conversely, substituting any seven- and
eleven-circulants and a singleton into a triangle is strong and carries the
product rotation of order 77. This proves the exact minimum order and
completeness of the min-order search universe.

## Count with the retained two-path sidecar

For a tournament A, write a=H(A) and

    x=F2(A)=sum_(empty != U != V(A)) H(A[U]) H(A[V(A) minus U]).

Thus x counts ordered spanning two-path covers. Let b,y be the corresponding
quantities for B. The inherited THM-3121
(`path-cover-walk-content-substitution-kernel`) gives

    H(C3[A,B,1]) = 3ab + xb + ay + xy.                 (1)

For completeness, every path splits uniquely into maximal runs in the three
modules. The singleton is used once, so each other block has one or two
runs. Run contents (1,1,1) have three possible cyclic starting letters;
(2,1,1), (1,2,1), and (2,2,1) each have one. The ordered covers supply the
run assignments. This proves (1) directly and also shows why H alone cannot
be substituted without F2.

The complete finite inputs in ../gaps/probe.json give all eight seven
orientations: six have a=175 and two have a=189, while all have x=1638.
All 32 eleven orientations have b>=92411 and y>=1571526, with both minima
simultaneously attained by exactly ten orientations. These finite minima
are independently audited, not assumed from Paley symmetry. Since every
coefficient in (1) is positive, the minimum is attained precisely by these
6*10 pairs. Either of two directed quotient triangles gives the same
count, hence 120 labelled minimizers. A concrete pair uses connection sets
{1,2,3} mod7 and {2,3,4,5,10} mod11.

## What this connection preserves and loses

Source: an actual cyclic action of composite order77. Target: its orbit
quotient decorated by internal path-cover data. The quotient preserves
strongness here because cross incidences are uniform. Bare orbit sizes and
H discard the run-count coordinate; F2 restores exactly what (1) needs.
The cheap hostile is deleting the singleton: the quotient becomes a
transitive pair and the tournament ceases to be strong. A bare group with
order divisible77 need not contain an element of order77, so that change
of quantifier is forbidden. The result prices the first possible strong
joint-symmetry carrier; it neither supplies arbitrary target H values nor
settles Chan–Pak completeness.
