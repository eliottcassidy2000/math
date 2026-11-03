# EC1: six signed polynomials for an edge between two retained ports

**PROVED for arbitrary finite tournament edge pairs.** This carrier records
the difference of two observers. It does not replace the source tournament
or the full eighteen-polynomial observer.

Use unordered spanning covers by nonempty directed paths. For each component
retain its ordered endpoint pair among ports a=0,b=1; use -1 for an unnamed
endpoint. Entirely unnamed components contribute the polynomial exponent.
This is the inherited eighteen-skeleton observer of
[GRADED_OBSERVER](../source_selection_closure_20261103/GRADED_OBSERVER.md)
and [GLUING](../coverage_machine_induction_20261103/gluing/NOTE.md).

Let Delta be the minus(a<-b) observer minus the plus(a->b) observer, with
all other arcs identical. Covers not using that edge cancel. In every
remaining cover, both ports lie consecutively in the same component.
Its endpoint skeleton can only be

    empty, a-start, b-start, a-end, b-end, b-to-a, a-to-b.

The last two coordinates are opposites: a component with both ports as
endpoints and using their edge is precisely that edge, leaving an arbitrary
cover of T-{a,b}. Hence only six polynomial coordinates are independent.
The free order used in the receipt is

    (empty, a-start, b-start, a-end, b-end, b-to-a),

and a-to-b is reconstructed as minus b-to-a. This proves lossless encoding
of every edge difference, without a restriction on the rest of the graph.

## Operations that really close

For each of the eight fixed-port uniform births, the new vertex has
arbitrary incidences with a,b and one common incidence with every other
old vertex. Both outputs still differ only on ab. Thus the six-coordinate
space closes. Inserting the new vertex in a cover can join at most one
anonymous component: a uniform anonymous endpoint cannot simultaneously
provide both an incoming and outgoing incidence. The exact polynomial
operator therefore needs only first derivatives.

The cyclic wrapper G(T)=T->c->d->T similarly preserves the edge-pair
relation and has the inherited second-order differential operator. Signed
restriction of that operator to the six coordinates is exact. Beware that
the inherited positive Counter implementation discards negative entries;
carrier6.py uses a separate signed implementation.

Converse swaps starts and ends and negates the b-to-a coordinate:

    (E,Sa,Sb,Ea,Eb,C) -> (E,Ea,Eb,Sa,Sb,-C).

Physical converse commutes with G after interchanging its two fresh names.
None of this identifies field Frobenius with tournament converse.
Dropping an old port or selecting another old vertex is outside the
six-coordinate closure claim. Such steps use the full source and full
eighteen-polynomial transport, as in golden_transport.py.

## Sharp pointwise linear lower bound, even for strong sources

Take these six actual five-vertex sources, their ab flips, and apply G
to both sides, retaining a=0,b=1:

    (2,0,3,7,15), (6,0,2,7,15), (10,0,3,6,15),
    (18,0,3,7,14), (2,4,1,7,15), (2,8,3,5,15).

The resulting twelve tournaments are strong and all have order7. The
determinant of their six free polynomial differences is

    -2 z^2 (z+1)^2 (z+3)
       (z^3+7z^2+10z+3)(z^3+7z^2+12z+2)
       (2z^3+9z^2+12z+3)^2 (z^4+5z^3+5z^2+3).

It is nonzero. Therefore a common pointwise Q(z)-linear encoding from
which the whole edge response is linearly recoverable needs at least six
coordinates. The construction supplies six. This does not bound arbitrary
nonlinear encodings, a scalar-only observer, or a faithful whole-source
carrier. Polynomial degrees still grow with size.

## H497 and controls

For its neutral reset the empty-skeleton difference is

    z^2(3z^3+28z^2+62z+28),

while the b-to-a coordinate is

    z(z^6+21z^5+140z^4+350z^3+301z^2+63z+1).

The latter is exactly the unordered cover polynomial of the transitive
seven-vertex base: its coefficients are Stirling numbers of the second
kind. Neutral total H does not make the signed rooted differences vanish;
the opposite root effects in RC1 are retained by these polynomials.

Independent permutation-and-cut enumeration checks all37 edge pairs on
orders2..4. Subset/partition enumeration checks H497 under all eight births,
literal G, and literal converse. Six iterated G steps compare full and
restricted signed transfer. carrier6.json contains every reduced operator,
the full rank matrix and the twelve literal strong witnesses.
