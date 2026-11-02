# Sharp weak determinant-gap stability at order six

Status: PROVED, independently audited (Session138 audit lane).

For every weighted tournament A of order six, with the notation of
SHARP_LOW_ORDER_STABILITY.md,

    delta=H-u
      >= tau+(37/48)V+(3/4)C4+C5+C6+3M33.           (1)

The coefficient37/48 is sharp. Equality with V>0 in the unreserved
inequality delta>=37V/48 occurs exactly for three fair ordinal pairs.

The exact order-six cycle expansion gives

    delta=V+tau+C4+C5+C6-M22+M23-M24+3M33+M222.

Because every 4-cycle has a complementary pair, M24<=C4/4. Hence

    delta-tau-37V/48
      >= (3/4)C4+C5+C6+3M33+I,
    I=M222+M23+(11/48)V-M22.

The triangle-variance lemma in TRIANGLE_VARIANCE.md implies M23>=P21(q),
the sum of products on a two-edge wedge and a vertex-disjoint edge.
The resulting lower polynomial

    I_low=M222+P21(q)+(11/48)V-M22

is affine in each q_e in [0,1/4]. At an endpoint undirected graph on six
vertices, with edge count m, matching counts m2,m3 and wedge-plus-edge count
p21, it becomes

    192 I_low=3m3+3p21+11m-12m2.

`graph6.cpp` exhausts all 2^15=32768 labeled simple graphs, with no filters.
The minimum is zero, attained at76 graphs; `graph6.json` is the receipt.
The independent audit uses a different edge-toggle pattern-count kernel
and checks the same complete universe. Multiaffine interpolation therefore
proves I_low>=0 for the full real q-box, establishing (1).

Three fair ordinal pairs have H=1,u=27/64,V=3/4, giving
(H-u)/V=(37/64)/(3/4)=37/48. Conversely, unreserved equality with V>0 in
(1) forces tau=0 and hence H=1. The H=1 classification and scalar minimum
then force three fair pairs.

A separate larger exact check (`weak6.py`, `weak6.jsonl`) certifies the
weaker displayed delta>=tau+37V/48 directly by degree-two Bernstein tensors
on all56 binary orientation orbits. It is corroborative, not required by
the small graph certificate proof above.

    g++ -O3 -std=c++17 graph6.cpp -o graph6
    ./graph6 > graph6.json
