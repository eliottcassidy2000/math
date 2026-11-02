# Why 56 occurs: regularity, cycle deficits and the next missing coordinate

**PROVED elementary identities; known five-cycle extremal bound explicitly
CITED below.** No novelty priority is claimed for the identities. Their use
here is to explain the regular-seven module response and locate its first
failure at order nine. Counts are raw directed cycles up to rotation, with
orientation retained. M33 counts unordered pairs of vertex-disjoint triangles.

For every regular tournament of odd order n,

    C5 + 2 M33 = n(n^2-1)(n-3)(5n^2-12n-41)/2880.       (1)

The values at n=3,5,7,9,11 are 0,2,56,384,1584. This is an all-order
identity, not a polynomial fitted to the three regular seven-vertex classes.

## Proof through arc loads

Write d=(n-1)/2 and let t_e count the directed triangles containing arc e.
Each vertex is in t=d(d+1)/2=(n^2-1)/8 triangles: summing the outdegrees
of its d out-neighbours gives d^2, of which d(d-1)/2 remain within that
out-neighbourhood. Thus the total number of triangles is c=nt/3.

Let A be the adjacency matrix and S=A-A^t=2A-J+I. Regularity gives S1=0.
For e=u->v, the common out- and in-neighbour counts both equal d-t_e;
the number of intermediate vertices u->w->v is t_e-1. Consequently

    (S^2)_(uv)=4t_e-n,    (S^2)_(uu)=-(n-1).

Using sum_e t_e=3c and n(n-1)/2 arcs gives

    tr(S^4)=32 sum_e t_e^2-n^4+n.                       (2)

A closed walk of length five in a tournament is a simple directed cycle:
otherwise splitting at a repeated vertex produces two shorter closed walks,
each of length at least three, which cannot total five. Hence tr(A^5)=5C5.
On the constant vector A has eigenvalue d; on its complement A=(S-I)/2.
Odd traces of the real skew matrix S vanish, so expansion gives

    tr(A^5)=((n-1)^5+1-n+10n(n-1)-5tr(S^4))/32.

Combining this with (2),

    C5+sum_e t_e^2=K(n):=n(n^2-1)(n^2+11)/160.          (3)

Two distinct triangles share zero, one or two vertices; sharing two means
sharing their common arc. Inclusion-exclusion therefore gives

    M33=binom(c,2)-n binom(t,2)+sum_e binom(t_e,2).

Equation(3) and sum_e t_e=3c yield

    C5+2M33=K(n)+c^2-c-nt^2,

which factors as (1). The argument includes n=1,3, with zero short-cycle
counts. QED.

## What the path-cover observer retains

In the inherited odd-cycle-collection formula (THM-002, `THM-002-ocf.md`),
a collection of disjoint odd cycles contributes 2 to the power of its
number of cycles. Applied to each colour class of r labelled slots, it gives

    E_r(T)=sum_collections 2^(#cycles) r^(n-sum_cycle(|cycle|-1)).

E_r counts ordered spanning r-path covers allowing empty slots. The first
three possible coefficients for regular n>=5 are therefore

    E_r = r^n + n(n^2-1)/12 r^(n-2)
          + n(n^2-1)(n-3)(5n^2-12n-41)/1440 r^(n-4)
          + terms of lower degree.                    (4)

The source is an actual regular tournament; the map is its cover polynomial;
the preserved predicate is the number of ordered path covers. The map loses
the distinction between one five-cycle and two disjoint triangles, since
both have deficit four. Arc triangle loads supply the proof coordinate.

For n=7 the only remaining term is linear. As E_r=sum_j binom(r,j)C_j,
changing only that linear coefficient changes C_1=H and no other C_j.
This proves the homogeneous-port exchange law in SYMMETRY.md structurally.
It also explains why substituting the circular and Paley seven modules can
have identical multiple-run responses but different symmetry and prime rules.

At n=9 there is an additional cubic coefficient alpha:

    E_r=r^9+60r^7+768r^5+alpha r^3+beta r.

It records 2C7+4M35+8M333, where M35 is the number of disjoint triangle/
five-cycle pairs and M333 the number of unordered disjoint triangle triples.
Three explicit circulants give (H,alpha,beta)

    (3267,1994,444), (3357,2030,498), (3159,1916,414).

Their C2 values are 31110+6alpha and differ. Thus the next sidecar is the
deficit-six cycle collection, not only H. In full nine-vertex homogeneous
contexts, replacing B by B' changes the path count by

    delta H * b1 + 6 delta alpha * (b2+b3).              (5)

Indeed r^3=binom(r,1)+6binom(r,2)+6binom(r,3).
Formula(5) is an exact next-order repair of the failed one-coordinate
seven-vertex rule; it retains the actual marked quotient coefficients b_j.

More generally, for two regular modules of the same odd order n>=7,
their E_r difference has degree at most n-6. Therefore their six highest
nonempty-run responses C_(n-5),...,C_n agree. Writing r^m as
sum_j j! S(m,j) binom(r,j), these shared values are explicitly

    C_j=j! [S(n,j)+a S(n-2,j)+b S(n-4,j)], j>n-6,

where a,b are the two coefficients in(4) and S(m,j) are Stirling numbers
of the second kind. At most (n-5)/2 variable monomial coefficients remain.
This is a rigorous all-order compression for this particular homogeneous
path-cover observer; it does not identify the underlying tournaments or
settle which coefficient vectors are realizable.

## Extremal boundary and attribution

The mean arc load is (n+1)/4. Applying Cauchy--Schwarz in (3) gives

    C5 <= n(n^2-1)(n-2)(n-3)/160,

with equality precisely when every arc has equal load, equivalently the
tournament is doubly regular. This bound is **not new**: Savchenko's
2024 paper proves it even for all odd-order tournaments. Its introduction
also locates earlier regular-tournament results. Our derivation recovers
the regular case as a byproduct. For n=1 mod4, integer loads have half-integer
mean, giving the further elementary correction -n(n-1)/8 to this bound.
The equality condition for this correction is that all arc loads are the
two nearest integers; no existence at every order is asserted.

Sources checked in this session:
- S. V. Savchenko, *On 5-cycles and strong 5-subtournaments in a tournament
  of odd order n*, abstract/introduction and stated bound:
  https://arxiv.org/html/2403.19555v1
- N. Komarov and J. Mackey, *On the number of 5-cycles in a tournament*,
  Theorem2.1 (a more general edge-score identity):
  https://arxiv.org/html/1410.6828v3

The fixed-balance law does not extend to every shared score sequence.
`regular_balance.json` records the first hostile in the complete small-class
scan. Literal cycle enumeration, pair enumeration, matrix multiplication,
and ordered path-cover convolution check different sides of the proof.
The all-order identity follows from the proof, not from these finite checks.
