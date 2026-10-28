# Recovering information from incidence and transport

2026-10-28. **PROVED elementary statements and applications; CITED classical
frameworks; FINITE-EXACT enumerations.** Classical Hamming codes, trace forms,
signed graphs and parity coding are not claimed as new. Local names P1–P6
below are packet labels, not reserved canon IDs.

## P1 The exact cost of a bit when linear information is already known

Let x range over ALL of F2^n. The observer already knows Ax, can ask for
individual coordinates of x, and must determine h x, where h is a binary
row. The minimum worst-case number of coordinate queries, even allowing
adaptivity, is

    min_{a in row(A)} weight(h+a).                         (1)

**Proof.** A minimizing a gives an algorithm: add the already known linear
combination a x to the queried parity (h+a)x. Conversely, at any leaf of
an exact decision tree, let S be the queried coordinates. If h is not in
span(row(A), {e_i:i in S}), there is a vector delta annihilated by A and
every e_i in S but with h delta=1. Both x and x+delta follow the same
adaptive transcript and need opposite answers. Thus h=a+s with s supported
in S, and |S| is at least (1). This also proves the lower bound for every
nonempty affine fibre Ax=y. QED.

This is a standard linear-algebra/query-complexity mechanism. Its concrete
use here is to distinguish a jointly held sign from a separately stored
sign distribution. Information is not created: a previously paid-for parity
replaces one coordinate query.

Let H_r have as columns every nonzero vector v of F2^r, once, so
n=2^r-1, with r>=2. Write P=sum_v x_v. For each nonzero c,

    c dot (H_r x) = P + sum_{v != 0, c dot v=0} x_v.       (2)

The row cH_r has weight 2^(r-1); its complement is the nonzero part of a
hyperplane, of size 2^(r-1)-1. Hence one syndrome check costs exactly
2^(r-1) coordinate queries without P, and exactly one fewer with P.

For deciding the complete predicate H_r x=0, worst-case query counts are
exactly n without P and n-1 with P. The upper bounds read all coordinates,
or all but one and infer the last from P. For the lower bound follow an
accepting input. Without P, any unread coordinate can be toggled to change
the syndrome. With P, any two unread distinct coordinates can be toggled
together, preserving P but changing the syndrome by two distinct columns.
This also allows adaptive algorithms; no fixed-query restriction is hidden.

### The original seven-prime hostile now has a joint separator

Use the inherited bank {2,5,7,19,29,59,89} with square-class coordinates
(2,3,5). A selection x gives squarefree n_x, so

    lambda(n_x)=mu(n_x)=(-1)^P,
    Q(sigma(n_x))=H_3 x.

For c=(1,0,0), the hyperplane's nonzero columns are vectors 2,4,6,
labelled by source primes 2,19,59. Formula (2) becomes

    v_2(sigma(n_x)) mod 2 = P + x_2 + x_19 + x_59.         (3)

Uniform laws on C_0={Hx=0} and C_1={Hx=(1,0,0)} have identical
three-coordinate marginals and identical separate lambda/mu laws, as
proved in the [earlier packet A3a](../incidence_attractors_20261028/THEORY.md).
But their JOINT law of (P,x_2,x_19,x_59) is disjoint under the parity
test (3). P jointly with any at-most-two coordinates still has identical
laws. The checker exhausts all 29 such coordinate sets; an independent
adaptive decision-tree optimizer obtains the sharp costs 4/3 and 7/6.

More generally, uniform H_r-cosets agree on all coordinate marginals of
order <=2^(r-1)-1. One proof: after removing that many columns, the remaining
at least 2^(r-1) distinct nonzero vectors cannot fit in a proper hyperplane,
so the projection fibres have equal cardinality. With P included jointly,
the uniform laws agree through order 2^(r-1)-2. Indeed the nonzero rows of
H_r have weight 2^(r-1), and their complements have weight 2^(r-1)-1;
no lower-support dual relation can couple P and the coordinates. Equation
(2) reaches the first possible distinguishing order. This is the precise
Vitali/alignment lesson in this example: the missing resource was joint
correlation, not another independently estimated marginal.

## P2 Fifteen explicit primes realize a larger arithmetic geometry

On square-class row basis (2,3,5,7), choose the following prime for each
nonzero vector. Vectors are binary integers, least significant bit first.
Each row has p+1=d t^2. The displayed p is the LEAST prime of that fixed
square class d, by checking every smaller positive t; no optimality over
all choices of row basis is asserted.

| v | d | p | t |
|---:|---:|---:|---:|
|1|2|7|2|
|2|3|2|1|
|3|6|5|1|
|4|5|19|2|
|5|10|89|3|
|6|15|59|2|
|7|30|29|1|
|8|7|251|6|
|9|14|13|1|
|10|21|83|2|
|11|42|41|1|
|12|35|139|2|
|13|70|5669|9|
|14|105|419|2|
|15|210|839|2|

For all 32768 subset products, multiplicativity gives

    sigma(n_x) square iff H_4 x=0.

Thus exactly 2048 accepted products realize the Hamming [15,11,3] code.
Its weight enumerator, independently checked by character averaging, is

    1 +35z^3 +105z^4 +168z^5 +280z^6 +435z^7
      +435z^8 +280z^9 +168z^10 +105z^11 +35z^12 +z^15.

The 35 weight-three words are precisely the lines {a,b,a+b} in PG(3,2).
Seven lines pass through each point: pair the other fourteen points by
b <-> a+b. Weight-three words do NOT exhaust the matroid's circuits;
the earlier Fano-circuit correction remains in force.

The arithmetic construction is unconditional for this finite bank. Extending
it to every rank with a prescribed square-class basis would require primes
among several quadratic sequences d t^2-1. A successful finite search is
not an infinitude theorem, and the zero class is exceptional: p+1 square
forces p=3. The graph obstruction persists because this matroid contains
the original nongraphic Fano restriction.

### A literal bridge between seven, fifteen and five

Choose a nondegenerate alternating bilinear form omega on F2^4. Among the
35 lines, exactly 15 are totally isotropic; each point lies on three of
them. These form the classical symplectic generalized quadrangle W(3,2),
often called the doily. There are 28 possible such forms on the marked
space, so arithmetic zero-sums alone do not specify the chosen doily.

The 35-line geometry has 56 spreads (partitions of the fifteen points into
five disjoint lines). For a fixed omega, exactly six are isotropic spreads.
These finite counts are exhausted in check.py. For the standard form
omega(a,b)=a0*b2+a2*b0+a1*b3+a3*b1, one resulting partition into five
square-total prime triples is

    {7,2,5}, {19,251,139}, {89,83,839},
    {59,41,5669}, {29,13,419}.                            (4)

Their divisor-sum square roots are respectively 12,840,2520,3780,420.
All source primes are distinct. Hence their products are coprime arithmetic
blocks as well as disjoint projective lines.

The link with the earlier six-label matching work is explicit. Let
V={even subsets of six labels}/complement, with intersection-parity pairing.
Each of the fifteen duads is a nonzero vector. Each perfect matching gives
three disjoint duads whose classes sum to zero and span a Lagrangian plane.
The fifteen matchings give exactly the isotropic lines. A one-factorization
of K6 consists of five matchings partitioning its fifteen edges, hence is
a spread. The independent audit enumerates the six factorizations and
checks the map using basis classes {i,5}, i=0,1,2,3. This basis has a
different matrix for omega from (4); a symplectic change of basis is needed
before identifying individual labels. The abstract incidence is preserved.

The same marked space describes two-qubit Pauli operators modulo scalar
phase: (a,b,c,d) maps to X^a Z^c tensor X^b Z^d. Commutation is exactly
omega=0. The prime bank thus realizes their 35 projective dependence
triples, while only 15 are commuting triples after omega is chosen.
Integer multiplication does not encode Pauli commutator signs or scalar
phases. This is an incidence bridge, not an algebra isomorphism. The Pauli
geometry is classical: [Planat–Saniga](https://arxiv.org/abs/quant-ph/0701211).

## P3 Subfield incidence and Frobenius reconstruct a symplectic form

Here is an all-size way to recover a missing choice from two structures
already held. Put q=2^m, m>=1, E=F_(q^2), and K=F_q. Regard E as a
2m-dimensional binary vector space, retaining its actual field structure.
Let S be its spread of one-dimensional K-subspaces xK.

**Claim.** There is exactly one nondegenerate alternating F2-bilinear form
on E that makes every member of S totally isotropic and is invariant under
absolute Frobenius x -> x^2. It is

    B(x,y)=Tr_(E/F2)(x y^q).                              (5)

**Proof of existence.** Since x^(q+1) lies in K and the extension E/K has
degree two, its absolute trace is zero. Hence B(x,x)=0. The trace pairing
is nondegenerate and y -> y^q is bijective, so B is nondegenerate.
For a,b in K, B(ax,bx)=Tr(ab x^(q+1))=0, proving isotropy. Finally
B(x^2,y^2)=Tr((xy^q)^2)=B(x,y).

**Classification before Frobenius.** Let C be any alternating binary form
vanishing on each K-line. For a in K, expanding C(a(x+y),x+y)=0 gives
C(ax,y)=C(x,ay). Choose a K-basis e,f of E. Its two coordinate lines
are isotropic. If ell(t)=C(e,tf), then

    C(ae+bf,ce+df)=ell(ad+bc).

Conversely every binary linear functional ell:K->F2 gives such a form.
There are q of them, including zero. The q forms

    B_c(x,y)=Tr_(E/F2)(c x y^q), c in K,

are distinct and have this property, so exhaust the possibilities. Every
nonzero c gives a nondegenerate form by the same trace-pairing argument.

**Frobenius selects the coefficient.** Trace invariance under squaring gives
B_c(x^2,y^2)=B_(sqrt(c))(x,y), where the square root is taken in K.
Distinctness of the forms forces c=sqrt(c), hence c in F2. Nondegeneracy
excludes c=0, leaving exactly c=1. QED.

Trace constructions and subfield spreads are classical finite geometry;
for the PG(3,q) symplectic-spread framework see
[Ball–Zieve, Symplectic spreads and permutation polynomials](https://web.mat.upc.edu/simeon.michael.ball/symplecticspreadpps.pdf).
The elementary all-dimension uniqueness argument here is
included in full; no literature-priority claim is made. Its use here is
as a concrete reconstruction rule for the arithmetic carrier.

At E=F16 this has a particularly sharp interpretation. Frobenius alone
preserves exactly two of the 28 nondegenerate alternating forms. The
F4-spread alone admits three. Their intersection contains one. In the
polynomial basis for alpha^4+alpha+1, forms are encoded by six upper-triangle
bits in order (01,02,03,12,13,23): the Frobenius-compatible masks are
12 and60, and the unique intersection is60. Form12 has no Frobenius-fixed
spread; form60 has two. The subfield spread is one of the latter.
These are FINITE-EXACT controls, not a proposed 28/56 numerology.

### Phase survives the reconstruction

For F8 use alpha^3+alpha+1 and normal basis vectors (3,5,7). For F16 use
alpha^4+alpha+1 and normal basis (8,12,15,10). Squaring rotates normal
coordinates by one place. Nonzero Frobenius orbit lengths are respectively
1,3,3 and 1,4,2,4,4. A vector's orbit length is its degree over F2; its
stabilizer in the cyclic Frobenius group has size r/degree.

The symplectic reconstruction does not select an absolute normal-basis
phase. Rotation of an unmarked normal basis still gives another equally
valid basis. Nor does its transport on prime labels preserve p, p+1,
aliquot height or natural order. It preserves the indicated binary
incidence; a named polynomial/normal basis is part of the arithmetic map.

## P4 Signed incidence is an affine holonomy calculation

The [previous relation packet D](../unifying_relations_20261028/THEORY.md)
already proves integral edge-sum reconstruction and its Smith factor2.
Here it sits in a slightly larger reusable construction.

For each directed edge u->v retain the equation

    x_v = s_e x_u + w_e,    s_e in {+1,-1}.               (6)

Difference labels use +1; square-sum labels use -1 and w_e=k_e^2.
On a spanning tree, propagation gives x_v=rho_v z+t_v, where z is the
root label and rho_v is the sign product on the root path. Each chord gives

    (rho_v-s_e rho_u) z = w_e-t_v+s_e t_u.                (7)

The coefficient is 0 or +/-2. A cycle with sign product +1 imposes only
compatibility. A cycle with sign product -1 pins the root, subject to an
even-numerator condition over Z and consistency with every other chord.
For all-negative edges, these are even and odd cycles respectively.
This explains a free root on a square-sum tree and a determined root on
an odd square-sum cycle. It also explains why reduction modulo2 loses the
distinction: +1 and -1 become the same coefficient.

Examples: the positive distinct triangle (5,20,44) has edge sums25,64,49;
its root is (25-64+49)/2=5. The square edge data (4,4,1) instead force
(1/2,7/2,1/2), passing rational compatibility but failing integer realization.
The even square-sum cycle (1,8,136,440) has sums9,144,576,441 and can be
shifted to (2,7,137,439), preserving every edge sum. Positivity and label
distinctness remain separate restrictions on the free root.

In the general vector version x_v=A_e x_u+b_e, with invertible matrices
over a field, tree propagation gives x_v=G_v z+t_v. The chords give

    (G_v-A_e G_u) z = A_e t_u+b_e-t_v.                   (8)

Stacking (8) gives an exact affine solution space. The root is recovered
precisely when it is nonempty and the stacked coefficient matrix has full
column rank. Over integers one must additionally enforce the corresponding
integral-image conditions, for example by Smith normal form. Noninvertible
edge relations require keeping their extra variables rather than silently
using inverse transport. Classical signed-graph rank is documented by
[Zaslavsky, Theorem3.6](https://people.math.binghamton.edu/zaslav/Tpapers/sggm.pdf).

LCA now has a specific role. It identifies shared prefixes in the chosen
tree, allowing root-to-vertex transports to be combined into the path
u->LCA(u,v)->v, with inverses and order retained. It speeds repeated chord
checks; it does not remove the checks or recover unrecorded branch signs.
For Berggren's noncommuting matrices, letter counts are insufficient.

For K3,3 the all-negative system is bipartite and retains one root
parameter; K5 and Petersen are nonbipartite and pin the rational root if
compatible. Their integer obstruction has a final Smith factor2. Petersen
spoke contraction still requires its one/two-run boundary information and
global connectivity for Hamiltonian questions. A rank calculation does not
turn the six disconnected Petersen 2-factors into Hamiltonian cycles.

## P5 Collatz shows when the next bit is not free

Use T_s(n)=n/2 for even n and (3n+s)/2 for odd n, s=+1 or-1.
For a legal parity word e_0,...,e_(k-1), put r=sum e_i and

    B_0=0,  B_(i+1)=3^(e_i) B_i + e_i 2^i.

Induction gives 2^k T_s^k(n)=3^r n+s B_k. The word is realized exactly on

    n = -s B_k (3^r)^(-1) mod 2^k.                      (9)

Necessity follows from the formula. For sufficiency, the usual parity
branch induction gives one residue for every word, with the two next
branches the two lifts modulo2^(k+1). This is the classical parity coding
of [Bernstein–Lagarias](https://websites.umich.edu/~lagarias/doc/bernstein.pdf),
not a new Collatz conjugacy.

The two starts n and n+2^k have the same first k parities but their
endpoints differ by 3^r, an odd number, so their next parities differ.
The next bit is not determined by the prefix. The starts n+2^k t have
unbounded ordinary height despite a fixed word. An independently proved
interval of width less than2^k can select at most one such start, but
proving the interval is work that a residue graph alone does not do.

A repeated word supplies only a candidate cycle
n=s B_k/(2^k-3^r); integrality, positivity, legality, and global orbit
questions remain. Reversing s changes the candidate equation, not a
proved count of all positive cycles. The proposed 'which edge starts first'
analogy therefore does not establish the single-cycle/triple-cycle claims.
The checker exhausts all words of lengths1..10 for both signs and verifies
their opposite next-bit lifts: 4092 words in total.

The recovered Antihydra/Mahler lane adds a different necessary coordinate:
minimum internal credit. After its marked scalar translation z=b+3 the
map is ceil(3z/2), but a legal machine block can fail before its positive
net gain is realized. Tree-prefix depth, affine transport, ordinary height
and minimum-prefix credit must not be collapsed into one scalar.

## P6 Degree bounds and geometric anchors remove other ambiguities

The [Session43 degree theorem](../session43_observers_phase_guards_20261023/FINITE_GUARDS_AND_DEGREE.md)
is inherited in full. For F=(P,Q) over C of degree <=D, assume exact
equivariance under (x,y)->(zeta x,zeta^-1 y), zeta of order m. Monomial
charges a-b lie in [-D,D]. If m>D+1 the required charge congruences
are equalities, giving full C*-equivariance. If det JF=c!=0 then
F=(xA(xy),yB(xy)), and (tAB)'=c forces A,B constant.

The sharp boundary is the shear (x+y^(m-1),y), degree m-1, determinant1.
Thus exact C5/C7/C15 symmetry gives linearity through degrees3/5/13.
Joint C5 and C7 gives effective modulus35 and degree33; C5 and C15 still
gives15, not75; C7 and C15 gives105 and degree103. This is a bounded-charge
reconstruction theorem, not a result for arbitrary planar Keller maps or
finite point samples. Planar JC remains outside the conclusion.

For geometric graphs, distances to two fixed anchors usually leave two
mirror positions. Subtracting squared-distance equations for three
noncollinear anchors gives two independent linear equations, hence a
unique point. With collinear anchors the reflection ambiguity survives.
This is another instance of recovering a choice through existing incidence
constraints. It requires the anchor coordinates and correct distances.

Hadwiger–Nelson colourability is a further predicate. The norm-seven
Eisenstein quotient T->K7 preserves its graph-cover incidence but loses
the metric lift and nonperiodic colourings; T itself is three-colourable.
K4 is already a cheap hostile to a metric-free lift: the Gram matrix of
three unit vectors from one vertex of a unit K4 has diagonal1 and
off-diagonal1/2, hence rank3 and no planar realization. The preserved
period shifts and Gram constraints cannot be inferred from a finite
colour quotient. Our new results give no plane chromatic bound.

The genuinely useful six-palette criterion from Session61 is this: for a
fixed ternary colouring q, keep the edges with equal ternary labels. A
binary lift exists iff that subgraph is bipartite. The graph depends on
q; testing one cycle space for every q is invalid. This is the same
transport/holonomy mechanism as (6), now over F2, with a different target.
