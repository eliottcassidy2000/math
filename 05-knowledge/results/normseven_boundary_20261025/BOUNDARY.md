# Boundary normalization with marked characters and arithmetic sheets

2026-10-25. **PROVED below; FINITE-EXACT controls in check.py and an
independent monoid-generation algorithm in audit_monoid.py.** Internal
adversarial audit, not an independent human/agent audit. No claim of a
whole-surface F1 descent, rational hourglass point, or Collatz convergence.

## 1. Inheritance, and the actual additional result

The closest proved mechanism is [KUMMER_AND_F1, K1–K4](../normseven_f1_20261025/KUMMER_AND_F1.md): complementary characters recover the full
twisted root algebra on the open torus. The chart cocycle and full
eigensheaf multiplication were already proved in [CHARACTERS, Theorem 1](../six_seven_gale_cube_20261025/CHARACTERS.md). We reuse them, not claim
them as new. The canonical hostile is the collapse of coordinate axes
under a monomial map. The corrected near miss is that a unimodular Laurent
change need not be regular where roots vanish. The least-used sidecar is
the *entire eight-character quotient algebra*, together with its integral
coefficient ideals, rather than three chosen character generators.

The earlier boundary gate [THM-2685, equivariant Kummer boundary parity](../../../01-canon/theorems/THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate.md)
detects ramification through divisor parity and retains unit squareclasses.
[THM-3397, torsor killing versus effective boundary valuations](../../../01-canon/theorems/THM-3397-torsor-killing-versus-effective-boundary-valuations.md)
prevents interpreting a finite cover as erasing poles. Our construction
instead adjoins integral functions omitted by the product algebra.

The additional result is an explicit conductor, its compatible seven-chart
atlas, and the residual arithmetic obstruction to lifting boundary points.

Concept board: character algebra / conductor / residual squareclasses /
monoid saturation / native source ownership / clocked affine packets.
Anchor: compatible norm-seven boundary charts. Niche: source-owned Collatz
histograms. Wildcard: zero-length tape blocks and exact native acceleration.

## 2. A general complementary-root algebra theorem

Let R be a UFD of characteristic zero, b_1,...,b_n nonzero, pairwise
nonassociate prime elements times units, and

    A = R[u_1,...,u_n]/(u_i^2-b_i).

Write V=F_2^n, u_z=product u_i^(z_i), and b^v=product b_i^(v_i).
The monic equations give the free R-basis {u_z:z in V}. For complementary
binary subspaces C,D with V=C direct_sum D, set

    A_C = direct_sum_(p in C) R u_p,
    A_D = direct_sum_(q in D) R u_q.

These are subalgebras because

    u_p u_q = b^(p AND q) u_(p XOR q).                 (1)

In characteristic zero they are the fixed algebras of the annihilator
subgroups of the sign deck. This identification is not asserted in
characteristic two. For every z, let z=p(z) XOR q(z) be its unique
decomposition and put e(z)=p(z) AND q(z). Let J be the set of coordinates
whose coordinate functional is nonzero on *both* C and D.

**B1 (image, conductor, and normalization).** Multiplication injects
A_C tensor_R A_D into A, with image

    B = direct_sum_z R b^e(z) u_z.                    (2)

Its conductor in A is

    (B:A) = (product_(i in J) u_i) A.                 (3)

If A is normal, A is the normalization of B in their common fraction
field. The ideal in (3) is principal on the normalized side; it need not
be principal as a B-ideal.

**Proof.** The tensor basis u_p tensor u_q maps to a nonzero scalar times
u_(p XOR q). The character labels are distinct, proving injectivity and
(2). The independent prime valuations show that the b_i have independent
squareclasses, so A is a domain and has fraction field of degree 2^n over
Frac(R). Inverting all b_i makes (2) an equality, hence B and A have the
same fraction field. Every u_i is integral over B by u_i^2=b_i. Normality
of A therefore gives the normalization statement.

For the conductor, a homogeneous element f_z u_z multiplies u_y into B
exactly when

    b^e(z XOR y) divides f_z b^(z AND y), for all y.   (4)

For an arbitrary sum of homogeneous elements, multiplication by each
fixed u_y still has distinct character labels, so (4) is necessary and
sufficient separately in every component. No cancellation of coefficients
from different characters is permitted. Since R is a UFD, the least
exponent required on b_i in f_z is

    c_i(z)=max_y max(0,e_i(z XOR y)-z_i y_i).          (5)

If i is not in J, e_i is identically zero. If i is in J, p_i and q_i are
independent nonzero linear functionals on V=C direct_sum D. If z_i=0,
choose z XOR y with p_i=q_i=1, obtaining c_i(z)=1. If z_i=1 and y_i=0,
the two components of z XOR y have opposite i-th bits, so e_i=0. If
z_i=y_i=1, subtraction in (5) makes the result nonpositive. Thus

    c_i(z) = 1-z_i for i in J, and 0 otherwise.        (6)

Multiplying u_(J XOR z) by product_(i in J)u_i has exactly coefficient
b^c(z) and character z. This proves (3), including minimality: lowering
any required exponent violates one of (4). QED.

The UFD and distinct-prime assumptions do real work. If b_1=b_2=t, the
two Kummer characters are no longer independent and the fraction-field
argument changes. If a coefficient is a unit, its divisibility condition
disappears. We do not transport (3) blindly through such specializations.

## 3. Apply the theorem to the two norm-seven quotients

Use seven forms in the inherited branch order

    (E, E+x, E-x, E+z, E-z, E+x+z, E-x-z).

Retain arbitrary constants d_i in k*, char(k)=0, and projective roots
s_i^2=d_i L_i. On chart L_j!=0, put

    b_i^(j)=d_i L_i/(d_j L_j),   u_i^(j)=s_i/s_j, i!=j.

The six b_i are distinct affine-linear primes times constants. The normal
root-cover result in [TORIFICATION_OBSTRUCTION, T1](../session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md)
applies; the three triple intersections give A1 points, not a failure of
normality. Constant twists preserve geometric normality.

The displayed complementary Fano markings are

    U=(1,2,4,5,6,3,7), D=(1,2,4,6,3,7,5).

Their character codes consist of the seven weight-four words obtained by
nonzero linear functionals on these labels, together with zero. Puncture
coordinate j to obtain complementary three-dimensional subspaces of
F_2^6. Every coordinate functional is nonzero on both codes, so J is all
six coordinates. The same is true for every complementary labelled Fano
pair, including the second inherited norm-seven pair.

**B2 (the seven-chart normalization).** If Y is the fibre product of the
two quotient surfaces over P^2, the original root surface X is its
normalization. On chart j the image algebra is (2) and its conductor is

    P_j A_j,    P_j=product_(i!=j)(s_i/s_j).

This concerns the normal finite root surfaces, not a presumed finite map
between their smooth resolutions. The tensor algebra is integral by (2).

The full quotient algebra has eight free character components. In
general three selected character functions do not generate it over the
base: multiplying two of those functions introduces the branch factors in
(1), and dividing them back out may be illegal. For example, in chart 0,
the U-code characters from masks 1 and 2 have overlapping root support;
their XOR character cannot be obtained from that product without its
branch correction. Relative to an independent three-character basis,
every monomial for this XOR uses those two generators an odd number of
times, so adding the third generator or even powers cannot remove the
unwanted factor. This is why the full eight-component algebra is kept.

## 4. Gluing preserves the old ports, including twists

Use *even subsets* S of {0,...,6} as the global character labels. Inherited
CHARACTERS Theorem 1, with the twists retained, gives

    e_(S,j)=product_(i in S)s_i / s_j^|S|,
    e_(S,k)=(b_k^(j))^(-|S|/2) e_(S,j).               (7)

The transition is a base unit on an overlap, since |S| is even. It changes
neither S nor any original port expressed in these marked functions. The
nonzero Fano characters have |S|=4, so their exponent is always -2.
The transitions satisfy the triple-overlap cocycle by cancellation of
ratios d_i L_i. There is no unrecorded square-root choice.

For the conductor,

    P_k=(u_k^(j))^(-7)P_j.                            (8)

Indeed P_j=(product_i s_i)/s_j^7. Consequently the conductor ideal on X is
the invertible sheaf O_X(-7H), generated locally by the product of the
ramifying roots, with H=O_X(1). It is the ideal of the reduced ramification
divisor: distinct root divisors share no codimension-one component and
each is generically reduced. Normality and the Cohen–Macaulay property
exclude hidden embedded components of this Cartier divisor. On this
seven-root canonical model K_X=H, so it is also O_X(-7K_X). The last
identification is special to seven roots; it is not an all-n formula.

The seven here counts branch coordinates. No additional identification
between a conductor degree and the Eisenstein norm of 2-omega is inferred
merely from their equal numerical values.

**B3 (determinant valuation).** The inclusion B_j -> A_j has determinant
a unit times product_(i!=j)b_i^16. Each of the independent pairs of bits
(p_i,q_i) equals (1,1) for 16 of the 64 characters, proving the assertion
directly from (2). In particular, at the generic point of each branch
line the R-module A/B has length 16. The trace discriminants differ by
the square of this determinant. In the squarefree root basis,

    disc(A/R)=64^64 product_i b_i^32,
    disc(B/R)=64^64 product_i b_i^64,                  (9)

up to the harmless sign from ordering bases. Trace is zero on nontrivial
characters and trace(u_z^2)=64 b^z; this proves (9) without a numerical
discriminant fit. These are base-divisor valuations, not a claim of 16
missing geometric sheets: the latter is given next.

## 5. Exactly which sheets and arithmetic data are absent?

Suppose r branch lines vanish at a geometric base point and their inertia
I has rank r. Assume I injects into each quotient deck. For the chosen
hourglass pair this holds at every branch stratum, including triples
012,034,056. Then the reduced geometric fibres have cardinalities

    X: 2^(6-r),  each quotient: 2^(3-r),
    Y: 2^(6-2r),  X over each point of Y: 2^r.         (10)

**Proof.** A fibre of the sign cover is the transitive G-set G/I. Its
quotient by H is G/(H+I). Since I intersects either kernel trivially and
the kernels are complementary, the joint image has 2^(6-2r) points and
the map from G/I has constant fibre 2^r. Alternatively evaluate the
characters on the 2^(6-r) choices of nonvanishing root signs. QED.

Equation (10) counts *reduced geometric points*, not scheme lengths. Both
X and Y remain finite flat of rank 64 over the base; their special fibres
have nilpotents. Thus r=1,2,3 require respectively one, two, three extra
binary sheet choices. At a triple, Y has one reduced point and X has
eight A1 points above it.

Over a nonclosed field the choices form a possibly nontrivial elementary
quadratic torsor. This is arithmetic information, not just sign notation.
For a rational base point, let V_0 be the (6-r)-dimensional character space
of its nonzero root coordinates. Characters from either quotient that
survive at this point form subspaces C_0,D_0 of dimension 3-r. Their sum
W_0 is direct, hence

    dim V_0/W_0 = r.                                 (11)

A rational point of each quotient trivializes the prescribed squareclasses
on W_0. The remaining root squareclass map V_0 -> k*/k*2 therefore factors
through V_0/W_0. Full lifting requires it to vanish there too. If all
classes vanish, the compatible character values extend and give 2^r
rational lifts; if any is nonzero, there is no rational lift over that
chosen product point. Character extension over F_2 proves sufficiency.

The r missing classes can be prescribed independently by changing the
unit twists d_i at the fixed point. Thus (11) is sharp, not an upper bound
which disappears after remembering that the charts come from lines.

**Explicit hostile, r=1.** Take (E,x,z)=(1,-1,3), chart 0, and

    d=(1,1,1/2,1/2,-1,2/3,-1).

The seven normalized squared roots are (1,0,1,2,2,2,1). Every surviving
nontrivial U or D character has value 2; vanished characters have value 0.
The full multiplication table (1) is satisfied over Q in both quotients,
but X has no rational point here because s_3/s_0 squares to 2. All 16
rational reduced product points have quadratic, rather than rational,
lifts. No change of Laurent generators removes this obstruction.

**Explicit hostile, r=3.** Take (E,x,z)=(0,0,1), chart 3, and

    d=(1,1,1,1,-2,3,-5).

The zero roots are 012; the three remaining nontrivial unit radicands
are 2,3,5. Each quotient has one rational reduced boundary point, while
the normalization above it has residue field Q(sqrt(2),sqrt(3),sqrt(5)),
of degree eight by the three independent prime valuations. Over Qbar it
splits into the eight points in (10).

## 6. A literal monoid/F1 model, and its exact stopping boundary

Before imposing the additive branch equations, replace the six radicands
by independent parameters t_i. The universal normalized algebra is
Z[u_1,...,u_n], with t_i=u_i^2. The product algebra is the monoid algebra
Z[M], where

    M = generated_by {2e_i} union C union D in N^n,
      = disjoint_union_z (z+2e(z)+2N^n).              (12)

Here binary vectors C,D are interpreted as nonnegative integer exponent
vectors. The equality follows either from (1) or directly by reducing
products within C and D and collecting even exponents.

**B4 (monoid saturation).** M has group completion Z^n, saturation N^n,
and conductor into N^n equal to 1_J+N^n.

**Proof.** It contains 2Z^n in its group completion and maps onto F_2^n,
so the group completion is Z^n. Every v in N^n has 2v in M; any vector
with a negative coordinate has no positive multiple in M. This proves
the saturation. Apply the coefficientwise calculation (4)–(6) to the
independent variables t_i, or use (12) to test translations by each parity
vector. The conductor is precisely the exponent set of
(product_(i in J)u_i)Z[u]. QED.

This gives an actual affine monoid scheme and its normalization in the
sense of [Deitmar, Schemes over F1, Sections 1–2](https://arxiv.org/html/math/0404185), with base extension Z[M]. It is a specific
universal multiplicative object, not an asserted F1 model of the full
surface. The latter is obtained by t_i=d_iL_i/(d_jL_j), which uses additive
linear relations. Normalization does not commute with arbitrary base
change. In this case B1–B2 prove the claimed normalization *after* the
specialization, using the distinct branch primes and actual normality.

This is a precise reroute: keep the monoid normalization as the reusable
multiplicative template and record the additive specialization separately.
The inherited geometric-genus obstruction to finite torification of the
whole seven-root surface remains in force.

## 7. The representation graph that retains the boundary

The faithful vertices are the 64 characters z. Multiplication by root i is
a labelled operator

    T_i u_z = b_i^(z_i) u_(z XOR e_i),
    T_i^2=b_i I,    T_i T_j=T_j T_i.                  (13)

This is the regular module of the root algebra, with labelled commuting
squares. At b_i=0 its i-th operator becomes square-zero. It cannot be
treated as the invertible sign-toggle available on the open torus.

The image B is a diagonal lattice in this module. Its conductor is the
largest A-stable submodule contained in B. Formula (5) is a finite
max-plus feasibility calculation on its multiplication rows, and the
closed formula (6) solves it at all sizes under B1's hypotheses.

An honest tournament probe clarifies what would be lost. Given a base
valuation v, multiplying u_x by the character joining x to z costs

    cost(x->z)=sum_i v(b_i) x_i(1-z_i).

The difference between the two directions is height(x)-height(z), where
height(x)=sum_i v(b_i)x_i. Keeping only the cheaper direction gives a
transitive preorder with ties. It loses the overlap factors and the port
labels that determine (2). For example e_1 and e_2 tie when their branch
valuations agree but act differently on every marked coordinate port.
We therefore retain both directed costs, all ties, and the commuting
squares. This representation-theoretic graph is useful precisely because
it does not force an uninformative tournament.

## 8. Boundaries and next questions

The construction closes the stated algebraic chart-compatibility problem;
it does not establish rational points or six-colour forcing. The exact
remaining arithmetic problem is the vanishing of the residual class map
in (11), with the original twists fixed. Local valuation parity alone is
not enough: the r=1 example is entirely a residual unit-squareclass issue.

The cyclic-root extension was pursued within this session and is proved in
[CYCLIC, C1](CYCLIC.md), including the composite-modulus correction and
the classical cusp. The next extension is to a normal non-UFD base, using
fractional ideal intersections and class-group data from THM-2685.
Repeated radicands, dependent branch inertia and coefficient characteristic
p remain the first hostile cases; no characteristic-p normalization claim
is inferred from a characteristic-zero character calculation.
