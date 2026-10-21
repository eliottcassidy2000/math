# Independent audit: explicit R5 raw circulation series

Session 23b, 2026-10-21. Result: **PASS** for the submitted addendum and all its explicit formulas and asymptotic claims. No correction to its mathematical conclusions is required. The proof below supplies independent coverage, nonoverlap, boundary, and polynomial certificates rather than trusting the producer's pulling-triangulation implementation or fitting coefficients.

Audited source: `session23b-flows-extra/ADDENDUM.md`, SHA-256 `ab0465a1d9fc9b1ed371205b5f12ab45aa1e95d24af79cba41671442cf5639a8`. The source `check.py`, both submitted receipts, and every entry of its manifest were also verified. The independently reviewed primary flow packet is an input, not re-audited here. The anti-fixed formula is accepted from the separate five-case proof in `session23b-pentagon/PROOF.md`; the present audit also checks it independently through weight 100.

## 1. Lattice and ray completeness

On the actual ten-arc tournament write x_i for i -> i+1 and y_i for i -> i+2. Balance gives

    x_i+y_i = x_(i-1)+y_(i-2).

Therefore the five integers x_i+y_i+y_(i-1) coincide. Their common value t gives the mutually inverse integral maps

    (x,y) -> (t,y),
    x_i = t-y_i-y_(i-1).

Consequently the circulation lattice is exactly Z^6 in coordinates (t,y), not a proper sublattice. Nonnegativity becomes y_i >= 0 and t >= y_i+y_(i-1); this is precisely the displayed cone after cyclically reindexing the five adjacent-pair inequalities. The literal arc total is 5t-sum(y_i). The checker independently verifies balance and nonnegativity for every displayed ray.

At t=0 the cone has only the origin. The slice t=1 is bounded because 0 <= y_i <= 1. An extreme point of this full-dimensional five-dimensional polytope has five linearly independent active inequalities. An independent Fraction-only Gaussian-elimination enumeration of all C(10,5)=252 active sets finds 122 nonsingular systems and exactly the twelve asserted vertices. Thus the displayed rays generate the entire cone. Every integer ray vector has coordinate gcd 1. In particular the all-half slice vertex produces primitive ray (2;1,1,1,1,1), not a spurious doubled lattice generator.

## 2. Independent full-fan certificate

The audit takes only the producer's 25 proposed six-ray lists as candidate data; it does not import or reproduce its recursive pulling routine.

For each ray matrix M, independent exact Gaussian elimination gives det(M)=+1 or -1. Each cone is full dimensional and unimodular in the saturated circulation lattice.

For each of the 300 unordered pairs of cones C_a, C_b, the audit constructs and exactly verifies a rational linear functional h such that:

- h(r)=0 on every shared ray;
- h(r)>0 on every nonshared ray of C_a;
- h(r)<0 on every nonshared ray of C_b.

Every certificate is stored in `normal.json`, `cones.pairwise_separators`. Because coefficients in the conical combinations are nonnegative, an element of both cones must have h=0 and can use only shared rays. Conversely every conical combination of shared rays is in both. Thus

    C_a intersect C_b = cone(shared rays).

This proves all pairwise intersections are the intended common faces and all full-dimensional interiors are disjoint. The separators are found by finite integer search in sums of dual-basis rows, and verified by exact rational arithmetic. No numerical optimization is used.

The 150 simplex-facet incidences have exactly 100 distinct five-ray facets:

- 50 facets have two incident cones, whose remaining rays lie strictly on opposite sides of the facet hyperplane;
- 50 facets have one incident cone, and each lies in at least one original supporting hyperplane of the pentagon cone.

The complete incidence and supporting-inequality data are stored in `normal.json`, `cones.facets`. There is no singly incident facet in the interior of the original cone, and no twice incident exterior facet.

These finite checks imply coverage. Let U be the finite union of the candidate cones, a closed subset of the original cone C. It contains a full-dimensional open set. If int(C) contained a point outside U, openness of the complement would permit a generic segment from a point inside a candidate cone to a point outside U, wholly in int(C), that avoids the finitely many faces of codimension at least two. Its first exit from U would cross the relative interior of a simplex facet. Such a facet cannot be an original boundary facet. It must therefore have two incident cones on opposite sides, so their union contains a neighborhood of the crossing point. This contradicts first exit. Hence int(C) is contained in U, and closedness gives C=U.

This proves coverage and conformity independently of the pulling-recursion claim.

## 3. Half-open boundary assignment

Set w=sum_i 2^i r_i. Direct exact checks verify that all original inequalities are strictly positive on w and all coordinates of M^(-1)w are nonzero for every simplex matrix M.

For x in C, x+epsilon*w remains in int(C). For sufficiently small positive epsilon it avoids every simplex facet hyperplane: on a hyperplane through x, its sign is fixed by the nonzero functional value at w; off such a hyperplane, its sign is unchanged when epsilon is sufficiently small. Thus it lies in exactly one simplex interior, by Section 2.

Write lambda=M^(-1)x and alpha=M^(-1)w. Membership in the uniquely assigned half-open simplex is exactly

    lambda_j >= 0 when alpha_j > 0,
    lambda_j > 0 when alpha_j < 0.

This proves the producer's sign convention, including points on the original boundary, the origin, rays, shared facets, and higher-codimension intersections. No lattice points are omitted or counted twice. Since M is unimodular, integral x has integral lambda. An open ray therefore starts at coefficient 1; a closed ray starts at 0. Its contribution is precisely

    q^(sum of open-ray weights) / product_j(1-q^(ray weight_j)).

As a hostile bounded test, the independent checker explicitly assigns every one of the 8,921 integral points of total weight at most 25 to these half-open cones. Each has exactly one owner. Of these, 8,896 lie on at least one simplex boundary; thus this control genuinely exercises boundary cases rather than merely generic interior points.

## 4. Numerator certificate and exact raw series

Independent recomputation of all dual coordinates, open rays, shifts, and determinants agrees with every submitted cell record. There are exactly two weight types. The fifteen cones of type (3,3,4,4,5,5) have shift polynomial

    A=1+q^3+3q^4+2q^6+4q^7+q^8+3q^10.

The ten cones of type (3,3,3,4,5,5) have shift polynomial

    B=2q^3+2q^6+2q^7+q^9+2q^10+q^13.

Dense exact polynomial arithmetic verifies

    A(1-q^3)+B(1-q^4)=P(1-q^5),
    P=1+2q^3+3q^4+q^5+3q^6+q^7+3q^8+2q^9+q^12.

Sections 1-3 and this polynomial identity prove, to all orders,

    R(q)=P/((1-q^3)^3(1-q^4)^2(1-q^5)).

The proof does not infer a numerator, denominator, or recurrence from finite coefficients. A separate direct enumeration of the original pentagon inequalities agrees with this formula through total weight 45, extending the producer's weight-40 check. The bound is complete because summing the five inequalities gives 2 sum(y_i)<=5t and hence weight>=5t/2.

Ordinary Burnside gives f=(R+4/(1-q^5)^2)/5, as a nontrivial rotation fixes exactly the two constant arc-orbit coordinates. The supplied independent anti-fixed formula is

    k(q)=(1+q^4+q^5+q^6+q^7+q^8+q^12)
         /((1-q^3)(1-q^5)(1-q^6)(1-q^8)).

The present checker independently enumerates its four-dimensional lattice domain through total weight 100. All coefficients agree. Ordinary and converse Burnside integrality/nonnegativity also pass through weight 45.

## 5. Denominators, recurrences, quasiperiods, and growth

Exact rational Euclidean polynomial arithmetic verifies coprimality of the displayed R and k numerators and denominators. It also independently reduces f and H=(f+k)/2. Their complete cyclotomic pole profiles are:

    series   Phi_1 Phi_2 Phi_3 Phi_4 Phi_5 Phi_6 Phi_8
    R           6     2     3     2     1     0     0
    f           6     2     3     2     2     0     0
    k           4     2     2     1     1     1     1
    H           6     2     3     2     2     1     1

Every stated pole is genuine; the reduced numerator and denominator are coprime. Consequently the minimal orders of eventual homogeneous constant-coefficient recurrences over Q are respectively 22, 26, 22, and 32. In particular, the order-22 raw denominator must not be assigned to the ordinary orbit series f. The submitted addendum makes no such erroneous identification.

The minimal quasiperiods are respectively 60, 60, 120, and 120; the degrees are respectively 5, 5, 3, and 5. These are exact quasipolynomials for nonnegative indices because all four reduced rational functions are proper. Partial fractions give an upper bound equal to the stated period, while the genuine primitive-root poles force divisibility by its full least common multiple. The expanded recurrence-denominator coefficient arrays are recorded in `normal.json`, `rational_consequences.recurrence_denominators`. Their convolution identities are also checked against the independent bounded coefficients.

At q=1 the raw pole has order six and P(1)=17. The fixed pole has order four and its numerator evaluates to 7. Thus

    R_m ~ 17 m^5/259200,
    f_m ~ 17 m^5/1296000,
    k_m ~ 7 m^3/4320,
    k_m/f_m ~ (2100/17) m^(-2).

No other pole has the top order: for R and f the largest nonidentity order is three, and for k it is two. Therefore these are full-sequence asymptotics, not just selected residue-class statements. The ratio in the addendum is correct.

This audit does not replace the separately assigned general reciprocity proof. The addendum's algebraic reciprocal-shift checks are consistent with its reduced numerators and denominators.

## 6. Reproducibility and scope

`check.py` uses the Python standard library only, exact Fraction arithmetic, and explicit `need` checks rather than optimization-sensitive assertions. Its normal and `python -O` receipts are byte-identical. Both modes of the producer's checker were also independently rerun; their outputs match the corresponding submitted receipts byte-for-byte. All source-manifest byte counts and SHA-256 values match.

The audit directory contains:

- `check.py`: independent certificate and controls;
- `normal.json`, `optimized.json`: the complete exact independent receipts;
- `producer-normal-rerun.json`, `producer-optimized-rerun.json`: verified producer reruns;
- `manifest.json`: hashes of this audit packet;
- this `REVIEW.md`.

The conclusion is limited to the explicit R5 cone, its ordinary and anti-fixed quotients, the stated rational/algebraic consequences, and the use of the pre-reviewed primary packet. No general cyclic-tournament formula or algebraic independence of primitive cycles is inferred.
