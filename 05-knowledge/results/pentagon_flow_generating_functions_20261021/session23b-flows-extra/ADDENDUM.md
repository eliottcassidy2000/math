# Explicit R5 flow series from a certified pentagon triangulation

Session 23b, 2026-10-21. This is a separate addendum to the frozen `session23b-flows/PROOF.md`. Its formulas are derived by exact finite geometry and elementary sums, not guessed from coefficients. The accompanying `check.py` is optimization-safe and uses exact SymPy rational arithmetic, never floating point. Independent review of the anti-fixed formula is being supplied separately; independent review of the raw triangulation should be read alongside this certificate before publication.

## 1. The three explicit weights

Let R(q) count all nonnegative integer circulations on the labelled cyclic tournament R5 with arcs i->i+1 and i->i+2 on Z/5, weighted by total edge use. Let f(q) count their Aut(R5)-orbits. Let k(q) count the ordinary flow-orbits fixed by converse. Then

    P(q)=1+2q^3+3q^4+q^5+3q^6+q^7+3q^8+2q^9+q^12,

    R(q)=P(q)/((1-q^3)^3(1-q^4)^2(1-q^5)),          (1)

    f(q)=(R(q)+4/(1-q^5)^2)/5.                     (2)

For the anti-fixed series,

    N(q)=1+q^4+q^5+q^6+q^7+q^8+q^12,

    k(q)=N(q)/((1-q^3)(1-q^5)(1-q^6)(1-q^8)).       (3)

The converse-merged flow-orbit series for this one underlying atom is (f+k)/2. These are different enumerators: R counts labelled flow assignments on a fixed tournament, f takes the actual order-five automorphism quotient, and k counts fixed ordinary orbits, not all flows on a self-converse underlying tournament.

The original degree-zero formula in the primary packet is consequently entirely explicit:

    A_ind(z,q)=1/(1-z^2-z^6/(1-q^3)^2),
    A_mir(z,q)=1/(1-z^2-z^6/(1-q^6)),

    B_0^flow(z,q)=1/2 * [
      A_ind(z,q)*(1+z+z^3/(1-q^3)+z^5 f(q))
      + A_mir(z,q)*(1+z+z^3/(1-q^3)+z^5 k(q)) ].     (4)

This uses the inherited complete classification of underlying degree-zero centers, not a classification of a new decorated metagraph.

## 2. The integral pentagon cone

In coordinates (t,y0,...,y4), the flow cone is

    y_i>=0,  t-y_i-y_(i+1)>=0  for i mod 5.         (5)

The corresponding short-edge coordinates are x_i=t-y_i-y_(i-1). Total use is 5t-sum y_i. This is an integral lattice isomorphism: t=x_i+y_i+y_(i-1) is an integer for every integral circulation, and all inverse coordinates are integral.

Intersect (5) with t=1. The vertices are 0, the five unit vectors, the five independent-pair vectors e_i+e_(i+2), and the all-half vector. A direct proof of completeness is useful. At a vertex, any fractional-coordinate component of the graph of tight adjacent-sum constraints must either be anchored by a zero/integral coordinate, forcing its whole alternating component integral, or have a nonzero perturbation alternating along that component. A path component has such a perturbation and is not a vertex. The only unanchored cyclic component is the full odd pentagon, which forces every coordinate to 1/2. All remaining vertices are 0/1 independent sets of C5, of sizes 0,1,2. The script independently exhausts all 252 possible sets of five active inequalities on t=1 and obtains exactly these twelve vertices.

Use this order for the primitive integer rays:

    r0=(1;0,0,0,0,0),
    r_(1+i)=(1;e_i),                  i=0,...,4,
    r_(6+i)=(1;e_i+e_(i+2)),          i=0,...,4,
    r11=(2;1,1,1,1,1).

Their total-use weights are respectively 5, five copies of 4, five copies of 3, and 5. The last ray is essential: it is the directed Hamiltonian cycle using only length-two arcs. Omitting the fractional vertex would omit genuine flows.

## 3. An exact 25-cone unimodular certificate

Pull the polytope vertices in the displayed order. The following recursive definition fixes the triangulation without a numerical geometric library:

- For a face with a single vertex, return that vertex.
- Otherwise take its least-index vertex v.
- Obtain the facets of the face as the maximal proper nonempty intersections of its vertex set with the ten zero-sets in (5).
- Triangulate each facet not containing v recursively, and adjoin v to every resulting simplex.

This is a triangulation by induction: projecting away from v sends the polytope minus v onto its facets not containing v; the recursive subdivisions agree on their common faces because they use the same vertex order. The same argument applies to every face. The supplied exact certificate checks all recursive facet codimensions.

Coning the resulting triangulation from the origin gives 25 full-dimensional simplicial cones. Every 6 by 6 ray determinant is +1 or -1 in the integral (t,y) lattice. Thus every cone is unimodular, even though the unsliced pentagon polytope has a half-integral vertex.

For disjoint lattice enumeration, use the explicitly checked generic interior vector

    w = sum_(i=0)^11 2^i r_i.

For a simplex with ordered ray matrix M, all coordinates alpha=M^(-1)w are nonzero. Keep the coefficient of ray r_i strictly positive exactly when alpha_i<0, and nonnegative otherwise. These half-open cones partition the whole cone: for each lattice point x, x+epsilon*w belongs to one simplex interior for sufficiently small epsilon>0, and these sign rules describe that simplex. Unimodularity therefore gives one numerator monomial q^shift, where shift is the sum of the weights of the open rays.

The complete 25 simplex lists, signed determinants, all alpha coordinates, open-ray sets, and shifts are stored in `normal.json` under `triangulation.cells`. Their generating functions group into only two denominator types:

    A(q)=1+q^3+3q^4+2q^6+4q^7+q^8+3q^10,
    B(q)=2q^3+2q^6+2q^7+q^9+2q^10+q^13.

The fifteen cones with weights (3,3,4,4,5,5) contribute A divided by their common denominator. The ten with weights (3,3,3,4,5,5) contribute B divided by theirs. Consequently

    R(q)= A(q)/((1-q^3)^2(1-q^4)^2(1-q^5)^2)
         +B(q)/((1-q^3)^3(1-q^4)(1-q^5)^2).

The exact polynomial identity

    A(q)(1-q^3)+B(q)(1-q^4)=P(q)(1-q^5)

proves (1). This is the short cancellation certificate; no truncation or fitted recurrence is used to prove it.

## 4. An independent elementary derivation of k(q)

All five anti-isomorphisms of R5 are conjugate under its rotations. Thus anti-coset Burnside reduces to the fixed flows of i->-i. In the y coordinates those have form (a,b,b,a,c), and

    k(q)=sum_{a,b,c>=0, t>=max(a+b,2b,a+c)}
                  q^(5t-2a-2b-c).                            (6)

Set t=2b+d and a=b+d-u. The bounds become

    b,d>=0,  0<=u<=b+d,  0<=c<=b+u,

and the exponent is 6b+3d+2u-c. Sum first over c, then u, then b,d. As an identity of rational functions the result is

    1/(1-q^-1) * {
      1/(1-q^2) * [1/((1-q^6)(1-q^3))
                         - q^2/((1-q^8)(1-q^5))]
      -q^-1/(1-q) * [1/((1-q^5)(1-q^3))
                         - q/((1-q^6)(1-q^4))] }.

All sums before simplification have nonnegative total exponents; the intermediate q^-1 is only a finite-geometric-sum expression. Simplifying yields (3). This derivation is independent of the 25-cone triangulation.

A nonidentity rotation fixes exactly constant short-edge and long-edge coordinates, with weight 5(x+y). Ordinary Burnside then proves (2).

## 5. Consequences and controls

Only q=1 has pole order six in R and f. Its leading coefficient gives

    [q^m]R(q) ~ 17 m^5/259200,
    [q^m]f(q) ~ 17 m^5/1296000.

Only q=1 has pole order four in k, giving

    [q^m]k(q) ~ 7 m^3/4320.

The nonidentity rotation terms in (2) have only order-two poles and cannot affect the leading term. The raw nonidentity roots have pole orders at most three; the anti-fixed nonidentity roots have pole orders at most two. Hence the fraction of ordinary R5 flow-orbits that are converse-fixed satisfies

    [q^m]k(q) / [q^m]f(q) ~ (2100/17) m^-2.          (7)

The denominators and numerators in (1) and (3) are coprime. Their numerator degrees are 12 while their denominator-weight sums are 22, giving reciprocal shift 10, equal to the ten arcs of R5. Exact checks verify R(q^-1)=q^10 R(q), and the same identity for f, k, and (f+k)/2. This is a cross-check, not a replacement for the cone enumeration or a proof of a general regular-tournament reciprocity theorem.

The script independently enumerates the original pentagon inequalities through total use 40 and the anti-fixed sum (6) through total use 80. It proves the bound t<=2m/5 from summing the five adjacent-sum inequalities, so these are complete slices. It compares the resulting coefficients with (1) and (3), and checks ordinary/converse Burnside integrality. Normal and optimized output files agree byte-for-byte. The all-order formulas come from the exact geometric certificate and geometric sums; the bounded comparisons are hostile controls only.

Scope stops here: one explicit R5 cone, actual automorphism and anti-coset quotients, and insertion into the existing underlying-degree-zero formula. No general R_n free-monoid claim, Euler-tour identification, or solution to an unrelated open problem is inferred.
