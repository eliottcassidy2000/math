# Why six is balanced and seven is the first surplus

**Status: PROVED in the stated smooth complete-intersection scope; the
six-line matrix identity below is an exact independent calculation.**
Adjunction, Koszul vanishing, and Gale duality are classical mechanisms.
This note packages their application to the current root models and claims
no historical priority for those mechanisms.

BG1 is also recovered in concurrent [Session60, G1](../session60_six_seven_geometry_20261025/thresholds/ROOT_THRESHOLDS.md); the two-parameter outer BG2 identity below is separate from its central Kummer map.

## 1. A common parameter behind several integers

Let d>=1, r>=2, and n>=d+2. Over an algebraically closed field of
characteristic zero, take n hyperplanes in P^d in linear general position.
Their forms L_i span the d+1-dimensional parameter space. Pull that linear
P^d back under coordinatewise r-th powers in P^(n-1), obtaining X.
Write c=n-d-1 for the number of independent relations and H=O_X(1).

**Theorem BG1.** X is smooth and integral, is a complete intersection of c
degree-r equations, and has

    dim X=d,  H^d=r^c,  K_X=(rc-n)H,
    rc-n=(r-1)n-r(d+1).

Its Hilbert series is (1-t^r)^c/(1-t)^n. For d>=2,
H^i(X,O_X)=0 for 0<i<d. In particular, when rc=n its canonical bundle is
trivial and its only nonzero H^i(O_X) occur at i=0,d.

**Proof.** Choose d+1 independent forms. The other root coordinates obey
monic degree-r equations over the polynomial ring in the chosen root
coordinates. These equations form a regular sequence and give degree r^c.
Over the generic parameter field the r-th root extensions have full degree:
valuation along each distinct numerator hyperplane detects its exponent
modulo r. This also proves geometric integrality. At a point at most d
hyperplanes vanish. General position makes their equations independent
local parameters. Taking their r-th roots and etale roots of the unit
equations leaves a regular local ring, proving smoothness. Adjunction gives
K_X=(rc-n)H. The Koszul resolution gives the Hilbert series and intermediate
cohomology vanishing. These arguments do not concern arbitrary singular
root arrangements. In the actual surface arrangement, Session59 separately
proves that its nodes admit crepant resolutions. QED.

For square roots r=2, the canonical coefficient is

    c-(d+1)=n-2d-2.

It is literally the **number of relations minus the number of parameter
coordinates**. Gale duality exchanges these two dimensions. Therefore the
quadratic trivial-canonical threshold is exactly the dimension-balanced
Gale threshold n=2d+2. Balanced dimensions alone do not specify an
identification of the two configurations.

| Dimension d | Balanced square-root model | First model with K=H |
|---:|---|---|
| 1 | n=4, genus-one intersection of two quadrics | n=5, canonical genus-five curve |
| 2 | n=6, K3 intersection of three quadrics | n=7, canonical surface, p_g=7 and K^2=16 |
| 3 | n=8, trivial-canonical intersection of four quadrics | n=9, canonical threefold with h^0(K)=9 and K^3=32 |

In general the first surplus model n=2d+3 has K=H,
h^0(K)=n and K^d=2^(d+2). The sections are the n coordinate linear forms;
there are no degree-one equations. Thus the seven and sixteen in the
hourglass are consequences of one codimension calculation, not independent
numerical coincidences.

**Hostile change of operation.** For fourth roots on a surface, n=4 already
has K=0: a smooth quartic K3 in P^3. Its parameter dimension is three and
relation dimension one, so it is not Gale-balanced. The equality of the
canonical and Gale thresholds is specific to square roots. For cubics the
first surface with K=H has n=5 and K^2=9. There is no operation-independent
meaning of the bare integers six and seven.

## 2. An actual Gale return, with coefficients retained

For u,v nonzero, take the six parameter forms

    E, E+x, E-x, E+z, E-z, E+u x+v z.

Their coefficient matrix and one relation matrix are

    A = [[1, 1, 1, 1, 1, 1],
         [0, 1,-1, 0, 0, u],
         [0, 0, 0, 1,-1, v]],

    C = [[-2,1,1,0,0,0],
         [-2,0,0,1,1,0],
         [u+v-1,-u,0,-v,0,1]].

They have rank three and A C^t=0. The columns of C are the Gale-associated
configuration. Put

    M = [[u,v,2],[u,0,0],[0,v,0]],
    p = (5,2,1,4,3,0),
    D = diag(-2,-u,u,-v,v,2).

**Theorem BG2.** M C=A_p D, where column i of A_p is column p(i) of A.
Since det M=2uv, the two labelled configurations are projectively
equivalent after exactly the displayed permutation and column rescalings.
The identity holds over every characteristic-zero field containing u,v.

**Proof.** Direct multiplication gives the six columns
(-2 A_5,-u A_2,u A_1,-v A_4,v A_3,2 A_0). Nonzero u,v make all scalars
and the determinant invertible. QED.

The actual outer-deletion intermediate has u=v=1. A complete permutation
check finds eight matroid-compatible permutations and two projective
matches for that configuration. The other matching exchanges x and z.
Those are finite statements about the displayed rational configuration,
not the moduli of all six-line arrangements.

## 3. Net discriminants: which K3 is on which side?

The six-root surface is the intersection of the three diagonal quadrics
given by the rows of C. The determinant of a member of their net is

    det diag(C^t lambda)=product_i (C_i dot lambda).

Consequently its discriminant sextic is the union of the six Gale lines.
BG2 identifies this branch arrangement with the original six-line
arrangement. The exact product of the six rescalings is

    product_i D_i=-4u^2 v^2.

Thus after the displayed base change and division of the cover coordinate
by 2uv, the net double plane has equation w^2=-product_i L_i. Over an
algebraic closure it is isomorphic to the original six-line double plane;
over Q(u,v) this displayed identification retains the **minus-one quadratic
twist**. No claim of nonisomorphism under every possible other map is made.
This calculation does not identify the degree-eight K3
intersection itself with either degree-two K3 double plane.

The six-root intersection also maps to the double plane of the original
six lines: the character selecting all six branch labels gives a degree-16
quotient of its degree-32 map to P^2. These are two explicit ways the same
six-line double-plane geometry enters. A polarization, transcendental
lattice, or derived equivalence is additional structure, not supplied by
the word “K3” or the equality of p_g.

The classical smooth-net correspondence is documented by Madonna–Nikulin,
[On a classical correspondence between K3 surfaces](https://arxiv.org/abs/math/0206158).
Here only the displayed determinant and field maps are asserted for the
singular special net; no smooth-net moduli theorem is imported into that
specialization. The standard complete-intersection K3 examples and
six-line double planes are recorded in Huybrechts,
[Lectures on K3 surfaces, Chapter 1, Example 1.3](https://www.math.uni-bonn.de/people/huybrech/K3Global.pdf).

## 4. Connection contract and controls

- Source: root equations, with every coefficient and root order retained.
- Target: canonical degree and the associated relation configuration.
- Map: regular-sequence adjunction, followed by C=ker A and explicit M,p,D.
- Preserved: dimension and the six-line discriminant arrangement over the
  stated field. Lifting the individual column rescalings to root coordinates
  requires adjoining their square roots; it is not asserted over Q(u,v).
- Lost by merely recording the numbers: branch incidence, integral square
  conditions, coefficient squareclasses, polarization and rational points.
- Sidecars: root order r, matrix A, permutation p, diagonal D, base field,
  and singularity discrepancies. The minus-one twist is an explicit
  arithmetic sidecar, not an unspecified scalar that may be discarded.
- Controls: change r to four; set u or v to zero and watch M become
  singular; independently enumerate all permutations at u=v=1.

Reproduce with `python3 -B 04-computation/six_seven_gale_cube_20261025/verify_budget_gale.py`.
The finite checker corroborates identities; the proofs above quantify
beyond its finite parameter grid.
