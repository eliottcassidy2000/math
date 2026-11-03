# A signed doubling lift of the 290 wheel

**PROVED O1–O4; independently audited PASS.**
[Independent audit](../audit_orientations/REPORT.md). Classical character
orthogonality and cyclotomic degree are CITED; the elementary transport and
scope proofs are supplied below. No historical-priority claim is made.

## Inheritance and the exact question

The closest proved mechanism is [PRIME-MRB](../../jacobian_eth_prime_wait_20261102/PRIME-MRB.md):
it already proves the290 prime-delay certificates and the MRB mean calculation.
The canonical hostile is treating increasing primes as independent uniform
samples of U(m)/{±1}. The corrected near miss is its56-object coincidence
with tournament order6: a290-edge graph on56 vertices cannot be a regular
Cayley graph because its average degree145/14 is not integral. The least-used
coordinate is the sign erased by inversion, now retained along an actual
multiplicative orbit. [The arithmetic-braids geometry note](../../arithmetic_braids_20260917_geometry.md)
already gives the exact seventh-root trace map and quadratic cycle reversal.

Our source is the multiplicative unit group of an odd conductor, our map is
doubling, and the target is its inversion quotient or its cyclotomic trace.
Neither map advances the increasing-prime clock. The cheapest decisive test
is whether half an orbit reaches -1: it does at145 and does not at7.

## O1. The two possible signed doubling geometries

Let N>1 be odd, U_N=(Z/NZ)^×, d=ord_N(2), and let negation J act on U_N.
Let V_+ and V_- be the complex functions satisfying f(-r)=+f(r) and
f(-r)=-f(r), respectively. Let Tf(r)=f(2r).

If -1 belongs to <2>, then d is even and 2^(d/2)=-1. Every unit orbit
has d elements, its inversion quotient has d/2, and

    charpoly(T|V_+) = (X^(d/2)-1)^(phi(N)/d),
    charpoly(T|V_-) = (X^(d/2)+1)^(phi(N)/d).             (O1)

If -1 does not belong to <2>, negation exchanges distinct doubling orbits.
The quotient retains period d, and both characteristic polynomials are

    (X^d-1)^(phi(N)/(2d)).                               (O2)

**Proof.** A unit r returns under doubling exactly when2^k=1, so every
orbit has length d. If2^k=-1, its square is1 and the unique element of
order2 in the cyclic group <2> is2^(d/2). Negation is then the halfway
shift. On even functions this closes a length-d/2 cycle; on odd functions
it closes with sign -1. Their shift matrices have the polynomials in(O1).
Otherwise an orbit and its negative are disjoint, and a function of either
parity is freely specified on one of the two length-d cycles, proving(O2).
This also proves all asserted divisibility and multiplicity statements. QED.

## O2. Why145 and28 meet exactly

Reduction modulo145 identifies U(290) with U(145): every residue modulo145
has a unique odd lift modulo290. Under this isomorphism doubling is
multiplication by147, the odd lift of2. Literal multiplication by2 modulo290
would leave the units and is not the operation being used.

Modulo5,2 has order4. Modulo29,

    2^7=12,  2^14=-1,  2^4=16 !=1,

so its order is28: the proper divisors not already excluded by2^14≠1
divide4. CRT gives ord_145(2)=28 and2^14=-1 modulo both5 and29.
Consequently the112 unit residues form four28-cycles. Folding r and-r
leaves four14-cycles, with

    charpoly(T|V_+)=(X^14-1)^4,
    charpoly(T|V_-)=(X^14+1)^4.                         (O3)

In the oriented representation T^14=J; on the signed sector T^14=-I.
This is an exact reason to distinguish the two half-orbit endpoints.
It is a finite permutation operator spectrum, not a tournament adjacency
spectrum, PRIMEGAME running time, or prime-arrival distribution.

At N=7,ord_7(2)=3 and-1 is outside<2>; the two oriented3-cycles fold to
one3-cycle and both sectors have polynomialX^3-1. At N=9,ord_9(2)=6 and
2^3=-1; the folded3-cycle has even polynomialX^3-1 and odd polynomial
X^3+1. Thus a common three-cycle picture does not determine its signed lift.
The F8 multiplicative clock has order7 and Frobenius has order3; they are
separate operations. Only the latter residue-doubling clock appears at N=7.

## O3. MRB convergence and an explicit signed tail bound

Put f(x)=exp(log(x)/x)-1 for x≥1. The inherited periodic-mask theorem says
that a periodic coefficient a_n of mean μ satisfies

    Σ_(n≤M) a_n f(n) = (μ/2)(log M)^2 + C + O(log² M/M).   (O4)

The constant and error may depend on the fixed period and coefficients.
The displayed error follows directly by writing a_n=μ+b_n. Bounded
partial sums of b_n give a convergent b-sum with tailO(f(M)); moreover
f(n)-log(n)/n=O(log²(n)/n²), whose tail isO(log² M/M).
The integral comparison for log(n)/n has errorO(log M/M). The leading term and
convergence criterion are inherited, not a new MRB formula.

For the unit-filtered alternating series, a_n=-1 on U290 and0 elsewhere:
its mean is-112/290=-56/145. Therefore

    Σ_(n≤M,gcd(n,290)=1) (-1)^n f(n)
       = -(28/145)(log M)^2 + C + O(log² M/M).           (O5)

For the two unsigned target residues the corresponding coefficient is
-1/290. If instead a_n=1_(n=1 mod290)-1_(n=-1 mod290), its period sum
is zero and the series converges. Its ordinary prefix sums are always0
or1. For every integer M≥3, summation by parts therefore gives the
particularly simple rigorous tail estimate

    |Σ_(n>M) a_n f(n)| ≤ f(M+1).                        (O6)

Indeed f is positive, decreasing for x>e, and tends to0. Partial sums
of the tail's coefficients are in{0,1} or{0,-1}; Abel summation bounds
the entire tail by its first f-value. At block cut M=290K, the tail is
nonnegative, since each subsequent pair is
f(290j+1)-f(290j+289)>0 for j≥K≥1. The signed sum is not absolutely
convergent: absolute values have mean2/290, giving a positive log² term.

The character formula makes the discarded information explicit:

    1_(n=1)-1_(n=-1) = (2/112) Σ_(chi(-1)=-1) chi(n)

on the290 residue classes, with characters extended by0 off the units.
[Dirichlet character orthogonality](https://dlmf.nist.gov/27.8#E6) proves it;
there are56 characters of each parity. The indicator difference is an
antisymmetric observable, not itself a multiplicative character. The quotient
by±1 erases it, whereas the lift inO2 transports it exactly. Any antisymmetric
periodic mask has zero mean, but zero-mean symmetric masks can also converge:
antisymmetry is sufficient, not necessary.

## O4. From the wheel to quadratic dynamics, with a firm boundary

Let ζ_N be a primitive Nth root of unity and t_r=ζ_N^r+ζ_N^(-r). Then

    D(t_r)=t_r²-2=t_(2r),
    t_r=t_s iff r=±s modN.                              (O7)

The second statement follows becauseζ_N^r andζ_N^(-r) are the two roots
of Z²-t_r Z+1. Thus the trace is an exact realization of the folded wheel,
with no further collisions among primitive exponents. Its period is the
least k>0 with2^k=±1 modN. In particular the145 trace set has56 points
in four exact14-cycles, and the7 trace set has three points in one3-cycle.

The sign has an exact differential manifestation as well. Differentiating
D^k(z+z^(-1))=z^(2^k)+z^(-2^k) gives

    (D^k)'(t_r)=2^k (z^(2^k)-z^(-2^k))/(z-z^(-1)),
    z=ζ_N^r.                                            (O7a)

At the least folded return k this is +2^k or-2^k according as the
oriented lift returns to+1 or-1. Therefore all four145 trace cycles have
multiplier-16384; the7 and9 trace cycles have multipliers+8 and-8.
This recovers return monodromy from the differential of the actual D map;
a bare cycle graph, or an arbitrarily chosen affine transport, lacks it.

The standard cyclotomic degree theorem gives degreeφ(N) forζ_N;
see [Conrad, Cyclotomic Extensions, Theorem2.5](https://kconrad.math.uconn.edu/blurbs/galoistheory/cyclotomic.pdf).
Sinceζ_N is nonreal and quadratic over the real trace field, t_r has
rational degreeφ(N)/2. These degrees are56 and3 respectively. The checker
also constructs the integral trace polynomials by

    Φ_N(Z)=Z^(phi(N)/2) R_N(Z+Z^(-1))

and checks D's exact periods by polynomial reduction. This exact computation
is an independent control, not a replacement for the degree theorem.

For N=7, R_7(t)=t³+t²-2t-1. Put L(t)=-t-1/2. The inherited identity is

    f_(-7/4)(L(t))-L(D²(t))=(t-1)R_7(t).               (O8)

So L transports the reversed chosen7-trace cycle to the parabolic3-cycle
of x²-7/4. This holds on that chosen three-point set; it is not a global
conjugacy of degree2 and degree4 maps. The quadratic cycle has multiplier1,
while the Chebyshev cycle has multiplier8, so derivatives are not preserved.
For R_9=t³-3t+1 the same reversal gives parameter-11/4 and multiplier19,
a decisive hostile to identifying all folded three-cycles with-7/4.

Finally the inherited chosen-cycle trace parameterη=2σ+1 satisfies

    c=-(η²+7)/4,       reversal η↦-η-1.

The unique fixed trace isη=-1/2, hence c=-29/16=-7/4-1/16. Its rational
cycle is -7/4→5/4→-1/4→-7/4, with multiplier35/8. This is the precise
28+1=29 relationship:16c=-28-1 at that fixed trace. It is an equation
in a marked quadratic family, not a deduction from the count28 alone.
There is no affine cycle conjugacy from the145 trace set to this rational
three-cycle: the exact periods14 and3 already disagree. The fact29 also
occurs in145 supplies a modulus to test, not an identification of dynamics.

## Connection ledger and reproduction

| Map | Preserved | Lost or required |
|---|---|---|
| U290→U145 | multiplicative group, unique odd lift | the native integer and increasing-prime clock |
| Unit→unit/± | folded doubling, exact period as inO1 | sign sector and antisymmetric MRB readout |
| Unit/±→cyclotomic trace | exact finite D dynamics | oriented lift; real derivatives of other maps |
| R7 trace→chosen-7/4 cycle | reversed three-point dynamics, algebraic degree | global map and multiplier |
| Fixed chosen trace→-29/16 | rational three-cycle and its reflection | no145-cycle correspondence |

Run check.py normally and with-O. Exact gates cover group orders, orbit
partitions, signed-polynomial relations, character parities, trace-polynomial
reductions and the two quadratic controls. Numerical MRB sums are exploratory
only; convergence and(O6) follow from the proof. Prime arrivals and odd-pair
rows are independently guarded in the companion machine lane.
