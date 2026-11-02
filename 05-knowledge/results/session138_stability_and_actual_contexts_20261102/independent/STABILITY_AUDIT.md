# Independent Session138 stability audit

Status: ACCEPT within the domains below; analytic proofs and finite exact
certificates are distinguished. No all-order optimal coefficient is claimed.

## H=1 locus and scalar optimization (all orders)

The supported-completion characterization, fractional-edge matching argument,
homogeneous contraction and transitive quotient are valid. Binary uniqueness
is most cleanly justified by the odd-cycle formula, not a local rerouting
claim whose exterior attachments have not been checked. The scalar derivative
is `C(1+S0)-1-epsilon D S0`, strictly negative when another positive pair
variance exists. Thus the advertised sharp H=1 constants and equality domains
hold. These are upper bounds on the global optimum until a global lower bound
is separately proved. The separate-edge determinant interpolation signs also
check directly: both u and v acquire `-q_e` times the pair-deleted coordinate.

## Independent degree-two Bernstein certificates

`audit_n5.py` reconstructs the complete signed coefficient tensors directly
from local degree-two Bernstein factors, avoiding producer power-basis
assembly and conversion. A single directed arc gives
`(2,2+epsilon,2+2epsilon)/4`; a reciprocal pair gives `(1,1,0)/4`; an absent
edge gives `(4,4,4)/4`. Tensor products and signed determinant sums therefore
have denominator `4^10`. All24 tensors match the producer's arrays after the
exact factor32 dictated by its denominator `2^(5+10)`. All coefficients are
nonnegative. Independent permutation action partitions all1024 orientations
into12 orbits, with sizes120,40,120,40,120,120,120,40,120,40,120,24.
A Python arbitrary-precision/object replay of the transitive class matches,
and24 Fraction evaluations with independent Gaussian determinants match.
A conservative bound on every direct arithmetic intermediate is3093299200.

This proves the universal order-five bound `Delta>=7V/8`. Its positive-V
zero set is exactly an ordinal sum of two fair pairs and one singleton.
The zero-support computation checks all relatively open cube faces, not
only coefficient-index corners: no positive-dimensional zero face survives.

`audit_n6.py` independently assembles Hamilton and determinant polynomials
by subset dynamic programming (rather than permutation enumeration), and
uses a blockwise tensor conversion. Its56 independent orientation orbits
partition all32768 labelled tournaments. All112 certificate minima, zero
counts and argmins agree with the producer. Every coefficient is nonnegative.
The conservative power-assembly/transform bound is18166579200 (<2^35), so
int64 overflow cannot affect the computation. Support sizes511,826,826,16
and absolute coefficient sums11520,1912,1610,480 provide additional receipts.
The order-six positive-V equality set is exactly an ordinal sum of three
fair pairs. The transitive orthant has1 zero for one sign and9 for the other;
the latter consists of the binary point and the cube of indices0/1 on the
three adjacent matching edges, which permits only the fully fair matching
as a nonbinary actual zero.

## All-order triangle reserve

Let tau be the total product-weight of directed triangles. The root
recurrence retains length-one old paths for variance and length-two old
paths for triangles; these are disjoint summands. In the weak recurrence
`2H_R-u_R>=H_R>=1`, including the empty remainder. In each signed strong
recurrence `2H_R-u_R+/-v_R>=H_R>=1` for nonempty R. Thus

    delta >= tau+c_n V,  n>=2;
    Delta >= tau+beta_n V, n>=4,

with the previously proved Session137 coefficients. Average the parent
triangle count plus triangles through the root; their sum is exactly tau.
The strong n>=4 domain is essential: binary C3 has Delta0 but tau1.
If a nearest binary rounding contains a directed triangle, its three chosen
weights are each at least1/2, hence Delta>=tau>=1/8 for n>=4.

`audit_triangle_reserve.py` independently uses exact Fraction subset
Hamilton DP and Gaussian determinants. All888 matrices (all ternary-edge
matrices through order4, and32 fixed-seed rational cases each order5..8),
3835 weak root bounds and3748 strong root bounds pass. Normal and optimized
Python receipts are identical. The semantic digest is
`85e7a7dfa4969b25240185351a6217b1e7d725ed24a9e26c9c01c17d72612672`.

## Exact order-five triangle/variance region

`audit_joint_n5.py` directly reconstructs all24 Bernstein tensors for
`Delta-tau-13V/16`, using the local-factor method above. All coefficients
are nonnegative and match the producer's minima/zero counts/argmins at
scale factor32. For V>0 equality is precisely a fair3-block ordinal fair2-
block, in either order. This also witnesses why simultaneously claiming
`Delta>=tau+7V/8` would be false: H3/2,u3/8,v1/16,tau1/4,V1 gives
`Delta-tau=13/16`.

The analytic corner `Delta>=3tau/4+7V/8` is independently checked as follows.
Let Ck count weighted directed k-cycles; M22 counts unordered disjoint
variance-edge pairs, M23 counts disjoint variance-edge/triangle pairs.
For n5:

    H=1+2tau+2C5,
    u=1-V+tau-C4+C5+M22-M23,
    v=M23-C5.

For a triple R, let Q_R be its total edge variance and s_R its cyclic
skew sum. Then `Q_R=2tau_R+1/4-s_R^2/4`. Since
`sum_e q_e Q_(e^c)=2M22`, write `S=sum_e q_e s_(e^c)^2` to obtain
`V/8-M22=-M23+S/8`. The two signed corner residuals are

    C4+2C5+sum_R tau_R(1/4-q_(R^c))+S/8,
    C4+tau/4+M23+S/8.

Every term is nonnegative. Together with the certified corner(1,13/16),
convexity and downward closure prove the exact real coefficient region
for `Delta>=a tau+b V`:

    a<=1, b<=7/8, a/4+b<=17/16.

Necessity is witnessed respectively by C3 with two ordinal singletons,
two fair pairs with one singleton, and fair3 ordinal fair2.
The corner(3/4,7/8) has a continuous equality family: any weighted3-block
with cyclic skew sum zero, ordinal a fair2-block. It must not be mistaken
for only the isolated equality sets at the other corners.

## Analytic order-six sharp joint bound

Let M24 count disjoint variance-edge/4-cycle pairs, M33 unordered disjoint
triangle pairs, and M222 unordered triples of disjoint variance edges.
At order6 the independently checked expansions are

    H=1+2tau+2C5+4M33,
    u=1-V+tau-C4+C5-C6+M22-M23+M24+M33-M222,
    v=-C6+M24+M33-M222.

For `S=sum_e q_e sum_(R subset e^c,|R|3) s_R^2`, the triangle identity
sums to `4M22=2M23+V-S/4`: each disjoint edge pair is counted four times.
Hence `M22=M23/2+V/4-S/16`. The two signed residuals after subtracting
`tau+3V/4` become

    C4+C5+M23/2+4M33+S/16,
    C4+C5+2C6+M23/2-2M24+2M33+2M222+S/16.

Finally `M24<=C4/4`, since a4-cycle's complementary edge has variance at
most1/4. Both expressions are nonnegative. This proves analytically
`Delta>=tau+3V/4`; three fair pairs show the variance coefficient is sharp.
The separate Bernstein audit above verifies the weaker sharp variance
bound and its equality locus independently of this cycle proof.

The order-six full coefficient region is therefore the rectangle
`a<=1,b<=3/4`: C3 with three ordinal singletons forces a<=1, and three fair
pairs force b<=3/4. The joint bound gives sufficiency by downward closure.
The same rectangle at order4 follows from the inherited beta4=3/4 triangle
reserve and its analogous witnesses. Order5's slanted facet is a real contrast.

Audit scripts now use active RuntimeError checks (not optimization-disabled
assertions) and resolve the producer folder relative to this directory. Both
n5 audit logs agree byte-for-byte with their Python -O replays.
