# Conway amplitude and golden carry memory

2026-10-13. **CITED:** the classical Pascal-word encoding of Conway's
sequence. **PROVED from that encoding:** its uniform Gaussian-parametrized
hump limit and exact limiting relative-maximum constant. **FINITE-EXACT:**
rational brackets for relative maxima through generation 1024 and binary
golden carries through length 14. **OPEN:** the long-time behaviour and
uniform error of the separate lagged continuum model, Q's sampling
assumptions, and transport to a universal Collatz bound. Self-audited.

The actual Conway sequence has

    lim sqrt(k) max_(2^k<=n<=2^(k+1)) (C(n)/n-1/2)
      = 0.1380149023990716... .                              (1)

The decimal evaluates an exact constant characterized below. This closes
the repo's open item for the *integer sequence*, by recovering classical
structure rather than extrapolating its continuum simulation. It does not
prove the continuum simulation has the same limit. The new Collatz
[source-refinement note](collatz_wavelet_energy_20261013d.md) uses this as a
control for retaining the original source and its denominator.

## 1. Classical inheritance and the exact object

Let `C(1)=C(2)=1` and
`C(n)=C(C(n-1))+C(n-C(n-1))`. The
[previous tree note](golden_nested_trees_20261013c.md) proves the elementary
charge identity: `e(n)=2C(n)-n` counts terminal leaves labelled 1 in the
computation tree. The other terminal label 2 contributes zero. Retaining
the index at each split gives total leaf mass n.

Use words `B_(k,0)=+`, `B_(k,k)=-`, and
`B_(k,j)=B_(k-1,j-1) B_(k-1,j)` for 0<j<k. The classical encoding says
that row k, concatenated left to right, is exactly the increment word of e
on `2^k+1,...,2^(k+1)`. This is explicitly recorded in Appendix B of
[Janvresse, de la Rue and Velenik, *Self-Similar Corrections to the Ergodic
Theorem for the Pascal-Adic Transformation*](https://www.unige.ch/math/folks/velenik/papers/JRV05a.pdf),
with credit there to Mallows. Their Theorem 2.1 resolves finer, centered
within-word fluctuations. The [Kubo–Vakil publisher record](https://doi.org/10.1016/0012-365X(94)00303-Z)
also explicitly credits their proof of Mallows's asymptotic-shape conjecture.
Thus the large-scale limiting shape is classical, not a new open problem.

The proof below is an elementary derivation from that cited encoding.
It does not need an uninspected formula from the Kubo–Vakil full text, or
any probabilistic model for a nested recurrence. Independent literal
recurrence checks authenticate the encoding through 2^18 here (through 2^20
in the preceding packet), and the max-plus DAG is checked against literal
relative maxima at generations4,10,17.

## 2. Boundary positions and charges determine the limit

Put `b_(k,j)=binom(k,j)`. The length and signed total of a word are

    length B_(k,j) = b_(k,j),
    total B_(k,j) = binom(k-1,j)-binom(k-1,j-1)
                   = (1-2j/k)b_(k,j).                       (2)

At the end of cell j, therefore,

    n_(k,j)=2^k+sum_(i=0)^j binom(k,i),
    e(n_(k,j))=binom(k-1,j).                                (3)

The total telescopes. With `z=(2j-k)/sqrt(k)`, Stirling's formula and the
binomial central limit give, uniformly on bounded z intervals,

    n_(k,j)/2^k -> 1+Phi(z),
    sqrt(k)e(n_(k,j))/2^k -> phi_N(z),                       (4)

where `phi_N(z)=exp(-z^2/2)/sqrt(2pi)` and Phi is its distribution
function. The subscript distinguishes normal density from the golden
ratio. Formula (4) alone leaves a gap: a word between two boundaries could
have a large hidden excursion. The following bound closes it.

### A uniform bound inside a Pascal word

For a word of length b and total S, subtract the straight line of slope S/b
from its piecewise linear prefix path; let D be the largest absolute
deviation. At a split `B_(k,j)=A B`, the left endpoint of this centered
bridge has height

    delta_(k,j)=2j(k-j)/(k^2(k-1)) b_(k,j)
                 <= b_(k,j)/(2(k-1)).                       (5)

This follows by substituting (2). It is nonnegative. Each child bridge is
added to a nonnegative linear tent of height at most delta, hence

    0 <= centered bridge,
    D_(k,j) <= delta_(k,j)+max(D_A,D_B).                     (6)

For `3k/8<=j<=5k/8`, follow any branch for `h=floor(k/8)` generations.
At each of these splits both child-length fractions lie between1/4 and3/4:
after d<=k/8 steps the index has fallen by at most d, and the level is k-d.
Consequently every depth-d child has length at most `(3/4)^d b_(k,j)`.
For k>=8 the denominators in (5) satisfy `k-d-1>=3k/4` during these h
splits. Applying (6) h times and bounding the residual bridge by twice
its length gives the explicit estimate

    D_(k,j) <= b_(k,j) [8/(3k)+2(3/4)^floor(k/8)]
             = O(b_(k,j)/k),                                (7)

uniformly in this central band. This estimate is deliberately coarse;
the classical Pascal-adic theorem describes the surviving smaller fractal.

On any fixed bounded z interval, cell lengths are O(2^k/sqrt(k)) and
their means in (2) are O(1/sqrt(k)). Equations (2) and (7) make the
*whole* within-cell excursion o(2^k/sqrt(k)); the boundary mesh in n/2^k
is O(1/sqrt(k)). Thus (4) extends uniformly across every compact subarc
of 1<x<2. For the two tails, the crude bound
`|prefix within cell|<=b_(k,j)`, together with monotonicity of binomial
coefficients away from the center and their local normal estimate,
bounds the scaled tail by O(exp(-R^2/2)) for |z|>=R, in the successive
limits k->infinity then R->infinity. Both endpoint values are zero.
This proves uniform convergence on the full closed arc:

    sqrt(k)(C(2^k x)/(2^k x)-1/2) ->
       g(x)=phi_N(Phi^(-1)(x-1))/(2x),   1<x<2,
    g(1)=g(2)=0.                                            (8)

Here the charge is linearly interpolated between integer n; dividing by
n afterwards preserves the location of an interval maximum at an endpoint.
So the maximum in (8) is also the integer-source maximum.

## 3. Exact amplitude constant and the continuum hostile

Maximizing (8) is equivalent to maximizing

    f(z)=phi_N(z)/(2(1+Phi(z))).                             (9)

Its derivative has the opposite sign to
`F(z)=z(1+Phi(z))+phi_N(z)`. Since `F'(z)=1+Phi(z)>0`, F has exactly one
zero z_*, lying in (-1,0). Both tails of f vanish, so this is the unique
global maximum. At its critical point,

    C_*=-z_*/2,
    z_*(1+Phi(z_*))+phi_N(z_*)=0.                            (10)

Numerically `z_*=-0.2760298047981433`, `C_*=0.13801490239907166`, and the
limiting peak is at `x_*=1+Phi(z_*)=1.3912625871346345`. These displayed
decimals are ordinary numerical evaluations; (10) is the exact answer.

The [max-plus checker](../../04-computation/experiments/conway_golden_memory_20261013d.py)
uses integer arithmetic to bracket the actual A_k, without enumerating
2^k sources. For a trial rational lambda, give each letter weight
`sign-lambda`. The word summary `(length,total,max prefix)` composes as
`max_AB=max(max_A,total_A+max_B)`. Subtracting `lambda*2^k` then decides
whether any source has `e(n)/n>lambda`. Bisection on a fixed rational grid
gives certified A_k brackets. The scaled decimals below display those
brackets after multiplication by sqrt(k).

| k | lower display for sqrt(k) A_k | upper display |
|---:|---:|---:|
| 20 | 0.1461362382 | 0.1461362404 |
| 29 | 0.1438076293 | 0.1438076319 |
| 128 | 0.1390872079 | 0.1390872132 |
| 400 | 0.1384044438 | 0.1384044532 |
| 1024 | 0.1381513774 | 0.1381513924 |

The authoritative exact rational endpoints are in the
[saved JSON](conway_golden_memory_20261013d.json). The k400 bracket already
separates the true sequence from the pinned continuum model's roughly 0.137.
The earlier inference “the limit, if positive, is at most0.137” is invalid
even for a model from a finite decreasing segment, and false for the actual
sequence by (10). The integer sequence's limit and shape are now closed;
the continuum map's own asymptotics and approximation error remain open.
The earlier fit0.1377 is a finite-window fit, not the exact constant.

### A rate for the true maximum

The same argument gives the quantitative asymptotic

    sqrt(k) A_k = C_* + O(1/k).                             (11)

Here the implied constant is not supplied as an effective finite-k
certificate. Use the continuity-corrected coordinate
`z_c=(2j+1-k)/sqrt(k)` in (3). Stirling's expansion and midpoint summation
give, uniformly for bounded z_c,

    sqrt(k) binom(k-1,j)/2^k = phi_N(z_c)+O(1/k),
    2^(-k) sum_(i<=j) binom(k,i) = Phi(z_c)+O(1/k).

For completeness, the midpoint error is O(1/k) because the mesh is
2/sqrt(k) and the normal density has integrable second derivative. The
local binomial error is bounded by a constant times
`k^(-1)(1+z^4)phi_N(z)` times that mesh on |z|<=k^(1/8); the omitted tails
are exponentially small. The cell endpoints therefore have scaled
relative charge `f(z_c)+O(1/k)`. A grid point lies within1/sqrt(k) of z_*;
since f'(z_*)=0, their maximum is `C_*+O(1/k)`.

Uniform convergence already confines all maximizing cells to a fixed
central band for large k. Within such a cell, the straight chord's
relative charge is bounded by its endpoint maximum (an affine numerator
over the positive affine source coordinate is monotone). Bound (7) adds
only O(1/k) after scaling by sqrt(k)/(2n). This proves (11), including the
possibility that the true maximizer is inside a cell. The smaller
Pascal-adic fractal can affect this correction without changing C_*.

### The formal nonlinear eigenprofile is explicit

The repo's formal diffusion operator was `2h^2(x^2 h')'`. The limiting
profile just recovered satisfies its self-similar equation exactly:

    g^2(x^2 g')' = -g/4,  g(1)=g(2)=0,  g>0 on (1,2).     (12)

This closes the separate *existence of a positive eigenprofile* question.
To see the mechanism, set u=xg. Equation (12) becomes

    u u''=-1/4.                                            (13)

Alternatively, from `x=1+Phi(z)` and `u=phi_N(z)/2`, direct differentiation
gives `u'=-z/2` and `u''=-1/(4u)`. Thus no statistical assumption is needed
to obtain the nonlinear equation from the deterministic binomial profile.

It is the unique positive solution of (13) continuous on [1,2], zero at
both ends and twice differentiable inside. Concavity gives one maximum
u_0, and the first integral is

    (u')^2 = (1/2) log(u_0/u).

The time from zero to the maximum is `u_0 sqrt(2pi)` on either side;
the interval has width one, so `u_0=1/(2sqrt(2pi))` and the maximum is at
x=3/2. This pins the whole solution. Dividing by the retained source x
moves the *relative* maximum to x_*=1.391262587... and gives (10).
For a profile normalized to maximum one, the eigenvalue is
`1/(4 C_*^2)` instead of1/4. The similarity form `h_k(x)~g(x)/sqrt(k)`
has `h_k-h_(k-1)~-g/(2k^(3/2))`, matching the operator's leading value.
This match and profile existence do not prove convergence of the lagged
scheme: g has singular endpoint slopes, and a uniform consistency and
stability estimate is still required.

## 4. Golden 11 as an exact shift with memory

The inherited [Lucas clock note, section1](collatz_lucas_monotile_discrepancy_20260930.md)
already proves the requested identities. Multiplication by the golden ratio
in the basis (1,phi) is

    M = [[0,1],[1,1]],  M^2=M+I,
    M^4=[[2,3],[3,5]],  trace 7, determinant 1,
    M^5=[[3,5],[5,8]],  trace 11, determinant -1.

Cayley–Hamilton gives `M^8+I=7M^4` and `M^10=11M^5+I`, hence
`phi^8+1=7phi^4` and `phi^10=1+11phi^5`. The decimal integer11 here is a
Lucas trace. The digit string11 below is a different object, tied to the
same quadratic relation rather than identified with that integer.

For a finite binary word read in base phi, retain its value as the exact
pair (a,b) representing a+b phi. Appending a bit d performs

    (a,b) -> (b+d,a+b).                                     (14)

This is literally a shift carrying a second coordinate. No irrational
rounding is needed. The local rewrite `011 -> 100` preserves the pair.
Choose the leftmost11, padding with a leading zero if necessary: its
preceding bit is zero, so the rewrite remains binary and decreases the
number of ones by one. It terminates at a word with no consecutive ones.
Overlaps require this context: `111 -> 1001`, not an untracked replacement
of every occurrence simultaneously.

The resulting finite binary normal form is unique. If two such words first
differ at their highest differing power h, the smaller side's largest
possible lower tail is bounded by
`phi^(h-1)+phi^(h-3)+...<phi^h`; it cannot cancel the leading difference.
Only finite words in nonnegative powers are being considered. The rewrite
preserves exact value, not the original spelling: retain a rewrite DAG
when the history itself is the certificate.

This gives a precise comparison with the new Collatz packet. Its address
carry has `B -> B` or `B -> 3B+2^j`; forgetting B and retaining only the
number of each branch loses the endpoint. Golden normalization can safely
forget spelling for a value query; a Collatz ROOT or energy query must
also retain the marked source, terminal relation and rounding carry.
The common object is an affine action plus an additive carry, with
different relations and different admissible quotients.

For Q, the [proved covariance-model calculation](golden_nested_trees_20261013c.md#6-incoming-q-covariance-work-a-genuine-two-step-pell-mechanism)
has the two-generation matrix `[[3,4],[2,3]]`, a Pell unit of determinant 1.
The golden matrix is a unit of determinant -1 on a different quadratic
lattice. Both admit exact integral memory and characteristic-polynomial
identities. This comparison transfers a representation method, not the
unproved white-sampling premise or a covariance value.

## 5. What this suggests for Collatz

The Conway result succeeds because the tree supplies all three needed
coordinates: cell mass, signed charge, and a bound on the bridge inside
each cell. Boundary values alone would have been insufficient. The
Collatz wavelet packet now has the corresponding exact source and boundary
rules, but still needs its analogue of (7), with actual correlations.

The useful helper question is therefore: can a depth-dependent Collatz
packet have a deterministic bridge bound, after its source-specific affine
part is subtracted? A generic Brownian assumption is unnecessary in the
Conway control. Conversely, the existence of a deterministic limiting
shape for Conway does not supply one for Collatz or Q. Test any proposed
bridge norm against both signs of the cross-term in the companion note.

```sh
python3 04-computation/experiments/conway_golden_memory_20261013d.py --output 05-knowledge/results/conway_golden_memory_20261013d.json
```

The run brackets nine Pascal generations through 1024, checks 32,766 binary
words through length 14 for exact carry and canonical consistency, and
compares literal Conway values through 2^18 with the independent word DAG.
Normal and optimized Python runs are required to agree. It does not
enumerate the 2^1024 sources individually.
