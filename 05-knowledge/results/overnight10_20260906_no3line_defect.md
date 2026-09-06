# An exponential bound from the diagonal excess defect

Status: **PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED**.
The [independent referee](overnight10_20260906_no3line_exponential_audit.md)
accepts the full proof and reconstructs the rook means, conditional ranges,
factorial deficits, and finite constants by separate exact paths.
No scarce theorem identifier is requested.
This is a probability theorem about random labels, not an extremal
nonexistence theorem. No external priority claim is made.

The parent research lane proposed replacing the unbounded local triple count
by the positive-part occupancy defect. This report independently derives the
concentration argument, an exact cycle-profile mean, uniform finite mean
bounds, and standalone literal controls.

## 1. Inheritance and the source-target contract

The closest proved mechanism is the uniformly controlled line-occupancy
expansion in
[the eighth diagonal-density theorem](overnight8_20260906_no3line_diagonal_density.md).
Its slope-one triple statistic has linear mean and variance and gives a
uniform polynomial probability bound. The newer
[two-diagonal calculation](overnight9_20260906_no3line_two_diagonals.md)
retains the shared-cell contribution to the second-order constant. This
remains distinct information even after its probability consequence is
strengthened below.

The canonical hostile is n=4: C8 and 2C4 have the same full triple mean, but
their probabilities of no three in line are 1/36 and 1/2. The corrected near
miss is the nonpositive formal type-weight functional in
[the cycle-defect note](overnight4_20260906_no3line_cycle_defect.md).
The least-used sidecar is the exact matching polynomial of the shore-labelled
cycle skeleton. The live board is: actual label permutations; line occupancy;
positive-part excess; factorial-moment truncation; conditional range length;
the distinction between selected-direction and full success.

The map sends a positively distributed saturated board to its occupancies
on the slope-one diagonals, then to their excess above two. It preserves
the zero event for triples in that direction exactly; full no-three-in-line
success implies this zero event. It discards all other slopes and much of the
cycle profile. Unlike the cubic triple count, its one-point changes are
bounded, and the two label permutations remain available as the concentration
source. No theorem is transferred from a Fourier, Laurent, Smith, or LRC
object. The probability estimate and its elementary exponential-moment lemma
are proved here, so no external permutation-concentration theorem is invoked.

## 2. Model and main theorem

Let G be any simple bipartite 2-regular graph on n vertices on each shore.
Independently label its shores by uniform bijections rho,sigma to
{0,...,n-1}. Its image B has exactly two points in every row and column.
All collinearity is over the integers. Let X count every nonaxis collinear
triple of B, and define

```text
D_d = {(r,c): c-r=d},       L_d=n-|d|,       1-n<=d<=n-1,
Y_d = |B intersect D_d|,
f(y) = (y-2)_+,
F = sum_d f(Y_d),           mu_G=E_G F,
alpha = 1-5/e^2 = 0.3233235838169366... .
```

Then F=0 if and only if no D_d contains three points. In particular
`{X=0} subset {F=0}`; the converse is false.

**Theorem.** Uniformly over every such skeleton G, for n>=4,

```text
2n/15 <= mu_G,
alpha*n-17 <= mu_G <= alpha*n+16,                         (1)

P_G(X=0) <= P_G(F=0)
         <= exp(-mu_G^2/[16(n-1)])
         <= exp(-b_n^2/[16(n-1)]),
b_n=max(2n/15, alpha*n-17).                              (2)
```

Consequently, at every n>=4,

```text
P_G(X=0) <= exp(-n^2/[900(n-1)]) <= exp(-n/900),           (3)
```

and

```text
limsup_(n->infinity) (1/n) log sup_G P_G(X=0)
 <= -(1-5/e^2)^2/16
 = -0.0065336337407642265... .                           (4)
```

The exact mean in Section 4 may be substituted into (2); the numerical rate
in (4) is the rate supplied by this defect and conditional-range argument,
not an optimality claim. The uniform mean bounds and the probability bounds
using b_n, (3), or (4) remain valid for mixtures over cycle types, including
the uniform distribution on all simple saturated boards. The bound using the
particular mu_G in (2) is conditional on that fixed skeleton; no concentration
claim around an arbitrary mixture's mean is inferred. Small positive success
probabilities, or successful boards at arbitrarily large n, are consistent
with this theorem.

## 3. The sharp swap bound and the concentration constant

Swapping two row labels moves at most four points: remove their old locations,
then insert their new locations. Each removal decreases F by zero or one,
and each insertion increases F by zero or one. The total decrease and total
increase therefore each lie in [0,4]; their difference has absolute value
at most four. The same argument applies to a column-label transposition.
This is a net-change bound of four, not eight.

The value four is attained already at n=5. Take the board with row i occupied
at columns i and i+1 modulo 5. Its defect is five. Swapping rows 0 and 2
gives defect one. Thus a smaller universal transposition bound cannot be
inserted without additional structure.

Expose the first n-1 values of rho, then the first n-1 values of sigma.
The final value of each permutation is determined. Let M_j be the conditional
expectation of F after j revealed labels; there are 2(n-1) nontrivial steps.
At a step, compare two possible labels a,b for the currently revealed vertex.
Swapping labels a,b in a completion gives a measure-preserving bijection
between the two uniform completion sets. All earlier reveals stay fixed,
and corresponding completed boards differ by one row or column transposition.
The unrevealed other permutation can be coupled identically. Thus the two
conditional expectations differ by at most four. Conditional on the past,
the possible values of M_j-M_(j-1) lie in an interval of **length** at most
four and have mean zero.

For completeness, if a mean-zero random variable W lies in an interval of
length h, then

```text
E exp(tW) <= exp(t^2 h^2/8).                             (5)
```

Indeed let K(t)=log E exp(tW). Its second derivative is the variance under
the exponentially tilted distribution. Subtracting the midpoint of the same
interval bounds this variance by h^2/4. Since K(0)=K'(0)=0, integration proves
(5), for either sign of t. Apply this conditional lemma successively with
h=4. For s>=0,

```text
E exp(-s(F-mu_G)) <= exp(4(n-1)s^2).
```

The exponential Markov inequality at F=0 and the choice
s=mu_G/[8(n-1)] give (2). Using 2n reveals would instead give the valid but
slightly weaker denominator 16n. Replacing conditional interval length four
by the enclosing symmetric interval [-4,4] would unnecessarily lose another
factor of four in the exponent.

## 4. Exact finite mean and alternating lower bounds

Write m_k(G) for the number of unordered k-edge matchings in G. For a grid
partial matching T of L cells and Y=|B intersect T|, inverse labels are
uniform injections on both shores, so for 0<=k<=n,

```text
E(Y)_k = (L)_k k! m_k(G)/(n)_k^2.                       (6)
```

This is non-induced containment; no further automorphism factor appears.
For k>L the numerator is zero. For k>n define E(Y)_k=0 directly, without
evaluating the quotient in (6).

For every nonnegative integer y,

```text
f(y)=sum_(k=3)^y (-1)^(k+1)(k-2) binom(y,k).             (7)
```

One derivation is f(y)=y-2+2*1_(y=0)+1_(y=1), followed by the two binomial
expansions of the indicator functions. More precisely, for y>=2 the partial
sum through K>=3 minus f(y) equals

```text
(-1)^(K+1) * ((K-2)y+2)/K * binom(y-2,K-1).             (8)
```

Use the convention that a binomial coefficient is zero if its lower index
exceeds its nonnegative upper index. For y=0,1 both sides of the required
truncation inequalities vanish directly. Thus every even truncation in (7)
is a lower bound, every odd truncation an upper bound.

The diagonal lengths consist of n once and every 1,...,n-1 twice; hence

```text
sum_d binom(L_d,k)=binom(n,k)+2binom(n,k+1).
```

Combining with (6)-(7) proves the exact all-profile mean

```text
mu_G = sum_(k=3)^n
 (-1)^(k+1)(k-2)(2n-k+1)/(k+1) * m_k(G)/(n)_k.          (9)
```

An even truncation of (9) is a rigorous finite lower bound. These formulas
require no asymptotic approximation or enumeration of label permutations.
If G has cycle components C_(2a), its matching polynomial is the product of

```text
R_(2a)(t) = sum_(k=0)^a
 [binom(2a-k,k)+binom(2a-k-1,k-1)] t^k.                 (10)
```

The second binomial is zero at k=0. Splitting matchings according to whether
they use the closing edge of a path proves (10). Polynomial multiplication
therefore evaluates (9) using only the cycle profile.

For n>=4 let c4 be the number of C4 components. The exact first coefficients
needed for the fourth-order lower bound are

```text
m_3 = 2n(n-2)(2n-5)/3,
m_4 = n(2n-5)(2n-6)(2n-7)/12 + c4.                    (11)
```

For example these follow directly from (10): with formal roots
u,v of z^2-z-t=0, the cycle polynomial is u^(2a)+v^(2a).
Since u=1+O(t), v=-t+O(t^2), the product over cycles agrees through degree
four with u^(2n)+c4*t^4; all other component corrections start at degree six.
The corresponding coefficients of the single cycle C_(2n), using (10),
give (11). This also explains the positive sign of the C4 correction.

The pointwise inequality f(y)>=binom(y,3)-2binom(y,4) now gives

```text
mu_G >= B4(n,c4)
 := (2n-5)(n^2+5n-11)/[15(n-1)(n-2)]
    -2(2n-3)c4/[5n(n-1)(n-2)(n-3)].                   (12)
```

Since c4<=n/2, put t=n-4>=0. Direct subtraction yields

```text
B4(n,n/2)-2n/15
 = [11t^3+48t^2+58t+12]/[15(n-1)(n-2)(n-3)] > 0.      (13)
```

This proves the first inequality in (1) and the finite exponential bound
(3). No cycle type is assumed to minimize the full mean. Indeed at n=4 the
means for 2C4 and C8 are 2/3 and 5/6, whereas at n=5 the means for C4+C6
and C10 are 29/25 and 17/15. One cannot infer a universal mean ordering
from a single matching coefficient.

## 5. Uniform approximation with an explicit finite error

This section proves the needed uniformity directly, avoiding any unstated
uniform-integrability step in a Poisson heuristic.

Fix a target matching of length L, set theta=L/n and lambda=2theta, and put

```text
delta_k=lambda^k-E(Y)_k.
```

For all n>=2 and k>=2,

```text
0 <= delta_k
 <= 2^k k(k-1) [theta^(k-1)/(2n)+theta^k/(4(n-1))].     (14)
```

Here is an elementary proof. Choose k ordered distinct rows of G uniformly,
then choose one of the two incident edges independently at each row. Let q_k
be the probability that all chosen columns are distinct. For k<=n,

```text
k! m_k(G) = (n)_k 2^k q_k,
E(Y)_k = 2^k (L)_k/(n)_k * q_k.                       (15)
```

For any selected pair of distinct rows the probability of a column collision
is exactly 1/[2(n-1)]: there are 2n ordered edge pairs sharing a column and
4(n)_2 equally likely ordered choices. A union bound gives

```text
0<=1-q_k<=k(k-1)/[4(n-1)].                            (16)
```

Also (L)_k/(n)_k<=theta^k. To bound its deficit, count ordered k-tuples
drawn from L positions with replacement. A specified collision pair gives
at most L^(k-1) such tuples, so

```text
0<=theta^k-(L)_k/(n)_k
 <=theta^k-(L)_k/n^k
 <=binom(k,2)*theta^(k-1)/n.                          (17)
```

Equations (15)-(17) imply (14) when k<=n. For k>n the factorial moment is
zero, and the first term alone in (14) bounds lambda^k because
k(k-1)/(2n)>=1 and theta<=1. The case theta=0 is direct. In particular
`E(Y)_k<=lambda^k<=2^k` for every k.

For a Poisson variable of mean lambda, (7) has expectation

```text
g(lambda)=lambda-2+(lambda+2)e^(-lambda).              (18)
```

All relevant factorial series are absolutely bounded by
sum_(k>=3)(k-2)2^k/k!, uniformly in G,n,L. Taking expectations and differences
is therefore justified. Odd k in (7) carry positive signs, even k negative.
Using (14) separately on those two classes gives the stronger one-sided bounds

```text
-cosh(lambda) H_n(lambda)
 <= E f(Y)-g(lambda)
 <= sinh(lambda) H_n(lambda),
H_n(lambda)=lambda^2/n+lambda^3/[4(n-1)].              (19)
```

For example, writing j=k-2, the odd sum uses
sum_(j odd) j lambda^j/j!=lambda*cosh(lambda); the even sum uses
lambda*sinh(lambda). This establishes (19) for the complete series,
including the terms k>n.

The discrete sum over the diagonal lengths is exactly n times the trapezoid
rule on [0,2] with mesh 2/n. Since g(0)=0, g is increasing, and
g''(lambda)=lambda*e^(-lambda)>=0,

```text
n*alpha <= sum_d g(2L_d/n) <= n*alpha+4e^(-2),
alpha=integral_0^2 g(lambda) d lambda=1-5e^(-2).        (20)
```

For any increasing h>=0 with h(0)=0, the same grid satisfies
sum_d h(2L_d/n)<=n*integral_0^2 h+h(2). Apply this to the four terms
in (19). Write S=sinh(2), C=cosh(2). The elementary integrals are

| integrand | integral from 0 to 2 |
|---|---:|
|lambda^2 cosh(lambda)|6S-4C|
|lambda^3 cosh(lambda)|20S-18C+6|
|lambda^2 sinh(lambda)|6C-4S-2|
|lambda^3 sinh(lambda)|20C-18S|

Consequently the lower error is at most

```text
C_n^-=(6S-4C)+4C/n
       +n(20S-18C+6)/[4(n-1)]+2C/(n-1),              (21)
```

and the upper error, in addition to 4e^(-2) from (20), is at most

```text
C_n^+=(6C-4S-2)+4S/n
       +n(20C-18S)/[4(n-1)]+2S/(n-1).                (22)
```

These are decreasing functions of n>=4; their maximum is at n=4.
There,

```text
C_4^-=(38S-25C+6)/3 =16.5886010737... <17,
C_4^++4e^(-2)=(38C-25S-6)/3+4e^(-2)
                             =15.9719831546... <16.  (23)
```

The strict rational enclosures follow, for instance, from the elementary
Taylor bound 7<e^2<37/5, as checked by rational arithmetic in the source.
This proves all of (1), including uniform `mu_G=alpha*n+O(1)`. Substitution
in the concentration inequality proves (4). For a sharper finite bound one
can use (12), or (20)-(21), rather than their rounded constants.

## 6. Exact controls, hostiles, and reproduction

The standard-library
[source](../../04-computation/overnight10_20260906_no3line_defect.py) and
[output](overnight10_20260906_no3line_defect.out)
import no repository or producer code. There are no floating-point tests
or assert statements. Both ordinary and optimized runs retain every gate.

The complete labelled-board universe at n=2,3,4,5 has respectively
1,6,90,2040 boards, totalling **2,137**. Every row and column has degree two;
each is generated exactly once by its row-neighbor pairs. Literal components
classify the cycle type, and the independently known orbit count is
`(n!)^2/product_a (2a)^(c_a)c_a!`. Within each type, uniform independent
shore labels give exactly the uniform board law used by the averages.

The source reconstructs the exact mean and every individual diagonal
factorial moment, checks all **41,918** row/column transpositions, and confirms
the sharp n5 change 5->1. Two complete n4 permutation reveal trees each
contain **425** conditional nodes; their exact rational conditional ranges
are at most two, within the universal bound four. The proof does not assume
that the small-tree maximum persists at larger n.

All **384** cycle profiles with shore size 2..18 verify the rook coefficients,
finite fourth-order bound, and all-k factorial deficit at selected lengths
including zero and full length, with k extending three beyond n. Pointwise
truncation controls cover y=0..80 and K=3..16. Rational Taylor enclosures
verify the rounded mean constants. These finite controls support the separate
analytic proofs; they are not the basis of any all-n extrapolation.

The success-event firewall is visible literally: among the 72 C8 boards at
n=4, 25 have F=0, but only two have X=0. Thus the selected direction leaves
substantial information out even at small order. Counting zero defects alone
is never used as an extremal no-three-in-line certificate.

```text
python 04-computation/overnight10_20260906_no3line_defect.py
python -O 04-computation/overnight10_20260906_no3line_defect.py

PASS 55520 always-active gates
normal and optimized transcripts: byte-identical, LF

source SHA256:
74ce43c6435da3c4662b5397c9e89ef86874452fcb1b445521a57b55816399d6
output SHA256:
b729b941bb7e4b8b2decb6486a1713554c59915ee31717982847495da7faa295
```

The missing coordinate in the earlier variance-only route was a bounded
local defect with the same selected zero event, not a stronger independence
assumption on geometric triples. The exact second-order moment work remains
useful structural information. Optimizing the exponential rate or exploiting
more directions without proportionately increasing the swap range is a
separate open task; this report makes no such optimization claim.

**Filing:** root integrated these audited artifacts in the tenth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
