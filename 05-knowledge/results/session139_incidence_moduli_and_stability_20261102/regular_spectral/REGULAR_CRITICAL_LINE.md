# A finite functional equation and critical circle do not determine weighted Hamiltonian stability

Status: PROVED by explicit linear algebra and incidence identities below;
exact coefficientwise verifier passes in normal and optimized Python.
Baseline source commit: `291b4d9757842c02aa451beb65acfade2848609b`.
This is a finite-matrix analogue, **not** a theorem about Riemann zeta zeros,
RH, GRH, or any special analytic L-function. No priority claim is made.

## 1. Every regular weighted tournament has this finite spectral symmetry

**Literature status.** The binary regular-tournament eigenvalue-line fact is
known: see Eschenbach, Hall, Hemasinha, Kirkland, Li, Shader, Stuart and
Weaver, [*Properties of Tournaments Among Well-Matched Players*, American
Mathematical Monthly 107 (December2000), pp.881ff](https://galton.uchicago.edu/~lekheng/meetings/mathofranking/ref/shader.pdf),
printed pp.886--887 (including Theorem5's regularity/normality equivalence).
The weighted statement and Möbius identity used below are derived directly;
no novelty claim or historical open-problem status is imported from that paper.

Let A be a real weighted tournament of order n>=2:
Aii=0, 0<=Aij<=1, A+A^T=J-I. Put K=A-A^T. Assume regularity
K1=0, equivalently A1=d1 with d=(n-1)/2. Then A^T1=d1 as well, and
R^n splits orthogonally into span(1) and its invariant complement E=1^perp.
On E,

    A|E=(-I+K|E)/2.

Since K|E is real skew-symmetric, every eigenvalue on E has the form

    lambda=-1/2+i mu/2,    mu real.                       (1)

Thus **all nonprincipal eigenvalues lie on Re(lambda)=-1/2**. This is
automatic for regular weighted tournaments; it is not a difficult
zero-location hypothesis. Real coefficients make their multiset invariant
under lambda -> conjugate(lambda)=-1-lambda.

Write

    P_A(z)=det(I-zA)=(1-dz) R_A(z),
    R_A(z)=product_(lambda on E)(1-z lambda).             (2)

There is no zero eigenvalue on E, so R_A has degree n-1. Its roots are the
poles z=1/lambda of the nonprincipal dynamical zeta 1/R_A(z). By (1),

    |z+1|=|(lambda+1)/lambda|=1.                         (3)

Consequently **all reduced-zeta poles lie on the circle centered at -1
with radius 1**. They are poles of a finite rational function, not zeros of
the classical Riemann zeta function.

The same involution gives the exact polynomial functional identity

    R_A(z)=(1+z)^(n-1) R_A(-z/(1+z)).                   (4)

Indeed the right side is product_(lambda on E)(1+z+z lambda), and
replacing lambda by -1-lambda in the multiset yields R_A(z).
This calculation is first made away from z=-1; clearing denominators
proves the identity as polynomials everywhere. The fractional-linear map
z -> -z/(1+z) is an involution. Equations (1), (3), and (4) are exact
finite critical-line/circle/functional-equation counterparts, but supply no
map to the analytic variable or zeros of Riemann zeta.

Finally the all-one marked resolvent and all-one walk moments are

    1^T(I-zA)^(-1)1=n/(1-dz),
    1^T A^k 1=n d^k.                                    (5)

They depend only on n throughout the regular class.

## 2. A whole regular isospectral circle at five vertices

Use vertices 1,2,3,4,5. For real t and x=(x1,x2,x3), define the skew
matrix K by

    Ki4=xi,  Ki5=-xi (i=1,2,3),  K45=t,
    Kji=-Kij; all remaining entries zero.

Impose

    x1+x2+x3=t,    x1^2+x2^2+x3^2=t^2,    |t|<=1.        (6)

Then K1=0. Also |xi|<=|t|<=1, so A(x,t)=(J-I+K)/2 is a regular
weighted tournament. At t=0 it is the fully fair tournament. For fixed
nonzero t the allowed x form a circle: the intersection of a plane and a
sphere, centered at (t/3,t/3,t/3).

K acts only on the span of (x,0,0), e4, e5. Its nonzero eigenvalues are

    ±i sqrt(2||x||^2+t^2)=±i sqrt(3) t

as an unordered pair; its other three eigenvalues are zero. This also
follows by multiplying K^3=-(2||x||^2+t^2)K. Since one zero direction
is 1, the adjacency spectrum is

    2,  -1/2,  -1/2,  (-1+i sqrt(3)t)/2,
                         (-1-i sqrt(3)t)/2.              (7)

Every matrix on the circle has the same full spectrum, the same all-one
marked resolvent, (1)--(5), and common determinant polynomial

    P_t(z)=(1-2z)(1+z/2)^2
               *[(1+z/2)^2+3t^2 z^2/4].                 (8)

In particular

    u=det(I+A)=3(1+3t^2)/16,
    v=-det A=-(1+3t^2)/8<0.                              (9)

## 3. The incidence coordinate varies despite all those symmetries

Use the Session138 quantities V=sum_edges Aij Aji, tau=total directed
triangle weight, and M23=sum_edges (Aij Aji) tau(A[e^c]). Directly,

    V=(10-3t^2)/4,    tau=(10+3t^2)/4=5-V.               (10)

For completeness, a triangle on three vertices of {1,2,3} has weight
1/4; a triangle on i,j,4 or i,j,5 has weight (1-xi xj)/4; and a
triangle on i,4,5 has weight (1+2t xi+xi^2)/4. These statements follow
from 4 tau_ijk=1+Kij Kjk+Kjk Kki+Kki Kij.

The ten complementary edge/triangle terms split into three types:

    16 M23 = sum_i(1+2t xi+xi^2) + (1-t^2)
              +2 sum_i(1-xi^2)(1-xj xk), {i,j,k}={1,2,3}.

By (6), sum_(i<j)xi xj=0. Also
sum_i xi^2 xj xk=(x1x2x3)sum_i xi=t x1x2x3. Therefore

    M23=5/8+t x1x2x3/8.                                 (11)

The weighted odd-cycle and determinant expansions at order five give

    H=1+2tau+2M23-2[z^5]P_t,
    [z^5]P_t=v.

Combining with (9)--(11),

    H(A(x,t))=15/2+9t^2/4+t x1x2x3/4,                   (12)
    Delta(A(x,t))=115/16+21t^2/16+t x1x2x3/4.            (13)

Thus neither the entire spectrum nor the critical-line/circle condition,
finite functional equation, regularity, and all-one marked resolvent
*together* determine H or Delta. The lost information is vertex-disjoint
edge/triangle incidence, not a missing zero-location condition.

## 4. Exact oscillation width on this isospectral circle

For t!=0 write si=xi/t. Then sum si=1 and sum si^2=1, so their pair sum
is zero. Put r=s1s2s3. The si are the three real roots of

    f(y)=y^3-y^2-r.

Its discriminant is

    disc(f)=-r(4+27r)>=0,

and hence -4/27<=r<=0. The endpoints occur exactly at permutations of

    r=0:      (s1,s2,s3)=(1,0,0),
    r=-4/27:  (s1,s2,s3)=(-1/3,2/3,2/3).

The circle is connected and r is continuous, so every intermediate value
is attained. At t=0 the statement is immediate without normalization.
Consequently the **exact** ranges, for fixed t, are

    H in [15/2+9t^2/4-t^4/27, 15/2+9t^2/4],            (14)
    Delta in [115/16+21t^2/16-t^4/27,
                                      115/16+21t^2/16]. (15)

This classifies the incidence oscillation on the specific 3+2-support
regular rank-two circle (6), **not the full five-vertex isospectral fibre**.

A rational parameterization of the normalized circle, apart from its
limiting point, is

    (s1,s2,s3)=(-q,1+q,q(1+q))/(1+q+q^2).

The product bound has the elementary exact certificate

    4(1+q+q^2)^3-27q^2(1+q)^2
        =[(q-1)(q+2)(2q+1)]^2>=0.                       (16)

This provides a second route to the sharp endpoints without relying on a
numerical root calculation.

## 5. Explicit conjugate extremizers

Let A_t be the upper endpoint x=(t,0,0), so K has the directed triangle
1->4->5->1 of skew-weight t and is zero elsewhere. Set

    w=(1,1,-2,0,0)^T,  Q=I-ww^T/3,  B_t=Q A_t Q^T.

Q is orthogonal and fixes 1. Its conjugate has x=(2t/3,-t/3,2t/3),
the lower endpoint of the circle. Thus

    H(A_t)=15/2+9t^2/4,
    H(B_t)=15/2+9t^2/4-t^4/27,
    Delta(B_t)-Delta(A_t)=-t^4/27.                      (17)

At t=1/4 both are strictly interior weighted tournaments, and

| invariant | A_(1/4) | B_(1/4) |
|---|---:|---:|
| V | 157/64 | 157/64 |
| tau | 163/64 | 163/64 |
| u | 57/256 | 57/256 |
| v | -19/128 | -19/128 |
| H | 489/64 | 52811/6912 |
| Delta | 1861/256 | 25123/3456 |

Their H and Delta difference is -1/6912. At t=1 the respective H values
are 39/4 and 1049/108, with difference -1/27.

For any modulus m>=2 coprime to 6, all the rational entries and polynomial
Hamilton sums at t=1/4 reduce canonically from Z[1/6] to Z/mZ. The common
characteristic polynomial remains common, while the H difference is the
unit -1/6912. Thus this obstruction survives modulo5, modulo11, and every
other such modulus. This is a statement about polynomial evaluation over
rings; positivity, absolute value, and real Delta inequalities do not
transport under reduction. Moduli divisible by2 or3 require denominator
clearing with a retained scale, not illicit inverses.

## 6. What this recovers from the older RH-adjacent work

Source pins remain the baseline commit above:
- `07-reflections/riemann-zeta-tournament.md` (March17; August25 correction)
  proposed complement/functional-equation analogies but no RH equivalence.
- `01-canon/theorems/THM-924-walk-reciprocity.md` (July16) gives exact
  symmetrization of adjacency spectrum about Re(lambda)=-1/2.
- `01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md`
  identifies the finite dynamical zeta and distinguishes adjacency from
  skew normalization (October23 correction).
- `01-canon/theorems/THM-499-the-spectral-reframe-boundary-H-is-finer-than-the-spectrum.md`
  (June13) locates binary spectral blindness at six; weighted digons move
  the first boundary to five, as established in the first-checkpoint note.
- `05-knowledge/results/session138_stability_and_actual_contexts_20261102/stability/SHARP_LOW_ORDER_STABILITY.md`
  supplies the weighted cycle expansion and exact five-vertex inequalities.

**Typed connection.** Source: a regular weighted tournament and its finite
closed-walk zeta. Target: exact five-vertex H and Delta. Map: the principal
factor deletion and Möbius transformation (2)--(4). Preserved: full spectral
and marked-root data; regularity. Lost: M23's vertex-exclusion incidence.
Repair sidecar: M23 suffices at order five. Cheapest decisive test: (17) at
t=1/4. The critical circle is a finite linear-algebra consequence, not a
Riemann zero set. No implication for classical RH/GRH follows.

This also does **not** rule out spectral sufficient lower bounds or a
spectral component in a stability proof. It disproves exact reconstruction
from the specified spectral data. The family is not a sharp Delta/V
extremizer, and it does not itself close any finite-stability order.

## 7. Reproduction and validation

    python3 -B verify_regular_critical_line.py > verify_regular_critical_line.out
    python3 -B -O verify_regular_critical_line.py > verify_regular_critical_line.optimized.out
    cmp verify_regular_critical_line.out verify_regular_critical_line.optimized.out

The standard-library verifier works coefficientwise over Q[t], with
Fraction arithmetic and no assertions, floating constants, SymPy, numerical
roots, or interpolation. It separately implements permutation Hamilton sums
and subset-DP Hamilton sums; determinants use permutation expansion over
principal subsets. It checks conjugacy, row sums, coefficient bounds, full
polynomial equality, H/M23/Delta identities, the reduced functional equation,
and the sharp product square certificate. Six rational circle parameters
supply extra all-t controls; the all-circle proof is the displayed algebra.

Mutation flags `--mutate-conjugation`, `--mutate-difference`, and
`--mutate-functional-sign` must all return nonzero even under `python3 -O`.
Their exact rejection records are in `mutations.json`.
