# Four original windows repair the next positive-amplitude obstruction

**Status: PROVED fixed-row representation; FINITE-EXACT and INDEPENDENTLY
AUDITED.** The [independent referee](continuing2_20260906_laurent_sparse_amplitude_audit.md)
accepts the proof and reconstructs its coefficients and determinant signs
with a different exact elimination path. This is a fixed-row representation theorem for the actual support
**(-27,1,15)**. Four positive real algebraic coefficients represent its full
carried doubled response at all four positive first phases, using the original
coupled derivative windows and amplitude powers at most three. This does not
give a general Laurent noncancellation theorem or a uniform certificate for
the finite-phase anchored model.

## 1. Inheritance and the precise gain

The nearest proved mechanism is the
[minimal coupled-window cone](creative_20260906_laurent_bridge.md),
whose sign input is **THM-4440**, `01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md`.
The literal source map is the proved
[A2B3 midpoint transport](overnight7_20260906_laurent_midpoint_transport.md)
and [alpha completion](overnight8_20260906_alpha_completion.md).
No new external theorem is invoked.

The canonical hostile is the audited
[h=3 all-integer-phase cone separator](continuing1_20260906_laurent_cone_separator.md).
Its repaired form allows square-root amplitude weights and real algebraic
coefficients on the selected positive branches; the
[h=3 positive repair](continuing1_20260906_laurent_real_power_boundary.md)
is independently audited. The next corrected near miss is the actual
[h=4 minimal-amplitude obstruction](continuing1_20260906_laurent_amplitude_h4_boundary.md):
the unique polynomial A(z) of degree at most three with
A(sqrt(r_i))=T(r_i)/J_0(r_i) has

    -344095/10^6 <= [z^2]A <= -344094/10^6 < 0.                 (1)

That refutes coefficient positivity of the single-window recipe. It does
not refute a mixed-window representation. The least-used sidecar is the
joint root-evaluation vector of several differently differentiated windows.
The board here consists of the complete carried row, the original Euler
scale, the finite first-row quotient, positive amplitude multiplication,
mixed derivative orders, and the coefficient field/branch restriction.

The gain below keeps the maximum amplitude degree three and changes the
window choices. It is an existence certificate with four atoms, not a claim
that four atoms or amplitude degree three are minimal.

## 2. Actual source, target, and normalization

Set h=4, x=1, q=5, g=14, k=9, n=24. All three support exponents are congruent
to one modulo 14, and counts (1,12,1) give a return at mass 14. Thus the first
support-return mass is exactly 14. The full first and doubled rows are

    P_raw(t) = sum_(j=0)^4 14! t^j /[(1+j)!(12-3j)!(1+2j)!],
    Q_raw(t) = sum_(j=-1)^8 28! t^j /[(2+j)!(24-3j)!(2+2j)!].  (2)

The lower carry j=-1 in Q_raw is present. For the Laurent polynomial with
coefficients (1,1,gamma), the literal constant terms are gamma P_raw(gamma^2)
and gamma^2 Q_raw(gamma^2). Keep this same coefficient normalization and put
t=-s, s>0, z=sqrt(s)>0; thus gamma=i z in this representative normalization.

The complete reversed alternating beta kernel and carrier are

    B(v) = v^5-13v^4+55v^3-84v^2+35v-1,
    K_B(s,u) = s^5 B(u^2/s),
    H(s,u) = (1+u)^14 K_B(s,u).

The displayed B has positive roots: the inherited unsigned path polynomial
has negative roots, and the reversal/alternating substitution changes their
sign. Hence K_B and H are real-rooted in u for s>0. Direct extraction gives

    [u^9]H = -s P_raw(-s),
    p(s) = P_raw(-s)/2002 = s^4-60s^3+84s^2-10s+1/11.

Let r_1<r_2<r_3<r_4 be the four positive roots of p and z_i=sqrt(r_i)>0.
Define the original basic windows and actual response by

    J_r(s) = s^-1 [u^(18-2r)] (D_u^r H(s,u))^2,  0<=r<9,
    T(s) = s Q_raw(-s).                                      (3)

All divisions by s are exact as polynomials. At every r_i the genuine
interior coefficient [u^9]H vanishes. The inherited window theorem gives
J_r(r_i)<0; the certificate independently verifies all 36 strict signs.
In particular J_0=s W_raw(-s) is the complete alpha-completed hit response.

The source is the literal charge fibre at masses 14 and 28. The target is
membership of its four-component response vector in a positive cone. The
map preserves every coefficient, the lower carry, the original phase, and
all four first-root values. Reducing a polynomial modulo p discards only
multiples of the actual first equation. Choosing positive square roots and
allowing real algebraic common coefficients loses Galois invariance across
the eight z-conjugates; Section 5 records that loss explicitly.

## 3. The four-atom identity

Write (9)_r=9!/(9-r)! and use the positive normalizations
G_r=J_r/(9)_r^2. These factors are stated explicitly so that they cannot
alter the original response normalization. Choose atoms

    (z_power, derivative_order) = (0,0), (1,7), (1,8), (3,2).

For i=1,...,4 define the exact real algebraic matrix and right side

    A_(i,l) = z_i^(d_l) J_(r_l)(r_i) /[(9)_(r_l)^2 J_0(r_i)],
    b_i = T(r_i)/J_0(r_i).

The first column of A is exactly one. Let A[l<-b] replace column l by b and
define beta_l=det(A[l<-b])/det(A), for l=0,1,2,3. The exact certificate in
Section 4 proves det(A)<0 and every replacement determinant negative.
Therefore all four beta_l are strictly positive. Cramer's rule gives the
same four coefficients at all four positive phases:

    T(r_i) = beta_0 J_0(r_i)
           + beta_1 z_i J_7(r_i)/32920473600
           + beta_2 z_i J_8(r_i)/131681894400
           + beta_3 z_i^3 J_2(r_i)/5184,     i=1,...,4.        (4)

The coefficients are exactly the algebraic determinant ratios above, not
rounded decimals. Their certified bounds are

| Coefficient | Lower bound | Upper bound |
|---|---:|---:|
| beta_0 | 984869/10^6 | 984870/10^6 |
| beta_1 | 421733/10^6 | 421734/10^6 |
| beta_2 | 4973466/10^6 | 4973467/10^6 |
| beta_3 | 5162/10^6 | 5163/10^6 |

Every summand in (4) is negative. Thus (4) is a new positive sign
certificate for this fixed actual row. The fact T(r_i)<0 was already
verified in the inherited work; that underlying sign is not claimed as a
new discovery. The new mechanism is a positive mixed-window representation
despite the negative coefficient in (1).

## 4. Exact certificate and independent reconstruction paths

These initial rational brackets are positive, disjoint, and have opposite
p endpoint signs. Four sign changes exhaust its degree, so they isolate all
four simple roots without relying on floating-point root approximations.

| Root | Left endpoint | Right endpoint |
|---|---:|---:|
| r_1 | 8419/849544 | 11993/1210189 |
| r_2 | 259526/2155711 | 291249/2419213 |
| r_3 | 11376199/8744207 | 1124341/864214 |
| r_4 | 120947883/2065060 | 91906768/1569213 |

The source makes 43 exact sign-preserving bisections in every bracket.
For each resulting rational interval [a,b], integer square roots give
rational endpoints enclosing sqrt(a),sqrt(b), checked by squaring. Rational
Horner arithmetic bounds each remainder J_r mod p and T mod p, then each
matrix entry. The determinant enclosure expands all 24 signed permutations;
it does not assume independent choices of the correlated root coordinates.
Interval inclusion is sufficient even though it deliberately forgets those
correlations. The following outer bounds are all rational:

| Determinant | Lower bound | Upper bound |
|---|---:|---:|
| det(A) | -371065/10^6 | -371064/10^6 |
| det(A[0<-b]) | -365451/10^6 | -365450/10^6 |
| det(A[1<-b]) | -156491/10^6 | -156490/10^6 |
| det(A[2<-b]) | -1845480/10^6 | -1845479/10^6 |
| det(A[3<-b]) | -1916/10^6 | -1915/10^6 |

The [exact JSON certificate](../../04-computation/continuing2_20260906_laurent_sparse_amplitude_certificate.json)
contains the complete polynomials, remainders, final rational root and square
root brackets, every matrix interval, determinant/numerator intervals, and
the sharper coefficient intervals from which these tables follow.

The [standalone producer](../../04-computation/continuing2_20260906_laurent_sparse_amplitude.py)
uses only Python's standard library. It enumerates both literal charge fibres
to reconstruct every coefficient of (2), builds every differentiated window
from H, and independently checks J_0 by complete beta-path convolution. Its
finite universe is this single support, all first/doubled channels, all nine
basic windows at all four first roots, and the five four-by-four determinants.
Numerical cone optimization only located the candidate atoms and is not a
proof input. No numerical infeasibility result is promoted.

## 5. Coefficient-field and full-eight-root boundaries

The primitive numerator of p is 11s^4-660s^3+924s^2-110s+1. Its monic
reduction modulo 3 is s^4+2s+2. It has no linear factor and none of the three
monic irreducible quadratic factors over F_3. The verifier exhausts this
small factor universe, proving p irreducible over Q. In Q(r_i),
Norm(r_i)=1/11 is not a rational square, so r_i is not a square in that field.
Consequently z_i has degree eight over Q and p(z^2) is irreducible.

Suppose a nonnegative real combination of atoms z^d J_r(z^2), with integer
d, agreed with T(z^2) on all eight roots. Subtracting its values at z_i and
-z_i gives twice its odd-power part. Every nonzero odd atom is strictly
negative at z_i>0, because its coefficient is nonnegative and J_r(r_i)<0.
Thus every odd coefficient would have to vanish. The positive odd terms in
(4) make (4) impossible as an all-eight-root identity.

Moreover the four coefficients in (4) cannot all be rational. A rational
Laurent-polynomial identity holding at even one z_i would, after clearing
an even power of z, be divisible by the irreducible p(z^2) and hence hold at
all eight roots, contradicting the odd-part argument. The same reasoning
excludes any rational positive representation containing a nonzero odd
atom. It does **not** exclude a possible all-even/integer-s certificate for
this h=4 row. Unlike the h=3 separator theorem, no all-integer h=4 no-go is
claimed here. No uniform formula in h is proved.

## 6. Connection to the incoming quadratic anchor

The incoming proved
[quadratic-anchor restriction](open_frontier_sep06_quadratic_anchor.md)
keeps e_1=13 and e_2=55, requires both specified contiguous interlacers,
and proves e_4 is bounded away from zero on its compact model class. It
thereby excludes every escaping original-phase obstruction and proves a
qualitative uniform negative tail. That result supplies no numerical threshold. The
[effective refinement](continuing2_20260906_effective_anchor.md) now proves
e4>1/100 and negativity at every original phase s>=118163898523. The
finite-phase sign problem remains open.

Here every beta coefficient is fixed at its genuine value
(e_1,e_2,e_3,e_4,e_5)=(13,55,84,35,1), the original phase is retained, and
the exact certificate concerns all four finite phases of that one row.
These are complementary uses of retained data: the incoming theorem
controls an entire model boundary through two coefficient anchors and both
interlacers; (4) supplies a finite-dimensional positive representation at
one actual coefficient point. Algebraic coefficients chosen from this
point's ordered roots provide no uniform finite-phase model certificate.
The next bounded structural question is which mixed-window selection rule
can make the relevant determinants have a uniform sign as the genuine
coefficient parameters vary.

## 7. Reproduction and frozen artifacts

From the repository root:

```text
python -B 04-computation/continuing2_20260906_laurent_sparse_amplitude.py --certificate 04-computation/continuing2_20260906_laurent_sparse_amplitude_certificate.json
python -B -O 04-computation/continuing2_20260906_laurent_sparse_amplitude.py
```

Both runs pass **439 always-active gates**, with byte-identical LF
[output](continuing2_20260906_laurent_sparse_amplitude.out).
No optimization, numerical root solver, or external algebra library enters
the verifier. The old h=3 no-go, its positive-branch repair, and the h=4
single-window obstruction remain unchanged.

```text
source      2afc7df14b92e70a700a9a06a17b34ee2331eabde89902dcf408f597dccb9e28
output      7b95d9102c18309ce245d9017f3456b0f75447117f1ec6ac363393c63f183777
certificate df8471af42220b5b26f7e4ff78dd9611d7999c0578e36603a6a7a33ba3c0afa2
```
