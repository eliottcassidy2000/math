# Independent audit of the outward gap-one content theorem

Date: 20 October 2026.

**Verdict: PASS.** The proposed full-content bound

    gcd(|2^(p+r)-3^p|, B(w)) <= 22^(4 p^(2/3))

is proved for the precise primitive mechanical outward gap-one family stated
in `../OUTWARD_GAP_ONE_CONTENT.md`, conditional only on its stated standard
Minkowski second-theorem input and the already audited moving-anchor carry
identity. The present audit independently checks those inputs at the point of
use. No repair to the mathematical argument is required. The proof is not an
all-length integer-cycle exclusion and does not concern arbitrary disjoint
pairs.

The reviewer changed no frozen producer, theorem, or receipt. The source and
frozen receipts were rerun unchanged; new files are confined to this audit
directory. The inherited root gate was independently checked on actual words,
not merely on the reduced exponent equations.

## 1. Inheritance and audit board

The closest proved mechanism is the two-one overlap lattice-content theorem,
`collatz_two_swap_overlap_lattice_content_20261020/OVERLAPPING_TWO_ONE_CONTENT.md`.
Its external-prime valuation witness does not apply to the present smooth
constant -4. The relevant near miss is substituting a generic smooth-S-unit
noncancellation lemma for this legal-range theorem. The canonical hostile is
p=3, r=1, u=9, with common root 3/2. The least-used necessary sidecar is the
place-by-place normalized valuation together with the actual physical outer
zeros.

The audit board was:

1. Physical cyclic `0110` patterns, including the shared-zero q=3 boundary.
2. Full prime-power carry content versus a clock-forgetting resultant.
3. Rank of `(m,c)` projections versus rank of full lattice rows.
4. Normalized valuations, field degree, and archimedean norms.
5. Exact lattice height, coefficient growth, and constant-polynomial cases.

The applicable method card is “Pay algebraic rank, then audit the physical
action”: the rank calculation alone is insufficient without the legal physical
range. No new general-purpose method card is claimed from this single audit.

## 2. Exact physical scope and modular gate

The hypotheses are p>=3, 1<=r<p, gcd(p,r)=1, q=p+r, with

    max(0,p-2r) <= u < min(r,p-r).

Two consecutive source ones occupy adjacent physical positions and have
residues u+r, u. A left move at residue u+r requires u+r>=p-r; a right move
at residue u requires u<r. A gap-one pair also has 0<=u<p-r. Their conjunction
is exactly the displayed interval.

For every nonempty case, r>p/3 and r+u<=p-1. If r<=p/2, nonemptiness forces
p-2r<r. If r>p/2, the first inequality is immediate. The integer upper bound
on u gives the second inequality.

The only numerically admitted primitive q=3 boundary is p=2,r=1,u=0. Its
source is `011`. The outward candidate edges are 0 and 2; both target the same
zero. After either first move, the other edge is `11`, so the second swap is
not actual. This is not a disjoint path. The theorem expressly excludes it,
and the independent census requires four distinct physical positions.

In a nonempty legal case q>=5, the four positions of `0110` are distinct.
Changing the two disjoint edges gives `1001` in either order. This verifies
both existence and commutation. Original-cut crossings are cyclic positions,
not a separate exception.

For every h dividing |D|, h>1, 2 and 3 are units. A Bezout combination of the
clock exponents gives the unique z satisfying z^p=2 and z^r=3/2 modulo h. The
identity

    (z-1)(1+z+...+z^(p-1)) = z^p-1 = 1

makes z-1 a unit over the full ring, including nonreduced prime-power rings.
The moving-anchor gate reduces as follows, without dividing by a nonunit:

    Q = 1 + (z-1)(z^u-z^(u+r)/2)
      = 1 + z^u(z-1)(1-3/4)
      = (4+z^u(z-1))/4.

Thus h divides B(w) if and only if h divides 4+z^u(z-1). This holds for every
prime-power divisor h of D, not just for primes. Reconstructing z modulo the
whole |D| gives exactly the claimed equality of integer gcds. The independent
checker instead reconstructs z from z^p=2, z^q=3 using a Bezout relation for
(p,q), then verifies the complete gcd directly against a reverse-scan physical
carry. This is a distinct clock parametrization from the producer.

## 3. Cleared polynomials and all degeneracies

The congruence lattice has determinant p because the coefficient of m in
m-br-cu modulo p is 1. For H=max(|b|,|c|,|m|),

    |a| = |m-br-cu|/p <= H(1+r+u)/p <= H.

Negating a row to make c>=0 is valid because z and z-1 are units modulo the
full g. With e=a-b+2c, the proposed clearing uses only 2 and 3, which are units
modulo g. Therefore g divides every F_v(z), including all prime-power depth.

The apparent clearing roots are absent:

- At X=0, the only possible cancellation has m=0 and e=b=0. Then a=-2c
  and 0=c(u-2p), which forces c=0 in the legal range and gives the zero row.
  Other signs of m leave one nonzero constant term.
- At X=1, c>0 leaves -V. If c=0, cancellation forces a=b=m=0.
- A possible leading cancellation with m<0 has -m=c and U=V. Unique prime
  factorization gives b=0, a=-2c, and hence u=2p-1. It is outside the legal
  range. The sign condition on c can only further restrict this impossible
  case.
- Consequently the degree is zero exactly when m=c=0. Two independent full
  rows cannot both be on that b-axis. Nonzero constant polynomials must be
  allowed, not silently discarded.

The independent implementation creates the rational polynomial in QQ[X] and
then uses exact denominator clearing, rather than copying the producer's
integer binomial-array construction. Explicit constant and
independent-full-row/dependent-projection controls are retained in its receipt.

## 4. Common-root proof: the potentially dangerous steps

Assume two nonconstant F_i share x. The previous section gives x!=0,1, so
both original rational multiplicative equations hold in K=Q(x).

### 4.1 Dependent projection is not dependent full row

If (m_i,c_i) are dependent, a nonzero integer combination kills both m and c.
Independence of the full rows ensures its remaining (a,b) is nonzero. Since
ap+br=0 and gcd(p,r)=1, it is (rt,-pt), t!=0. The same combination of the
multiplicative equations yields

    1 = 2^(qt) 3^(-pt),

which is impossible. Negative combination coefficients are legitimate because
x and x-1 are nonzero. There is no hidden sign left: the combined c is zero.

### 4.2 Uniformity is over every prime above a rational prime

Suppose the projected 2-by-2 matrix is invertible. At each prime P of K over
ell, use v_P/e_P, so the rational number ell has valuation 1. The right sides
of both row equations have the same valuations at every P above ell. Solving
the same invertible rational linear system therefore makes the valuations of
x and x-1 equal at all these P. This is the missing premise that would make
an indiscriminate norm argument fail; it is genuinely present here.

For ell outside 2 and 3, the two right-hand valuations vanish, forcing both
unknown valuations to be zero. At 2 and 3 the row equations are

    m_i alpha_2 + c_i beta_2 = a_i-b_i+2c_i,
    m_i alpha_3 + c_i beta_3 = b_i.

Multiplying by p and q respectively, adding, and substituting the lattice
relation gives

    m_i(p alpha_2+q alpha_3-1)
      + c_i(p beta_2+q beta_3-(2p-u)) = 0.

Invertibility now proves both identities (4.1) in the candidate, without a
field-degree or ramification assumption.

### 4.3 The valuation pattern really is exhaustive

For a nonarchimedean valuation, the possibilities for x and x-1 obey:

- alpha<0 implies beta=alpha;
- alpha>0 implies beta=0;
- alpha=0 implies beta>=0.

If alpha_2<0, the first identity forces alpha_3>0, hence beta_3=0, while the
second identity would have negative left side and positive right side 2p-u.
The reversed case works with coefficients p and q exchanged. Thus neither
alpha is negative. Both positive would make both beta zero; both zero would
violate the first identity. Exactly one alpha is positive and the other is
zero. Its beta is zero; the other beta is nonnegative.

All finite valuations of x and x-1 are therefore nonnegative, so both are
algebraic integers. No unproved rationality of x, splitting of 2 or 3, or
normality of K is being used.

### 4.4 Degree and norm accounting

Let d=[K:Q]. For y=x-1, the exponent of ell in |Norm(y)| is

    sum_(P|ell) f_P v_P(y)
      = beta_ell sum_(P|ell) f_P e_P
      = d beta_ell.

The equality sum e_P f_P=d holds for every rational prime in a number field.
There is no missing residue degree or ramification factor. The exponents are
integers even though the displayed normalized beta values can be rational.
All other finite valuations vanish, so the candidate's exact norm formula
follows. The absolute value handles a possible negative norm.

### 4.5 Eliminating c survives torsion

For M_1=2/x^p, M_2=3/(2x^r), M_3=-4/[x^u(x-1)], the row relations become
M_1^a_i M_2^b_i M_3^c_i=1. Eliminating c gives exponents

    A = a_1 c_2-a_2 c_1,
    B = b_1 c_2-b_2 c_1.

They are not both zero because

    pA+rB = m_1 c_2-m_2 c_1 != 0.

This also covers a zero c-coordinate on one row. Both c-coordinates cannot
vanish in the independent-projection case.

If alpha_2>0, the exact values are alpha_2=1/p, alpha_3=0, beta_2=0,
beta_3=(2p-u)/q>1. At every prime over 3, v(M_1)=0 and v(M_2)=1, forcing
B=0 and A!=0. Thus M_1 is a root of unity. Under every embedding of K into C
it has absolute value 1, so every conjugate of x has modulus 2^(1/p). Hence

    |Norm(x-1)| <= (1+2^(1/p))^d <= (1+sqrt(2))^d < 3^d,

whereas the finite norm formula gives 3^(d beta_3)>3^d.

If alpha_3>0, the exact values are alpha_3=1/q, alpha_2=0, beta_3=0,
beta_2=2-u/p>q/p>4/3. At every prime over 2, v(M_1)=1 and v(M_2)=-1, forcing
A=B!=0. Therefore 3/x^q is a root of unity, and every conjugate of x has
modulus 3^(1/q). Consequently

    |Norm(x-1)| <= (1+3^(1/q))^d
      <= (1+cuberoot(3))^d < (5/2)^d < 2^(4d/3),

contrary to the finite norm formula. The two strict numerical comparisons
have exact rational certificates: (3/2)^3>3 and (5/2)^3<16. In particular no
floating-point logarithm or incorrectly rounded bound is load-bearing.

This proves coprimality for any two independent full rows, not merely for
primitive rows or the rows picked by the finite search. Possible roots of
unity in the eliminated equations have not been incorrectly replaced by 1.

## 5. Height and integer-content conclusion

For the candidate coefficient envelope, each summand is the exponential of
a convex piecewise-linear function of (a,b,c), hence is convex. A convex
function on the box |a|,|b|<=H, 0<=c<=H is bounded by the maximum over its
vertices. At each vertex its two terms are A^H+B^H with A,B positive and
H>=1, so this is at most (A+B)^H. The eight values at H=1 are exactly

    4, 10, 7, 5, 7, 22, 4, 14.

Thus the bound 22^H is valid uniformly, including H>1 and all signs. A crude
separate maximization of the two terms would lose this useful constant.

Minkowski's second theorem was checked against [Barvinok's lattice notes,
Theorem 22.5, printed page 46](https://sites.lsa.umich.edu/barvinok/wp-content/uploads/sites/1434/2025/05/latticenotes669.pdf).
For the cube of volume 8 and this determinant-p lattice it gives
lambda_1 lambda_2 lambda_3<=p. Since lambda_3>=sqrt(lambda_1 lambda_2),
there are independent vectors with H_1 H_2<=p^(2/3). The proof does not
incorrectly claim that both H_i separately are at most p^(1/3).

The resultant of the two cleared polynomials is a nonzero integer. Its
integral Bezout identity evaluated at z modulo g proves g divides it. This
is full integer divisibility, not radical divisibility. If one polynomial is
a nonzero constant, its resultant is the appropriate positive power of that
constant; g divisibility still follows. Two constants are excluded above.

The Sylvester-Hadamard estimate uses d_i<=2H_i and ||F_i||_2<=||F_i||_1,
giving exactly

    g <= |Res(F_1,F_2)|
      <= 22^(H_1 d_2+H_2 d_1)
      <= 22^(4 H_1 H_2)
      <= 22^(4 p^(2/3)).

No squarefreeness, prime-only reduction, or one-sided clock-sign argument is
hidden in the final inequality.

## 6. Independent exact computation

The independent script is `audit_outward_gap_one.py`. It imports no producer
or repository mathematics. It uses only Python exact integers and rationals,
plus SymPy 1.14.0 for QQ denominator clearing, polynomial gcds, and resultants.
The producer instead builds binomial coefficient arrays and computes explicit
Sylvester determinants by Bareiss elimination.

The independent physical universe is all coprime 2<=p<=120, 1<=r<p, with
source positions built by inverse floor jumps ceil(jq/p)-1. It enumerates
four-distinct-position cyclic `0110` patterns, compares their residue labels
with the exact legal interval, executes both orders, evaluates the physical
carry by a reverse suffix scan, and verifies the whole-modulus gcd equality.
For p<=20 it additionally closes the actual branch map on the rational point
B/D using Fraction arithmetic.

Results:

- 4,385 primitive source shapes.
- 58,535 unordered outward pairs, hence 117,070 ordered actual paths.
- 4,384 ordered paths cross the original cut.
- 30,277 unordered positive-clock and 28,258 unordered negative-clock paths.
- 7,367 short-pair resultants, one for every legal parameter through p=60.
- 1,795 arbitrary independent row pairs of sup-height at most 2, covering
  every legal parameter through p=12.
- All polynomial gcds and resultants pass; every whole g divides its resultant.
- No integral output in this finite universe.
- Maximum observed g is 4217 at (p,r,u)=(116,59,24).
- 675,812 explicit checks, all active under `python -O`.

The independent normal and optimized receipts are byte-identical. The
producer's normal and optimized reruns also match its frozen receipt exactly.
Their source hash is checked within the emitted receipts. These are finite
supporting computations, not substitutes for the all-parameter proof.

### Prime-power and resultant-surplus controls

The receipt retains actual square factors, including:

    (p,r,u)=(26,17,5): g=49=7^2,
    (p,r,u)=(29,14,4): g=25=5^2.

The full g, not merely 7 or 5, is tested against the common modular root and
the exact resultant. Dependent-projection and constant controls include the
first of these square-depth cases.

The resultant can still have irrelevant clock factors. For example, at
(p,r,u)=(3,2,0), the actual content is 1, while a dependent-projection pair
has resultant -5 and a constant-row pair has resultant 5. Since D=5,
`gcd(|D|,|Res|)` is 5, not 1. The proof correctly uses only g divides Res.

### Illegal smooth common-root hostile

At p=3,r=1,u=9, the independent checker reconstructs

    F_1=2X-3,
    F_2=4X^2-8X+3=(2X-1)(2X-3),

checks their common root 3/2, and finds resultant zero. The full row vectors
(1,0,1) and (0,2,0) are independent, so this really refutes an unrestricted
smooth-root lemma. The failure is at the first valuation-pattern step:
2p-u=-3 is negative and no longer excludes a negative 2-adic valuation.
This witness must remain outside the legal universe and inside the audit
record. It is not a counterexample to the scoped theorem.

## 7. Reproduction and exact stopping point

Run from the packet directory:

    python3 audit/audit_outward_gap_one.py > audit/receipt.replay.json
    python3 -O audit/audit_outward_gap_one.py > audit/receipt.optimized.replay.json
    cmp audit/receipt.json audit/receipt.replay.json
    cmp audit/receipt.json audit/receipt.optimized.replay.json
    python3 verify_outward_gap_one.py > audit/producer.replay.json
    python3 -O verify_outward_gap_one.py > audit/producer.optimized.replay.json
    cmp receipt.json audit/producer.replay.json
    cmp receipt.json audit/producer.optimized.replay.json

`SHA256SUMS` binds the audited candidate, producer, original receipts, and this
audit's source and receipts. The only wording clarification worth making in
a future edit is to spell out “2<=p<=60” in the producer census, since its
1,101 shapes intentionally include the p=2 boundary although the theorem
starts at p=3. The subsequent boundary sentence already makes the intended
universe clear; this is not a proof defect.

Direct progress: the new two-prime valuation-and-norm noncancellation argument
passes independent audit, and its actual full-content consequence is sound.
The useful new connection is from physical legality to simultaneous normalized
valuation constraints, then to a contradiction between finite and infinite
norm bounds. The destroyed coordinate in a bare resultant is still the chosen
clock root, witnessed by the surplus-factor controls.

No all-length integral exclusion, signed finite-descent transfer, p>=2744
reduction, or 904-path closing census is asserted or checked by this audit.
Those remain separate work with separate dependencies and receipts. No change
to the global Collatz status follows.
