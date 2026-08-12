---
id: THM-3335
title: "Square-triangular Pell/Markov selector, square-leg Pythagorean compiler, skew-EW gate, and seven-clock LRC benchmark"
status: >
  PROVED + VERIFIED-EXACT + CITED CANNONBALL CLASSIFICATION /
  INDEPENDENT HOSTILE AUDIT PASSED.  The
  square-triangular solutions are exactly the consecutive-parameter primitive
  Pythagorean triples whose even leg is an even square.  A half-Hadamard map
  identifies the same orbit with the fixed-two Markov branch and gives a
  two-state exact sequence compiler.  The orbit selects precisely the
  arithmetic skew-Ehlich--Wojtas candidate orders whose offset by two is an
  even square, but supplies no design-existence theorem.  Separately, square universal
  tournament arc count and skew-EW attainment are compatible only at order
  two.  In seven labelled LRC one-tail planes the same consecutive spinor has
  exact maxima 1/j for j=7,...,13; j=6 fails sharply.  Six planes have the
  same Gaussian norm, Kelvin defect, and signed hull owner but different
  maxima, proving that the labelled clock placement is an essential sidecar.
  No LRC(14), planar-JC, or FC(3) frontier is closed.
source: codex-2026-08-12-pythagorean-pell-synthesis
depends_on:
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-476-skew-ew-square-law
  - THM-1920-the-spectral-insertion-response
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-928-two-scale-certificate
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
related:
  - THM-1880-the-a-b-functional-frame-chebyshev-pell-companions
  - THM-2142-the-half-angle-bridge-ab-monoid-is-the-ctu-cyclotomic-skeleton
  - THM-900-guy-square-triangular-bothclean
  - THM-2057-scaled-zeta-core-one-tail-closure
  - THM-3010-ballot-column-newton-ratios-and-metallic-alternation
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
  - THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors
  - HYP-3075-lrc14-hurwitz-markov-pell-cannonball-carrier
script: 04-computation/square_triangular_pell_markov_pythagorean_thm3335.py
output: 05-knowledge/results/square_triangular_pell_markov_pythagorean_thm3335.out
script_sha256: 9f543e1f319405ace4d33de5d950a1947b919469eed3e71ca06fcb33fb69f9b7
output_sha256: 02011cccf205a9f88965f1988591a42d63d0c327a824d1c6da2611ddf42b6ba2
hash_basis: LF-normalized bytes
---

# THM-3335 -- the square-triangular selector

**PROVED + VERIFIED-EXACT + CITED CANNONBALL CLASSIFICATION /
INDEPENDENT HOSTILE AUDIT PASSED.**

This is a repository synthesis and proof interface.  No literature-priority or
global-novelty claim is made.  The claims about square-triangular numbers and
Pell equations are elementary; the new payload here is their typed assembly
with the current Pythagorean, Markov, tournament, LRC, and sequence interfaces,
together with the exact loss boundaries.

## 1. Inheritance pass and notation

The closest proved mechanism is [THM-3333, the Gaussian-square Pythagorean
light cone](THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md).
For `m>n>0` it uses

```text
Phi(m,n)=(m^2-n^2, 2mn, m^2+n^2).                       (1)
```

Its exact side difference is

```text
C-B=(m-n)^2.                                             (2)
```

Thus the unit side difference `C-B=1` is not a vague triangular analogy: it
forces consecutive Euclid parameters `m=n+1`.

Three canonical hostiles show why both conditions below are load-bearing:

```text
(5,12,13)=Phi(3,2):     C-B=1, but B is not a square;
(63,16,65)=Phi(8,1):    B is a square, but C-B=49;
(15,36,39)=3(5,12,13): scaling can manufacture a square even leg, but
                        destroys the unit gap.            (3)
```

The corrected near miss is [THM-1880, the a/b Chebyshev--Pell
frame](THM-1880-the-a-b-functional-frame-chebyshev-pell-companions.md): its
first coupled recurrence transposed the even and odd outputs.  THM-2142 had
recorded the correct crossed formula, and MISTAKE-367 now repairs the original
truth surface.  The defining closed forms used here were never affected.
THM-900's Guy identity is also inherited only in its proved `G/Q` scope; its
historical both-clean search verdict was superseded by MISTAKE-151 and
THM-928(C).

The least-used relevant sidecar is
[HYP-3075, the Hurwitz--Markov--Pell/cannonball carrier](../../05-knowledge/hypotheses/HYP-3075-lrc14-hurwitz-markov-pell-cannonball-carrier.md).
It observed the scalar coincidence `29*169=4901=70^2+1` but correctly left it
as address-only.  Sections 2--4 identify the full infinite compiler that
contains that row.

Write

```text
T_n=n(n+1)/2.                                             (4)
```

For selector depth `k>=0`, the symbols `n_k,q_k,x_k,s_k` below are reserved
for the square-triangular index, square root, odd Pell coordinate, and even
Pell coordinate respectively.  Tournament orders are denoted `R` or `W`, not
`n`, to avoid conflating the two order maps in Section 6.

The matrix compiler includes the useful base state `k=0`, where
`(n,q,x,s)=(0,0,1,0)`.  That state is a degenerate null point, not a positive
Pythagorean triple; the positive equivalence theorem therefore begins at
`n,q>0`.

## 2. Five exact faces of one selector

### Theorem 2.1 (lossless equivalence)

For integers `n,q>0`, put

```text
x=2n+1,                 s=2q.                            (5)
```

The following are equivalent:

1. `T_n=q^2`.
2. `x^2-8q^2=1`, equivalently `x^2-2s^2=1`.
3. The consecutive-spinor image is

   ```text
   Phi(n+1,n)=(x,s^2,s^2+1),                             (6)
   ```

   a primitive Pythagorean triple whose even leg is an even square and whose
   hypotenuse exceeds that leg by one.
4. With

   ```text
   a=x-s,                    b=x+s,                      (7)
   ```

   the pair lies on the positive fixed-two Markov branch seeded by `(1,1)`:

   ```text
   a^2+b^2+4=6ab,            ab=s^2+1.                  (8)
   ```

5. The integer

   ```text
   W=s^2+2=4q^2+2                                             (9)
   ```

   satisfies the THM-476 arithmetic skew-Ehlich--Wojtas gate and a square
   offset:

   ```text
   2W-3=x^2,                 W-2=s^2 with s even.        (10)
   ```

The inverse maps are

```text
n=(x-1)/2,       q=s/2,       x=(a+b)/2,       s=(b-a)/2,
W=ab+1.                                                     (11)
```

The phrase **arithmetic skew-EW candidate** in (5) means only that the proved
necessary square gate from THM-476 is met.  It does not assert that a skew-EW
tournament or matrix exists.

*Proof.*  From (4)--(5),

```text
x^2-8q^2=(2n+1)^2-4n(n+1)=1,                             (12)
```

which proves (1) iff (2).  Equation (1) gives

```text
Phi(n+1,n)=(2n+1,2n(n+1),2n^2+2n+1)
           =(x,4q^2,4q^2+1)=(x,s^2,s^2+1).              (13)
```

Consecutive integers are coprime and of opposite parity, so this Euclid triple
is primitive.  Conversely, `x^2+s^4=(s^2+1)^2` is exactly
`x^2-2s^2=1`, recovering (2).

For (4), the half-sum/difference substitution (7) gives

```text
a^2+b^2-6ab+4=-4(x^2-2s^2-1),
ab=x^2-s^2=s^2+1.                                       (14)
```

It therefore carries the Pell conic to the fixed-two Markov conic and is
invertible on this parity branch.  The asserted orbit is exhaustive, not only
generated: for any positive solution with `q>1`, the inverse-unit step

```text
(x,q) -> (x',q')=(3x-8q,3q-x)                           (14a)
```

again solves `x'^2-8q'^2=1`.  Since

```text
8q/3 < x < 3q                                           (14b)
```

(the upper bound uses `q>1`), both new coordinates are positive and
`q'<q`.  When `q=1`, necessarily `x=3`, and (14a) gives `(1,0)` directly.
Iteration therefore reaches the seed `(1,0)`; reversing the unique steps
gives (19), and the half-Hadamard map gives the entire positive fixed-two
Markov branch exactly once.  Finally (9)--(10) are immediate from
`x^2=2s^2+1`; conversely (10) with even `s` returns (2).  QED

### Corollary 2.2 (intrinsic triple characterization)

Every positive integral Pythagorean triple with an **ordered even leg** `B`
and `C-B=1` has the form

```text
(A,B,C)=(2n+1,2n(n+1),2n(n+1)+1).                       (15)
```

It is automatically primitive.  Its even leg is a square iff `T_n` is a
square, in which case that leg is the square of an even integer.

Indeed, `A^2+B^2=(B+1)^2` gives `A^2=2B+1`; hence `A` is odd,
`A=2n+1`, and `B=(A^2-1)/2=2n(n+1)`.  Consecutive side and hypotenuse force
their common divisor, and hence the triple's common divisor, to be one.

## 3. The ambient square-even-leg locus and its sharp intersection

The selector is much thinner than the locus of primitive Pythagorean triples
with a square even leg.

### Theorem 3.1 (two valuation sheets)

Let `m>n>0`, `gcd(m,n)=1`, and let `m,n` have opposite parity.  Then

```text
2mn is a square
iff
(m,n)=(u^2,2v^2) or (m,n)=(2u^2,v^2)                    (16)
```

for coprime positive integers `u,v`.

*Proof.*  Every odd prime divides at most one of `m,n`; its exponent in
`2mn` must therefore be even.  The odd parts of both parameters are squares.
Exactly one parameter is even.  If its 2-adic valuation is `e`, the exponent
of two in `2mn` is `e+1`, so `e` is odd.  The even parameter is consequently
twice a square, while the odd parameter is a square.  The converse and the
coprimality statement are immediate.  QED

Intersecting (16) with `m-n=1` gives the alternating Pell sheets

```text
u^2-2v^2=1                    or
2u^2-v^2=1.                                             (17)
```

This explains the controls in (3): square leg and unit side difference are
independent selectors, while unit side difference itself forces primitivity.
The scaled control shows why a square leg seen after arbitrary scaling does
not recover the selector.

### Corollary 3.2 (a sparse selector on the Berggren parabolic spine)

Let `U` be THM-3334's parabolic Berggren matrix

```text
U=[1 -2 2; 2 -1 2; 2 -2 3].                            (17a)
```

For every `r>=1`, its parameter intertwining law gives

```text
U Phi(r+1,r)=Phi(r+2,r+1),
Phi(r+1,r)=U^(r-1)(3,4,5).                              (17b)
```

Consequently the square-triangular triples are exactly the sparse nodes

```text
Phi(n_k+1,n_k)=U^(n_k-1)(3,4,5),
n_k-1=0,7,48,287,1680,...                               (17c)
```

selected from the full unit-gap spine by the square-even-leg predicate.
This is an operation-level bridge: the Berggren step preserves primitivity
and `C-B=1`, while the Pell state (19) decides which spine positions have a
square even leg.

It is not a new Farey path on the selected nodes.  With fixed cusp
`s=(1,1)` and `d_r=(r+1,r)`,

```text
det(s,d_r)=-1,                 det(d_r,d_t)=t-r,
<Phi(d_r),Phi(d_t)>_L=2(t-r)^2.                          (17d)
```

Thus all selected vertices lie in one fixed-cusp Farey fan, but successive
positive selector depths have gaps `7,41,239,...` and are not Farey adjacent.
The retained sidecar is Pell depth or Berggren position; an unweighted Farey
distance through the common cusp collapses their unbounded determinant and
Lorentz separation.

## 4. Markov mutation and an exact closed-form sequence compiler

Put

```text
lambda=3+2 sqrt(2),              lambda^(-1)=3-2 sqrt(2),
x_k=(lambda^k+lambda^(-k))/2,
q_k=(lambda^k-lambda^(-k))/(4 sqrt(2)),
n_k=(x_k-1)/2.                                            (18)
```

These are integers.  Starting from `(x_0,q_0)=(1,0)`, they obey

```text
[x_(k+1)]   [3 8][x_k]
[q_(k+1)] = [1 3][q_k].                                  (19)
```

Thus the `k`th row is computable using `O(log k)` exact matrix
multiplications by binary powering.  This is an operation-count statement;
the integers themselves have `Theta(k)` bits, so no constant-bit-complexity
claim is intended.

The first positive rows are

```text
k  n_k   q_k   x_k   s_k   B_k=s_k^2   C_k=B_k+1   W_k=B_k+2
1    1     1     3     2          4            5            6
2    8     6    17    12        144          145          146
3   49    35    99    70       4900         4901         4902
4  288   204   577   408     166464       166465       166466. (20)
```

The scalar recurrences are

```text
x_(k+2)=6x_(k+1)-x_k,
q_(k+2)=6q_(k+1)-q_k,
n_(k+2)=6n_(k+1)-n_k+2,

Y_(k+2)=34Y_(k+1)-Y_k+2,             Y_k=q_k^2,
B_(k+2)=34B_(k+1)-B_k+8,
C_(k+2)=34C_(k+1)-C_k-24,
W_(k+2)=34W_(k+1)-W_k-56.                              (21)
```

Let

```text
D_6(z)=1-6z+z^2,                 D_34(z)=1-34z+z^2.      (22)
```

The ordinary generating functions, all starting at `k=0`, are

```text
X(z)=sum x_k z^k=(1-3z)/D_6(z),
Q(z)=sum q_k z^k=z/D_6(z),
N(z)=sum n_k z^k=z(1+z)/((1-z)D_6(z)),

Y(z)=sum q_k^2 z^k=z(1+z)/((1-z)D_34(z)),
B(z)=4Y(z),
C(z)=sum C_k z^k=(1-30z+5z^2)/((1-z)D_34(z)),
W(z)=sum W_k z^k=2(1-32z+3z^2)/((1-z)D_34(z)).           (23)
```

In particular,

```text
W_k=(lambda^(2k)+lambda^(-2k)+14)/8.                    (24)
```

These are structural rational functions derived from (19), not fitted
interpolants.  The verifier independently expands each rational function and
compares it with the matrix orbit.

### The fixed-two Markov clock

Let the ordinary Pell numbers be

```text
P_0=0,       P_1=1,       P_(r+2)=2P_(r+1)+P_r.          (25)
```

Then for `k>=1`,

```text
a_k=x_k-2q_k=P_(2k-1),
b_k=x_k+2q_k=P_(2k+1),
a_k b_k=P_(2k)^2+1=4q_k^2+1=C_k.                        (26)
```

The Markov Vieta mutation and the Pell unit are conjugate under the
half-Hadamard map:

```text
(a,b) -> (b,6b-a)

(x,s)=((a+b)/2,(b-a)/2)
      -> (3x+4s,2x+3s).                                 (27)
```

Thus the Markov seed/mutation word, not the scalar `ab`, is the faithful
arithmetic address.

### The corrected THM-1880 evaluator

For the defining THM-1880 polynomials

```text
E_r(X)=((X+1)^r+(X-1)^r)/2,
O_r(X)=((X+1)^r-(X-1)^r)/2,                              (28)
```

one exact specialization compiles all three sequence channels:

```text
E_(2k)(sqrt(2))=x_k,
O_(2k)(sqrt(2))=2 sqrt(2) q_k,
E_(2k+1)(sqrt(2))=sqrt(2) P_(2k+1).                     (29)
```

The correct coupled recurrence, now repaired by MISTAKE-367, is

```text
E_r=X E_(r-1)+O_(r-1),
O_r=E_(r-1)+X O_(r-1).                                  (30)
```

This makes (29) a useful parity-typed evaluator.  It does not identify an
arbitrary tournament with a Pell orbit: `E_r` is the transitive skew
characteristic polynomial, whereas every tournament merely shares the
universal arc coefficient from THM-1920.

## 5. The cannonball row and even/odd square channels

At `k=3`, (20) and (26) give

```text
(a,b)=(29,169),       (x,s^2,C)=(99,4900,4901),
29*169=4901,          4900=70^2.                         (31)
```

This turns HYP-3075's visible `29,70,169` coincidence into one addressed
Markov/Pell/Pythagorean row.  It also meets the exceptional cannonball identity

```text
sum_(r=1)^24 r^2=70^2.                                  (32)
```

The requested even/odd square split is literal:

```text
sum_(r=1)^m (2r)^2     =2m(m+1)(2m+1)/3,
sum_(r=1)^m (2r-1)^2   =m(4m^2-1)/3.                    (33)
```

At `m=12`,

```text
70^2=2600+2300,
99^2+(2600+2300)^2=4901^2.                              (34)
```

This is an exceptional intersection, not a general orbit law.  The classical
cannonball classification, imported through [Bennett's primary-source record](../../05-knowledge/reference/CORE-PAPERS-PYTHAGOREAN.md),
proves that the positive solutions of

```text
sum_(r=1)^N r^2=s^2                                      (34a)
```

are exactly `(N,s)=(1,1),(24,70)`.  Since every positive selector coordinate
`s_k=2q_k` is even, (31)--(32) give the **global** selector/cannonball
intersection

```text
(k,N,s_k)=(3,24,70).                                     (34b)
```

The exact companion independently reproduces the cannonball solutions through
`N=100000`.  In contrast, the question whether the square-triangular root
`q_k` is itself triangular is only **FINITE-EXACT** here: the scan through
`k=30` finds `q=1=T_1` and `q=6=T_3`, giving square-triangular values `1` and
`36`.  Its all-range Diophantine classification remains open in this theorem.

## 6. Two tournament-order maps: one search gate and one no-go

THM-1920 proves that every tournament on `R` vertices has exactly

```text
T_(R-1)=binom(R,2)                                       (35)
```

arcs.  Therefore the square-triangular indices give the **square universal
arc-count orders**

```text
R_k=n_k+1=2,9,50,289,1682,... .                         (36)
```

This scalar property is universal at fixed order and forgets every orientation
and isomorphism coordinate.

The distinct skew-EW order map from (9) is

```text
W_k=4T_(n_k)+2=2,6,146,4902,166466,... .                 (37)
```

For every odd `x`, a general THM-476 arithmetic candidate

```text
W=(x^2+3)/2                                              (38)
```

has the consecutive-spinor Pythagorean header

```text
(x,W-2,W-1)=Phi((x+1)/2,(x-1)/2).                       (39)
```

The square-triangular selector is exactly the subfamily in which `W-2` is an
even square.  Thus `W=146` is a concrete sparse search address, not a promised
construction.  The inherited hostile `W=86 -> (13,84,85)` passes the square
gate but the THM-476 multiplier-symmetric two-circulant search is empty; even
that negative is scoped to its ansatz.

### Theorem 6.1 (square arcs exclude nontrivial skew-EW attainment)

If a tournament order `R` has square universal arc count and an EW-attaining
skew tournament at that order, then `R=2`.

*Proof.*  Write `R=n+1` and `T_n=q^2`.  THM-476's necessary condition gives

```text
2R-3=2n-1=r^2                                           (40)
```

for an odd integer `r`.  Hence

```text
8q^2=(r^2+1)(r^2+3).                                    (41)
```

Set

```text
A=(r^2+1)/2,                B=(r^2+3)/4.                 (42)
```

Because an odd square is `1 mod 8`, both are integers.  Moreover
`2B-A=1`, so `gcd(A,B)=1`, while (41) says `AB=q^2`.  Therefore `B=v^2` for
some integer `v`.  Now

```text
(2v-r)(2v+r)=3,                                         (43)
```

and positivity forces the two factors to be `1,3`.  Thus `v=r=1`, `n=1`, and
`R=2`.  QED

This is a genuine structural exclusion, not a bounded census.  It also shows
why the two visually similar tournament maps (36) and (37) must not be merged.

## 7. A seven-rung consecutive-spinor clock ladder

The strongest LRC-facing use is a benchmark that exposes exactly what the
Gaussian/Kelvin state forgets.

For `7<=j<=13` and `n>=1`, put `a=n+1` and define the thirteen-speed row

```text
S_(j,n)=a({1,...,13}\{j}) union {(j+1)a-1}.               (44)
```

Equivalently, in the labelled rank-two plane with columns

```text
c_i=(i,delta_(ij)),                                      (45)
```

the parameter direction is the same consecutive spinor

```text
d_n=(n+1,n).                                             (46)
```

Choose

```text
kappa_j=2                 for j odd,
kappa_j=3                 for j=8,10,
kappa_12=5.                                               (47)
```

### Theorem 7.1 (exact seven-clock ladder)

For every `7<=j<=13` and `n>=1`,

```text
M(S_(j,n))=1/j,                                          (48)
```

witnessed at

```text
t=kappa_j/(ja).                                          (49)
```

Every row is primitive.

*Proof.*  In (47),

```text
gcd(kappa_j,j)=1,             2<=kappa_j<=j/2.           (50)
```

Because `2j>13`, the only multiple of `j` among `1,...,13` is the deleted
index `j`.  Hence every core speed `ra` has at (49) distance at least `1/j`;
some `r in {1,...,j-1}` is the inverse of `kappa_j` up to sign, so the core
minimum is exactly `1/j`.

The tail is

```text
(j+1)a-1=ja+n.                                          (51)
```

Its centered numerator at (49) is `kappa_j n`, and

```text
kappa_j n >= n+1=a,
2 kappa_j n < j(n+1)=ja.                                (52)
```

Thus its distance is at least `1/j` and below `1/2`.  This proves the lower
bound in (48).

For the upper bound, the row contains

```text
{a,2a,...,(j-1)a}.                                      (53)
```

Sort the `j` points `0,{at},...,{(j-1)at}` cyclically on the circle.  Their
`j` cyclic gaps sum to one, so one gap has length at most `1/j`.  Let its
endpoints be `{pat}` and `{qat}` with `p!=q`, and put
`r=|p-q| in {1,...,j-1}`.  Their circular separation represents
`(p-q)at mod 1`, so `||rat||=||(p-q)at||<=1/j`.  Hence `M<=1/j`, proving
equality.  Finally
`gcd(a,(j+1)a-1)=1`, so the row is primitive.  QED

The threshold `j=7` is sharp for this theorem.  At `j=6,n=1`, the undeleted
core still contains `12a`, which vanishes at every 6-clock AP extremizer, and
the exact value is

```text
M({2,4,6,8,10,14,16,18,20,22,24,26,13})=2/23 != 1/6.   (54)
```

The independent exact lower-envelope and pair-sum engines agree.  At `n=1`,
the lower rungs `j=2,3,4,5,6` have values

```text
2/17, 2/17, 2/19, 2/19, 2/23.                           (55)
```

Equation (48) is an explicit family within the general safe-killer principle
already visible in THM-757; the new content is the unified seven-rung
consecutive-spinor, owner, and Kelvin placement.

### Kelvin gate versus actual phase

For (45)--(46), the Gaussian norm is

```text
C(d_n)=2n^2+2n+1.                                       (56)
```

The determinant maximum and unique signed axial owner are

```text
D_j(d_n)=13n,       owner=-c_13,       7<=j<=12;
D_13(d_n)=12n,      owner=-c_12.                         (57)
```

Indeed, for `i!=j`, `|det(d_n,c_i)|=in`, while the special column has value
`|1-(j-1)n|`; (57) follows immediately.  The THM-2056/3333 Kelvin defect
`F=C-91D` is therefore

```text
F_j(n)=2n^2-1181n+1,             7<=j<=12,
F_13(n)=2n^2-1090n+1.                                  (58)
```

The first nonnegative integer values occur at

```text
n=591                    for j<=12,
n=545                    for j=13.                       (59)
```

Thus the square-triangular selector meets the Kelvin residual in exactly

```text
n=1,8,49,288                                                (60)
```

on every rung.  All four are nevertheless closed by the explicit clock (49).
This is a clean false-positive calibration for a sufficient geometric gate.

More sharply, fix any `n`.  The six rows `j=7,...,12` have identical

```text
Phi(d_n), C(d_n), D(d_n), signed owner, and F(d_n),       (61)
```

as a **compressed Gaussian/Kelvin scalar tuple**, but their full labelled
column configurations and polar polygons differ.  Their true maxima are

```text
1/7,1/8,1/9,1/10,1/11,1/12.                             (62)
```

There is no contradiction with THM-3333's fixed-owner scalar recurrence: that
recurrence is inside one labelled owner plane.  Equations (61)--(62) prove
that it cannot be transported **across** labelled planes as a complete phase
state.  The destroyed coordinate is the non-hull tail/clock placement `j`.

Finally, each selector-coordinate sequence

```text
n_k,q_k,x_k,s_k,B_k,C_k,W_k,a_k,b_k                    (62a)
```

(indexed by `k`, with positive terms only) grows term-to-term by a ratio at
least `7/3`.  For `x_k,q_k` this follows directly from (19); `n_k=(x_k-1)/2`,
the squared/offset coordinates grow faster, and `a_k,b_k` obey the
six-recurrence in (27) with first ratio five.  This statement does **not**
apply to the undecimated Pell sequence `P_r` in (25).  Any thirteen consecutive
positive terms from one sequence in (62a) are therefore LRC(14)-safe by
THM-928(A').  Such a selector-coordinate packet cannot inhabit the hard
frontier.  Its honest use is as a mixed coordinate address, clock benchmark,
or hostile control.

## 8. Two cross-frontier loss ledgers

### Factorial/Gaussian torus quotient

[THM-3300](THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go.md)
identifies the factorial class with the `U(1)^n`-invariant Gaussian
subalgebra.  In one complex coordinate that quotient retains

```text
C=|z|^2                                                   (63)
```

but kills the weight-two coordinate

```text
z^2=A+iB                                                  (64)
```

needed by THM-3333's represented-sum/Farey carrier.  The exact hostile is

```text
u=(8,1),          u'=(7,4),          v=(1,0),
C(u)=C(u')=65,
|det(u,v)|=1,     |det(u',v)|=4,
<Phi(u),Phi(v)>_L=2,
<Phi(u'),Phi(v)>_L=32.                                   (65)
```

So the direct norm-only port to the factorial algebra is impossible.  A
weight-two `(A,B)` sidecar outside that invariant quotient would be required;
this observation proves nothing about FC(3).

### Planar-JC metallic alternation

THM-3010's full silver recurrence has characteristic roots `1+-sqrt(2)` of
product `-1`, and hence the alternating Cassini law

```text
P_(r-1)P_(r+1)-P_r^2=(-1)^r.                             (66)
```

The square-triangular selector keeps only even indices:

```text
q_k=P_(2k)/2,
q_(k+1)=6q_k-q_(k-1),
q_(k-1)q_(k+1)-q_k^2=-1.                                (67)
```

Its characteristic roots are `3+-2sqrt(2)`, of product `+1`; the alternating
norm sign has become constant.  Consequently every defined Newton ratio is
`q_k^2/(q_k^2-1)>1` rather than maximally alternating.  Even-power decimation
therefore destroys precisely the parity sidecar used by the metallic circuit
stratification.  This is a stopping boundary, not a planar-JC bridge.

## 9. Connection contract and honest frontier

| Source | Target / map | Preserved predicate | Destroyed information | Sidecar / decisive control |
|---|---|---|---|---|
| `T_n=q^2` | `Phi(n+1,n)` | Pell equation, primitive triple, square even leg, unit side difference | none; inverse is (11) | hostile pair in (3) tests both selector keys |
| Pell point `(x,s)` | fixed-two Markov pair `(x-s,x+s)` | conic, mutation, product/hypotenuse | branch orientation if pair is unordered | retain seed and mutation depth |
| selector row | skew-EW order `W=s^2+2` | THM-476 necessary square gate | every matrix/tournament sign and Gram coordinate | `W=86` necessary-only hostile |
| tournament order `R` | arc count `T_(R-1)` | universal coefficient of `char_S` | all orientations and isomorphism data | THM-1880 polynomial only on transitive tower |
| LRC plane `(j,n)` | `(C,F,owner)` | determinant gate inside the labelled plane | tail placement and clock/phase | six-plane hostile (61)--(62) |
| Gaussian spinor | factorial norm quotient | radial norm | represented sum, angle, determinant shell | weight-two `(A,B)`; norm-65 hostile |
| full silver orbit | even Pell subsequence | magnitude and norm-`+1` recurrence | original index parity / alternating Cassini sign | retain parity before circuit transfer |

Direct LRC progress is an exact benchmark and no-go: the selector exposes four
Kelvin false positives per rung and proves which clock coordinate repairs them,
but every row is already safe and no open ledger row closes.  Tournament
progress is a sparse arithmetic search family plus the global square-arcs/EW
incompatibility theorem, not a new EW design.  The main reusable artifact is
the two-state exact compiler (19), with the fixed-two Markov seed/mutation word
as its address sidecar.

## 10. Verification and scope

Run

```bash
python3 04-computation/square_triangular_pell_markov_pythagorean_thm3335.py
python3 -O 04-computation/square_triangular_pell_markov_pythagorean_thm3335.py
```

Both modes reproduce the stored transcript byte for byte.  The verifier uses:

- an independent direct square-triangular scan through `n=100000`, compared
  with the Pell recurrence;
- every primitive Euclid parameter pair with `m<=500`, finding `134`
  square-even-leg triples and checking the two valuation sheets;
- selector depths through `k=30`, fourteen Markov mutations, six rational
  generating functions, and exact arithmetic in `Q(sqrt(2))`;
- the square-pyramidal scan through `100000`, corroborating the cited global
  classification, plus the exact even/odd split;
- the seven LRC rungs through `n=200`, with independent pair-sum and full
  lower-envelope exact-max engines at the sharp `n=1` boundary;
- positive controls `j=7,...,13`, lower-rung hostiles `j=2,...,6`, the
  norm-`65` represented-sum collision, and the `W=86` necessary-only control.

Equations (5)--(17), (19), (21), (26)--(30), (35)--(43), (44)--(53), and
(56)--(67) are proved for all integers in their displayed domains.  The
bounded scans in Section 5 are explicitly FINITE-EXACT and make no all-range
claim.  The skew-EW statements use THM-476 only as a necessary gate except for
Theorem 6.1's exclusion.  The Kelvin gate remains sufficient, not necessary.
No tournament orientation is manufactured from the symmetric Farey or Pell
relations, and no result here proves LRC(14), FC(3), HFC(3), or planar JC.
