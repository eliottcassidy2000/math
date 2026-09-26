# Prime brackets, quadratic dynamics, and the arithmetic information in theta

**Status:** INHERITED PROVED mechanisms explicitly credited below; PROVED
elementary fibre bijection and finite-rational-chart obstruction; CITED
classical theta identities; FINITE-EXACT controls. Collatz, HYP-9127 and
HYP-9130 remain OPEN. No novelty or literature-priority claim is made.

The useful synthesis has two directions. Gaussian squaring really connects
quadratic dynamics to primitive Pythagorean triples, and preserves an entire
fixed-hypotenuse fibre. But finitely many rational formulas cannot encode
Collatz as a degree-two dynamical system while retaining unbounded input
information. Theta gives a separate arithmetic weighting; its fourth power,
rather than its square, sees every positive integer.

## 1. Inheritance and the working board

Closest recovered sources:

- [Prime/square decoder, section2](decoder_prime_square_20260925.md) and
  [odd-square pairing transitions](procgen_brackets_20260924_pairings_transitions.md)
  already classify the prime set {2,3,11}.
- [Arithmetic braids geometry, sections3,5–6](arithmetic_braids_20260917_geometry.md)
  already proves the small quadratic critical graphs, signed Gaussian
  squaring, fresh odd-leg prime support, denominator growth, and the folded
  angle map. These are recovered results, not new discoveries here.
- [THM-3334, Gaussian collision torsors](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
  retains prime-choice coordinates on equal-hypotenuse fibres.
  [THM-3341, square-hypotenuse transplant](../../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md)
  already distinguishes a Gaussian-square transplant from a tree homomorphism.
- [THM-3756, odd-square ordinal descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md)
  supplies the complete two-coordinate tree chart; [THM-3819, Chebyshev–Pell
  return cocycle](../../01-canon/theorems/THM-3819-chebyshev-pell-berggren-return-conductor-staircase.md)
  retains state-dependent return times rather than identifying unrelated clocks.
- [HYP-9127, cube-swap theta](../hypotheses/HYP-9127-cube-swap-cubic-theta.md)
  and [HYP-9130, cube Pade zero estimate](../hypotheses/HYP-9130-cube-theta-2adic-zero-estimate.md)
  are explicit unresolved p-adic questions, not consequences of the classical
  Jacobi theta identities.

The canonical hostile is bounded real motion with unbounded denominator.
The corrected near miss is treating Gaussian angle doubling as a Berggren
branch operation. The least-used coordinates are the sign sheet, primitive
content, denominator, and the degree of a proposed encoding.

| Object | Operation worth testing | Necessary retained information |
|---|---|---|
| Odd-square brackets | Multiply a number inside its own bracket | Both bracket boundaries |
| Modulo30 wheel | Project primes or a quadratic orbit | Prime factors and residue modulo4 |
| Chebyshev coordinate | Square the Gaussian phase | Sign sheet and denominator |
| Primitive triples | Square the hypotenuse | Odd/even leg convention; ancestry separately |
| Theta coefficients | Project to primitive representations | Multiplicity and square-divisor content |
| Collatz encoding | Use finitely many rational charts | Rational-map degree and source identity |

## 2. What {2,3,11} means exactly

Use B_m=((2m-1)^2,(2m+1)^2], with1 treated separately. An integer n in
B_m has a nontrivial integer multiple in B_m exactly when2n is at most the
upper boundary. For m>=3,

    2((2m-1)^2+1)-(2m+1)^2=4m^2-12m+3>0.

Thus only the first two brackets need inspection. Their integer exceptions
are2,3,4,10,11,12; their prime exceptions are exactly2,3,11. The full prime
and multiplier list is

    (p,k)=(2,2),(2,3),(2,4),(3,2),(3,3),(11,2).

This recovers the examples4,6,8 below9;6 below9;22 below25. A strict upper
boundary deletes the equality3*3=9 but does not change the prime set.
The bracket is essential: with consecutive integer-square brackets instead
of odd-square brackets, the analogous closed-upper prime set is only{2}.
The assertion2p<(p+2)^2, by contrast, is true for every positive p and cannot
single out three primes.

Modulo30 is the wheel for the three distinct smallest primes2,3,5. Its
unit classes are1,7,11,13,17,19,23,29. Every prime greater than5 lies there,
but49 is already a composite survivor. Each unit class contains composites,
for example products(r+30a)(1+30b) with positive a,b.

For f(x)=x^2-2, every wheel unit maps to17 or29 modulo30, and both residues
are fixed. This is a genuine finite quotient, but11 maps to119=7*17.
Even two primes in the same wheel class can have different sum-of-two-squares
types:13 and43 agree modulo30, whereas r2(13)=8 and r2(43)=0. Modulo30
forgets modulo4 because30 is not divisible by4.

## 3. The small quadratic graphs and the exact larger map

The graphs, with their signs specified, are

    x^2:    0->0, -1->1->1;
    x^2-1:  1->0->-1->0;
    x^2-2:  1->-1->-1, 0->-2->2->2.

For any of these integer-coefficient monic quadratics, a nonintegral reduced
rational a/b has next denominator b^2: gcd(a^2+c b^2,b)=1. Induction gives
b^(2^j), so such a point is never preperiodic. The remaining integer escape
check gives rational preperiodic sets{-1,0,1} for c=0,-1, and
{-2,-1,0,1,2} for c=-2. This statement concerns the starting points of these
three fixed maps, separately from the inherited classification of rational
parameters with finite critical orbit.

Let J(z)=z+z^(-1). Direct expansion gives

    J(z^2)=J(z)^2-2.

For z on the unit circle this is2cos(theta)->2cos(2theta), the classical
Chebyshev formula [DLMF18.5.1](https://dlmf.nist.gov/18.5.E1).
The quotient identifies z and z^(-1); the inverse phase is not retained.

For a signed primitive Pythagorean triple(A,B,C), with A odd and B even,
z=(A+iB)/C gives the actual integer lift

    (A,B,C)->(A^2-B^2,2AB,C^2).

The scalar x=2A/C follows x^2-2 exactly. Its real values stay in[-2,2],
but its reduced denominator is C^(2^j). This inherited example separates
bounded real position from bounded arithmetic height. Canonical positive
triples instead use D(a,b,c)=(|a^2-b^2|,2ab,c^2), so their x coordinate
follows |x^2-2|. Restoring the sign sheet is necessary for the unfurled map.

## 4. A complete hypotenuse-fibre bijection

Let P_c be the set of positive primitive triples(a,b,c), ordered so that
a is odd and b even. For every integer c>1, Gaussian squaring induces a
bijection

    D:P_c -> P_(c^2).

This includes empty fibres. Forward primitivity follows because an odd
prime dividing a^2-b^2 and2ab would divide both a,b;2 cannot divide the odd
first coordinate. For surjectivity take any(A,B,c^2) in the target. Its
unique Euclid parameters u>v>0 satisfy

    A=u^2-v^2, B=2uv, c^2=u^2+v^2,
    gcd(u,v)=1, u-v odd.

Thus(u,v,c) is itself a primitive right triangle. Choose a to be its odd
leg and b its even leg. Then D(a,b,c)=(A,B,c^2). Uniqueness of the Euclid
parameters proves injectivity as well. Repetition gives bijections between
P_c and P_(c^(2^j)). Repeated prime exponents cause no exception.

This bijection is compatible with the Gaussian prime-choice fibre modulo
units and global conjugation: squaring doubles valuations but retains which
prime over each split rational prime was selected. It is not a map of
Berggren ancestry. In the inherited U,A,D matrix convention,

    D(3,4,5)=(7,24,25)=U^2(3,4,5),
    D(U(3,4,5))=D(5,12,13)=(119,120,169)=A^2(3,4,5).

The latter lies outside the U subtree. No fixed appended tree word can
implement D on all vertices. The positive result preserves the norm equation,
primitivity and fibre identity; tree address and signs are separate data.
Romik's [primary paper](https://arxiv.org/abs/math/0406512) supplies the
different, piecewise three-branch circle dynamics whose rational expansions
terminate in the Pythagorean tree. Angle doubling is not that descent map.

The exact census independently checks792 primitive triples through c=5000,
and proves finite surjectivity by enumerating every square-hypotenuse target
with c^2<=250000 separately. The all-c proof is the Euclid argument above.

## 5. A finite-rational-chart obstruction to a Collatz quadratic model

**Proposition.** Let T be shortcut Collatz. Let s:N+->S be any colouring by
a finite set, and let each R_s be a rational function over C, finite at the
integers where it is used. Set H(n)=R_(s(n))(n). Let F be a rational map of
degree d>=2. If

    H(T(n))=F(H(n)) for every positive integer n,

then every R_s used on infinitely many integers is constant. In particular
H has finite range. The statement permits arbitrary, nonperiodic colouring.

**Proof.** For each infinite colour class choose one input parity and target
colour t occurring together for infinitely many inputs. On those inputs T
is one nonconstant affine map A. Therefore

    R_t o A = F o R_s

as rational identities, since equality on infinitely many points forces the
cleared numerator polynomial to vanish. Degrees give deg(R_t)=d deg(R_s),
with degree0 for constants. The target class is infinite because A is
injective. Each vertex in the finite graph of infinite classes therefore
has an outgoing edge with this degree relation. Following edges reaches a
cycle, where deg(R)=d^ell deg(R) forces degree0. Propagating backwards gives
degree0 at the starting class. Finite colour classes supply only finitely
many additional output values. QED.

For one chart and F(x)=x^2-2, even inputs already force
R(x/2)=R(x)^2-2. The only constant possibilities are2 and-1. This rules out
a rational one-coordinate explanation of Collatz by Chebyshev dynamics;
it does not rule out nonrational maps, infinitely many charts, variable
return times, or models retaining an unbounded history.

The finite/closed-chart condition matters. On a finite unrolled halving
chain, R_j(x)=F^j(2^j x) satisfies R_(j+1)(x/2)=F(R_j(x)). Its degree2^j
grows, so it cannot close into finitely many recurrent rational charts.
The probe includes this positive control and an exhaustive small Mobius
coefficient census. The proof, rather than that census, establishes the
universal obstruction. Root independently audited the graph/degree argument.

## 6. Theta coefficients recover the fibres, with a primitive-content filter

Write vartheta(q)=theta_3(0,q)=sum_(m in Z) q^(m^2). The standard definition
is [DLMF20.2.3](https://dlmf.nist.gov/20.2.E3). For t>0,

    vartheta(e^(-pi t))^2=sum_(n>=0) r2(n)e^(-pi n t).

Here r2 counts ordered signed representations. Jacobi's identity is
r2(n)=4 sum_(d|n) chi_4(d), where chi_4 is0 on even integers and +/-1 on
1/3 modulo4 [DLMF27.13.5–7](https://dlmf.nist.gov/27.13).
Consequently, for Re(s)>1, absolute convergence and the gamma integral give

    integral_0^infinity (vartheta(e^(-pi t))^2-1)t^(s-1)dt
      = 4 Gamma(s) pi^(-s) zeta(s) beta(s).

The subtraction of1 and the convergence half-plane are essential. This is
an analytic identity, not a primality test for a recurrence.

Removing common integer content by Mobius inversion gives, for c>1,

    #P_c = (1/8) sum_(d|c) mu(d) r2((c/d)^2).

The factor8 accounts for signs and order; primitive axes occur only at c=1,
which is excluded. Euler factors of the divisor formula yield

    #P_c=2^(omega(c)-1) if every prime divisor of c is1 modulo4;
    #P_c=0 otherwise.

The condition also excludes even c. Thus the theta formula independently
confirms #P_c=#P_(c^2). Coefficients alone do not retain the Berggren address;
the explicit bijection in section4 supplies that additional structure.

There is an instructive p-adic hostile. For a product c of r distinct primes
that are1 modulo4, the primitive theta coefficient is8#P_c=2^(r+2).
All summands count positive objects, yet its2-adic valuation is unbounded.
Positivity in the real field does not provide the nondivisible leading
coefficient required by HYP-9130.

Three different expressions must not be conflated:
vartheta(q) is a square-exponent Jacobi series; vartheta(q)^3 counts three
squares; the repo's cubic series Psi_3(q)=sum_(k>=1) q^(k^3) is neither.
HYP-9127 concerns Psi_3 at the rational2-adic argument2^10/3^9, and HYP-9130
asks for low-height Pade forms with a controlled leading valuation. The
classical Mellin identity supplies no such zero estimate.

## 7. A full-support theta carrier for the actual first-descent question

The square theta weight misses integers such as3; its cube misses7. Among
the positive integer powers vartheta^d, the first with nonzero coefficient
at every nonnegative integer is d=4. The failures for d=1,2,3 are2,3,7.
Jacobi's four-square formula gives, for n>0,

    r4(n)=8 sigma_1(n) if n is odd,
    r4(n)=24 sum_(d|n,d odd) d if n is even,

so every weight is positive. The formula is cited on the same
[DLMF page](https://dlmf.nist.gov/27.13).

Fix a real0<q<1 and let B_k be the positive integers n>=2 with T^j(n)>=n
for every1<=j<=k. Define

    A_k(q)=sum_(n in B_k) r4(n) q^n.

**Exact equivalence:** A_k(q)->0 if and only if every positive integer n>1
eventually descends below itself, hence if and only if Collatz holds.
Indeed the sets B_k decrease; their weights are summable, being bounded by
vartheta(q)^4. If all starts descend, dominated convergence gives zero.
If a start n survives forever, A_k(q)>=r4(n)q^n>0 for every k. Strong
induction converts universal first descent to convergence to1.

This is the root lane's full-support atomic criterion with explicit divisor
weights. It supplies no decay theorem. It does locate the required new
estimate: control A_k at a fixed real q while retaining actual parity/carry
and source height. A Haar-density bound or theta-square weighted estimate
cannot replace this because either can assign zero mass to an individual
positive source. Nor may the positivity argument be moved to Q_2.

## 8. Reproduction and next bounded question

Run `python -X utf8 04-computation/experiments/crossroads_crossing_20260926_arithmetic.py`
and with `-O`. The [output](crossroads_crossing_20260926_arithmetic.out) is
exact stdout. Checks cover the first100 odd-square brackets; a complete
mod30 wheel;1520 rational denominator gates; small Mobius charts and an
unrolled positive control;792 primitive triples;80 independently enumerated
square-hypotenuse targets; all r2 coefficients through2000; primitive theta
projection through c=150; and r4 convolution through300. No disabled
assertions or floating-point mathematical decisions are used.

The strongest new construction is the whole-fibre squaring bijection. The
strongest rejected Collatz route is finite rational-chart transport into
quadratic dynamics. The constructive unresolved target is a source-sensitive
decay estimate for A_k(q), possibly exploiting its divisor weights. Any
proposed use of a quadratic or tree model should first pass the degree,
ancestry, sign, and full-support tests above.
