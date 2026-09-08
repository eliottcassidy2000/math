# Independent audit of the DG filtration and recovered cubic floor

**Status: INDEPENDENT ANALYTIC + EXACT SOURCE + REPLAY AUDIT PASS.**
I read the complete [primary note](planar_jc48_sep08_dg_quadratic.md)
and standalone source, recovered the load-bearing planar theorems in full,
and replayed the producer normally and under optimization. Both outputs
agree byte for byte with the frozen output. The result concerns the
specified surface and a source-degree condition on one output-pencil
member. It does not assert a new low-degree planar theorem or close JC(2).

## 1. The all-degree global filtration is complete

On the dense chart of the fixed DG surface,
`z=x^2+1/t`. For fixed `n`, multiplying a polynomial of t-degree at most
`n` by `(z-x^2)^n` gives a unique polynomial `Q(x,z)` of z-degree at most
`n`. This is an invertible change between the two displayed vector
spaces before imposing global regularity; no rational numerator is
silently admitted.

At the generic point of the actual boundary `D`, the parameters are
`r=1/x` and `z=b/(1+r^2 b)`. In particular `z` has residue `b`, which
is transcendental in the residue field `C(b)`. Therefore in

    r^(2n) Q(1/r,z)/(r^2 z-1)^n

the denominator is a unit, and a nonzero leading x-coefficient of `Q`
cannot vanish generically along `D`. This proves the necessary bound
`deg_x Q<=2n`, with no cancellation between distinct leading r-orders.

Conversely the entire bidegree box gives the functions

    x^a t^(n-j)(1+x^2 t)^j,
    0<=a<=2n, 0<=j<=n.

Their second-chart expressions are exactly

    (-1)^n r^(2n-a) b^j(1+r^2 b)^(n-j).

These are polynomials everywhere on that chart, including its locus
`1+r^2 b=0`. Thus generic-boundary regularity was used only for a
necessary degree bound; the converse restores full chart regularity.
Independence follows from the distinct numerator monomials `x^a z^j`.
Consequently the complete dimension is `(2n+1)(n+1)` for every `n>=0`.
For `n=0` this says that global functions independent of t are constants.

The multiplication assertion has its stated linear-span meaning.
Every exponent pair in the box for `m+n` splits into pairs in the boxes
for `m,n`: choose the first x-exponent `min(a,2m)` and the first
z-exponent `min(j,m)`. Both remaining exponents have the required bounds,
and the denominator powers add. This proves `span(L_m L_n)=L_(m+n)`,
including zero indices. Since every global function restricts to an
ordinary polynomial on the dense affine chart, these spaces exhaust
`O(W)`; injectivity of restriction follows from irreducibility.

## 2. The full quadratic layer and its boundary value

I independently separated the b-powers of the second-chart expression of
`F=A(x)t^2+B(x)t+C(x)`:

    b^2: r^8 A(1/r),
    b:   2r^6 A(1/r)-r^4 B(1/r),
    1:   r^4 A(1/r)-r^2 B(1/r)+C(1/r).

The first forces `deg A<=8`. The second forces `deg B<=6`,
`B_6=2A_8` and `B_5=2A_7`. The last forces
`C_4=A_8`, `C_3=A_7`, `C_2=B_4-A_6`, `C_1=B_3-A_5`,
and no higher coefficients of C. There are exactly nine free A
coefficients, five free lower B coefficients and one C constant.
The actual value on D is

    A_8 b^2+(2A_6-B_4)b+(A_4-B_2+C_0).

It retains genuine boundary separators; a collapsed-boundary argument
would not prove this result. The independent incoming linear-layer
classification, read from its tracked blob, agrees with the `n=1`
six-dimensional specialization.

## 3. The recovered planar theorem was checked, not extrapolated

The exact dependency is
[THM-2071, quadratic-fiber square/parity rigidity](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md),
together with
[THM-2063, one-fiber-linear Keller pairs](../../01-canon/theorems/THM-2063-one-fiber-linear-planar-keller-pairs.md).
I read their complete analytic proofs. Guardrail 60 in
`01-canon/ACTIVE-GUARDRAILS.md` correctly distinguishes this source-fibre
degree from generic covering degree and from JC(2).

For the quadratic theorem, minimizing the mate's fibre degree modulo
`C[P]` gives `q_n^2=c A^n`. Positive even `n` is removed by the honest
polynomial target shear whose top coefficient is a scalar multiple of
`A^(n/2)`; zero residual degree contradicts the Jacobian. Odd `n` therefore
forces `A=U^2`, after scalar absorption. In the rational centered coordinate
`z=Uy+B/(2U)`, write `P=z^2+D`. The exact Jacobian operator is

    J(P,-)=U(D' partial_z-2z partial_x),

where the x-derivative holds z fixed. Its even kernel is precisely
`C[P]`, since independence of the powers of P forces the rational
coefficient derivatives to vanish. Removing that even part is thus an
actual polynomial target shear, not merely a rational operation.

For an odd part of degree `2r+1`, the triangular recurrence makes its
coefficients polynomials in D. Their leading terms give, at a putative
finite pole of `h=B/(2U)`, the common highest-pole coefficient

    c sum_(k=0)^r (-1)^k binom(r+1/2,k)
      =c (-1)^r binom(2r,r)/4^r !=0.

All lower terms have strictly smaller pole order because `C` is polynomial
and `D=C-h^2`. The original constant fibre coefficient is still polynomial
after the target shear, so this noncancellation proves that h has no finite
pole. Then D and the recurrence coefficients are polynomial; the final
equation `D'a_0(D)=kappa/U` forces U constant. Comparing degrees gives
`(r+1)deg D=1`, hence `r=0` and `deg D=1`. The displayed polynomial inverse
in the theorem completes the proof for every mate degree. No finite
recurrence sample is used for that unbounded conclusion.

The lower-degree theorem is the same honest leading-coefficient descent:
subtract scalar powers of the first component until the complement
depends only on x, and use the nonzero constant product `-A R'`.
If the first component has fibre degree zero, the two derivative factors
are units directly. Both cases have explicit polynomial inverses.

## 4. The all-degree cubic theorem also applies

A further bounded inheritance search recovered
[THM-2118, all-degree cubic boundary-flux coprimality and source-fibre closure](../../01-canon/theorems/THM-2118-all-degree-cubic-faber-boundary-flux-coprimality.md).
The earlier THM-2084 and THM-2110 finite reduced-degree bounds explicitly
route to this complete theorem. I read THM-2118 in full, including its
uniform pole cleanup, and its load-bearing
[THM-2102, power-free weighted-face descent](../../01-canon/theorems/THM-2102-power-free-weight-face-and-first-defect-descent.md).
The targeted correction search found MISTAKE-244, which repairs the old
Newton-face/DvdK identification by precisely this direct weighted
argument; it does not retract either theorem. Their applicability is
to one cubic coordinate with an unrestricted polynomial mate.

The cubic coefficient law reduces the mate modulo honest `C[P]` shears,
excludes degrees divisible by three and makes the leading cubic
coefficient a cube. In the resulting rational depressed variable, the
full Faber combination has two equations `Phi'=0` and `R'=kappa/U`.
I checked the translated identity

    e_n=[w^n](1+3w+s w^2)^(n/3),
    phi_n=3[w^(n+1)](1+3w+s w^2)^(n/3),

its first-order boundary-flux relation and the polynomial ODE. Away
from its two singular points, a common zero would be a repeated root
whose lowest ODE order cannot cancel. At `s=0,9/4` the two binomial
evaluations are nonzero since `3` does not divide n. Thus this is an
all-degree coprimality proof, independently of its finite gcd referees.

The remaining valuation arguments are essential. If the centering pole
has order H and the depressed quadratic coefficient has pole order rho,
then `rho>2H` gives a unique highest top-flux term with a strictly
positive gap to every lower Faber representative. For `rho<2H`, the
top boundary value has order nH and its finite alternating-binomial
coefficient is nonzero; all lower representatives have smaller order.
At equality, polynomiality of the original constant fibre coefficient
and constancy of the full first flux force the forbidden common zero.
The possible zero leading constant-cubic coefficient on that boundary
is retained. Once the center is polynomial, the remaining depressed
pole is excluded by the unique top pure power in the first flux for
odd n, and by the unique top pure power in the original boundary value
for even n. Hence the depressed coefficients are polynomial, and
`R'=kappa/U` forces the scale U to be constant.

For the resulting honest depressed cubic `z^3+p(x)z+q(x)`, the positive
weight `(6,max(3 deg p,2 deg q))` exposes `z^3` and at least one other
term when p or q is nonconstant. Any proper power would have to be a
cube of a linear polynomial in z. Its leading coefficient is constant
and the absent z-squared coefficient forces the linear polynomial's
constant term to vanish, contradicting the extra exposed term. If both
p and q are constant, the Jacobian equation is directly impossible.

I also checked the cited weighted descent rather than assuming its
conclusion. The weighted Euler/UFD lemma makes any commuting leading
mate form a scalar power of a power-free leading component. Subtracting
that full polynomial power strictly lowers the mate's weighted degree.
At the terminal nonzero bracket the degrees sum to the sum of the source
weights, forcing the first component to have no mixed monomial and one
linear axis term. This gives the polynomial inverse. No DvdK or
generic-cover-degree inference enters. The recovered cubic theorem is
therefore accepted as a load-bearing dependency of the strengthened
DG conclusion, with no new source computation required.

## 5. The global inverse contradiction has the required scope

An invertible constant output transformation can put any nonconstant
pencil member of source-t degree at most three first. Restriction to U0
preserves polynomiality and changes the nonzero constant Jacobian only
by that target determinant. The checked canonical theorems therefore
make the restricted pair a polynomial automorphism, without a degree
bound on the other component.

Its inverse expresses the original x as a polynomial in the two global
functions. Such an expression is regular on W and agrees with x on the
dense U0, hence equals x as a rational function. But `x=1/r` has a genuine
simple pole on the actual nonempty divisor D. This is the contradiction.
It does not require assuming the proposed map is finite or identifies W
with a Keller envelope. An affine source change also leaves x polynomial
in the new source coordinates, so the stated all-linear-directions
extension is valid. It does not identify those different filtrations
with the fixed-chart `L_n`.

The affine-plane control `(t,-x)` has Jacobian one; its second component
has the exact boundary pole `-1/r`. This isolates global regularity as
the missing predicate. The argument excludes the declared low-degree
pencil class on W while leaving the higher-degree global pair question
outside its conclusion.

## 6. Source and replay acceptance

The standalone source uses always-active checks. I inspected its complete
Laurent matrix construction, symbolic fifteen-parameter family, boundary
restriction and multiplication controls. The finite matrix universes for
`n=0,1,2,3` have respective dimensions `1,10,27,52`, ranks `0,4,12,24`
and kernels `1,6,15,28`. The degree boxes are justified by the analytic
numerator bound: expanding its basis gives source x-degree at most `4n`.
Thus these are full bounded kernels, not pruned ansatz tests. Their
finite range does not replace the all-n numerator proof or the recovered
quadratic and cubic theorems. The frozen producer's output labels the
quadratic recovery; the later cubic extension is an analytic inherited
consequence, not an extra claim about those 122 finite gates.

I ran both commands independently and compared each output with the frozen
file. All three are the same **409 bytes**, with **122 always-active gates**.

```bash
python3 -B 04-computation/planar_jc48_sep08_dg_quadratic.py
python3 -B -O 04-computation/planar_jc48_sep08_dg_quadratic.py
```

Frozen SHA-256:

    source b170cb02793e5eea5423aa1a838ee296406e993c18b26c97a033d15a01f9d568
    output 98dbafefbba323c2401b5af303de6ca4c5b2fa8ac37cb28b04acb7baca0761a8

Semantic output digest:

    5f756f8279fc34f19893d021f63f5a1a4feb64b1c071140fb5736cb9947abc6d

No source correction was required. The initial quadratic conclusion was
strengthened by the recovered cubic theorem, with frozen source and
output unchanged. The conclusion is a
recovered canonical obstruction transported to a precisely specified
global surface ring; no literature-priority claim is made.
