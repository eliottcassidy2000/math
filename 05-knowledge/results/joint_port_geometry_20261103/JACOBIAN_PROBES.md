# A mixed response recovers the missing polynomial

**PROVED elementary reconstruction and image criterion; CITED kernel theorem.**
Local labels JP1--3. A constant-Jacobian mate is supplied throughout.
These statements neither construct such a mate nor prove a polynomial
inverse for an arbitrary Keller map.

## JP1. A split exact two-query observer

Let k have characteristic zero, V=k[x,y], D=Jac(P,-), and let B be supplied
with DB=kappa in k*, a nonzero constant. For any unknown Q define

`R1=DQ, R2=D(BQ)`.

The derivation rule gives the complete inverse

`Q=(R2-B R1)/kappa`.                                       (J1)

The image consists exactly of pairs satisfying the single polynomial
differential equation

`C(R1,R2):=D R2-B D R1-2 kappa R1=0`.                       (J2)

Indeed applying D to (J1) gives DQ=R1+C/kappa. If C=0, the first response
is right and the product rule gives the second. Conversely every actual
response has C=0. More strongly, the linear map

`(R1,R2) -> (Q=(R2-B R1)/kappa, C(R1,R2))`

is an automorphism of V direct-sum V, with inverse

`R1=DQ-C/kappa, R2=B R1+kappa Q`.

Thus `0 -> V -> V^2 -> V -> 0`, with maps
`Q -> (DQ,D(BQ))` and `(R1,R2) -> C`, is split exact. This is a complete
image description for this supplied two-query measurement operation.

## JP2. Which extra multiplication probe suffices?

Replace B by any fixed f in V. The product rule gives

`D(fQ)-f DQ=(Df)Q`.                                        (J3)

The pair Q -> (DQ,D(fQ)) is injective **if and only if Df is nonzero**.
If it is nonzero, the domain property of V proves injectivity from (J3).
If Df=0, the two responses both vanish on the nonzero polynomial Q=1.
This is a precise selector threshold, not a search for accidental matching
values.

When a=Df is nonunit, response inversion has the explicit image guard
`R2-f R1 in aV`, followed by `D((R2-f R1)/a)=R1`.
For example f=B^2 gives a=2 kappa B. The response pair (0,1) fails
divisibility; a rational quotient is not a polynomial source.

The existing [THM-2230, planar response-fiber theorem](../../../01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md)
identifies ker D=k[P] using the Cheng--McKay--Wang theorem, as extended in
[Moskowicz, Theorem2.3](https://arxiv.org/html/1802.04685v2).
Consequently, under this supplied-mate hypothesis, JP2 can equivalently be
stated as f not in k[P]. We use that cited kernel identification; the
product-rule reconstruction itself needs only the displayed derivation and
the value of DB. The old rank-one tournament response loss and integral
guards are recovered in [JACOBIAN_FIBRES.md](../spectrum_descent_jacobian_20261102/JACOBIAN_FIBRES.md).

## JP3. Why finite exact closures require a boundary condition

On V, multiplication M_B satisfies `[D,M_B]=kappa I`. On any nonzero
finite-dimensional characteristic-zero space the same identity is
impossible: taking traces gives0=kappa dim(V). Thus there is no nonzero
finite-dimensional representation preserving these two operators and
their exact relation.

For the degree cap on k[B], normalize D to d/dB and let M truncate the
top-degree product. On degrees0,...,d the exact identity becomes

`[D,M]=I-(d+1)E_top`.                                      (J4)

Below the top degree the usual product rule survives; on B^d, multiplication
leaves the cap and is set to zero, giving commutator -d. The correction is
rank one, with an essential coefficient. It cannot be discarded as small.

The [September torsion audit](../planar_jc48_sep06_torsion_audit.md) already
uses a response connection with `[nabla,M_P]=I` to exclude nonzero finite
connection-stable torsion. This is a recovered mechanism. Here it supplies
a concrete warning about transporting finite observer constructions from
tournaments: the tournament cap works because its allowed births cannot
bring discarded terminal states back into the measured sector. The
polynomial multiplication/derivative pair does allow return across a
degree boundary, producing (J4).

[jacobian_observer.py](jacobian_observer.py) checks (J1)--(J4) symbolically
for three supplied triangular-composition mates and degree caps1--12.
The proofs above establish their stated general scope.
