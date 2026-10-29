# A compatible CRT Frobenius lift, with the observation selector retained

2026-10-29. **PROVED elementary all-degree identity; FINITE-EXACT application
to an inherited native scalar response.** Same-agent alternative-path
verification. This does not supply a vertex automorphism, an attachment
operation, or a Selmer map.

Incoming078be2b84d proves an actual periodic family of first complements
of the interval tournament has a35-state joint-port evaluator and a
minimal scalar recurrence of order18 over Q. Its ternary scalar response
has the cyclic module

    A = F3[X]/((X+1)^2 g(X)^2),
    g=X^5+X^4-X^3-X^2-1.                         (F1)

The two factors give F3[delta]/delta^2 and F243[epsilon]/epsilon^2.
The field generator alpha has order121, field Frobenius has order5,
and the full multiplication clock T=X has order726=2*3*121. See
[PERIODIC_COMPLEMENT E1--E5](../overnight_rank_return_20261029/PERIODIC_COMPLEMENT.md).
Its main recurrence checker was replayed normally and with -O here;
the general bounded-defect transfer checker was also replayed.

The new point is an explicit compatible lift of Frobenius to this full
response module. It preserves the nilpotent coordinate, but the scalar
observation must be transported as well.

## F1. A general two-jet identity

Let K=Fp(alpha) have degree d, and let alpha have multiplicative order m.
In K[epsilon]/epsilon^2 put beta=alpha+epsilon. It has order mp. The map

    sigma(a+b*epsilon)=a^p+alpha^(p-1)*b^p*epsilon      (F2)

is an algebra automorphism, sigma^d=1, and its restriction to K is
Frobenius. Indeed iterating multiplies epsilon by
alpha^((p-1)(1+p+...+p^(d-1)))=1. Its exact order is d because its
restriction already has that order.

Let u be the unique residue modulo mp satisfying

    u=p mod m,           u=1 mod p.                    (F3)

Then u is a unit, and beta^u=alpha^p+alpha^(p-1)*epsilon. Therefore

    sigma * multiplication(beta) * sigma^-1
       = multiplication(beta)^u.                     (F4)

This proves a clock-conjugating lift without killing the infinitesimal
coordinate. In contrast, the algebra's literal pth-power map sends
epsilon to zero and is not an automorphism. The identity is classical
finite-field/dual-number algebra applied to the present declared clock;
no general novelty claim is made for the underlying construction.

For the F243 factor, F3 gives u=124 modulo363. For BOTH factors in F1,
the semisimple order is242 and the nilpotent order is3. The compatible
CRT solution is

    u=487 mod726,   u=3 mod242,   u=1 mod3.             (F5)

Thus sigma:X->X^487 is an order-five automorphism of A and
sigma T sigma^-1=T^487. It acts trivially on the (X+1)^2 factor. On the
degree-five factor it sends alpha->alpha^3 and epsilon->alpha^2*epsilon.
The checker verifies every basis-product identity, invertibility,
sigma^5=1, and this literal twelve-dimensional matrix conjugacy.
Naive cubing instead has rank6, losing both nilpotent layers.

## F2. The actual response supplies the missing selector

Let b_r be the inherited response, and let H_ij=b_(i+j) modulo3,
0<=i,j<12. The inherited Hankel certificate gives det(H)=1. Hence

    y(r)=(b_r,...,b_(r+11))^T = H T^r e_0              (F6)

is an invertible chart from polynomial coordinates to the actual future
observation window. Store Phi=H sigma H^-1. It fixes the actual initial
window and sends y(r) to y(487r). All726 orbit states are checked.

However the original scalar selector ell=(1,0,...,0) is not fixed:

    b_2=0 mod3,             b_248=2 mod3,
    487*2=248 mod726.                                 (F7)

The exact transported selector is ell'=ell Phi^-1, so ell' Phi=ell.
This gives the commuting source/observation diagram. The receipt also
stores the differently directed pullback ell Phi; these two selectors
are not silently interchanged. Merely retaining a period or the abstract
field does not retain the original observation.

The native source here is a particular cyclic scalar-response module,
with a supplied polynomial chart and named initial/terminal maps. The
map is an algebra automorphism conjugating its multiplication clock;
the preserved predicate is the scalar response AFTER its selector is
transported. The lost information under naive cubing is six nilpotent
coordinates. Its required sidecars are H and the selector. No action on
the original tournament's vertices, all joint-port entries, or the
complement/re-marking language has been constructed.

## Connections and limits

This is a concrete 5/11/12 connection: degree5 Frobenius, order11^2 field
phase, and a12-dimensional observable module. The equality of periods
was already a hostile in the incoming work: forgetting only the five
F243 nilpotent coordinates leaves dimension7 and still has period726.
The present CRT lift repairs the clock action while explicitly recording
the observation change. It does not make the smaller quotient faithful.

The same mechanism explains why field Frobenius and full precision
transport must be distinguished in the F64/F256 and Ramanujan12/13
threads. Their actual operators and minimal polynomials differ; no
identification by matching integers is claimed. The Zagier/class-invariant
trace comparison remains in [ARITHMETIC_INTEGRATION](ARITHMETIC_INTEGRATION.md).

Compared with [TERNARY_EXTENSION_OBSTRUCTION](TERNARY_EXTENSION_OBSTRUCTION.md),
this is a successful lift on a different, explicitly supplied carrier.
There the native group relations prohibit every radical-reading first
step; here a CRT exponent gives an invertible clock lift, but fails to
preserve a fixed selector. Both tests keep the source map and the joint
operation laws. Neither supplies an arithmetic Sha finiteness argument.

Reproduce with `python -B periodic_frobenius_lift.py` and its -O variant.
The179 active checks include the full orbit and144 basis-product checks.
The JSON stores the incoming canonical source hash, every matrix, both
selector directions, and the literal mismatch. Normal/optimized receipts
agree. Next ask whether a legal native operation realizes Phi together
with the transported selector; the algebra construction alone does not.
