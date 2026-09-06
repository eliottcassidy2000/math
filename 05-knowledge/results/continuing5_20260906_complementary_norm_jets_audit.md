# Independent audit: complementary norms and the exact carried jet

**Status: PASS — analytic proof and independent exact controls. No repair requested.** The general leading-coefficient identity and its exact-order iff are accepted for all h>=1 and1<=r<=h. The corollary gives seven exact boundary diagonals0<=h-r<=6 at arbitrary height. It does not prove an all-height actual response-sign or endpoint theorem.

I read the full producer proof and source, the incoming sharp endpoint jet in `third_20260906_trace.md`, and the reflection/regular-family proof in `continuing4_20260906_regular_duality.md`. The reflected identities and their complete residual certificates were independently audited previously. The present proof uses only positivity of the regular norm polynomial from that certificate, so it does not need the separately cited real-root supplier or the primitive support-return gcd condition.

## Analytic check, including the if-and-only-if boundary

Use exactly the producer's monic p, full Laurent q and characteristic norm c_(h,h). Put H=h-r, z=2r+1, delta=x+r, and

```
a0=(2h+1)! H!/[(3H)!(2r+1)!],
A=(-1)^(r-1)(2h+1)! H!(r-1)!/(3h)!,
B=2(2H)!(2r)!/(6h+3)!.
```

I independently differentiated the falling products: the unique zero in the constant p coefficient occurs at index H, leaving H! times(-1)^(r-1)(r-1)!; the unique zero in the carry occurs at index2H, leaving(2H)!(2r)! with derivative2. Thus A is nonzero of sign(-1)^(r-1), B is positive, and the regular complement constant a0 is positive.

At delta=0, p=t^r Phi with Phi(0)=a0. These two factors are coprime even when Phi itself has repeated roots. Formal lifting separates the two invariant blocks. On putting delta=epsilon^r, the r small roots have t_i=epsilon v_i(epsilon), with nonzero distinct leading values solving

`a0 v^r+A=0`, hence `product v_i(0)=(-1)^r A/a0`.

The only leading small-root contribution is the inverse carry. Its coefficient is B delta+O(delta^2), so its value has epsilon order r-1. Every nonnegative channel e<2r is divisible by delta and has order at least r+e; every channel e>=2r has order at least e>=2r. This strict separation remains true for r=1, when the carry has nonzero order-zero response. Therefore the characteristic norm of the small block is

`B^r a0/A *delta^(r-1)+O(delta^r)`.

This initially follows to higher Puiseux order. The canceled full response is already polynomial in the quotient; the lifted block and its determinant lie in C[[delta]]. Integer exponents then justify precisely O(delta^r). The argument does not set the vanishing raw carry to zero before cancellation.

At the complementary block, t is a unit and the raw specialized q equals

`t^(2r) Psi(t)/(4h+2)!`.

Consequently its characteristic norm at delta=0 is

`a0^(2r) N_H(2r+1)/((4h+2)!)^H`.

Here the even power removes the sign in product rho=(-1)^H a0. This computation is a determinant identity with multiplicities, requiring neither simple complement roots nor N_H nonzero. Multiplying the two blocks yields

```
c_(h,h)(-r+delta)
 = [B^r a0^(2r+1) N_H(2r+1)/(A ((4h+2)!)^H)]
     *delta^(r-1)+O(delta^r).
```

Every prefactor apart from N_H is nonzero. Thus the order equals r-1 **if and only if** N_H(2r+1) is nonzero. If it vanishes, the order is at least r, with infinity allowed for an identically zero polynomial. No unsupported nonvanishing is inserted in the general theorem.

For H=0, Phi=Psi=N_0=1 and a0=1; this recovers the inherited all-height x=-h jet. For H=1,...,6, the audited regular norm is the product of positive paired factors and a polynomial with positive rational coefficients for every positive real z. Therefore the leading coefficient is nonzero of sign(-1)^(r-1) on all seven diagonals. H>=7 remains conditional on its own complementary nonvanishing. The six additional diagonals are an unbounded consequence of retaining the complement, not additional finite endpoint certifications.

## Independent exact path

The referee imports no producer code and does not use either previous norm certificate as an oracle for the new old-norm computations. For h=1,2,3,4 it constructs the literal rational p and t*q at every integer x from1 to2h^2+1, then computes their Sylvester resultant. Since p is monic,

`Res(p,tq)/p(0)=(-1)^h product_(p(rho)=0) q(rho)=c_(h,h)`.

This convention is checked on linear and quadratic controls. Positive interpolation parameters keep p(0) nonzero, avoiding singular evaluation. The inherited quotient degree bound deg_x c_(h,h)<=2h^2 makes these64 values sufficient to reconstruct the entire norm polynomials by Newton forward differences. Four additional noninteger resultants at x=1/2 agree independently.

The ten boundary pairs1<=r<=h<=4 then have every lower jet coefficient checked zero and their leading coefficients compared exactly with the proposed formula. Complementary N_H values in this comparison are computed by fresh rational Sylvester resultants of the regular Phi,Psi, rather than taken from coefficient tables. Both complete reflected rows are also checked at every pair, with the inverse-carry position preserved. The h1,r1 characteristic norm is1/90720; the response is-1/90720, while naive specialization at its zero-root block would give zero. That is the retained normalization hostile.

Direct product differentiation checks all78 pairs1<=r<=h<=12 against A,B,a0, without substituting the factorial expressions as their own oracle. These are finite controls of the displayed analytic formulas, not evidence used to extrapolate their quantifiers.

Finally, the referee reads the pinned previously audited regular certificate and checks every residual **norm** coefficient for H1..6 positive, with its complete degree and paired-factor reconstruction. It compares those data with each independently computed small complementary norm. This is a dependency replay, separate from the fresh old-norm path; it does not claim a second new proof of every high-degree regular certificate.

## Freeze and promotion basis

```
python ../../04-computation/continuing5_20260906_complementary_norm_jets_audit.py
python -O ../../04-computation/continuing5_20260906_complementary_norm_jets_audit.py
```

Both pass706 always-active exact gates, with byte-identical actual LF output. The source configures LF and capture rejects carriage returns instead of normalizing. No repository or Git change was made.

- Referee source SHA256:`c5a3f723ae99bef4ec825411ef1ca9a48e7a5afc8f9bf9e0581b76d4a6cecb61`.
- Referee output SHA256:`9330b4c4526f3709d18f8a2c02e9b0604a47329527e5e2c6a8e8061fbe71de6a`.
- Reviewed producer source:`d1be69ab693d9cd57c2b6e269634681d6e87fcb425588e6489310d93991ec6a5`.
- Reviewed candidate proof:`467868c809a1e78b3b812997316ec67c99acfe67de2508d608ead623771f2a01`.
- Inherited regular certificate:`0d5a65f03fc4f4295f3db38bca8609375cfea8805f21499a28e9a5e0d9a1ccd4`.

The promotion basis is the complete all-parameter two-block proof, the exact-order iff with its vanishing boundary, and a distinct rational-resultant reconstruction of all low-height norms and jets. Parent owns filing and status promotion. The source of the six new diagonals is the already-audited positive complementary norm, with no additional real-root or coprimality assumption.

Filed checkpoint provenance: the [raw-byte manifest](continuing5_20260906_manifest.json)
pins the final report, source and output. Reviewed candidate-report hashes
above refer to the pre-promotion bytes. Source-location defaults and local
links were made portable where necessary; all emitted outputs were replayed
as raw bytes. The independent audit supplies the stated promotion basis.
