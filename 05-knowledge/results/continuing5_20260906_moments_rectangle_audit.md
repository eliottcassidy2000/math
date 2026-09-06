# Independent audit: the coefficient rectangle is closed at every original phase

**Status: PASS — full analytic proof and complete independent exact certificate audit. No mathematical repair requested.** The producer may be promoted from candidate to PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED within its stated model scope.

Audited producer: `continuing5_20260906_moments_rectangle.md`, source, transcript and complete JSON in this directory. The theorem is: for every real83<=x<=85,34<=y<=36,0<=z<=5, every positive zero s of the original first polynomial has **the full carried Q(-s)<0**. If B has five nonnegative roots in this x,y rectangle, its fifth elementary coefficient automatically satisfies0<=z<540/109<5. Neither C nor D interlacing is needed. The full two-anchor model outside this rectangle remains OPEN; no neighboring coefficient is asserted to be an actual factorial row.

## 1. Root geometry and all weak boundaries

I checked the first-critical-point argument independently. A monic quintic B with five nonnegative roots has four nonnegative derivative roots, counted with multiplicity. B'(0)=y>0 excludes a derivative root at zero. The first derivative root eta therefore lies between the first two roots of B, including coincidence. The sign of B between its first two ordered roots is nonnegative: there are four nonpositive factors and one nonnegative factor; if the interval collapses, B=0 there. This remains valid with a zero first root or repeated roots.

The maximum derivative at8/25 in the rectangle is -146524/78125<0. Consequently eta lies strictly between0 and8/25. With F=B+z,

`z<=F(eta)<36 eta-(327/5)eta^2<=540/109<5`.

Strictness comes from eta>0 and eta<13 in the term eta^5-13eta^4; it does not assume simple beta roots. Nonnegative roots give z>=0 independently. No interlacer, truncated Hankel positivity, or z=0 admissibility hypothesis enters this proof. In particular it does not silently import the false general assertion that every fixed-x,y admissible fibre reaches zero.

## 2. Original phases are exhausted

The literal first row is

`P(-s)=2002[z s^4-(12/7)y s^3+x s^2-10s+1/11]`.

All56 endpoint-corner signs were independently reconstructed directly from the original O-star-beta row. Affinity in x,y,z extends each sign to the whole closed prism. Three disjoint sign changes place roots in(1/102,1/100), (1/9,13/100), (1,8/5). The negative sign at19/2 and a positive quartic leading coefficient supply a fourth root beyond19/2 when z>0. Degree exhaustion then proves all roots are positive and simple, with exactly one in each window. When z=0 the cubic has exactly the first three simple positive roots. Thus the escaping phase is covered, and no repeated first-row phase is hidden at a coefficient boundary.

The degree argument applies to first-row phases, not to the roots of B. Confusing those two polynomials would lose precisely the weak beta boundaries this theorem includes.

## 3. Full response and independent certificate method

The referee reconstructs both ordinary convolutions and both Laurent shifts numerically in x,y,z. In particular it retains O^2+t^-1 E^2, beta^2+2t C_raw D_raw, exponents -1 through8, and q_-1=28. It does not import the producer or its sparse polynomial routines.

For raw Q, the separate variable degrees are at most2. Exact tensor interpolation from x,y,z in{0,1,2} therefore reconstructs every monomial coefficient, which is compared with the complete producer JSON. All raw monomial coefficients are positive. This independently gives the first-branch signed monomial bound

`-728550005046322718853208807/1704830652000000000000000 < -427`

for sQ(-s), without eliminating z or dividing by a small root.

For the remaining branches, the referee first polarizes Q's quadratic dependence on z from its values at0,1,-1, then substitutes the exact original-root z expression. At the nine x,y points in{83,84,85} times{34,35,36}, it reconstructs all Laurent coefficients of h. Every negative s power cancels, and every surviving exponent is between0 and8. These are complete identities, not numerical samples: the eliminated coefficients have separate x,y degrees at most2, since the raw expression is quadratic and the substituted z is affine in x,y. Thus the nine-point tensor grid also proves the cancellations identically. The recovered h agrees with every coefficient in the producer's manifest, and27 further literal same-root evaluations check `sQ(-s)=-(14/11)h(s)` directly.

The transformations use separate polynomial convolution and Horner translation. Instead of the producer's monomial-to-Bernstein formula, the referee obtains each degree-two Bernstein vector from its values at0,1/2,1:

`(b0,b1,b2)=(f(0),2f(1/2)-(f(0)+f(1))/2,f(1))`.

Applying this independently in x and y reconstructs all243 tensor coefficients. They agree exactly with the manifest and are all strictly positive. Their minima are

- [1/9,13/100]:10396002758034286225411/24500000000000000;
- [1,8/5]:1992478168568/49;
- [19/2,infinity):165789872820/49.

The compact transforms are(1+u)^8 h((a+bu)/(1+u)); the unbounded transform is h(19/2+u). Each coefficient in u is a strictly positive Bernstein polynomial on the full closed rectangle. This proves positivity for every u>=0. The finite right endpoint is included by the positive leading coefficient, so no projective endpoint is omitted. At every retained original root, s>0 then gives Q(-s)<0.

## 4. Controls and stopping scope

For every x,y in the rectangle, z=1 has B signs -,+,-,+,-,+ at0,1/10,1,3,5,7. The referee checks all24 corner signs; coefficient affinity and degree exhaustion prove five simple positive beta roots throughout this section. This is genuine two-dimensional nonvacuity of beta geometry, not a claim that both interlacers hold throughout it.

The actual central support(-27,1,15) is verified by a separate charge-solving enumeration: with c positive15 terms, the negative27 count is a=(mass+14c)/28. Integral nonnegative compositions at masses14 and28 recover all five and ten channels through binomial products. This checks the ordinary core and doubled lower carry together.

Both stated hostiles were independently retained. Formal moment coefficients are reconstructed by the finite geometric series for the reciprocal denominator, and their5x5 determinants use the permutation formula rather than the producer's elimination. The degree-seven positive-response surrogate has its original first equation exactly zero and positive full Q, while both degree-eight Hankel determinants are negative. At x84,y35,z6 the original first polynomial changes sign in(16693/2000,41733/5000); all nine coefficients of the transformed h are negative there, proving positive response at its original root. The D moment determinant is exactly-25668. These hostiles show why neither root constraints nor the fifth-coefficient cap can be discarded outside the stated prism.

The proof includes all x,y endpoints, z=0 and5, and the whole unbounded original-phase branch. Actual beta discriminant boundaries lie inside the strict z cap and need no extra limiting argument. No universal free-x,y or all-support Laurent noncancellation claim follows.

## 5. Reproduction, pins and promotion basis

```
python ../../04-computation/continuing5_20260906_moments_rectangle_audit.py
python -O ../../04-computation/continuing5_20260906_moments_rectangle_audit.py
```

Both runs pass566 always-active exact gates with byte-identical actual LF stdout. The source configures LF, and capture rejects carriage returns instead of normalizing them. The full proof was read; finite computations certify complete polynomial identities and positivity coefficients using proved degree bounds, rather than purporting to sample a continuous theorem.

Audited producer pins:

- Source:`da6236b1828b05bf5d30fe09d97e95005db36fe1baaa390c1b422aba85f13a25`.
- Output:`4aab430109b707da57cd2654ee0138d1b141002a1d67d33164cbe66ba0bb4a28`.
- JSON:`1e86411634a8f62f98a512d741d4b8df8ba5a069fb262d92bafe64e40569a4f8`.
- Candidate report:`fcfa0824a6494444998f29d0e3245aea3f7c5d05b364dedd0d68979d383c5cb3`.

Independent referee pins:

- Source:`83febb8f0dd76e89711a6d5f880104ec6bcb9895769fa09e0dae7deba2096d64`.
- Output:`d3d9e185bf98fd3d5a9e918dcea0d7e812f741856a96cfac51bcf20ab6d7cb56`.

Promotion is supported by the all-parameter proof, complete independent coefficient reconstruction, actual factorial-channel controls and explicit domain hostiles. This audit leaves producer bytes unchanged and makes no repository or Git changes. Parent owns filing and status promotion.

Filed checkpoint provenance: the [raw-byte manifest](continuing5_20260906_manifest.json)
pins the final report, source and output. Reviewed candidate-report hashes
above refer to the pre-promotion bytes. Source-location defaults and local
links were made portable where necessary; all emitted outputs were replayed
as raw bytes. The independent audit supplies the stated promotion basis.
