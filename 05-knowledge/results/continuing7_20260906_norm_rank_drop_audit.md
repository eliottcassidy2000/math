# Independent audit of the nearby-endpoint norm identity

**PASS: the all-parameter modular identity, rank-one unbounded family, and new carried-boundary consequence are accepted.** Two wording clarifications were applied before this final audit: the new quotient is smaller only when m<H, and the iff concerns the modular norm; rational nonvanishing is a one-way consequence. The mathematical formulas and producer source were unchanged.

The reviewed source is `continuing7_20260906_norm_rank_drop.py`, with its same-stem report and certificate. No producer implementation is imported. Parent owns repository filing and status promotion.

## 1. The actual norm and every residual rank

For positive integers H,z, set S=z+2H and use the monic falling-factorial rows

`Phi(t)=sum_(k=0)^H (H)_k(S)_(2k)t^(H-k)/(3k)!`,

`Psi(t)=Phi_(2H,2z)(t)`.

The actual characteristic norm is

`N_H(z)=(-1)^H Res(Phi,Psi)`.

Let p>6H be prime, and0<=m<2H with2S congruent2m+1 modulo p. Then the accepted identity is

`Psi modp=t^(2H-m)R_(H,m)(t)`,

where

`R_(H,m)(t)=sum_(k=0)^m (2H)_k(2m+1)_(2k)t^(m-k)/(3k)!`.

Both Phi and R have unit constant terms; Phi in fact has every coefficient a unit. The norm becomes

`N_H(z) modp=(-1)^H Phi_bar(0)^(2H-m)Res(R_(H,m),Phi_bar)`.

The residual resultant is a degree-m quotient determinant, with value1 for m0. This decreases the dimension when m<H. For H<=m<2H the exact identity remains valid, but there is no dimension reduction. The modular norm is nonzero iff this residual determinant is nonzero. That is sufficient for rational coprimality, without a converse and without a real-sign conclusion.

## 2. Independent proof of factorization and signs

The standard residue of S is(p+2m+1)/2. Since2m+1<4H<p, this lies between3H and p. The last2k integers ending at S avoid zero modulo p for every k<=H. All height and denominator factorials in Phi have index at most3H<p. This proves the unit statement, including its constant term.

In Psi, the upper product reduces to(2m+1)_(2k). For k<=m, all its factors are nonzero. For k>=m+1, the factor zero occurs. No factorial denominator meets p because3k<=6H<p. Its leading2H-m zero block and its entire residual factor are therefore exactly as displayed. In particular no singular denominator has been discarded.

Now use resultant multiplicativity. Since

`Res(Phi,t)=(-1)^H Phi(0)` and `Res(Phi,R)=(-1)^(Hm)Res(R,Phi)`,

the combined sign in the characteristic norm is

`(-1)^(H+H(2H-m)+Hm)=(-1)^H`.

The finite-field identity is thus correctly oriented. It is about the actual norm, not only a coefficient pattern. The unit constant term of Phi excludes the zero-root block from every common-root obstruction, so the gcd obstruction is exactly the one between Phi and R.

This argument specializes to the old monomial theorem when m0 and p>6H. It does not replace the old square-prime extension at4H<p<=6H, whose denominator ledger differs.

## 3. The exceptional rank-one value is never zero over Q

For m1, the residual is t+2H. Its sole obstruction is the auxiliary rational value

`A_H=Phi_(H,S=3/2)(-2H)`.

The half-integer endpoint need not correspond to positive z. It is a polynomial coefficient specialization used only modulo p. All its denominators have prime factors at most3H, or2, so reduction is legal at p>6H. Thus, within these rank-one endpoint hypotheses, the modular norm is a unit exactly when p does not divide the reduced numerator of A_H.

Write

`a_k=(H)_k(3/2)_(2k)/[(3k)!(2H)^k]`.

All terms are positive: after the first positive pair, the falling half-integers occur in pairs of negative factors. Moreover a0=1,a1=1/16, and for1<=k<H,

`a_(k+1)/a_k=(H-k)/(2H) * (2k-3/2)(2k-1/2)/[(3k+1)(3k+2)(3k+3)]`.

The first factor is less than1/2. The denominator minus the second numerator is

`27k^3+50k^2+37k+21/4>0`.

Hence the ratio is positive and less than1/2. The finite alternating sum therefore satisfies

`15/16<=A_H/(-2H)^H<1`.

The lower equality occurs only at H1; for H>=2 the tail can be paired into strictly positive differences. In particular A_H never vanishes. This proves the exceptional-prime set finite for each H; it does not infer the sign of N_H(z) at a positive integer parameter.

For every H there are primes p>6H outside that finite set. Choosing z=(p+3)/2-2H gives a positive integer, and all nonnegative translates z+lp have the same reduction. Since p is odd, infinitely many such translates have odd z. The inherited norm-jet iff then gives exact carried-boundary orders with unbounded complementary height H. This argument uses prime infinitude and the nonzero rational value, not a finite sample of primes.

## 4. Exact hostile and genuine new consumer

At H1, A_1=-15/8 has no numerator prime divisor greater than6, so there is no eligible rank-one exceptional prime. At H2,

`Phi_(S=3/2)(t)=t^2+t/4+1/640`,

`A_2=9601/640`.

The number9601 is prime. With z4798 and p9601, the original positive-integer rows obey the endpoint hypotheses, their modular gcd is exactly t+4, and their norm is zero modulo9601. Thus this is a genuine minimal-height failure of automatic rank-one units, not merely a failed proof template.

The same rational rows have gcd1 and norm20135 modulo28807. Consequently their rational norm is nonzero. A failed endpoint prime does not imply rational cancellation. Periodic lifts can also carry the bad reduction to odd z, but the displayed good-prime value20135 is claimed only for the displayed z4798 rows.

The exact additional carried example is

`(h,r)=(11,4), H7,z9,p43, N_7(9)=16 mod43`.

The dictionary gives4h-1=43, while the previous endpoint4h+1=45 has only prime divisors3 and5, both below the old threshold4H. The inherited `continuing5_20260906_complementary_norm_jets.md` therefore supplies exact order3 at x=-4. Its full leading-jet normalization and inverse-carry proof remain dependencies. This unit does not determine the real leading-jet sign.

## 5. Independent exact computation

The referee constructs coefficients by successive lowering recurrences rather than direct factorials. It computes full and residual resultants through Euclidean polynomial division, separately checks their monic gcds, and compares every field of the frozen certificate. The producer instead uses full quotient multiplication determinants.

- All172 producer cases are reconstructed from the stated H/m/prime selection, including all periodic rank-one lifts.
- A separate complete residual-rank range H1..10,m0..2H-1 supplies110 controls. This explicitly checks m>=H as identities without a dimension decrease.
- Twenty literal rational Sylvester computations cover every m0..2H-1 at H1..4. Their reductions check the full/residual determinant identity and both orientation signs by a third path.
- All80 exact A_H values and their normalized alternating bounds are checked, including every consecutive ratio and the lower equality only at H1.
- An independent prime sieve certifies9601 and28807; actual polynomial gcds and resultants reproduce the hostile, its good-prime repair, and the new boundary(11,4).

All6743 gates remain active under optimization. Normal and optimized output is byte-identical actual LF. No producer implementation or numerical root approximation is used.

## 6. Reproduction and reviewed pins

```
python continuing7_20260906_norm_rank_drop_audit.py
python -O continuing7_20260906_norm_rank_drop_audit.py
```

Filed layout prefers the adjacent or sibling-results certificate. The outside fallback is `continuing7_20260906_norm_rank_drop_certificate.json`. The certificate is used for comparison, not as a proof oracle.

Reviewed producer pins:

- Report:`69d1bed387e2a123d439947090badacfd16a81f1e7cfe5cc5398be7fccda77e6`.
- Source:`cf78d34048d1920d258a5537e7ddda6639999800b43c9e4d1b623f52dbbd5792`.
- Output:`94baac5cd4c4664a60d12897b4a1b18b8c96dd15b0d0cc34bd4668204aad5caa`.
- Certificate:`609fd8a81d29d139922980b0197614877056bdcde63c80c1d863cad6ff87f35a`.

Independent referee pins:

- Source:`3b9ec6b4d56f72fb82d728ee98f8be02871b2e17877ef3a4be9bf7cf862ae4f3`.
- Normal and optimized output:`240928154896412f97d61fa775c1d5f58000b385ee494579abccc9370b0577fa`.

Promotion basis: the coefficient factorization, exact norm orientation, legal local-ring reduction, all-rank scope, exceptional-value proof, unbounded-family quantifiers and actual boundary consumer all pass. The two phrasing clarifications are reflected in the reviewed report; no mathematical source repair was required. Parent may promote to independently audited status while retaining this candidate-report pin as lineage.
