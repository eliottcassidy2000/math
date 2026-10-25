# Independent audit: uniform rational-observer conductor bound

## Verdict

PASS. The statement in session54-rational-cm-observer/THEOREM.md is valid with its stated nonconstant, reduced rational-map and finite-value hypotheses. All imaginary quadratic orders, including nonmaximal orders, are covered. The explicit N_d is a uniform necessary divisor bound, not a sufficient CM condition and not equivalent to the global square-subgroup condition. No novelty or optimality claim is needed.

Audited producer manifest SHA256: 7f7e26f77074700b5b3adb40b575d99723e2faa64e8e8ad81bbad87fee58d223.

## Independent mathematical review

1. Write t=2cos(2pi/n), F=Q(t), j=R(t), L=Q(j). The nonzero polynomial A(X)−jB(X) has degree at most d and vanishes at t. Nonzero follows from nonconstancy, even when j is rational. The pole exclusion is necessary. Thus [F:L]<=d; leading-coefficient cancellation only improves the bound.
2. F is totally real and abelian. Every conjugate of j is real. Proper invertible ideal classes of the CM order parameterize its conjugates injectively, with complex conjugation inverting classes. Hence the class group has exponent at most2. The generalized-dihedral ring class group consequently becomes elementary abelian, as does its quotient corresponding to L. No maximal-order restriction or complete list of exponent-two class groups is used.
3. The quotient G/N=Gal(L/Q) kills all squares of G=Gal(F/Q). Therefore |G²|<=|N|=[F:L]<=d. This proves the central observer statement independent of its coefficients.
4. For every prime power q|n, unit reduction is surjective by lifting at the same prime and using CRT at other primes. It descends to quotienting by signs, and remains surjective on square subgroups. Odd q=p^a gives a cyclic quotient of order p^(a−1)(p−1)/2, whose square subgroup has order p^(a−1)(p−1)/gcd(4,p−1). In particular, replacing that denominator by2 is incorrect for p=5.
5. For a>=3, U(2^a)=<-1>×<5>, with the latter factor of order2^(a−2). Quotienting by signs and then squaring gives2^(a−3), including a=3. Smaller powers2 and4 are automatically included in e2(d)>=3. All odd primes satisfy p<=4d+1, and all exponent maxima are finite. Thus n|N_d follows for all n, not from a scan cutoff.
6. The n5 positive control uses exactly t=2cos(2pi/5), whose polynomial is t²+t−1: divide Phi5(z)=0 by z² and use z²+z^−2=t²−2. For a=85995,b=−52515, elimination gives (J−b)²+a(J−b)−a²=J²+191025J−121287375. Its discriminant is nonsquare. Thus both conjugate affine values are roots of the cited H−15, proving a genuine degree-one CM example and testing normalization.
7. Constant8000 is CM independently of n and cannot obey a finite uniform conductor restriction. R=1/X at n4 is undefined; unreduced X²/X must be canceled first. R=X/2 at n3 gives the rational noninteger−1/2 and is not CM despite passing the group bound. Finally15|N1 but |G15²|=2>1 explicitly separates the divisor bound from the stronger group condition.

## Independent exact computation

The audit does not merely rerun the producer's local formula. independent.py forms quotient classes {a,−a} explicitly and squares their representatives. Every admissible odd prime power q obeys q<=6d, since q<=4dp/(p−1); every admissible two-power obeys q<=8d. Therefore direct enumeration of prime powers q<=8d followed by an LCM independently reconstructs the entire local bound for d1 through8.

The resulting N_d are:

- d1:120
- d2:240
- d3:65520
- d4:2227680
- d5:122522400
- d6:122522400
- d7:3553149600
- d8:7106299200

The audit also derives the global square-group size by CRT: take the product of local unit-square-group sizes and divide by2 exactly when −1 is a square modulo n (for n>2). It agrees with direct quotient enumeration for every3<=n<=500. Enumerating every divisor of each now-proved finite N_d then gives complete global square-group survivor counts8,16,30,36,46,69,71,83. The maximum is840 for d6,d7,d8. This does not contradict the producer's explicitly bounded counts68,70,82 through500; the additional conductor is840. These are necessary-group survivors, not a classification of CM observer values.

Normal and Python -O output is identical. Eight independent hostile runs reject a changed N3, an omitted prime5 factor, a false bounded survivor count and a doubled affine coefficient (wrong normalization). The producer's own wrong-gcd, wrong-N3 and wrong-affine mutations also reject in both modes.

## Portability and provenance

The ten-file text-only producer packet was copied to an unrelated temporary directory. Its exact-file-set verifier, checker and mutation suite all replayed normally and under -O, matching the frozen outputs exactly. Corrupt, missing, extra and nested-extra files were rejected in both modes. See closed-replay.json. No third-party PDF or article text is included or required at runtime.

The primary theoretical inputs were already independently inspected during the associated classification audit:

- Cremona–Sutherland, Computing the endomorphism ring of an elliptic curve over a number field, sections2.2–2.3. Used only for arbitrary-order ring class structure and class inversion, not catalogue completeness. https://lucant.org/papers/2023/230123-Cremona.pdf . Inspected PDF SHA2562037856cfc69c7db0703ef3a18e34631767eea9171d6a05f41ba090cca470d7a.
- Bilu–Luca–Pizarro-Madariaga, Rational Products of Singular Moduli, opening and Table2. Used for algebraic integrality and the explicit H−15 polynomial, not global class-number lists. https://arxiv.org/pdf/1410.1806 . Inspected PDF SHA256976e516695aaadcfa0b7f9be4f92a23f035fb3f1b6c58a503021d4680b951363.

References and inspected hashes are provenance, not redistributed article contents. All earlier classification and audit packets were preserved. No repository edits were made.
