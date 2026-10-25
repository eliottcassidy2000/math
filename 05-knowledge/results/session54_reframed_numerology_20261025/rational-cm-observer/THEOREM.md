# A uniform finite-conductor obstruction for rational cyclotomic CM observers

Status: **PROVED from standard cited CM theory**, with an elementary exact supplement. This is a reusable formulation of the mechanism in the tangent-Legendre classification; no claim of literature novelty or bound optimality is made.

## Theorem

Let R=A/B in Q(X) be a nonconstant rational function written with relatively prime A,B in Q[X]. Let its map degree be

    d=max(deg A,deg B)>=1.

For n>=3 put t_n=2cos(2pi/n) and F_n=Q(t_n). Assume B(t_n)!=0. If j=R(t_n) is a singular modulus (the j-invariant of an elliptic curve with complex multiplication by any imaginary quadratic order, not necessarily maximal), then

    |Gal(F_n/Q)^2| <= d.                              (1)

The superscript denotes the subgroup of squares. In particular, if p^a divides n for an odd prime p, then

    p^(a-1)(p-1)/gcd(4,p-1) <= d.                     (2)

If 2^a divides n and a>=3, then

    2^(a-3) <= d.                                    (3)

Define e_2(d)=3+floor(log_2 d). For each odd prime p<=4d+1, let e_p(d) be the largest positive integer a satisfying (2), or zero if no such a exists. Set

    N_d = 2^e_2(d) product_{odd primes p<=4d+1} p^e_p(d).

Then n divides N_d. The bound is independent of the coefficients of R. It is a necessary condition, not a characterization of CM or even of (1).

## Proof

### A. The low-degree observer loses at most d degrees

The real cyclotomic field F=F_n is Galois and abelian over Q, with group

    G=(Z/nZ)^*/{+1,-1}.

Let L=Q(j), a subfield of F. The element t_n satisfies A(X)-jB(X)=0 over L. This polynomial has degree at most d and is nonzero: otherwise A and B would be proportional, making R constant. Therefore

    [F:L] <= d.

Reduced form fixes the genuine rational-map degree and avoids treating a removable zero of an unreduced denominator as a pole. The nonvanishing condition ensures j is a finite value.

### B. A totally real singular-modulus field is elementary abelian

Because L is contained in the totally real field F, every conjugate of j is real. For the imaginary quadratic order O associated to j, standard complex multiplication theory gives a bijection from proper invertible ideal classes a of O to the conjugates j(a), with

    conjugate(j(a))=j(a^-1).

The class parametrization is injective. Thus reality of every j(a) forces a=a^-1 for every class, so Cl(O) has exponent at most two.

The ring class field H=K(j), K=Q(sqrt(Delta)), has Gal(H/K)=Cl(O), and Gal(H/Q) is its semidirect product with complex conjugation acting by inversion. In the exponent-two case inversion is trivial, and H/Q has elementary abelian group. Its subfield L/Q consequently also has elementary abelian Galois group. These statements apply to proper invertible ideal classes of arbitrary imaginary quadratic orders, so no fundamental-discriminant restriction is present.

Let N=Gal(F/L). Restriction gives G/N=Gal(L/Q), whose exponent is at most two. Hence g^2 belongs to N for every g in G, and

    |G^2| <= |N| = [F:L] <= d.

This proves (1).

### C. Prime-power quotients bound the conductor

For any prime power q dividing n, reduction U(n)=(Z/nZ)^* -> U(q) is surjective. To see this, lift a unit modulo q to a unit modulo the full power of its prime dividing n, and prescribe units at every other prime power via the Chinese remainder theorem. Reduction sends the subgroup {+1,-1} onto the corresponding subgroup modulo q, so it induces a surjection G_n -> G_q. A surjective group homomorphism maps the square subgroup onto the square subgroup. Therefore |G_q^2|<=|G_n^2|<=d.

For odd q=p^a, U(q) is cyclic of order p^(a-1)(p-1). Dividing by {+1,-1} gives a cyclic group of half that order. Its square subgroup has order

    p^(a-1)(p-1)/gcd(4,p-1),

which proves (2). For q=2^a, a>=3, the standard unit-group decomposition is C_2 x C_{2^(a-2)}, with -1 generating the first factor. Hence G_q is cyclic of order 2^(a-2), with square subgroup of order 2^(a-3), proving (3).

Equation (2) implies (p-1)/4<=d and therefore p<=4d+1. It also bounds every odd exponent a. Equation (3) bounds the 2-exponent by e_2(d). Every prime-power factor of n thus divides the corresponding factor of N_d, proving n|N_d. QED.

## Exact modest-degree examples

The certificate computes d=1,...,8 using integer arithmetic and proves maximality of each included exponent by checking the next exponent fails the relevant local inequality. In particular,

    N_1=120,
    N_2=240,
    N_3=65520=2^4 * 3^2 * 5 * 7 * 13.

The full output includes prime factorizations for all eight degrees. A unit-group check through n=500 is only an independent finite sanity test; the all-n conclusion comes from the proof above, not that cutoff.

For d=3 this gives a sharper coarse divisor bound than the valid 144-product superset used in the separately frozen tangent-Legendre classification. It changes neither that proof nor its checked 30-conductor universe.

## Positive and hostile controls

1. **Genuine d=1 CM control.** At n=5, t_n^2+t_n-1=0. The affine observer R(X)=85995X-52515 has value satisfying H_-15(J)=J^2+191025J-121287375, hence is CM. Exact substitution proves the polynomial identity. This checks that the denominator gcd(4,p-1) in (2) matters: replacing it by 2 would incorrectly reject p=5 for d=1. The polynomial H_-15 is pinned to the primary source below.
2. **Not sufficient for CM.** At n=3, R(X)=X/2 has degree one and yields j=-1/2. The group condition holds, and n divides N_1, but this value is not an algebraic integer and therefore is not CM.
3. **Divisor condition weaker than (1).** The conductor n=15 divides N_1=120, yet |G_15^2|=2>1. The independent prime-power bounds deliberately lose correlations between different prime factors.
4. **Constant maps must be excluded.** The constant function R=8000 has a CM value for every n, so no finite uniform conductor bound exists for constant observers. Arbitrarily large prime n already defeat such a conclusion. The degree-zero case is outside the theorem.
5. **Poles must be excluded.** For n=4, t_n=0 and R(X)=1/X has no finite value. Applying the CM statement there would be undefined, not a counterexample. Unreduced expressions such as X^2/X must first be canceled if interpreted as rational maps.

These are algebraic controls, not numerical CM recognition. The general observer need not preserve the tangent family's positivity: the affine n5 example has a negative conjugate, which is fully compatible with its CM order.

## Scope and provenance

This theorem concerns singular moduli as outputs of bounded-degree rational functions of real cyclotomic generators. It does not assert that every surviving conductor admits CM for a specified observer. It does not give an optimal N_d, enumerate all singular moduli, use an idoneal-number completeness claim, or imply a modular expression for the MRB constant or defect.

CM inputs: Cremona–Sutherland, *Computing the endomorphism ring of an elliptic curve over a number field*, Section 2.3, https://lucant.org/papers/2023/230123-Cremona.pdf . Used for the ring class field, inversion action, and reality criterion only. No completeness assertion about abelian-discriminant tables is imported.

Integrality and H_-15 control: Bilu–Luca–Pizarro-Madariaga, *Rational Products of Singular Moduli*, opening and Table 2, https://arxiv.org/pdf/1410.1806 . No global class-number catalogue is used.

SOURCES.json records precise bibliographic data and hashes of the privately inspected source copies, not runtime dependencies. No third-party PDFs are included. The elementary cyclotomic-field and unit-group facts used above are explicitly derived to the needed extent. Neither this supplement nor its verification mutates the frozen classification packets.
