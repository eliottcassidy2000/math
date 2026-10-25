# The tangent-Legendre family has exactly one CM member

Status: **PROVED subject to independent audit of this new synthesis**, using standard cited CM results and a small exact finite certificate. The all-integer reduction is proved analytically; it is not extrapolation from a numerical scan.

## Statement and boundary

For an integer m>4 put k_m=tan(pi/m), lambda_m=k_m^2, and

    j_m=256(1-lambda_m+lambda_m^2)^3/[lambda_m^2(1-lambda_m)^2].

Then j_m is a singular modulus (equivalently, the Legendre elliptic curve has complex multiplication) **if and only if m=8**. At m=8, j_m=8000 and the CM order has discriminant -8, as already proved in Session52. The novelty here is uniqueness within the full integer tangent family, not discovery of that CM curve.

This says nothing about a modular expression for the actual MRB constant C or defect D=A(8)-C. The map from branch-jump amplitude A(m) to tan(pi/m) forgets the alternating functional. No theta-tail identity, Ogg-prime criterion, LRC consequence, or idoneal-list completeness is asserted.

Inheritance: Session52 `session52-mrb-supersingular/RESEARCH.md`, its `INHERITED_BOLZA_EXCERPT.md`, and THM-4012's Bolza elliptic quotient. Closest proved mechanism is m8 -> j8000; hostile controls m6 and m12 have rational nonintegral j. Corrected near miss is an approximate MRB defect, not equality with a nome power. Least-used sidecar is the Galois square subgroup of the real cyclotomic field.

Concept board: branch amplitude; Legendre parameter; real cyclotomic Galois group; CM class-group inversion; algebraic integrality; dominant singular modulus.

## 1. The cubic carrier

Set n=m/gcd(m,2), z=exp(4*pi*i/m), and t=z+z^-1=2cos(4*pi/m). Then z is a primitive n-th root, n>=3, and F=Q(t) is the maximal real subfield of Q(z), with

    Gal(F/Q)=G=(Z/nZ)^*/{+1,-1}.
    j=64(10+3t)^3/[(2-t)^2(2+t)].                    (1)

The cubic equation obtained by clearing (1) implies [F:Q(j)]<=3. It is a nonzero polynomial (the numerator and denominator are not proportional). Every conjugate t lies strictly between -2 and 2. Moreover

    64(10+3t)^3 - 1728(2-t)^2(2+t) =256(9t+14)^2.

Thus all conjugates of j are >1728: equality would require t=-14/9, impossible because t is an algebraic integer. This also makes every root smaller than the trace of its minimal polynomial when its degree exceeds one.

## 2. CM forces a small square subgroup

**CITED CM input.** For a singular modulus j of order O and imaginary quadratic field K, H=K(j) is its ring class field. Its Galois group over K is Cl(O), and complex conjugation acts by inversion. The degree [Q(j):Q] is the class number h(O); conjugates are indexed bijectively by proper ideal classes. See Cremona–Sutherland, Section 2.3, and Bilu–Luca–Pizarro-Madariaga, beginning and Proposition 2.5.

Here Q(j) is a subfield of the totally real abelian field F. All its conjugates are real, so the identity conjugate(j(a))=j(a^-1) and injectivity of the class parametrization force a=a^-1 for every ideal class. Therefore Cl(O) has exponent at most two. The ring class Galois description then implies Q(j)/Q is elementary abelian of exponent two. (Equivalently, H/Q is the direct product of this class group with complex conjugation.)

Let N=Gal(F/Q(j)). Since G/N has exponent two, its subgroup of squares G^2 is contained in N. Consequently

    |G^2| <= |N|=[F:Q(j)] <=3.                       (2)

This is not a catalogue claim about exponent-two quadratic class groups.

## 3. A genuinely finite conductor universe

Write U=(Z/nZ)^*. Its square subgroup maps onto G^2 with kernel of size at most two. Thus |U^2| belongs to {1,2,3,4,6}. For an odd prime power p^a exactly dividing n,

    |((Z/p^aZ)^*)^2| = p^(a-1)(p-1)/2.

For 2^a, the corresponding factor is 1 for a<=3 and 2^(a-3) for a>=3. Chinese remaindering multiplies these factors. It follows that

    n=2^a 3^b 5^c 7^d 13^e,
    0<=a<=5, 0<=b<=2, c,d,e in {0,1}.

The apparent prime 11 is excluded because its square-group factor 5 cannot divide any of {1,2,3,4,6}. Larger primes have factor at least eight. This gives 144 products, including 1 and 2 which are discarded, not an unbounded search assumption.

The exact group screen (2) leaves precisely 30 conductors:

    3,4,5,6,7,8,9,10,12,13,14,15,16,18,20,21,24,26,
    28,30,36,40,42,48,56,60,72,84,120,168.

For EACH of these 30, exact elimination shows [Q(j):Q]=phi(n)/2=[F:Q]. The certificate first reconstructs the real cyclotomic polynomial Psi_n from Phi_n, then computes the resultant with the cubic (1), removes repeated factors, and verifies irreducibility. It also independently substitutes j(t) modulo Psi_n. Therefore N is trivial for each remaining conductor. Equation (2) strengthens to G^2={1}, leaving

    n in {3,4,5,6,8,10,12,24}.

Algebraic integrality, necessary for CM, eliminates n=3,5,6,10. The remaining minimal polynomials are:

    n=4:  J-8000.
    n=8:  J^2-159616J+282300416.
    n=12: J^2-840064J+1593413632.
    n=24: J^4-13839104J^3+365636968448J^2
          -2380381457219584J+4567157201705107456.

## 4. Three small discriminant exclusions

**CITED bound.** Bilu–Luca–Pizarro-Madariaga Proposition 2.1 states

    ||j(tau)|-exp(2*pi*Im(tau))| <=2079

on the standard fundamental domain. Every imaginary quadratic order of discriminant Delta has a principal reduced form with a=1, whose j-value is a conjugate and whose imaginary part is sqrt(|Delta|)/2. Hence, if all conjugates are positive with trace T,

    exp(pi*sqrt(|Delta|))-2079 < T.                  (3)

The certificate verifies these strict lower bounds using only rational arithmetic: pi>157/50, rational lower bounds 3.87<sqrt(15), 4.35<sqrt(19), 5.65<sqrt(32), and the first 100 nonnegative terms of exp's Taylor series.

- n=8 has T=159616, so (3) gives |Delta|<15. Every negative discriminant in that range has class number one, whereas this j has degree two.
- n=12 has T=840064, so |Delta|<19. The sole class-number-two discriminant in this range is -15. Its polynomial J^2+191025J-121287375 (Bilu et al., Table 2) is coprime to the displayed n12 polynomial. Hence impossible.
- n=24 has T=13839104, so |Delta|<32. Every negative discriminant in this range has class number at most three, whereas this j has degree four.

The small class-number statements are independently checked by enumerating ALL primitive reduced forms, with a<=sqrt(|Delta|/3) and the standard boundary convention. No global class-number table, idoneal catalogue, or GRH enters.

Only n=4 remains, corresponding to m=8. The inherited AGM/lattice argument gives its CM positively. QED.

## 5. Exact requested box and reproducibility

For m=5..24, j is integral only at m=8,15,16,21,24. Integrality is not sufficient: m15 has cyclic-quartic Galois obstruction; m21 has a degree-six obstruction; m16 and m24 fail the discriminant tests above. The all-m theorem also needs and excludes m48, which was outside the requested initial box.

Hostile controls: j6=21952/9 and j12=140608/3. Positive control: j8=8000. There are no numerical near-integer tests. The script uses explicit exceptions, not Python assert, so checks remain active under -O.

    python check.py > results.normal.json
    python -O check.py > results.optimized.json
    cmp results.normal.json results.optimized.json

Dependencies: Python standard library and SymPy. Exact rational/integer polynomial arithmetic only. The certificate's finite universe is justified in Section 3, rather than chosen by a convenient cutoff.

## Sources and scope

1. J. E. Cremona and A. V. Sutherland, *Computing the endomorphism ring of an elliptic curve over a number field*, Section 2.3, https://lucant.org/papers/2023/230123-Cremona.pdf . Used ONLY for CM class-field structure, inversion, and reality criterion; no statement about completeness of abelian-discriminant lists is imported.
2. Y. Bilu, F. Luca, A. Pizarro-Madariaga, *Rational Products of Singular Moduli*, https://arxiv.org/pdf/1410.1806 , opening degree/integrality fact, Proposition 2.1 (2079 bound), Proposition 2.5 (reduced forms and conjugates), Table 2 (H_-15). No global list completeness is required.
3. Session52 `session52-mrb-supersingular/RESEARCH.md`: positive CM(-8) identification and MRB/Bolza lineage. That proof uses DLMF 19.8.5 and 23.15 for elliptic-function conventions.

The connection map preserves a Legendre elliptic isomorphism class and exposes a cubic cyclotomic carrier. It destroys the analytic MRB sum, so those data must be retained as a separate sidecar for any future MRB claim. The general reusable mechanism is: a low-degree rational observer of a real cyclotomic generator, if CM, forces a bounded Galois square subgroup. Finite exact rank/degree checks can then upgrade a local numerical pattern to an all-parameter theorem.
