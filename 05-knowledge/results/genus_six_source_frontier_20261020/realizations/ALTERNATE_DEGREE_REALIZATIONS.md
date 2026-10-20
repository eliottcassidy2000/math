# Alternate degree presentations at conductor twelve

20 October 2026. **PROVED by explicit identities and degree-separated module bases; exact symbolic checks supplied.**

## Results

Over every characteristic-zero field k, both candidate alternate characteristic sequences are realized by immersed primitive polynomial subalgebras of k[t]:

- (4,10,5), with degree semigroup <4,5>;
- (6,9,4), with degree semigroup <4,6,9>.

Both have codimension six and numerical-semigroup conductor twelve. These are polynomial realizations, not deductions of realizability from abstract characteristic axioms. They do not establish any new Keller-source conductor bound or any Jacobian-conjecture claim.

Here “immersed” means gcd(p',q')=1, or a nonvanishing derivative pair. It does not mean injectivity or a locally closed scheme embedding. The alternate presentations refer to the same numerical-semigroup type; no isomorphism between the particular displayed algebras is asserted.

The inheritance is the full degree-separated free-module certificate used in the frozen conductor-ten source-exclusion packet. The hostile is the already-realized alternate (6,9,2), which prevents declaring an alternate characteristic sequence unrealizable merely because its smallest generator is absent from the original coordinate pair. The anchor here is exact realization; the niche is a parameter family retaining immersion; the wildcard is the degeneration to (6,9,2). The live concepts are cancellation remainders, characteristic gcd drops, integral power bases, distinct residue degrees, and derivative common zeros.

## 1. The (4,10,5) presentation

Set

    p = t^4+t,
    q = t^10+(5/2)t^7+(15/8)t^4+(5/16)t,
    r = q^2-p^5-(5/64)p^2.

Direct expansion gives

    r = t^5/64+5t^2/256.                         (Q)

Since p and r have coprime degrees, k(p,q)=k(t). The justification is the field-tower divisibility argument: [k(t):k(p,q)] divides the degree of every nonconstant polynomial in k[p,q], hence divides both 4 and 5.

The integral power basis 1,q,q^2,q^3 over k[p] changes unimodularly to

    1,q,r,qr,

because r is monic of degree two in q with coefficients in k[p]. Their degrees are

    0,10,5,15,

which represent the four distinct residue classes modulo four. No leading-term cancellation is possible in a nonzero k[p]-linear combination. Thus the exact degree set is

    4N union (10+4N) union (5+4N) union (15+4N) = <4,5>.

The missing degrees are 1,2,3,6,7,11, so the codimension is six and conductor is twelve. Equivalently, the Apéry sum is (0+10+5+15-0-1-2-3)/4=6.

The immersion certificate is especially short:

    q' mod p' = -5/8.

Consequently gcd(p',q')=1 over k, and also after any extension of k.

A one-parameter version, valid for b != 0, is

    p_b = t^4+b t,
    q_b = t^10+(5/2)b t^7+(15/8)b^2 t^4+(5/16)b^3 t,
    r_b = q_b^2-p_b^5-(5/64)b^4 p_b^2
        = (b^5/64)t^5+(5b^6/256)t^2.

Here q_b' mod p_b' = -5b^3/8. At b=0 the algebra is k[t^4,t^10] and is not primitive; that boundary must not be included.

## 2. The (6,9,4) presentation

Set

    p = t^6+t^4,
    q = t^9+(3/2)t^7+(3/8)t^5-(5/32)t^3-(3/128)t,
    r = q^2-p^3+(3/16)p^2-(7/1024)p.

Direct expansion gives

    r = t^4/2048+9t^2/16384.                    (S)

The field index [k(t):k(p,q)] divides 6,9,4 and is therefore one. The integral power basis 1,q,...,q^5 over k[p] changes unimodularly to

    1,q,r,qr,r^2,qr^2.

Its degrees are

    0,9,4,13,8,17,

which are distinct modulo six. Hence the exact degree set is

    6N union (9+6N) union (4+6N) union (13+6N)
       union (8+6N) union (17+6N) = <4,6,9>.

For example, reduce the exponent of 9 modulo two using 2*9=3*6, and the exponent of 4 modulo three using 3*4=2*6, to see the displayed union is the whole monoid. The gaps are 1,2,3,5,7,11. The codimension is six, also given by (0+9+4+13+8+17-15)/6=6, and the conductor is twelve.

For immersion,

    p' = 2t^3(3t^2+2),
    q'(0) = -3/128,
    q' mod (3t^2+2) = -27/128.

Neither zero nor a root of 3t^2+2 is a common derivative zero. Thus gcd(p',q')=1, including after extending k.

## 3. A two-parameter family and its exact boundaries

Let a,b belong to k and put D=a^2-4b. Define

    p = t^6+a t^4+b t^2,
    q = t^9+(3a/2)t^7+(3b/2+3a^2/8)t^5
          +a(-5a^2+36b)t^3/32+3(16b^2-a^4)t/128,
    A = -3aD/16,
    B = D^2(7a^2+8b)/1024.

The exact identity is

    q^2-p^3-Ap^2-Bp
       = D^3 t^2(9a^2+8a t^2-4b)/16384.        (F)

If aD != 0, the remainder has degree four and the same module proof gives the primitive degree semigroup <4,6,9>, with codimension six. The explicit example above is a=1,b=0.

The derivative condition is also exact. Write z=t^2 and

    G(z)=3z^2+2az+b,
    Q_1(z)=q'(t).

Then p'=2tG(t^2), and

    Q_1(0) = -3D(a^2+4b)/128,
    Q_1(z) mod G(z) = 3D(-a^2+12az+12b)/128,
    Res_z(Q_1,G) = 6561D^4/16384.

Thus the family is immersed if and only if D(a^2+4b) != 0. On the genus-six stratum aD != 0, the only non-immersed divisor is a^2+4b=0; the common derivative zero there is t=0.

The two other boundaries explain why the parameter restrictions matter:

- If a=0 and b != 0, formula (F) becomes r=b^4 t^2/64. Hence t^2 belongs to the algebra, and

      k[p,q] = k[t^2, t^9+(3b/2)t^5+(3b^2/8)t].

  The even/odd direct sum over k[t^2] gives exactly <2,9> and codimension four. This is the inherited alternate (6,9,2) mechanism. For b=4 it is the earlier f=t^6+4t^2, g=t^9+6t^5+6t.

- If D=0, put h=t^3+(a/2)t. Then p=h^2 and q=h^3. Consequently k(p,q)=k(h), of index three in k(t), and the algebra is not primitive. Its raw degree set must not be interpreted as a numerical semigroup of k[t].

When aD != 0 and a^2+4b=0, polynomial realizability remains true, but immersion fails. Therefore realizability and immersion should continue to be recorded as separate properties.

## 4. How the sextic family is forced inside the parity ansatz

This is a small classification within the displayed even/odd ansatz, not a normal-form classification of all polynomial degree-(6,9) pairs.

Start with p=t^6+a t^4+b t^2 and a monic odd q of degree nine. Suppose scalar constants A,B can be chosen so that r=q^2-p^3-Ap^2-Bp has degree at most four. The coefficients in degrees sixteen and fourteen force

    q = t^9+(3a/2)t^7+(3b/2+3a^2/8)t^5+u t^3+v t.

The degree-twelve coefficient gives

    A=a^3/8-3ab/2+2u.

The degree-ten coefficient forces

    v=7a^4/128-9a^2b/16+au/2+3b^2/8.

After these eliminations, the degree-eight coefficient is exactly

    D(5a^3-36ab+32u)/128.                       (E)

For D != 0, equation (E)=0 forces u=a(-5a^2+36b)/32 and yields exactly the family of Section 3. The degree-six term then determines B. Thus the family is exhaustive in this normalized parity ansatz with D != 0.

For D=0, equation (E) vanishes automatically, but this is not a hidden primitive family. With h=t^3+(a/2)t and c=u-a^3/8, the full formulas become

    p=h^2, q=h^3+c h, A=2c, B=c^2, r=0.

The field is k(h), not k(t). This supplies a sharp obstruction for the exceptional denominator branch, rather than discarding it by division.

The construction therefore consists of two controlled cancellations: truncate the formal 3/2 power at infinity, then tune the residual square-root component. Equation (E) is the inexpensive decisive test. It distinguishes a genuine degree-four cancellation from a common cubic-composition collapse.

## 5. Why the module statements give the whole algebra

For any monic p of degree n with k(p,q)=k(t), the element q is integral over k[p] because k[t] is a finite k[p]-module. Its minimal polynomial over k(p) has degree n by the field tower. Its coefficients are in the integrally closed ring k[p]. Consequently 1,q,...,q^(n-1) is a free k[p]-basis of k[p,q].

Replacing q^2 by r=q^2-H(p), and its powers by r^i or q r^i, is triangular with diagonal entries one. Thus each basis used above spans exactly k[p,q], rather than merely a smaller subalgebra with the same fraction field.

In a degree-separated free basis b_i, all nonzero summands f_i(p)b_i have degrees n deg(f_i)+deg(b_i) in different residue classes. This prevents cancellation between summands. The codimension is therefore the number of omitted degrees, not just a lower bound from displayed elements.

## 6. Reproduction and scope

Run, from this directory,

    python3 verify_realizations.py > results.json
    python3 -O verify_realizations.py > results-optimized.json
    cmp results.json results-optimized.json

The checker verifies the concrete identities; both parameter identities; derivative remainders, resultants and bad divisors; both module-transform determinants; primitive degree gcds; exact residue-separated monoids; their gap sets; and the parameter-degeneration controls. It uses explicit exceptions rather than bare assertions.

The symbolic identities prove families over characteristic zero. The finite semigroup checks cover a declared bound and are supplemental to the all-degree residue-reduction proofs above. No literature theorem is needed for realizability, and no external-source novelty claim is made. The construction says nothing by itself about allowing a polynomial x-dependence with nonzero constant planar Jacobian.
