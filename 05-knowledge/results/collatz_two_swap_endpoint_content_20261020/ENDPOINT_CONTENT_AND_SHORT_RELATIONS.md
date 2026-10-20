# Exact endpoint-chain content and a subexponential common-root bound

20 October 2026. Status: **PROVED exact arithmetic identities and content bounds, independently audited with exact controls.** These are two specified overlapping two-edit families. No assertion about all two-edit shapes follows from these identities alone.

## 1. Shapes, actual paths, and theorem

Let p>=3, q=p+r, p/2<r<p, and gcd(p,r)=1. The clock D=2^q-3^p may have either sign. Use the zero-intercept lower mechanical word and the residue labelling a_j=(-jr mod p) from the [two-edit sparse-form theorem](../collatz_two_swap_norm_barrier_20261020/symbolic/TWO_SWAP_SPARSE_FORMS_AND_BARRIER.md).

The endpoint-plus path moves the one of residue p-1-r right, then the one of residue p-1 right. These are a legal overlapping 110 -> 101 -> 011 chain. The endpoint-minus path moves the one of residue r left, then the one of residue 0 left, giving a legal 011 -> 101 -> 110 chain. Crossing the original cut is allowed and does not change the carry gcd.

For each final word w, define its actual content

    g=gcd(|D|,B(w)).

**Exact content theorem.** The respective contents are

    g_plus =gcd(2*13^p-10^p, 3*13^r-2*10^r),           (1.1)
    g_minus=gcd(9^p-2*5^p, 2*9^r-3*5^r).              (1.2)

These are equalities of integers, including every prime-power depth, not merely necessary divisibility gates.

**Short-relation theorem.** Put N=ceil(sqrt(p)). Then

    g_plus <52^N,    g_minus <36^N.                    (1.3)

In particular, for either chain and every p>=25,

    g/3^p <exp(-37p/300)<exp(-p/9).                    (1.4)

The exact norm modulus can have asymptotic base 3 on nearby shapes, but the actual content retains the common clock root and obeys (1.3). There is no conflict: a resultant discards that additional equation.

The coprime small boundary p=2,r=1,q=3 is legal but outside the strict displayed band. Directly, its endpoint-plus output is 101 and its endpoint-minus output is 110, with fixed points -7 and -5. Both formulas (1.1)-(1.2) give content 1 there. Any later integer-exclusion corollary must retain this genuine negative-cycle exception.

## 2. The modular root and the discarded clock relation

For any divisor h of |D|, h is coprime to 6. If h>1, choose integers u,v with pu+rv=1 and set

    z=2^u(3/2)^v modulo h.

All divisions and negative powers are in the unit group of Z/hZ. The clock relation 2^(p+r)=3^p gives

    z^p=2,    z^r=3/2.                                (2.1)

Moreover z is the unique unit satisfying these two equations: a common solution must equal (z^p)^u(z^r)^v. The element z-1 is also a unit, because

    (z-1)(1+z+...+z^(p-1))=z^p-1=1.

The inherited moving-anchor carry identity says that h divides the final carry if and only if its canonical sparse polynomial vanishes at this z. Every omitted factor is a power of 2 or 3, a power of z, or z-1, and hence is a unit modulo h. This is valid over the whole possibly non-squarefree modulus h.

The two integral sparse representatives are

    H_plus(Z)=3-Z^(p-1)+Z^(p-r)-Z^(p-r-1),
    H_minus(Z)=3-Z-Z^(r+1)+Z^r.

For clarity, H_minus is the integral representative 2Q_raw, not the un-cleared rational edit expression. This distinction matters in the next calculation.

Using both equations in (2.1),

    z H_plus(z)
      =3z-2+(z-1)z^(p-r)
      =3z-2+(4/3)(z-1)
      =(13z-10)/3,

and

    H_minus(z)=3-z+(1-z)(3/2)=(9-5z)/2.               (2.2)

Thus the actual-content equations are equivalent to

    plus: 13z=10,    minus: 5z=9.                     (2.3)

## 3. Necessity and sufficiency: why the gcd formulas are exact

### Plus chain

Apply (2.3) modulo the actual g. It forces 13 to be a unit: a prime 13 dividing g would make 10 zero modulo 13. It also forces 5 to be a unit, since z is a unit and 13z=10. Hence z=10/13 and the two clock equations give

    10^p=2*13^p,
    2*10^r=3*13^r modulo g.

This proves that g divides the right side of (1.1).

Conversely let G be the gcd on the right side of (1.1). If G=1 there is nothing to show. The first gcd argument is 1 modulo 3 and nonzero modulo 5 and 13; the second is odd. Consequently G is coprime to 2,3,5,13. Define z=10/13 modulo G. The two gcd equations give z^p=2 and z^r=3/2. Raising them to the r-th and p-th powers shows

    2^r=(3/2)^p,

so G divides D. By Bezout, z is exactly the unique clock root of Section 2. Equation 13z=10 now gives H_plus(z)=0, which is equivalent to G dividing the carry. Thus G divides g, proving (1.1).

### Minus chain

Equation 5z=9 modulo the actual g forces 5 to be a unit, so z=9/5. The two clock equations yield (1.2) as a necessary divisibility gate.

For the reverse implication let G be its right side. The first argument 9^p-2*5^p is odd and is nonzero modulo 3 and 5. Thus G is coprime to 2,3,5. Set z=9/5 modulo G. Both clock equations hold, they imply G divides D, and Bezout identifies z with the clock root. Equation 5z=9 makes H_minus(z)=0. Hence G divides the carry, proving (1.2).

No radicalization, squarefree reduction, or bound on multiplicative order occurs in either argument.

## 4. A short simultaneous exponent relation

Let N=ceil(sqrt(p)), so N<p for p>=3. The N+1 residues

    R_j=j r mod p,    j=0,...,N,

are distinct, since gcd(p,r)=1. Put them into the N half-open bins indexed by floor(N R_j/p). Two, say R_i and R_j with i>j, share a bin. Define

    b=i-j,
    a=-(floor(ir/p)-floor(jr/p)),
    m=ap+br=R_i-R_j.

Then

    1<=b<=N,    -N<=a<=0,    0<|m|<p/N<=N.           (4.1)

The bound on a follows from 0<r<p and i,j in {0,...,N}. Distinctness of the residues makes m nonzero.

For either actual clock root, (2.1) implies

    z^m=(z^p)^a(z^r)^b=2^(a-b)3^b.                  (4.2)

Now substitute the rational value of z. This is the information that the norm construction did not retain.

### Plus bound

For z=10/13, clearing units in (4.2) gives an integer divisible by g:

    E=2^(b-a)10^m-3^b13^m,             if m>0,
    E=2^(b-a)13^(-m)-3^b10^(-m),       if m<0.

Neither integer is zero: in each case precisely one positive term is divisible by 13. Because b-a<=2N and |m|<=N, their two terms are bounded respectively by

    m>0: 40^N and 39^N,
    m<0: 52^N and 30^N.

The absolute difference of two positive terms is strictly less than their larger upper bound. Therefore 0<|E|<52^N, and g<=|E| proves the plus bound in (1.3).

### Minus bound

For z=9/5 the corresponding integers are

    E=2^(b-a)9^m-3^b5^m,               if m>0,
    E=2^(b-a)5^(-m)-3^b9^(-m),         if m<0.

Exactly one term is divisible by 5, so E is nonzero. The term bounds are

    m>0: 36^N and 15^N,
    m<0: 20^N and 27^N.

Hence g<=|E|<36^N, proving the second bound.

The proof is elementary simultaneous approximation on a finite set of residues. It is a bound on actual gcd content, not a hypothesis that those residues are statistically distributed.

## 5. A convenient exponential corollary

For p>=25,

    N<=sqrt(p)+1<=p/5+p/25=6p/25.

The exact positive Taylor series gives exp(4)>52, hence log52<4. Also

    log3=2 artanh(1/2)
         >2(1/2+(1/2)^3/3)=13/12.

Using the weaker plus bound for both families,

    log(g/3^p)<N log52-p log3
                <(24/25-13/12)p
                =-37p/300<-p/9.

This proves (1.4). If one of these words were integral, then |D|=g and the multiplicative logarithmic form would satisfy

    |2^q 3^(-p)-1|<exp(-p/9).

For p>=25 this is less than 1/2, and therefore

    |q log2-p log3|<2 exp(-p/9).

This is an input for a separately certified Diophantine closure, not a claim here that these inequalities alone exclude every p.

## 6. Exact checks, prime-power controls, and limitations

`verify_endpoint_content.py` directly generates both actual paths for every coprime shape with 3<=p<=160 and p/2<r<p. It compares the positional edit result's carry gcd with the exact formulas, checks the modular root reconstruction, constructs the short relation, and checks the nonzero integer bound. It uses the affine carry recurrence, not the gcd formula, to evaluate each edited word.

The main universe has 3,902 shapes and 7,804 paths, including 158 original-cut crossings and 663 negative-clock shapes. The completed run has 58,985 explicit checks, with byte-identical normal and optimized receipts. No floating-point arithmetic controls a decision.

Two additional exact controls retain nontrivial prime powers:

| family | p | r | q | actual g |
|---|---:|---:|---:|---:|
| plus | 418 | 237 | 655 | 11849=17^2*41 |
| minus | 203 | 165 | 368 | 169=13^2 |

For the legal plus norm barrier at p=28,r=17, actual content is 1 despite the norm exceeding 3^28. At the first tuned near-critical member p=18,r=13, actual content is 17. The exact formulas preserve this distinction.

Reproduce from this directory:

    python3 verify_endpoint_content.py > receipt.replay.json
    python3 -O verify_endpoint_content.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json

The separate independent audit checks direct binary-word carries, rotation transport, both prime-power controls, the reverse gcd implication, and all four short-integer cases. Its files are in `audit/`.

The mechanism is specific: two endpoint polynomials become linear once the clock root is retained. Most other two-edit displacement patterns do not become either equation in (2.3). Nothing here supplies a uniform actual-content bound for those remaining primitive shapes.
