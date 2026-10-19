# All two-cycles in the cofactor deformation, with a sharp core bound

19 October 2026. PROVED, with independent audit PASS. Let F_t(cp)=(c-1)p+2tc-1 on positive even integers, where c is the exact dyadic factor and p is odd, and t is a positive integer. This is an auxiliary arithmetic family; only t=1 on the odd-prime domain equals the aliquot map.

## 1. Exact classification

For any two distinct dyadic cores c=2^a<d=2^b, set

    D=c+d-1,
    g=2^gcd(a+1,b+1)-1.

There is a nontrivial two-cycle with these exact cores if and only if

    D/g divides 2t-1.

When it exists it is unique for those cores. Put

    j=(2t-1)(2d-1)/D,
    j'=4t-2-j.

The cycle consists of

    n=c[2t(d-1)+j],
    m=d[2t(c-1)+j'].

Both bracketed quantities are positive odd integers.

**Proof.** Write n=cp and m=dq. The two exact equations

    (c-1)p+2tc-1=dq,
    (d-1)q+2td-1=cp

have determinant c+d-1>0. Solving gives

    p=2t(d-1)+(2t-1)(2d-1)/D,
    q=2t(c-1)+(2t-1)(2c-1)/D.

The two final fractions sum to 4t-2. Thus either fraction is integral exactly when both are; being quotients of odd integers, their integral values are odd. Positivity is immediate. The resulting points therefore have exactly the specified cores and genuinely satisfy the map equations.

Now

    gcd(D,2d-1)=gcd(2c-1,2d-1)
               =2^gcd(a+1,b+1)-1=g,

using the elementary gcd identity for Mersenne numbers. Hence integrality is exactly D/g | 2t-1.

If both points had the same core c, they would follow the same affine map of slope (c-1)/c strictly between zero and one. Its second iterate has the same unique fixed point, so there is no nontrivial same-core two-cycle. Distinct cores also ensure distinct points. This proves completeness. ∎

## 2. Finite search bounds for every fixed parameter

Every nontrivial two-cycle satisfies

    c<=2t-2,
    d<=8(t-1)^2.

Thus the preceding criterion is a finite exact classification for each supplied t, even though orbit values and dyadic exponents initially had no bound. In particular t=1 has no nontrivial two-cycles.

To prove the bounds, c<d gives

    2t-1 < j <4t-2.

Because d is divisible by c, reducing jD=(2t-1)(2d-1) modulo c gives

    j ≡ 2t-1 (mod c).

Hence the strictly positive integer j-(2t-1) is a multiple of c and is at most 2t-2 (j is odd). This proves the bound on c.

Rearranging the same equation gives

    d=[j(c-1)+(2t-1)]/[(4t-2)-j].

The denominator is a positive integer. Using j<=4t-3 and c<=2t-2 gives

    d<=(4t-3)(2t-3)+(2t-1)=8(t-1)^2.

No numerical cycle census is used in these all-parameter bounds.

## 3. The quadratic bound is attained infinitely often

For every dyadic c>=2 choose

    t=c/2+1,       d=2c².

Then j=2c+1=4t-3 and j'=1. The classification gives the explicit two-cycle

    c(2c³+4c²+c-1)  <->  2c²(c²+c-1).

Its larger core is exactly

    d=2c²=8(t-1)^2.

Thus the universal bound on the larger core is sharp for infinitely many parameters. For c=2, t=2 this is 66 <-> 40. For c=4, t=3 this is 780 <-> 608. For c=8, t=5 it is 10296 <-> 9088.

This infinite family varies the map parameter. It does not give infinitely many two-cycles for one fixed map; the proved finite bounds rule that out. It also explains why the two-cycle exclusion at t=1 is special rather than a consequence of having contracting branch slopes.

## 4. Why this is a useful transfer

The forward finite-odd-part reduction controls growth episodes by a valuation guard. The two-cycle theorem instead solves the actual affine cycle equations and retains exact dyadic cores. A Mersenne gcd then supplies a complete divisibility criterion, and a residue trap supplies a finite parameter-dependent bound. This is a concrete connection among local affine dynamics, inverse divisibility, and exponent lattices. It does not assume that formal core assignments are realized: the positive odd solutions certify that realization explicitly.
