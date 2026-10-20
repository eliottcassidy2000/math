# Fixed source separation: full content and a cube-root spacing obstruction

20 October 2026. Status: **PROVED, with cited Minkowski and Matveev inputs, independent exact controls, and a second full-proof audit.** The full-content theorem, quantitative spacing consequence, and fixed-separation finite-weight reduction are established in the scopes below.

## 1. Precise scope and conclusions

Let p>=6, 1<=r<p, gcd(p,r)=1, and q=p+r. The source is the lower zero-intercept mechanical word

    c_i=floor((i+1)p/q)-floor(ip/q), 0<=i<q.

Label its ones in increasing lifted order by

    s_j=j+ceil(jr/p)-1,  a_j=-jr mod p,
    s_(j+p)=s_j+q,      0<=a_j<p.

Take two distinct source labels at minimum cyclic label separation h, where

    3<=h<=floor(p/2).

Choose an orientation so the labels are j and j-h, and move each exactly once by an actual cyclic adjacent unequal-bit swap, in directions d_1,d_2 in {+1,-1}. The first direction refers to j. These two moves have disjoint affected bit sets and commute whenever they are legal. Write w for the output and

    B(w)=sum_(i:w_i=1) 2^i 3^(number of ones after i),
    D=2^q-3^p,  g=gcd(|D|,B(w)),  C_h=34*4^h.

The branch fixed point B(w)/D is integral exactly when g=|D|. It has the sign of D. The results below retain both signs and every prime-power divisor of g.

**Full-content theorem.** Every legal path just specified satisfies

    g <= C_h^(4 p^(2/3)).                              (1.1)

**Quantitative spacing consequence.** If p>=64^3, a legal output with

    3<=h<=p^(1/3)/8                                   (1.2)

cannot be integral. Equivalently, in this primitive two-distinct-source universe, any integral output with h>=3 at such a weight must have h>p^(1/3)/8.

**Effective fixed-h reduction.** Put

    N_h=ceil(720(18+7h)/881),
    P_h=max(100,N_h^3).

Every possible integral output in this family has p<P_h. For h>=3, N_h>=32, so P_h=N_h^3. This is an explicit O(h^3) bound on the remaining weights, not a proof that all those finitely many residuals are empty.

These are not universal two-edit exclusions. They do not include h=1 or h=2, repeated source words, overlapping moves, or q>=2p. Those require separate results. They do not prove the unrestricted rank-four noncancellation claim discussed in Section 10.

## 2. Inheritance and loss ledger

The closest proved mechanism is the external-denominator-prime rank-three argument of the frozen overlap-content packet, Sections 3--6. The source-label transport below supplies the additional monomial identity needed to apply that mechanism to a new family. The distance-two family is separate; this proof does not depend on its smooth-gate analysis or classification.

The source-to-target map is:

    physical labelled edits -> the actual common clock root -> a rational one-variable gate
    -> two short lattice relations -> a nonzero integer resultant.

The first two steps preserve exact carry divisibility, including prime powers. The resultant step preserves necessary divisibility but loses the converse: g divides the resultant, whereas gcd(|D|,resultant) need not equal g. The complete clock root and the original carry remain the sidecar needed to recover actual content. Section 9 includes an explicit hostile to the false converse.

The concept board is source-label transport, exact physical legality, denominator valuations, successive-minimum products, and coefficient dependence on the separation. The decisive new point is that h>=3 removes every smooth denominator by an elementary consecutive-powers argument. The coefficient dependence is retained rather than hidden in an unspecified constant.

## 3. Exact legality and the four rational gates

Put u=a_j and v=a_(j-h). Then

    v=u+hr-tp,  t=floor((u+hr)/p),  0<=t<=h.          (3.1)

An individual right move at residue a is legal iff a<r. An individual left move is legal iff a>=p-r. Thus the complete test here is

    d_1=+1: u<r;       d_1=-1: u>=p-r;
    d_2=+1: v<r;       d_2=-1: v>=p-r.               (3.2)

Because both cyclic label distances are at least three, neither moving source one is a neighbor of the other. A blocked move cannot be unlocked by moving the other chosen one: a source one can block only an immediately neighboring label. If two initially legal swaps shared a target zero, the chosen source ones would also be neighboring labels. A source position shared by the two edges would likewise require neighboring labels. Consequently (3.2) is necessary and sufficient, the affected sets are disjoint, and both orders are actual with the same output. Original-cut crossings cause no exception.

At h=p/2, reversing the orientation selects the same unordered pair. The algebra allows both presentations, with directions reassigned to their labels; it asserts the same content bound twice. A physical enumeration must deduplicate them. No extra physical path is inferred from this orientation choice.

For M=|D|>1, choose integers A,B with pA+rB=1 and define

    z=2^A (3/2)^B mod M.

The numbers 2 and 3 are units modulo M. Since (3/2)^p=2^r modulo M, this is the unique simultaneous root

    z^p=2,  z^r=3/2.                                 (3.3)

Uniqueness follows from z=z^(pA+rB). Also z-1 is a unit because

    (z-1)(1+z+...+z^(p-1))=z^p-1=1.

All of this is over the full modulus, not its radical. If M=1 the content bound is immediate; in the present range this convention is harmless.

Let alpha(+1)=1 and alpha(-1)=-1/2. The frozen moving-anchor identity, Section 3 of `TWO_SWAP_SPARSE_FORMS_AND_BARRIER.md`, gives the exact common-root carry gate

    Q(z)=1+(z-1)[alpha(d_1)z^u+alpha(d_2)z^v].       (3.4)

Carry divisibility by any divisor of M is equivalent to Q(z)=0 modulo that divisor. Every factor transported by the moving cut is a unit. In particular, this is not an assumption that the source anchor stays fixed.

Transport (3.1) using both clock equations:

    z^v/z^u=(z^r)^h/(z^p)^t=3^h/2^(h+t).

Put T=2^(h+t) and W=3^h. The four gates z^u(z-1)=K are exactly

    ++: K=-T/(T+W),
    --: K= 2T/(T+W),
    +-: K=-2T/(2T-W),
    -+: K=-T/(W-T/2).                                (3.5)

The two mixed denominators can have either sign. Normalize each fraction as

    K=epsilon*2^delta/J,  epsilon in {+1,-1}, J>0.

Explicitly:

    ++: J=T+W,      delta=h+t,   epsilon=-1;
    --: J=T+W,      delta=h+t+1, epsilon=+1;
    +-: J=|2T-W|,   delta=h+t+1, epsilon=-sign(2T-W);
    -+: J=|W-T/2|,  delta=h+t,   epsilon=-sign(W-T/2).

Every denominator is odd and coprime to 3, and no further cancellation with its dyadic numerator occurs. Therefore

    h<=delta<=2h+1,    J<=2*4^h+3^h.                 (3.6)

Writing X=z^u(z-1), the original coefficient in (3.4) is -epsilon J/2^delta. Hence

    2^delta Q(z)=-epsilon[JX-epsilon*2^delta].

The multiplying factor is a unit modulo M; importantly, this argument does not divide by J modulo M. Thus the exact identity is

    g=gcd(M,J z^u(z-1)-epsilon*2^delta).              (3.7)

Every prime dividing J is immediately excluded from g by the right side, so gcd(g,6J)=1. Only now do we use the divided equation modulo g. This order of operations preserves complete prime-power depth.

## 4. No smooth denominator when h>=3

The same-sign denominators T+3^h are greater than one. A mixed denominator equal to one would give

    |3^h-2^n|=1,

with n=h+t+1 in the +- case and n=h+t-1 in the -+ case. In particular n>=2 when h>=3.

Here is a complete elementary check. If 2^n-3^h=1 and n>=3, reduction modulo 8 is impossible: 3^h is 1 or 3 modulo 8, so -3^h is 7 or 5, never 1. For n=1 or 2, direct substitution leaves only (n,h)=(2,1).

If 3^h-2^n=1 and n>=3, reduction modulo 8 makes h even. Then

    (3^(h/2)-1)(3^(h/2)+1)=2^n.

The factors are positive even integers differing by two and with gcd two. Since their product is a power of two, the smaller factor must be 2 and the larger 4. Therefore h=2,n=3. For n=1 the only solution is h=1, while n=2 would require 3^h=5, impossible.

Thus every solution has h<=2. In particular

    J>1 and gcd(J,6)=1 for all h>=3.                 (4.1)

Choose any prime ell dividing J; it is outside {2,3}. No Catalan or Mihailescu theorem is needed. The restriction is genuine: at h=2,t=0,+- the gate is K=8, and at h=2,t=2,-+ it is K=-16. Those smooth gates need a separate mechanism.

## 5. Rank-three relations, full-modulus polynomials, and degeneracies

Use the determinant-p lattice

    L={(b,c,m) in Z^3 : m-br-cu is divisible by p}.

Its determinant is p because the congruence map onto Z/pZ has coefficient one on m. For a nonzero vector define

    a=(m-br-cu)/p,  H=max(|b|,|c|,|m|)>=1.

Since 0<=u<p and 1<=r<p,

    |a|<=H(1+r+u)/p<2H.                              (5.1)

The relaxed bound |a|<=2H is sufficient. Negate the vector if necessary so c>=0. From the two clock equations and (3.5), modulo the whole g,

    z^m(z-1)^c=epsilon^c 2^(a-b+delta c)3^b J^(-c).

Put e=a-b+delta c and

    U=2^max(-e,0)3^max(-b,0)J^c,
    V=epsilon^c 2^max(e,0)3^max(b,0),
    F_v(X)=U X^max(m,0)(X-1)^c - V X^max(-m,0).      (5.2)

All denominators cleared are units modulo g, so F_v(z)=0 modulo g, not merely modulo each prime dividing it.

For completeness, every nonzero vector gives a nonzero polynomial with neither 0 nor 1 as a root. If m!=0, evaluation at 0 leaves one nonzero term. If m=0,c>0, cancellation there would force the rational right side above to be (-1)^c, contradicting its negative ell-adic valuation. If c=0, such cancellation forces a=b=0 by unique factorization, hence m=0 and the zero vector. Evaluation at 1 leaves -V when c>0, and the c=0 case is the same unique-factorization argument.

The only possible nontrivial leading cancellation has m<0, -m=c>0, and U=V. That is again impossible at ell. Thus F_v is constant exactly when m=c=0, and

    deg F_v<=|m|+c<=2H.                              (5.3)

Two independent vectors cannot both give constants: m=c=0 confines them to a single line.

## 6. Exact coefficient bound and characteristic-zero noncancellation

The binomial theorem gives

    ||F_v||_1 <= 2^(max(-e,0)+c)3^max(-b,0)J^c
                 +2^max(e,0)3^max(b,0).              (6.1)

As a function of (a,b,c), each summand is an exponential of a convex piecewise-linear function, so their sum is convex. Its maximum on the relaxed box

    |a|<=2H, |b|<=H, 0<=c<=H

is at a vertex. At H=1 and c=0 the four values are 7,11,11,7. At c=1, delta>=3 implies e>=0 at all vertices. For b=-1 the expression is 6J+2^(delta+a+1), maximized at a=2. For b=+1 it is 2J+3*2^(delta+a-1), which is at most the preceding maximum. Therefore the exact maximum on this box is

    C=6J+2^(delta+3).                                (6.2)

For integer H>=1 a vertex value is A^H+B^H<=(A+B)^H. Hence ||F_v||_1<=C^H. Using (3.6),

    C <=6(2*4^h+3^h)+2^(2h+4)
      =28*4^h+6*3^h
      <=34*4^h=C_h.                                  (6.3)

This proves the explicit coefficient bound. It is uniform in r,u,t and the two directions, but not independent of h.

Now take any two independent lattice vectors. Their polynomials are coprime over Q. Indeed, a nonzero constant cannot share a root with anything. Otherwise a common algebraic root x would satisfy x!=0,1 and both multiplicative identities used in (5.2).

If the projected rows (m_i,c_i) are dependent, an integer linear combination kills m and c. Independence of the full rows leaves a nonzero pair (a,b) with ap+br=0, so (a,b)=(rt,-pt) for some nonzero integer t. The same combination of the two root identities forces

    1=2^((p+r)t)3^(-pt),

contradicting unique factorization.

If the projected rows are independent, choose a valuation nu on Q(x) above ell, normalized by nu(ell)=1, and write e_ell=nu(J)>0. The two equations give

    m_i nu(x)+c_i[nu(x-1)+e_ell]=0, i=1,2.

Independence forces nu(x)=0 and nu(x-1)=-e_ell<0. This contradicts the ultrametric inequality nu(x-1)>=min(nu(x),nu(1))=0. The contradiction uses the full valuation e_ell; it does not replace J by its radical. This establishes noncancellation for every admissible gate.

## 7. Short independent relations imply the full-content theorem

Let lambda_1<=lambda_2<=lambda_3 be the successive minima of L for the cube [-1,1]^3. The cited Minkowski second theorem gives

    lambda_1 lambda_2 lambda_3<=det L=p.

The cube's volume eight cancels the theorem's factor 2^3. Its integer lattice minima are attained. Choose independent vectors of heights H_1=lambda_1,H_2=lambda_2. Since lambda_3>=sqrt(lambda_1 lambda_2),

    H_1 H_2<=p^(2/3).                                (7.1)

This is a product bound; it is not an assertion that both heights are individually O(p^(1/3)). The exact source is Barvinok, *Math 669: Combinatorics, Geometry and Complexity of Integer Points*, Theorem 22.5, printed page 46. The primary PDF was checked for this audit: https://sites.lsa.umich.edu/barvinok/wp-content/uploads/sites/1434/2025/05/latticenotes669.pdf .

By Section 6, R=Res(F_1,F_2) is a nonzero integer. The integral resultant Bezout identity evaluated at their common modular root gives g|R, including each entire prime power dividing g. If one polynomial is constant, the resultant is a positive power of that constant and the same divisibility holds directly. The two-constant case was excluded in Section 5.

Sylvester--Hadamard, (5.3), (6.3), and (7.1) yield

    g<=|R|
      <=||F_1||_1^deg(F_2) ||F_2||_1^deg(F_1)
      <=C_h^(4H_1H_2)
      <=C_h^(4p^(2/3)).

This proves (1.1) for the actual gcd, not merely for its radical and not for the raw sparse norm.

## 8. Exact arithmetic-gap transfer and the two quantitative corollaries

The inherited absolute-gap lemma, Sections 3--4 of `ENDPOINT_TWO_EDIT_CLASSIFICATION.md`, says:

    If gcd(p,q)=1, p>=100 and p<q<2p, then
    0<|2^q/3^p-1|<exp(-p/9) is impossible.            (8.1)

Its proof depends on the displayed arithmetic hypotheses, not on endpoint edits. It applies the real Matveev estimate to U=2^q/3^p-1, not directly to L=q log2-p log3. It first gives p<10^12, then uses |L|<2|U| for both signs, Legendre's criterion, and exact rational intervals for all eighteen relevant convergents. The primary cited estimate is Bugeaud--Mignotte--Siksek, *Classical and modular approaches to exponential Diophantine equations I: Fibonacci and Lucas perfect powers*, Theorem 9.4, real case, printed page 16: https://arxiv.org/pdf/math/0403046 . The theorem's real-case constant and the multiplicative-minus-one form were checked directly. This packet reruns the frozen independent large-range certificate in normal and optimized Python, with identical receipts. No signed-descent or global Collatz-convergence assumption is used here.

The rational positive-series log bounds certified in `certify_spacing_constants.py` are

    log34<18/5,  log4<7/5,  log3>109/100.             (8.2)

From (1.1),

    log(g/3^p)/p <=4(log34+h log4)/p^(1/3)-log3.

If p>=64^3 and h<=p^(1/3)/8, (8.2) gives strictly

    log(g/3^p)/p
      <4*(18/5)/64+4*(1/8)*(7/5)-109/100
      =-33/200<-1/9.                                (8.3)

An integral output would have g=|D| and would therefore violate (8.1). Its domain p>=100 is automatic because p>=64^3. This proves the spacing consequence.

For fixed h, (8.2) instead gives

    log C_h <(18+7h)/5,
    log3-1/9>881/900.

Thus p^(1/3)>=720(18+7h)/881 is sufficient for the strict bound g/3^p<exp(-p/9). Combining with the explicit p>=100 requirement in (8.1) proves the cutoff P_h=max(100,N_h^3) stated in Section 1. For h>=3, N_h>=32, so the maximum is redundant but has not been silently dropped from the transfer. Examples are N_3=32 and cutoff 32^3, N_4=38 and cutoff 38^3, and N_8=61 and cutoff 61^3.

For a fixed h there are only finitely many p<P_h, primitive r, source pairs, and direction choices. Direct physical carry computation therefore gives an effective finite decision problem. This packet does not execute or assert emptiness of that complete residual. A finite height reduction is not an all-height exclusion.

## 9. Independent physical, prime-power, and resultant controls

`verify_fixed_separation.py` imports no repository mathematics. It constructs source bits directly, discovers their one positions by scanning the actual word, tracks the chosen source identities, and attempts each swap against the current bits in both orders. It then compares the result with the residue-legality test. In particular, illegal presentations are tested, not silently discarded by the proposed formula.

Its exhaustive universe is every primitive shape 6<=p<=40, every 3<=h<=floor(p/2), and all four direction choices. It also tests every primitive r at p=43,61,89 for h in {3,4,5,7,11,floor(p/2)}. The algebraic gate table and relaxed-box constants are checked for every 3<=h<=64 and 0<=t<=h. Actual carries use two independent formulas. The exact clock, raw carry gate, normalized gate, and full gcd all agree.

Through p=18, it enumerates exact first and second sup-norm lattice minima and checks cleared polynomials and their resultants. Unlike the inherited Sylvester/Bareiss implementation, these resultants use a rational Euclidean recurrence. Additional controls cover a prime-cube content, equal projected rows, and a nonzero constant polynomial.

The completed physical run records:

- 670 primitive shapes and 6,562 shape/separation pairs;
- 328,268 legal oriented presentations and 658,996 illegal presentations;
- 650,572 distinct ordered physical paths, including 25,908 cut crossings;
- 262,287 legal positive-clock and 65,981 legal negative-clock presentations;
- 5,964 h=p/2 oriented presentations and exactly 5,964 physical ordered paths, correctly removing their duplicate orientation;
- 7,474 exhaustive short-pair resultant checks, plus three explicit special controls;
- no integral output in this finite universe.

The largest tested content is g=56,461 at p=89,r=79,h=7,u=17,v=36 with directions --. A prime-cube control is

    p=19,r=14,h=3,u=2,v=6,t=2,directions +-,
    K=-64/37,  g=125.

At this last shape, vectors (0,1,2) and (19,1,2) are independent but have equal (m,c) projections; they still give a nonzero resultant divisible by 125. Pairing (19,0,0) with (0,1,2) instead supplies a constant polynomial D. The resulting gcd(|D|,R)=|D| is strictly greater than the actual g=125. This demonstrates why resultant divisibility must not be promoted to an iff content gate.

The h=2 smooth rows K=8 and K=-16 are retained as explicit algebraic hostiles to removing h>=3. The p=6,r=1,h=3 invalid moves are retained as physical legality controls. At p=6,r=5,h=3, sixteen oriented presentations give sixteen ordered paths, not thirty-two.

Reproduce from this directory:

    python3 verify_fixed_separation.py > receipt.replay.json
    python3 -O verify_fixed_separation.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json
    python3 certify_spacing_constants.py > spacing_constants.replay.json
    cmp spacing_constants.json spacing_constants.replay.json

`REPRODUCE.md` gives the frozen inherited absolute-gap replay command. All assertions are explicit runtime checks and remain active under Python optimization. `DEPENDENCIES.json` pins the actual source hashes and distinguishes proved inputs from the candidate read for provenance.

## 10. Exact remaining boundary

For arbitrary disjoint edits, the faithful equation is

    1+(X-1)(alpha Y+beta W)=0,
    X=z, Y=z^u, W=z^v.

Without a bounded source-label transport, the congruence lattice in four coordinates gives bivariate equations after one eliminates W. Independent rank-four lattice rows do not automatically become independent univariate equations, and the present external-denominator valuation proof does not establish their universal coprimality. Clearing denominators must still retain X,Y,W,X-1!=0; spurious boundary components are not actual clock roots.

Every pair does have a minimum label separation h<=p/2, so (1.1) is a valid numerical bound throughout the stated h>=3 family. Its exponent, however, is O((h+1)p^(2/3)). It is sublinear in p only when h=o(p^(1/3)), and the explicit spacing exclusion uses the narrower sufficient range (1.2). At h comparable to p the displayed estimate is far too large to prove nonintegrality.

The proven advancement is therefore a separation-sensitive full-content mechanism, a cube-root necessary spacing condition at large weights, and an effective O(h^3) fixed-h reduction. Unrestricted rank-four noncancellation, complete finite-residual emptiness, and arbitrary two-edit exclusion remain open.
