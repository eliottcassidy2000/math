# Subexponential actual content for every primitive mechanical two-one overlap

20 October 2026. Status: **PROVED, with the cited Minkowski theorem and an independent exact audit.** The core valuation, lattice, and coefficient estimates have passed independent mathematical review and exact checks. A separate finite-weight consequence is stated at the end, not an all-length exclusion.

## 1. Scope and theorem

Let

    p>=2,  1<=r<p,  gcd(p,r)=1,  q=p+r,
    c_i=floor((i+1)p/q)-floor(ip/q).

Label source ones by positions and residues

    s_j=j+ceil(jr/p)-1 modulo q,
    a_j=-jr modulo p.

We consider two actual adjacent swaps whose affected edges overlap, and which move two distinct labelled ones in the same direction:

* **Plus:** the local path 110 -> 101 -> 011. For a parameter u, move residue u right first, then residue u+r right.
* **Minus:** the local path 011 -> 101 -> 110. Move residue u+r left first, then residue u left.

The complete legal parameter ranges, including both clock signs, are

    plus:  0<=u<min(r,p-r),
    minus: max(0,p-2r)<=u<p-r.                         (1.1)

They are not interchangeable when r<p/2. Swaps across the original cut are included.

For the final word w, define

    B(w)=sum_{i:w_i=1} 2^i 3^(number of ones after i),
    D=2^q-3^p,    g=gcd(|D|,B(w)).

**Theorem.** Every legal path in (1.1) satisfies

    g_plus <=38^(4 p^(2/3)),
    g_minus<=46^(4 p^(2/3))<exp(16 p^(2/3)).           (1.2)

These bounds concern the full actual gcd, including every prime-power depth. They are uniform in r, u, and the sign of D. In particular, actual content is subexponential in p throughout these two families.

The result is not a comparable bound for the original sparse norm. It retains the simultaneous clock equations discarded by that norm. It does not cover arbitrary disjoint edits, repeated mechanical sources, or a single labelled one moving twice through 100 -> 010 -> 001. The latter operation cannot occur in this k=1 source universe because its zeros are isolated, but it exists for other mechanical shapes and must not be conflated with a two-one overlap.

## 2. Exact legality and inherited common-root equations

Consecutive source ones have residues u+r followed by u precisely when their source gap is one, with 0<=u<p-r. A first right move is legal exactly for residue less than r. Therefore moving the leading one u right unlocks the following one u+r exactly when u<r, yielding the plus range.

A first left move is legal exactly for residue at least p-r. Moving the preceding one u+r left unlocks u exactly when u+r>=p-r, yielding the minus range. Reversing either order is blocked.

For q>=4, distinct adjacent edges sharing a bit have a three-bit union. Checking the possible actual two-swap paths gives 100 <-> 001, 110 <-> 011, or undo. The first type is absent because mechanical k=1 sources have no adjacent zeros. Thus (1.1) captures every nonundo overlapping path in this stated source universe. The q=3 boundary is checked directly and has the two declared two-one paths.

For any divisor h>1 of |D|, 2 and 3 are units modulo h. Choose integers v,t with pv+rt=1 and define

    z=2^v(3/2)^t modulo h.

The clock gives the unique simultaneous root

    z^p=2,  z^r=3/2.                                 (2.1)

Also z-1 is a unit, since (z-1)(1+z+...+z^(p-1))=1. The [moving-anchor carry identity](../collatz_two_swap_norm_barrier_20261020/symbolic/TWO_SWAP_SPARSE_FORMS_AND_BARRIER.md) is valid over this whole possibly non-squarefree modulus. It gives

    plus:  Q_raw=1+(z-1)(z^u+z^(u+r)),
    minus: Q_raw=1-(z-1)(z^u+z^(u+r))/2.

Carry vanishing is equivalent to Q_raw=0; all transported factors are units. With (2.1), the two equations become

    z^u(z-1)=K,
    K=-2/5 for plus,  K=4/5 for minus.                (2.2)

Before dividing by 5, these equations read 2+5z^u(z-1)=0 and 4-5z^u(z-1)=0. Thus 5 cannot divide the actual g. Consequently g is coprime to 30 and (2.2) is legitimate modulo its full prime powers.

For g=1 the desired bound is immediate; hence all following modular reasoning can assume g>1.

## 3. The exponent lattice and cleared polynomials

Write

    K=epsilon*2^delta/5,
    (epsilon,delta)=(-1,1) for plus, (+1,2) for minus.

Consider the rank-three lattice

    L={ (b,c,m) in Z^3 : m-br-cu is divisible by p }.  (3.1)

The use of c as a lattice coordinate here is unrelated to the source bits c_i. The lattice has determinant p: the congruence map to Z/pZ is onto because the coefficient of m is one. Equivalently, columns (1,0,r), (0,1,u), (0,0,p) form a basis.

For a lattice vector define

    a=(m-br-cu)/p,
    H=max(|b|,|c|,|m|).

Every legal u has r+u<=p-1. Hence for a nonzero vector, H>=1 and

    |a|<=H(1+r+u)/p<=H.                              (3.2)

The three modular equations imply

    z^m(z-1)^c=2^(a-b)3^b K^c.                       (3.3)

We may negate the whole vector to arrange c>=0. Put

    e=a-b+delta*c,
    U=2^max(-e,0) 3^max(-b,0) 5^c,
    V=epsilon^c 2^max(e,0) 3^max(b,0).

Then the integral polynomial

    F_v(X)=U X^max(m,0)(X-1)^c - V X^max(-m,0)         (3.4)

vanishes at the actual z modulo the full g. Denominator clearing used only powers of 2, 3 and 5, all units modulo g.

### Degeneracy bookkeeping

For every nonzero lattice vector, F_v is not the zero polynomial. Neither 0 nor 1 is a root:

* A possible evaluation at 0 can cancel only when m=0 and the rational right side of (3.3) equals (-1)^c. If c>0, its 5-adic valuation is -c, impossible. If c=0, cancellation would require a=b=0, hence the zero vector.
* At 1, c>0 leaves the nonzero second term. For c=0, cancellation would require 2^(a-b)3^b=1, again forcing a=b=m=0.

The degree is zero exactly when m=c=0. If m<0 and c>0, a possible leading cancellation requires -m=c and U=V. The latter would give 5-adic valuation -c=0 for the right side of (3.3), impossible. All other nonconstant cases have unequal term degrees.

Two independent lattice vectors cannot both give constant polynomials: m=c=0 confines them to the same b-axis. This matters because the formal resultant of two constants is not the required integral Bezout certificate.

## 4. Universal noncancellation from the prime 5

**Coprimality lemma.** The two polynomials obtained from any two linearly independent vectors of L are coprime over Q.

If one is a nonzero constant the assertion is immediate. Otherwise suppose they share an algebraic root x. Section 3 gives x!=0,1, so both rational equations (3.3) hold at x.

Take a nonarchimedean valuation nu on Q(x) above 5, normalized by nu(5)=1. Let

    X=nu(x),  Y=nu(x-1).

Because the other rational primes in (3.3) are 2 and 3, the two equations give

    m_i X+c_i(Y+1)=0,      i=1,2.                    (4.1)

If the two projected pairs (m_i,c_i) are independent, (4.1) forces X=0 and Y=-1. This contradicts

    nu(x-1)>=min(nu(x),nu(1))=0.

If the projected pairs are dependent, some nonzero integer linear combination kills both m and c. Since the full lattice vectors are independent, the resulting (a,b) is not (0,0). The relation ap+br=0 and gcd(p,r)=1 give

    a=rt, b=-pt,    t a nonzero integer.

Taking the same combination of (3.3) gives

    1=2^((p+r)t)3^(-pt),

which is impossible by unique factorization. Thus no common algebraic root exists in either case.

The prime 5 is used as a characteristic-zero valuation witness to noncancellation. We are not asserting that 5 divides the actual content; it was explicitly excluded from that content in Section 2.

## 5. Two short independent vectors, not merely one

Let lambda_1<=lambda_2<=lambda_3 be the successive minima of L for the closed cube [-1,1]^3. Since L is an integer lattice, these minima are positive integers, attained by independent lattice vectors.

Minkowski's second theorem gives

    lambda_1 lambda_2 lambda_3 <= det L=p.             (5.1)

Here the cube has volume 8, exactly cancelling the theorem's dimension-three factor 2^3. The precise theorem, including its upper-bound direction, is available in [Barvinok's lattice notes, Theorem 22.5, printed page 46](https://sites.lsa.umich.edu/barvinok/wp-content/uploads/sites/1434/2025/05/latticenotes669.pdf).

Choose two independent vectors with H_1=lambda_1 and H_2=lambda_2. Since lambda_3>=sqrt(lambda_1 lambda_2), (5.1) implies

    H_1 H_2<=p^(2/3).                                (5.2)

One cannot replace this argument with an assertion that both vectors individually have size O(p^(1/3)): lattices can have unequal successive minima. The product bound (5.2) is precisely what the resultant estimate needs.

## 6. Explicit coefficient and resultant bounds

For F_v in (3.4), write ||F_v||_1 for the sum of absolute coefficients. Its degree satisfies

    d_v<=|m|+c<=2H.                                  (6.1)

By the binomial theorem,

    ||F_v||_1
      <= 2^(max(-e,0)+c) 3^max(-b,0) 5^c
         +2^max(e,0) 3^max(b,0)
      =T_1+T_2.                                     (6.2)

Using |a|,|b|,c<=H and delta>=1:

* If b>=0, max(-e,0)+c<=2H, so T_1<=20^H.
* If b<=0, write t=-b. Then max(-e,0)+c=max(-a-t-(delta-1)c,c)<=H, so T_1<=30^H.
* If b>=0 and e>=0, T_2=2^(a+delta*c)(3/2)^b<=(3*2^delta)^H. If e<0, T_2<=3^H.
* If b<=0, max(e,0)<= (delta+2)H, so T_2<=2^((delta+2)H).

Therefore, for every integer H>=1,

    ||F_v||_1<=30^H+2^((delta+2)H)
               <=(30+2^(delta+2))^H
               =38^H for plus, 46^H for minus.       (6.3)

Apply this to the two short independent vectors of Section 5. By the coprimality lemma,

    R=Res(F_1,F_2)

is a nonzero integer. At least one polynomial has positive degree. The integral resultant Bezout identity gives R=A(X)F_1(X)+B(X)F_2(X) for integer polynomials A,B. Evaluating at the actual common modular root yields

    g divides R.                                    (6.4)

This is integer divisibility, not just divisibility by the radical of g. If one polynomial is a constant, its resultant with the other is the appropriate positive power of that constant, and the same divisibility holds directly.

Hadamard's inequality on the Sylvester matrix, followed by ||F||_2<=||F||_1, gives

    |R|<=||F_1||_1^d_2 ||F_2||_1^d_1
         <=C^(H_1 d_2+H_2 d_1)
         <=C^(4H_1H_2)
         <=C^(4p^(2/3)),                             (6.5)

where C=38 or 46. Together with (6.4), this proves (1.2). Finally 46<exp(4), so the uniform upper bound is strictly below exp(16 p^(2/3)).

## 7. Exact controls and why the resultant is not an iff gate

The producer `verify_overlap_lattice.py` uses no floating-point arithmetic or symbolic package. It checks every coprime shape with 2<=p<=40 and 1<=r<p, exhaustively compares the legal parameter ranges with actual edge-by-edge overlapping paths, and evaluates carries by two methods. It finds exact first and second sup-norm minima by integer shell enumeration. Cleared polynomials use explicit binomial coefficients; resultants use fraction-free Bareiss determinants of the Sylvester matrix.

The main finite universe has:

* 489 source shapes and 6,546 actual paths;
* 3,273 paths in each direction;
* 3,272 paths below r=p/2, so the noncritical legality correction is tested;
* 2,262 positive-clock and 4,284 negative-clock paths;
* 489 original-cut crossing paths.

Including explicit controls, the producer performs 470,590 checks, with byte-identical normal and optimized receipts.

### Partial-factor and dependent-projection hostiles

At p=5,r=3,u=1 in the plus family, take independent relation vectors

    (b,c,m)=(0,1,1), (5,1,1).

Their (m,c) projections are equal, not independent. The polynomials are

    5X(X-1)+2,    640X(X-1)+243,

and their nonzero resultant is 4225=25*13^2. Taking instead the first vector (5,0,0) gives the constant polynomial 13 and resultant 169 with the quadratic. Neither degenerate case may be dropped.

For this actual word, D=13 and g=1, whereas gcd(D,R)=13 in both controls. Thus resultant divisibility is necessary, not an exact formula for actual content. The simultaneous clock root remains relevant even after using two short relations.

### Prime-power depth controls

A nonendpoint plus chain at (p,r,u)=(18,5,2) has

    g=6647=17^2*23.

The producer additionally checks the endpoint-derived controls

    plus  (418,237,180): g=11849=17^2*41;
    minus (203,165,0):   g=169=13^2.

Full g, not just its radical, divides the nonzero resultant in every case. The p=2,r=1 negative cycle is retained: the plus path 011 -> 110 -> 101 has fixed point -7, and the minus path 011 -> 101 -> 110 has fixed point -5. The finite universe has no other integral endpoint or nonendpoint overlap output. This last observation is finite evidence only.

The independent audit additionally checks 175,670 legal overlap paths through p<=120, all 4,759 lattices with p<=36 and every 0<=u<p-r, and 9,518 resultants. Its methods and receipts are in `audit/`.

Reproduce from this directory:

    python3 verify_overlap_lattice.py > receipt.replay.json
    python3 -O verify_overlap_lattice.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json

## 8. Separate finite-weight consequence, not a full classification

The core theorem above needs Minkowski's second theorem, but not Matveev. For this section only, import the cited Matveev bound and exact eighteen-convergent absolute-gap certificate from the [audited endpoint closure](../collatz_two_swap_endpoint_exclusion_20261020/ENDPOINT_TWO_EDIT_CLASSIFICATION.md), Sections 3--4. Those sections require only a hypothetical integral word with gcd(p,q)=1, p<q<2p, and

    |D|/3^p<exp(-p/9).

A positive rational Taylor sum certifies exp(23/6)>46, hence log46<23/6. For p>=4096, p^(1/3)>=16. The content theorem therefore gives

    log(g/3^p)
      <=4 log46 p^(2/3)-p log3
      <(23/24-13/12)p
      =-p/8<-p/9.                                   (8.1)

If the word were integral, |D|=g, so the hypothesis of the imported absolute-gap certificate would hold. Its Matveev step gives p<10^12; its Legendre step and complete rational convergent check exclude every such p>=100 for both signs. Therefore:

**Finite-reduction corollary.** Any integer cycle among these primitive k=1 two-one overlap families has p<=4095.

This does not assert that the remaining all-clock residual has been checked. The source universe still includes every coprime 1<=r<p and every parameter in (1.1), not just the critical clocks and not just the two endpoint values. The two genuine p=2 negative paths remain inside it. Other overlapping operation types and arbitrary two edits remain outside this theorem's scope.

The unreduced finite universe 2<=p<4096 contains exactly 5,097,971 source shapes and 6,959,303,782 ordered paths. `count_unchecked_residual.py` verifies these counts both directly from the legality intervals and by Euler phi plus Mobius inversion of weighted coprime residue sums. It tests no carries. Its count receipt is expressly not an integrality census or a completed exclusion of that residual.
