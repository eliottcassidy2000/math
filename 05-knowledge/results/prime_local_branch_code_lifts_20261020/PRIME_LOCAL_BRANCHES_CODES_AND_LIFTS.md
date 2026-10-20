# Prime-local branch integrality, Pascal codes, and prime-power lift clocks

Date: 20 October 2026. Status: **PROVED scoped algebraic results; CITED classical code-distance and Frobenius inputs; FINITE-EXACT controls.** No Collatz convergence result, unrestricted mixed-prime symmetry correspondence, or new classical coding-theory priority is claimed.

## 1. Inheritance and the connection contract

The branch-integrality mechanism is inherited from `05-knowledge/results/switching_branch_integrality_20261020/BRANCH_INTEGRALITY.md`: a p-regular Burnside sum has only p-power denominator, while its conjugacy-class expression can separately be p-integral. Their intersection is integral. This report extends the concrete carrier to oriented gain-switching classes over every prime and prime power, with an explicit sign-aware fixed-space formula.

The other two lanes are independent, precise uses of prime-power structure: Pascal rows as repeated-root cyclic codes over a field, and Jordan clocks over the ring Z/p^s. They are not the same object and are not implicitly interchanged. The concept board is invariant quotient / p-regular averaging / local denominator / repeated root / valuation depth / common witness. The anchor is all-prime branch integrality; the niche is a classical code identification with a closed distance staircase; the wildcard is the first prime-power Jordan layer where reduction loses period information.

The Collatz connection is only the *method* of retaining the actual prime-power obstruction. The map there is a parity word to the valuations of its carry and clock; here it is a finite group action to a rational restricted Burnside sum. These maps do not identify their objects. A finite code distance, an integral orbit contribution, or a finite-ring period supplies no sign-controlled descent on positive integers.

## 2. Oriented gain switching and the all-prime theorem

Let R be either F_(p^s) or Z/p^s, with p prime and s>=1. Write Q=|R|=p^s. On the complete graph with vertices 0,...,n-1, an oriented gain assignment is a collection a_ij in R with a_ji=-a_ij and a_ii=0. Switching by potentials u_i replaces

    a_ij by a_ij+u_j-u_i.

Let X_R(n) be the quotient of all assignments by all such potential differences. Fixing a_0i=0 identifies it with R^((n-1)(n-2)/2). The symmetric group acts by permuting vertices, with signs retained when edge orientation reverses. This is an ordinary finite S_n-set, not an unlabelled family of graphs with a silently chosen orientation.

For a partition mu=(l_1,...,l_k) of n, let

    z_mu=product_j j^(m_j) m_j!,
    e(mu)=sum_i floor((l_i-1)/2)+sum_(i<j) gcd(l_i,l_j)-k+1.

Define the p-regular branch and complementary branch by

    B_(p,s)(n) = (1/n!) sum_(g in S_n, p does not divide ord(g)) |X_R(n)^g|,
    C_(p,s)(n) = |X_R(n)/S_n|-B_(p,s)(n).

**Theorem 1.** For every n>=1,

    B_(p,s)(n)=sum_(mu partition n, p does not divide any part) Q^(e(mu))/z_mu.    (1)

Both B and C are nonnegative integers except when p=2,n=2, in which case B=C=1/2 for every s. In particular the formula and integrality hold over both rings of size p^s named above, although those rings need not be isomorphic and their p-singular fixed sets need not agree.

### 2.1 Proof of the fixed-space formula

Let E be the R-module of oriented edge assignments and V=R^n the vertex module. The sequence

    0 -> R -> V -> E -> X_R(n) -> 0

has constants as its first map and potential differences as its second. For a p-regular permutation g, its order is a unit in R. Averaging over its cyclic group makes the fixed-point functor exact.

If g has k vertex cycles, V^g has rank k, and the invariant potential-difference module has rank k-1. Within a vertex cycle of odd length l there are (l-1)/2 free edge-orbit parameters. Within a cycle of even length l there are (l-2)/2 such parameters; the antipodal orbit reverses orientation after half a turn and requires a=-a. When p is odd, 2 is invertible, forcing a=0. When p=2, p-regularity forbids even cycles entirely. Thus in every allowed case the within-cycle rank is floor((l-1)/2). Between cycles of lengths l_i,l_j there are gcd(l_i,l_j) free parameters, with no orientation-reversal constraint.

Therefore E^g is free of rank the first two sums defining e(mu). Exactness gives X^g as its quotient by a free invariant potential module of rank k-1, hence |X^g|=Q^e. One can obtain the cardinality directly from exactness without needing to assert that every quotient of free modules is free. Summing by cycle type yields (1).

### 2.2 The only possible denominator is a power of p

For any finite group action G on X, orbit-stabilizer double counting gives

    (1/|G|) sum_(g p-regular) |X^g|
      =sum_([x] in X/G) #{p-regular h in G_x}/|G_x|.

Write |H|=p^a m with gcd(p,m)=1. The p-regular elements of H are exactly the roots of h^m=1. Frobenius's root-count divisibility theorem makes their number divisible by m. Each ratio therefore belongs to Z[1/p], and so does their sum. This is an inherited general lemma, not a new claim about each stabilizer ratio being integral.

The imported theorem is explicitly stated on printed page 437 in Brusyanskaya–Klyachko–Vasil'ev, *What do Frobenius's, Solomon's, and Iwasaki's theorems on divisibility in groups have in common?*, Pacific J. Math. 302 (2019), 437–452, [primary PDF](https://msp.org/pjm/2019/302-2/pjm-v302-n2-p03-s.pdf), DOI 10.2140/pjm.2019.302.437. The cited input is gcd(|H|,m)-divisibility of the number of roots, not a normal-complement theorem.

### 2.3 Each cycle-type term is p-integral

Because the permitted l_i are prime to p,

    v_p(z_mu)=sum_j v_p(m_j!) <= v_p(k!).

Also all pairwise gcds are at least one, so

    e(mu) >= sum_i floor((l_i-1)/2)+binom(k-1,2).                    (2)

For odd p, k=1 or 2 gives v_p(k!)=0 and e>=0. For k>=3,

    v_p(k!) <= (k-1)/(p-1) <= (k-1)/2 <= binom(k-1,2) <= e(mu).

For p=2, every l_i is odd, and (2) becomes

    e(mu)>=(n-k)/2+binom(k-1,2).

For k=1 the denominator is odd. For k=2 and n!=2, n>=4 gives e>=1=v_2(2!). For k=3, e>=1=v_2(3!). For k>=4, binom(k-1,2)>=k-1>=v_2(k!). The only failure is mu=(1,1), n=2.

Since s e>=e, every nonexceptional summand Q^e/z_mu is p-integral. Thus B belongs to Z[1/p] intersect Z_(p)=Z. Subtracting from the integral full orbit count proves integrality of C. Both branches are nonnegative because they sum nonnegative fixed-point counts over complementary subsets of permutations. At p=2,n=2 the quotient X has one element and the p-regular subset consists only of the identity, giving 1/2 in each branch. ∎

### 2.4 Limits and a prime-power distinction

At p=2,s=1, this linear gain carrier is graph switching, not the affine tournament-switching torsor. Its regular branch equals the inherited N_odd, but its complementary branch must not be called N_lev. Already n=3 gives graph-switching total 2 and complement 1, whereas the tournament-switching total is 1 and N_lev=0.

**Affine extension, also PROVED.** Define an affine gain assignment by a_ji=1-a_ij for distinct vertices, with the same potential switching. Diagonal entries are not part of the data and satisfy no affine relation. Its quotient A_R(n) is a torsor for X_R(n). Permutations act affinely. Every p-regular permutation has a fixed point on this torsor: average any point over its cyclic orbit using the inverse of its order in R. This affine average is well defined because its coefficients sum to 1. Once a fixed point is chosen, the fixed set is a torsor for X_R(n)^g, so its cardinality is again Q^e. Therefore the regular branch is exactly (1) for A_R(n) too, and the complementary branch is integral with precisely the same p=2,n=2 exception. At R=F_2 this affine carrier is exactly tournament switching, so its complementary branch does recover N_lev. For odd p, shifting every off-diagonal gain by 1/2 gives an equivariant identification with antisymmetric gains, so both full orbit counts agree. When 2 is not invertible, no such shift is available.

The proof gives numerical integrality, not an involution, matching, or intrinsic decomposition of the orbit set into two subsets. Individual stabilizer terms can be fractions; C_p contributes 1/p. The inherited n=8 switching obstruction to pairing an unadorned carrier is not removed.

For a concrete ring/field distinction, X_(F4)(3) has four full orbits, but X_(Z/4)(3) has three. Both regular branches are 2. On this one-dimensional quotient the action is sign: it is trivial in characteristic two, but acts by x→-x modulo 4. Thus equal regular branches over equal-size alphabets do not identify their singular branches.

The theorem works for prime-power gain alphabets of a fixed characteristic prime because the p-regular averaging denominator is invertible. One must not substitute an arbitrary mixed modulus while retaining the same formula: a permutation prime to one prime may still have order noninvertible at another. CRT can decompose the underlying assignment module, but simultaneous permutation actions and fixed-set conditions must still be retained. Independent one-prime matches do not automatically yield a common mixed-prime match.

## 3. Pascal rows are a classical repeated-root code

Let P=p^m, m>=1, and 0<=F<=P. Over F_p define

    C_F=span { (binom(j,h) mod p)_(j=0)^(P-1) : 0<=h<F }.

**Theorem 2, classical identification and specialized consequence.** The dimension is F, cyclic shifts preserve C_F, and under the polynomial coordinate map it is the ideal

    C_F = <(x-1)^(P-F)> in F_p[x]/(x^P-1)=F_p[x]/((x-1)^P).        (3)

The Euclidean dual is C_(P-F). For F>0 its minimum Hamming distance is

    d(C_F)=min_(0<=h<F) product_(i=0)^(m-1) (p-h_i),               (4)

where h_i are the base-p digits of h. If F<P and L is the largest p-power not exceeding F, then

    d(C_F) = P/(pL) * (p-floor(F/L)+1).                            (5)

At F=P the distance is 1. The zero code F=0 is excluded from minimum-distance assertions.

**Attribution.** Equation (4) is not a new coding-theory theorem. It is a specialization of the repeated-root cyclic-code theory of Castagnoli, Massey, Schoeller, and von Seemann, *On Repeated-Root Cyclic Codes*, IEEE Trans. Inf. Theory 37(2) (1991), 337–342, [primary PDF](https://www.isiweb.ee.ethz.ch/archive/massey_pub/pdf/BI433.pdf), Theorem 1 and its preceding base-p weight formulas. An equally direct source is Camps, López, Matthews, and Sarmiento, *Polar Decreasing Monomial-Cartesian Codes*, Theorem 3.9, [author PDF](https://personal.math.vt.edu/gmatthews/Polar_decreasing.pdf), [DOI](https://doi.org/10.1109/TIT.2020.3047624). The present note supplies the explicit Pascal identification and elementary derivation below, with no priority claim.

### Proof of the identification and distance

The first F columns form an upper triangular matrix with diagonal 1, so the dimension is F. Pascal's identity sends the h-th row under j→j+1 to row h plus row h-1. It also holds across j=P-1→0 because binom(P,h)=0 mod p for 0<h<P. Thus the code is cyclic. Every ideal of F_p[x]/((x-1)^P) is generated by a unique power of x-1; its dimension fixes that power as P-F.

To prove duality directly, take h<F and k<P-F, so h+k<=P-2. By the base-p binomial digit identity, the inner product factors into m one-digit sums. At least one digit satisfies h_i+k_i<p-1, since otherwise h+k>=P-1. For that digit, binom(t,h_i)binom(t,k_i) is a polynomial of degree less than p-1. Its sum over F_p is zero, including the constant case. Hence the rows are orthogonal; dimensions F and P-F prove equality of the dual codes.

The same digit identity identifies C_F with polynomial evaluation on F_p^m using the basis product_i binom(x_i,h_i), h<F. The index set is downward closed under coordinatewise comparison. Triangular change from these binomial polynomials to monomials gives precisely the span of x_0^(h_0)...x_(m-1)^(h_(m-1)) with h<F.

For completeness, a nonzero polynomial of individual degrees below p with lexicographically leading exponent h has at least product_i(p-h_i) nonzero evaluations. Induct on the number of variables: its leading coefficient in the first variable has at least the corresponding product on the remaining grid, and above each such point a degree-h_0 polynomial vanishes at no more than h_0 of the p first coordinates. This proves the lower bound in (4). Equality is attained for any chosen h by product_i product_(a in A_i)(x_i-a), where A_i has h_i distinct field elements: all its monomials are in the downward-closed index set and it has exactly the asserted number of nonzero evaluations. Taking a minimizing h proves (4).

For (5), write L=p^r and F=aL+b with 1<=a<=p-1 and 0<=b<L. All digits above r vanish in h<F and contribute p^(m-r-1). The value h=aL-1 has digit r equal to a-1 and every lower digit equal to p-1, giving the factor p-a+1. If digit r is smaller, no lower product improves this; if it is a, then b>0 and the lower digits cannot all be p-1. Their product is at least 2, so the total low factor is at least 2(p-a)>=p-a+1. This proves the staircase. ∎

No absolute-difference dynamics is identified with addition at odd primes. For example |0-1| mod3=1 whereas 0-1 mod3=2; the binary identity |a-b|=a+b for a,b in {0,1} is characteristic-two-specific. A Pascal array shared with an integer recurrence is not itself a lossless dynamical map.

## 4. Jordan clocks over Z/p^s: exact order and lost lift data

Let N be the single nilpotent Jordan block of size F, with (Nx)_i=x_(i+1) and x_F=0. Set U=I+N over R=Z/p^s. Assume F>=2, and let m be the least integer with F<=p^m.

**Theorem 3.** The order of U is

    ord(U)=p^(s+m-1).                                              (6)

For F=1 the order is 1 instead.

**Proof.** For 1<=k<=p^a,

    v_p(binom(p^a,k))=a-v_p(k),

where the endpoint k=p^a gives valuation 0. For k<p^a this follows from binom(p^a,k)=(p^a/k)binom(p^a-1,k-1), whose second factor is a p-unit by its base-p digits. In U^(p^a)=sum_k binom(p^a,k)N^k, the relevant k satisfy 1<=k<F. Their largest p-adic valuation is m-1, attained at k=p^(m-1). Taking a=s+m-1 makes every coefficient divisible by p^s. Taking the immediately smaller exponent a-1 leaves the k=p^(m-1) coefficient with valuation s-1. Different N^k occupy different superdiagonals, so it cannot cancel. The upper bound is a p-power, making this lower test sufficient for exact order. ∎

The result is classical unipotent arithmetic, presented with proof; no novelty assertion is made. It is a statement about Z/p^s, not the field F_(p^s), where U has order p^m.

### First carry-sensitive layer: F=p+1, modulus p^2

**Theorem 4.** For U on (Z/p^2)^(p+1), the numbers of cycles of exact lengths 1,p,p^2,p^3 are respectively

    p^2,
    p^(p+1)-p,
    p^(2p-1)-p^p,
    p^(2p-2)(p-1).                                                (7)

No other lengths occur. A point x=(x_0,...,x_p) is fixed by U^p exactly when

    x_i=0 mod p for 2<=i<=p-1,
    x_p+p x_1=0 mod p^2.                                         (8)

The empty middle condition at p=2 is intentional.

**Proof.** The order from (6) is p^3. U fixes exactly p^2 points because x_1=...=x_p=0. In U^p-I the coefficients of N,...,N^(p-1) are p times units and the coefficient of N^p is 1. Solving its coordinates backwards first forces x_2,...,x_p to be multiples of p; the first coordinate then reduces to x_p+p x_1=0 mod p^2. Conversely those conditions annihilate every coordinate. Thus U^p has p^2 choices each for x_0,x_1, p choices for each of x_2,...,x_(p-1), and x_p determined: p^(p+2) fixed points.

In U^(p^2)-I only the N^p coefficient survives modulo p^2, and it is p times a unit. Thus x_p is any multiple of p and all other p coordinates are free, giving p^(2p+1) fixed points. U^(p^3) fixes all p^(2p+2) points. Successive differences of these fixed counts divided by the relevant period give (7). ∎

For p=3 the complete cycle census is lengths 1,3,9,27 with counts 9,78,216,162. The states (0,1,0,0) and (0,1,0,6) have identical reductions modulo 3 but periods 9 and 3 modulo 9. Reduction therefore loses lift-dependent cancellation data. Maximum matrix order also does not determine the period of each state.

### Mixed modulus: the correct CRT composition

If M=product p_i^(s_i), CRT conjugates U on (Z/M)^F to the product of its local systems. For a tuple of local states of periods t_i, the global period is lcm(t_i). A product of selected local cycles of lengths t_i contains product(t_i)/lcm(t_i) global cycles. This follows because the product set has product(t_i) points and every global orbit in it has the indicated lcm length. For this unipotent CRT family, each t_i is a power of its distinct characteristic prime p_i by Theorem 3 (or is 1). They are therefore pairwise coprime: product(t_i)=lcm(t_i), each selected tuple of local cycles gives exactly one global cycle, and the total global cycle count is the product of the local cycle counts. The full period distribution still needs the local lengths. For arbitrary products of permutation systems with overlapping period primes, the more general factor product(t_i)/lcm(t_i) can exceed 1; that contrast is not an obstruction in the present CRT family.

## 5. Exact controls and independent replay

`check_bridges.py` is standard-library-only and imports no repository mathematics. The recovered `controls.json` is reproduced byte for byte by a fresh run. It records:

- 38,130 enumerated Pascal-code words in six stated (p,P,F) universes, and 6,156 dual inner products.
- Complete p-regular branch tables through n=16 for p=2,3,5,7,11, and 99,102 direct fixed-state tests for the nine small (p,n) cases printed in the receipt.
- 660 matrix-order cases, p in {2,3,5,7,11}, 1<=s<=4, 2<=F<=34.
- Complete orbits at modulus 4,F=3 (64 states) and modulus 9,F=4 (6561 states), plus the identical-reduction/different-period hostile.

`independent_bridge_controls.py` separately uses modular row reduction for p-regular fixed spaces, direct ring state enumeration for p^s gains, explicit cyclic-code ideal generators, and direct local/global orbit traversal for CRT. It also independently enumerates both linear and affine gain actions, keeping their complementary branches separate. Its finite ranges and counts are printed in `independent-controls.json`. The general theorems are proved above; the controls do not replace the all-order arguments.

Reproduce from this directory:

    python3 check_bridges.py > controls.replay.json
    python3 -O check_bridges.py > controls.optimized.replay.json
    cmp controls.json controls.replay.json
    cmp controls.json controls.optimized.replay.json
    python3 independent_bridge_controls.py > independent-controls.replay.json
    cmp independent-controls.json independent-controls.replay.json

## 6. What transfers and what does not

- Collatz word -> carry/clock valuations preserves the exact integer-cycle predicate if the common word and full prime-power depth are retained. Passing to prime support loses depth; testing different primes on different words loses simultaneity.
- Gain quotient -> restricted Burnside sum preserves the specified weighted orbit contribution. Frobenius removes prime-to-p denominators; fixed-space rank pays the remaining p-part. The result is a number, not a canonical carrier-level pairing.
- Pascal coefficient vector -> digit-polynomial evaluation is an actual code isomorphism, preserving weight after the fixed coordinate identification. Characteristic and digit order are retained. This map does not transport an unrelated integer recurrence.
- A ring state -> its field reduction preserves a factor trajectory but not exact period. The missing coordinate is its p-adic lift; (8) displays the first cancellation explicitly.
- Mixed-modulus state -> CRT tuple is lossless. Local lengths determine individual global periods; because their prime supports are disjoint here, total cycle counts multiply even without retaining those lengths.

Next targets are a decorated combinatorial carrier explaining the regular-branch integers, explicit p-singular fixed-space formulas over F_(p^s) versus Z/p^s, and full state-period stratification beyond F=p+1 at modulus p^2. None can be inferred from the common word “prime power.”
