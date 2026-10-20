# All-order integrality of the odd and even-level switching branches

Status: PROVED using the cited Frobenius root-count divisibility theorem. This establishes the numerical integrality assertion in HYP-9172(b) for every n >= 3, including all positive n divisible by 4. It does not supply the requested structural combinatorial pairing. The n = 2 exception is sharp: both branches are 1/2.

## Inherited objects and exact question

THM-479, `01-canon/theorems/THM-479-level-holonomy-and-branch-split.md`, defines

    N_odd(n) = sum_{mu partition of n, all parts odd} 2^(o_2(mu)-k(mu)+1)/z_mu,
    N_lev(n) = A049313(n) - N_odd(n).

Here k(mu) is the number of parts, z_mu = product_j j^(m_j) m_j!, and

    o_2(mu) = sum_i floor(l_i/2) + sum_{i<j} gcd(l_i,l_j).

The total is the integer number of isomorphism classes of tournament switching classes. THM-4531, `01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md`, Theorem B, already provides a specific pairing for n = 2 modulo 4, n >= 6. The present argument is different and applies to every order.

## Theorem

For every integer n >= 1 except n = 2, both N_odd(n) and N_lev(n) are integers. At n = 2 their values are both 1/2.

### Step 1: the only possible denominator is a power of two

For any finite group H, let a(H) be the number of its elements of odd order, including the identity, and put rho(H) = a(H)/|H|. Write |H| = 2^s m with m odd. By Lagrange's theorem, the elements of odd order are exactly the solutions of h^m = 1. Frobenius's root-count theorem says their number is divisible by m. Therefore

    rho(H) belongs to Z[1/2].

Apply the orbit-stabilizer double count to the S_n-action on labelled tournament switching classes Q. Every orbit with stabilizer H contributes a(H)/|H| to the restricted Burnside sum over odd-order permutations. Thus

    N_odd(n) = sum_{[q] in Q/S_n} rho(Stab(q)) belongs to Z[1/2].

This equality uses THM-479's odd-branch fixed-class formula; it does not presume integrality of either branch.

### Step 2: each odd-partition summand has an odd denominator

Let mu have k odd parts. Since every gcd(l_i,l_j) is at least 1,

    e(mu) := o_2(mu)-k+1
           >= (n-k)/2 + binomial(k-1,2).

The odd factors j in z_mu contribute no 2-adic valuation. Moreover product_j m_j! divides k! (the quotient is a multinomial coefficient), so

    v_2(z_mu) = sum_j v_2(m_j!) <= v_2(k!).

We check e(mu) >= v_2(z_mu):

- k = 1: z_mu is odd and e = (n-1)/2 >= 0.
- k = 2: if n != 2, then n >= 4 and e >= (n-2)/2 >= 1 >= v_2(2!).
- k = 3: e >= 1 = v_2(3!).
- k >= 4: binomial(k-1,2) >= k-1 >= v_2(k!).

The sole exceptional partition is mu = (1,1), at n = 2. For all other orders, every summand 2^e/z_mu lies in Z_(2), the rationals with odd denominator. Hence their sum does also.

### Step 3: intersect the local conclusions

The intersection Z[1/2] intersect Z_(2) is Z. Therefore N_odd(n) is an integer for n != 2. Subtract it from the integral orbit count A049313(n) to obtain integrality of N_lev(n). At n = 2, direct substitution gives N_odd = N_lev = 1/2. This proves the theorem.

## Cited dependency, checked from a primary research source

Frobenius's theorem is: in a finite group H, the number of solutions of x^m = 1 is divisible by gcd(m,|H|). A modern primary research paper reproduces this statement explicitly on printed page 437 (PDF page 2), and proves a generalization containing it:

Elena K. Brusyanskaya, Anton A. Klyachko, Andrey V. Vasil'ev, *What do Frobenius's, Solomon's, and Iwasaki's theorems on divisibility in groups have in common?*, Pacific Journal of Mathematics 302 (2019), 437-452, DOI 10.2140/pjm.2019.302.437.

https://msp.org/pjm/2019/302-2/pjm-v302-n2-p03-s.pdf

This is the root-count divisibility theorem of Frobenius (1895), not the unrelated normal-complement theorem for Frobenius groups and not the stronger Frobenius conjecture about exactly m roots.

## Scope, hostiles and remaining question

- The proof is uniform; no finite computation is load-bearing.
- It does not say rho(H) is integral individually: C_2 gives 1/2, A_4 gives 3/4, PSL(2,7) gives 5/8.
- The two local arguments genuinely need each other. Partition summands can have nontrivial odd denominators, while individual stabilizer contributions can have dyadic denominators.
- At n = 8 the known 19 classes with Sylow 2-subgroup C_2 still forbid a fixed-point-free pairing contained within that Sylow stratum. The present proof does not remove that obstruction or define a pairing mixing the strata.
- HYP-9172(b) should therefore be split: branch integrality PROVED; an explicit combinatorial pairing explaining it at n = 0 modulo 4 remains OPEN.
- The abstract method generalizes: for any finite group action, the odd-order restricted Burnside average lies in Z[1/2]. Any independent 2-integral class-sum formula upgrades it to an integer.

## General prime-local restricted-Burnside criterion

Let a finite group G act on a finite set X and let p be any prime. Define

    B_p(X,G) = (1/|G|) sum_{g in G, p does not divide ord(g)} |X^g|.

For each stabilizer H, Frobenius's theorem applied to the p-prime part m of |H| implies that its proportion of p-regular elements lies in Z[1/p]. Orbit-stabilizer double counting therefore gives B_p(X,G) in Z[1/p].

If, additionally, every p-regular conjugacy-class representative g with |X^g| != 0 satisfies

    v_p(|X^g|) >= v_p(|C_G(g)|),

then the conjugacy-class formula

    B_p(X,G) = sum_{[g], g p-regular} |X^g| / |C_G(g)|

is termwise p-integral (zero terms cause no issue). Hence B_p(X,G) belongs to Z[1/p] intersect Z_(p) = Z. The switching-class application is p = 2, with |X^g| = 2^e and |C_{S_n}(g)| = z_mu. This is a sufficient criterion, not a necessary one.
