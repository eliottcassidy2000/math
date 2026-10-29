# Arithmetic boundary recurrences with the multiplier source retained

2026-10-29. **PROVED here:** exact arity recurrences and their minimal orders for uniform nonzero multiplier sources; a weighted-source annihilator; a minimal deterministic divisor carrier; the explicit 36-modulus source obstruction. **INHERITED:** the floor/boundary identity from [Arithmetic braids II](../../arithmetic_braids2_20260917_floor_reciprocity.md), pinned at `d3741fe6ed02f6072dbe6bfda2f3d04193d75ce0`. **FINITE-EXACT:** the accompanying independent residue, divisor-state, polynomial and matrix checks. No external-priority claim, Collatz conclusion, or tournament equivalence is made.

## 1. Recovering a less-used boundary rather than another analogy

The September floor-reciprocity note already identifies the correction in a product-floor sum as the number of nonzero-coordinate tuples whose product vanishes modulo n. It proves the all-arity identity but does not give the arity recurrence of that correction. Its hostile is the distinction between a reduced composite ring and a field: the cube defect vanishes at 6, while the product defect does not.

Here the retained sidecar is the **whole truncated prime-valuation vector of a partial product**, not merely whether that product is zero. The concept board is:

* floor boundary counts and exact zero-product incidence;
* CRT decomposition and local nilpotent depth;
* rational generating functions and minimal scalar recurrence order;
* arbitrary source weights versus named future multipliers;
* short scalar certificates versus a carrier closed under every future operation.

The resulting connection to recent shared-port research is an exact information distinction. Uniform repeated averaging has a short scalar recurrence. Naming the next multiplier can reveal a mode invisible to every term of that recurrence.

## 2. The exact uniform formula

Let n>=2 have prime factorization n=product p^a, and let

    Z_r(n) = #{(x_1,...,x_r) in {1,...,n-1}^r : n divides product x_i}.

The empty product is 1, so Z_0(n)=0. For each prime power p^a||n define the polynomial

    P_(p,a)(r) = sum_(t=0)^(a-1) r^(overline t)/(t! p^t).

Thus P_(p,a)(r)=sum binom(r+t-1,t)/p^t for positive integer r. The rising-factorial definition also specifies r=0 without ambiguous negative binomial notation.

For S a nonempty subset of the distinct prime divisors, put

    lambda_S = n product_(p in S) (p-1)/p,
    P_S(r) = product_(p^a||n, p in S) P_(p,a)(r),
    D_S = sum_(p^a||n, p in S) (a-1).

Then for every r>=0,

    Z_r(n) = (n-1)^r + sum_(S nonempty) (-1)^|S| lambda_S^r P_S(r).       (1)

Every lambda_S is an integer. This formula is not obtained by assuming that the nonzero-coordinate source factors independently under CRT. That assumption would be false: a globally nonzero residue can have zero components.

**Proof.** First allow all residues, including zero in every coordinate. At p^a the number of r-tuples whose product is NOT zero is

    phi(p^a)^r sum_(t=0)^(a-1) binom(r+t-1,t)/p^t.

Indeed each valuation s<a has phi(p^a)/p^s residue choices; total valuation t<a has binom(r+t-1,t) weak compositions. The formula at r=0 is also 1. Therefore the local zero-product count is p^(ar)-phi(p^a)^r P_(p,a)(r). The all-residue global zero count is their product, by CRT. Every tuple with at least one globally zero coordinate has zero product; exactly n^r-(n-1)^r such tuples exist. Subtract them and expand the product. The all-empty subset term n^r cancels, leaving (1). QED.

For example,

    Z_r(12) = 11^r - 8^r - 6^r(1+r/2) + 4^r(1+r/2).

At r=2 this gives 17, the inherited table's boundary count. The sequence at r=0,...,4 is 0,0,17,439,7425.

## 3. An exact minimal recurrence, not just an annihilator

Write E for the shift (Eu)_r=u_(r+1). If n is composite, the minimal monic constant-coefficient rational recurrence polynomial of Z_r(n), even allowing a recurrence only eventually, is

    Q_n(E) = (E-(n-1)) product_(S nonempty) (E-lambda_S)^(D_S+1).        (2)

Let w=omega(n) be the number of distinct prime factors, and D=Omega(n)-omega(n). Its order is

    degree Q_n = 2^w + 2^(w-1) D.                                    (3)

Thus Z_r(n), as a series in arity, is rational with denominator

    (1-(n-1)z) product_(S nonempty) (1-lambda_S z)^(D_S+1),

and that denominator is reduced. Arithmetic complexity and integer bit length are still charged; this does not factor n for free.

**Proof of minimality.** The lambda_S are pairwise distinct. If two subset products product(1-1/p) agree, cancel common primes and choose the largest prime q in the remaining symmetric difference. In a cross multiplication, q divides one denominator but none of the opposite numerator factors p-1, whose primes satisfy p<q. This is impossible. Moreover lambda_S<=n-n/p<n-1 for composite n, choosing any p in S. Each P_S has degree D_S and nonzero leading coefficient, explicitly

    product_(p^a||n, p in S) 1/[p^(a-1)(a-1)!].

For distinct nonzero bases lambda, the exponential-polynomial sequences lambda^r r^j are linearly independent, and a degree-d polynomial times lambda^r has minimal shift factor (E-lambda)^(d+1). One elementary proof applies successive difference operators E-mu to remove the other bases, then uses that E-lambda lowers the polynomial degree exactly one. This proves (2), including eventual minimality. Summing D_S over all nonempty subsets counts each a-1 exactly 2^(w-1) times, giving (3). QED.

For prime n the two terms of (1) cancel identically: Z_r(n)=0 and its minimal recurrence order is zero. The composite hypothesis in (2) is essential.

Examples: orders at n=4,6,12,36,60,72 are respectively 3,4,6,8,12,10. In particular, squarefree composite n has order 2^omega(n), whereas prime-power depth adds multiplicity rather than new spectral roots.

## 4. The floor sums inherit a different leading root

For r>=1 put

    F_r(n) = sum_(1<=x_1,...,x_r<n) floor(product x_i / n),
    alpha = n(n-1)/2.

The inherited antipodal boundary lemma gives

    F_r(n) = alpha^r/n - (n-1)^r/2 + Z_r(n)/2
           = alpha^r/n + (1/2) sum_(S nonempty) (-1)^|S| lambda_S^r P_S(r). (4)

For composite n its exact minimal recurrence is obtained from (2) by replacing the root n-1 by alpha. Since alpha>n-1 for n>2, no new cancellation occurs. The order is again (3). This statement is for the tail r>=1; the empty-product floor at r=0 does not obey the displayed antipodal formula.

For prime n>2, the floor recurrence has order TWO, even though the boundary count vanishes: F_r=alpha^r/n-(n-1)^r/2. Thus at 11 its roots are 55 and 10; at 13 they are 78 and 12. At 12 the roots are 66,8,6,4 with multiplicities 1,1,2,2. These 11/12/13 distinctions arise from exact arithmetic, not a numerical analogy. For n=2 every floor sum vanishes.

## 5. Arbitrary iid source weights: the valuation algebra

Assign a rational weight w_x to each residue x modulo n, including zero if desired. Weights can be signed for the algebraic statement; nonnegative weights are needed for probability language. Keep the same source distribution at every step, independently. Let v_p(x) be capped at a, and work in

    A_n = Q[X_p : p|n] / (X_p^(a+1)-X_p^a for each p^a||n).

The source polynomial is

    W = sum_(x mod n) w_x product_p X_p^min(a,v_p(x)).

A coefficient vector in the monomial basis retains weighted mass at every divisor/gcd state. Multiplication by W is exactly the next-step update. The coefficient of product X_p^a is the zero-product mass. Arbitrary initial partial-product masses can be retained in an initial polynomial V, giving the observable coefficient in V W^r.

The factorization X^(a+1)-X^a=X^a(X-1) yields

    A_n ≅ product_(S subset primes) Q[X_p:p in S]/(X_p^a:p in S).

Coordinates outside S are set to 1. In component S the constant term of W is

    lambda_S(w) = sum_(x : p does not divide x for every p in S) w_x.

In particular lambda_empty is the total source mass. The remaining term lies in the nilpotent ideal generated by X_p,p in S, whose (D_S+1)-st power is zero. Consequently every scalar observation of V W^r is annihilated by

    product_(S subset primes) (E-lambda_S(w))^(D_S+1),                 (5)

where D_empty=0. Its order is bounded by (3), independent of how noninvariant the source is. Equal eigenvalues or special source weights can reduce the actual order; (5) is not asserted minimal for arbitrary weights. For a varying source sequence W_1,W_2,..., retain the exact product V W_1...W_r; no fixed-coefficient recurrence is promised.

This algebra keeps CRT dependence correctly. It never replaces an arbitrary source by independent local weights. Nor does it identify source mass with an orbit-averaged mass.

## 6. Nine necessary states but an eight-order scalar recurrence

Under the menu "multiply by a named nonzero residue; ask whether the product is zero", a partial product is represented exactly by gcd(product,n), equivalently its capped valuation vector. For composite n there are tau(n)=product(a+1) reachable states, and this deterministic carrier is minimal.

**Proof.** Every proper divisor is reachable by itself, while n is reachable as the product of two proper factors. Two distinct valuation vectors u,v have, after swapping their roles if necessary, a coordinate u_p<v_p. Choose a completing multiplier whose p-valuation is a-v_p, and whose other valuations saturate v. It is a proper divisor of n because a-v_p<a at the distinguished coordinate. It makes v zero but does not make u zero. Thus any two states can be distinguished by one legal named multiplier. QED.

The exact dimension gap is also explicit. Write b_p=a_p-1. Expanding product(a_p+1)=product(2+b_p) gives

    tau(n)-degree Q_n
      = sum_(S subset primes, |S|>=2) 2^(w-|S|) product_(p in S) b_p.     (6)

Thus the uniform scalar recurrence has the full divisor-state dimension exactly when **at most one prime factor is repeated**. A genuine gap requires joint nilpotent depth at two distinct primes. The smallest such modulus is 2²3²=36. This is a precise reason that taking the square of the first squarefree composite creates a new compression boundary; the numerical occurrence of 6 alone supplies no theorem about other subjects.

At n=36 the deterministic state count is 9 but the minimal uniform-source scalar recurrence order is 8. Here is an exact hostile with nonnegative source masses, each of total mass 8:

    A: 2 copies of 3, 3 copies of 4, 1 copy of 6, 2 copies of 18;
    B: 3 copies of 2, 2 copies of 9, 3 copies of 12.

For every r>=0 they have the SAME zero-product total after appending r independently uniform COUNTED multipliers from {1,...,35}. They therefore also agree under the normalized uniform probability source. Yet appending the single named multiplier 2 gives zero mass 2 for A and 0 for B, a probability difference of 1/4 after dividing by 8.

An all-order certificate is small. Order gcd states as

    (1,2,3,4,6,9,12,18,36),
    h=A-B=(0,-3,2,3,1,-2,-3,2,0).

Let M_(d,e)=#{1<=x<36 : gcd(dx,36)=e}. Direct integer multiplication gives h M=12h, and h has zero coordinate at 36. Hence h M^r has zero final acceptance coordinate for every r. The checker verifies the matrix equation, not just a finite prefix. The named multiplier 2 does not preserve that invisible mode.

This is the precise connection to recent carrier results: an eight-order recurrence is an exact compressed scalar calculation, but it does not replace the nine-state representation when later multipliers can be named. No contradiction exists between recurrence efficiency and faithful future-operation memory.

## 7. A complete one-step probe bank by classical Möbius inversion

The hidden initial source masses can nevertheless be reconstructed without long continuation words. Suppose an arbitrary initial mass m_g is assigned to every divisor state g|n, and its total mass is known. For each divisor d>1, append the named multiplier n/d and measure the resulting zero-product mass. This is a legal nonzero multiplier, and its reading is

    P_d = sum_(g|n, d|g) m_g.

Set P_1 equal to the known total mass. Classical divisor-poset Möbius inversion gives the exact inverse

    m_g = sum_(d|n, g|d) mu(d/g) P_d.                                  (7)

Indeed gcd(g*n/d,n)=n iff d|g, so these are precisely the upper-divisor zeta sums. Substituting them into (7) leaves the coefficient sum_(e|h/g)mu(e), which equals 1 for h=g and 0 otherwise.

Thus **tau(n)-1 one-step scalar probes plus the known normalization reconstruct every initial divisor-state mass**. This is optimal among fixed banks of linear scalar observations on unrestricted normalized real masses: the normalized simplex has affine dimension tau(n)-1, so fewer linear observations have a nonzero kernel direction; a sufficiently small displacement of an interior positive mass vector along that direction remains a pair of distinct normalized nonnegative sources with identical readings. The statement concerns arbitrary distributions, not adaptive identification of a single discrete state, for which a different query count may suffice.

At 36, eight named probes plus normalization recover all nine masses even though the entire infinite uniform-continuation sequence leaves the displayed signed mode invisible. This uses classical Möbius inversion; the application and hostile comparison are the contribution here, not a new inversion theorem.

## 8. Verification and stopping boundary

Run `python check_zero_product.py` and `python -O check_zero_product.py`. Checks use exact integers and Fractions and remain active under optimization. They compare:

* direct tuple enumeration and floor sums against (1),(4) for small full universes;
* divisor-state dynamic programming against (1) for n<=100 and arities through 12;
* characteristic-polynomial recurrences for composite n<=100;
* full scalar Hankel rank against the proved minimal order for selected moduli;
* arbitrary rational/signed weighted sources against (5), including sources at zero;
* the all-order 36-modulus eigenvector certificate and its named-multiplier separation;
* the one-step distinguishing multiplier for every pair of divisor states at composite n<=100;
* the complete normalized one-step Möbius probe reconstruction for positive and signed initial masses at n<=100.

The proofs cover all stated sizes; finite tests are controls, not substitutes. We do not assert that this aggregate recurrence reconstructs prime factors from a short prefix, preserves arithmetic height, yields Collatz descent, or enumerates unrestricted tournaments. The bridge is the shared algebraic mechanism: local nilpotent depth creates repeated scalar spectral roots, while named future operations can still expose hidden joint data.
