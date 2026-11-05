# Opposed apexes, odd subsets, and a forbidden-pattern barrier

Session 197, 2026-11-05. **All-order physical lower bound and construction-class obstruction. The full `{7,21}` Hamiltonian-path spectrum remains open.**

This is a separate strengthening of the session's regular-language entropy packet. That packet uses Busch's universal lower bound and gives the Fibonacci controller exponent `3 log(phi)/log(5)`. Here an explicit family of actual odd cycles improves the exponent to `log(phi)/log(2)<7/10`, and rules out every fixed forbidden binary pattern in a specified opposed-apex carrier. It does not require Busch's quantitative theorem. No claim of literature priority is made.

## 1. Inherited input and physical carrier

The inherited odd-cycle collection formula is

    H(T) = sum_Gamma 2^|Gamma|,                       (1)

where Gamma ranges over collections of pairwise vertex-disjoint directed odd cycles of length at least three, including the empty collection; cycles are counted modulo rotation. The [S196 symmetry packet, equation (1)](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/symmetry/PROOF.md) explicitly states this inherited Grinberg--Stanley identity and its source. The external reference is Darij Grinberg and Richard P. Stanley, [*The Redei--Berge symmetric function of a directed graph*](https://arxiv.org/abs/2307.05569), whose tournament power-sum expansion gives (1). This packet uses the established identity; it does not prove it anew.

Only the empty-collection and single-cycle terms will be needed:

    H(T) >= 1 + 2 * (# directed odd cycles in T).     (2)

Fix two apexes a->b. A word w on {1,2} adds an increasing chain x_1->...->x_m, with the opposed incidences

    1: a->x_i->b,       2: b->x_i->a.

This is the same physical convention as [S191 two-apex transfer](https://github.com/eliottcassidy2000/math/blob/8d04f8557ea50c985521f039d5f58aa335617c9a/05-knowledge/results/session191_tournament_growth_20261104/two_apex/TWO_APEX.md). It is the restriction to two opposed letters, rather than arbitrary four-letter words, that enables the stronger bound.

A word containing no 2 gives the transitive order a,chain,b and H=1. Every word containing a 2 is strong. Indeed b reaches a through that 2-vertex, and a->b already exists. The apex component reaches the first chain vertex and is reached from the last chain vertex, regardless of whether either endpoint letter is 1 or 2. The increasing chain then joins every vertex to that component.

## 2. A physical odd-subset injection

**Theorem 1.** Suppose w contains a 2, has length m, and has u letters equal to 1. Then

    H(T(w)) >= 2^m-2^u+1, if u>=1,
    H(T(w)) >= 2^m+1,     if u=0.                    (3)

In particular,

    H(T(w)) >= 2^(m-1)+1.                            (4)

**Proof.** Choose any odd-cardinality subset R of chain vertices containing at least one 2-vertex. Divide it into its 1-vertices S and 2-vertices U. The literal cyclic vertex word

    a, increasing(S), b, increasing(U), a            (5)

is directed. Each nonempty S-run begins at an out-neighbor of a and ends at an in-neighbor of b. If S is empty, use the actual edge a->b. The U-run is nonempty, begins at an out-neighbor of b, and ends at an in-neighbor of a. Each run is internally increasing, so all its edges exist irrespective of how the two subsets interleave in the original chain.

The cycle has length |R|+2, which is odd and at least three. Different subsets R have different vertex sets, so their cycles are distinct even modulo cyclic rotation.

There are 2^(m-1) odd subsets of all m chain vertices. If u>=1, precisely 2^(u-1) contain only 1-vertices; if u=0 none do. Thus (5) constructs respectively

    2^(m-1)-2^(u-1),       or       2^(m-1)

distinct directed odd cycles. Apply (2) to get (3). Finally u<=m-1 because a 2 is present, giving (4). QED.

All these cycles share both apexes. We deliberately do **not** multiply their weights as though they were vertex-disjoint. Only the separate single-cycle terms in (1) are used.

The uniform exponential base two in (4) cannot be increased. The permitted family `1^(m-1)2` is nearly transitive and has H=2^m+1, an inherited scalar family checked independently here. Thus for any gamma>2 and fixed positive c, H>=c gamma^m eventually fails in this carrier. No pointwise sharpness of the lower estimate (3) is asserted.

### The leading coefficient is asymptotically sharp too

For integers a,b>=1, a direct one-apex calculation gives

    H(T(1^a 2 1^b))=(2^a+1)(2^b+1).                 (5a)

To check it, delete the unique 2-vertex r. The remainder is transitive in the order `a_apex, all chain 1-vertices, b_apex`. Relative to that order, the out-neighbor indicator of r is `1 0^a 1^b 0`. Every directed cycle contains r. A directed odd cycle through r consists of an increasing remainder path from an out-neighbor to a later in-neighbor, with an even number of internal selected vertices.

There are exactly three endpoint categories: the initial out-neighbor to a vertex in the first zero block, contributing 2^(a-1) cycles; a vertex in the last one block to the terminal zero, contributing 2^(b-1); and the first vertex to the last, contributing 2^(a+b-1). No two cycles can be disjoint because all contain r. Equation (1) therefore gives (5a).

Taking a=b=t makes m=2t+1 and gives

    H(T(1^t 2 1^t))/2^(m-1)=(1+2^(-t))^2 -> 1.

Thus the leading constant in (4) is asymptotically sharp, even inside the no-22 language. This is a transparent one-apex calculation used as a boundary control, not a claim of a newly discovered scalar family.

## 3. Transfer from word growth to scalar growth

Let L be any family of binary words, not necessarily regular, with at most `C(1+m)^d lambda^m` words of length m, where lambda>=1. All nontrivial outputs in this opposed carrier are strong, while all-1 words contribute only one. By (4), an output with 1<H<=X has

    m <= 1 + log_2(X-1).                             (6)

Thus its distinct strong generator values B satisfy

    # {b in B:b<=X} = O(poly(log X) X^rho),
    rho = log(lambda)/log(2).                        (7)

The exact degree of the harmless logarithmic factor is not needed.

For completeness, the multiplicative extension from the session's entropy packet is elementary. If (7) holds and s>rho, then `sum_(b in B) b^(-s)` converges by dyadic summation. Since b>=3, the Euler product

    Z_B(s) = product_(b in B) (1-b^(-s))^(-1)

converges too. Expanding finite partial products and taking a monotone limit gives

    sum_(n in <B>) n^(-s) <= Z_B(s),
    #(<B> intersect [1,X]) <= Z_B(s) X^s.             (8)

The product is indexed by **distinct integer generator values**, not ordered construction words. Different factorizations yielding the same integer only overcount.

Ordered joins multiply H, and transitive singleton padding contributes one. Therefore, for every epsilon>0, the full ordered-join and padding closure of L has at most

    O_epsilon(X^(rho+epsilon))                       (9)

distinct counts through X. If lambda<2, choose rho+epsilon<1 to obtain natural density zero. Moreover (8) is finite at s=1, whereas the reciprocal sum of the primes diverges. Hence infinitely many odd primes are omitted. This does not claim that omitted primes are globally impossible or globally attained.

## 4. Every fixed forbidden binary pattern is an obstruction

**Theorem 2.** Fix any nonempty binary pattern p of length r. Take any family of opposed-apex words avoiding p, then permit arbitrary ordered joins and transitive padding. Its distinct scalar support has natural density zero and omits infinitely many odd primes.

One explicit bound is

    O_epsilon(X^(rho_p+epsilon)),
    rho_p = log(2^r-1)/(r log(2)) < 1.               (10)

**Proof.** Divide a length-m avoiding word into disjoint full r-blocks and a final shorter block. Every full block must differ from p. There are at most

    (2^r-1)^floor(m/r) 2^(m mod r)
      <= C_r ((2^r-1)^(1/r))^m

such words. This is only an upper bound: it ignores forbidden copies crossing the chosen block boundaries, which can only lower the true inventory. Apply (9) with lambda=(2^r-1)^(1/r)<2. QED.

The boundary r=1 is included. For p=2, the language has only all-1 words and no nontrivial strong output; its count set is just {1}. For p=1, it has all-2 words, an exponentially growing scalar ray but only one word at each length. Neither boundary requires an invented strong witness.

This theorem rules out any fixed forbidden local pattern, however long, in this carrier. It is not a theorem about arbitrary four-letter two-apex words, changing the apex set, or general tournament constructions. Allowing all binary patterns is necessary to escape this particular fixed-pattern obstruction, not sufficient to prove coverage.

## 5. Fibonacci controller: a stronger exponent and explicit constants

Let L avoid the pattern 22. Its length-m word inventory is F_(m+2), since the two-state transition matrix is

    [[1,1],[1,0]].

Its growth constant is phi=(1+sqrt(5))/2. Equations (7)-(9) give

    rho_*=log(phi)/log(2)=0.6942419136306174... .

In particular the product closure has O(X^(7/10)) distinct values. The comparison can be certified exactly: `phi^10=(123+55sqrt(5))/2<128`, since `5*55^2<133^2`.

A quantitative bound also follows. For s>rho_* set z=2^(-s)<1/phi. The strong words are all valid words except the all-1 word at each length, so

    sum_(m>=0) (F_(m+2)-1) z^m
      = z / ((1-z)(1-z-z^2)).                       (11)

Choose one physical word for each distinct strong H-value. Using (4) and then overcounting by every strong word gives

    sum_(b in B) b^(-s)
      <= 2^s sum_(m>=0)(F_(m+2)-1)z^m
       = 1/((1-z)(1-z-z^2)).                        (12)

Consequently, for the whole product support M,

    # (M intersect [1,X]) <= exp(C(s)) X^s,
    C(s)=1/((1-3^(-s))(1-z)(1-z-z^2)).              (13)

The same threshold is convergence of the weighted Fibonacci matrix z[[1,1],[1,0]], whose spectral radius is z phi. This is a direct bridge from a finite-state physical constructor to an arithmetic Euler product. The constant is not claimed sharp, especially near s=rho_*.

### Physical gap and escape, without Busch

If a valid word had H=11, (4) would force m<=4. Direct enumeration of the no-22 words through length four gives the set

    {1,3,5,9,13,15,17,29,33}.

There is no 11, and multiplying positive integer generator values cannot create the prime 11. This proves an all-order extra omission for the entire Fibonacci/product class using only the new bound and exact finite graphs.

Allowing the forbidden transition gives the actual word `122`, with outgoing rows

    [6,24,26,17,1], H=11.

The allowed word `212` gives H=13. Thus crossing the one-apex H=13 obstruction does not entail realizing the smaller H=11. These small graph counts are independently checked by the verifier.

### Independently bounded scalar image through 5001

The new bound gives m<=13 whenever 1<H<=5001. There are 1595 no-22 words through that length, of which 1581 contain a 2. They give exactly 157 distinct strong generator values at most 5001 and 414 values after arbitrary multiplication, including one.

This agrees with the separate entropy packet's weaker Busch-based length-14 census. Here the completeness cutoff depends only on the explicit odd-cycle injection. `bounded_image.json` records a literal physical word for every generator value and a product factorization for every value in the closure. This is an all-order bounded image **only for this specified class**.

## 6. Verification

Run `python3 verify.py` and `python3 -O verify.py`; compare stdout and both JSON outputs. There are no external dependencies and no removable assertions.

The verifier includes independent adjacency/terminal-removal DP and the inherited physical nonnegative slot transfer. It checks the actual cyclic edges in (5) for every binary word through length six and every qualifying odd subset, verifies the exact odd-subset inventory, compares the resulting lower bound with independent graph counts, checks finite forbidden-pattern inventories, and verifies the base-two sharpness family. It also repeats the existing 1530 general fixed-apex transfer/SCC controls and the exact Fibonacci image, using the stronger length cutoff.

The all-order conclusions come from the cycle injection, elementary word counts, and Euler-product proof. They do not extrapolate from the finite tests. Only the established OCF identity is imported for the quantitative lower bound.
