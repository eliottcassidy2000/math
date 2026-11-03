# Path-cover polynomial, odd-cycle packings, and 2-adic gates

**Status.** Exact consequences of the odd-cycle collection formula (OCF) in the repository's `06-writeups/research-document-for-mathematicians.md`, and elementary path-cover bijections. These results give arithmetic constraints and explicit families of strong tournaments, not a characterization of the global Hamilton-path spectrum. No novelty claim is made.

Let T be a tournament on n>=1 vertices and h(T) its Hamilton-path count. Set h(empty)=1. Define F_T(r) for nonnegative integers r as the number of ordered r-tuples of directed paths (empty paths allowed) covering the vertices of T exactly once. Also define F_empty(r)=1. Then

    F_T(r) = sum_{k=0}^n binom(r,k) C_k(T)
           = sum_{f:V(T)->[r]} product_{i=1}^r h(T[f^{-1}(i)]).

For nonempty T, F_T(0)=0, F_T(1)=h(T), and F_T(2)=2h(T)+C_2(T), so this is the full polynomial extension of the S158 two-slot state. The formula displays a polynomial in r of degree n.

## Cycle-packing expansion

We use the OCF as an input: h(U) is the sum over all collections of pairwise vertex-disjoint directed odd cycles of U of 2^(number of cycles); cycles have length at least 3, and the empty collection contributes 1. In the expression above for F_T(r), a given collection K of such cycles can be assigned to slots in r^c ways, c=|K|; each of the u=n−sum_{C in K}|C| uncovered vertices can also choose a slot independently. Therefore

    F_T(r) = sum_K 2^c r^(u+c)
           = sum_K 2^|K| r^[n−sum_{C in K}(|C|−1)].          (1)

Every exponent has the same parity as n, so F_T(−r)=(−1)^n F_T(r). All coefficients in the power basis are nonnegative integers. In particular, the r^(n−2) coefficient is twice the number of directed triangles; the r^(n−4) coefficient is twice the number of directed 5-cycles plus four times the number of disjoint pairs of directed triangles. Each directed cycle is counted once up to cyclic rotation, not both orientations unless both exist (impossible in a tournament).

A second self-contained derivation of reciprocity uses d(pi), the number of backward consecutive edges in a vertex permutation pi:

    F_T(r) = sum_pi binom(n+r−1−d(pi),n).

A slot assignment is a weakly increasing sequence of slot labels along pi, strictly increasing at backward edges. Reversing pi changes d into n−1−d; the generalized-binomial negative-argument identity proves reciprocity.

## Two-slot congruences and sharpness

At r=2 each packing contributes exactly 2^(u+2c). Consequently, with cycles counted up to cyclic rotation:

* F_T(2) is divisible by 4 for every n>=2.
* For even n>=4, F_T(2) = 8 times the number of directed (n−1)-cycles mod 16. In particular it is divisible by 8.
* For odd n>=5, F_T(2) = 4 times the number of directed Hamilton cycles mod 16.
* At n=3 the same last formula has an additional 8 from the empty packing: F_T(2)=8+4*(number of directed 3-cycles).

For the proofs modulo 16 only packings of exponent u+2c<4 matter. The empty packing contributes 2^n, absent for n>=4. Any nonempty packing with exponent at most 3 is one cycle with u<=1; this cycle is spanning when n is odd, or misses one vertex when n is even. This establishes the displayed cases and handles the n=3 exception.

The bounds are best possible for every order. Define X_m (m>=3) by orienting i->j for i<j except m->1. Every directed cycle uses the reversed edge and its remaining vertices occur in increasing order, so on each odd vertex subset containing {1,m} there is exactly one such cycle; all odd cycles share two vertices. Equation (1) yields

    F_{X_m}(r) = r^m + 2 sum_{j odd>=1} binom(m−2,j) r^(m−j−1),
    F_{X_m}(2) = 2^m + 2(3^(m−2)−1).

For odd m>=3 this is 4 mod 8. For every even n>=4, let T=X_{n−1} ordinal-sum a singleton. The slotwise ordinal-sum bijection gives F_T(2)=2 F_{X_{n−1}}(2)=8 mod 16. Thus universal divisibility by 4 in odd orders and by 8 in even orders cannot in general be improved.

More generally, for s>=1, (1) expresses F_T(2^s) as a positive sum of powers 2^[s*u+(s+1)c]. In particular the minimum of that exponent over directed odd-cycle packings gives a certified lower bound on its 2-adic valuation. Equality additionally requires an odd number of minimizing packings; no equality is claimed in general.

## Exact arithmetic operations and a strong gate

For ordinal sum X⊕Y, slotwise restriction and concatenation give

    F_{X⊕Y}(r) = F_X(r) F_Y(r).

Let G(T) be the 3-cycle substitution with T and two singletons. It is strong for every nonempty T, and S158 gives h(G(T))=h(T)+F_T(2). Counting r-slot covers according to whether the two singleton vertices occupy distinct or the same slot yields

    F_{G(T)}(r) = r^2 F_T(r) + r[F_T(r+1)−F_T(r−1)].       (2)

Indeed distinct slots give r(r−1)F_T(r). Identical slots give r[F_T(r)+F_T(r+1)−F_T(r−1)]: the last term corrects the case where that slot contains neither T-vertex. These hold at all nonnegative integers r, with F_T(0)=0 for nonempty T, and hence as polynomial identities. The first repeated wrappers of a singleton have polynomials r^3+2r, r^5+8r^3+6r, and r^7+18r^5+74r^3+30r, giving strong Hamilton-path counts 3,15,123.

Thus multiplication and the finite-difference operator (2) form an exact arithmetic-circuit calculus on the entire path-cover profile. Its mod-4 obstruction is immediate: h(G(T))=h(T) mod 4 when n>=2. If n is even and >=4, h(G(T))=h(T) mod 8. For an odd-order X_m above, the difference is 4 mod 8; this shows the latter stronger congruence is special to the even-order input condition.

## Independent finite control

`python3 check.py` and `python3 -O check.py` enumerate all 1,099 labeled tournaments of orders 1 through 5. They count h on subsets by endpoint DP, separately enumerate directed odd cycles and disjoint cycle collections, and verify (1) against direct r-slot subset partition counting at r=1,2,3,4. They check parity of the power exponents and both mod-16 statements in their stated ranges. These computations do not replace the OCF or the general bijections.

## Smith normal form of the sample lattice

Fix n>=1 and let d=ceil(n/2), ε=1 if n is odd and ε=2 if n is even. Every F_T(r) lies in the free abelian group generated by r^ε,r^(ε+2),...,r^(ε+2d−2). Consider the d-by-d evaluation matrix M with rows i=1,...,d and columns k=0,...,d−1, M[i,k]=i^(ε+2k). Its Smith invariant factors are

    ε=1: 1!, 3!, ..., (2d−1)!,
    ε=2: 2!/2, 4!/2, ..., (2d)!/2.

Proof. The central-factorial polynomials

    g_k(r) = r^ε product_{j=1}^k (r²−j²),    k=0,...,d−1,

are a unimodular change of integer basis. They vanish at the first k sample points. At row i they evaluate, respectively, as

    ε=1: g_k(i) = (2k+1)! binom(i+k,2k+1),
    ε=2: g_k(i) = ((2k+2)!/2) [binom(i+k,2k+1)+2 binom(i+k,2k+2)].

The bracketed quotients (also for the first displayed line) form an integer unit lower triangular matrix; integer row operations therefore diagonalize M to the stated factorial sequence. Each term divides the next, so this is its Smith form.

Consequently the sample vectors of parity-restricted integer polynomials form a full-rank lattice with exactly the stated invariant factors. Membership gives necessary and sufficient congruence conditions for *integer-polynomial interpolation within this parity class*, not for a tournament path-cover profile. A genuine tournament profile satisfies stronger restrictions: its leading coefficient is 1, every other coefficient is even and nonnegative, and (1) requires those coefficients to arise from compatible cycle packings. The Smith calculation must not be used as a sufficiency claim for tournament realization.
