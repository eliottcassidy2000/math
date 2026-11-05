# Odd-prime radix constraints for mixed fixed-port periods

## Exact scope and theorem

Start with any nonempty finite tournament and k distinct physical ports. At a right birth, the newborn beats a selected subset of those ports and loses to every other old vertex. At a left birth, it loses to the selected ports and beats every other old vertex. Keep exactly the original ports forever: no port is replaced, no newborn becomes a retained port, and no forgotten old vertex is re-exposed.

Fix a nonempty word W of these births and repeat it. Let h_n count directed Hamilton paths after n periods. The theorem permits an arbitrary seed, arbitrary fixed masks, either or both directions, and any finite k.

**Theorem.** Let p be an odd prime. Suppose the sequence (h_n) is not eventually constant and, for all sufficiently large n, satisfies the recurrence

`h_(n+1)=p*h_n+c`.

Then, with n still indexing periods from the original seed,

`h_n=A*p^n+beta` eventually,

where A is a positive integer divisible by `(p-1)!`, and beta is an odd integer. Moreover:

1. `p<=k+1`;
2. one birth direction occurs exactly once in W;
3. `c=-(p-1)*beta`, so `c=p-1 modulo 2(p-1)`.

Thus no such period has zero offset. At p=3, offsets must be2 modulo4. At p=5 they must be4 modulo8, and at p=7 they must be6 modulo12. These are necessary conditions, not existence claims for every permitted offset.

The statement is not about port-changing periods, adaptive choices of words, internal arc edits, or graph rebuilding. It does not assert a universal scalar-spectrum obstruction.

## Uniform finite interface and spectral reduction

We recall the needed all-order mechanism from the separately audited fixed-mixed interface theorem. For an unordered spanning path cover, write c_0 for its component count and P for its total port endpoint occurrences, counting a singleton port twice. Retain covers with `c_0-P<=1`, and record each component's endpoint port labels or a nonport symbol. There are finitely many such shapes because `c_0<=2k+1`.

A birth inserts its vertex between optional predecessor and successor components. For a right birth, any successor starts at a selected old port; for a left birth, any predecessor ends at a selected old port. Deleting an internal newborn increases component count by one while restoring at least one port endpoint. Endpoint or singleton deletion also cannot increase `c_0-P`. Therefore the retained covers are exactly closed under every birth, by a deletion/insertion bijection, and the Hamilton count is a linear readout of this finite nonnegative integral interface.

For a right transfer use the potential

`2k-P+(number of nonport starts)`.

It lies between0 and2k+1 and never decreases. Equality occurs precisely when the newborn extends a nonport end; the shape then stays unchanged. Hence its diagonal is the number d of nonport ends, with `0<=d<=k+1`, independently of the mask. A right-only word of length r has eigenvalues among `d^r`. The finite characteristic polynomial is an annihilator independent of the initial graph. Left-only words have the analogous statement using ends and starts interchanged.

Opposite-side births commute physically, fixing the original graph and all ports: the left newborn beats the right newborn in either chronological order. Thus their cover matrices commute. If W has r right and l left births, its period eigenvalues lie among `d^r e^l` when both sides occur, or among the corresponding one-sided bases when a side is absent. This follows by simultaneous triangularization of two commuting complex matrices, or by the two uniform annihilators.

An eventual nonconstant affine recurrence has the displayed form with A>0: a negative amplitude contradicts positive Hamilton counts for large n. Also c is an integer since consecutive counts are integers. The period recurrence forces p to be an eigenvalue. Primality then forces at least one of r,l to equal1; if the other direction is absent, the period length itself is1. It also gives p<=k+1. This reasoning is valid after an arbitrary finite transient.

## The factorial amplitude

For a fixed rooted tournament, let a_j count root skeletons with j root occurrences and no consecutive roots. Replacing the root by a transitive module of m vertices gives

`H_m=sum_j a_j j! S(m,j)`.

Here each skeleton lifts by a surjection of the ordered clone vertices onto its j ordered runs. The j-1 nonempty fragments between consecutive roots can be freely permuted, keeping any outer fragments fixed. Their vertex sets are disjoint and nonempty, so this action is free. Consequently `(j-1)!` divides a_j.

If only one birth occurs per period, the newborns form a transitive clone module. Positivity excludes all run indices above p under growth O(p^n). The p^n coefficient is therefore a_p, divisible by `(p-1)!`. This argument starts at n>=1 and includes any eventual recurrence.

Now assume both directions occur. Reverse every arc if necessary, so the left subword consists of one birth and the right subword has length r>=1. Commuting opposite sides shows that after n periods the n left-born vertices form one transitive clone module over the graph containing the seed and the n right subwords. Write a_j(n) for its skeleton coefficients. Then

`h_n=sum_(j=1)^(k+1) a_j(n) j! S(n,j)`.

The bound is uniform: each of the j-1 internal fragments ends at a distinct selected old port, since the left root loses only to those ports. In particular there are only finitely many possible j even though the right-side source grows. All coefficients are nonnegative integers and satisfy `(j-1)! | a_j(n)`.

For each fixed q, the count H_q(n) with q left clones and n right subwords is a right-only continuation from a fixed physical source. Its uniform recurrence has bases d^r. The triangular Stirling system at q=1,...,j expresses a_j(n) as a fixed rational combination of the H_q(n). Hence a_j(n), apart from a finite zero-eigenvalue transient, is an exponential-polynomial sequence with bases d^r. This recurrence does not depend on the growing order of the graph.

Positivity of the run sum and h_n=O(p^n) gives

`a_j(n)=O((p/j)^n)`.

For each of the finitely many j>p, integrality forces eventual zero. For j=p, the sequence is bounded, and its only possible bounded nonzero base is1. Thus a_p(n) is eventually a constant a, divisible by `(p-1)!`.

To identify the p^n coefficient exactly, expand

`j! S(n,j)=sum_(i=1)^j (-1)^(j-i) binom(j,i) i^n`.

Every resulting exponential base is `i*d^r`. Distinct positive exponential bases, with their polynomial factors, are linearly independent on every tail. Since p is prime, a p^n term has only two possibilities:

- i=p and d=1. Since j<=p eventually, this is precisely the term j=p and contributes a.
- i=1 and d^r=p. This requires r=1 and d=p. For j>=2, the bound on a_j(n) excludes a p^n component. Thus only a_1(n) can contribute.

In the second case a_1(n)=H_1(n) is itself an ordinary fixed-root transitive-clone sequence: there is one fixed left vertex and one right birth per repetition. Its growth O(p^n) excludes skeleton indices above p, and its p^n coefficient b is divisible by `(p-1)!` by the free-fragment action. No polynomial multiple of p^n appears in this literal clone formula.

Therefore A=a when r>1 and A=a+b when r=1. Both are integers divisible by `(p-1)!`. This also covers a_p=0 or b=0; nonconstant growth ensures their total is positive. All statements were made on tails, so no initial transient changes the conclusion.

Finally p is odd and `(p-1)!` is even. Every tournament has an odd Hamilton count. Subtracting the even integer A*p^n shows that beta is an odd integer. The claimed congruence for c follows.

## Physical examples and nonclaims

Two strong physical p=3 examples show that either offset sign is possible.

- Rows `[2,12,145,149,99,143,175,19]`, root7, fixed ports `[0,1,4]`, repeated `R(mask7,drop3)`: counts `60*3^m-1`, m>=1, returning `3H+2`.
- Rows `[30,64,10,66,174,335,21,367,95]`, root0, fixed ports `[1,2,3,4]`, repeated `R(mask15,drop4)`: counts `72*3^m+1`, m>=1, returning `3H-2`.

The source run vectors are respectively `(a0,a1,a2,a3)=(59,179,180,60)` and `(73,217,216,72)`. The clone formula proves both identities at every m. With n counting births from the displayed seeds, their amplitudes are180 and216, both divisible by2!, and their constants are-1 and+1.

A genuinely mixed +2 example inserts a universal-source birth `L(mask0)` before each negative-clock right birth, retaining all original ports. It has the same counts as that negative clock, with a growing transitive source prefix. It is not strong after a nonempty period, which is allowed by the theorem.

No physical p=5 or p=7 affine clock is claimed here. The congruence examples for those primes illustrate necessary arithmetic restrictions only.

## Verification

The mathematical proof is all-order and all-k. The standalone verifier checks free-fragment factorial divisibility by literal skeleton enumeration on all5,405 rooted labelled tournaments through order5, compares those coefficients with a separate subset recurrence, and checks the two concrete clock vectors and repeated physical births. Its arithmetic controls verify the prime factorization step over a finite grid, as a debugging aid rather than the proof of that step. Normal and optimized Python receipts agree.

The separately frozen foundation is `fixed-mixed-interface/PROOF.md`, SHA256 `95bc75e7e42d9212aa57c850bca9580958e38a88b0f63530bcaef3cf7c33b767`. This corollary does not alter that snapshot. The odd-prime extension was suggested during the independent ribbon review and is proved here with its own eventual and parity cases explicit.
