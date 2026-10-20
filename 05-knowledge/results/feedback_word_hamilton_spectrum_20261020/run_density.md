# Addendum: positive Stirling expansion and a density-zero barrier for bounded-run constructions

**Status: PROVED, independently reproducible matrix controls.** This extends the fixed-feedback word automaton in `proofs.md`, which is frozen separately for audit. No novelty claim is made. The global odd spectrum remains OPEN.

## 1. Commuting idempotents behind a homogeneous run

Fix a feedback set of size `k`, a feedback-vertex permutation, and a letter `a`. Let `P_j(a)` be the matrix for assigning the current chain vertex specifically to gap `j`, for `j=0,...,k`. It is a partial deterministic transition with nonnegative integer entries.

The matrix acts only on the state of gap `j`. If an empty gap fails its first-vertex condition, its transition is zero. Otherwise it becomes nonempty with the appropriate right-end status. Every already nonempty gap is reset to that same right-end status. Consequently

`P_j(a)^2=P_j(a)` and `P_i(a)P_j(a)=P_j(a)P_i(a)`.

The letter matrix is

`M_a=sum_{j=0}^k P_j(a)`.

These identities also hold on the block diagonal sum over feedback permutations.

**Spectral consequence.** The commuting idempotents are simultaneously diagonalizable over the rationals. Thus `M_a` is diagonalizable, its eigenvalues lie in `{0,1,...,k+1}`, and

`product_{j=0}^{k+1}(M_a-jI)=0`.                              (1)

Every fixed-context family `P a^r S` therefore has an exact expression as a rational linear combination of `0^r,1^r,...,(k+1)^r`, with no polynomial-in-`r` factors. Here `0^0=1` in the matrix-power expansion.

The more useful statement keeps **positivity**, rather than merely the eigenvalues.

## 2. Exact positive Stirling expansion

Let `Surj(r,t)=t! S(r,t)` count surjections from an `r`-element set onto a fixed `t`-element set. Set `Surj(0,0)=1`, and `Surj(r,0)=0` for `r>=1`.

Expanding `(sum P_j)^r`, group all sequences of chosen gaps by their support. Idempotence and commutativity show that every sequence with support `J` has product `product_{j in J}P_j`. The number of such sequences is `Surj(r,|J|)`. Hence

`M_a^r=sum_{J subseteq {0,...,k}} Surj(r,|J|) product_{j in J}P_j(a)`.    (2)

For fixed prefix and suffix words, let `u_P` and `v_S` be their nonnegative integer initial/final vectors. Define

`c_t=sum_{|J|=t} u_P (product_{j in J}P_j(a)) v_S`.

Every `c_t` is a **nonnegative integer**, and the ordinary path count is exactly

`H(P a^r S)=sum_{t=0}^{k+1} c_t Surj(r,t)`.                  (3)

This is stronger than an arbitrary exponential-polynomial fit. It has no negative coefficient cancellation in its native basis.

For example, with two feedback vertices `x_0 -> x_1` and homogeneous letter `a=2` (so `x_1 -> v -> x_0`), the count for `r>=1` is `2^r+1=3 Surj(r,1)+Surj(r,2)`.

**Run-growth boundary.** For `r>=1`, either all `c_t=0` for `t>=2`, in which case the count is constant, or the count grows exponentially. In particular a single fixed-type run cannot generate a nonconstant arithmetic progression or any nonconstant polynomial in its length. The spectral description proves this asymptotically; (3) gives an effective bound without cancellation assumptions.

## 3. A logarithmic compression lemma for long runs

Let a represented tournament have positive count `h`. Consider one nonempty run of identical letters of length `r`. The prefix and suffix, including their current lengths, are fixed while applying (3).

Suppose some `c_t>0` for `t>=2`, and `r>=k+1`. To exhibit surjections onto `t` fixed gaps, assign the first `t` positions to these gaps in a fixed order and assign the remaining `r-t` positions to either of two chosen gaps. Therefore

`h >= Surj(r,t) >= 2^(r-t) >= 2^(r-k-1)`.

It follows that

`r <= k+1+floor(log_2 h)`                                   (4)

whenever any `c_t` with `t>=2` is positive.

If instead

`r > k+1+floor(log_2 h)`,

all those coefficients vanish, so (3) says that replacing the entire run by **one** copy of its letter preserves `h` exactly. The argument needs no decision about which automaton states are actually reachable beyond the nonnegative integer formula.

Applying this to each run in turn preserves the count at every step. Thus every H-value witnessed with `k` feedback vertices and at most `R` orientation runs also has such a witness whose every run has length at most `k+1+floor(log_2 h)`.

This is value-preserving compression. The density conclusion below does not require that the compression preserve a separate property such as strong connectivity; it counts the larger unrestricted presentation class. Vertex deletion also preserves the property of having no two disjoint triangles.

## 4. Fixed feedback plus bounded runs is a density-zero numerical grammar

For fixed `k,R`, let `A_(k,R)` be the set of ordinary path counts of tournaments admitting a feedback set of size `k` such that the ordered transitive remainder has at most `R` maximal constant-letter runs. The orientation among the feedback vertices is unrestricted.

For `X>=1`, put `L=k+1+floor(log_2 X)`. There are at most `2^(k(k-1)/2)` feedback tournaments, at most `(2^k)^r` letter patterns with `r` runs, and at most `L^r` compressed run-length tuples. Hence

`|A_(k,R) intersect [1,X]|`

` <= 2^(k(k-1)/2) sum_{r=0}^R (2^k L)^r`

` = O_(k,R)((log X)^R)`.                                  (5)

Allowing `k` to range over any fixed finite bound only changes the constant. Therefore every bounded-feedback, bounded-run grammar has natural density zero among the positive integers, and also among the positive odd integers. It cannot realize all positive odds except a finite forbidden set.

This result is not a claim that fixed feedback number alone is inadequate. **Unbounded run complexity is deliberately left open.** In particular the four-letter `k=2` automaton is still a meaningful candidate for an unbounded word construction.

### Multiplicative closure does not repair the prime deficit

Let `A` be any finite union of classes covered by (5), with a fixed upper bound on `k,R`. Close `A` under finite products and adjoin the value one (order joins and source/sink padding).

If a prime `p` belongs to this multiplicative closure, one of its factors must equal `p`, so `p` was already in `A`. But `A` has at most a polylogarithmic number of values up to `X`, whereas the number of primes up to `X` is at least a constant multiple of `X/log X`. An elementary proof of the needed lower bound suffices: the central binomial coefficient is at least `4^n/(2n+1)`, and its prime-power contribution for any prime `p<=2n` is at most `2n` (each floor-difference in its p-adic valuation is zero or one). Hence

`4^n/(2n+1) <= binom(2n,n) <= (2n)^pi(2n)`.

Taking logarithms gives the asserted lower bound. Thus the multiplicative closure also misses infinitely many odd primes, and cannot discharge the prime lane in THM-4094.

**Interpretation.** The earlier finite-bank obstruction extends to an infinite bank of arbitrarily large tournaments whenever both feedback size and run-template complexity are bounded. Merely inserting longer homogeneous transitive blocks cannot supply the missing prime atoms. A viable completeness construction must escape at least one of those bounds or use a genuinely different operation.

## 5. Exact probe boundary

Two small algebraic probes were run on the two-feedback, four-letter representation. They test scalar affine append identities, not tournament coverage:

- Prefix and suffix words of length at most three, every pair of distinct letters, and middle words through length three: no candidate with distinct offsets for `H(PwaS)=2H(PwS)+c_a` and `H(PwbS)=2H(PwS)+c_b` survives.
- Prefix/suffix length at most two, append blocks of length one or two, slopes `2,3,4,5,9`, and middle block words through depth two: no distinct-offset candidate survives.

Because every exact identity would have to survive these tests, the negative results rule out only these specified short-context schemes. They do not prove that longer blocks, higher-dimensional invariants, other slopes, or a non-affine inverse compiler fail. Scripts and raw outputs are preserved separately as `probe_affine.*` and `probe_blocks.*`.

## 6. The intercyclic lane fits k=2, but bounded runs are not automatic

A primary-source theorem strengthens the naive feedback bound from three to two:

Bang-Jensen and Bessy, [*Small cycle transversals in tournaments*](https://www.lirmm.fr/~bessy/publis/SmallTransversTournament.pdf), preprint dated 16 January 2013, Theorem 1.8, proves that every tournament with no two vertex-disjoint directed cycles has a feedback vertex set of size at most two, and that this is sharp. In tournaments, the condition is equivalent to having no two vertex-disjoint directed triangles.

Thus the whole triangle-packing-one lane is contained in the four-letter gap automaton, not merely the eight-letter `k=3` automaton. This is a **CITED** structural reduction. The paper's theorem gives no bounded-run condition, and none is inferred here. Arbitrarily alternating one-apex words already warn against such an inference for an arbitrary chosen feedback presentation. No completeness or new missing value for the full intercyclic class is asserted by this addendum.
