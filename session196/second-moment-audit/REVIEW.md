# Independent audit: S196 variance expansion and planted-spine construction

Date: 2026-11-05. Scope: `session196/sequence-structure/PROOF.md`,
`verify.py`, its saved receipt, and the named inherited sources.

## Outcome

**PASS on the intended all-order mathematics and implementation.** The
independent exact checks and fresh derivation below confirm the second-moment
formula, orbit weights, all-order asymptotic argument, and deterministic
construction. The parent/producer were notified immediately of a significant
inheritance omission, which was repaired in the proof and README. The exact
pair formula, even-run proof, ordinary generating series, and leading `2/n`
rate are inherited results, not new S196 results.

One quantifier precision edit was requested and applied: the truncation order
`J` in Lemma 3 and Theorem 4 is now explicitly a nonnegative **integer**.
This is the intended use in the proof and verifier. A real-valued interpretation of `J>=0` would
not justify the displayed geometric-tail exponent `J+1`.

No dirty checkout was edited. All audit artifacts and producer replays are
in this separate directory. The snapshots are identified by
`SOURCE_SHA256SUMS`; exact verifier outputs were unchanged by the prose fixes.

## Independent derivation through Parseval

This gives an independent route to each positive summand, rather than simply
repeating the packet's fixed-spine inclusion-exclusion.

Let the fair independent orientation variables be `X_ij` for `i<j`, extended
antisymmetrically. The Hamiltonian count has the multilinear expansion

    H = 2^(1-n) sum_pi product_(consecutive a,b in pi) (1+X_ab).

A selected edge support can occur in this expansion only if its nontrivial
components are disjoint paths. For a path component of odd edge length,
reversing the ordering changes the sign and cancels its coefficient. For a
support with `2j` edges and `r` positive even-length path components, its
Walsh coefficient has absolute value

    2^(1-n) 2^r (n-2j)!.

There are exactly

    N(n,j,r) = n! binom(j-1,r-1)
               / [(n-2j-r)! 2^r r!]

such supports: choose an ordered composition of `j` into `r` half-lengths,
assign an ordered list of distinct vertices to these paths, then divide by
the two orientations of each path and the ordering of components. Unused
vertices are isolated. The allowable range is exactly
`1<=r<=min(j,n-2j)`.

Orthogonality of distinct Walsh characters gives

    E[H^2] = (n!/2^(n-1))^2
             + sum_(j,r) N(n,j,r)
                          [2^r (n-2j)! / 2^(n-1)]^2.

Using `(n-2j)!/[(n-2j-r)!r!] = binom(n-2j,r)` produces precisely

    E[H^2] = n! W_n / 4^(n-1)

with the packet's positive double sum. This argument also explains why no
negative variance terms arise, even though individual inverse-power
asymptotic coefficients may be negative. It agrees with the inherited
Walsh-energy formulation in THM-848; no priority claim attaches to it.

## Conventions and boundary cases

- `H` counts directed Hamiltonian vertex **words**. A path and its reversed
  word are not identified. A reversed word is not a directed path in the
  same nontrivial tournament. Ordered pairs are not unordered pairs.
- The probabilistic law is uniform over **labelled** tournaments, with each
  orientation an independent fair bit. Uniform isomorphism classes are a
  different law.
- For `n=1`, `W_1=1`, `H=1`, both moments equal 1, and both orbit masses equal 1.
- For `n=2`, `W_2=2`, both labelled tournaments have `H=1`, and both moments and
  orbit masses again equal 1. In the second orbit identity the factor is
  `2^(-1)W_2=1`, so negative exponents introduce no exception.
- `W_0=1` is a formal-series convention. The moment theorem is correctly
  restricted to `n>=1`; substituting `n=0` in `n!/2^(n-1)` would be wrong.
- The automorphism action on words is free, so `|Aut(T)|` divides `H(T)`.
  Its diagonal action on ordered pairs is also free. The orbit identities
  require `H/|Aut|` and `H^2/|Aut|`, not `H` or `(H/|Aut|)^2`.
- Explicitly, for the two order-three classes, `(H,|Aut|)=(1,1),(3,3)`.
  Their weighted first mass is 2 and second mass is 4. This is a useful
  hostile case for accidentally using uniform unlabelled probabilities.

## Tail and asymptotic proof audit

Set `m=n-2j`. Cancel the last `r` denominator factors against `binom(m,r)`:

    binom(m,r)/(n)_(2j)
      <= 1/[r! (n)_(2j-r)] <= 1/[r! (n)_j].

The first inequality follows since each numerator factor in
`m(m-1)...(m-r+1)` is no larger than the corresponding positive factor in
`(m+1)...(m+r)`. The second uses `r<=j`. Dropping `1/r!` and enlarging the
range to all `1<=r<=j` gives the layer bound

    sum_r t_(j,r) <= 2*3^(j-1)/(n)_j <= (2/3)(6/n)^j.

The final step holds uniformly over all allowed `j`, since
`j<=floor((n-1)/2)`. For integer `J`, summing the geometric series proves the
stated tail. For each fixed integer `J`, the remaining finite rational
functions have Taylor expansions at `1/n=0`, giving the stated
`O(n^(-J-1))` error. The implicit constant may depend on `J`.

The displayed coefficients through `n^-6`, and the receipt's `n^-7`
coefficient, are correct. Independent polynomial long division also gives

    [x^0,...,x^10] R_(1/x) =
    0, 2, 0, -14/3, -58/3, -916/15, -6986/45,
    -13010/63, 344948/315, 1079114/81, 189551482/2025.

These are additional finite symbolic controls, not an assertion that the
full asymptotic series converges. The `48/n^2` elementary error bound follows
from `24/[n^2(1-6/n)]<=48/n^2` for `n>=12`.

## Conditional-expectation audit

Conditioning on a specified Hamiltonian word gives the exact mean
`W_n/2^(n-1)`. Equivalently, choosing a uniform planted word and then a
compatible random tournament size-biases the original law by `H`, hence the
mean is `E[H^2]/E[H]`. Permutation invariance makes every fixed planted word
have that same mean.

Every required edge in a Hamiltonian word concerns a distinct unordered
pair, so the product of fixed/unfixed edge probabilities is exactly its
conditional probability. The weighted endpoint recurrence is valid, and
the two possible orientations of any free pair average to the current
mean. A deterministic larger-mean choice cannot decrease it. This proves
the all-order guarantee, including `n=1,2` where there are no decisions.

There are `O(n^2)` decisions, each requiring two `O(n^2 2^n)` subset DPs,
with `O(n 2^n)` storage. The stated arithmetic-operation qualification is
correct. The 0/1/2 implementation consistently scales every mean by
`2^(n-1)`, regardless of how many edges have already been fixed.

The construction proves a large-count witness at each vertex order. It does
not prescribe an arbitrary odd `H` and does not settle the `{7,21}` support
problem or a best-known maximum-H asymptotic.

## Exact independent controls

`independent_audit.py` imports none of the producer code. It checks:

1. Every labelled tournament on `n=1,...,6`, counting Hamiltonian paths by
   accumulating every word into all compatible completions, rather than
   endpoint DP. Both moments agree exactly.
2. Actual relabelling orbits and free-action divisibility through `n=6`.
   Class totals are `1,1,2,4,12,56`; both weighted masses agree.
3. Every binary selected-spine-edge string through `n=18`, rejecting any
   odd-length run, independently equals the double sum.
4. An odd-block recurrence independently evaluates the formal composition
   coefficients through `n=80`.
5. Rational polynomial long division obtains the asymptotic coefficients
   through `n^-10`.
6. All truncations `J` for every `n=12,...,120`, totaling 3,624 exact rational
   tail inequalities.
7. Every recorded conditional-expectation decision through `n=10`, using a
   separate recursive memoized weighted count. For `n<=6`, each branch mean
   is additionally recomputed as the literal average of the earlier
   exhaustive tournament census restricted to compatible completions.

Independent witness counts are `1,1,3,5,15,43,175,653,3003,15281`.

Both `python3 independent_audit.py` and `python3 -O independent_audit.py`
passed with byte-identical output. Both normal and optimized executions of
the frozen **producer** verifier also passed with byte-identical output and
receipts. The producer's replay receipt matches its original supplied receipt.
No removable `assert` statements are used in the independent checker.

## Precise inherited overlap and external source scope

The following repository sources were inspected directly:

- `01-canon/theorems/THM-589-cvH-even-run-closed-form-and-2-over-n-rate.md`:
  exact same signed-adjacency permutation sum, positive even-run formula,
  second moment, and leading `2/n` concentration rate.
- `01-canon/theorems/THM-589-metagraph-H-variance-closed-form-nud-succession.md`:
  same moment formula, weighted odd compositions, and kernel
  `z(1+z^2)/(1-z^2)`. Its 2026-10-14 correction explicitly retracts transfer
  of the unweighted two-term recurrence to weighted `W_n`.
- `05-knowledge/results/oeis_transfer_atlas_20261014h.md`:
  inherited run identity plus a general compiler and explicit seven-step
  polynomial-coefficient recurrence for `W_n`.
- `01-canon/theorems/THM-848-hamiltonian-path-walsh-needle-stalk.md`:
  exact even-path-forest Walsh support/coefficients and level-energy EGF.
- `01-canon/MISTAKES.md`, MISTAKE-151 concerning THM-204/THM-848:
  the tempting monomial level-energy extrapolation fails first at `n=8`.
  The current packet avoids this error and gives `R_8=131/560`.

The two THM-589 filenames matter because their numeric theorem IDs collide.
The revised packet now cites both exact slugs and states the inheritance.
The identified incremental comparison is the rigorous uniform tail bound,
arbitrary-order expansion, and planted-spine algorithm. This is not a
repository-wide or literature-wide novelty determination.

Primary external sources inspected on 2026-11-05:

- Noga Alon, *The Maximum Number of Hamiltonian Paths in Tournaments*,
  <https://web.math.princeton.edu/~nalon/PDFS/hamilton.pdf>, introduction:
  confirms the classical first-moment lower bound and its attribution to
  Szele.
- N.C. Wormald, *Tournaments with many Hamilton cycles*,
  <https://arxiv.org/abs/2610.16245>, abstract: states a substantially stronger
  asymptotic lower bound and says the same statements hold for paths.
- Alexander Divoux and Shouda Wang, *The maximum number of Hamiltonian cycles
  in tournaments*, <https://arxiv.org/abs/2610.23755>, abstract: states an
  asymptotic extremal result with an implication for Hamiltonian paths.

The latter two proofs were **not** independently audited and are not used
as mathematical dependencies. They warrant the explicit scope statement
that this packet gives a simple deterministic guarantee, not a best-known
extremal bound. No literature-priority claim is made.
