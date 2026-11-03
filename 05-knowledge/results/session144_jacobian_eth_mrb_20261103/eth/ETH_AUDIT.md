# ETH import: scope, exact interfaces, and reduction budgets

Session 144, 2026-11-03. Repository baseline: `origin/main` at `5e1cedf`.

**Status:** primary-source statement audit and independent finite-exact interface
checks. The large imported algorithm is **CITED / EXTERNAL**, not newly proved
or built here. The reduction-budget and correction-sidecar lemmas below have
elementary proofs independent of the imported algorithm. Their algorithmic
applications are **CONDITIONAL** on the cited runtime theorem.

## 1. What changed since the prior import

The inherited entry is
`05-knowledge/results/external_preprints_mub6_local_smoothing_eth_merge_20261027.md`,
section 3. Its explicit claim that nothing had been re-derived remains accurate
for that session. Its October 27 search finding is not a current review status.

The organizers' [October 26 statement](https://ethdiscussion.blog/2026/10/26/eth-refutation-examination-and-discussion/)
distinguishes organizing an examination from certifying the result. Their
[October 30 replication thread](https://ethdiscussion.blog/2026/10/30/lean-replications/)
contains independent reports of successful Lean checks. [Carlos Toledo's own
audit](https://carlostoledo.co/eth/) reports both Lean and nanoda acceptance and
separately identifies the interpretation of machine models and input encodings.
These are public primary reports, not executions performed by this session.
We should neither erase that positive evidence nor relabel it as our own full
verification. The review forum still calls the result a claimed refutation.

The current upstream commit observed through the GitHub connector is
`4eaee87548250b2682affda7d05d947c05554bae`. `SOURCES.json` records pinned URLs,
exact byte lengths, SHA-256 hashes and Git blob hashes for every source read.
No source is inferred from an unverified secondary news headline.

**Source-materialization correction.** The first temporary manifest hashed
local copies containing one extra trailing LF introduced by this session's
patch-writing wrapper. At 01:54:41 UTC that single added byte was removed from
each of the eight copies. This was not an upstream revision. `SOURCES.json`
retains the initial local lengths/hashes alongside the corrected source hashes,
explains the exact writer error and correction, and distinguishes independently
verified source bytes from temporary materializations.
The independent session audit subsequently re-downloaded all eight files from
their commit-pinned raw URLs and confirmed byte-for-byte agreement, lengths,
SHA-256 hashes and Git blob hashes. The manifest identifies the verification
script and receipt by name and SHA-256.

## 2. Three bounds that must not be merged

The [October 25 manuscript, Theorem 1.1](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/main_full_version.tex)
claims, for each fixed width `k>=2`, deterministic exact counting with

1. general input: `poly_k(n+m) 2^{O_k(n loglog n / sqrt(log n))}`;
2. fixed density `m<=Dn`: `poly_{k,D}(n+m) 2^{O_{k,D}(n/sqrt(log n))}`;
3. its earlier rational interpolation route: the weaker exponent
   `rho(n)=n(loglog n)^2/sqrt(log n)`.

The [formalization scope note](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/docs/378.md)
explicitly identifies its theorem with bound 3, expressed using encoding length
`L`, and polynomial visited tape space. It does not claim to formally certify
the two sharper exponents. All three exponents are `o(n)` as `n` grows with
`k` fixed. Thus the weaker bound already has the headline ETH consequence;
the gap matters for quantitative transfers and practical crossover claims.

Here `n` is the number of declared Boolean variables, `m` the number of clauses,
and `L` the length of the actual encoded input. They are not interchangeable.
For unnormalized inputs with repeated clauses, time must include reading the
input. After removing repetitions, fixed-width CNFs have polynomially many
possible distinct clauses. A pure `2^{o(n)}` statement is therefore consistent
with the displayed polynomial input factor, but that factor must not disappear
when discussing enormous explicitly repeated inputs.

## 3. Definition and quantifier audit

Read the [self-contained Comparator target](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/ComparatorChallenges/SubexponentialSAT.lean)
before reading the proof implementation. Its order is

`for every fixed k>=2, there exists one finite machine and constants,
 for every n and every width-at-most-k CNF, there exists a bounded correct run`.

The machine is chosen **before** `n` and the formula. This is uniform in input
size, not a family of unrelated per-size advice machines. It does not supply
a useful bound uniform in growing `k`, nor an efficient compiler whose cost is
controlled as a function of `k`. The usual ETH/SETH consequences need neither.

Its primitive machine has finitely many states and tapes. A step reads one
symbol from one selected tape and either halts or performs one tape move/write.
There is no free large-integer operation, unbounded alphabet, semantic SAT
oracle, or density promise in this target. The real constants appear in the
bound, not as machine instructions.

The [solution endpoint](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/OAI/Algorithms/SubexponentialSAT/SqrtLogSpace/Final.lean)
contains `exact_count_and_sat` and `exact_count_eth`. It returns the exact
assignment count as a binary word of length at most `n+1`, the SAT bit, and a
polynomial visited-space bound for every prefix of the run. This is stronger
than a nondeterministic existence assertion about a satisfying assignment.

The Comparator challenge's terminal `sorry` is its requested statement, not
the implementation proof. The corresponding JSON points at solution module
`OAI.Algorithms.SubexponentialSAT.Main` and permits only `propext`,
`Quot.sound`, and `Classical.choice`; `enable_nanoda` is false as shipped.
Independent reports enabling another kernel are separate evidence. This
session did not compile either the solution or its full dependencies.

### Unary encoding is polynomial here

The [encoding source](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/lean/OAI/Algorithms/SubexponentialSAT/Encoding.lean)
uses unary numeric fields, but literal codes are below `2n`, not exponential
integers. Its proved `formulaBitLength_upper` states

`L <= n+2+[4+(8n+8)k]m`.

One public replication report left this cost unresolved. It is explicitly
bounded in the source. The endpoint's `exact_count_eth` transports the encoding
factor to a polynomial in `n+m+1`. Canonicalization just sorts and deduplicates
literal and clause lists; it preserves SAT and the number of assignments,
including unused variables. It does not preserve the multiplicity-weighted
Max-SAT objective, which requires its own input handling.

Independently unfolding the frames gives the exact length identity

`L = n+2+4m + sum_{all literal codes c}(4c+8)`.

Since `c<2n`, this also proves a slightly sharper upper bound
`L <= n+2+[4+(8n+4)k]m`. This elementary strengthening concerns encoding only.
The script checks both the exact length and the published weaker inequality
on all 1,434 CNFs in the declared finite universe.

## 4. Consequences with their directions exposed

Let `A_k` denote the asserted deterministic exact counter for width `k`.
The following implications are logical deductions from that assertion:

- `A_3` decides SAT by testing whether the count is positive, so refutes ETH.
  It counts exactly, so also refutes deterministic `#ETH`.
- `A_k` for **every fixed** `k` refutes usual SETH: fix, for example, an
  exponent improvement of `1/2`; each width eventually admits a better bound.
  Constants and thresholds may depend on the width.
- A deterministic solver is also a randomized solver that uses no coins.
  The corresponding bounded-error randomized hypotheses therefore fail.
- Testing whether the count is zero decides UNSAT deterministically, hence
  also nondeterministically. This refutes usual NSETH.
- Exact 3-SAT decision solves the satisfiable-versus-constant-gap promise, so
  the deterministic and randomized Gap-ETH versions fail.
- It does **not** give polynomial-time SAT, settle `P=NP`, give a practical
  SAT implementation, or prove subexponential time in the original input
  length for every problem in NP.
- Fixed-width guarantees do not, just by changing the symbol `k`, cover
  unrestricted-width CNF-SAT. Nor does refuting SETH refute OVC: the known
  implication `SETH => OVC` cannot be reversed or contraposed incorrectly.

Conditional lower bounds are not invalid theorems. If a paper proved
`ETH => lower_bound`, negating ETH removes the assumption's force but does not
negate `lower_bound`. The Boolean valuation
`SAT_algorithm=true, ETH=false, separate_target_lower_bound=true` satisfies both
implications. The lower bound here concerns a separate target problem, not
the SAT lower bound directly contradicted by the asserted algorithm.
The script enumerates that exact logical hostile. To obtain an algorithm for
a target one needs a reduction **from the target to** the fast SAT problem,
with a useful size bound, not merely an old SAT-hardness reduction to the target.

## 5. A quantitative reduction-budget lemma

**PROVED conditional composition lemma.** Fix a width and suppose a solver has
time `poly(L) 2^{C f_p(N)}`, where

`f_p(N)=N(loglog N)^p/sqrt(log N)` for sufficiently large `N`.

Suppose an input of original scale `b` is transformed in `2^{o(b)}` time into
at most `2^{o(b)}` formulas, all with `L<=b^{O(1)}` and `N<=g(b)`, where `g`
is polynomially bounded. If `f_p(g(b))=o(b)`, applying the solver and a
`2^{o(b)}` aggregation procedure costs `2^{o(b)}`. Count preservation must be
proved separately if the desired answer is a count.

**Proof.** For large arguments `f_p` is increasing. The logarithm of the total
cost is bounded by the reduction/aggregation cost in the exponent, the log
number of calls, `O(log b)`, and `C f_p(g(b))`. Every term is `o(b)`.
Finitely many small instances can be handled separately. This proves a
sufficient transfer criterion, not a lower bound when it fails.

For an explicit size growth
`g(b)=Theta(b (log b)^a/(loglog b)^c)`, with `a>=0`, the transported exponential
factor has exponent divided by `b` of order

`(log b)^(a-1/2) (loglog b)^(p-c)`.

Thus it tends to zero for `a<1/2`, or for `a=1/2` and `c>p`. At the equality
boundary it is merely exponential in `b`; beyond it this upper-bound transfer
does not beat that scale. For the formal bound `p=2`, a useful sufficient
budget is

`N=o(b sqrt(log b)/(loglog b)^2)`.

For the sharp **sparse** paper bound, use `p=0`, provided the reduction also
guarantees a fixed clause/variable density. Ordinary gate-by-gate Tseitin
encodings of a bounded-fan-in Boolean circuit have this property and unique
gate extensions. For a witness with `b` bits and a verifier with `S(b)` gates,
the relevant variable cost is `N=O(b+S(b))`, not just `b`.

This repairs the inherited reflection's overly narrow suggestion that only
linear encodings can benefit: mildly superlinear encodings can too. Conversely
`N=Theta(b^2)` transports to an exponent of order
`b^2(loglog b)^p/sqrt(log b)`, not `o(b)`. A generic polynomial-time reduction
is not enough. This statement compares guarantees, not actual runtimes.

The prior position-by-vertex Hamiltonian-path encoding uses quadratic size.
It therefore supplies no improvement over the repository's `2^v poly(v)`
subset DP by this transfer. That is a diagnosis of that encoding, not a proof
that no better reduction or algorithm exists. Similar bookkeeping is required
before speeding a bounded Jacobian coefficient/collision search. A faster
finite search does not provide an exhaustive classification or a bound beyond
which witnesses cannot occur.

## 6. Exact correction sidecars: a genuine mechanism, not an analogy

The [October 30 exposition by Ryan O'Donnell](https://www.cs.cmu.edu/~odonnell/some_eth_false_exposition.pdf)
extracts a low-degree matrix family with easy block structure at selected
points plus one common low-rank correction. The correction is evaluated,
not thrown away. Hermite data then recover the permanent polynomial. This
suggests a concrete lesson for classification work: compressed coordinates
can be useful only when their lost information is reconstructed at affordable
cost. It is not evidence that the same construction applies to Jacobian maps.

Here is a completely explicit audit of the identity. On the eight vertices
of the three-cube let `X_1,X_2,X_3` flip coordinates and let

`S=X_1+X_2+X_3`, `J=X_1 X_2 X_3`, `K=X_1+X_2/2+X_3/3`.

The polynomial `p(z)=-z` agrees with parity on the two central Hamming levels.
Its error matrix `E=J+S` has rank two: in the common Walsh basis it vanishes
on weights one and two and is nonzero on weights zero and three. Define

`A(T)=J+TK`, `B(T)=-S+TK`, so `A(T)=B(T)+E`.

At `T=i`, `B(i)` omits `X_i`, hence splits into the two blocks with coordinate
`i` fixed. This gives the source's exact structure-plus-correction pattern.
Nevertheless `per(A(0))=1`, while `per(B(0))=81`. Small rank is not small error
in an exact count. The script verifies all matrix identities and computes
the permanent polynomial by subset DP and independently at five points by
enumerating all `8!` permutations.

### Why the number of jets matters

`P(T)=per(A(T))` has degree at most eight. Three Taylor coefficients at each
of `T=1,2,3` recover it exactly by a rational nine-by-nine linear solve.
If only two coefficients at `T=3` are retained, the nonzero polynomial

`V(T)=(T-1)^3(T-2)^3(T-3)^2`

is invisible to every retained jet but changes `P(0)` by `V(0)=72`. Thus even
lossless local computations do not imply global recovery unless their total
multiplicity exceeds the degree. This is an explicit hostile to that overclaim,
not a counterexample to the paper, whose interpolation budget includes the
needed derivatives.

## 7. A bounded-integer Keller-collision interface

The following is an elementary, deliberately restricted bridge proposed in
this session. It is not a reduction of all complex Jacobian counterexamples
to a finite census.

**PROVED interface.** Fix `r>=1`, `D>=1`, and positive bit bounds `B,H`.
Use canonical signed binary encodings for the coefficients of a polynomial
map `F:Z^r -> Z^r` of total degree at most `D`, and for two integer points
`u,v`. There is an effectively constructible bounded-fan-in Boolean circuit
that accepts exactly the labelled triples satisfying

`det JF is a nonzero constant`, `F(u)=F(v)`, and `u != v`.

For fixed `r,D`, the circuit has size polynomial in `B+H`. Its gate-by-gate
3-CNF encoding has `O(W+G)` variables and clauses and a unique extension for
each input tuple, where `G` is the circuit's gate count and

`W = r binom(r+D,D) B + 2rH`

is the number of free input bits. Consequently exact SAT counting counts
labelled ordered triples, not equivalence classes of polynomial maps.

**Proof and size certificate.** Set `L=r(D-1)`. Every variable in `det JF`
has degree at most `L`. Evaluate the determinant exactly at every point of
`{0,...,L}^r`, and require every value to equal its nonzero value at the origin.
A polynomial with degree at most `L` in each variable that vanishes on this
grid is zero: apply the univariate root bound one coordinate at a time.
Thus these finitely many equalities are equivalent to the polynomial identity,
over characteristic zero. Also evaluate `F` at both unknown integer points,
compare the outputs and test inequality of the input vectors.

These are integer computations with explicitly bounded intermediate bit lengths.
Write `M=binom(r+D,D)` and `T=(L+1)^r`. At a collision point, an output has
`O(B+DH+log(M+1))` bits. At a grid point, a derivative entry has
`O(B+log(DM+1)+(D-1)log(max(1,L)))` bits; determinants add a factor of `r`
and `log(r!)` to that bound. Expanded determinant evaluation uses at most
`r!` products per grid point. Schoolbook Boolean arithmetic therefore gives
an explicit polynomial circuit for fixed `r,D`, with no real-number oracle
or unit-cost multiplication assumption. If `r,D` vary, the monomial count,
grid size `T`, determinant expansion, and arithmetic widths must remain in
the cost ledger; there is no uniform polynomial-size claim in them.

For Tseitin conversion use equivalences between each gate variable and its
gate operation, and require the output gate to be true. Each input assignment
then has precisely one gate extension; ordinary one-way implications would
not suffice to preserve counts. Two's-complement coefficient/point encodings
avoid multiple encodings of zero.

**CONDITIONAL algorithmic conclusion.** With `N=O(W+G)`, the imported general
formal bound gives `poly(N)2^{O(N(loglog N)^2/sqrt(log N))}` for this bounded
search. Because the circuit encoding is sparse, the sharper paper claim
would give `poly(N)2^{O(N/sqrt(log N))}`. Both are subexponential in the
encoded variable scale `N`. They are subexponential in the original free-bit
scale `W` only when the appropriate section 5 budget is met. Merely knowing
`G=poly(B+H)` does not establish that.

**Boundary.** The search is already finite and decidable without ETH. Its
acceleration gives no bounds on degree, coefficient height, or collision
height. Maps with nonintegral complex coefficients or without a bounded
integer collision witness are outside this universe. Even a complete tuple
count weights a map by its ordered collision pairs. No conjugacy, stable
equivalence, minimal dimension, or classification normal form is inferred.
The finite test checks that the `r=2,L=2` full interpolation matrix has rank
nine, whereas the undersized grid `{0,1}^2` has rank four and misses `x(x-1)`.

## 8. Reproduction and exact limits

Run from this directory:

```sh
python3 verify_eth_interfaces.py > verification.json
python3 -O verify_eth_interfaces.py > verification.optimized.json
cmp verification.json verification.optimized.json
```

The script uses only the Python standard library and exact rational arithmetic.
Checks are not disabled by `-O`. Receipt universe: every distinct-clause CNF
with `n=0..3`, `k=0..3`, and at most two clauses; one fully specified cube
replacement family; five independent permanent evaluations; exact Hermite
reconstruction and its missing-jet hostile; Boolean implication enumeration.
The Keller interface also has the complete small interpolation-grid test
described above; no large Boolean solver was constructed.

**Not done:** full Lean build or axiom audit, implementation/execution of the
ETH counting algorithm, independent proof of the tensor-contraction runtime,
an unbounded-width uniform algorithm, a faster repository Hamiltonian counter,
or a Jacobian classification. Source text lives only in the temporary audit
directory; publish this report, script, receipts and source manifest, not the
downloaded source tree.
