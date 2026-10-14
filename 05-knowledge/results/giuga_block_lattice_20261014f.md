# Signed prime completions, Giuga index jumps, and a genuine mod-eight lattice

**PROVED scoped:** an exact finite two-prime completion compiler, index-jump
obstructions for proper Giuga multiples, and the star-matrix identities below.
**CITED:** the Giuga/Carmichael characterization, classical neighbour
constructions, and van der Blij's theorem; the recent star-lattice congruence
is credited to Mavrotheris. **FINITE-EXACT:** the declared compiler and matrix
controls. Giuga's conjecture, universal Collatz, rational distances to all four
square corners, plane colouring and planar Jacobian remain **OPEN** here.
No historical novelty claim is made for Egyptian-fraction factorization.

Artifacts: [checker](../../04-computation/experiments/giuga_block_lattice_20261014f.py),
[JSON](giuga_block_lattice_20261014f.json),
[one-prime/CRT interfaces](giuga_korselt_defects_20261014f.md),
[Collatz source guards](collatz_successor_guards_20261014f.md), and
[finite witnesses and continuum boundary](continuum_finite_witness_boundary_20261014f.md).

## 1. Inheritance, portfolio and exact meanings

The [September Giuga circuit audit](collatz_mod6_20260921_grand_circuit_typing.md)
already rederived the classical equivalence: a composite counterexample to
Giuga's primality conjecture is simultaneously a Giuga and a Carmichael number.
The [successor-block packet](carmichael_successor_blocks_20261014b.md) already
gave complete one/two-prime Carmichael extensions of 561, 1105 and 1729.
The canonical opposed witnesses are 30 (Giuga, not Carmichael) and 561
(Carmichael, not Giuga). There is no known distinguished integer that is both,
and the conjecture does not assert uniqueness of such an integer.

Our anchor is complete arithmetic completion of a retained core. The niche is
the integral quadratic form behind signed reciprocal-prime sums. The wildcard
is Guy D19 and continuum-dependent constructions. The live board is **core /
sign / integer index / prime-address guards / determinant / finite witness**.
The least-used coordinate is the integer index, which changes under adjoining
factors even when a modular condition looks unchanged.

The relevant [META-PATTERNS](../../00-navigation/META-PATTERNS.md) are
"Search the statement before the method", "Type every analogy and every
implication", and "Audit saturation and basis covariance before naming a
lattice bridge". We apply them by retaining the factor labels and proving the
determinant; an unlabelled isometry class is not claimed to recover the integer.

For squarefree `N>1` define

\[
S(N)=\sum_{p\mid N}\frac1p,\qquad
\kappa_\mu(N)=S(N)+\frac\mu N,\quad \mu\in\{-1,1\}.
\]

The negative sign and a composite integer with integral index give a Giuga
number. The positive sign with index one gives a primary pseudoperfect number
(PPN); other integral positive-sign indices are a wider class. Squarefreeness
and primality of the retained labels are essential.

In the companion's cofactor convention `p|(N/p-alpha)`, `alpha=-mu`.
The index-one neighbouring constructions `N(N-1)` and `N(N+1)` are classical;
[Sondow--Tsukerman, Theorem 8](https://arxiv.org/abs/1401.0322) is the attribution.
The companion gives the complete one-prime consumer for arbitrary squarefree
cores, with the exceptional dual-sign source 2 retained.

## 2. Every two-prime completion reduces to finite divisor lists

Fix any squarefree core `N>1`. Seek distinct new primes `q<r`, neither dividing
`N`, for which `Nqr` has a chosen signed integral index. No bound on these new
primes will be assumed. Let `q0<r0` be the two smallest primes not dividing `N`.
For any completion with index `k`,

\[
S(N)<k\le S(N)+\frac1{q_0}+\frac1{r_0}
                  +\frac\mu{Nq_0r_0}.                         \tag{1}
\]

The lower inequality follows from `N(q+r)+mu>0`. For the upper one, the
function `1/q+1/r+mu/(Nqr)` decreases in each variable: after clearing positive
denominators the relevant coefficients are `Nr+mu>0` and `Nq+mu>0`.
Thus (1) gives a finite, computable list of integer indices.
In fact it gives **at most one**: its width is at most
`1/2+1/3+1/12=11/12<1`, using `N>=2` and `q0>=2,r0>=3`.

For each such `k`, put

\[
A=N(k-S(N))=kN-\sum_{p\mid N}N/p>0.
\]

Clearing the target reciprocal equation gives

\[
Aqr=N(q+r)+\mu,
\qquad
\boxed{(Aq-N)(Ar-N)=N^2+\mu A.}                 \tag{2}
\]

Both factors on the left are positive: for example
`Aq-N=(Nq+mu)/r>0`. They are strictly ordered because `A>0` and `q<r`.
Consequently the following is a complete terminating classifier:

1. Enumerate the indices in (1), and form `D=N^2+mu A`.
2. Enumerate positive divisor pairs `u<v`, `uv=D`.
3. Require `A|(u+N)` and `A|(v+N)`, then recover
   `q=(u+N)/A`, `r=(v+N)/A`.
4. Authenticate primality, distinctness and coprimality with the core.

Every completion appears. Conversely each accepted pair satisfies the exact
reciprocal identity and its local prime conditions. Factoring by trial division
is a terminating algorithm, with no efficiency or polynomial-time assertion.
This generalizes the familiar index-one factorization into an explicitly
bounded consumer for any retained core and allowed index.

Some complete catalog entries are:

| Core | Target | All pairs of new primes | Resulting integers |
|---:|---|---|---|
|2|Giuga|`(3,5)`|30|
|2|positive index one|`(3,7)`|42|
|6|Giuga|`(7,41)`, `(11,13)`|1722, 858|
|6|positive index one|`(7,43)`|1806|
|561|Giuga|`(2,59)`|66198|
|1806|Giuga|`(3041,4447)`|24423128562|

The index interval contains no integer for many empty rows; others are rejected
by the finite factor lists and native prime guards. These are complete
two-prime catalogs for the named cores, not catalogs of all signed integers.

## 3. The two ways to complete 561 pay different obligations

The inherited unique two-prime **Carmichael** extension is

\[
561\cdot401\cdot3041=684106401.
\]

The unique two-prime **Giuga** extension is instead

\[
561\cdot2\cdot59=66198.
\]

The first fails the Giuga cofactor condition; the second fails Korselt and is
even. Thus neither is a counterexample to Giuga's conjecture. This supplies an
actual common source and two different completion maps, rather than identifying
their outputs because both are prime-factor repairs.

The [CRT companion](giuga_korselt_defects_20261014f.md) combines the old-prime
obligations into one residue class when `gcd(N,lambda(N))=1`. The two-prime
compiler additionally discharges every new prime's cofactor equation. The
Korselt periods of those new primes are still a separate final check.

## 4. A Giuga multiple must pay an integer index jump

Suppose `N` and `NQ` are Giuga, where `Q>1` is squarefree and coprime to `N`.
Writing its new prime factors as `q_1<...<q_t`, the exact change is

\[
m=\kappa_{-1}(NQ)-\kappa_{-1}(N)
 =\sum_{i=1}^t\frac1{q_i}+\frac{Q-1}{NQ}.        \tag{3}
\]

It is a strictly positive integer. In particular, a fixed-index Giuga family
is an antichain under divisibility. Extending an index-one example to another
Giuga integer would require a higher-index target; it cannot simply reproduce
the same index-one construction.

Let `ell_1<...<ell_t` be the smallest primes not dividing the retained core.
Then

\[
1\le m<\frac1N+\sum_{i=1}^t\frac1{\ell_i}.       \tag{4}
\]

Whenever the right side is at most one, **every** block with that many new
primes is excluded, regardless of their sizes. Exact rational sums give:

| Giuga core | Minimum number of new primes in any proper Giuga multiple |
|---:|---:|
|30|56|
|858|55|
|1722|55|
|66198|54|
|2214408306|54|

The first permitted count is a necessary lower bound, not an existence claim.
For the first four rows the reciprocal threshold is first crossed at prime277;
the retained support explains the differing counts. This is the classical
reciprocal-capacity mechanism specialized to a source, not a claimed improvement
over the far larger published bounds for Giuga--Carmichael intersections.

This is a decisive limit to the analogy with the two-prime repair of561:
adjoining two primes cannot repair an already Giuga core into a larger Giuga
integer. Signed transitions and changing the retained core are different moves.

## 5. A literal graph-shaped lattice and its eightfold constraint

There is an actual arithmetic-to-quadratic-form map here. For retained distinct
primes `p_1,...,p_r` with integral signed index `k`, form the symmetric matrix

\[
B=\begin{pmatrix}
k&1&\cdots&1\\
1&p_1&&0\\
\vdots&&\ddots&\\
1&0&&p_r
\end{pmatrix}.
\]

It is the Gram matrix of an integral lattice on a weighted star. Its Schur
complement is `k-S(N)=mu/N`. Hence

\[
\det B=\mu,\qquad \operatorname{signature}(B)=r+\mu.           \tag{5}
\]

The prime labels and distinguished central vertex recover `N` and the signed
equation. Forgetting them and retaining only an isometry class is not claimed
to be lossless. An arbitrary rational reciprocal approximation does not produce
an integral unimodular lattice: integrality and determinant are the gates.

When `N` is odd, the central basis vector `e_0` is characteristic, because
`B e_0` agrees with the diagonal of `B` modulo two. Its norm is `k`.
[Van der Blij's theorem](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/vdb.pdf)
therefore yields

\[
\boxed{r+\mu\equiv k\pmod8.}                                  \tag{6}
\]

This star-lattice application is explicitly credited to
[Mavrotheris, *Characteristic-Vector Congruences for Weak Primary Pseudoperfect
and Giuga numbers* (2026 preprint)](https://zenodo.org/records/22085830).
The elementary matrix and characteristic-vector steps are displayed above;
the non-elementary input used is precisely van der Blij's congruence.

For an odd Giuga number, `r=k+1 mod8`. If its index is one, its number of prime
factors is therefore `2 mod8`; index one must not be assumed for an arbitrary
unknown Giuga number. For two odd Giuga integers related by `N -> NQ`, (6)
and (3) further imply `t=m mod8`. Equation (4) and this congruence can be
applied together to a proposed block. Neither gives a bound on every possible
number of new primes.

For the known even example30, the star has determinant `-1`, signature2,
and characteristic vector `(0,1,1,1)` of norm10: the difference is exactly8.
By contrast, using prime leaves3,5 and an arbitrary central coefficient1
gives determinant7, signature3, and characteristic norm1. The congruence fails
because the lattice is not unimodular. A graph shape or dimension alone cannot
replace this hypothesis.

The relevant eightfold mechanism is now an explicit signature theorem. Bott
periodicity belongs to a different specified functor; the occurrence of eight
does not supply a map from this lattice into a Collatz terminal certificate.

## 6. What transfers to the other proof targets

**A genuine, scoped index-to-Collatz implication.** Suppose a hypothetical
composite `N` satisfies both Giuga and Korselt, and its Giuga index `k` is odd.
Korselt makes `N` odd. Its number `r` of prime factors is even: reducing
`kN=sum(N/p)-1` modulo two already gives `k=r-1 mod2`. This is the parity
shadow of (6); the full signature theorem is not needed for this step.
If any prime factor is `1 mod4`, Korselt forces `N=1 mod4`. Otherwise all prime
factors are `3 mod4`, and their even count again gives `N=1 mod4`. Thus

\[
U(N)=\frac{3N+1}{2^{v_2(3N+1)}}<N.
\]

In particular this applies to index one, without asserting that every unknown
Giuga number has index one. A joint source outside this immediate-descent rule
must have **even index, an odd number of factors, and all prime factors
`3 mod4`**. This pays one actual step, conditional on the stated arithmetic
source; it proves neither the nonexistence of that source nor its ROOT suffix.

The [source-guard companion](collatz_successor_guards_20261014f.md) extracts a
literal descending Collatz prefix from the complete Korselt conductor. It also
shows that partial old-prime conditions coexist with every finite legal valuation
word along an infinite prime progression. Better prime-distribution estimates
could bound prime supply there; they would not authenticate the deliberately
unpaid new-prime condition. The full distinction survives assuming GRH.

Guy D19 is a rational-distance section, not a name for Giuga's conjecture.
The [finite-witness companion](continuum_finite_witness_boundary_20261014f.md)
retains its four square-corner equations and shows why algebraic realization
does not supply rational distances. The same note distinguishes CH independence,
forcing-invariance of finite witnesses, and the possible proof strength of
additional axioms. No continuum axiom is presented as a Collatz, HN or planar-JC
proof. Rational/irrational and algebraic/transcendental are different partitions;
transcendental real numbers are still real numbers.

Incoming work changed the choice of mechanism: the repaired collector note's
ensemble law cannot replace a source trajectory (MISTAKE-607); the HN experiments
needed a rigid piece in addition to a forced radius. Two stale unconditional
cycle-entry sentences in that collector note were qualified in this checkpoint
to match its existing correction. The pre-audit Fermat/circulant note is treated
as an idea source, not a proved dependency.

The next precise arithmetic task is an adaptive completion search that retains
the old-prime CRT class, signed integer index, and newly introduced prime guards
simultaneously. One- and two-prime completion is now finite and exact. Letting the count
grow still needs a coverage or stopping argument; no unbounded closure follows
from the compiler.

## 7. Reproduction and audit scope

```text
python -B 04-computation/experiments/giuga_block_lattice_20261014f.py
python -B -O 04-computation/experiments/giuga_block_lattice_20261014f.py
```

Both runs give **7,970 checks**, identical LF JSON, SHA256
`0f4fd5177bd192408a12351a3e09a67489e7820d37bf37a9d61112953fdf7d70`.
The independent local-prime oracle uses every squarefree core2..150 and every
pair of distinct new primes at most31. All returned completions, including
larger primes, are independently replayed. The eleven named cores use the
proved unbounded-prime classifier. Matrix controls use exact rational
determinants, characteristic vectors, Schur complements, signed examples and
the nonunimodular hostile. No optional assertion disables a check.

The two-prime completeness proof, index-jump bound, star determinant and odd
signature implication were independently audited by the Collatz lane, with an
optimized replay matching the saved output. The companion packets received
independent cross-audits; their exact scopes and reproduction commands are
retained in the linked notes.
