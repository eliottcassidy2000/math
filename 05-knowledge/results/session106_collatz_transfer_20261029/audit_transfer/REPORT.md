# Independent adversarial proof audit: transfer, spectral boundary, observability

2026-10-29. **PASS within stated scope. No substantive proof defect found.**
No producer files changed. This review does not prove Collatz convergence, supply a
source-tail estimate, verify literature priority, or audit the separate descent,
reverse, and history packages.

## Findings first

1. The six-coordinate source-interval carrier is valid for every finite positive
   valuation word and odd integer multiplier b>=3. It is a growing exact carrier,
   not a finite-state convergence certificate.
2. Ordinary native killed ℓ¹ operator norms remain one at every depth even if
   Collatz terminates. Full-support atomic extinction is equivalent to **strong**
   stability instead. The notes keep these quantifiers separate.
3. The explicit inverse chain supplies all exact finite stopping clocks and hence
   unbounded rank of all-source survival responses. This is not a complexity or
   nonrationality theorem for an aggregated scalar series.
4. Independent checks pass 22,134 explicit conditions, identically with normal
   Python and `-O`. They include exhaustive *excluded* residue tests, b=9 controls,
   genuinely substochastic randomized kernels, and native rank witnesses over
   GF(2) through dimension 100. No producer modules are imported or rerun.

## Exact cylinders and interval filtering

The guard must be `bt+1+2bAq = 2^a mod 2^(a+1)`, not mere divisibility by
2^a. Dividing by two yields a congruence modulo 2^a with odd invertible
coefficient bA. Thus precisely one s in [0,2^a) is valid; the stated child
coordinates are integral and the child endpoint is odd. Its least source obeys
1<=r+Ms<M2^a. No hidden floor, translation, or rational-source inclusion enters.
Induction from (1,2,1,1) proves both directions of the exact-cylinder statement:
every selected source realizes every earlier guard, and every source realizing
the word lies in that cylinder. The final affine residue formula is also valid;
it is not being used to replace the intermediate guards without proof.

For every interval [L,U], congruence preimage is exactly the stated integer
ceil/floor interval in q'. Negative numerators in the ceiling are permitted.
Survival requires comparison with the immutable original source. Its coefficient
is 2A'-M', not an endpoint-to-endpoint drift. Dividing by a negative coefficient
reverses the inequality and gives the stated upper floor; coefficient zero has
the correct all-or-nothing rule. Old prefix constraints are retained by the
interval preimage, so successive intersections prove all-prefix closure.

For finite N, retained source sets are disjoint and nonempty, bounding cell count
by the number of initial odd sources. On a finite cell endpoints increase with q,
so `floor(log2(b*y_max+1))` bounds all next valuations, including the equality
case where the numerator itself is a power of two. The producer's bit-length
range includes that case. Finite-depth termination requires no conjecture.

The generating mass formula has exponent equal to **source**, not endpoint.
Countable unfiltered child sums equal their parent coefficientwise because each
positive integer has a unique finite valuation; there is no missing infinite-
valuation integer. Formal sums are locally finite coefficientwise. For real
0<z<1 all terms are nonnegative, so regrouping is valid. Survival indicators
are decreasing in depth. Their full-support weighted limit is exactly the mass
of sources never descending, by dominated/monotone continuity. If any such source
exists its unchanged atom prevents decay. Conversely, pointwise extinction and
summability imply decay. Strong induction then makes all-source first descent
equivalent to termination for b=3. No finite census proves this limit.

Hostile checks: the b=5 native cycle 13,33,83 has exact guards (1,1,5), and its
minimum immutable source 13 survives forever, so its positive atom blocks an
all-b contraction. The rational fixed points 19/101 and 31/101 reclose after
appending 1 to 37/35 and 221/175, respectively. This is carry-sensitive algebraic
reclosure, not native advancement of an integer periodic orbit. The note makes
that distinction explicitly.

## Spectral and Abel claims

The column killed operator acts on ℓ¹(Y), not on an un-killed stationary space.
At depth k, Q^k applied to a source atom has mass one exactly when τ>k. The
initial terminal mass has τ=0 and is absent already at k=0. The identity for A_k
therefore has no off-by-one or missing-terminal defect.

Full-support A_k→0 implies extinction of each point mass, hence termination of
all sources; extinction of all point masses implies strong stability by finite-
support approximation and the contraction norm. The converse from strong
stability is immediate. Signed input vectors introduce no issue: approximate in
ℓ¹ and use the norm bound, or dominate by their absolute values.

For every k>=1 the changing source 2^(k+1)-1 has all first k valuations exactly
one and endpoints >1, including the last endpoint 2*3^k-1. Thus each Q^k has
norm one. The spectral-radius formula applies to this bounded operator and yields
radius one; k=0 is the identity. There is no assertion of a single infinite orbit.

For positive substochastic Q,R, L>=(1-alpha)Q gives L^k>=(1-alpha)^k Q^k by
positivity and induction, without commutation. The transfer bound requires a
strictly smaller rate than 1-alpha, or its explicitly stated ratio limit. The
one-state fixed-point/refresh-to-terminal example proves the boundary cannot
be weakened to equality. A stationary mixing gap on another finite carrier is
not this killed-survival estimate.

The Neumann series exists for every z<1 regardless of termination. Nonnegative
Tonelli interchange and the finite/infinite geometric sums give
D_z=sum w(x)(1-z^τ(x)); initial terminal sources contribute zero. Dominated
convergence at z↑1 gives precisely the nonterminating mass, not finite mean
lifetime. Bounded resolvent norm would instead require finite weighted mean.
The one-atom strict tests follow because that atom's contribution is a lower
bound on D_z. A strict inequality `w(x)(1-z^m)>d` excludes τ>=m, including
infinity. Equality is insufficient and the note does not use it.

A nonnegative summable supersolution dominates every finite partial Neumann sum:
iterate the inequality, retaining the nonnegative remainder. Monotone passage
gives the claimed mass upper bound. Constructing such a family with improved
mass remains the arithmetic obligation, not a result of the certificate format.
Finite-head bounds correctly charge unresolved sources their full atoms and
retain the exact telescoping tail. They cannot newly certify an unresolved
source by throwing away its own mass.

The lazy-clock factorization and normalized Abel identity hold algebraically.
For state-dependent strictly positive acceptance, each finite native edge has
almost surely finite waiting time; a countable intersection gives preservation
of the infinite path, but no uniform physical-time tail follows. Zero acceptance
can trap a source. The aggregation counterexample correctly distinguishes native
survival mass 1/2 from resampled survival 2^-k.

## Inverse chain and rank

For odd y not divisible by three, the proposed exponent pairs ensure integrality
and oddness. If both candidate predecessors were divisible by three, both
numerators would vanish modulo nine; multiplication by four makes this
impossible. At least one candidate is therefore admissible. All candidates with
a>=2 exceed y except y=1,a=2, and the explicit first step excludes that exception.
Because y is odd, `3x+1=2^a*y` proves the *exact* valuation. The factor 32/3 bounds
each step, so the displayed global size bound follows by induction.

Strictly increasing constructed sources cannot hit 1 prematurely when followed
back along this inverse chain. Thus τ(x_j)=j, not just τ(x_j)<=j. This is a
constructive terminating family, independent of the Collatz conjecture.

Rows i=0,...,N-1 and columns j=1,...,N give entries 1[j>i]. Under these indices
the matrix is upper triangular with diagonal one, hence determinant one over
every field. Factorization through d-dimensional source states forces d>=N,
even with time-dependent linear output functionals. Nonzero source-column
weights preserve rank over their coefficient field. The statement does not
prohibit unbounded integer coordinates, nonlinear encodings, scalar averaging,
or a concise formula for a terminating rank. The deterministic finite-state
variant likewise needs distinct states for distinct complete stopping responses.

## Independent receipt scope and reproduction

`check_independent.py` uses only standard-library arithmetic and explicit
exceptions, with no `assert` and no producer imports. It checks:

* 117 words of lengths 1..3 on letters 1..3 at b=3,5,9, enumerating every odd
  residue in each full cylinder modulus (both acceptance and rejection);
* 468 all-prefix interval cuts in final source coordinates, 585 independently
  enumerated next guards, and 1,755 interval pullbacks;
* literal rational generating masses at nonempty prefixes and both hostiles;
* all 36 two-state half-integral substochastic kernels and every pair thereof
  for noncommutative mixed-kernel domination; analytical 2x2 resolvents and lazy
  identities, plus every two-state deterministic map with initial terminal mass;
* 100 inverse levels constructed by exhaustive small-exponent search rather
  than residue cases, with native replay, size bounds, and the actual response
  matrix's rank over GF(2); 64 changing-source norm-one controls.

```sh
python3 session106_collatz_transfer/audit_transfer/check_independent.py > session106_collatz_transfer/audit_transfer/checks.normal.json
python3 -O session106_collatz_transfer/audit_transfer/check_independent.py > session106_collatz_transfer/audit_transfer/checks.optimized.json
cmp session106_collatz_transfer/audit_transfer/checks.normal.json session106_collatz_transfer/audit_transfer/checks.optimized.json
```

Both executions passed; receipts match byte-for-byte. Finite checks are not
extrapolated to universal results: those rest on the proofs reviewed above.

### Reviewed input fingerprints

SHA-256, relative to `session106_collatz_transfer/`:

```
6b4b583f52675a20fe98d5610713d83ece131999826c660e014541baa97975b9 recurrence/SOURCE_INTERVAL_TRANSFER.md
8f3501aef71107b43d870c54336a8c82052e921baa991d84f974c7325fd0164e recurrence/check_recurrence.py
5582f87245fea9c62547e3bb8320bbce8141f1a0ec67b8742db08502f919f153 spectral/NATIVE_SPECTRAL_BOUNDARY.md
42010a4f53c5c3316ae77e29d7ecd1f7e35a2bce1769be14e621c1caec51d21f spectral/check_spectral.py
cf276a23f387802c74dc9c6234a962400fcebcf2ace9b7522ad95710709a4e72 observability/SOURCE_RESPONSE_RANK.md
6f331b02987bce3a7b8ca77602dc6f106e07e7585ede5c274bff93f0cef8a52f observability/check_observability.py
```

No blocked verification in this assigned scope. Remote publication is left to
the parent session; this local workspace itself is not a Git checkout.
