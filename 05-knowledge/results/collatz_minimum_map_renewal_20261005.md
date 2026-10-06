# The minimum map: renewal structure of hard excursions along running minima, coverage decay as a killed renewal, and five structural rhymes

2026-10-05, opus session `opus-2026-10-05-S10` (minimum-map renewal). Owner's
seed: "work the renewal structure of hard excursions along running minima;
merge in the ideas of five papers and their connections with unexplored past
work in this repo; be loose and free in your creativity." The papers:
Hoelscher, semicomplete sequences and the Pell constant (arXiv:2102.07083);
Liu, Wang, Zhang, Hankel transforms and (α,β) Somos-4 sequences
(arXiv:2608.08703); Fox and Schildkraut, nearly consecutive sequences without
long arithmetic progressions (arXiv:2608.20533); Riordan and Scott, a short
proof of the Erdős–Sós conjecture (arXiv:2609.15893); Tóth, ceiling Mills
functions (arXiv:1801.08014).

**Status.** PROVED (elementary, with Terras's cylinder structure as the only
imported input): Theorem RN (renewal of first-descent blocks), Theorem CD
(every excess-rate family has natural density zero, with an explicit rigorous
rate), Proposition KR (the killed-renewal exponent). FINITE-EXACT: the block
laws and Hankel probe. VERIFIED: the census below `2^22` (renewal tests,
exponents, coalescence spectrum, renewal density). HEURISTIC and typed: the
five rhymes of section 6; nothing there is a theorem about Collatz. **OPEN:**
universal entry; the asymptotic coalescence spectrum. Not a canon promotion.

## 0. Answers

- **The object.** Let `m(x)` be the first value below `x` on the odd orbit of
  `x` (the endpoint of its first-descent segment). The chain of running minima
  of `n` is the `m`-orbit of `n`; it ends at 1. The functional graph of `m` is
  a tree with ROOT 1, internal nodes the units, leaves the odd multiples of
  three (never images). The excess-rate families `F_q` are the sources whose
  whole `m`-orbit uses `q`-admissible blocks.
- **Renewal (Theorem RN).** Under the uniform law on odd residues modulo
  `2^K`, the types `(l, A)` of the successive first-descent blocks are
  independent and identically distributed with the cylinder law
  `p(l,A) = N(l,A)/2^A` of the expense note, as long as they fit in `K` bits;
  the actual running minima coincide with the coefficient descents except at
  the least representatives, at most one per word. Census below `2^22`: the
  laws of blocks 1 to 5 of sources in `[2^21, 2^22)` are within total variation
  `0.001` to `0.013` of the cylinder law, and the law of block 2 given the type
  of block 1 is within `0.002` to `0.009` of the unconditional law.
- **Scaling.** `E[l] = 3.4926`, `E[A] = 2E[l]` (Wald), `E[excess] = (2 -
  log_2 3) E[l] = 1.4496` bits per block. Sources in `[2^21, 2^22)` have on
  average `14.62` running minima (renewal `14.83`) and `53.3` odd steps
  (renewal `51.8`). The renewal density is flat in log scale: `0.665` to
  `0.696` visits per dyadic block per source, against `1/E[excess] = 0.690`.
- **Coverage decays as a power (Theorem CD, Proposition KR).** The proportion
  of `F_q` in `[X, 2X)` decays like `X^(-c_q)` where `c_q` is the root of
  `sum_(good types) p(l,A) 2^(c excess) = 1`. Predicted and fitted on
  `2^10 .. 2^21`: `q=3`: `0.2275` vs `0.2250`; `q=8`: `0.1559` vs `0.1577`;
  `q=16`: `0.0751` vs `0.0727`; `q=32`: `0.0527` vs `0.0507`; `q=100`:
  `0.0066` vs `0.0062`. The naive exponent `-log_2(1-p_q)/E[excess]`
  overshoots (`0.298` at `q=3`) because hard blocks have small excess, so
  survivors need fewer blocks. Every bounded-rate family, every two-tier
  family with a fixed threshold, and every finite-bank family therefore has
  natural density zero. This closes the question left by the expense note.
- **Coalescence.** Chains coalesce into the `m`-tree. The visit probability
  `pi(x)` (share of sources whose chain passes through `x`) is `0.938` at 5
  (the trunk-entry share of THM-4517, 0.93796, is the same event up to the
  source 3), `0.434` at 23, `0.254` at 11, `0.224` at 61, `0.191` at 7,
  `0.179` at 19, `0.176` at 37, then the three hard excursions `47, 31, 91`
  at `0.095, 0.093, 0.088`. The ten most visited minima carry 21% of all chain
  visits, the top thousand 55%, the top ten thousand 73%: a positive-density
  renewal in log scale concentrated on a sparse skeleton of attractive minima.
- **Rhymes.** The first-descent words form a Kraft-complete prefix code
  (Hoelscher's completeness, with the never-descended mass as the deficit);
  the block-length law has no finite Hankel order, is not a positive-measure
  moment sequence (`H_4 < 0`) and has no Somos-4 structure, while the
  integer-slope control has `H_k = 2^(-k(2k-1))` exactly, a degenerate
  Somos-4 (the Liu–Wang–Zhang class begins at slope 2 and is already lost at
  the convergent 3/2); the extremal rising words are Sturmian words of
  slope `log_2 3` (Fox–Schildkraut's small-step, structure-avoiding shape);
  Erdős–Sós is the contrast case where density forces a tree, while here
  density never forces an atom; Tóth's floor-to-ceiling variation is the
  upper/lower semiconvergent duality of the expense note. Section 6 types
  each with map, preserved predicate, lost data, and a cheap probe.

## 1. Inheritance and board

Closest proved mechanisms: Theorem T of
[the expense note](collatz_expense_diophantine_20261005.md) (the cylinder law
of first-descent types), Theorem D and Proposition R of
[the weak-reset family](collatz_weak_reset_family_20261005.md), Terras's
theorem that the first `K` parity bits are a bijection with residues modulo
`2^K` (CITED; the shadow theorem's forward half), THM-4512 (at most one
exceptional member per contracting cylinder), THM-4517 (trunk-entry
spectrum, mac-mini), and the thread's earlier renewal results (S22, Lemma R'
and Theorem C on the Fourier side). Canonical hostile: the excursion `47 -> 23`,
a running minimum of 9.5% of all sources with rate 31. Corrected near miss:
the naive survival exponent `-log_2(1-p_q)/E[excess]`, which ignores that
killed blocks are the short-drop ones. Least-used sidecar: the renewal density
in log scale, which is flat while the visit measure is wildly non-uniform.
Board: **minimum map / block type / renewal / killed renewal / visit measure /
coalescence skeleton / Kraft completeness.**

**Repository sweep (one Explore agent, 70 tool calls, six themes).** The
exact hit is THM-4514 Theorem 3 (proved, audited): a tight spine block of
length `L >= 2` exists iff `L = 1 + floor(b beta)` with
`beta = log_2 3/(log_2 3 - 1)`, and such a block has `ceil(L log_3 2)` odd
steps. Since `1/alpha + 1/beta = 1` for `alpha = log_2 3`, the Beatty
sequences `floor(b alpha)` and `floor(b beta)` are complementary, and the
primitive cone depths of Lemma T4b (`{l alpha} < alpha - 1`, equivalently
`floor(l alpha) = floor((l-1) alpha) + 2`) correspond exactly to the tight
spine lengths through `L = floor((l-1) alpha) + 2 = floor(l alpha)`, the
total valuation of the primitive word. Lemma T4b is THM-4514 Theorem 3 in
backward coordinates; its sufficiency direction is therefore proved twice,
by the explicit word and by that theorem. The other findings are typed in
section 6.

## 2. Theorem RN: renewal of first-descent blocks

Let `b = (b_1, b_2, ...)` be the parity vector of an odd source under the
shortcut map `T` (`b_1 = 1`). The coefficient first-descent time is the least
`j` with `3^(k_j) < 2^j`, `k_j` the number of ones among the first `j` bits;
the block is the prefix of length `j`, with type `(l, A) = (k_j, j)`. The next
block is read from the remaining bits.

**Theorem RN (PROVED).** Under the uniform (Bernoulli 1/2) law on parity
vectors, the successive block types are independent and identically
distributed with law `p(l,A) = N(l,A)/2^A`. For odd integers below
`X = 2^K`, uniformly chosen, the first `k` blocks have exactly this joint law on
the event that their total length `sum (A_i + 1)` is at most `K`, and the
actual running minima agree with the coefficient descents except when a
running minimum is the least representative of its cylinder, which happens for
at most `X^(h*+o(1))` sources below `X`.

*Proof.* The coefficient first-descent time is a stopping time of the i.i.d.
bit sequence, so the block decomposition is a renewal process and the
post-block bits are again i.i.d. (strong Markov property of the shift). Terras:
the first `K` bits of the parity vector are determined by, and determine,
`n mod 2^K`, so among the `2^(K-1)` odd residues the bits are uniform. Actual
descent along a block word is equivalent to coefficient descent unless the
start is at most the positive rational cycle point of the word (Proposition
R), and by THM-4512 at most one member of each cylinder below `X` is in that
position; the number of first-descent words with `A < K` is at most the
rising-prefix count, which is `2^(h* K + o(K))`. QED.

Census (sources in `[2^21, 2^22)`, 1,048,576 of them): total variation to the
cylinder law `0.0010, 0.0030, 0.0044, 0.0072, 0.0129, 0.0259` for blocks 1 to
6 (the drift is the finite size: block `k` starts near `2^(21 - 1.45 k)` where
the cylinders of long words are no longer fully represented); block 2 given
block 1 of type `(1,2), (1,3), (2,4), (5,8)`: `0.0022, 0.0034, 0.0048, 0.0088`.
Mean lengths of blocks 1 to 6: `3.49, 3.48, 3.48, 3.48, 3.48, 3.49` against
`E[l] = 3.4926`.

## 3. Theorem CD and Proposition KR: coverage decays as a killed renewal

**Theorem CD (PROVED).** For every `q`, the natural density of `F_q` is zero;
more precisely, for every integer `k`,

    dens(F_q ∩ [1, 2^K)) <= (1 - p_q)^k + P(A_1 + ... + A_k + k > K) + 2^(-(1-h*)K + o(K)),

where `p_q` is the exact tail of Theorem T. Choosing `k = K/(E[A] + 1 + eps)`
gives a rate `2^(-c' K)` with `c' = -log_2(1 - p_q)/(E[A] + 1 + eps)` (about
`0.054` at `q = 3`). The same holds for the two-tier families with a fixed
threshold and for the families relative to a finite bank, since a fixed table
only shortens the chain by a bounded number of blocks.

*Proof.* A member of `F_q` has every block admissible, in particular its first
`k` coefficient blocks when they fit in `K` bits and are actual descents; by
Theorem RN those blocks are i.i.d. with the cylinder law, so the probability
that all are good is `(1 - p_q)^k`; the remaining sources either have blocks
that do not fit, or are exceptional members. Weak law of large numbers for the
renewal gives the second term for `k` below `K/(E[A]+1)`. QED.

**Proposition KR (killed renewal; PROVED as a renewal statement, heuristic as
applied beyond the first `K` bits).** Let `S(L)` be the probability that a
renewal chain with i.i.d. block types, killed at the first block of rate
above `q`, survives until the accumulated excess reaches `L`. Then
`S(L) = 2^(-c_q L + o(L))` with `c_q` the unique root of

    sum_(types with rate <= q) p(l,A) 2^(c (A - l log_2 3)) = 1.

The chain from a source `n` needs excess `log_2 n - O(1)` to reach the bottom,
so the renewal model predicts `dens(F_q ∩ [X,2X)) ~ X^(-c_q)`.

*Proof of the renewal statement.* `S` satisfies `S(L) = sum_(good t) p(t)
S(L - excess(t))` for `L > 0`; this is a defective renewal equation whose
exponential rate is the Cramér root. QED.

| `q` | `p_q` | naive `-log_2(1-p_q)/E[excess]` | killed-renewal `c_q` | fitted on `2^10..2^21` |
|---:|---:|---:|---:|---:|
| 3 | 0.25903 | 0.2984 | 0.2275 | 0.2250 |
| 8 | 0.17631 | 0.1930 | 0.1559 | 0.1577 |
| 16 | 0.08008 | 0.0831 | 0.0751 | 0.0727 |
| 32 | 0.05545 | 0.0568 | 0.0527 | 0.0507 |
| 100 | 0.00668 | 0.0067 | 0.0066 | 0.0062 |

Coverage of `F_3` per dyadic block: `0.205` at `2^10`, `0.129` at `2^13`,
`0.081` at `2^16`, `0.050` at `2^19`. The agreement of the killed-renewal
exponents with the fits to within `0.003` at every rate is the quantitative
content of "renewal structure of hard excursions": the hazard per block is
the Diophantine tail of Theorem T, and survivors are selected for large drops.

Consequence for the family program: the run-block tree (inside `F_3`), the
excess-rate families, the two-tier families and the bank-relative families
are all density-zero sets, with explicit exponents. Coverage proportions
reported below `2^20` are finite-size values on their way to zero.

## 4. Coalescence: the visit measure of the minimum map

Define `pi(x)` as the proportion of sources below `X` whose chain passes
through `x`. Below `2^22`:

| `x` | `pi(x)` | own block | own rate |
|---:|---:|---:|---:|
| 5 | 0.9382 | (1,4) | 1 |
| 23 | 0.4335 | (3,7) | 2 |
| 11 | 0.2537 | (3,6) | 3 |
| 61 | 0.2241 | (1,3) | 1 |
| 7 | 0.1914 | (4,7) | 7 |
| 19 | 0.1792 | (2,4) | 3 |
| 37 | 0.1757 | (1,4) | 1 |
| 47 | 0.0948 | (34,55) | 31 |
| 31 | 0.0928 | (35,56) | 67 |
| 91 | 0.0878 | (28,45) | 46 |
| 59 | 0.0760 | (4,8) | 3 |
| 49 | 0.0703 | (1,2) | 3 |
| 43 | 0.0690 | (3,5) | 13 |
| 65 | 0.0682 | (1,2) | 3 |
| 29 | 0.0591 | (1,3) | 1 |

`pi(5) = 0.9382` is the trunk-entry share through 5 of THM-4517 (`0.93796`
to `2^30`): "5 is a running minimum" and "the orbit enters the trunk through
5" differ only for the source 3. Shares of all chain visits: top 10 minima
21%, top 100 38%, top 1,000 55%, top 10,000 73%. Among `x` in `[2^8, 2^16)`,
63% are visited at least 100 times and 0.6% at least 10,000 times. In log
scale the renewal measure is flat (`0.69` visits per dyadic block per source),
so the concentration is entirely in which integers at a given scale are
attractive: `pi(x)` divided by the uniform expectation `2/(E[excess] ln 2 x)`
is `2.3` at 5, `5.0` at 23, `6.8` at 61, `2.3` at 47, and `0.66` at 7. The
in-degree of `m` among sources below `2^22` has mean `1.5` on units (exactly
two thirds of the sources are units and every source has one image), 36% of
units receive no first descent from below `2^22`, and the maximum is 34 at
`x = 17711`, the Fibonacci number `F_22`. A fair test kills the curiosity:
within factor-two windows the other unit Fibonacci numbers have in-degree
ranks `0.36, 0.11, 0.46, 0.23, 0.65, 0.71, 0.30, 0.25`; only 17711 is
exceptional (rank 1.000), a single outlier, not a Fibonacci effect.

The asymptotic visit measure `pi` on the integers is the open object of this
note: it is the `m`-subtree mass, satisfies `pi(x) = sum_(m(x')=x) pi(x')`
plus the starting mass, and its support on the popular hard excursions is
exactly why inheritance dominates the chain rate.

## 5. What this settles

- Lemma T4b's Beatty criterion is THM-4514 Theorem 3 in backward coordinates
  (section 1), which closes that obligation a second time.
- The renewal structure is proved and quantitatively exact: i.i.d. blocks
  with the Theorem T law, flat renewal density `1/E[excess]`, killed-renewal
  coverage exponents matching the census to `0.003`.
- Every proved family of the program has density zero (Theorem CD). The
  program's honest asymptotic is the exponent `c_q`, not a coverage
  proportion. Universal coverage needs an argument of a different kind.
- The next admissible target is the visit measure `pi`: its concentration on
  a sparse skeleton is what makes inheritance dominate, and the skeleton is
  determined by the local minima of popular excursions (47 and 91 inside 27's
  excursion). A proof that `pi` has a limit, or an exact formula for `pi` at
  small `x` in terms of basin densities (THM-4517's objects), is the open
  question.

## 6. Five structural rhymes, typed

Each row: the paper's mechanism, the Collatz object it rhymes with, what the
map preserves and loses, and the cheap probe run here. These are analogies of
shape, as the owner's standing instruction asks; none is a theorem about
Collatz unless marked.

| Paper | Mechanism there | Rhyme here | Preserved / lost | Probe and verdict |
|---|---|---|---|---|
| Hoelscher, semicomplete sequences and the Pell constant | a sequence is complete up to a threshold built from fractional sums; only three semicomplete arithmetic sequences; a generating-function identity for the density of solvable negative Pell equations (odd-period continued fractions) | the first-descent words form a prefix code with Kraft sum `sum N(l,A)/2^A = 1`; truncating at length `L` leaves the deficit `P(sigma > L)`: `0.148, 0.065, 0.019, 0.0030, 0.00062` at `L = 5, 10, 20, 40, 60` (exact law), the "semicompleteness threshold" of the code | preserved: completeness as a Kraft identity and the Diophantine density flavour (period parity there, upper versus lower semiconvergents here); lost: no Pell structure, no hypercube count | PROVED as a Kraft identity (Theorem RN); the Pell side is a rhyme of two continued-fraction densities, nothing more |
| Liu, Wang, Zhang, Hankel transforms and Somos-4 | Hankel determinants of a generating function obey `S_n S_(n-4) = a S_(n-1) S_(n-3) + b S_(n-2)^2` when the function has a J-fraction of elliptic type; Sulanke–Xin transform; periodic Hankel determinants | the block-length law `c_l = 1/2, 1/8, 1/8, 3/64, 7/128, 3/128, 15/1024, 85/4096, ...` and its Hankel transform `H_k = 1/2, 3/64, 1/1024, -171/2^24, ...` | preserved: the question "is the law's generating function of finite Hankel order"; lost: algebraicity, since the lattice paths run under a line of irrational slope `log_2 3` | FINITE-EXACT: no `H_k` vanishes through `k = 10`, `H_4 < 0` (not a positive-measure moment sequence, no positive J-fraction), and the `(a,b)` fitted from `H_5, H_6` fails beyond; the integer-slope control (threshold `2j`) gives `H_k = 2^(-k(2k-1))` exactly, i.e. Somos-4 with `(alpha,beta) = (2^-12, 0)`, and the convergents `3/2, 8/5, 19/12, 65/41` already show sign changes: the elliptic/quadratic class lives at slope 2 and nowhere on the way to `log_2 3` |
| Fox and Schildkraut, nearly consecutive AP-free sequences | steps 1 or 2, length `2^k/k^2`, no `k`-term arithmetic progression; the construction balances local freedom against global structure | rising words have letters mostly 1 and 2 and must stay under the line `A_j < j log_2 3`; the extremal rising word of each admissible length is the lower mechanical (Sturmian) word of slope `log_2 3` (Lemma T4b of the measurement-independence note), aperiodic with minimal complexity | preserved: small steps plus a global constraint give exponential but entropy-reduced families (`2^(h* A)`, `h* = 0.95`); lost: AP-avoidance is not the Collatz constraint, and Beatty sequences do contain long progressions | HEURISTIC rhyme only; recorded so that nobody reads Sturmian aperiodicity as AP-freeness |
| Riordan and Scott, Erdős–Sós | average degree above `k - 1` forces every `k`-edge tree; extremal graphs determined | the minimum-map tree has infinite in-degree at every unit in the limit (every unit has infinitely many predecessors with a smaller first descent), so every finite tree embeds trivially; the repository's standing obstruction is the opposite direction: density of the backward tree never forces an atom | preserved: "density forces structure" as a slogan; lost: the Collatz statement is universal over integers, not existential over copies | contrast, not transfer; the in-degree spectrum of `m` (mean 1.5 on units, 36% in-degree zero below `2^22`, maximum 34) is the finite shadow of the infinite branching |
| Tóth, ceiling Mills functions | `ceil(B^(c^n))` prime for all `n`, dual to Mills' `floor(A^(3^n))`; nested intervals built from short-interval primes encode infinitely many choices in one constant | the upper and lower semiconvergents of `log_2 3`: descents from above use `A_min(l) = ceil(l log_2 3)` (Theorem Q), rising cones from below use `floor(l log_2 3)` (Lemma T4b); `E_inf` (always-rising parity vectors) is a nested-cylinder construction of 2-adic points, like Mills' constant | preserved: floor/ceiling duality and the nested-interval existence of a point with an infinite property; lost: for Collatz the question is whether such a point is a positive integer, the analogue of asking whether Mills' constant is rational | analogy; it names the two sides of the Diophantine structure in one picture |

**Unexplored repository work surfaced by the sweep (pointers, typed).**

- Pell. The 2026-08-23 reflection on cube packets and Pell returns
  (`07-reflections/incoming-agent-signals-cube-packet-pell-return-and-nodal-closures-20260823.md`,
  monad-claudebox, not mac-mini) and THM-3819 (Pell numbers `h^2 - 2g^2 = -1`,
  return times `7, 41, 239, 1393`) are first-return structures with a
  deterministic return, the contrast case to the independent renewals of
  Theorem RN; the Pell constant, Stevenhagen's conjecture and period parity
  have zero hits, and a `0.580574` in an LRC output is not the Pell constant
  `0.5805776` (numerology trap). Two next steps of that reflection were never
  taken: a general modular orbit sieve beyond THM-3858's mod 9, and the owner-
  labelled Pell prefix packets. Brown's completeness criterion appears only in
  two drafts.
- Hankel. THM-3288 (LRC, Hankel rank 15 against denominator degree 14);
  Collatz: Proposition 2 of the depth-layers note of 2026-09-27 (eventual
  periodicity iff the Hankel determinants of `2^(d_k(n))` vanish, Kronecker;
  audit owed), the theta lane's 2-adic Hankel irrationality (Bézivin), and the
  unproved `±1` shifted Hankel pattern of the parity skeleton `sum w^(2^j)`
  (checked to `n = 100`, no literature search). The present probe sharpens to:
  `H_4 < 0`, so the block-length law is not the moment sequence of a positive
  measure and admits no positive J-fraction. The rational-slope control of
  section 7 shows where the Liu–Wang–Zhang class begins: at integer slope 2
  the law is Catalan-like and `H_k = 2^(-k(2k-1))` exactly (hexagonal
  exponents), a Somos-4 relation with `(alpha, beta) = (2^-12, 0)`; at every
  convergent of `log_2 3` below 2 (`3/2, 8/5, 19/12, 65/41`) and at the limit
  the determinants change sign and no Somos-4 relation fits.
- Beatty and AP structure. THM-4514 (above), THM-4505's window theorem along
  the Beatty word of `log_2 3`, the Lonely Runner AP axis (LEM-010/012/013,
  THM-1120, THM-1171, THM-1158) and the Christoffel work (THM-536, THM-778 to
  THM-853, THM-4422). Alon–Zaks and nearly consecutive sequences have zero
  hits; nothing measures AP content inside a valuation word, and the obvious
  measurement (longest progression among partial sums) is dominated by runs of
  ones, hence by the minimal-sum record words, so it was not pursued.
- Trees. Erdős–Sós and Sumner's conjecture are greenfield (`PROBLEM-LEDGER.md`
  marks Sumner and Seymour as zero work); the nearest results are THM-4526
  (largest oriented graphs in every tournament, citing Linial–Saks–Sós) and
  THM-4533 (parity of embedding counts). THM-4523's rigidity says the
  minimum-map forest, unlike an abstract tree, recovers its labels.
- Mills. `05-knowledge/results/seam_mills_20260925.md` already gives the exact
  correspondence between prime-shell chains and nested-interval constants for
  `floor(A^(c^n))` (cubic children of 2: `11, 13, 17, 19, 23`); no ceiling
  variant exists. Tóth's version has shells `((p-1)^c, p^c]`, so the ceiling
  children of `p` are the floor children of `p - 1`: the two trees are one
  tree with shifted parent labels, and the endpoint `(p-1)^c + 1` is divisible
  by `p` exactly for odd `c`, so even `c` can keep a prime endpoint (`p = 3`,
  `c = 2`, endpoint 5). A two-line remark, recorded, not computed.
- Renewal in the thread. Lemma R' and Theorem C (five-mirrors note, Fourier
  renewal over ceiling levels, conditional on square-root cancellation),
  THM-4519 (mac-mini, exact Lemma R'', Theorem C' conditional on HYP-9166),
  THM-1770 (GMC2, parts B–D retracted: distinct atoms can cancel), THM-4495
  (ballot factor `k^(-3/2)` in the no-descent count, the polynomial factor
  seen in the slow convergence of `log_2 N(A)/A` to `h*`), THM-4512, THM-4517,
  THM-4530. The present note's renewal is the arithmetic one (blocks of the
  parity vector), not the Fourier one; the two meet at the stopping-time law.

## 7. Reproduction

[Script](../../04-computation/experiments/collatz_minimum_map_renewal_20261005.py),
[output](collatz_minimum_map_renewal_20261005.out),
[JSON](collatz_minimum_map_renewal_20261005.json):

```text
python3 04-computation/experiments/collatz_minimum_map_renewal_20261005.py --census-bits 22 --json 05-knowledge/results/collatz_minimum_map_renewal_20261005.json
python3 -O 04-computation/experiments/collatz_minimum_map_renewal_20261005.py --census-bits 22
```

The DP (exact fractions) gives the block law to length 200; numba computes the
chains of all odd sources below `2^22` (two seconds), the block types of the
first six blocks, the visit counts and the in-degrees; the Hankel determinants
are exact (sympy over fractions). Hostiles: Wald's identity `E[A] = 2E[l]` and
the excess identity are checked; the renewal total-variation bounds, the
nonvanishing of the Hankel determinants and the failure of the Somos-4 fit are
asserted; the fitted exponents are reported against both predictions.
