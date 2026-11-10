# Factorial-gauge sequence engines

Status: PROVED elementary implementation identities, FINITE-EXACT benchmarks;
no literature priority, polynomial complexity, or fastest-known claim.

## Inheritance and boundaries

The source is Davis/Burnside. Fixed-point stripping and the integer kernel
K[m,t] are inherited from THM-514, `01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md`.
The future-response quotient is inherited from
`05-knowledge/results/tournament_future_boundary_20261020/README.md`, §§2–3.
The Gaussian solution, pole expansion and non-P-recursiveness are inherited
from THM-4845, `01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md`.
The old `04-computation/a000568_int_v3.py` already avoids rational arithmetic
using a global denominator budget. Our implementation instead scales each
state by its own degree factorial, retaining the newer finite-horizon response
compression. This is a representation/compiler improvement, not a new
Burnside formula. No closed core kernel is claimed.

Anchor: exact efficient A000568 prefixes. Niche: composite-modulus residues.
Wildcard: factorial scaling as an integral divided-power implementation.
Concept board: future response, integral state gauge, core/tail separation,
precision filtration, information lost by counting.
Canonical hostile: equal degree and cycle count need not imply equal gcd
responses. Corrected near miss: truncating high or low exponents does not
recover a full integer without an absolute bound. Least-used sidecar here is
the degree factorial, which turns rational response transitions into integers.
No transfer from sequence counts to an individual tournament or LRC is made.

## 1. Integral future-response recurrence

Fix a maximum degree N. Process the odd lengths in descending order. A history
mu is a multiset of processed cycles; d=|mu| and

    e(mu) = sum_(a in mu) (a-1)/2 + sum_(unordered pairs a,b in mu) gcd(a,b),
    z(mu) = product_l l^(m_l) m_l!.

For the remaining lengths j <= N-d, retain R_j=sum_(a in mu) gcd(a,j).
At a state (d,R), define

    W(d,R) = sum_(histories in state) (d!/z(mu)) 2^e(mu).

Every summand is an integer: d!/z(mu) is a conjugacy-class size in S_d.
Initially W(0,0)=1. Adding k cycles of the next length l changes degree to
D=d+lk, responses to R'_j=R_j+k gcd(l,j), and exponent by

    Delta = k(l-1)/2 + l k(k-1)/2 + k R_l.

The exact transition is

    W(D,R') += W(d,R) C(d,l,k) 2^Delta,
    C(d,l,k) = (d+lk)! / (d! l^k k!).

C is an integer: choose the lk new labels and then give them a permutation of
cycle type l^k, so C=binom(d+lk,lk)(lk)!/(l^k k!). Algebraically it is also the
ratio of new and old class sizes, since length l has not previously occurred.
Equal (d,R) histories have identical Delta and C for each possible insertion,
so adding their integer weights preserves all continuations. Coordinates
j>N-D can be discarded because those lengths can never be inserted. Terminal
states have empty response and W(n,empty)=n! a_n. Exact division yields all
coefficients 0..N simultaneously. This proves the code, including n=0.

This is a diagonal change of basis on the inherited rational transfer
operator, W_d=d! w_d. The operator becomes integral. It uses exactly the same
state graph as the inherited quotient, except zero states may be removed in
residue mode. It does not establish a polynomial bound on state count.

## 2. Safe residue arithmetic, including noninvertible denominators

For any positive integer M, run the integer recurrence in Z/(M N!)Z.
Every transition is integral, so its terminal value is congruent to n!a_n.
Because n! divides M N!, the canonical residue is divisible by n!. Division
by n! and reduction modulo M recover a_n modulo M. This works even when M is
composite or shares prime factors with n!. Dividing a generic Burnside
fraction modulo M is invalid in exactly that hostile case.

For M=2^b, only L=b+v2(N!) bits are needed, where
v2(N!)=N-popcount(N). Run modulo 2^L. At terminal degree n, let
v=v2(n!) and u=n!/2^v, odd. The canonical residue w is divisible by 2^v, and

    a_n mod 2^b = (w/2^v) u^(-1) mod 2^b.

Indeed the discarded multiple of 2^L becomes a multiple of 2^(L-v), hence a
multiple of 2^b. A transition with Delta>=L is zero and may be skipped. A
state whose accumulated residue is zero may also be removed: all subsequent
maps are integer-linear. Skipping these terms is exact for the residue only,
not for the full integer. This engine complements the parent's few-cycle
precision theorem; it neither assumes nor reproves that theorem.

## 3. Integral core-prefix implementation (inherited formula)

Enumerate nondecreasing odd parts >=3 only. For core mu of mass m, t cycles,
exponent e and class size c=m!/z(mu), accumulate K[m,t]+=c 2^e.
If appending length r produces multiplicity k of that length, update

    c' = c (m+r)!/(m! r k),
    e' = e+(r-1)/2+sum_(s in mu) gcd(r,s).

Although the multiplicative factor alone need not be integral, c' is a
conjugacy-class size, so the full numerator is exactly divisible. The code
checks this. An empty core is allowed. Finally,

    n! a_n = sum_(m,t) binom(n,m) K[m,t]
               2^(binom(n-m,2)+(n-m)t).

Each nonunit core is visited once for the entire prefix, rather than once per
possible count of fixed points. This implementation remains partition-bound;
its core-node count through N equals the number of odd partitions of N
(by filling unused mass with ones). At N=80 the compressed integer response
engine is faster, so core stripping is not presented as the best engine.

## 4. Reproduction and evidence

From repository root:

    python3 05-knowledge/results/session239_sequence_engines_20261110/burnside/benchmark.py --orders 40 60 80 --repeats 3
    python3 05-knowledge/results/session239_sequence_engines_20261110/burnside/validate_extended.py

`benchmark.json` records Python/platform, median wall times, state counts,
sequence hashes and final integers. Four paths agree: direct Fraction Davis,
inherited Fraction response quotient, integer response quotient, integer
core-prefix DFS. Direct known-prefix controls run through n=15. Extended
validation checks the inherited external reference through its full available
range and compares against historical compiled GMP v4 (single term rather
than prefix). Timings are machine-specific; no cross-machine comparison.

Proof consequences: arbitrary-composite-modulus evaluations are now safe and
exact without modular division by a nonunit; low-bit experiments can be
turned into falsifiable conjectures without calculating all quadratic-size
output bits. These computations alone do not prove observed congruences.
The finite-state recurrence still grows with N, consistent with the inherited
non-P-recursiveness theorem.

### Measured comparison against a stronger historical implementation

`validate_extended.py` retrieves the unchanged historical C/GMP source with
`git show HEAD:04-computation/a000568_gmp_enum_v4.c`, compiles with gcc -O3,
and runs one thread. Median of three process-inclusive single-term times:

| N | historical C/GMP single term | new Python entire prefix |
|---|---:|---:|
| 80 | 0.0191 s | 0.0502 s |
| 100 | 0.1171 s | 0.1267 s |
| 120 | 0.6977 s | 0.3194 s |

The historical compiled program wins at 80 and slightly at 100 despite its
partition enumeration. The response quotient crosses over in this sample at
120 and also provides every preceding term. This is not an exhaustive ranking
of all historical implementations (batch CRT, multi-threaded variants, and
GMP optimizations were not exhaustively benchmarked).
The 16-bit prefix at N=120 takes 0.0904 s, with 1062 peak states and 37530
nonvanishing transitions versus 6776 states and 202287 transitions for full
integers. Timings include the Python integer engine once in this extended
comparison; the dedicated benchmark uses medians of three.

`scale.py` separately compares the new and inherited compressed engines at
N=120,160,200, checking every coefficient, with source/output saved together.
Its large-order timings are single-run observations, not medians. None of
these experiments asserts complexity improvement in N for the state graph;
the full-precision state graph is the same inherited graph.
