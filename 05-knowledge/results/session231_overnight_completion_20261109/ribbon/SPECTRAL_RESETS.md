# Short fixed ribbon blocks cannot support exact doubling orbits

Status: FINITE-EXACT matrix certificates, with PROVED all-time consequence for each tested block; independently audited ([audit](../ribbon_audit/AUDIT.md)). No all-block impossibility theorem is claimed.

## Input and controls
Inherited exact integer ribbon representations from THM-4880, directory `04-computation/experiments/ribbon_spectrum_20261104`: span2 has dimension14 and span3 dimension62. Span2 matrices are rebuilt by the inherited exact module-invariance code `intrep.build(2)`; span3 matrices are read from `span3_intrep.txt`. The scan checks det(M_B-aI) modulo the prime1,000,003. The original span2 producer filters modular-singular candidates through exact determinant tests, so its empty `found` list alone establishes only absence of exact eigenvalues. The independent audit replay explicitly counts modular singulars and certifies zero in every listed universe; consequently all listed determinants are nonzero modulo that prime and over the rationals. Products at these lengths fit int64; the largest entries are far below2^63. The script uses modulus after each elimination step, with products bounded by10^12. No floating point is used.

`search_resets.py`: all5,460 nonempty blocks of lengths1..6 over all four span2 letters; radices2,4,8. No eigenvalue candidates.

`search_span3.py`: all584 nonempty blocks of lengths1..3 over all eight span3 letters; radices2,4,8. No eigenvalue candidates.

Span2 letter indices equal bit masks:0 empty,1={1},2={2},3={1,2}. The span3 stored matrix order is masks[0,1,2,4,3,5,6,7]; index-based span3 words must use that mask list. The exhaustive scan is unchanged by this permutation.

`search_long_resets.py`: all29,523 blocks of lengths1..9 over the three nonempty span2 letters1,2,3; radices2,4,8. No modular singular candidates. This larger restricted-alphabet scan does not cover blocks containing letter0 beyond length6.

## Exact implication
Let a>1, let M be the matrix of a fixed block, and suppose det(M-aI) is nonzero. No state vector y and scalar d can satisfy an unbounded scalar orbit

    h M^t y = a^t c - d/(a-1), c != 0, for every t>=0.

Proof: let p be the characteristic polynomial of M. Cayley–Hamilton implies the scalar sequence satisfies p(E)s=0. Applying it to the proposed expression gives

    c p(a) a^t - (d/(a-1))p(1) = 0

for every t. Since a>1 and c!=0, p(a)=0, contradiction.

Therefore none of the scanned blocks can be repeatedly applied on an invariant family with nonconstant, growing exact radix2,4,8 affine recurrence. This holds even if all prefixes are strong. It does NOT exclude scalar-degenerate repeated branches, variable-length blocks, switching between different blocks, or longer unscanned blocks. In a proposed two-block strong-state digit family, any branch whose repeated scalar recurrence is nonconstant must evade the scan's length/span bounds.

## Hostile probe against an easy structural explanation
A tempting explanation is that every block's characteristic polynomial has odd lowest nonzero coefficient. It is false. Span2 word

    [3,1,2,3,3,3,1,2,1]

has characteristic polynomial

    x^11 (x^3 - 971x^2 - 1302x + 96).

Its stable matrix rank is3 over Q and1 over F2. Hence an odd-constant-term argument cannot promote the finite eigenvalue exclusion to all lengths. The exact polynomial is independently computable with SymPy from the inherited matrices. `parity_probe.py` found this with deterministic seed231 on sample50; the probe is discovery only, and the displayed exact characteristic polynomial is the decisive witness.

## Reproduction
Set MATH_REPO to a complete checkout containing the inherited matrix sources, then run with PYTHONDONTWRITEBYTECODE=1:

    python3 search_resets.py
    python3 search_span3.py
    python3 search_long_resets.py
    python3 parity_probe.py

The all-time statements above are exact implications of finite certificates, not finite samples of trajectories.

The independent audit also verifies the stored quotient matrices through exact intertwiners with the inherited profile representation, certifies all integer overflow bounds, detects an eigenvalue3 positive control, and checks the hostile polynomial and stable ranks exactly. See its `check.py` and `RESULTS.jsonl` for the complete replay.
