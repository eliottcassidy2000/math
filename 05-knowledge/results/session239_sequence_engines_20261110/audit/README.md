# Independent audit of Session239

Status: PROVED audit of stated identities; FINITE-EXACT independent controls.
No commits or edits outside this audit directory were made by the auditor.

## Verdict

The factorial-scaled response recurrence, its two residue modes, the few-cycle
precision bound, both proposed congruences with their stated lower-order
boundaries, and the value-refined fixed-permutation SCC recurrence pass.
No fastest-known or polynomial-time full-count claim is warranted.

Reproduce from repository root:

    python3 05-knowledge/results/session239_sequence_engines_20261110/audit/verify.py
    python3 -O 05-knowledge/results/session239_sequence_engines_20261110/audit/verify.py

Both runs pass. The checks use explicit exceptions, never optimized-away
assertions. The matching saved outputs are `verify.out` and
`verify.optimized.out`.

## Independent arithmetic controls

The auditor constructs odd-cycle multiplicity profiles independently and uses

    sum_(i,j) gcd(c_i,c_j) = sum_d phi(d) (#cycles divisible by d)^2

rather than either engine's pairwise insertion formula. The fixed-point
exponent is half this quantity minus half the number of cycles. Class sizes
n!/z are evaluated exactly as integers, followed by one terminal division.
All 46 terms n=0..45 agree. There are 16,857 checked profiles and 8,971 odd-order
profile-bound checks. The inherited independent reference n=0..76 also agrees.

Horizon consistency is checked at 0,1,2,3,4,7,12,31,64,120. At each horizon,
all prefix terms are checked modulo 1,2,4,6,9,72,1001 and modulo powers of two
with 1,2,3,7,32 bits. Five invalid-parameter cases are rejected. This probes
empty/singleton universes, noninvertible denominators, vanishing residues,
large factorial valuations, and coordinate truncation across horizons.

The integral multiplier proof in `../burnside/PROOF.md` is correct:
C(d,l,k)=binom(d+lk,lk)(lk)!/(l^k k!). It requires that this length has not
previously been inserted, exactly as the code ensures. A canonical residue
modulo M N! remains divisible by n!, because n! divides the modulus. For
binary residues, dividing by 2^v2(n!) leaves at least the requested precision,
and the odd factorial part is invertible. Removing zero states and shifts
at least the working precision is therefore safe. These are residue-only
prunings and do not provide the full integer.

## Precision theorem and congruences

For odd n, an odd-part partition has odd length r. If G is its cross-gcd sum
and m_j its multiplicities, its contribution to U_n=a(n)/2^((n-1)/2) has
valuation

    delta = G-(r-1)/2-sum_j v2(m_j!).

Since r!/product(m_j!) is an integer, the last sum is at most v2(r!). Since
every pairwise gcd is at least one,

    delta >= (r-1)^2/2-v2(r!).

Thus only cycle counts with this bound below b can affect b bits. The weaker
bound (r-1)(r-3)/2 follows from v2(r!)<=r-1 and makes r=O(sqrt(b)) explicit.
For fixed b, enumerating the first r-1 parts determines the last, so the
stated O(n^(rmax-1)) candidate bound holds. This is a candidate bound, not
an assertion of constant-cost integer arithmetic or a polynomial-time
algorithm for the full quadratic-bit count. Use modular odd-denominator
arithmetic after extracting each term's powers of two.

For mod8, all r>=5 terms vanish: r=5 already has lower bound5, and the weaker
bound is at least4 for r>=5. For r=3:

* Three distinct parts have delta=G-1. They survive iff all three pairwise
  gcds are1, when delta=2 and their odd denominator contributes4 modulo8.
* Parts (x,x,y), x!=y, have delta=x+2gcd(x,y)-2. It is odd and below3 iff
  x=1. For n>=5 this gives (1,1,n-2), contributing 2/(n-2).
* Three equal parts x>=3 have delta=3x-2>=7 and vanish.

The single n-cycle contributes 1/n. Hence the parent's formula

    U_n = n^(-1)+2(n-2)^(-1)+4 C_n (mod8), odd n>=5,

holds, with C_n counting unordered distinct positive odd pairwise-coprime
triples of sum n. Reducing modulo4 gives U_n=-n (mod4). At n3, direct U3=1
verifies the mod4 statement but the displayed mod8 expression would give5;
the exceptional all-one triple explains this failure. At n1, even the mod4
formula fails, so its stated n>=3 boundary is necessary.

Independent exact integer tests cover mod4 at every odd n=3..119 and mod8
at every odd n=5..119. THM-305,
`01-canon/theorems/THM-305-two-adic-valuation-tournament-count.md`, already
proves the underlying valuation; this is a precision/congruence refinement,
not a new valuation theorem. Targeted relevant-file searches found no earlier
mod4/mod8 statement, but this is not a comprehensive priority search.

## Value-refined SCC recurrence

For a fixed permutation sigma, every invariant tournament has a unique first
strong component. Automorphisms preserve the total order on tournament SCCs,
so the component is a union of whole sigma-cycles. Choosing b_d cycles among
m_d contributes product binom(m_d,b_d). The remaining tournament is arbitrary
and fixed by the restricted permutation. All cross-arcs have the forced
forward direction. Every Hamiltonian path traverses components in this order,
so H is multiplicative and the correct coefficient algebra uses [h][k]=[hk].
Subtracting proper nonempty first-component profiles therefore proves the
recurrence. Ordinary polynomial multiplication would be incorrect.

The audit separately enumerates every labelled tournament through order5,
counts paths by permutations (not subset DP), recognizes strong connectivity
by forward/reverse set closure (not the worker's bit-reach routine), and
checks fixedness by relabeling arc sets. All18 nonempty cycle types pass the
SCC recurrence. The auditor reuses only the worker's permutation constructor
and subprofile iterator, not its orientation generator, path counter,
strong-connectivity counter, convolution, or scalar recurrence. Fixed H is
also divisible by permutation order in every tested case. The general reason
is the free action of a permutation on ordered Hamiltonian vertex sequences.

## Benchmark scope

The three-repeat n40,60,80 comparison against the inherited Fraction response
implementation is a valid same-task prefix comparison. Integer-core versus
naive direct Fraction is a valid stated baseline but not the best competing
engine. The inherited v3 output aborts at n2 and is not a valid speed baseline.
Historical machine timings cannot establish current speed ratios.

The new extended benchmark appropriately adds the historical compiled GMP
v4 evaluator: one term versus an entire prefix, explicitly different tasks.
Recorded single-term/full-prefix times are approximately .019/.050s at80,
.117/.127s at100, and .698/.319s at120. These support the reported observations
only, not universal dominance. This auditor inspected the benchmark code and
outputs, but did not independently rerun its compiled-GMP timing experiment.
Full-prefix integer response versus inherited Fraction is the clean claim.
