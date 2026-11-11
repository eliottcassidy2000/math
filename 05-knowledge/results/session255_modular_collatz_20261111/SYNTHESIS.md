# Central tables become exact edge atlases and finite machines

**Session255, 11 November 2026.** New proofs and constructions are separated
from inherited identities and finite controls. This session follows the user's
central-table, mod-eight and prime/composite ideas into concrete finite tools.
It does not claim a proof of Collatz, Chan–Pak, JC(2), general LRC, or the
Hadwiger–Nelson target, and makes no literature-priority claim.

## 1. The central picture is retained, with its exact convention

For the nonzero-index multiplication table modulo odd n, the four central
products are the complementary pair 4^(-1), −4^(-1), each repeated twice.
Thus the user's formulas are exactly the inherited central blocks:

    n=4N+1: [[3N+1,N],[N,3N+1]],
    n=4N−1: [[N,3N−1],[3N−1,N]].

Each central row sums to n as an INTEGER. The natural augmented display
therefore retains the two copies of n as row sums beside the four products.
This is a proposed reading of the six-number description, not an assertion
that a four-cell block has six entries or that its products include zero.

The older [paired-prime note](../paired_prime_collatz_20261003.md) already
connects these centers to adjacent-source shortcut pairings and odd
multiplication a*b=a+b+2ab. Those identities are credited, not rediscovered.
The present single-source raw-edge atlas uses different source labels from
that shortcut pairing; the apparent sign reversal is resolved explicitly.

## 2. Every signed Collatz edge has one table code

Fix s=+1 or −1 and use the ordinary map C_s(x)=3x+s on odd x, x/2 on even x.
Assign source x the table modulus **m=4x+s**.

- Odd x: its central values x and 3x+s directly give the odd arrow.
  These are the classes 5 mod8 for the plus map and 3 mod8 for the minus map.
- Even x=2h: the center is 2h,6h+s. Reverse the pair, subtract s before
  division by6, and obtain h→2h. Reverse that arrow for the ordinary even
  edge 2h→h. Class1 subtracts1; class7 adds1.

The parity flag, decoder and orientation give a bijection between all sources
and all graph edges. “Half” is the partition by source parity, not a claim
about equal edge counts inside each orbit. This is the user's four-family
picture with its exact map and missing coordinates supplied.

There is also an infinite inverse-edge compiler. For a positive odd target y
not divisible by3, choose k0=1 or2 so 2^k0 y≡s mod3. All accelerated odd
predecessors are

    x_j=(2^k0 4^j y−s)/3,   j≥0.

Both the sources and their table codes obey the SAME recurrence z→4z+s.
Modulo any odd L their exact least period is

    ord_(3L/gcd(y,L))(4).

A finite list of this length, together with j=qP+r, reconstructs every actual
integer predecessor and its full division exponent. The infinite object is
compressed without discarding height.

[Dictionary, convention comparison and proof](collatz_dictionary/README.md).

## 3. The prime/composite split has an exact information geometry

Let p be the least prime divisor of n. The largest centered interval of row
labels whose multiplication channels are all reversible has width **p−1**.
It fills the entire nonzero table exactly when n is prime. Every odd modulus
has the same perfect central seed; composites are distinguished by how soon
that reversible seed hits a singular row.

Multiplication by a modulo n loses exactly log2(gcd(a,n)) bits for a uniform
input with known multiplier. An optimal reconstruction tag has gcd(a,n)
values. At prime powers it records the missing HIGH digits; recording low
digits again does not repair the loss, as n=9,a=3 shows.

Under scaling by units, the nonzero rows have tau(n)−1 defect types. Primes
have exactly one, the smallest possible alphabet. Capped valuation vectors
form a minimal tau(n)-state machine tracking kernel sizes of arbitrarily long
multiplication words. This is a precise interpretation of maximal transfer
and minimal prime complexity, not a raw description-length claim.

There is a sharp boundary: when gcd(n,6)=1, any deterministic observer that
handles BOTH affine sheets 3x±1 and detects zero needs all n residue states.
Their generated maps contain translation by2, so every residue can be
separated by a future word. Multiplicative information compresses to divisor
states; full two-sheet affine location does not.

[Five laws, optimal tags and lower bound](central_table/README.md).

## 4. An infinite ternary tower has one spectral formula

The right infinite carrier is a finite odd-modulus register driven by a
complete binary parity tape. At every finite depth this is an exact CRT
model, not an assumption that one integer orbit behaves randomly.

The crucial product rule is **choose the shared branch first, then combine
prime-power registers**. Modulo15, source0 has only the two outputs0 and8.
Letting the mod3 and mod5 registers choose different branches fabricates
outputs3 and5. This small example identifies an actual lost coordinate.

Write A_q(u) for the weighted shortcut register: even steps have weight1,
odd steps weightu. For every a≥1 and either global sign,

    det(I−z A_(3^a)(u))
      = (1−(1+u)z) product_(j=1)^a (1−z^(2·3^(j−1))).

An exponentially growing family of state graphs has a linear-size factored
description. More strongly, for every q=3^a m with gcd(m,6)=1, all additional
spectral factors are explicit even-branch cycles. They do not depend on the
odd-step weight or sign. Odd transitions lower ternary Fourier conductor;
even transitions preserve and rotate it. This filtration is the structural
reason, and it gives a second proof besides fixed-point counting.

A finite parity-language controller can be included as well; its synchronized
matrix product has rational accepted-path series. The exact all-depth
factorization and its controller extension are in the proof packet.

The arithmetic boundary is explicit. A length-k word with s odd letters has
an ordered numerator B and rational periodic source

    x=B/(2^k−3^s).

It is a positive integer cycle only when this quotient is integral and
positive. Word100 gives x=1/5: every ternary depth accepts it, but modulo5
rejects it. Thus the closed ternary spectrum does not secretly prove an
integer-cycle theorem. The denominator and ordinary order remain essential.

[All-depth formulas, finite controller and realization receipt](finite_infinite/PROOF.md).

## 5. A native Chan–Pak repair search becomes finite without bounding size

The transfer to recent tournament work is an actual compiler, not an analogy
between trajectories. In a common positive d-coordinate frame, compare two
actual source states x≥y. For any nonnegative operation matrix, the support
of their difference evolves exactly by Boolean reachability. Therefore:

**Exact equality of one final coordinate plus any modulo-m response predicate
is recognized by at most 2^d m^d states.** If a legal operation word exists,
one of length at most2^d m^d−1 exists. No bound on the growing integer state
or graph size is needed. The starting graphs and common operation word are
retained as the realization data.

For the current O/Z/J four-coordinate frame the equality observer minimizes
to eleven states, improving the witness bound to11m^4−1. The recent actual
401 repair permits precisely O* count-neutral continuations. Asking also
for U≡1 becomes a finite orbit test. Modulo7 it fails despite a unit source,
because the source lies in the wrong doubling coset; modulo5 one O succeeds.
This distinguishes divisibility failure from phase failure.

Fixed signed-Collatz words have a complementary exact memory law: only their
3-primary component can collapse. After r odd multiplications, the kernel
size modulo m is3^min(r,v3(m)). Composite moduli coprime to3 remain reversible.
The prime/composite distinction is real, but must be attached to the actual
operation rather than used as a blanket claim about all composite channels.

[Finite repair theorem, eleven-state certificate and actual401 controls](cross_frontiers/PROOF.md).

## 6. Connections that survive, and next constructive targets

The reusable move is **retain the shared operation and the exact unpaid
coordinate, then compress everything else**. Here that coordinate is a row
kernel tag, an edge orientation, an inverse-ladder height, a common CRT bit,
a rational denominator, or a native positive difference support.

- HN's geometric conditional multiplication and modular CRT both multiply
  fibres only after fixing the same boundary state. Summing first can produce
  false compatibility; no map equating colourings with numbers is asserted.
- Session254's prime-power substitution interfaces and the ternary register
  both form compatible towers. This session computes the new ternary layer
  explicitly; it does not identify Chan–Pak curvature with Collatz state.
- The finite positive repair compiler applies when a lawful matrix action and
  comparable actual sources exist. It does not cover signed cancellation,
  arbitrary graph resets, JC's unpaid source debts, or LRC phase completion.

Three focused next targets emerge: attach uniform height inequalities to
inverse-ladder records; classify denominator-divisibility word families
using the shared register; and query actual positive repair libraries under
complete downstream interfaces. These are concrete remaining tasks, not
renamed claims that a major conjecture is now one finite check away.

## Evidence and publication

Four runnable engines check only this session's new constructions. Inherited
claims and outside conjectures were not broadly re-audited. The proofs supply
all-height/all-depth assertions; bounded tests are labelled separately.
See README.md for commands, SOURCE_PINS.json for inherited sources,
MANIFEST.json and SHA256SUMS for the saved packet.
