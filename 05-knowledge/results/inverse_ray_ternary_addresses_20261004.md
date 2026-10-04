# Inverse rays with ternary addresses and certificates that need not expand

2026-10-04. **PROVED** elementary parameterizations, address lifts, and
first-hit certificate rules below; **FINITE-EXACT** for the declared controls.
**OPEN:** whether the resulting rooted grammar contains every positive odd
integer. No novelty is claimed for inverse Collatz trees or their ternary
odometer. The contribution is an explicit, tested certificate representation
and its modular evaluator, with the relevant losses and guards retained.

The practical result is a node carrying `(parent, source row, block height)`.
It defines an actual integer, an exact edge to its parent, and a route to 1.
It can answer binary and ternary residue queries and accept further inverse
nodes without expanding that integer. The finite demonstration includes a
21-odd-step certificate whose source has more than `10^102` binary digits.
That is a symbolic integer with a proved route, not an expanded computation.

## 1. Inheritance: what is already proved

The incoming [sixth-clock note](sixth_clock_branches_20261004.md), section 2,
recovers the three rows with first nonroot sources 21, 5, 85 and distinguishes
ordinary Collatz, single-halving H(n)=(3n+1)/2, and accelerated odd U(n).
Its section 7 already proposes attaching certified target pointers to inverse
rays. The corrected convention is recorded in `01-canon/MISTAKES.md`.

The [September 17 braid note](arithmetic_braids_20260917_collatz.md), section 2,
already proves the full inverse fibre, R(n)=4n+1, and
`v3(R^t(n)-n)=v3(t)`. Its
[inverse-completion companion](arithmetic_braids2_20260917_inverse_completion.md)
already completes arbitrary finite valuation words to prescribed targets and
proves their full finite ternary support. The
[reverse-tree study](collatz_mod6_20260917_reverse_tree_pieces.md), sections 1
and 5, already gives the mod-9 branch table, one-digit precision cost per
generation, and minimal-child survival laws. These are inherited results,
not new consequences of finding the number 63 again.

Closest mechanism: guarded inverse rays and a retained terminal certificate.
Hostiles: the root self-loop, a row-3 leaf, and two equal low-precision target
addresses with different deeper continuations. Corrected near miss: a full
residue cycle does not prove coverage of fixed integers. Least-used sidecar:
the quotient of an unbounded block height by a finite ternary address period.
The live board is hub / row / height / ternary carry / first-hit rank / size.

## 2. Three channels of a complete inverse fibre

Use U(n)=(3n+1)/2^v2(3n+1) on positive odd integers. Fix positive odd u with
3 not dividing u. Write r=0,1,2 for source residue modulo 3; these mean the
odd residue rows 3,1,5 modulo 6, respectively.

There is exactly one exponent kappa_r(u) in {1,...,6} satisfying

    2^kappa u = 1+3r mod9.                              (1)

Indeed 2 has order 6 modulo 9 and runs through all its units. The exact table is

| u mod9 | kappa_0 | kappa_1 | kappa_2 |
|---:|---:|---:|---:|
| 1 | 6 | 2 | 4 |
| 2 | 5 | 1 | 3 |
| 4 | 4 | 6 | 2 |
| 5 | 1 | 3 | 5 |
| 7 | 2 | 4 | 6 |
| 8 | 3 | 5 | 1 |

**PROVED full-channel parameterization.** All odd positive sources in row r
whose U-image is u are exactly

    n_b=(2^(kappa_r(u)+6b)*u-1)/3, b>=0.               (2)

They are integral and odd; `3n_b+1=2^(kappa+6b)u` has exactly that valuation
because u is odd. Conversely that identity for any predecessor fixes its
exponent modulo 6 by (1), giving its unique b. Thus

    n_(b+1)=T(n_b), T(n)=64n+21=R^3(n).               (3)

Rows r=1,2 can be used as new inverse targets. Row r=0 cannot: U never
has an image divisible by 3. A row-zero source still has a completed forward
route whenever u does; being an inverse leaf is not a forward failure.

For u=1, the r=0 and r=2 first members are 21 and 5. The r=1 first member
is 1, whose self-return must be removed for a first-hit certificate. Its
next member is 85. Since H(n_b)=2^(kappa+6b-1)u, this gives the requested
ordered H exponents 5,3,7 in source rows 3,5,1 modulo 6. No sign change is used.

## 3. Ternary address lifting uses one carry subtraction

For every integer n and positive t, (3) gives

    T^t(n)-n=(64^t-1)(3n+1)/3,
    v3(T^t(n)-n)=1+v3(t).                            (4)

The factor 63=64-1 has valuation 2 at 3; the division by 3 leaves the
one-digit shift in (4). This is the relevant retained prime-power depth.
It is not a claim about a new prime divisor of a different sequence.

For a>=1 and d in {0,1,2}, the more precise congruence is

    n_(b+d*3^(a-1)) = n_b+d*3^a mod3^(a+1).          (5)

**Proof.** The binomial expansion of `64=1+7*3^2` gives
`64^(d*3^(a-1))=1+d*3^(a+1) mod3^(a+2)`; this also follows by successive
cubing, keeping the first nonzero ternary coefficient. Substitute in (4)
and use `3n_b+1=1 mod3`. The d=0 case is immediate. In particular, the
coefficient of the new digit is 1, independently of u and the current b.

It follows that, for fixed u,r,a, the map

    b mod3^(a-1) -> n_b mod3^a                       (6)

is a bijection onto the residues congruent to r modulo 3. Its period is
exactly `3^(a-1)` by (4). Equivalently, the inverse exponent has period
`6*3^(a-1)=2*3^a`, matching the inherited order of 2 modulo `3^(a+1)`.
The extra power of 3 is necessary because forming a source divides by 3.

For a desired source address s modulo `3^A`, set r=s mod3 and start b=0
at precision a=1. Suppose b already matches s modulo `3^a`. Compute n_b
modulo `3^(a+1)` and put

    d=((s-n_b)/3^a) mod3,
    b <- b+d*3^(a-1).                               (7)

Equation (5) proves that this chooses the next digit without searching the
whole clock. Continue to a=A-1. The result b_0 in `[0,3^(A-1))` describes
the complete requested address family

    b=b_0+3^(A-1)t, t>=0.                           (8)

For a first-hit family over root 1, exclude t=0 exactly when r=1,b_0=0.
At precision A this shifts the first allowed original R-index from 0 to
`3^A`, recovering the all-depth root offset in the incoming note.

The division precision is explicit: to obtain n modulo `3^a`, one needs
u modulo `3^(a+1)` (at least mod9 for kappa). Choosing the *next* digit in
(7), modulo `3^(a+1)`, therefore needs u modulo `3^(a+2)`. More inverse
generations and more digits within one fibre are different operations.

Among row-zero block indices the depth distribution is exact: within a
complete sufficiently long ternary period, a fraction `2/3^d` has
`v3(n_b)=d`, for d>=1. This follows by counting the addresses in (6), or
by lifting the unique zero address one digit at a time. In the root row,

    n_b=(64^(b+1)-1)/3,
    v3(n_b)=1+v3(b+1).                              (9)

These are frequencies in the **block index**, not densities of source
integers. A fixed ray grows exponentially and has natural density zero.

## 4. A canonical certificate that can stay compressed

Start with the root certificate for integer 1. A nonroot node stores

    (parent certificate, row r, nonnegative block height b).

Its integer semantics is (2), applied to the parent's integer. Require
the parent not to be divisible by 3, and reject the node whose parent is
root 1 and whose fields are r=1,b=0. No other positive first-hit exception
is needed. The implementation validates these structural guards before
accepting supplied nodes, including nodes constructed outside its helper.

**PROVED first-hit invariant.** Every valid node represents a positive odd
integer n>1 whose exact U-image is the parent's integer. Its first-hit odd
rank equals its depth above root. Its ordinary Collatz first-hit rank is

    sum_edges (kappa_r(parent)+6b+1).                 (10)

Induct on depth. Integrality, positivity, oddness, and exact valuation follow
from (1), (2). The equation n=1 would force parent=1 and exponent=2,
precisely the rejected node. The parent already has a first-hit route;
prepending the exact odd edge therefore adds one odd step. Before arriving
at its odd target, the ordinary edge passes through even integers, so it
cannot visit 1 early unless the target itself is 1. It adds exactly k+1
ordinary steps. This proves (10) and also rules out repeated nonroot states.

Distinct valid rooted codes represent distinct integers. Applying U to an
equal source recovers the same parent and exact exponent; the row and
block height are then uniquely recovered from (1), (2). First-hit depth
prevents equality between different depths. Induction gives identical codes.

Conversely, the finite first-hit route of any convergent odd positive integer
decodes into this grammar by reading its exact valuations backwards. The
image is therefore **exactly the positive odd basin of 1**, with unique
first-hit codes. This is a representation theorem for that basin, not a
proof that it contains every positive odd integer. The forward encoder has
a declared cap; no totality of that search is assumed.

### Evaluating residues without expanding the source

At ternary precision a>=1, request the parent's residue modulo `3^(a+1)`;
derive kappa from it modulo 9, and evaluate

    n mod3^a = ((2^k*u-1) mod3^(a+1))/3.             (11)

The numerator's canonical residue is divisible by 3. At binary precision H,
3 is invertible, so

    n mod2^H = (2^k*u-1)*3^(-1) mod2^H.             (12)

If k>=H, the power of two in (12) is already zero. Thus neither projection
requires the source integer. At depth D a ternary query of precision a
ultimately uses root precision a+D; the increasing modulus is retained,
not hidden inside a claim of a fixed finite-state machine. Exact exponents,
both first-hit ranks, and further permitted nodes are computable from this
finite representation. Modular exponentiation takes a number of squarings
controlled by the bit length of k, not by the integer `2^k`.

The controls start with the certificate for 5 and append 20 alternating
internal rows with heights `10^100+i*10^30`, i=0,...,19. This has odd rank 21;
its source has more than `10^102` binary digits. Binary residues modulo
`2^64`, ternary residues modulo `3^12`, and every stored edge are checked
without constructing the source. Appending a row-zero node gives a completed
leaf of rank 22. Structural induction proves these are actual integer routes;
finite residue tests alone would not establish that conclusion.

For completeness the evaluator also provides rigorous bit-length bounds.
If the parent has L binary digits and the edge exponent is k, then

    k+L-2 <= bitlength(n) <= k+L-1.

The valid positive inputs have k+L>=3; these inequalities follow from
`2^(L-1)<=u<2^L` and (2). Propagating the two bounds separately avoids any
logarithmic rounding. The huge internal certificate's bounds differ by
only 21 bits while their magnitude exceeds `10^102`.

## 5. Concrete hostiles and a completed family through 7

One generation's target residue does not decide every later generation.
The certified hubs 11 and 29 are equal modulo 9. Taking the smallest
admissible exponent at each inverse step gives

    11 -> 7 -> 9,     29 -> 19 -> 25.

The first children share row 1, but the second is a leaf in the first path
and remains internal in the second. The missing coordinate is another
ternary target digit. This is the inherited precision-consumption boundary,
shown here without using the root self-loop as the example.

A finite source address also loses ordinary height: b and
`b+3^(a-1)` give the same residue modulo `3^a`, but different actual integers.
The discarded quotient t in (8) must stay in a certificate that identifies
one integer. Similarly, attaching a row-zero child is valid completion;
attempting to use it as an inverse target is invalid. The verifier rejects
that attempt, the root self-return, and tampered rows or negative heights.

There is an exact family through the preceding session's structural miss 7:

    n_j=(22*4^j-1)/3, j>=0,
    7,29,117,469,... -> 11 -> 17 -> 13 -> 5 -> 1.     (13)

The displayed initial list consists of alternative sources, each mapping
directly to 11. Every member has exactly five odd steps and `16+2j`
ordinary steps to 1. The entire fibre lies outside the
[universal twice-repeated favorable atlas](periodic_chart_separation_20261004.md).
For j>=1 its first valuation is at least 3, whereas a favorable growing
word must start with valuation 1. For j=0, source 7 has first descent at
time four and exact initial valuations (1,1,2,3); no twice-repeated nominal
word can match its first three entries at that time.

This is not new positive-density coverage relative to the old bank. In fact
the whole cylinder `7 mod128` already belongs to the frozen swap-lift bank
(rows indexed 25,89,153,217,281). It has first-descent time four: for
`n=7+128t`, the first three odd nodes are
`11+192t,17+288t,13+216t`, all greater than n; the fourth is
`oddpart(5+81t)<n`. The first three valuations are always (1,1,2).
Thus being missed by one representation does not imply novelty relative
to an older certificate bank. Equation (13) adds a compact completed suffix
representation, not a claimed new density increment.

Finally, for a fixed hub and fixed address precision, (8) gives only
`log X/(3^(a-1)*log64)+O(1)` sources up to X. Full ternary support and a
short symbolic certificate coexist with sparse ordinary-size coverage.
The construction does not classify unknown negative cycles. Reducing the
integer coefficient 2 into the sextic residue field F64 would make it zero,
destroying the exponent clock; its separate unit t+1 is not the same map.
Any bridge to that 63-cycle needs an explicitly marked clock map, not the
shared word "binary" or the factor 63 in (4).

## Reproduction and exact scope

Run from the repository root:

```text
python -X utf8 -B 04-computation/experiments/inverse_ray_ternary_addresses_20261004.py
python -X utf8 -B -O 04-computation/experiments/inverse_ray_ternary_addresses_20261004.py
```

The [script](../../04-computation/experiments/inverse_ray_ternary_addresses_20261004.py)
and [output](inverse_ray_ternary_addresses_20261004.out) retain:

- all 14,760 hub-residue / row / block-address cases at precisions 1--4,
  comparing power evaluation, independent affine orbit iteration, and the
  one-subtraction inverse; every digit lift is checked separately;
- four certified hubs, precisions 1--8, four requested addresses per depth,
  and three lifts per address; 300 root-leaf valuation controls;
- all 1,536 valid codes through depth five with rows 0,1,2 and block heights
  0,1, stopping at leaves and removing the root self-loop; independent literal
  forward replays, source/code round trips, and binary/ternary projections;
- the two enormous symbolic codes and their exact structural guards;
- 51 members of (13), 100 members of the inherited `7 mod128` cylinder,
  and the precision, quotient, leaf, and tampering hostiles above.

Normal and optimized checks use the same explicit exceptions. No floating
arithmetic, probabilistic primality test, or assumed global convergence is
used. The improved object stores a checked route and its address operations;
finding such an object for an arbitrary input remains the global obligation.
