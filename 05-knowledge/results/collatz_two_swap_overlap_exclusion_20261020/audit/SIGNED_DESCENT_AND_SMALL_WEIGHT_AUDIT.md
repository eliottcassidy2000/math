# Independent audit of signed descent and the small-weight overlap closure

20 October 2026. Status: **VERIFIED proof application and FINITE-EXACT independent computations.** This packet audits the new finite descent dependency and the entire reduced overlap residual. The inherited all-weight content estimate and the resulting bound `p<4096` remain separately identified dependencies.

## Verdict and scope

Both signed descent universes through magnitude `N=2^24` pass independent checked-integer verification. The product estimate, strict Legendre inequality, primitive-cycle bookkeeping, exact continued-fraction list, and direct overlap candidate checks also pass. No false positive, missed known negative cycle, integer overflow, or ignored computation failure was found.

Together with the frozen bound `p<4096` in the [overlap-content theorem](../../collatz_two_swap_overlap_lattice_content_20261020/OVERLAPPING_TWO_ONE_CONTENT.md), these inputs exclude every integral primitive `k=1` overlapping two-distinct-one output at `p>=3`. The two `p=2,q=3` paths remain genuine negative-cycle examples, with original-cut fixed points `-7` and `-5`.

This is a reduction followed by 906 direct checks, **not** an enumeration of all 6,959,303,782 paths in the unreduced `p<4096` universe. The separately reviewed Section 6 corollary also excludes every other nonundo overlapping path at any mechanical shape, including nonprimitive sources. It does not exclude arbitrary two edits, disjoint edits, arbitrary nonmechanical sources, or general signed Collatz cycles.

## 1. Independent finite signed descent

The independently written `audit_signed_descent.cpp` imports no producer implementation. Unlike the producer's accelerated odd map and unsigned 128-bit arithmetic, it uses every step of the original half map and checked unsigned 64-bit arithmetic:

- Positive: even `n -> n/2`; odd `n -> (3n+1)/2`.
- Negative magnitude: even `n -> n/2`; odd `n -> (3n-1)/2`.

For each sign it visits **every** start `1<=n<=16,777,216` in increasing order. It skips only the explicitly verified minimum of each known cycle. Every other start must reach a strictly smaller positive integer. Strong induction therefore proves convergence of every tested start to an explicitly verified known cycle.

The base cycles are checked edge by edge, with positive magnitudes, distinct entries, and the minimum first:

- Positive: `(1,2)`.
- Negative: `(1)`; `(5,7,10)`; `(17,25,37,55,82,41,61,91,136,68,34)`.

The negative cycles are magnitudes of the signed integer cycles. No positive/negative sign conversion is silently applied to the update rule.

### Failure behavior

Before `3*n+1`, the code requires `n <= (UINT64_MAX-1)/3`. Before `3*n-1`, it requires `n <= UINT64_MAX/3`. The latter check is conservative: it may reject an otherwise representable subtraction boundary, but can never accept overflowing arithmetic. Even steps cannot overflow. Zero magnitude is rejected. The sum of step counts is separately guarded.

The per-start step cap is one million. Exceeding the cap is a fatal nonzero exit, never a skipped start or an assumed descent. Both cycle errors and incomplete coverage also fail. No `assert` is used.

Internal hostile controls require rejection of zero and an overflowing odd input for each sign. Exact large safe-boundary controls check both odd formulas. A separate compilation with `-DAUDIT_STEP_LIMIT=1` exits nonzero with `FAIL: step limit at start 3`, demonstrating that a failed descent is not silently accepted.

### Results

| Map | Known minima | Certified strict descents | Half-steps | Maximum descent length | First maximizing start | Largest reached half-step value |
|---|---:|---:|---:|---:|---:|---:|
| Positive | 1 | 16,777,215 | 58,576,417 | 287 | 13,421,671 | 30,171,305,459,816 |
| Negative magnitude | 3 | 16,777,213 | 58,604,988 | 273 | 5,960,769 | 15,270,716,514,700 |

The positive maximum value is first reached from `6,631,675`; the negative one from `8,612,097`. Both maximum values are exactly half the producer's recorded maximum unhalved numerators. They are far below the 64-bit limit; width correctness is nevertheless enforced at every transition, not inferred retrospectively from these maxima.

The optimized build and the unoptimized build with UndefinedBehaviorSanitizer produce byte-identical receipts. Compiler: Debian g++ 14.2.0-19. This is an exhaustive finite-exact computation, not a global Collatz convergence claim.

### Producer inspection

The producer's accelerated map also passes inspection and replay. Every odd step computes a checked positive even numerator and removes all powers of two, so its next argument is again positive and odd. Termination at a known minimum is valid even if that minimum exceeds the original start, because those cycles were independently checked before the induction. Every other termination is a strict descent. Even starting values are covered by immediate halving. The loop covers all `N/2` odd starts, and its explicit coverage test distinguishes skipped known minima from tested starts.

Its trajectory arithmetic is checked unsigned 128-bit. Only the documented trace checksum wraps modulo `2^64`; that checksum is not used for a mathematical decision. All aggregate counters are bounded or explicitly guarded. The reproduced producer receipt is byte-identical to its frozen receipt.

## 2. The product inequalities

Let an integer cycle have `q` half-steps, `p>0` odd entries, and minimum positive magnitude `M`. Denote its odd magnitudes by `m_1,...,m_p`. Multiplying the step ratios gives

    2^q = 3^p product_j (1+s/(3m_j)),

where `s=+1` for positive cycles and `s=-1` for negative-magnitude cycles. Put `alpha=log3/log2`.

For positive cycles, `0<log(1+t)<t` yields

    0 < q/p-alpha < 1/(3M log2).

For negative cycles, `0<-log(1-t)<t/(1-t)` for `0<t<1` yields

    0 < alpha-q/p < 1/((3M-1) log2).

The second inequality must use `3M-1`; substituting the positive-sign denominator `3M` would be unjustified. The sign of the resulting clock agrees with the sign of the integer branch fixed point.

The finite descent computation shows that any cycle other than the four explicitly verified cycles has `M>N`. Since `M` is integral, `M>=N+1`, hence `3M-1>=3N+2`. The first positive term of the exact series for `log2` is `2/3`, and the remaining terms are positive. Therefore `log2>2/3`, and both displayed errors are strictly less than `1/(2N)`.

For every integer `p<4096`,

    p^2 <= 4095^2 < 4096^2 = N,

so the error is strictly less than `1/(2p^2)`. When `gcd(p,q)=1`, Legendre's criterion forces `q/p` to be an ordinary convergent of `alpha`. The known cycles are handled separately; none is subjected to the false premise `M>N`.

## 3. Branch parity, primitive periods, and known-cycle controls

For a bit word of weight `p>0`, the carry `B` is positive and the clock `D=2^q-3^p` is odd. The fixed point is `B/D`, so zero cannot occur and its sign is the sign of `D`.

If `B/D` is integral, its parity equals the initial bit: modulo 2, the carry is that bit and the clock is odd. The cyclically shifted branch identity then gives the next integer fixed point for the rotated word. Iteration proves the full parity itinerary. Equivalently, this follows directly from the affine composition identity and the successive parity congruences. The independent script additionally checks every parity and closure on each actual integral output.

A periodic integer orbit has a least positive period dividing every return time. If the given length-`q` word made `d>1` traversals of a shorter cycle, then `d` would divide both `q` and its odd count `p`. Here `gcd(p,q)=1`, so `d=1`. Thus known cycles cannot generate omitted repeated-weight shapes.

The four known cycles have `(p,q)` respectively

    (1,2), (1,1), (2,3), (7,11).

The first two are outside `p>=2` and `p<q<2p`. The `(2,3)` shape is already on the convergent list. The `(7,11)` shape is **not** on that list and must be included separately. All six legal overlapping paths at `(7,11)` have actual content `1` against `|D|=139`; none gives the `-17` cycle. This explicitly closes the exceptional known-cycle branch.

## 4. Independent exact continued fractions

The independent Python script uses rational arithmetic only, with 80 positive-series terms each for

    log2 = log((1+1/3)/(1-1/3)),
    log(3/2) = log((1+1/5)/(1-1/5)),
    log3 = log2 + log(3/2).

This differs from the producer's direct `log3` series at `t=1/2`. For `0<t<1`, after `K` terms, the omitted tail is positive and less than

    2 t^(2K+1) / ((2K+1)(1-t^2)).

The ratio enclosure is `[L3/U2,U3/L2]`. At every continued-fraction step, both endpoints have the same integer part. Subtracting it and reversing the interval endpoints under reciprocation therefore certifies each exact digit, without any floating-point guess.

The prefix through the first denominator at least 4096 is

    1,1,1,2,2,3,1,5,2,23.

The complete in-scope convergent shape list `(p,q)` is

    (2,3), (5,8), (12,19), (41,65),
    (53,84), (306,485), (665,1054).

The next convergent is `q/p=24727/15601`, already beyond the weight bound. Denominators of subsequent ordinary convergents increase, so there is no omitted denominator below 4096. Semiconvergents are not an additional possibility under the strict Legendre criterion.

## 5. Independent direct overlap census

The independent script generates the zero-intercept source directly from its floor differences, labels each occupied position, and enumerates all legal first cyclic edges and adjacent second cyclic edges. It explicitly rejects moving the same labelled one twice. It does **not** generate its test universe from the producer's residue formulas. The second movement direction is checked to agree with the first.

For each resulting word, the carry is computed by a reverse scan and independently by the positional power sum, using Python's arbitrary-precision integers. Integrality is tested by the actual gcd `gcd(B,|D|)`, not a resultant or modular surrogate.

| p | q | Plus paths | Minus paths | Total |
|---:|---:|---:|---:|---:|
| 2 | 3 | 1 | 1 | 2 |
| 5 | 8 | 2 | 2 | 4 |
| 7 | 11 | 3 | 3 | 6 |
| 12 | 19 | 5 | 5 | 10 |
| 41 | 65 | 17 | 17 | 34 |
| 53 | 84 | 22 | 22 | 44 |
| 306 | 485 | 127 | 127 | 254 |
| 665 | 1054 | 276 | 276 | 552 |
| **Total** | | **453** | **453** | **906** |

Eight paths cross the original cut. The only integral outputs are

- Plus at `(p,q)=(2,3)`: source `011`, edges `(2,1)`, output `101`, `B=7`, `D=-1`, orbit `(-7,-10,-5)`.
- Minus at `(p,q)=(2,3)`: source `011`, edges `(0,1)`, output `110`, `B=5`, `D=-1`, orbit `(-5,-7,-10)`.

After independent generation, the source position of the first moved one is used to recover the producer's residue parameter solely for record comparison. Every one of the 906 full records agrees, including word, actual edges, carry, clock, and actual content. The shared canonical trace is

    10721e4caf578f6475281c6a5085e27a955b1fe6c5c0b71d943fee50e437f505.

The independent audit performs 2,793 explicit checks and yields byte-identical normal and optimized Python receipts. Producer source hashes and normal/optimized receipt equality are also checked. The producer's candidate script was separately replayed with a byte-identical receipt.

## 6. Reproduction

From this audit directory:

    g++ -std=c++17 -O2 -Wall -Wextra -Wconversion -pedantic audit_signed_descent.cpp -o /tmp/audit_signed_descent
    /tmp/audit_signed_descent > descent.replay.txt
    cmp descent.receipt.txt descent.replay.txt
    g++ -std=c++17 -O0 -fsanitize=undefined -fno-sanitize-recover=all audit_signed_descent.cpp -o /tmp/audit_signed_descent_unoptimized
    /tmp/audit_signed_descent_unoptimized > descent.unoptimized.replay.txt
    cmp descent.receipt.txt descent.unoptimized.replay.txt
    python3 audit_small_cf_and_overlaps.py > candidates.replay.json
    python3 -O audit_small_cf_and_overlaps.py > candidates.optimized.replay.json
    cmp candidates.receipt.json candidates.replay.json
    cmp candidates.receipt.json candidates.optimized.replay.json

The deliberate cap-failure control is

    g++ -std=c++17 -O2 -DAUDIT_STEP_LIMIT=1 audit_signed_descent.cpp -o /tmp/audit_signed_descent_hostile
    ! /tmp/audit_signed_descent_hostile

A nonzero exit and `FAIL: step limit at start 3` are required in that last control. The normal verifier's success marker alone is insufficient if its process exited nonzero; all recorded successful runs completed with exit status zero.

## 7. Dependency boundary

The mathematical closure uses the already audited actual-content theorem plus its `p<4096` consequence, which in turn retains its stated Minkowski and Matveev dependencies. This audit does not silently replace those dependencies by a new global Collatz verification threshold.

The new dependency is only the fully executed finite signed descent through `2^24`, followed by elementary cycle products, Legendre's criterion, rational continued fractions, and 906 exact actual-word checks. Every computation and the inequality chain needed for this reduction has been independently checked here.


## 8. Broader overlapping-path corollary and final proof review

The [final proof](../OVERLAPPING_TWO_SWAP_CLASSIFICATION.md), including its Section 6 extension, has been reviewed. The broader corollary is justified, not merely extrapolated from the primitive `k=1` theorem.

1. Two different overlapping edges have a three-position union when `q>=3`. Inspecting actual bit exchanges leaves only a single labelled one moving twice or a two-distinct-one chain. At `q=3` this remains a valid three-position cyclic block. At `q=2` both directed edge labels describe the same physical transposition; any second actual swap is undo.
2. A mechanical single-one chain starts with `100` or `001`, so a gap between consecutive cyclic ones has length at least three. Mechanical gap lengths are `floor(q/p)` and `ceil(q/p)`, hence `q>2p`. A positive integer cycle has `2^q=product(3+1/a_i)<=4^p`, so `q<=2p`. A negative integer branch fixed point has `D<0`, so `q<p log_2(3)<2p`. Both contradict the gap condition.
3. A two-distinct-one chain starts with `110` or `011`, so it contains a gap of length one. The same mechanical gap law gives `q<2p`. Every primitive two-one case is therefore inside the newly classified `k=1` band. The equality case `q=2p` is alternating and has no nonundo overlapping path.
4. For `d=gcd(p,q)>=2`, the source is the `d`-fold repetition of the reduced source. The inherited repetition factor `H>=7` divides its carry and clock and is coprime to 6. After rotating the cut outside the affected three-bit block, the local carry difference is a 2,3-unit times `3` or `5`. Rotation preserves the full gcd with `H`, so the same obstruction holds at the original cut. Thus the edited carry cannot retain all of `H`. The factor `5` is not presumed invertible modulo `H`; `H>5` is the sufficient point.
5. The exceptional mirrored edits when `d=2` have edge starts separated by `q_0>=2` on a circle of length `2q_0`. Their affected pairs are disjoint, even for `q_0=2`, so they are not counterexamples to the overlap corollary.
6. Both exceptional `p=2,q=3` outputs are rotations of the mechanical source `011`. Their cyclic adjacent-swap necklace distance from that source is zero, not two. The corollary therefore correctly says that no integral output has exact necklace distance two **via an overlapping path**; the proof does not assert that primitive disjoint paths are impossible.

The repetition input was checked directly against the [current proved repetition packet](../../collatz_two_swap_repetition_classification_20261020/TWO_EDIT_REPETITION_CLASSIFICATION.md), Sections 2--3, and its independent audit. No use is made of an unproved primitive disjoint-edit statement.

As additional finite hostile controls, `audit_overlapping_path_types.py` independently enumerates all 2,016 nonconstant mechanical shapes of lengths 2 through 64 and their 74,494 ordered overlapping actual paths. There are 44,706 physical-edge undo paths, 14,894 same-one nonundo paths, and 14,894 two-one nonundo paths; 11,676 nonundo paths have a nonprimitive source. The gap implications and all actual content tests pass; the only integral nonundo outputs are the two declared `p=2,q=3` words. It performs 104,286 checks and has byte-identical normal/optimized receipts. This last computation is a hostile control, not a replacement for the all-length proof.

Reproduce these additional controls with

    python3 audit_overlapping_path_types.py > path-types.replay.json
    python3 -O audit_overlapping_path_types.py > path-types.optimized.replay.json
    cmp path-types.receipt.json path-types.replay.json
    cmp path-types.receipt.json path-types.optimized.replay.json

**Final audit verdict:** the primitive theorem and the broader overlapping-path corollary may both be promoted to **PROVED**, with their stated inherited analytic dependencies and finite-exact computational dependencies retained. No substantive mathematical gap remains in the reviewed claim. Primitive disjoint edits and the general Collatz conjecture remain outside scope.
