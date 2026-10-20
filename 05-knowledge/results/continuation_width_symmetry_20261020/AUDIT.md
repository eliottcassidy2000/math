# Proof audit and reproduction

**Status:** same-agent proof review with separate exact representations.
This is not an independent-agent review. All logical claims are scoped in
README.md; this packet is not a numbered canon promotion.

## Logical pressure points checked

1. **Descending lower bound.** Histories use lengths 3p above the cut, probes
   use primes p below it, and horizon 4 sum(p) admits every tested continuation.
   The response for p includes the common cycle count. That count must be
   retained: the indicator alone is not the formula. The fixed-degree
   rank argument also fixes cardinality, so its extra `2^(r^2)` is a
   constant. Positive historical/continuation factors do not change rank.
2. **Largest-prime upper bound.** The small kernel for an old cycle is its
   entire p-smooth part, not its radical. Powers must survive. Newly
   processed cycles in the current p-batch require a second colour and
   retain their full cofactor a/p. The argument covers intermediate batch
   positions. Actual committed degree remains a separate coordinate.
3. **Uniform analytic estimate.** The Rankin parameter varies with N but
   stays in [1/4,1/2), where every stated constant is uniform. The stronger
   bound uses classical `pi(t)=O(t/log t)`; the weaker log-log version
   needs only the integer Euler-product majorant. Neither uses a conjectural
   estimate for smooth numbers.
4. **Arithmetic cost.** The state bound alone would not prove a bit-cost
   bound. The explicit common denominator N! and polynomial numerator-bit
   estimates supply the additional step, including signed clone removal.
5. **Group support.** Apply regularity to actual whole-group orbits.
   Taking a union of the individual support bounds would be invalid.
6. **Skeleton.** Mark full orbits, not just endpoints. The orbit-size charge
   d and external backward-edge count e satisfy different inequalities.
   Quiet triangle contiguity follows from constant cross-orbit outdegree
   and the displayed order, and implies its module property.
7. **Action and ports.** Triangle modules are actual individual G-orbits.
   Therefore they cannot be exchanged. Independent rotations split off;
   pinning a vertex fixes that triangle's phase. The exceptional group H
   retains any synchronization. The delta=6 example prevents assuming this
   product for every triangle orbit.
8. **Counting.** The fixed-D skeleton series is a surjective upper bound;
   it is not an exact quotient count. The zero, three and four slack
   normal forms are separately proved injective, using strong-component
   order and recoverable extreme-score vertices.
9. **Feedback and group order.** No odd-order classification theorem is
   used. Involutions are impossible, and the elementary odd-factorial
   induction supplies the estimate. Equality is checked at every step.
10. **Collatz.** Affine evaluation by itself does not authenticate the word.
    The actual-source tests check exact valuations and terminal values.
    Refinement means splitting quiet runs without adding steps or changing
    the integer clock. No claim covers arbitrary source-changing refinement.

## Main exact run

`verify.py --bound 100 --large-bound 200` performs **360,949 always-active
checks**. Its normal and optimized Python outputs are byte-identical.

- All 33,868 labelled tournaments on 0 through 6 vertices. Automorphisms
  are enumerated among degree-preserving permutations, then checked on
  every pair. Minimum feedback uses subset DP; all orders independently
  verify it through five vertices.
- Common support, orbit cost, skeleton size, complete labelled arc
  reconstruction, the zero-slack classification, and its fixed-point-free
  mass. Larger ordered-component controls are also included.
- Six explicit prime-subset banks, with exact positive LDL pivots on
  common-mass/cardinality cells. Smooth-kernel gcd identities are checked
  at every prime batch for odd horizons through 83.
- Exact A000568 computations through 30,60,100 using both prime and
  descending orders; direct partition sums through 30. Signed core and
  fixed-point-free core series agree under both orders through 24.
- Prime-order extension through 200, agreeing with all earlier coefficients
  through 100. Its final 5,616-digit coefficient agrees with the frozen
  preceding packet's descending-order A200. This last comparison does not
  independently check every intermediate coefficient 101 through 199.
- All 21,845 exponent words of lengths 0 through 7 over {1,2,3,4}.
  There are 24,024 authenticated prefixes: the first 24 odd steps from
  every positive odd source at most 2001, including continuation at 1.
  Of the 1,040 supplied block shapes, 926 meet the sufficient contraction
  threshold. The other 114 are not declared non-descending.
- Thirty actual-source quiet-clock hostiles and 32 rotation-support
  hostiles with two binary defects.

| Horizon | Descending peak | Prime-batch peak | Descending transitions | Prime-batch transitions |
|---:|---:|---:|---:|---:|
| 30 | 72 | 72 | 1,271 | 1,340 |
| 60 | 420 | 420 | 12,468 | 12,597 |
| 100 | 2,953 | 1,896 | 90,219 | 80,277 |
| 200 | 96,836 | 40,504 | 2,994,814 | 1,480,961 |

The two smallest cases explicitly show that the new order need not improve
every finite transition count. The all-order separation is asymptotic.
The order-200 descending statistics are inherited from the preceding
packet, not rerun for this table.

## Separate representations

`independent_controls.py` uses a smooth-kernel **multiset** DP through 60,
instead of the response-vector DP, and an explicit subset-zeta Gram
factorization instead of LDL. It independently tallies all small classes
by `sum_labelled |Aut|/n!`.

The whole-group product, full permutation cycle index, and every choice of
zero, one or two fixed ports are checked for every labelled tournament
through order six. The group-order equality boundary and the exact
delta=3/4 class series are checked in that census. Generated larger
delta=3 normal forms and the synchronized-triangle hostile test the next
scope boundary. The exact check count is printed in independent-receipt.json.

Normal and optimized outputs must agree byte for byte. This second
implementation remains work by the same author.

The completed separate run has **1,257,908 always-active checks**, including
75 generated first-slack examples. Its normal and optimized outputs are
byte-identical. Together the two frozen receipts record **1,618,857 checks**.

## Commands

From this directory, with Python 3.10+ and only the standard library:

```text
PYTHONINTMAXSTRDIGITS=0 python3 verify.py --bound 100 --large-bound 200 > /tmp/width-main.json
PYTHONINTMAXSTRDIGITS=0 python3 -O verify.py --bound 100 --large-bound 200 > /tmp/width-main-O.json
cmp /tmp/width-main.json /tmp/width-main-O.json
cmp /tmp/width-main.json receipt.json
PYTHONINTMAXSTRDIGITS=0 python3 independent_controls.py > /tmp/width-separate.json
PYTHONINTMAXSTRDIGITS=0 python3 -O independent_controls.py > /tmp/width-separate-O.json
cmp /tmp/width-separate.json /tmp/width-separate-O.json
cmp /tmp/width-separate.json independent-receipt.json
```

The inherited checker and its frozen receipt are identified by SHA256 in
the manifest. All new files except the manifest itself are hashed there.
