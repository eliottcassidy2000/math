# Five-letter failure renewal without changing the marked source

**PROVED:** exact source-cylinder and affine-moment recurrences for consecutive
expanding returns of the total controller in [INTEGER_RETURNS](INTEGER_RETURNS.md).
**OPEN:** a pointwise descent certificate for every positive integer.

## 1. A finite failure alphabet with exact phases

Let `B` consist of the five words in (6) of INTEGER_RETURNS. They are
prefix-free: the queue identifies the end of each block. A concatenation
`w_1...w_k`, each `w_i in B`, remains the actual chronological valuation word
of one source. Its affine data compose by

\[
 J_{uv}=J_u+J_v,\quad D_{uv}=D_u+D_v,\quad
 C_{uv}=3^{J_v}C_u+2^{D_u}C_v.                         \tag{7}
\]

The final endpoint and guard must be evaluated with this total `C,J,D`, not by
substituting independent residue representatives at the intermediate returns.
As a port-preserving implementation, retain `(original source,J,D,C,current
endpoint,occupied external pins)`; this calculation creates no external port
identifications and authorizes no reuse of a pin.

## 2. Three exact masses

On odd 2-adic Haar sources, the cylinder of an exact valuation word has mass
`2^(-D)`. For one failure block set

\[
 p=\sum_{w\in B}2^{-D_w}=\frac{1059}{2048},\qquad
 s=\sum_{w\in B}\frac{3^{J_w}}{4^{D_w}}
   =\frac{3235233}{4194304},\qquad
 b=\sum_{w\in B}\frac{C_w}{4^{D_w}}
   =\frac{1270035}{4194304}.                            \tag{8}
\]

For `k` consecutive failure blocks, the probability, slope-weighted mass, and
carry-weighted mass are respectively

\[
 p^k,\qquad s^k,\qquad
 b_k=b\frac{s^k-p^k}{s-p}.                              \tag{9}
\]

Proof: cylinder masses multiply under word concatenation, and (7) gives
`b_k=s b_(k-1)+b p^(k-1)`, with `b_0=0`. This is a sum over exact source
cylinders, not a probabilistic replacement of the marked integer. The checker
enumerates all `5^k` word sequences through `k=5`, recomputes their guards and
affine maps, and realizes each on an integer in its original source cylinder.

The incoming [exact failure transducer](../session97_overnight_arithmetic_carriers_20261029/controllers/FAILURE_TRANSDUCER.md)
independently gives the product mass `p^k` and the guarded 25 two-block charts.
The additional quantities here are the slope and carry moments `s^k,b_k`;
the unweighted mass law is not counted as a separate discovery.

## 3. Source refinement pays three failure blocks per retained bit

For any fixed odd cylinder `I_K` modulo `2^(K+1)`, and any nonnegative quantity
`X`, `E[X|I_K]<=2^K E[X]`. With `k=3K+h`, (9) gives

\[
 E[\Lambda; k\text{ failures}\mid I_K]
 \le(2s^3)^K s^h,
\qquad
 E[B; k\text{ failures}\mid I_K]
 \le\frac b{s-p}(2s^3)^K s^h.                          \tag{10}
\]

The exact rational `2s^3` is strictly below one. Each failure block has
`J in {1,4,7}` and `D in {1,6,11}`, so `J<=7k`, `D<=11k`, and the physical
clock is at most `18k` on this alphabet. We retain these clocks rather than
renaming `k` as time. The probability version uses `p` and is even smaller.

## 4. An actual exceptional point and the integer boundary

The infinite word `111...` is legal on the odd 2-adic integer `-1`; each
return fixes it. All positive sources `2^L-1` share its first `L-1` odd
valuations and hence arbitrarily long strings of failing one-step blocks.
Thus no bound on failure length independent of source height is possible.
This is a hostile control for both uniformity and sign. The universally
terminating queue proved in INTEGER_RETURNS is a different stopping problem.

More generally a periodic concatenation of failure words has slope greater
than one and positive real carry; its affine fixed point is negative rational.
This excludes periodic all-failure histories for positive integers. It does
not exclude a nonperiodic infinite failure history.

The newer [all-period arithmetic certificate](../session97_overnight_arithmetic_carriers_20261029/controllers/PERIODIC_EXCEPTION_CERTIFICATE.md)
is stronger: fixed points of compositions lie in the convex hull of the five
individual fixed points, `[-16177/139,-1]`. Exact source guards on its 58 odd
integers leave only the periodic integer `-1`; the only integer sources with
eventually periodic infinite failure histories are `-1,-3`. Thus a positive
exception must have an aperiodic, strictly increasing, unbounded return
history. This is a genuine arithmetic exclusion, with the aperiodic case
still open.

## 5. The remaining sufficient targets

Either of the following would be substantial progress; neither is supplied
by (10):

1. A positive-integer height bound for the first escape from this five-letter
   failure alphabet, with an endpoint debt that can be paid across later
   descending returns. Escape alone need not undo earlier growth.
2. A well-founded potential on the marked source and current return state that
   pays the factor four at every queue restart, perhaps with a finite memory
   of the five failure phases and their actual accumulated carry.

Since the singleton marked integer has Haar measure zero, taking `K` to
infinity in (10), while allowing the tested horizon `3K+h` to grow, supplies
neither assertion. The new result closes integer **entry** into the auxiliary
return structure and identifies every expanding branch; it does not close
the all-integer **descent** gap.
