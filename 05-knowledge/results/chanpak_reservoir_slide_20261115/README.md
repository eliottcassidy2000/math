# Unit slides, armed fans, and target-selected reservoirs

The exact task is to separate a **local repair** of an already
nearly correct target from a **placement mechanism** that makes
the target nearly correct. The two do not need to be the same
operation. In fact the fixed-core results below show why requiring
one binary tape to do both is too restrictive.

## 1. One core now gives consecutive windows

An [eight-vertex marked core](UNIT_SLIDE.md), of ordinary
Hamilton count 45, has a 12-bit transitive slide with exact count
step 2. It preserves the whole two-sided context response. This
replaces the earlier [two-core 145/39 carry](../chanpak_window_descent_20261115/GROWING_WINDOWS.md)
when the only job is to repair a count by twos. At height \(k\) it
provides \(k+1\) consecutive odd counts on \(11+12k\) vertices,
all strong; small common padding gives a window at every order
\(n\ge 11\). Independent literal DP gives the first pair
**11,273,423 and 11,273,425** on 23 vertices. The local decoder
returns the position itself, not just the residue class.

The [finite slide cocycle](SLIDE_COCYCLE.md) computes the gcd of
the actual integer translations at a fixed returning interface,
over **all** tile lengths. The earlier 13-profile component has
index 10; the new 45-profile component has index 2 and returns
of both signs. A pair of finite walks realizes any even net
local translation there, by concatenation with nonnegative
multiplicities. This does not keep the absolute baseline fixed
while the length changes.

A [binary delay versus active carry theorem](SLIDE_ENTROPY.md)
limits this approach. If a return component can transmit both
binary input choices freely forever, it is a pure power-of-two
delay in its two unsigned channels and its count displacement
vanishes. Every nonzero repair component loses some binary
addresses, so a fixed or finite bank of such guarded sliding
cores alone has sublinear-power count capacity. For the sampled
45-core unit component a certified upper exponent is \(5/12\)
up to the logarithmic factor from the positions. The structural
lesson is not to demand a full binary center address inside
the device meant only to repair the last few units.

## 2. Both signs can be armed, but repeat-stable trees do not regenerate

A literal [count-59 one-parent tree](TWO_DIRECTION_RETURN.md) has
two distinct arms that are ready to repeat, forward \(-1\) and
reverse \(+1\). Their **readiness** can survive arbitrary mixing,
yet either move changes the other arm's actual unit defect. A
count-5 tree even has both opposite orientations at one site,
and distinct forward unit sites. These are genuine positive
answers to simultaneous availability, not binary recursion.

For arbitrary strong seeds and arbitrary interleaving of
repeat-stable moves in both directions, no new eligible vertex
can be created, and each direction can switch its used root
at most once. Finite seed banks in this grammar produce only
\(O((\log X)^4)\) counts through \(X\). A separate
[raw-unit SCC reduction](RAW_UNIT_COMPONENTS.md) shows that for
a reducible tournament of count greater than one at least one
of the two directional unit menus is empty, even without the
repeat-stability test. These theorems do not rule out non-leaf
preparation, an unstable unit move, or rebuilding the source
from the target.

## 3. Choose an exposed reservoir, not independent coordinates

[Wrapping one internal apex](RESERVOIR_CHARGE.md) makes a physical
reservoir whose ordinary count \(g\) varies while its exposed
start data stay fixed. In the appropriate mirror gluing the
output is an exact **shifted product** rather than a loose
congruence. All-marked tags change the shift through their
two-cover count; two literal six-vertex tags have the same 37
Hamilton paths but 352 and 376 two-covers. This is a small
example of changing the placement response without changing the
multiplier. It gives both a strong prime-valued example (4409)
and a strong \(7p\)-valued example (133). These are existence
examples, not a prime distribution theorem.

The explicit 37-count reservoir solves the residue debt for
every unit class modulo 37 at each unit-slide height not
divisible by 3, and then imposes the real interval test.
For example the prescribed count **11,721,300,267** decodes to
an actual reservoir and slide position. Its adjacent odd
control fails in that fixed bank. The [two-bank norm-shell
decoder](NORM_SHELL.md) explains another aspect of the geometry:
for symmetric reservoirs the target determines the sum of the
two slide positions by rounding, and a square discriminant
decides whether the two positions exist. This is a useful
bridge to square-sum and Berggren-style inversion: lowering or
selecting the enclosing shell is not enough unless the residual
has its exact bounded parameterization.

## Proof direction

The next theorem needed for the prime and \(7p\) least-missing
lanes is **reservoir supply**. Given the target, select a real
multiplier and marked profile that both pays the residue debt
and leaves an integer quotient inside the unit window. The
unit window removes all gaps inside that quotient interval;
it does not magically center it. Shifted products give an
exact alternative decoder, but a shifted prime or a failed
square test can rule out a chosen charge. The useful next
search is for independently adjustable charges at the same
multiplier, or for new selective contacts that widen the
readout without turning all free bits into a delay line.
Chan–Pak remains open.

The notes distinguish literal full subset-DP checks from
all-length proofs using exact ordinal transfer or finite graph
closure. Each replay has a frozen text receipt.
