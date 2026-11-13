# Why a finite dyadic atlas cannot renew the exact ternary target

Let a ready chart, in its transitive parameter, have full count
`H_i(k)=a_i 2^k+b_i`, with positive integer a_i, integer b_i, and k≥0.
The assertions below concern only this scalar projection. In particular
success in the scalar projection is not a certificate for entry. The
fixed-clock-shift cycle obstruction was already recorded in
[the two-entry funnel transport](../../two_entry_funnels_20261103/TRANSPORT.md).
The extension here allows arbitrary clock reindexing, sporadic handoffs,
and mixed affine rules.

## Infinite handoffs are rigid, and their directed graph has no cycle

Suppose H_j(l)=3H_i(k)+2 for infinitely many distinct k, with l allowed
to depend on k and with even arbitrarily sparse physical subdomains.
Then l−k takes finitely many integer values: after dividing by 2^k,
a_j 2^(l−k)=3a_i+(3b_i+2−b_j)2^(−k), which is bounded away from zero
and infinity for sufficiently large k. Along one infinite subsequence
l−k=s is constant. The right-hand side is constant for at most one
distinct k unless b_j=3b_i+2. It follows that

    a_j 2^s = 3a_i,         b_j = 3b_i+2.

Conversely those two equalities give a scalar handoff for every k with
k+s≥0. They do not prove any physical edge. The odd part of a_j is
three times the odd part of a_i, even when s<0. Thus v_3(a_i) strictly
increases by one across every infinite handoff, so a finite bank has
no directed cycle of such edges. This applies even if a chart uses
unbounded gaps in the allowed parameter set. Each other ordered pair
of charts has at most finitely many accidental handoffs, by the
bounded-shift argument. There are finitely many pairs in a finite
bank.

An **entire indefinitely iterated** orbit cannot stay in a finite
bank either, without any Diophantine finiteness theorem. Indeed its
counts H_{t+1}=3H_t+2 strictly increase for positive H_0. It can use
only finitely many accidental edges (each has finitely many input
counts). After the last accidental edge each remaining edge raises
the v_3 of its chart's leading coefficient. This is impossible in a
finite bank. When a value belongs to several charts one may choose
arbitrary chart labels along the trajectory: every resulting edge is
still either one of the finitely many accidental equalities or
raises v_3.

A stronger but nonconstructive statement concerns even sparse
revisits with intervening families outside the bank. For fixed
positive H_0 write D=H_0+1 and H_t+1=D 3^t. For chart i,
`D 3^t - a_i 2^k=b_i+1`. If b_i+1=0, comparing 3-valuations fixes
t. Otherwise divide by b_i+1. The resulting two nonzero terms sum
to 1 and are rational S-units for a fixed finite S (primes dividing
2·3·D·a_i·|b_i+1|). The standard rational S-unit equation theorem
gives finitely many pairs; it is a cited input only for this
stronger conclusion. The elementary finite-bank proof does not use
it.

## Mixed exact affine handoffs

The same bounded-shift proof works for a fixed positive rational slope
q and rational offset d, replacing 3 and 2 by q and d. On an edge that
is realized at infinitely many scalar inputs, the equations are

    a_j 2^s = q a_i,       b_j = q b_i + d.

Multiplying along a cycle forces the product of slopes to be a power
of two, and the composed offset map fixes the initial intercept.
A finite bank cannot support an unbounded monotone sequence that keeps
using a slope divisible by 3 with no compensating inverse 3-valuation:
each ordered pair and each rule has only finitely many accidental
inputs, and the remaining steps propagate the 3-adic valuation of the
positive dyadic coefficient exactly. The rule bank is assumed finite,
and slopes with negative 3-valuation are expressly outside the latter
claim. This is a scalar obstruction, not an impossibility theorem for
a graph family with a genuine mixed dyadic/ternary clock.

## Payoff for the mod-45 question

The inherited concrete entry and exit are `45·2^k+7 → 135·2^k+23`.
They pass the scalar compatibility test: v_3 of their leading
coefficients changes from 2 to 3. A next all-tail ready chart must
have odd leading part 405 and offset 71; a dyadic reindexing cannot
change that odd part. The next one thereafter needs odd leading
part 1215 and offset 215. Thus mod 45 is the first symptom of a
graded obstruction. Even finding a new physical chart at the next
target will not turn a finite collection of single dyadic rays into
an all-time exact 3H+2 recursion.

An escape can instead use an unbounded core/compiler with genuine
ternary scale, use a chart with more than two observed modes, or
change the arithmetic return map (e.g. preserve a branch choice
without insisting on H ↦ 3H+2 at every handoff). The inherited
`uniform_return_metagraph_20261111/reset/ITERABLE_RETURN.md` is an
example of the last distinction: its actual renewable move is
4H−33, not a repeatable exact ternary-target handoff.

The comparison with the merged metagraph is exact at the level of
proof discipline, not an identification of the two systems. A chart
count is a projection of the marked state; a merged node is a
projection of a framed tiling. A handoff or gluing must commute with
the retained observations before one may iterate on the quotient.
In the dyadic case we can already rule out a whole projected
transition system, which rules out every physical lift of that
specific system. In the metagraph case the compatible lower cards
determine a physical lift, as proved in the companion note.
