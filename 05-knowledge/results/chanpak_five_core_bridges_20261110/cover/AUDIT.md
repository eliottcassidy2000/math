# Independent cross-lane audit

2026-11-10. This audit separates proved reductions from their finite controls.
It does not assert all-height Chan–Pak coverage. It reads other lanes without
changing their files.

## Five-module interface and central switch

The rank decomposition in `../algebra/FIVE_MODULES.md` is valid over every
field. In particular, im(a) is contained in ker(b), the chosen W3,Wb,Wa,W0
are a direct decomposition, and lifts through b are independent modulo its
kernel. The five multiplicities follow and are nonnegative. The natural
matrix carrier is the length-three indecomposable, while the four smaller
ones annihilate the central element c=ab. This establishes the claimed
connection to the intrinsic five-module classification, without identifying
the five tournament cores with those modules.

Independent literal path enumeration for B95 masks12,20 gives respectively

    (h,d,L,i,o,s,e)=(351,95,1276,911,312,71,18),
                  (353,95,1276,911,312,71,18).

Coefficient-wise enumeration also finds identical full polynomial I(t),O(t)
channels, not merely equality at t=+/-1. The joint channel difference is

    J20(t)-J12(t)=-2t^2-6t^3+t^4+3t^5+2t^6.

Thus the claimed central difference holds polynomially for this physical
pair. Its +1 and -1 values are -2 and4, giving total count difference2.
This audit supplements the inherited Session229 endpoint propagation with
the current response translation; it does not claim a new unit-pair theorem.

## The all-additions readiness obstruction

`../alternatives/NOTE.md`, Proposition5, correctly invokes the strong
deletion condition. The inherited READY.md Section1 explicitly states that
an indegree-one ready root z with q>s has U-z strong. The fixed-start
injection therefore applies to S-u inside U-u even though S-u need not itself
be strong. The insertion rule is valid: after the last current vertex
beating an incoming outside vertex, that outside vertex beats the next one.
The outside-vertex schedule can depend only on the current vertex set, so
restricting the resulting path to the original vertices recovers its input.

The indegree-two case also closes: omitting the ready root leaves a selected
prefix source; omitting a port gives a one-apex six-chain with at most64
paths; retaining root and ports gives exactly47,49,39,55,59 at order7.
None contains a strong H73 copy with the stated order. The conclusion is
about readiness in strong supertournaments, as the note states.

## Predeclared bounded strong-extension audit

Independently reproduce the COMPLETE strong one-vertex extension frontiers
of the two named H141 sources through target427, without importing the
other lane's enumeration implementation. Use the even-path weighted cut
polynomial, independent of its two-path start/end convolution. Compare
every retained labelled path and literal child, including frontier closure;
verify the two-ear cost controls separately. The expected layer sizes are
25,4,0 and28,4,0, as reported by the originating lane. They are controls,
not assumptions used to prune the independent enumeration.

The all-order implication additionally needs the written strong-ear argument:
a larger strong tournament has either a one-vertex strong enlargement or
a cyclic two-vertex ear of S; the latter has H(S)+Q2(S)>= (|S|+2)H(S).
Here that is at least1410, already above427. This proof obligation is
separate from matching the finite frontiers.

**Audit outcome: PASS.** The final
[BIRTHS.md](../alternatives/BIRTHS.md) supplies that argument correctly.
With no mixed-contact outside vertex, outside vertices have uniform contact
with S. Strongness forces an edge from the class dominated by S into the
class dominating S, producing the specified cyclic two-ear. Its two
Hamilton-order cases give H+Q2 exactly. Ordinary path insertion ensures
that any induced subgraph of a target of count at most427 also has count
at most427. Since the ear would have count at least1410, every bounded
strong extension is reached through one-vertex strong extensions. The
empty computed frontier is therefore a closure certificate, not a fixed
order cutoff.

The last SCC reduction is also valid. A strong source lies in one strong
component C; H(C)>=141 and divides the total. Every proper divisor of425
is at most85, and every proper divisor of427 is at most61, so H(C) equals
the requested total in either case. Hence the strong-frontier obstruction
extends to arbitrary, possibly nonstrong, supertournaments at these two
specific totals. It must not be generalized to every total at most427:
for example a target with a proper divisor141 has a different SCC option.

[audit_frontier.py](audit_frontier.py) imports no repository implementation.
It uses weighted EVEN directed paths for the cut polynomial and a
first-vertex Hamilton recurrence. In contrast, the originating
[frontier.py](../alternatives/frontier.py) uses
[two_birth.py](../alternatives/two_birth.py)'s two-path start/end convolution.
All35,840 labelled cuts agree, including every retained child row and its
attachment history;126 literal controls pass. The independent frontiers
are exactly25,4,0 and28,4,0, with two-ear counts3237 and3069 and common
lower bound1410. See [frontier_audit.json](frontier_audit.json).

Normal and optimized outputs of this checker, the four cover proof checkers,
and the independent support-prefix checker are byte-identical; their JSON
objects match all retained receipts. The intentional frontier `--fail`
control exits nonzero. [lane_checks.json](lane_checks.json) records source
and output hashes. Classical module classification priority remains with
the cited algebra sources; inherited tournament laws remain credited to
their earlier packets.
