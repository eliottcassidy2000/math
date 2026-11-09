---
id: THM-4960
title: "Unbounded consecutive odd Hamilton-path counts in one physical interface: k+1 values on 12k+13 strong-tournament vertices"
status: >
  PROVED, independently audited: for every integer k >= 0, one physical
  deletion interface supports all k+1 consecutive odd counts D_k-2j,
  0 <= j <= k, where D_k=(245051264*4096^k+1602510*k+1118023)/24843.
  These full-interface witnesses are strong and have 12k+13 vertices. A later
  independently audited deletion refinement gives strong scalar interval
  witnesses on 12k+9 vertices, but not a common full interface. The proof uses an exact
  twelve-bit word identity, not finite-k extrapolation. The core is saturated,
  with h=45 and gcd(h,q)=45. Unrestricted Chan-Pak completeness and arbitrary
  target placement remain OPEN. No literature-priority or globally shortest
  interval-constructor claim is made.
session: Session231-overnight-20261109
source: 05-knowledge/results/session231_overnight_completion_20261109/fibers/SHORT_INTERVALS.md
verification:
  - 05-knowledge/results/session231_overnight_completion_20261109/fibers/verify_short12_unit.py
  - 05-knowledge/results/session231_overnight_completion_20261109/fibers/verify_short12_independent.py
  - 05-knowledge/results/session231_overnight_completion_20261109/fiber_audit/short12_audit.py
  - 05-knowledge/results/session231_overnight_completion_20261109/fiber_audit/SHORT12_THEOREM_AUDIT.md
---

# Twelve-bit moving-core consecutive-odd intervals

## Inheritance and precise advance

Session229 `dynamic/PROOF.md` supplies the seven-coordinate ordinal law,
endpoint-sensitive context law, and fixed-width physical switches. Session231
uses those inherited identities to obtain **unbounded consecutive width** by
moving one nontrivial core rather than repeating it. Earlier fixed-width
small-core censuses and growing product-spacing obstructions do not imply
this result. The all-k proof below does not depend on any complete census.

## Construction and proof

A marked deletion tournament has coordinates (h,q,A,G,T,S,E): h counts its
Hamilton paths; q counts ordered spanning two-path covers, allowing empty
paths; A/G constrain the second path's start / first path's end to the marked
cut; T imposes both constraints; S/E count marked starts/ends of full paths.
A singleton marked by b has state (1,2,b,b,0,b,b).

For nonempty ordinal blocks, write their states as (h,q,A,G,T,S,E) and
(j,r,B,F,U,V,W). Direct path decomposition gives

    (hj, qr, Ar+hB, Fq+jG, AF+Tj+hU+EV, Sj, hW).

Use the eight-vertex core with outgoing-neighbor masks

    C=[4,1,18,7,75,31,175,63], cut=124,
    state(C)=(45,1440,878,577,332,25,14),
    X=101000111000, Y=111011100110.

Let L_(k,j)=s1 => s0 => X^j => C => Y^(k-j) => s1. A word denotes
singletons in forward ordinal order, with its bits recording the cut. Close
L by vertices v,r with L=>v=>r and r beating exactly its cut vertices.
The path count is 2A+h-T. The selected first vertex is reached by r and
reaches all later deletion vertices; every deletion vertex reaches r through
v. Thus the closure is strong.

The literal core state and ordinal law give

    state(s0=>C=>Y=>s1)
      =(45,23592960,7535881,16564034,5290654,0,45),
    state(s0=>X=>C=>s1)
      =(45,23592960,7535881,16564034,5290656,0,45).

These differ only by 2 in T. Each move has a left boundary ending in 0 and
right boundary beginning in 1, since X ends in 0 and Y begins in 1. The
endpoint-sensitive context law therefore preserves this exact difference
at every adjacent core position. It follows that T_(k,j)=T_(k,0)+2j and the
closure counts are D_k-2j. There are 12k+8+3+2=12k+13 vertices.

For completeness, set J_k=(4096^k-1)/4095. The common interface is

    h=45, q=11520*4096^k,
    A=7516*4096^k+343260*J_k+45,
    G=2353+9532800*J_k+5760*4096^k,
    S=E=45.

The source proof derives the baseline by exact ordinal composition:

    D_0=9909,
    D_(k+1)=4096*D_k-264150*k-184225,
    D_k=(245051264*4096^k+1602510*k+1118023)/24843.

The recurrence proves integrality and oddness; the graph construction proves
positivity. The first two intervals are {9909} and {40403037,40403039}.

## Validation and boundary

Independent literal permutation-and-split counts recover the core. Separate
seven-coordinate, six-state path-cover, and standalone audit implementations
agree. The audit directly counts the 13-vertex base graph and checks all
positions through k=8 plus k=20,100. These computations validate the finite
identity; the context argument is what proves every k.

Arbitrary common singleton contexts retain the step two. General tournament
contexts multiply it by their Hamilton counts. The deletion core remains
saturated; this is not an unsaturated-fiber result. Its h is fixed rather than
a growing product of repeated nontrivial cores.

The unpadded baselines grow exponentially while width is only 2k. No overlap
or all-target selection theorem is supplied. Chan-Pak completeness remains
open; see the Session231 packet for the exact restricted-language selector
and the obstruction to merely iterating grafts at the same closing root.


## Audited smaller scalar witnesses

The later Session231 `rerooted/PROOF.md`, independently verified in
`rerooted_audit/AUDIT.md`, strips the four protective singleton/closure
vertices. Take B_(k,j)=X^j => C => Y^(k-j) and add just the marked apex r.
For all k>=0, 0<=j<=k, these graphs are strong on 12k+9 vertices and have
counts Z_k-2j, where

    Z_k=(13617152*4096^k+1602510*k+717259)/24843,
    Z_0=577,
    Z_(k+1)=4096*Z_k-264150*k-118165.

The first three upper endpoints are 577, 2245227 and 9196067477. The
shared bare A coordinate is A_B=878*4096^k+171630*J_k. Its rooted state is

    (D,E,U,d)=(Z_k-2j,45-e_j,A_B+45,45),
    e_j=0 for j<k, and e_k=14.

Thus this smaller family is not a common full rooted-state fiber: the last
witness changes E. The original 12k+13 family retains that stronger interface
property. Both conclusions are all-order, and neither selects arbitrary
targets. Strongness of the bare graphs uses the actual mixed cut and strong
core C at an empty-word endpoint, and X's first selected / Y's last
unselected bits otherwise. The independent audit checks the finite word
identity and directly counts the physical 9- and 21-vertex bare controls.
