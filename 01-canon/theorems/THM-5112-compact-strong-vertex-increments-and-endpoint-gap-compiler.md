---
id: THM-5112
title: "Compact strong one-vertex increments below five times the half-gap, and an endpoint-only gap compiler"
status: "PROVED, independently audited. Sharp cost within the declared constructor; global gap existence is inherited. Chan–Pak remains OPEN."
session: Session243-prime-gap-bridges-20261110
source: 05-knowledge/results/session243_prime_gap_bridges_20261110/constructions/GAP_UNIVERSALITY.md
audit: 05-knowledge/results/session243_prime_gap_bridges_20261110/construction_audit/AUDIT.md
---

# Exact increments with a physical deletion and a small absolute count

For every positive integer N there are strong tournaments L_N,R_N, with
R_N minus one designated vertex isomorphic to L_N, such that

    H(R_N)-H(L_N)=2N,
    2N+1 <= H(L_N) < H(R_N) <= 5N,
    |R_N| <= ceil(log2 N)+4.

For N>=2 the orders are exactly ceil(log2 N)+3 and ceil(log2 N)+4.
For N=1 they are3 and4, with counts3 and5. The5N bound is attained in this
constructor precisely at N=1 or N=2^k+1 for k>=1. No global optimum is claimed.

**Proof.** At rooted edge state (D,E,U,d)=(1,0,2,1), the inherited O/Z moves
append a sink or a fixed-root bridge, with

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D+U-E,d,2U-d,d).

Every U>=2 is obtained uniquely from2 using U->2U or U->2U-1; the inverse
strictly decreases, using U/2 when even and (U+1)/2 when odd. For N>=2
build response U=N and append O, obtaining (D,0,2N,1), with1<=D<=N because
root-splitting injects Hamilton paths into distinguished two-path covers.
For N=1 use the edge directly. The suffixes Z and OZ give D+2N and D+4N.
Both graphs are strong: the old root reaches all old vertices, all old
nonroot vertices reach the last bridge, and it reaches the root. Deleting
the penultimate newborn O vertex from the latter graph leaves the former.
The inverse word length gives the order bound. Tracking U-D proves equality
D=N exactly for a nonempty all-Z word, whose U is2^k+1.

## A universal five-letter endpoint comparison

Let J perform the same bridge as Z and mark the new vertex. Its state is
(D+U-E,D-E,D+U,D). For every physical rooted source,

    H(A.JOJOZ)=9D-3E+9U,
    H(A.JOOOZ)=9D-E+9U.

Thus their difference is exactly2E. With a strong source both final graphs
are strong and have the same order |A|+5. THM-5001's endpoint selector
supplies E=N, and increasing its near-transitive base order while preserving
those endpoint digits supplies unbounded distinct pairs.

This is a different compiler, not a new universal gap-existence theorem:
Session199 already gives infinitely many same-order strong single-edge pairs,
and Session221 gives an explicit logarithmic-order single-edge constructor.
The new compact one-vertex family and the endpoint-only word identity retain
different physical information.

The endpoint compiler's lower count is always divisible by3 and exceeds3
for these strong sources. It therefore supplies no prime pair at all.
Nor does the compact family place either count at an arbitrary target.
In particular it never yields a +2 extension from5 to the forbidden7.

[Full proof, equality boundary, exact inherited paths and replays](../../05-knowledge/results/session243_prime_gap_bridges_20261110/constructions/GAP_UNIVERSALITY.md).
Finite controls independently count the literal graphs and verify deletion,
strongness and the coefficient rows. Gap selection is not Chan–Pak coverage,
prime selection or a universal count-decreasing inverse.
