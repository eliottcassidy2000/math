---
id: THM-5000
title: "Coprime physical moving interfaces interlace dilated counts into consecutive odd windows"
status: PROVED, independently audited. Chan-Pak completeness and arbitrary target placement remain OPEN.
session: Session233-composition-20261109
source: 05-knowledge/results/session233_composition_20261109/graft/PROOF.md
audit: 05-knowledge/results/session233_composition_20261109/graft_audit/AUDIT.md
---

# Coprime interface composition

[Proof and certificates](../../05-knowledge/results/session233_composition_20261109/graft/PROOF.md) · [Independent audit](../../05-knowledge/results/session233_composition_20261109/graft_audit/AUDIT.md).

If gcd(a,b)=1, K>=b-1, L>=a-1 and C=(a-1)(b-1), then every integer in
[C,aK+bL-C] has the form ai+bj with 0<=i<=K and 0<=j<=L. The proof gives
an explicit bounded residue selector; gcd alone is not the claim.

Compose the inherited Session231 h=37 and h=45 full-interface moving-core
families. Their pure-T variations multiply by the opposite deletion count.
After reversing the second position index, the scalar rectangle is exactly
B_(K,L)+2(45i+37j), for 0<=i<=K,0<=j<=L. All witnesses are actual strong
tournaments, not abstract admissible state vectors.

Consequently, for every K>=36 and L>=44 one full physical interface supports
45K+37L-3167 consecutive odd Hamilton-path counts on 1188K+12L+23 vertices.
For K=36 this is 37L-1547 counts on 12L+42791 vertices. The construction
preserves every interface coordinate except T and remains saturated.

This composes distinct coprime cores; repeated copies of one core retain the
old divisor obstruction. The baselines grow exponentially and the widths
linearly. No arbitrary-target selector, interval-overlap theorem, unrestricted
coverage cutoff, or global optimality is asserted. Dependencies include
THM-4960-twelve-bit-moving-core-consecutive-odd-intervals.md and the audited
Session231 UNIT_INTERVALS.md, not automorphism omissions as completeness.

## Stronger short-loop composition

The later [weighted proof](../../05-knowledge/results/session233_composition_20261109/graft/WEIGHTED_PROOF.md),
with [independent supplement audit](../../05-knowledge/results/session233_composition_20261109/graft_audit/SUPPLEMENT_AUDIT.md),
reuses short loops whose individual scalar steps are not two.

- An eight-vertex h13 deletion core has a12-bit full-interface loop with T drift−10.
- An eight-vertex h29 deletion core has a28-bit loop with T drift−6.
- The inherited h45 core has a12-bit loop with T drift+2.

Compose protected families with respective movement ranges i=0,...,12,
j=0,...,28, and k=0,...,L. Reversing the final index, the exact scalar offsets
from the explicit baseline B_L are

    2(6525i+1755j+377t), 0<=t<=L.

The377 fixed representatives6525i+1755j cover every residue modulo377=13·29
exactly once. Their maximum is127440. An explicit CRT selector therefore gives
all odd counts in [B_L+254128, B_L+2(377L+376)].

**For every L>=337, one physical full interface supports377L−126687 consecutive
odd Hamilton-path counts on12L+962 strong-tournament vertices.**

At L=1000 this gives250313 consecutive odd counts on12962 vertices. The
asymptotic width/vertex ratio is377/12, versus1/12 for the inherited single
h45 core. This377-fold comparison is between these constructions, not a
claim of global optimality. The endpoints are explicit very large integers;
the L1000 lower endpoint has12966 bits.

Pair constructions and a nested alternative are also proved in the supplement.
The interval baseline still grows exponentially. The precise missing step is
an arbitrary-target placement or overlap theorem; scalar gcd, local endpoint
universality and arbitrarily long intervals do not provide it.
