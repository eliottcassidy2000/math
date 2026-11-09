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
