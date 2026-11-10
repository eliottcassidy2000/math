---
id: THM-5108
title: Five-prime physical descent and exact join-bank cost
status: PROVED all exponents in a specified smooth sector and fixed graph bank; independently audited. Sector closure is an application of inherited atom reduction, not new unrestricted coverage.
source: 05-knowledge/results/session242_congruent_descent_20261110/descent/PROOF.md
session: Session242, 2026-11-10
depends_on:
  - THM-1862-order-join-reduction-principle
  - THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness
---

# THM-5108: a real failure-preserving descent with exact physical cost

Fix the ten explicitly certified graph seeds of counts
3,5,31,41,35,49,63,217,287,343 and respective orders
3,4,7,6,7,8,8,10,10,11 in the linked proof. Every

    N=3^a 5^b 31^d 41^e 7^c, a,b,d,e,c≥0, N∉{7,21},

is realized by an ordinal sum of these graphs. This is a concrete application
of the earlier prime/7p atom mechanism, with full physical witnesses for the
additional suppliers 217 and 287, not a new global coverage theorem.

For c=0 remove ordinary prime blocks. For even c≥2 remove 49; for odd c≥3
remove 343. At c=1 use a 35,217,287 block when its ordinary prime is present;
otherwise a≥2 and a 63 block is available. Every quotient stays admissible
and strictly decreases. Ordinal join lifts *every* smaller count witness,
without a special rooted-state condition. Thus Fail(N) implies Fail(N/b)
for the chosen block b, contradicting the witnessed terminal target 1.
Unlike vertex deletion, this descends a genuinely hypothetical missing target.

Let C(a,b,d,e,c) be minimum total order among sums of these specified graphs,
using the empty sum at 1. Set F(negative)=F(1)=∞ and

    F(3q)=11q,
    F(3q+2)=11q+8,
    F(3q+1)=11q+5 for q≥1.

Then exactly

    C=3a+4b+7d+6e + min [2y+3x+3z+4w+F(c−y−x−z−w)],

where 0≤y≤floor(a/2), 0≤x≤b, 0≤z≤d, 0≤w≤e. These are the
multiplicities of blocks 63,35,217,287. The remaining seven powers use
49/343, whose cost F minimizes 8v+11t subject to 2v+3t=c.
For K=floor(a/2)+b+d+e, every candidate residual is at least 2 when c≥K+2,
so the exact eventual period is

    C(a,b,d,e,c+3)=C(a,b,d,e,c)+11, c≥K+2.

This optimizes the specified bank, not all tournament realizations. The
threshold is sufficient, not asserted sharp. Sound descent need not optimize
cost: at 7^6 three 49 blocks cost 24 vertices, while two 343 blocks cost 22.
Blind arithmetic division is unsafe: 35/5=7 and 63/3=21 leave the admissible
set. Removing the 217 supplier makes 217 impossible inside the remaining bank.

[Full proof, graph masks, scope and reproduction](../../05-knowledge/results/session242_congruent_descent_20261110/descent/PROOF.md),
[independent audit](../../05-knowledge/results/session242_congruent_descent_20261110/audit/REVIEW.md).
Ten seed graphs, eight actual joins, 1,053 cost comparisons and 2,560 period
controls are checked exactly; the unbounded conclusions are proved separately.
The prior atom reduction, smooth mechanisms and attainable values are credited.
No elliptic-to-tournament map, strong-output closure or Chan–Pak proof follows.
