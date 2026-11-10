---
id: THM-5099
title: "Value-refined fixed-permutation SCC recurrence and the sharp strong seven-symmetry gap"
status: "PROVED from inherited cyclic profile and path-action theorems plus independently replayed exact eight-circulant input; separately audited. Chan-Pak completeness OPEN."
session: Session239, 2026-11-10
source: 05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/PROOF.md
audit: 05-knowledge/results/session239_sequence_engines_20261110/audit/strong_gap_audit.md
---

# THM-5099: retain path values, then exclude an infinite symmetry sector

[Full proof and exact inputs](../../05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/PROOF.md)
· [Recurrence audit](../../05-knowledge/results/session239_sequence_engines_20261110/audit/README.md)
· [Gap audit](../../05-knowledge/results/session239_sequence_engines_20261110/audit/strong_gap_audit.md).

For a permutation profile m, let F_m count its fixed tournaments by their
Hamiltonian-path values in the monoid algebra [a][b]=[ab]; let G_m count the
strong fixed tournaments. Then

    G_m = F_m - sum_(0<k<m) product_d binom(m_d,k_d) G_k F_(m-k).

The unique first SCC is a union of whole permutation cycles. Weighted
Burnside gives exact ordinary class histograms. For fixed target h, projecting
to divisors of h is a legitimate rank-tau(h) ring quotient, so nondivisor
states may be pruned. This does not remove the costly identity-profile input.

Let gamma_p=min_R H(C3[R,1,1]) over p-circulants, p an odd prime, and mu_p
the inherited cyclic minimum of THM-5094
(`cyclic-orbit-profile-cost-and-sharp-sevenfold-path-floor`). If a strong T has
an order-p automorphism, either |T|=p or

    H(T) >= min(gamma_p,mu_p^2+2p).

With at least two moving p-orbits, the strict cyclic profile bound applies.
With exactly one, it is a module; a triangle through its vertex in the strong
quotient lifts to C3[R,1,1]. Induced-tournament path counts are monotone under
vertex insertion. The inherited two-path-cover kernel gives

    H(C3[R,1,1])=3H(R)+F2(R).

All eight seven-circulants have F2=1638; six have H175 and two have H189.
The audit independently counted all9! orders for every inflated tournament.
Thus a strong tournament with an order-seven automorphism either has exactly
seven vertices and H in{175,189}, or

    H(T) >= 2163.

The latter threshold is sharp on nine vertices at C3[Q7,1,1]. It is not an
exhaustive equality classification. The small-prime analogues are15 for p3
and135 for p5. Existing full-group module costs from THM-4987 are credited in
the proof; the new reduction handles arbitrary order-seven actions.

A further [audited corollary](../../05-knowledge/results/session239_sequence_engines_20261110/audit/palette_and_asymmetry_audit.md) is that every strong realization of203 or217 is asymmetric.
The inherited cited Busch bound in THM-4744
(`interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve`)
forces order<=10. The automorphism order divides the path count; factor7 is
excluded by this gap and factors29 or31 by the embedding into S_n. No
existence assertion is included.

Strongness is essential: transitive singleton padding around Q7 retains
H175 at arbitrarily large order. The group must act on the actual tournament;
H/|Aut(T)| is a different statistic. This removes an infinite strong symmetry
sector but does not settle Chan-Pak or construct asymmetric realizers.
