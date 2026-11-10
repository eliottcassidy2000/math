---
id: THM-5091
title: "Minimal modular tail memory and exact joint seed images"
status: "PROVED; independently audited. Chan–Pak completeness remains OPEN."
session: Session237, 2026-11-09
source: 05-knowledge/results/session237_chanpak_synthesis_20261109/README.md
proof: 05-knowledge/results/session237_chanpak_synthesis_20261109/arithmetic_states/PROOF.md
audit: 05-knowledge/results/session237_chanpak_synthesis_20261109/arithmetic_states_audit/AUDIT.md
---

# THM-5091: Minimal modular tail memory and exact joint seed images

Let Q(0)=0, Q(2a)=Q(a), Q(2a+1)=Q(a)+a+(a mod2). For every M>=1,
(m mod2M,Q(m)/2 modM) attains every pair, arbitrarily far out. This retains
both coordinates of the inherited Session206 block-pumping construction.

For the physical rooted-edge tail H(m)=1+2m-Q(m), the normalized full-count
observer (H(m)-1)/2 modM has exactly 2M^2 states in its minimal deterministic
MSB-first output automaton for M>=2, and one state for M=1. Reachability is
proved by pumping, and every distinct state pair is distinguished by the
empty suffix, 0, or 1.

For physical seeds s_i=(h_i,R_i,e_i,d_i) sharing one tail parameter, the exact
joint residues (H_i(m)-1)/2 modM are

    union_(b=0,1) { A(b)+j R-z d : j,z in Z/M },
    A_i(b)=(h_i-1+(R_i-e_i)b)/2.

These are synchronized-parameter images, not a theorem about finite-bank
integer coverage, density, or independence between separate image sets.
The O/Z law is inherited from Sessions195–206 and THM-4888/5085, not new.

For the physical apex endpoint E(m)=m-Q(m) from THM-5093, the exact minimal
state count for E modM is1 if M=1, 2M^2 for odd M>=3, and M^2/2 for even M.
At even M the state is (m modM,Q modM), with Q restricted to even residues.
At odd M the parity lift m mod2M is necessary. This gives fourfold reduction
relative to the NORMALIZED full-count observer at even M.

For avoidance of ambiguity, actual H modM has2M^2 states at odd M>=3,
M^2/2 at even M>=4, and one state at M=1,2. It is not the same observable
as (H-1)/2 modM. The proof supplement is
[ENDPOINT_OBSERVER](../../05-knowledge/results/session237_chanpak_synthesis_20261109/arithmetic_states/ENDPOINT_OBSERVER.md),
independently audited in the
[endpoint audit](../../05-knowledge/results/session237_chanpak_synthesis_20261109/endpoint_slack_audit/AUDIT.md).

[Full proof](../../05-knowledge/results/session237_chanpak_synthesis_20261109/arithmetic_states/PROOF.md) · [Independent audit](../../05-knowledge/results/session237_chanpak_synthesis_20261109/arithmetic_states_audit/AUDIT.md).
