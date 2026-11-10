---
id: THM-5093
title: "Endpoint half-bound equality and binary pair energy"
status: "PROVED; independently audited. Chan–Pak completeness remains OPEN."
session: Session237, 2026-11-09
source: 05-knowledge/results/session237_chanpak_synthesis_20261109/README.md
proof: 05-knowledge/results/session237_chanpak_synthesis_20261109/endpoint_slack/PROOF.md
audit: 05-knowledge/results/session237_chanpak_synthesis_20261109/endpoint_slack_audit/AUDIT.md
---

# THM-5093: Endpoint half-bound equality and binary pair energy

For a nonsink root a with positive endpoint count e_a, H(T)=2e_a+1 iff
T is a transitive prefix followed by the transitive tournament R_n with its
extreme arc reversed, rooted at the penultimate vertex (n>=3; all roots are
equivalent when n=3). Thus e_a=2^(n-3). The positive-endpoint hypothesis
excludes the many trivial (H,e)=(1,0) equalities.

For the explicit apex-chain tournament in the proof, let m encode the apex's
outneighbors in reversed chain order. With the inherited digit-pair energy
Q(m), its root endpoint E and full count H satisfy

    E=m-Q(m),  H=1+2m-Q(m),  H-2E=1+Q(m).

The pair contribution is2 at adjacent one-bits and2^(d-1) at distance d>=2.
This gives a physical interpretation of the same arithmetic statistic used
in fixed-port modular states. It retains the endpoint and bit order, not the
general rooted response.

If f(N) is the least full H of a strong tournament with endpoint count N,
then f(N)>=2N+1, with equality iff N is a power of2. Existence of f(N) is
THM-5001; the half bound and insertion injection are THM-4876. Its audit D
had the equality classification only through order9; Session231's related
join-equality classification had an extra cover-defect condition.
No formula for f(N) off powers of2 or Chan–Pak completeness is asserted.

[Full proof](../../05-knowledge/results/session237_chanpak_synthesis_20261109/endpoint_slack/PROOF.md) · [Independent audit](../../05-knowledge/results/session237_chanpak_synthesis_20261109/endpoint_slack_audit/AUDIT.md).
