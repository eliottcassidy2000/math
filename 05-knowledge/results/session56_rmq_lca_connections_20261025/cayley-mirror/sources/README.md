# Session31: tree moments and full-kernel reconstruction

Status: PROVED analytic identities and obstructions; exact rational numerical certificates; classical provenance explicitly recorded. No claim that the MRB defect is explained or evaluated.

Read `TREE_MOMENT_BRIDGE.md` for the proofs.

Key results:

- The tree series R(z)=−W₀(−z)/z is the Cauchy/Stieltjes transform of the probability density ρ(x(v))=v/π, where x(v)=exp(v cot v)sin(v)/v.
- Its entire admissible Mellin half-plane is explicit: M(s)=(s+1)^(s−1)/Γ(s+1), Re s>−1. This supplies strict Hankel positivity, log-convexity, coefficient bounds, and Gauss/Padé certificates.
- The complete MRB analytic kernel is reconstructed by w^(1/w)=w²/[Γ(1/w)M(1/w−1)] for Re w>0. This is a full-object nonlinear re-expression, not a q⁴ defect theorem.
- Fixed finite positive quadrature destroys the required negative-moment endpoint and makes the reconstructed MRB terms diverge. A continuous logarithmic endpoint model repairs this with an absolutely summable correction, whose leading error is −(log k)²/(2k²).
- At the actual hardest integer mirror input z=log3/3, an entirely rational 64-node certificate gives bracket width less than 3.258520×10^−8. At the separate rational benchmark z=1/3, 32 nodes give width less than 4.591×10^−18.

Reproduce:

    python verify_tree_moments.py
    python verify_mirror_endpoint.py
    sha256sum -c SHA256SUMS

The first script writes `verification.json`; the second writes `mirror_endpoint_verification.json`. The 64-node exact rational solve took about 266 seconds in this run. Saved `.out` files retain the concise completed-run summaries. The proof plus exact fractions certifies the enclosures. mpmath quadrature and Lambert-W values are only diagnostic cross-checks.

`SOURCES.json` records verified primary-source URLs, with no redistributed third-party full texts. `SHA256SUMS` hashes this packet's source and result files. No Git or repository mutation was performed in this lane.
