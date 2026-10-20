# Inverse trace lattice and finite-sheet selector

20 October 2026. **Scoped local proofs; independent review pending.**

The [proof](TRACE_HYPERBOLIC_PACKET.md) advances the previous inverse-trace frontier in three concrete ways:

- Maximal cancellation of m−1 trace moments forces every Smith exponent of an (m−1)-dimensional monomial trace lattice to equal −a, where a is the first trace pole order.
- In the local quartic A4 control, the rank-two trace form is hyperbolic and has a sharp mandatory lattice-repair index a. In all ranks over the stated complete DVR, the sharp index is ceil(a(m−1)/2).
- The coefficient ratio −c_(m+1)/c_m approximates the unique finite sheet to order a. Newton iteration improves that to order (2^(n+1)−1)a.

Planar generic degree four is already excluded by Domrina (2000), as confirmed by the primary publisher page and the repository’s earlier full-text source check. These quartic formulas are local-method controls; higher-degree sectors are the relevant open application.

The triangle/circle interpretation is exact finite quadrature: a three-root packet matches the first two holomorphic circle moments, and its cubic moment detects the pole. Its bilinear Gram form is not Hermitian-positive.

An integral-Gram hostile shows why a repaired trace lattice alone cannot prove properness. The remaining target is multiplication closure retaining the original source coordinates, or an actual common bound/gluing of finite-sheet selectors and the companion coordinate. No unconditional trace regularity, new global degree exclusion, or JC(2) resolution is claimed.

Reproduce:

    python3 verify_trace_hyperbolic_packet.py
    python3 -O verify_trace_hyperbolic_packet.py

The normal and optimized runs give identical stdout and JSON. The checker uses Newton identities and companion-matrix traces independently, along with truncated-series Newton iteration. It covers 72 root packets, 12 nontrivial quartics, 48 Newton residual checks, all-rank index controls, and explicit mixed-characteristic and integral-Gram hostiles.
