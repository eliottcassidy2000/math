# LRC14 Minkowski/Circuit/Ising/PDE Compiler

codex-2026-06-27-S264

The useful merge is not a bigger metaphor.  It is a proof compiler.

Minkowski handles the primal geometry-of-numbers side: relation lattice,
covolume, convex body, and successive minima.  Circuit complexity handles the
proof architecture side: how many gates, how much depth, which reductions are
uniform, and which missing sidecar gate makes a scalar quotient illegal.  The
Ising model handles transfer/free-energy and, through Lee-Yang, the root
geometry of a whole distribution.  PDE weak form language handles the dual
operator side: mass, stiffness, boundary data, and zero modes.  De Moivre's
quintic contributes the exact finite-depth fold:

```text
x = z - a/z
x^5 + 5*a*x^3 + 5*a^2*x = z^5 - a^5/z^5.
```

The S264 script verifies that identity as a Laurent-polynomial cancellation,
then ranks the proof carriers in a tournament.  The leading retention path is:

```text
finite_address_packet
-> endpoint_phi_gap_operator
-> pde_weak_form_stiffness
-> lee_yang_ising_transfer
-> minkowski_successive_minima
-> circuit_complexity_audit
-> de_moivre_quintic_fold
-> jacobi_theta_signed_tail
-> raw_scalar_p0
```

The ranking is not a theorem.  It is a guardrail against forgetting.  Raw
`p0`, raw volume, raw Ising energy, raw circuit size, or a pretty quintic fold
is not proof currency unless it exits through endpoint owners, boundary
conditions, low-height wall deletion, root strata, and finite-address or
observer-gluing packet data.

The PDE lookback was especially clarifying.  The repo already had Galerkin
stiffness/mass matrices in the tournament matrix atlas, exact endpoint-cover
circuit positivity in HYP-2108, the `Phi` gap functional in HYP-2112, the
normal-fan/Cech component-bound route in HYP-3101, and the `phi4`/Lee-Yang
root-curve route in HYP-3108/HYP-3109.  This means "PDE" should not be treated
as a decorative analogy; it is the language for keeping boundary conditions and
zero modes alive after discretization.

The next concrete move is small: pick a handful of HYP-2963/HYP-3107 residual
rows and annotate them with route type, mass/stiffness/boundary data if
available, relation-lattice low-height wall status, circuit-depth stage, and
finite-address exit.  If that ledger starts separating rows that scalar
measure or raw root counts merge, the compiler is earning its keep.
