# Session227: physical constructions toward Chan–Pak

**COMPLETE within stated scopes; independently audited. Chan–Pak completeness remains OPEN.**

The target is actual tournament realization of every positive odd Hamilton-path count except 7 and 21. This session obtains the following scoped all-order results. It does not prove the target or claim a new unrestricted finite coverage cutoff.

## Results

1. **The smallest escape from O/Z/J state histories.** A rooted graft with at least two gadget vertices leaves all three terminal faces exactly when both input roots are nonsinks. A source-rooted edge already suffices, despite rank three. Its image plane has an integer inverse, and iteration gives explicit strong cubic-recurrence families and root-start multiplication. This separates full states, not scalar spectra. [Proof](boundary_escape/PROOF.md), [audit](escape_audit/AUDIT.md).

2. **Exact physical-state hulls, with a strict realization gap.** The closed cone has four facets; the affine integer lattice has index eight. Their intersection still contains the nonphysical point (1,0,6,1). Invertible grafts transport these exact hulls, and every order at least four has strong full-rank gadgets of both deletion-ratio signs. The exact scalar affine lattice is H+gcd(2H,Q)Z. Graph-image recognition needs a marked module: two explicit strong six-vertex graphs share state (37,7,43,9) and a physical numerical inverse, but only one is a graft image. [Proof](graft_images/PROOF.md), [audit](image_audit/AUDIT.md).

3. **Multi-root constructive calculus and sharp limits.** A root subset gives an exact positive path-cover expansion of sharp degree min(|R|,|A−R|)+1, an SCC strongness criterion, and the sharp bound H≥2^|B|+1. Interface size r needs exactly r+1 rational linear gadget coordinates. Fixed k-type block menus give only O((log X)^k) values; even arbitrary outer ordered sums give exp(O((log X)^(k/(k+1)))) values. Recursive adaptive grafting is outside that bound. Explicit strong families include −5·3^m+4·12^m+2·33^m, and binomial(2m,m) nonisomorphic graphs of order 7m+5 with the same entire rooted state and count 3·15^m+12·384^m. [Proof](multiroot/PROOF.md), [audit](construction_audit/AUDIT.md).

4. **Physical scalar switches and their propagation law.** Two literal word pairs change D while keeping E,U,d fixed. Even arbitrary mixtures of all four words have at most two scalar values in any final sidecar fiber, at every positive iteration depth. The D-column of any common graft-containing context is itself a physical rooted state. Consequently an existing pure-D gap is multiplied by an odd positive Hamilton count, preserving its exact 2-adic valuation; its U and d differences cannot vanish after the graft. This limits one proposed digit-building mechanism, not all adaptive constructions. [Proof](adaptive_arithmetic/PROOF.md), [audit](arithmetic_audit/AUDIT.md).

5. **The bounded-intercept obstruction survives finite graft menus.** An exact optimal ratio κ_B controls U_after/D_after for every nonempty rooted gadget. Splitting off the final O-run extends the logarithmic bounded-intercept theorem and the inherited residue-gap line obstruction to any fixed finite gadget menu and arbitrary finite physical seed bank. The common missing odd residue remains essential; no unrestricted finite-line obstruction is claimed. [Proof](menu_boundary/PROOF.md), [audit](menu_audit/AUDIT.md).

The exact cone does not characterize physical integer states. The graph collision does not rule out state-based existence proofs. Fixed-menu sparsity does not bound recursive adaptive source construction. The scalar-switch theorem concerns a common continuation, not branch-dependent choices. These distinctions are load-bearing.

## Verification and inheritance

[Verification scope](VERIFICATION.md) separates written proofs from finite controls. The portable [aggregate checker](verify_all.py) replays each producer and independent checker in both normal and optimized Python, comparing exact committed receipts. [REPLAY.json](REPLAY.json) records the run count and source bindings; `sha256sum -c SHA256SUMS` verifies final packet bytes. No external data or packages are required.

The work continues [Session224's birth atlas](../session224_first_occurrence_20261107/README.md), [Session225's graft law](../session225_graft_recovery_20261108/PROOF.md), and [Session226's inverse, cone bounds, and line obstruction](../session226_structural_continuation_20261108/README.md). Individual endpoint inequalities, the O/Z/J source grammar, and general path-cover substitution are explicitly credited rather than reclaimed. No literature-priority assertion is made.

Anchor: physical predecessor constructions. Niche: multi-root interfaces and retained cover coordinates. Wildcard: exact graph-image and scalar-sensitivity failures. The concept board was literal tournaments, four-state profiles, graph modules, positive cones, root-start complexity, and height-sensitive scalar existence.

[Remaining proof obligations and tested stopping reasons](NEXT_STEPS.md) identify what an actual all-odd existence proof still needs. The session's short transfer search is FINITE-EXACT evidence for its stated universe, not a proof against longer radix constructions.
