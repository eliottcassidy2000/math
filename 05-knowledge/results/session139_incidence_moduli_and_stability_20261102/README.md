# Session139: tensor parity, sharp stability and arithmetic incidence

**Final research frontier — independently accepted, 2026-11-02.**

## Main result: exact stability from140, not closure of the whole gap

For weighted A, set `V=sum Aij*Aji`, `Delta=H-det(I+A)-abs(det A)` and
`R=floor(n/2)`. The optimal coefficients are exactly

```text
kappa_(2R)   = (4/R)[1-(3/4)^R-(1/4)^R],
kappa_(2R+1) = (4/R)[1-(3/4)^R].
```

They are proved at **n=3..9 and every n>=140**. This session closes1360
previously unresolved orders140..1499. **10..139 (130 orders) remains OPEN.**
Positive-variance equality in the new range is maximal ordinal fair pairs,
with one singleton in odd order; n=3 instead has Delta=V for every weighted
triangle. Threshold140 is sufficient, not asserted optimal.

The actual proof chain is

**odd Pfaffian parity -> odd-coset tensor norm -> sharp matching product
margin -> retained nonlinear reserve near1 -> n140.**

The tensor lemma has an independently audited second proof using
Eisenstein integers, F4 and Bernstein coefficients. It controls the energy
of whole Pfaffian groups; it does not revive a false termwise lower bound.
[New cutoff](stability/PROOF_N140.md) · [Independent audit](independent_threshold140/AUDIT.md)
· [Final source pin](independent_threshold140/FINAL_STATUS_REVIEW.md)
· [Tensor/product proof](tensor_matching/PROOF.md)
· [Arithmetic second proof](stability/EISENSTEIN_BERNSTEIN_TENSOR_PROOF.md).

A [separate all-modulus tensor appendix](tensor_matching_moduli/APPENDIX.md)
proves the sharp residue-coset factor `distance(rho,MZ)^2`; it does not upgrade
the ordinary tournament Pfaffian application beyond its parity hypotheses.

## Concurrent proof routes and arithmetic precision

The concurrently published [shared-baseline n276 proof](../stability_priority_moduli_20261102/STABILITY_N276.md)
combines a sharper greedy degree bound with this packet's per-edge geometric-tail
mechanism. It and the independent [n314](../stability_748_guarded_encodings_20261102/STABILITY_314.md)
and [n466](../stability_priority_spectrum_20261102/STABILITY_466.md) routes remain
credited, separately scoped proofs; the new tensor-based n140 route is now the
strongest accepted range. Their packets and receipts are not overwritten.

The [concurrent guarded-arithmetic work](../stability_priority_moduli_20261102/FRACTRAN_PRECISION.md)
proves exact FRACTRAN moduli and the minimal reciprocal clock observer. Its
THM-4852 precision correction and provenance entries are retained unchanged.
[Reconciliation record](RECONCILIATION.md).

## Other accepted connections, with their boundaries intact

- **Every modulus:** exact factorial-transform context functions and clocks;
  context content detects universal divisibility, with zero states exactly
  at odd moduli. Actual H5/H11 projector factors create content5 from two
  content-one factors. [Arithmetic](all_moduli/PROOF.md) · [Content](all_moduli/CONTEXT_CONTENT.md).
- **Actual finite prime states:** exact counts5,49,685,322101 at3,5,7,11,
  with weighted/binary context simulators and precise source-padding cycles.
  General odd-prime cycle formulas describe the formal space only.
  [Simulation](all_moduli/WEIGHTED_BINARY_SIMULATION.md) · [Padding](all_moduli/PRIME_PADDING_ORBITS.md).
- **Spectral information loss:** weighted first-blindness is n5; M23 restores
  exact H there. Even regular critical-circle data leave an exact t^4/27
  oscillation on the specified rank-two circle. [Boundary](riemann_five_eleven/WEIGHTED_ZETA_BOUNDARY.md)
  · [Regular family](regular_spectral/REGULAR_CRITICAL_LINE.md).
- **Integral versus mean data:** four-face conservation gives the sharp binary
  support floor n-3; weighted second moments require M22. [Cocycle](cocycle/PROOF.md).
- **Recovered RH/5/11 history:** corrected and typed, not a proof of RH/GRH or
  a map from the golden discriminant to the order-five7/8 stability constant.
  [Historical recovery](riemann_five_eleven/HISTORICAL_RECOVERY.md).

The spectral/RH-adjacent work supplies rigorous sidecar and reconstruction
boundaries, **not the dimension-cutoff improvement**. It does not rule out
spectral sufficient inequalities. The arithmetic5/11 projector, clock and
affine55 actions are separate precise maps, not shared-prime-label inference.

[**Synthesis**](SYNTHESIS.md) · [**Concept board**](BOARD.md) · [**Verification**](VERIFICATION.md)
· [Final scope](CHECKPOINT4.md) · [Next-session handoff](FINAL_HANDOFF.md)

All-prime/full-composite contextual realization, exact integer-spectrum
completeness, general scalar zero-rigidity and LRC(14) remain OPEN.
