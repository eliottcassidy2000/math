# Odd eighths: uniform formulas and complete family classifications

**Session 115, 30 October 2026. PROVED scoped results, with independent audits.** The major conjectures remain open. Classical Walsh, Lucas, finite-field and transfer mechanisms retain their credit; no literature-priority or formal-certification claim is made.

Read [the synthesis and formula atlas](SYNTHESIS.md) for the explicit maps, observer shifts and transport boundaries. [Checks](CHECKS.md) give reproduction commands and finite universes.

## Complete classifications and uniform formulas

1. **[Norm gauges and Boolean phases](quadratic/ODD_EIGHTH_PHASES.md).** All sixteen coherent relative Fano gauges have linear-coordinate orbit sizes 1,7,7,1. Their cubic sector gives the odd-eighth Walsh profiles. In every binary dimension d>=3, a Walsh coefficient of magnitude 2^d−2 characterizes exactly a single-point defect of an affine phase. The chosen base table and character label matter; this is not an invariant of abstract octonion isomorphism.

2. **[Affine groups and centralizer weights](switching/ODD_EIGHTH_GROUPS.md).** Every odd linear complement H in GL₃(2) has order 1,3,7 or21. The even-order proportions of F₂³⋊H are respectively 7/8,5/8,1/8 and3/8: all four fractions are values of one observable on one completely classified family. In every binary dimension the odd-order proportion is the average of 2^(−dim Fix(h)). A second theorem covers every finite group with Sylow D₈: its even-order proportion is (2c+1)/8, where c counts involution classes. No arbitrary switching-group realization is assumed.

3. **[Word-address transport](unitary/WORD_UNITARY_CLASSIFICATION.md).** At depth A, changing an odd multiplier q to p preserves products, real values, positivity and permutation-pullback structure iff p≡q modulo 2^(A−1). The first failure for distinct multipliers is v₂(p−q)+2, with an exact signed half-character defect. Thus 3→5 first fails at the odd eighths, even though linear unitarity survives.

4. **[Harmonic-tail families](lrc/HARMONIC_COMPLETION.md).** For the almost-AP core C_q={1,...,q−1,q+1}, harmonic completion has a complete optimum-preservation threshold, complete maximizing-point census, Möbius floor counts and physical boundary weights. The high-chamber spectrum is exact; the lower chamber is not claimed classified. At q=7 the weight is w₈(j)=88/105+(32/105)χ₂(j), retaining the quadratic character modulo8. The critical owner change gives 8/49, not the old two-owner prediction 8/7.

5. **[Every one-reset Collatz controller](collatz/GENERAL_RESET.md).** For every a>1 and L>=0, w=1^L a has an explicit maximal ordinary all-layer mean over all positive entry words, including all lifting stalls. Every maximizing entry is classified by a cyclic-intersection modulus and a complete mixed-radix no-borrow rule. Lucas binomial support is an exact special case, not a description of every radix. [The a=2 family](collatz/ONE_RESET_CLASSIFICATION.md) has a uniform strict clock bound below1/2. [Constant-letter controllers](collatz/CONSTANT_CONTROLLERS.md) have the sharp clock bound 5/8, with equality precisely at the all-one controller and all-one entries. No universal bound for arbitrary controllers is claimed.

## Explicit cross-field maps and a non-transfer theorem

**Group order to a Boolean measurement.** For odd-order h, let N_h be its binary norm operator and set f_h(v)=1[N_hv≠0]. This is exactly the even-order indicator of the affine element (v,h). Group even-order mass is therefore an average of named zero-character Boolean readings. Identity h supplies the OR phase, in the unique 7/8 cubic relative-gauge class. Normalizing that gauge shifts the observer from character0 to111; forgetting this shift changes the reading.

**A sharp obstruction to transporting physical weights.** On odd eighths, χ₂(j)=√2 cos(2πj/8). The word-address unitary U₍₃→₅₎ sends it to −iχ₋₂(j), so a positive LRC endpoint weight becomes nonreal. Positivity-preserving transport requires p≡q modulo4; fixing this particular weight requires the stronger congruence modulo8. Shared phase labels or Hilbert norms alone cannot transfer the weighted conclusion.

A [concurrent conditional-energy decoder](../odd_eighths_uniform_families_20261030/PHASE_BRIDGE.md), preserved and credited in the final merge, gives a different valid route from a phase-conditioned Collatz energy ensemble to the LRC weight. It retains clock parity and physical point multiplicities; it is not the failed multiplier-changing unitary. See the [synthesis addendum](SYNTHESIS.md#7-concurrent-merge-a-valid-conditional-energy-decoder).

The Collatz connection to Lucas/Frobenius is separately exact: maximizing every return layer becomes a no-borrow digit predicate, with completion-word data retained. Its parameter-Haar 3/8, orbit means, clock rates and LRC point weights are different measures and are never silently identified.

## Honest frontier

The packet does not settle Collatz convergence, general LRC, unrestricted switching enumeration or a general cross-field equivalence. It closes the named families and specifies precisely where their quotients lose information. [Independent audits](CHECKS.md) distinguish written all-parameter proofs from finite computation.
