# Fejer Kernel Wild Session

**Instance:** codex-2026-05-30-fejer-wild
**Date:** 2026-05-30
**Status:** exploratory; computationally reproducible in `04-computation/fejer_kernel_wild_session.py`

## Why This Session

The existing Fejer thread is already strong:

- interval circulant eigenvalue magnitudes are exactly sampled Fejer kernels;
- additive energy equals spectral fourth moment / IPR;
- Paley is spectrally flat while the interval is spectrally concentrated;
- the remaining proof gap is the passage from spectral concentration to spatial cycle localization and then to larger OCF packings.

This pass treats that gap as the object, not as background.

## Hard Checks

For `S={1,...,(p-1)/2}` in `Z_p`, the script verifies three identical objects:

1. `Q_k = |sum_{s in S} exp(2*pi*i*k*s/p)|^2`.
2. `Q_k = sin^2(pi*m*k/p) / sin^2(pi*k/p)`.
3. `Q_k` is the discrete Fourier transform of the difference autocorrelation
   `a(d)=|S cap (S-d)|`.

So the cleanest Fejer object is not directly a cycle object. It is the
autocorrelation shadow of the connection set. Cycles enter after orientation
and distinctness constraints are imposed.

## The Useful Distinction

There are two localization profiles that should not be conflated:

- **Difference localization:** `a(d)=|S cap (S-d)|`. This is the pure Fejer /
  additive-energy profile.
- **Pinned cycle localization:** `J_k(0,v)`, the number of simple directed
  `k`-cycles passing through both `0` and `v`.

For `k=3`, `J_3` is an oriented wedge fold of the same connection-set geometry:
one counts third vertices `w` making either `0 -> v -> w -> 0` or
`0 -> w -> v -> 0`. For the interval this becomes the known linear profile
`J_3(0,v)=min(v,p-v)`, while Paley is flat in the `p=3 mod 4` cases.

This suggests the proof should factor as:

```text
Fejer autocorrelation
  -> oriented wedge localization
  -> pinned simple-cycle localization
  -> overlap multiplicity distribution
  -> Omega independence advantage
```

The old phrase "spectral concentration implies spatial localization" is true
morally, but the missing lemma probably lives at the second arrow, not at the
first.

## The New Sharp Pattern

For every full orientation cube tested through `p=23`, the variance of the
pinned triangle profile is exactly the additive-energy axis:

```text
Var_{v != 0} J_3(0,v) = E(S)/(p-1) - (p^2 - 2p + 5)/16
```

The fitted slope equals `1/(p-1)` to roundoff, the intercept matches the
displayed formula, and the affine fit has roundoff-scale error. Since IPR is
already an affine transform of additive energy, this means:

```text
Fejer/IPR concentration <-> additive energy <-> J_3 localization variance.
```

So `J_3` variance is not an independent predictor. It is the first spatial
incarnation of the same Fourier/autocorrelation quantity.

## Observed H Pattern

For `p=7,11,13`, the script also compares that common
energy/localization axis with Hamiltonian path count.

The sign-reversal story survives:

- At `p=7`, Paley wins `H` despite minimal spectral concentration.
- At `p=11`, Paley still wins `H`, and `IPR` is only weakly anti-correlated with
  `H`.
- At `p=13`, interval-type sets maximize both `IPR` and `H`.

The new lens is that the phase transition is not "does spectral concentration
produce triangle localization?" It does. The transition is "when does that
localization beat Paley's larger low-order cycle count inside the OCF
hard-core partition function?"

## Mutation Picture

Single-pair flips away from the interval act like finite-difference probes of
the rearrangement inequality. They usually reduce:

- spectral IPR,
- the `J_3` gradient/variance,
- and, once the interval phase has turned on, `H`.

This is potentially more proof-friendly than global optimization. A local
optimality theorem might be easier to express as "every signed-pair flip lowers
the oriented wedge concentration" and then bootstrap to higher cycle packings.

## Wild Tangents

1. **Pinned-exclusion expansion.** Start from pinned closed-walk counts, which
   are directly spectral, then subtract repeated-vertex diagrams. For `k=3`
   there are no correction diagrams. For `k=5`, the first corrections are small
   enough to classify. This could become the missing bridge.

2. **Fejer as a diffraction pattern; Omega as crystallization.** The interval
   has a coherent diffraction peak. The OCF advantage is not "more cycles" but
   "cycles arranged into a packable crystal." Paley is better mixed and therefore
   worse at high-fugacity hard-core packing.

3. **Autocorrelation-to-wedge operator.** Define an operator `W` taking a
   connection-set autocorrelation profile to the oriented third-vertex profile.
   For intervals, `W` preserves triangularity. For Paley difference sets, `W`
   collapses to a constant. This may be a tiny theorem.

4. **Heat-flow analogy.** Normalized `Q_k` can be read as a probability
   distribution on Fourier modes. Paley is high-entropy; interval has effective
   participation ratio near `3`. The phase transition is a low-temperature
   ordering transition in the hard-core gas on odd cycles.

5. **Side-lobe arithmetic.** The interval's top Fejer mode gets the headline,
   but the side lobes decide finite primes. The `p=7,11` exceptions may be
   exactly the regime where side-lobe corrections feed `alpha_1` more than
   disjoint packing.

## Next Concrete Moves

- Derive the exact `J_3` wedge formula for arbitrary connection sets and add it
  as a lemma candidate.
- For `J_5`, enumerate the repeated-vertex correction diagrams and express the
  pinned simple-cycle profile as pinned walk profile minus corrections.
- Use `J_3` variance, `J_5` variance, and additive energy as a three-feature
  predictor for `H` across `p=17,19` sampled orientations.
- Try a local theorem: every single signed-pair flip away from an interval lowers
  the lexicographically sorted `J_3` profile in majorization order.
