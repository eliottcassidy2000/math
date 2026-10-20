# Scoped correction: a critical limiting slope is not a finite integer-cycle slope

Target: `01-canon/theorems/THM-4717-sturmian-circles-of-the-collatz-map-explicit-three-smooth-ladder-golden-circle-and-the-critical-slope.md`.

The file literally states, in its title and in the last sentence of Statement 2:

> the critical slope is exactly the slope of every integer `3x+1` cycle (THM-4633).

Its related-entry gloss also says integer cycles have slope at the convergents of log_3 2 without the conditional denominator/minimum qualification.

The exact-equality claim is false for every finite nonzero positive cycle: its slope p/q is rational whereas log_3 2 is irrational. Concretely, the shortcut cycle 1→2→1 has q=2,p=1 and slope 1/2. The cited THM-4633 actually proves a strict approximation inequality, with an additional hypothesis for the convergent conclusion.

Correct replacement under the file's binary-word convention: for a positive cycle with q shortcut steps, p odd steps, and minimum odd member m,

    2^q / 3^p = product over odd cycle members n of (1+1/(3n)),
    0 < q log 2 − p log 3 < p/(3m),
    log 2 / log(3+1/m) <= p/q < log_3 2.

The non-strict lower bound can be equality for the trivial cycle. In particular a family with m tending to infinity has p/q tending to log_3 2 from below. The convergence statement is conditional on such a family existing. The real-series failure at the critical parameter is a separate claim; this correction does not audit or endorse the topology and pointwise convergence assertions in THM-4717.

Recommended minimal repair:

- Replace the exact-cycle-slope phrase in the title and Statement 2 with: “the critical limiting slope approached from below by positive cycles whose minimum tends to infinity.”
- Replace the related-entry gloss for THM-4633 with: “positive cycles have slopes below log_3 2, with the quantitative minimum-dependent approximation and conditional convergent criterion stated there.”
- Record the rational-versus-limiting-parameter error without retracting unrelated ladder identities.

The new local rotation-defect theorem is independent of THM-4717 and uses only elementary carry identities.
