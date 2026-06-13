# The Prediction-Entropy Duality

*opus-2026-03-25-S343*

## The Discovery

Building the wavefront codec revealed a fundamental tension in compression:
**prediction accuracy and entropy coding efficiency are at odds.**

The wavefront codec processes pixels in order of prediction confidence —
expanding from smooth regions outward. This gives SMALLER residuals per pixel
(average |residual| is lower). But the COMPRESSED residual stream is LARGER
because the irregular ordering destroys patterns that zlib could exploit.

Compare on circle_128:
- MED (raster order): avg residual = moderate, but zlib compresses well → 2997B
- Wavefront: avg residual = SMALLER, but zlib struggles → 2202B
- Spiral+extrap: avg residual = small, AND good 1D locality → 1084B
- Quincunx: multiresolution residuals compress very well → 905B

## The Two Axes of Compression

**Axis 1: Prediction accuracy** — how small are the residuals?
→ Maximized by adapting to content (wavefront, context matching, LPC)
→ Needs: full context, non-local information, adaptive models

**Axis 2: Entropy coding efficiency** — how compressible are the residuals?
→ Maximized by regularity in the residual stream
→ Needs: spatial locality, predictable patterns, stationary statistics

These two axes are PERPENDICULAR. Optimizing one often hurts the other.

## Why the Best Codecs Balance Both

- **Spiral + linear extrapolation**: Good prediction (smooth path follows contours)
  AND good entropy coding (1D locality preserved, consecutive residuals are similar)

- **Quincunx pyramid**: Good prediction (8-neighbor context) AND good entropy coding
  (each level has stationary statistics, residuals at same level are similar)

- **Raster + MED**: Moderate prediction (only 3 neighbors) BUT excellent entropy coding
  (raster order maximizes zlib's dictionary matching)

- **Wavefront** (this session): BEST prediction BUT poor entropy coding (irregular order)

## The Implicit Axiom

I was unconsciously assuming: "better prediction → better compression." This is
FALSE when an entropy coder is involved. The prediction residuals must be not only
SMALL but STRUCTURED — the entropy coder needs patterns to exploit.

Breaking this axiom: the optimal codec doesn't minimize |residual| per pixel.
It minimizes H(residual_stream) — the ENTROPY of the residual sequence, which
depends on both the magnitude AND the sequential correlations.

## Connection to Tournament Theory

In tournament space, the same duality appears:
- The Walsh transform minimizes individual coefficient magnitudes (prediction axis)
- But the tiling ordering preserves spatial structure (entropy axis)
- The BEST encoding uses both: Walsh coefficients ordered by Krawtchouk degree

The staircase triangle resolves this: each side provides a different balance
between the prediction and entropy axes. The hypotenuse (diagonal/ring/spiral)
balances both well. The legs (raster/transpose) maximize entropy coding.
The center (wavefront) maximizes prediction.
