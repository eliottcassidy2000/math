# Stern depth packets and the prime Kloosterman bound

**Primary reference pin, checked 2026-08-25.** This sidecar records the sole
external estimate imported by
[THM-4061](../../01-canon/theorems/THM-4061-stern-depth-modular-hyperbola-and-prime-apex-balance.md).

## Weil -- *On Some Exponential Sums*

- **Primary / freshness:** André Weil,
  [*On Some Exponential Sums*](https://pmc.ncbi.nlm.nih.gov/articles/PMC1079093/),
  *Proceedings of the National Academy of Sciences* **34** (1948), 204--207,
  [DOI 10.1073/pnas.34.5.204](https://doi.org/10.1073/pnas.34.5.204).
  **PUBLISHED / stable.**
- **Imported role:** for an odd prime `p` and nonzero residues `h,k`, the
  complete Kloosterman sum

  ```text
  sum_(a in F_p^*) exp(2*pi*i(ha+k/a)/p)
  ```

  has absolute value at most `2sqrt(p)`.
- **Repo consumer:** THM-4061 completes the representative-parity function
  on `F_p`, applies this pointwise bound, and derives
  `S(p)=O(sqrt(p)log^2 p)` for the THM-4059 depth packet.
- **Does not prove:** THM-4061's modular-hyperbola identity or Fourier
  coefficients, cancellation between distinct Kloosterman sums, a composite
  or prime-power packet asymptotic, a zero-packet classification, an Euler
  product, any Khinchin constant statement, or LRC(14).
