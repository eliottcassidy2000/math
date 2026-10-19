# Independent audit of the polynomial-bit contact theorem

**Verdict: PASS.** No proof-critical error was found in
`square_sum_contact_output_height_20261019.md`. The audit treats the previously
proved non-torsion and signed recurrence results as inherited theorems, while
checking their hypotheses and their use in the new height argument.

## Checks completed

- **Integral model.** Direct symbolic expansion verifies the change
  `xi=9x+3a2`, `eta=27y` and both coefficients in equation (4). The displayed
  bounds `12C²`, `108C³`, `513C⁴`, and `5400C⁶` are valid under `0<A<2C`.
  Hence `R=6C` satisfies all four weighted inequalities. The model is
  nonsingular, and its anchor has rational coordinate heights `O(L)`.

- **Division polynomials.** The normalization removes all divisions from
  recurrences (7)–(8). In an odd recurrence the four even-index factors
  contribute exactly `(2eta)^4=16G²`; the even recurrence has no residual
  factor. Weighted degrees and the induction base norms were independently
  checked for indices 1 through 7. Their norms are
  `1, 1, 22, 88, 7072, 123628, 51453588`, all within (9).

- **Uniform induction.** The even products have sum of squared child indices
  `n²+6`, and the odd products have `n²+3`. Every child is smaller than `n`
  once `n>=8`. The norm cost is at most `log 288`, and inequality (10),
  including its endpoint `n=8`, supplies more than the required slack after
  multiplication by 4. Weighted evaluation therefore proves (11) uniformly.

- **Heights and exceptions.** Integer division-polynomial evaluations bound
  the numerator and denominator of the multiple before reduction. Equation
  (13) follows by clearing the denominator in the integral cubic equation.
  Translation by the fixed anchor, model changes, and the ordinary inverse
  chart each have fixed degree, so their heights remain `O(q²(L+log(|q|+1)))`.
  Identity, inverse, doubling, and zero-denominator cases are accounted for.
  The inherited successful arcs avoid the omitted inverse-chart points.

- **Output and complexity.** The budget satisfies `N=O(L)`. Substituting
  `|q|<=2N` gives `O(L³)` bits for coordinates, the primitive output, and the
  rational contact certificate. Reducing after each rational operation
  bounds intermediate integer sizes as well. A constant number of
  `O(L³)`-bit values suffices; iterative schoolbook arithmetic and the stated
  unamortized Euclidean bound give `O(L^10)` time and `O(L³)` working space.

## Source and scope

The cited [MIT 18.783 Lecture 5](https://ocw.mit.edu/courses/18-783-elliptic-curves-fall-2025/mit18_783_f25_lec05.pdf)
was checked directly: section 5.5 supplies the recurrences, Lemma 5.20 the
parity normalization, and Theorem 5.21 the multiplication formula. Its
characteristic and nonsingularity requirements hold over `Q`; Lemma 5.24
also confirms the numerator-denominator coprimality behind the absence of
division-polynomial zeros at a non-torsion point.

This audit certifies the uniform proof, not an optimal exponent or constant.
It does not extend the conclusion to prescribed spread directions, narrow
windows, or arbitrary forbidden sets. The companion implementation was not
present when this audit began and is not certified by this note.
