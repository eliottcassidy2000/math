# Chan–Pak: exact window viability and its Hardy–Littlewood shape

**Universal exact placement remains open.** This packet separates the
number of presentations, the joint residue image of each window, and
the positive support of exact tournament counts. It also classifies
the full slope-preserving transport algebra of one physical carrier.

1. **[Exact window generating function](WINDOW_GENERATING_FUNCTION.md).**
   A four-rule parity recursion counts the true pair checksum on
   each prescribed sum and window. Applying the exact quadratic
   readout makes coefficient positivity equivalent to a physical hit
   in that window. There are arbitrarily high windows passing mod 4
   and mod 7 separately but not the joint test mod 28. Thus
   local marginals within a proposed window cannot simply be
   multiplied.
2. **[Real density and normalized local factors](LOCAL_DENSITY.md).**
   The unrestricted cyclic carrier has ordered-pair mass
   (X/5760) log X + KX + O(sqrt(X) log X), with an explicit
   Gamma constant. For each fixed modulus the same two main
   terms divide exactly equally among compatible odd residues,
   with power-saving relative error. The real mass is uniform over
   *logarithmic* aspect ratios, not concentrated near a tangent:
   such near-balanced windows have only O(X^(3/4)) exact support.
   **[A tournament countermodel](MARKED_WORD_COUNTERMODEL.md)**
   shows why this cannot prove positivity by itself. Marking each
   10 in a transitive-spine word gives uniform local factors
   and X log X presentations, while the true support omits at
   least a relative 1−log 2 of the odds in lower density.
3. **[Complete transport spectrum](ALL_WIDTH_TRANSPORTS.md).**
   Widths are exactly 12k. A unique periodic ternary moment
   word describes 4·32^(k−1) matched charts and gives their
   exact Laurent displacement polynomial. For k>=9 it
   specifies both finite fringes and every central even signed
   shift. The [primitive law](PRIMITIVE_REPAIR_LANGUAGE.md)
   gives 4·28^(k−1) primitive charts and an exact free
   factorization. There are infinitely many truly primitive
   +2 charts, and a [specific 36-bit primitive return](SUFFIX_CHART_SWITCHES.md)
   is explicit alongside the inherited 12-bit unit slide.
4. **[Exact maximal windows](WINDOW_ORBIT_SERIES.md).**
   Arbitrary mixing of the two explicit +2 codes gives
   consecutive odd counts with both outer addresses still free.
   The canonical endpoint recurrence and geometric window kernel
   give an exact target-viability series. The primitive
   +2 codes are suffix-free and extend the same window law to
   the entire infinite unit alphabet. Their limiting presentation
   length PGF is (1−rho)^2/(1−rho t)^2, with rho obtained
   from the coefficient of loss +2 in the primitive series.
   The [full width tariff](SHIFT_MULTIPLICITY_AND_TARIFF.md)
   proves that even all allowed long transports cover a
   density-zero set above logarithmic-in-logarithm depth.
   There are many unit charts but they are not many independent
   estimates at any prescribed absolute target.

The most promising change is to put the coarse selector in the
unbalanced part of the actual representation geometry and make
the correction act there, with canonical chart multiplicities
or a genuine bound on the remaining Fourier contribution.
Extending the balanced tangent or adding more congruence tests
cannot supply exact absolute positivity by itself.

## Reproduction

Run the six `replay_*.py` scripts using Python 3, and compare the
matching receipts. The direct hyperbola replay is C++17:
`g++ -O3 -std=c++17 replay_hyperbola.cpp -o /tmp/vshape && /tmp/vshape`.
