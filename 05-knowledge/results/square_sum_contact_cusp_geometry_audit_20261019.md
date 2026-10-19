# Independent geometric audit of the cusp lower bound

**Result: PASS.** This audit checks the proof in
`square_sum_contact_cusp_lower_bound_20261019.md`, with the period convention
and upper bound in the recurrence note. It does not certify the companion
computational checker or claim a sharper first-success asymptotic.

- **Complete real curve and both charts.** Since `K>1` and `|f(u)|<=1`, the
  normalized coordinate `z/(1+u²)` never vanishes. Compactifying the parameter
  line therefore gives two real circles, including the points over infinity.
  The differential `−du/z` is nonsingular there. Its absolute value in the
  angular coordinate is exactly the one used in the proof. Both accepted
  positive-spread charts lie in the enlarged set `0<f<1/2` on these circles.

- **Cluster geometry.** Solving `−sin(4θ)<=1/2` gives exactly the two stated
  closed arcs modulo `π`. Each contains a middle interval of width `π/4`
  with `K−f>=1`, and two end intervals of total width `π/12` with
  `K−f>=1/2`. Thus its invariant length is at most `C₀`. The angular shift
  by `π/2` is fixed-point-free, preserves the differential on either circle,
  and partitions a full circuit into equal invariant lengths. It is an
  exact half-period. From any base in one arc, the other arc is consequently
  within `C₀` of the half-period displacement, while the original arc is
  within `C₀` of zero.

- **Translating step and explicit constants.** The displayed doubling
  formula puts `2P` on the negative-`Y` half of the identity component.
  Integrating `dX/(2Y)` from infinity along that half gives `α(K)`, with no
  missing factor of two. Substitution `r=sqrt(X+1)` at `K=1` gives
  `α₀=sqrt(2) log(1+sqrt(2))`. After `v=X−K`, restricting to `v>=1/4`
  and comparing `v+δ` and `v+2+δ` gives inequality (10) in the claimed
  direction. All rational comparisons in (11) are valid; in particular
  `3500*294 > 341*3009`. Thus `α(K)>C₀` throughout `0<δ<=1/1000`.

- **Even and odd indices.** For a nonzero integer `n`, the unreduced
  displacement `|n|α` exceeds `C₀`. If it is below `Ω/2−C₀`, it lies in
  the forbidden gap before any wrapping is possible. Larger displacements
  already satisfy the desired lower bound, so arbitrary winding causes no
  loophole. For odd indices the two stated exact inverse parameters have
  negative spread: the first has `u<−1` for `1<t<4/3`, and the second has
  `0<u<1` for `t>sqrt(3)`. Admissibility excludes the intervening interval
  because `f` is increasing on `[4/3,sqrt(3)]` and starts above `1/2`.
  These bases are themselves unsuccessful, and `|s+2n|>=2|n|−1` gives
  exactly the claimed one-unit loss.

- **Uniform cusp rate.** The integration endpoint and changing factors
  defining `α(K)` stay away from zero near `K=1`, giving
  `α(K)=α₀+O(δ)`. Together with
  `Ω=(1/sqrt(2))log(1/δ)+O(1)`, this yields the lower coefficient
  `1/[2 log(1+sqrt(2))]` with an absolute additive constant. The same
  estimate applies at both anchor cusps because its bounds depend only on
  `δ=f(t)`. Combining it with the established signed-index upper bound
  proves the asserted `Θ(log(1/f(t)))` order for this prescribed orbit.

No gap or counterexample was found in the requested geometric and
asymptotic scope. The proof correctly does not assert a first-success
asymptotic equivalent, a bound for other translating points, or a bit
complexity result.
