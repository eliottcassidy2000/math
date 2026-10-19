# Independent proof-critical audit: bounded rational-section atlas

Date: 2026-10-19. Scope: the staged theorem note, its inherited recurrence
statement, and the optional generic independence argument. No main-file edits,
pushes, or new broad census were performed.

**Verdict: PASS. No proof-critical correction required.** The universal
polynomial numerical-height / linear-output-bit conclusion follows from the
stated hypotheses. The optional independence proposition also passes.

## Algebra and closed cusp interval

The independent script `square_sum_bounded_section_audit_20261019.py` reconstructs the
rational-function group law without importing the companion implementation.
It passed **35 exact checks**, including:

- `P,Q` on the cubic; the formulas for `Q+P`, both coordinates of `S=2(Q+P)`,
  and the inverse quartic coordinates of `Q+S=3Q+2P`.
- The full polynomial square identity (10), using the printed
  coefficient list for `Z`, independently of the group-law route.
- All four shifted factors, every printed coefficient of `-(N-tE)/h` and
  `4N-5E`, and both exact rational tail bounds.
- The cutoff/cap spread values; numerator/denominator degree 11, coprimality,
  coefficient norms 760 and 446, and constant 776516.
- Both leg-swap identities and the displayed derivative.

The sign proof is uniform over `0<h<=1/16`, not inferred from probes:
`A,C<0`, `B,D>0` give `E<0`; the positive displacement polynomial gives
`u>t`; the negative-tail bound gives `4N-5E>12`, hence `u<5/4`. Monotonicity
then gives the strict changed spread and admissibility. The square root is
nonzero because its square exceeds `1/2`. There are no omitted denominator
zeros in the direct chart.

## Global coverage and complexity

The swap involution exchanges the normalized legs. Its minimum selects the
low admissible cone, and a reduced parameter satisfies `n^2+m^2=delta C`
with `delta` equal to 1 or 2. This handles the odd/odd swapped representation.

Above the cusp cutoff, `f0>1/10` and `(2+f0)/f0<25<32`, so `J<=5`.
The inherited recurrence bound is exactly
`2 floor(100(4+5)/7)=256`, with 256 nonzero signed even candidates.
Its guaranteed point is in the ordinary positive-parameter, positive-root
chart and has changed spread. Skipping unusable inverse charts is therefore
safe; the elliptic additions themselves must retain their usual identity,
inverse, and doubling branches, as the note and implementation do.

The cusp estimate `C'<=776516(2C)^11` is correct. For the complement, a fixed
finite rational-arithmetic decision tree suffices even at exceptional
specializations: each valid operation multiplies rational heights by at most
a fixed power/factor. Thus every operand and the returned triple/root have
`O(log C)` bits. This proves an effective polynomial numerical bound, not
merely a search-index bound. The conservative cubic bit-time and linear
working-bit claims are valid with ordinary arithmetic and Euclidean reduction.
Primitive normalization preserves spread, and changed spread ensures a
new unordered primitive triple. The five-vertex graph consequence has the
claimed numerical bound.

## Optional independence proposition

At `t=1`, the special cubic is nodal, with smooth reductions `P0=(0,1)`,
`Q0=(3,4)`, `T=P0+Q0=(-1,0)`. The displayed multiplicative coordinate has
`zeta(P0)=-(1+sqrt(2))^2`, `zeta(Q0)=(1+sqrt(2))^-2`, and `zeta(T)=-1`.
Consequently reduction of a relation forces `a=b` by absolute values.

Independently checked: `h^2 X(S)->1`, `h^3 Y(S)->1`, and
`(-X(S)/Y(S))/h->-1`. Thus `S` is a nonzero point in the formal reduction
kernel. Over `Q[[h]]`, the linear term of multiplication is `n v`; every
nonzero integer is a unit, so no nonzero kernel point is torsion. This
excludes the remaining relation. No specialization-wise independence is
implied or used.

For the standard facts, verified the author-hosted primary exposition:
[J. S. Milne, Elliptic Curves](https://www.jmilne.org/math/Books/ectext6.pdf),
Chapter II §3, printed pp. 57–58 (smooth nodal group), §4, pp. 61–62
(smooth reduction subgroup/homomorphism), and Remark 4.3, p. 64 (formal
law and reduction kernel). The same local algebra applies to the complete
DVR `Q[[h]]`; its characteristic-zero leading-term consequence is explicitly
proved above. Adding this citation to the optional proposition would improve
traceability, but it is not a mathematical repair.
