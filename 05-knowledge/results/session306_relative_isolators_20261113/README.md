# Session 306: relative isolators and a physical Collatz edge

**Status:** Proved scoped results on both problems; neither conjecture is closed.

1. [An ordinary-minus Collatz edge by one arc flip](MATCHED_COLLATZ.md). A
   single adjacent slide in a strong tournament from the crossing matching
   chart realizes 3H−1=2H′ exactly, for every parameter in an infinite
   family, with no change of order. Comparing the 3^b, 2^b, and constant
   coefficients classifies all uniform Collatz edges in the two crossing and
   nested matching charts, including zero outer widths. The only edges are
   this minus family and its reflection. Its output has unequal outer widths,
   so the edge cannot reenter its own input class. A return must change charts
   more substantially: no fixed signed Collatz macro of length ≥2 maps one of
   these chart families to either other family, even if intermediate counts
   leave the chart.
2. [Marked symmetry breaking](ROOTED_ISOLATORS.md) and [the exact visit
   hierarchy](RELATIVE.md). The attached FMCAD paper's ExactlyOne is safe
   on the actual port-preserving extension orbits, not on unrestricted
   child isomorphism classes. Universal unit fixing does not carry over to
   the bits of a fixed attachment cube. For uniform module substitution, an
   exterior budget of m vertices sees exactly the first floor(m/2)+1
   normalized ordered path-cover coordinates: no more are needed, and cyclic
   contexts recover every one. This locates precisely why a two-visit
   neutral preparation can cease to be neutral after a change of core.
3. [The infinite ordinary-Collatz survivor language](COLLATZ_LANGUAGE.md)
   and [an escape constraint](ESCAPE_BUFFER.md). The exact noncontracting
   parity-prefix language is not context-free, and its 2-adic survivor set
   has Hausdorff dimension h₂(log₃2). No regular safe sublanguage attains
   that exponent. Negative periodic shadows are dense but do not provide
   positive counterexamples. The escape note states an explicit logarithmic
   buffer for a hypothetical nonrepeating positive orbit and explains how
   the repository's inherited summability result connects the exact
   survivor set to the divergence half of ordinary Collatz.

The shared lesson is operational: a count, an isomorphism representative,
and a safe parity prefix each discard different information. Prove that the
data you keep support the next operation before identifying two states.
An exact physical edge is useful evidence about this problem, but it is
not a closed compiler for Session 305's reverse blocks.

## Replay

Run `python verify_matched.py`, `python verify_relative.py`, and
`python verify_rooted.py` next to this file. Also run with `python -O`;
the verifiers require only the standard library. The attachment itself is
linked from the mathematical notes; the PDF is not copied into this packet.
