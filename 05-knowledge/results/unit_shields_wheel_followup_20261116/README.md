# Restricted two-apex census, exact unit shields and wheel follow-up

- [Spine 14, second apex at most three wins](apex/CENSUS.md):
  15,400,960 oriented tapes, 1,138,644 distinct Hamilton counts.
  The previous first hole at 256,521 fills; the next is 729,421.
  This census is complete under the stated restriction only.
- [Composable one-arc shield](apex/SHIELD.md): a balanced
  pair of twelve-site two-apex blocks has gaps ±2; arbitrary
  concatenation adds these gaps exactly. Mirrored exteriors also
  preserve a +2 move. [The module formula](SHELL.md) explains the
  earlier constant guard as a literal factor-one case: the gap is
  2H(Z) and the baseline reads P_Z(2), P_Z(3).
- [A different physical seam](bridge/SHIELD.md) uses a primitive
  linear projector equal to one against an entire growing module.
  It has exactly two strong right rays for its first projector
  and [three more distinct projectors](bridge/ADDENDUM.md).
  [Two-contact cancellation in the binary reader](READER.md)
  gives Θ(√H) distinct +2 anchors on selected scales without swapping
  seeds; it does not retain the next-read slope.
- [Wheel follow-up](wheel/WHEEL_FOLLOWUP.md): an exact slice proves
  the corrected nesting criterion with a separate dyadic-only
  alternative. Least periods for q=11,13,17,19,23,29,31 are
  352,416,1088,1216,1472,1856,1984 respectively. W(120) cannot
  serve both signs of one coefficient, and literal period 6424
  does not cover at all. [A sieve bound](CUTOFF.md) proves that
  the least ambient period is 64q for every prime q≥52,000.
- [Synthesis](SYNTHESIS.md): constant-guard pumping obstruction,
  phase projection and exact landing, 8 and 9 / the square-sum
  seam, Golden/Fibonacci comparisons, fences, and the role of
  Catalan's 8 and 9 in the missing primitive clock.

The strongest new reusable construction here is an exact
single-arc signed gap after an arbitrary word of balanced
atoms. The gap is programmable but grows at most linearly in
the number of these atoms; the count scale grows exponentially.
Finding compatible arms at an arbitrary prescribed absolute
target, and Chan–Pak, remain open. These results do not turn
the earlier O(N^(2/5)) bracket into an exact target theorem.
