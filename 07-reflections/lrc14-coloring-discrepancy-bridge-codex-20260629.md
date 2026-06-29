# LRC14 Coloring/Discrepancy Bridge - Codex 2026-06-29

The useful thing from going back through the coloring work is not just that
LRC can be phrased as a coloring problem.  We already knew that from the
distance-graph/circular-coloring notes.  The sharper lesson is that every
coloring is a quotient, and every quotient needs guardrails about what it is
allowed to forget.

The older threads say this in several dialects:

- the centered edge variable `s_e=A_e-1/2` is the pair-first color bit;
- the W-polynomial only multiplies cleanly before edge colors are collapsed;
- the distance-graph LRC model makes a time into a regular circular coloring;
- the colored discrepancy bound says CRT placement must be counted by color;
- the Haar zipper cocycle says margins alone are not enough;
- the Paris-Harrington miniature says a color needs extension rank.

The AP84 tail is now a perfect place to make that precise.  HYP-3438 gives a
mod-`35` gate law:

```text
outer gate iff 7 does not divide m
inner gate iff 5 does not divide m
```

HYP-3456 gives the floor correction:

```text
N(m)=floor(12m/35)+d[(m-1) mod 35].
```

The new audit shows the gate color and the floor-correction color are not the
same quotient:

```text
both_outer_inner -> d values [1,2]
d=1 -> gate buckets both_outer_inner, inner_only_5_gate, outer_only_7_gate
```

The pair `(gate_bucket,d)` is clean for the local AP84 target, but the
endpoint phase still is not periodic: `m=1` and `m=36` have the same residue
and correction color, while HYP-3457 makes `m=1` a mixed endpoint transient
and HYP-3454 makes `m=36` rank-one `E/E`.

That is the main synthesis:

```text
gate color + floor color + height color + endpoint phase
+ branch mirror + incident C3/Qsqrt(-7) color + Haar zipper cocycle
```

is the legal color packet for the AP84 splice.  Anything smaller has to prove
the missing coordinate is reconstructible, or route it to HYP-3438/HYP-3453,
HYP-3455, owner-current, two-adic descent, exact-period/state-lift, or SPEC
debt.

This changes how I would try the next proof step.  Instead of proving a scalar
component-count statement, build a color-legality matrix for the HYP-3438
survivor gates.  Rows are gates/components.  Columns are the seven packet
colors.  The theorem target is:

```text
the local-to-global AP84 splice is a homomorphism for the color-packet product.
```

That is exactly the Haar product lesson, but now finite and arithmetically
named.  If the product is closed, the AP84 bridge can be spliced into
HYP-3439 without illegal quotient loss.  If not, the first missing color names
the next residual debt instead of letting it dissolve into "discrepancy."

Tournament Analysis used quotient/color carriers as vertices, not runners or
arcs.  The path ranked the full labelled color-packet theorem first, then
residue gate plus floor word, Haar zipper repair, endpoint phase, incident
`C3/Qsqrt(-7)` router, branch-mask discrepancy, distance-graph coloring, raw
mod-`35` gate color, and raw scalar escape count.  The tournament is
transitive with no directed `3`-cycles.

The challenged assumption was that balanced coloring classes are safe.  They
are not.  Balance helps only after the quotient has declared which LRC
predicate it preserves and which colors it has intentionally forgotten.
