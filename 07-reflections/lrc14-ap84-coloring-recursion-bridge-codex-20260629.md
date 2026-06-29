# LRC14 AP84 Coloring-Recursion Bridge, codex-2026-06-29

The prompt asked to reconnect with old coloring work and explore freely.  The
useful return route was HYP-2247: a coloring should not be counted without its
extension-rank sidecar.  In Paris-Harrington this becomes proof-theoretic
growth.  In the AP84 tail it becomes something much smaller but much more
usable: a finite `35`-state coloring automaton with an endpoint-rank subcolor.

For the canonical family

```text
S_m = {1,2,...,11,13,84m}
```

HYP-3456 already gives the number of high `E:84m` gaps meeting the left
HYP-3431 corridor `C1=[8/49,6/35]`:

```text
N(m)=floor((504m-6)/70)-floor((96m-13)/14).
```

Subtracting the carrier clock `floor(12m/35)` gives the residue color
`d_r=N(r)-floor(12r/35)` on `r=1..35`.  That was the visible coloring.  The
missing HYP-2247 coordinate is the address/rank of the canonical endpoint
witness from HYP-3433/HYP-3454:

```text
a_m = ceil(48m/7).
```

If `L_m=floor((96m-13)/14)+1` is the first high-gap index meeting `C1`, then
the selected witness rank is `a_m-L_m+1`.  Write `m=7q+s`, `0<=s<7`.  Then

```text
a_m = 48q + ceil(48s/7)
L_m = 48q + floor((96s-13)/14) + 1
rank_C1(a_m) = ceil(48s/7) - floor((96s-13)/14).
```

The seven residues give

```text
s:    0 1 2 3 4 5 6
rank: 1 2 2 2 2 2 2
```

So the endpoint-rank subcolor is exact:

```text
rank_C1(a_m)=1 iff 7|m, otherwise 2.
```

That is the new little hinge.  The mod-35 word counts all high-gap escapes,
but the mod-7 subcolor tells where the homogeneous endpoint trace sits inside
those escapes.  This is exactly the kind of coordinate the older coloring
sessions warned us not to forget.

## What Changed

Before this pass, the AP84 packet had endpoint intervals, escape counts, and
finite transients.  After this pass, it has a coloring-recursion state:

```text
phase color        mixed for m=1..4, pure for m>=5
boundary color     d_r in {0,1,2} on residues mod 35
endpoint subcolor  rank_C1(a_m) in {1,2} on residues mod 7
```

The mixed phase is the finite bad branch.  It appears in the first period at
`m=1..4` and is gone after the outer extension `m -> m+35`.  Unlike true
Paris-Harrington, this does not grow out of reach; it collapses to a small
automaton.  But the proof-design lesson is the same: a quotient that keeps only
escape count and forgets endpoint rank is too coarse.

The incoming S311 coordination note sharpened this rather than distracting from
it.  S311's no-free-slider rule says the finite transient packet is usable only
because the mixed endpoints have explicit addresses.  HYP-3458 is the tail
analogue: once the transient windows are closed, the pure period still needs an
explicit endpoint-rank address so the `35`-state count cannot silently slide
between first-gap and second-gap witnesses.

## Connections Back To Colorings

HYP-2247 split the problem into a side-choice filter and a bad-child extension
rank.  Here the side choice is the corridor/gap color, and the extension rank
is whether the selected endpoint remains first, second, mixed, or dead under
period extension.

HYP-2243's outer-extension usability theorem also fits cleanly.  The AP84 tail
is usable because extension by `35` preserves the boundary word up to the
shift `N(m+35)-N(m)=12`, preserves the endpoint-rank subcolor, and kills the
finite mixed phase.  That is a concrete finite version of an embedding theorem:
the old witness embeds into the next period with a known address.

The older Borel/embedded-maximality work says a constructive witness must keep
the address that computes it.  Here that address is literally `ceil(48m/7)`,
and the address rank is the color that prevents the escape clock from becoming
just a scalar.

## Pull Toward LRC14

The AP84 side of HYP-3439 now wants to be stated as a packet splice:

```text
HYP-3457 finite transients
+ HYP-3454 endpoint interval
+ HYP-3456 floor count
+ HYP-3458 coloring-recursion sidecar
=> AP-tail rank-5 descent packet
```

For noncanonical rows, the same pattern suggests a search target.  Do not ask
only whether a low-rank gate or component exists.  Color each candidate by:

```text
component/gate type
endpoint-owner address
minimal cover rank
outer-extension child count
whether the address survives, shifts, or dies
```

Then prove that every bad color either has finite height, routes into the
HYP-3455 seven-owner gluing clause, or reconstructs the strict HYP-3437
rescue cut.  That is the AP84 finite automaton idea exported to the noisy bank.

## Creative Extensions

One possible general theorem shape:

```text
If an AP-tail family has a finite residue coloring whose bad phase support is
bounded, whose endpoint-address rank is periodic, and whose period extension
is monotone in escape count, then the family reduces to finitely many
transients plus one symbolic endpoint-rank lemma.
```

For LRC14, this would turn a growing family of AP-tail witnesses into three
small obligations: exact corridor identity, exact endpoint rank, and finite
bad-color extinction.

The tournament version should use proof carriers as vertices:

```text
crt_state -> endpoint_rank_subcolor -> bad_phase_derivative
          -> outer_extension_shift -> boundary_word -> raw_count.
```

This is a useful guardrail because the raw count is last.  It is allowed to be
a shadow only after the endpoint address and extension derivative have been
kept.

## Next Hooks

1. Add a HYP-3438/HYP-3453 version of this scout: color survivor gates by
   endpoint-rank, cover-rank, and child survival under the natural AP/collar
   extension.
2. Convert the AP84 splice into a theorem interface with explicit inputs:
   HYP-3431 fixed corridor, HYP-3454 interval containment, HYP-3456 floor count,
   HYP-3457 finite transients, and HYP-3458 endpoint-rank sidecar.
3. Test whether the `5 x 7` CRT grid is the cleanest carrier for the canonical
   gate law: `7|m` controls outer `(13,7)` gates and the endpoint rank, while
   `5|m` controls inner `(11,5)` gates.
4. Revisit HYP-2247's baby coloring code with AP84-style sidecars: instead of
   raw bad coloring counts, store the first surviving address and its child
   rank.  This may give a portable "coloring recursion compiler" for later
   LRC packets.

The main proof movement is modest but real: the AP84 clock is no longer just a
periodic count.  It has a finite recursive state and a symbolic rank lemma,
which makes it much easier to state what a legal quotient must remember.
