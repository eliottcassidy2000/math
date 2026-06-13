---
source: codex-2026-06-03-S590
status: synthesis plus exact orbit-sheaf audit
tags: [LRC, rigidity, orbit, sheaf, monodromy, n14, tournaments, projection-defect, HYP-2127]
---

# Orbit-sheaf monodromy: the next rigidity layer

The prompt asked for rigidity types beyond local/global and to see everything
as orbits.  After pulling the latest work, the clean move was not to add one
more taxonomy beside HYP-2126.  It was to glue the two HYP-2126s together.

One HYP-2126 says rigidity is a type of orbit behavior: static unit orbit,
dynamic doubling orbit, spectral blocks, isostatic constraints.  The other says
rigidity is labelled fixed-point propagation: source, endpoint, fold,
denominator, pincer, middle state.  Those are not competing frames.  They are
base and stalk.

```text
orbit quotient = base
proof labels   = stalk
certificate    = section
rigidity       = section glues
defect         = monodromy or mixed fibers
```

That is HYP-2127.

Rebasing over HYP-2128 sharpened the two-face reading.  HYP-2128 is the
additive `+2`/triangular face: `2n-1` is the odd-square root of `8*C(n,2)+1`.
HYP-2127 is the multiplicative `x2`/boundary face: at even `n`, doubling exits
the unit witness sheaf into gcd-2 residues.

## The new seam: doubling is a boundary map

HYP-2124 says the AP witness set is exactly the unit clock orbit:

```text
W_n = (Z/n)^*.
```

HYP-2126 says doubling fragments this orbit when `n` is even.  S590 sharpens
that.  At even `n`, doubling is not merely a bad action on the unit orbit.  It
is not an action on the unit orbit at all.

It is a boundary morphism:

```text
(Z/n)^* -> gcd-2 residues.
```

For `n=14`:

```text
1  -> 2
3  -> 6
5  -> 10
9  -> 4
11 -> 8
13 -> 12
```

Every image has two total lifts under `x -> 2x`: one unit lift and one gcd-2
nonunit lift.  Example:

```text
2 has lifts 1 and 8 mod 14.
```

That is the crisp n=14 picture: the unit witness orbit is rigid, but the
doubling transport crosses a boundary where each sheet splits into unit/nonunit
lift ambiguity.  The apex `7` is not part of that six-residue image; it is the
self-antipodal chart that must be glued separately.

## The point-set warning

The old slogan

```text
vertex-transitive trienerment <=> regular polygon point-set
```

was too strong.  HYP-2125 gave the correction; S590 repeats it in sheaf
language.

Through `N<=18`, dihedral vertex-transitive point-sets split:

```text
31 regular polygons
52 imprimitive bracelets
```

All nonregular examples in the audit have gap period `2`.  So local vertex
profiles can glue with either trivial monodromy or order-2 block monodromy:

```text
regular polygon  = cyclic trivial gluing
bracelet         = dihedral two-block gluing
```

This is exactly the user's "sometimes outside, sometimes whole mesh" intuition:
the outside cycle is one sheaf chart; the mesh/block structure is the gluing
data.

## The tournament warning

The rooted tournament quotient audit is the direct projection-defect version.
At `n=6`, the full rooted class has `296` values.  Weaker lenses collapse it:

```text
unrooted_plus_score   196 values
split_no_cross         36 values
unrooted               56 values
delete_parent          12 values
score_sequence         22 values
root_score              6 values
```

Score-aware lenses preserve source and sink purity, but they do not preserve
all rigidity labels.  For instance, `split_no_cross` has mixed fixed-root and
orbit-size labels on `8` fibers and mixed parent labels on `12` fibers.

This is the same lesson as HYP-1977 and S535: unmarked quotients can be
beautiful and still forget the decisive observer-coupled data.  A quotient is
proof-safe only when the stalk labels are pure on every fiber.

## Rigidity menu

I now want to name the types this way:

```text
local       stabilizer-fixed basepoint or endpoint
global      transitive transport of a label
static      symmetry action preserves the witness orbit
dynamic     proof operator preserves the witness orbit
boundary    operator exits the witness orbit into a labelled stratum
quotient    forgetful fibers remain label-pure
gluing      local certificates agree on overlaps
monodromy   loops return the same certificate
spectral    character blocks isolate defect modes
isostatic   active constraints match degrees of freedom
pincer      witness and blocker fronts meet or export a labelled core
automaton   L/M/R side changes must pass through middle
```

The n=14 problem looks less like "find more symmetry" and more like:

```text
prove every boundary/gluing defect is owned.
```

Owned means endpoint-owned, denominator-shielded, pincer-routed, CRT-split, or
positive-measure escaping.

## The creative hypothesis

For `n=14`, split the certificate sheaf into:

```text
unit sheet:       six AP witness clocks, already rigid
gcd-2 boundary:  six doubled residues, each with unit/nonunit lift ambiguity
apex chart:      the self-antipodal residue 7
```

The odd `q=7` block should behave like solved unit dynamics.  The proof
pressure should concentrate in the 2-block and the apex chart.  The conjectural
terminal statement:

```text
no boundary monodromy class survives endpoint-owner peeling,
denominator shields, pincer closure, and CRT splitting.
```

This is close to the S579 apex-lift certificate sheaf, but S590 adds a sharper
boundary map: the exact way the unit clocks leave the sheaf under doubling.
Together with HYP-2128, this says the `n=14` worry surface has an additive
modular shell `2n-1=27` and a multiplicative boundary seam from unit clocks to
gcd-2 lifts.

## Assumption challenge

For this lens, vertices are not runners.  The useful vertices are orbit sheets,
boundary residues, quotient fibers, local certificates, endpoint owners, and
proof obligations.  The LRC predicate survives only when the labels that make a
local certificate glue are retained.

That is the sober version of the abstract picture:

```text
regular polygon = trivial monodromy;
bracelet        = order-2 monodromy;
nonabelian VT   = relator monodromy;
n=14 doubling   = boundary monodromy.
```
