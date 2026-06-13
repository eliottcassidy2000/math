---
id: HYP-2127
status: OPEN synthesis with exact S590 orbit-sheaf audits; bridged to HYP-2134 orbit-functor atlas; proposes boundary/gluing rigidity for n=14
source: codex-2026-06-03-S590
related: [HYP-2134, HYP-2133, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2128, HYP-2126, HYP-2125, HYP-2124, HYP-2123, HYP-2122, HYP-2121, HYP-2101, HYP-2020, HYP-1977, THM-401, THM-400, THM-381]
---

# HYP-2127: rigidity is an orbit sheaf with monodromy defects at forgotten labels

## Claim

Rigidity in the LRC/tournament interface is best modelled as a sheaf over orbit
quotients:

```text
base       = an orbit quotient, such as unit clocks, rooted tournament classes,
             dihedral local profiles, endpoint charts, or proof-obligation cells
stalk      = labels needed by the predicate
restriction= transport/deletion/doubling/CRT/endpoint/pincer maps
section    = a local certificate with all required labels
rigid      = the section glues with trivial or controlled monodromy
defective  = forgotten fibers create monodromy or mixed predicate labels
```

The short slogan is:

```text
rigidity = a section over an orbit quotient that glues without losing the
predicate labels;
failure = monodromy or projection defect in a forgotten fiber.
```

This fuses the two HYP-2126 strands.  One strand says to see every object as an
orbit and compare rigidity types on that orbit.  The other says proof-relevant
rigidity is labelled fixed-point propagation.  HYP-2127 says these are the same
problem: a labelled orbit section either glues, hits a boundary, or exposes a
projection defect.

HYP-2134 complements this sheaf framing by classifying the same AP witness
payload under explicit orbit functors: static units, doubling, reflection, CRT,
quotient transport, monodromy, stiffness, and source deletion.  HYP-2133 is the
labelled fixed-point cascade side of the same story, while HYP-2131 is the
prime-2 localization atlas that names the bad place for the same defect.
HYP-2132/THM-401 pins the additive pair-sum modulus `2n-1` that this sheaf uses
as its odd shell.

## Types Of Rigidity

```text
local rigidity        stabilizer fixes a labelled basepoint or endpoint
global rigidity       a transitive action transports that label everywhere
static rigidity       a symmetry group preserves the witness orbit
dynamic rigidity      an operator, such as doubling, preserves the witness orbit
boundary rigidity     an operator leaves the orbit but lands in a labelled stratum
quotient rigidity     fibers of a forgetful map are label-pure
gluing rigidity       local certificates agree on overlaps
monodromy rigidity    transport around a loop returns the same certificate
spectral rigidity     character blocks separate principal and defect modes
isostatic rigidity    active constraints match orbit degrees of freedom
pincer rigidity       opposing certificate fronts meet or export a labelled core
automaton rigidity    L/M/R states force every side change through a middle chart
```

The user's local/global split is therefore the first two entries in a larger
orbit-sheaf menu.

## Evidence From S590

`04-computation/lrc_orbit_sheaf_monodromy_s590.py` gives three exact toy audits.

### 1. Doubling Is A Boundary Morphism At Even n

For AP, HYP-2124 proves the lonely clock residues are exactly the unit orbit:

```text
W_n = (Z/n)^*.
```

S590 tests the doubling map on this witness sheaf for `5 <= n <= 18`.

For odd `n`, doubling is an automorphism of the unit sheaf:

```text
n=5,7,9,11,13,15,17: kept_in_unit_sheaf = |(Z/n)^*|
```

For even `n`, no unit witness remains in the unit sheaf after doubling:

```text
n=6,8,10,12,14,16,18: kept_in_unit_sheaf = 0,
                      boundary_count = phi(n),
                      image_gcd_hist = {2: phi(n)}.
```

At `n=14` the boundary map is especially clean:

```text
1  ->  2,  lifts {1 unit, 8 nonunit}
3  ->  6,  lifts {3 unit, 10 nonunit}
5  -> 10,  lifts {5 unit, 12 nonunit}
9  ->  4,  lifts {9 unit, 2 nonunit}
11 ->  8,  lifts {11 unit, 4 nonunit}
13 -> 12,  lifts {13 unit, 6 nonunit}
```

So HYP-2126's "fragmentation" sharpens to:

```text
doubling is not an internal dynamical symmetry at even n;
it is a boundary morphism from the unit witness sheaf to the gcd-2 stratum.
```

For `n=14`, every doubled witness has exactly one unit lift and one gcd-2
nonunit lift.  That lift ambiguity is the 2-block seam in sheaf language.

### 2. Dihedral Point-Sets Have Bracelet Monodromy

The dihedral point-set audit through `N <= 18` reproduces HYP-2125:

```text
dihedral_VT = 83
regular     = 31
bracelet    = 52
gap_period_hist = {1: 31, 2: 52}
```

The first non-regular object is:

```text
N=6, P=(0,1,3,4), gaps=(1,2,1,2).
```

There is no period greater than `2` in the audited dihedral-VT range.  Thus the
dihedral local sheaf glues in exactly two ways here:

```text
trivial cyclic monodromy       regular polygon
order-2 block monodromy        imprimitive bracelet
```

### 3. Rooted Tournament Quotients Mix Labels

The rooted tournament audit through `n=6` treats quotient lenses as base
objects and asks whether labels are pure on fibers.

At `n=6`:

```text
full_rooted           unique=296 collisions=  0 max_fiber=  1
unrooted_plus_score   unique=196 collisions=100 max_fiber=  4
split_no_cross        unique= 36 collisions=260 max_fiber= 64
unrooted              unique= 56 collisions=240 max_fiber=  6
delete_parent         unique= 12 collisions=284 max_fiber= 32
score_sequence        unique= 22 collisions=274 max_fiber= 72
root_score            unique=  6 collisions=290 max_fiber= 88
```

Source and sink labels are pure for score-aware lenses, but the finer rigidity
labels are not.  For example:

```text
unrooted_plus_score: parent mixed on 74 fibers
split_no_cross:      fixed/orbit_size mixed on 8 fibers, parent mixed on 12
root_score:          fixed/orbit_size mixed on 4 fibers, parent mixed on 6
```

This is the quotient-sheaf warning.  A quotient may preserve the obvious
source score but still lose the fixed-root, parent, orbit-size, cross-edge, or
endpoint-owner labels that a proof needs.

Tournament Analysis over rigidity lenses uses defect tuple

```text
(boundary_loss, projection_collisions, max_fiber, label_mix)
```

and is transitive in S590:

```text
cyclic_polygon_spine
-> source_root_sheaf
-> static_unit_sheaf
-> dihedral_bracelet_monodromy
-> unmarked_shadow
-> split_profile_no_cross
-> root_score_shadow
-> doubling_even_boundary_n14
```

The last placement is deliberate: the chosen switch treats boundary loss as a
harder defect than quotient collision, because a proof riding the doubling
operator has left its witness sheaf entirely.

## Creative Subhypotheses

### H7: Boundary Monodromy At Even n

For even `n`, the doubling map is a boundary morphism:

```text
(Z/n)^* -> {r mod n : gcd(r,n)=2}.
```

For `n=2q` with `q` odd, every doubled unit residue has one unit lift and one
gcd-2 nonunit lift.  The n=14 residual lives in that lift ambiguity, not in the
unit orbit itself.

### H8: Gluing Purity Criterion

A tournament or LRC quotient is proof-safe only if every required predicate
label is pure on its fibers.  The label list should include:

```text
observer-source state, fixed-root status, parent/deletion class, augmentation
grade, endpoint owner, pair denominator, shield blocker, pincer side, L/M/R
middle state, CRT lane, and p-adic depth.
```

If any one of these mixes, the quotient is a cache layer, not a certificate
layer.

### H9: CRT Sheaf Splitting For n=2q

For `n=2q`, the orbit sheaf should split into:

```text
q-block: solved odd/unit dynamics
2-block: boundary lift ambiguity
```

The proof of n=14 should first show the `q=7` block glues, then prove every
2-block boundary monodromy either transports to a labelled pincer certificate
or opens positive measure.

### H10: Monodromy-To-Escape

Any nontrivial monodromy loop in the certificate sheaf must carry an endpoint
owner or private pincer pivot.  Therefore failed gluing is not a counterexample
core by itself; it should route to:

```text
failed gluing -> cover circuit -> endpoint owner/private pivot -> peel ->
positive measure.
```

This is the HYP-2101 apex-lift certificate-sheaf dream, but with the S590
doubling boundary map as the explicit `n=14` seam.

## n=14 Proof Reading

The AP unit clocks are statically rigid:

```text
{1,3,5,9,11,13}/14.
```

Doubling sends them to the six gcd-2 boundary residues:

```text
{2,4,6,8,10,12}/14.
```

The missing nonunit residue `7/14` is the self-antipodal apex.  Thus the n=14
proof should not ask whether the unit orbit is rigid; it is.  It should ask:

```text
Can every gcd-2 boundary lift and the apex chart be glued to a labelled
endpoint/pincer certificate or a positive-measure escape?
```

Equivalently:

```text
No nontrivial boundary monodromy class survives after CRT splitting,
endpoint-owner peeling, denominator shields, and pincer closure.
```

After rebasing over HYP-2128, the two faces are clearer.  HYP-2128 gives the
additive `+2`/triangular face: `2n-1` is the odd-square root of `8*C(n,2)+1`.
HYP-2127 gives the multiplicative `x2`/boundary face: at even `n`, doubling
leaves the unit witness sheaf for gcd-2 boundary residues.  For `n=14`, the
proof should use both: the `2n-1=27` modular shell is the additive face, while
the gcd-2 lift ambiguity is the dynamical boundary face.

## Assumption Challenge

Do not assume the vertices are runners.  In this hypothesis the useful vertices
can be:

```text
unit-clock sheets, doubled boundary residues, CRT blocks, rooted tournament
fibers, local trienerment profiles, dihedral block systems, endpoint owners,
pair denominators, pincer middle states, or proof obligations.
```

The quotient preserves the LRC predicate only if it preserves the labels needed
to glue a certificate.  It destroys fine phase geometry and raw runner identity
on purpose; when that destruction mixes predicate labels, the quotient must be
lifted.

## See

`04-computation/lrc_orbit_sheaf_monodromy_s590.py`,
`05-knowledge/results/lrc_orbit_sheaf_monodromy_s590.out`,
`07-reflections/lrc-orbit-sheaf-monodromy-rigidity-s590.md`,
HYP-2134, HYP-2133, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2128, HYP-2126, HYP-2125, HYP-2124, HYP-2123, HYP-2122, HYP-2101,
HYP-2020, HYP-1977, THM-401, THM-400, THM-381.
