# LRC(14): the proof object is a ramified coset-cover system

**Date:** 2026-07-17  
**Scope:** structural synthesis after the independently replayed exact
common-scale H6 closure through scale 31.  This is a research map, not a claim
that LRC(14), the H5 bank, or the global `n=12` sporadic branch is closed.

The long common-scale computation has stopped looking like a list of unrelated
scale accidents.  The same object keeps reappearing under different
projections.  It is not fundamentally a tournament on runners, and it is not
fundamentally a set cover on individual sheets.  It is a **coloured,
owner-covariant coset-cover constraint system**, with a prime-power flag on the
sheet circle and one globally shared unit choice at each provider.

## 1. The object before quotienting

Fix a common scale `c`.  A scalar row consists of six distinct labels
`a_i in F_13^*` and six effective orders `D_i | c`, subject to the hereditary
leave-one-out lcm grammar.  The unit fibre at provider `i` is

```text
U_i = (Z/D_i Z)^*                       (a singleton when D_i=1).
```

For an owner label `y`, a provider and a unit determine a literal sheet mask

```text
M_y(i,u_i) subset Z/cZ.
```

All six owners see the **same** unit word `u=(u_1,...,u_6)`.  Define the owner
obligation

```text
O_y = {u in Product_i U_i : Union_i M_y(i,u_i) = Z/cZ}.
```

The common-sheet question is therefore

```text
Intersection_{y in support} O_y != empty.                         (1)
```

This formulation separates four logically different losses which had often
been compressed into one computation:

| layer | retained predicate | information discarded |
|---|---|---|
| scalar capacity | `sum_i |M_y(i,u_i)| >= c` | all overlap and unit correlation |
| quotient flag/coset cover | coverage in fibres of `Z/c -> Z/d` | affine positions inside a fibre |
| owner-local reachable bank | `O_y != empty` | the requirement that all owners use one unit word |
| owner-obligation nerve | intersections of the `O_y` | exact metric-component placement |
| literal component recursion | the actual common-sheet/metric terminal | almost nothing |

The first three layers are a **pre-nerve**.  If one `O_y` is empty, the global
intersection dies before pairwise or higher owner compatibility matters.  If
all are nonempty, their intersection nerve—not their cardinality ranking—is
the next object.

This also clarifies scope.  The system in (1) is the proof object of the
primitive proper AP-centred common-sheet H6 reduction.  It is not by itself a
replacement for the non-AP-centred, deep-sheet, H5, or metric-lift branches.

## 2. Boolean group rings explain the two recurring proof styles

Encode a sheet mask by

```text
m_y(i,u_i) = Sum_{s in M_y(i,u_i)} X^s
```

in the Boolean group algebra of `Z/cZ`.  Coverage says that the Boolean sum of
the six mask polynomials is the all-ones polynomial.  Two regimes follow.

### Tight regime: characters see exact partitions

If scalar capacity is exactly `c`, coverage forces the masks to be disjoint,
so the Boolean identity lifts to a coefficientwise identity over the integers.
Every nontrivial additive character must then annihilate the sum.  This is the
mechanism behind the scale-21 cubic-character obstruction: the unit-pair
coefficients in the Eisenstein character ring cannot sum to zero.  Scale seven
uses the same philosophy one layer earlier, through row products and a
quadratic power sum.

### Slack regime: quotient fibres measure overlap debt

When scalar capacity exceeds `c`, characters alone can hide collisions.  Push
the masks through a quotient `Z/c -> Z/d` and retain, for each fibre, either
the exact incidence pattern or its capacity saturated at the fibre size.
Coverage is then a coloured capacitated hypergraph problem.

- At `c=24`, the mod-three sheet classes and cubic label cosets expose the
  exact five-sheet debt at the two failed owners.
- At `c=25`, five owner-local sheet cosets turn the two order-five masks into
  thick fibres and the order-25 masks into rigid needles; a four-incidence tax
  gives the sharp bounds 22 and 21.  The name of the self coset must be chosen
  owner by owner.
- At `c=27`, the nested flags are `Z/27 -> Z/3` and `Z/27 -> Z/9`.
  Order-three anchors close 420 of the 450 scalar rows through nine-sheet
  fibres; order-nine anchors close the remaining 30 through three-sheet
  fibres, while order-27 masks are point needles.  A sound upper relaxation
  lets each nonanchor needle optimize independently, so the information loss
  can only enlarge the union.  Independently, the saturated `Z/27 -> Z/9`
  flag treats orders 1, 3, and 9 as unions of whole three-point fibres and
  order-27 masks as transversals; it makes every order-27 owner impossible.
  Three literal-bank implementations agree on every exact headline census.
- At `c=30`, no single prime-power chain is the right carrier.  The first
  flag retains orders dividing six and reduces 54,050 scalar rows to 120.
  Every survivor has lost all of those anchors.  The proof must rotate to the
  incomparable quotient retaining orders dividing ten; its 720 owner bounds
  are all at most 28.  Thus the flag calculus is a **cover of divisor
  directions**, not necessarily a nested chain.

Thus “character proof” and “nerve proof” are not competing methods.  They are
the zero-slack and positive-slack faces of the same group-ring cover system.

## 3. The Fourier dual is the same ramification data

THM-1091 makes the finite-sheet duality exact.  If a mask on `Z/cZ` is
invariant under translation by `D`, then its Fourier transform is supported
on the annihilator

```text
{r in Z/cZ : (c/D) divides r}.
```

Whole quotient fibres in the primal picture are sparse frequency support in
the dual picture.  An exact-capacity cover is a partition, so every nonzero
character sum cancels.  With slack, the cover error

```text
e(t) = sum_i 1_(M_i)(t) - 1 >= 0
```

has a Parseval energy equal to the total nonzero character energy.  The
zero-mode capacity and the nonzero-mode cancellation are therefore not two
different invariants; they are the two halves of one group-ring identity.

This explains why the complementary flags at `c=30` work.  Period-six
anchors live on modes divisible by five; period-ten anchors live on modes
divisible by three.  The residual of the first carrier is not expected to be
visible to a refinement of the same carrier: it has selected configurations
whose surviving information lies in a transverse annihilator.  Switching
from six to ten is a change of chart on the same ramified object.

THM-1092 supplies the continuous analogue.  Character orthogonality selects
the resonance lattice `sum_i a_i n_i=0`, and a dyadic-shell argument proves
the full-support reciprocal series converges absolutely.  Thus the continuous
danger-intersection moments and the finite-sheet cover identities share an
honest resonance/annihilator language, not merely a numerical analogy.  What
is still missing is quantitative cancellation: absolute convergence licenses
regrouping, but does not make the higher-support moments small.

## 4. How the carrier sharpened across the scale history

The historical sequence is now structurally legible.

- `c=1,2,3`: literal components and affine matching codes were still needed;
  the metric recursion was the faithful terminal object.
- `c=4`: an edge-coloured `K6` and its affine unit words appeared before the
  remaining component recursion.
- `c=5`: owner triples, not pair edges, became the first terminal hypergraph.
- `c=6`: owner obligations were affine four-flats with an octahedral nerve.
- `c=7`: the correct object changed to a multiplicative Latin row; product
  characters and square sums killed the hard cosets.
- `c=8,9,10`: the terminal objects were respectively a `K_{3,3}` nerve, sparse
  `3K_2/2K_2` nerves, and a projective triangular prism.
- `c=11,12`: owner orthogonality became exact; the c=12 terminal quotient is
  now small enough for a kernel-checked statement.
- `c=14,...,22`: scalar gates plus owner-local reachable banks usually killed
  a row before a global nerve was born.  The important invariant was the
  labelled vector of maximum unions, not an ordering of its entries.
- prime scales: directed Cayley graphs, spectra, and multiplicative closure
  replaced sheet enumeration at the scalar layer.
- `c=21,24,25,27`: ramification made the missing structure explicit—gcd
  strata, multiplicative character classes, and quotient flags on sheets.
- `c=28`: the four thick fibres of `Z/28 -> Z/4` dominate the transverse
  seven-adic geometry.  Even after every nonanchor toothpick is allowed to
  choose its best contribution independently, at most two of six owner
  obligations remain live.  The full `Z/4 x Z/7` flag does not improve the
  terminal threshold on any scalar survivor.
- `c=30`: a nonnested two-chart atlas is necessary.  `Z/6` finds the thin
  120-row residual and `Z/10` kills it uniformly; the exact DP shows the
  relaxed terminal bounds lose at most two sheets and never cross threshold.

The recursion is from large literal geometry to a smaller carrier which still
preserves the current predicate.  No single carrier is universally faithful.
Changing carrier when the current obstruction becomes blind is part of the
proof, not a failure of uniformity.

## 5. What the tournament program actually taught us

Tournament analysis was useful because it repeatedly forced the question
“what are the vertices, what is the switch, and what survives the quotient?”
Its strongest conclusion here is negative but precise.

The completed owner tournaments at scales 15, 20, 21, 22, 24, 25, 27, and 30 are
overwhelmingly or universally transitive under natural lexicographic gauges.
That stability does **not** mean the proof object is ordered.  It means that
ranking owners by `(feasible, maximum union, capacity, bank size)` erases the
absolute coverage threshold and the mask intersections.  Rows with different
obstructions collapse to the same score word `(0,1,2,3,4,5)`.

The proof-facing pair relations have instead been partial or non-tournament
objects:

- directed forbidden-ratio Cayley graphs;
- disjointness and intersection graphs of owner obligations;
- projective cycles and prisms;
- the weighted `3K_4` missing-owner graph at scale 24;
- character-labelled ratio graphs.

Forcing these relations to complete tournaments is lossy.  The useful
tournament output is therefore an **audit sidecar**: score histogram, cycles,
SCCs, edge flips, ties, and Hamiltonian-path count certify that two
implementations used the same gauge.  The Hamiltonian tie path is a canonical
serialization device, not the obstruction.

Three tournament-related experiments remain worth doing, with the threshold
data kept beside them:

1. **Forced-edge completion spectrum.**  Orient only comparisons invariant
   under every admissible unit word, enumerate the completions of the remaining
   partial relation, and ask which cycle/SCC fingerprints survive every
   completion.  This distinguishes proof-bearing edges from tie-path artefacts.
2. **Switch graph on witnesses.**  Use local unit words as vertices and join
   two when one provider-unit switch transforms one into the other.  Colour an
   edge by the owners gained and lost.  Directed cycles can then encode genuine
   owner correlation, unlike a scalar ranking tournament.
3. **Mode/obligation vertices.**  At tight scales use nontrivial character
   modes; at slack scales use quotient fibres or owner obligations.  Orient by
   signed character contribution or one-sided overlap debt, but retain the
   magnitude and zero/tie label.  The binary orientation alone is not enough.

At scale 30 the quotient flags `{6,10,15}` themselves form a useful diagnostic
tournament, but all 720 completions are transitive.  The proof is the absolute
fact `U_10<30`, not the edge `10 -> 6` or the score word.  This is the cleanest
case yet where tournament vertices can be proof methods while the actual
carrier remains owner-coloured fibres with numerical capacities.

These challenge the old assumption that vertices must be runners or arcs.
Gaps, sheet fibres, residues, wall crossings, Fourier modes, and proof
obligations are all legitimate vertices when they preserve the predicate being
tested.

## 6. The missing general theorem

The scale-by-scale closures have exposed the shape of a uniform statement but
have not proved it.  The desired theorem should be recursive in the
prime-power flag of `c`.

1. Factor `c` and enumerate only hereditary valuation-colour words.
2. Apply the exact ratio-cardinality table to eliminate scalar-impossible
   supports.
3. Choose a small **atlas of quotient flags** on which selected low-order
   masks are unions of fibres and the remaining masks have bounded transversal
   multiplicity.  The flags need not form a chain; switch charts when one
   relaxation selects an anchor-free residual, as at scale 30.
4. If total capacity is tight, pass to character identities in the integer
   group ring.
5. If there is slack, saturate fibre capacities and prove an overlap-debt
   bound.
6. If every owner survives locally, form the nerve of the globally covariant
   unit obligations.
7. Invoke literal metric components only after all smaller carriers fail.

Part of this infrastructure now exists.  `LRCPreNerveProjection` formalizes
the one-way global-to-owner implication; `LRCRamifiedCosetCover` formalizes
anchor/nonanchor upper bounds, saturated fibre flags, strict-cardinality
deficits, and complete-fibre presentations; THM-1091 supplies the mathematical
finite-character identity.  The remaining formal suppliers are sharper:

- an owner-covariance lemma which explicitly transports sheet classes and
  prevents the scale-25 “global self coset” mistake;
- a kernel-checked character obstruction/energy lemma for tight Boolean
  covers;
- a multi-flag atlas theorem which permits nonnested residual switching and
  composes its local deficits;
- certificate-import interfaces connecting complete raw-bank enumeration to
  the generic Lean consumers without hiding it behind `native_decide`; and
- a quantitative chart-selection theorem, ideally dual in fibre capacity and
  annihilator support, that closes many valuation grammars uniformly.

The common-scale H6 bank is only one face of LRC(14), but it has now supplied a
candidate underlying language.  The next real gain is proving enough of this
ramified coset-cover calculus that many
scales close from their valuation grammar at once, then transporting the same
language to the smooth H5 and non-AP/deep-sheet residuals.
