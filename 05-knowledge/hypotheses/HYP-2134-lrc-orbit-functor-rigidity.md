---
id: HYP-2134
status: OPEN orbit-functor atlas + finite AP-clock audit; S590 isolates n=14 as statically rigid, reflectively paired, and dynamically broken in one 2-adic bit
source: codex-2026-06-03-S590
related: [HYP-2141, HYP-2140, HYP-2138, HYP-2137, HYP-2136, HYP-2135, HYP-2133, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2128, HYP-2127, HYP-2126, HYP-2125, HYP-2124, HYP-2123, HYP-2122, HYP-2121, HYP-2120, HYP-2118, HYP-2117, HYP-2116, HYP-2113, HYP-2112, HYP-1783, THM-401, THM-400, THM-385, THM-381]
---

# HYP-2134: LRC rigidity should be classified by orbit functors, not by one local/global axis

## Claim

The useful question is not only whether an object is locally or globally rigid.
It is:

```text
which orbit functor is acting,
which LRC predicate labels it preserves,
and which labels it forgets or permutes.
```

The same AP witness set

```text
U_n = (Z/n)* = {j : AP is lonely at j/n}
```

has several incompatible rigidity structures:

```text
static multiplicative  (Z/n)*        transitive for every n
dynamical doubling     <x -> 2x>     internal iff n is odd
reflective dihedral    j <-> -j      pairs antipodal pincer jaws
CRT factor             n -> prod p^a localizes defects by prime block
quotient transport     X -> Q(X)     empty/rigid/leaking/transport fibers
monodromy              loop in Q     trivial or label-permuting lift
isostatic              constraint graph critically rigid or leaking
spectral character     Fourier/Dirichlet blocks carry defect mass
fiber stiffness        boundary labels resist flips more than inner labels
source functor         delete/add marked source exactly
```

Thus "rigidity" is an action-specific property, not a scalar amount of
symmetry.  A local fixed point becomes proof-relevant only when it is a section
for the functor and its labels are natural under the action.  Global rigidity
is the natural cascade of those labels through the functor.  Projection defect
is the failure of that cascade after forgetting a label or after looping with
nontrivial monodromy.

This refines the slogan

```text
vertex-transitive trienerment <=> regular polygon point-set.
```

The slogan is true for the cyclic primitive functor.  General
vertex-transitive trienerment means local rooted profile plus the relevant
cascade law: cyclic rotation, dihedral bracelet, nonabelian Cayley relator,
source-root extension, CRT channel, or labelled LRC proof automaton.

## Evidence From S590

`04-computation/orbit_rigidity_functor_atlas_s590.py` audits the AP unit clock
under several actions.

At `n=14`:

```text
unit witnesses:              [1, 3, 5, 9, 11, 13]
static unit orbit:           [[1, 3, 5, 9, 11, 13]]
doubling x2 on units:        [[1], [3], [5], [9], [11], [13]]
reflection/antipodal pairs:  [[1, 13], [3, 11], [5, 9]]
CRT signature:
  mod 2: unit point [1], x2 collapses to 0
  mod 7: x2 cycles [[1, 2, 4], [3, 5, 6]]
```

So `n=14` is not "nonrigid."  It is:

```text
statically rigid        one unit orbit
reflectively rigid      three pincer pairs
dynamically broken      six singleton doubling fragments
CRT-localized           the defect is exactly the mod-2 collapse
```

The finite table for `n in {7,8,9,10,12,14,15,16,18,21,27,28,30}` shows the
same pattern: full unit multiplication is one orbit on the AP witness clock,
while even `n` makes doubling exit the unit fiber and fragment the dynamic
proof route.  Odd `n` keeps doubling internal, though not necessarily a single
cycle.

S590's Tournament Analysis uses proof lenses as vertices, not runners or arcs.
The pair observable is:

```text
(payload_retention, orbit_connectedness,
 defect_localization, proof_maturity, -cost)
```

The resulting transitive path is:

```text
labelled_source_cascade
> pincer_denominator_endpoint
> static_unit_clock
> CRT_factor_projection
> dihedral_reflection_pincer
> quotient_fiber_transport
> cyclic_polygon_action
> dynamical_doubling_action
> spectral_character_blocks
> dihedral_bracelet_action
> isostatic_constraint_graph
> raw_unmarked_shadow
```

There are zero directed 3-cycles.  The important flip is that cost-only
ranking would choose much coarser lenses; the proof ranking prefers lenses that
keep the observer/source, denominator, endpoint, pincer, or CRT labels.

## Rigidity Dictionary

```text
local fixed point
  source observer, straddle-pinch clock, endpoint owner, denominator shield

global cascade
  functorial transport of the retained label through isomorphism/action

projection defect
  mixed predicate values inside a forgotten fiber

monodromy defect
  a closed loop in the quotient lifts to a nontrivial permutation of labels

static rigidity
  full unit action keeps the AP witness set finite and transitive

dynamical rigidity
  doubling acts internally on the witness fiber; fails exactly at even n

reflective rigidity
  negation pairs witness clocks into antipodal pincer jaws

factor rigidity
  CRT splits the proof into prime-power channels and localizes the leak

fiber stiffness
  boundary-label flips are costly relative to inner flips in a fixed class
```

This also clarifies multiplication versus addition.  Multiplication by a unit
is a symmetry of the AP clock.  Multiplication by `2` is a symmetry only when
`n` is odd; at even `n`, it becomes a projection/collapse.  Addition supplies
fold and denominator gates (`a+b=D`), so additive rigidity is not the same as
multiplicative rigidity.  It is observer-coupled and must carry denominator
and endpoint labels.

## New Falsifiable Hypotheses

**H7 Monodromy rigidity.** Any residual LRC counterexample core must contain a
nontrivial label-monodromy loop in an observer-marked quotient.  If monodromy
is trivial, source, denominator, pincer, and endpoint labels should cascade to
a witness, a positive `Phi` escape, or a smaller labelled core.

**H8 Dihedral pincer rigidity.** Even-`n` doubling fragments are still paired
by `j <-> -j`.  A labelled pincer proof should discharge each antipodal pair
unless a CRT/endpoint owner loop exports a middle-state core.

**H9 Fiber-stiffness criterion.** Hard buckets are maximally stiff in the
observer fiber, not necessarily maximally symmetric.  The ratio of boundary
label flips to inner flips should predict which quotient fibers leak safe and
unsafe states.

**H10 Character leakage.** For `n=2q`, all non-principal orbit-defect mass is
carried by the 2-adic character block; the `q`-block is already prime/odd and
should be certifiable by existing clocks.

## Proof Route For n=14

Treat the AP/`V*` boundary as a functor-atlas problem:

1. Use HYP-2124 static unit rigidity to reduce to the six unit witnesses.
2. Pair them by HYP-2123/HYP-2134 dihedral pincers:
   `(1,13)`, `(3,11)`, `(5,9)`.
3. Use CRT to split `14=2*7`.  The mod-7 coordinate has ordinary odd-prime
   doubling cycles; the mod-2 coordinate is the only collapse.
4. Push source and denominator labels through HYP-2127's rigidity closure:
   observer-source, denominator shields, endpoint owners, pincer ledgers, and
   `L/M/R` middle states.
5. Any residual core should be audited specifically for monodromy or fiber
   stiffness.  A raw unmarked tournament class, score sequence, or balanced
   energy bucket is too coarse unless it preserves those labels.

The desired terminal theorem is:

```text
Every n=14 residual either has trivial label monodromy and closes, or has a
nontrivial monodromy/stiffness witness that can be isolated in the 2-block and
discharged by endpoint/Phi positivity.
```

## Assumption Challenge

Do not assume the orbit vertices are runners.  The useful vertex sets here are
proof lenses, witness residues, CRT factors, antipodal jaws, quotient fibers,
source-root states, endpoint-owner clauses, character blocks, constraint
nodes, and monodromy loops.

The chosen quotient preserves the predicate "this labelled local fixed point
can be transported to an LRC certificate."  It destroys exact runner phases and
unlabelled tournament detail on purpose.  If that destruction creates mixed
safe/unsafe fibers, the quotient must be lifted.

## See

`04-computation/orbit_rigidity_functor_atlas_s590.py`,
`05-knowledge/results/orbit_rigidity_functor_atlas_s590.out`,
`07-reflections/lrc-orbit-functor-rigidity-s590.md`,
HYP-2141, HYP-2140, HYP-2138, HYP-2137, HYP-2136, HYP-2135, HYP-2133, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2128, HYP-2127, HYP-2126, HYP-2125, HYP-2124, HYP-2123, HYP-2122,
HYP-2121, HYP-2120, HYP-2118, HYP-2117, HYP-2116, HYP-2113,
HYP-2112, HYP-1783, THM-401, THM-400, THM-385, THM-381.
