# LRC(14): free refinement versus global choice

The long history in this repository has not merely accumulated lemmas. It has
repeatedly corrected what we thought the underlying object was. The newest
prime-seven and two-sheet results make one common distinction visible:

> some dynamics only refine the description for free; progress requires
> quotienting those dynamics out and then making a global choice among the
> genuinely interacting pieces.

This is why two initially plausible opening targets now have different honest
statuses. A uniform good denominator `q<=25` is false: divisor loading can
kill every prescribed bounded denominator ladder. Uniform emptiness of the n=12
sporadic branch remains open: the branch is finite in principle and sharply
sheet-stratified, but its remaining deep covers are not captured by a bounded
denominator, raw height, scalar measure, or one local phase packet.

**2026-07-15 update.** THM-815/816 close both Hamming-four alternatives left
by THM-810, so every AP-centred residue-preserving packet through radius four
is uniformly loose.  THM-815 has independent component-discrepancy and
collar/doubling certificates. The same exact component/comb discrepancy makes
the scale-one radii five and six finite-decidable, and loses coercivity at
seven remaining combs; arbitrary higher-radius deck ramification is still
unclassified. THM-820 makes the Hamming-five finiteness operational: after
intersecting with THM-815, only the doubling box
`(146,292,584,1168,2336)` and an exceptional top-four box with
`x<=146,v<=1986,max<=7944` remain. Its tournament-identical centre-9/centre-22
cycles show that the global choice must act on residual intervals with legal
continuations, not on the collar graph.
THM-817 decomposes every disconnected return into signed maximum-speed cells
and sharpens the deep selector to `2c_E N_R+2W-2g`; an exact family with
`N_R=Theta(B)` proves that global choice must retain those satellites rather
than assume them away.  The refined common rule is: quotient free refinement,
then recurse on the full residual interval union with labelled comb and
endpoint-owner incidence.

## What each historical viewpoint really retained

| viewpoint | exact gain | information it loses |
|---|---|---|
| rational good periods | explicit witness and signed residue blocker deck | divisor loading can kill any fixed bounded ladder |
| covering/divisor pins | routes non-covering families and forces sheet structure | does not identify a compatible global time |
| capped envelope / safe measure | quantitative room and per-core finite tails | raw component count grows under scale |
| binding denominator `p/(13s)` | exposes shallow versus deep sheet fibres | a residue snapshot does not encode persistence over the loose set |
| folded diamond | lossless two-odd-runner predicate | scalar diamond measure forgets component incidence |
| phase pigeonhole | height-free safe mass and returns | an arbitrary anchor can be locally trapped |
| prime-seven tokens | exact sheet coverage polynomial | instantaneous tokens omit endpoint schedule and metric base |
| centered Beatty/continued fractions | exact owner word, ties, and parity cocycle | the word alone omits token fibre and core incidence |
| metagraph / `HP(T)/Aut(T)` | intrinsic finite tiling fibres | iso nodes and masks omit owner assignment and continuation |
| wall-event tournament | chronological Hamiltonian path | finest resolution counts free same-owner subdivision as complexity |
| ownership hypergraph | simultaneous compatibility and transversal deficit | needs metric endpoints and recursion sidecars |
| dyadic deletion tree | exact descent through imprimitive two-sheet seams | does not choose the terminal safe component |

The recurring error pattern is now unmistakable. Raw height, component count,
wall count, and phase-cell depth can all be inflated by a free operation, while
a prescribed denominator ladder can be erased by divisor loading. Dilation
proliferates components; a very fast owner inflates wall count. In THM-789's
exact row, refining one phase cell adds coordinates while staying inside the
same bad diamond. These are resolution variables, not compactness coordinates.

## The prime-seven correction

THM-784 makes the failure exact. On `(5/16,7/20)`, the seven slow owners
`{1,2,3,4,5,8,10}` already form the complete token rainbow. The fastest owner
`560N+1` contributes `21N` consecutive covered walls without changing the
slow state. The event tournament simply becomes a longer transitive order.

THM-788 supplies the first quotient. Write a run as

```text
E_0, V_1, E_1, ..., V_A, E_A.
```

Each `E_j` is a fastest-owner block absorbing empty periods; each `V_j` is an
ordered zero-sum visitor packet. THM-794 proves that this quotient is still
noncompact: for `F=49H+1`, the full eight-owner packet repeats through `H-1`
active periods and `8H-8` genuine owner switches even though
`ceil(F/(F-7))=2`. The packet returns by common diagonal sheet translation,
so its raw token displacement is zero in the reduced deck `F_7^8/Delta` while
the metric base advances. Consequently `A`, switch count, and the proposed
universal extent bound are not proof costs. THM-788's exact conditional
conversions remain valid once some later invariant controls the normalized
return dynamics.

This also clarifies the role of Tournament Analysis. Chronological comparison
gives a transitive tournament with one Hamiltonian path, zero cycles, and
singleton SCCs. Those fingerprints are correct but nearly content-free. The
proof carrier is the decorated path: owner labels, absorbed block lengths,
absolute slow scale, zero-sum visitor hyperedges, ordered collision state, and
reduced return holonomy. Changing vertices from runners to wall events was not
enough; maximal free blocks and central packet cycles must both be contracted.

The remaining prime-seven theorem is therefore:

1. recognize legal packet loops in the centered-Beatty collision transducer
   and quotient their diagonal deck holonomy;
2. prove an exit from every remaining normalized SCC, retaining varying
   indices and possible order flips—MISTAKE-148 rules out the fixed-index
   de-phase shortcut;
3. lift that exit to metric time and intersect it with the actual core-safe
   component.

## The two-sheet correction

THM-789 strengthens the local metric substrate substantially. Symmetrization
doubles the guaranteed return mass and improves the component-width floor by
a factor of four. The pointwise thickness tax

```text
||wt||+(w/B)(phi_U(t)-1/13)<=2/13
```

turns safe depth into an explicit odd-runner penalty, while the erosion

```text
E_U subset H_(x,y) minus R_U
```

is the exact form of what hypothetical tightness must maintain.

THM-797 adds a genuinely global arithmetic selector. An odd divisor `q` of
one exception exposes a balanced-residue acceptance shell on every deep unit
class. At `q=13`, containment forces the folded ten-core support to be full,
or to omit only the other exception's class. Thus small-support and misaligned
branches die uniformly. Its sharp aligned survivor also explains the residual:
the deepest thirteenths can all be trapped while `7/33` escapes. The object is
the whole owner-labelled family of deep components, not a single prime grid.

But the exact core

```text
U={1,2,3,5,7,8,9,10,11,12}
```

shows the limitation of every local upgrade. At `t_0=4/17`, the entire natural
return set `(-1/858,1/858)` is contained in the folded diamond for `(13,9)`.
Same-cell symmetrization, extra x/y/a/b phase coordinates, and the ordinary
Lipschitz interval all remain trapped. At `14/19`, however, the same core is
deep and outside the diamond.

So the missing two-sheet theorem is not “make the packet wider.” It is:

> choose globally among the deep components and prove that at least one has
> positive escape margin after its full return erosion.

Here a useful tournament may take deep components as vertices and orient by
maximum escape margin, with left endpoint as a tie gauge. The tournament is
again transitive; its purpose is to retain alternatives, not manufacture
cycles. Its bare isomorphism class is useless unless component intervals,
owners, and margins remain attached.

## A unified underlying object

Both residuals are faces of one owner-labelled metric incidence skew product:

```text
metric base component / deep component
    -> exact endpoint or phase schedule
        -> owner-labelled sheet/token fibre
            -> coverage observation
                -> recursive quotient or next-component choice.
```

The object must retain:

- normalized scale and one absolute metric unit;
- strict/open versus closed endpoint conventions;
- owner labels and simultaneous-event blocks;
- sheet/token assignment and inverse steps;
- component locations, widths, and escape margins;
- the free-refinement equivalence relation;
- the next proof obligation after a transition.

Every successful quotient in the repository preserves a named predicate and
states what it destroys. Every recurrent overclaim skipped one of those two
sentences.

## Recursive program from here

The current frontier has four honest branches:

1. **Shallow n=12:** THM-795/800/804/806/810/815/816 close the full
   AP-centred Hamming-one through Hamming-four stars. Exhaust THM-815/820's
   finite scale-one Hamming-five boxes and Hamming-six tree, classify the
   arbitrary-scale five-/six-owner deck interfaces, and replace the failed
   discrepancy coefficient at the seven-comb wall by a decreasing
   height/owner potential.
2. **Deep two-sheet:** prove global erosion noncontainment by selecting a deep
   component, or force the dyadic seam tree into the existing finite bases.
3. **Prime-seven r=8:** quotient zero reduced-holonomy packet loops, prove an
   exit between normalized collision SCCs, and then use core incidence; do not
   count raw walls, active periods, or raw switches.
4. **Higher sheets / scale-normal families:** retain colour ownership and
   ramification through descent, rather than collapsing to residue counts.

The latest concurrent metagraph work sharpens this diagnosis. THM-785/787 give
exact `C3`/`E4` flow coordinates, while the proved blue-parity theorem currently
filed as `THM-790-blue-parity-law-proved.md` derives the half-tiling count, the
mod-16 parity split, and the fact that the transitive pipe drains through the
two path legs. These are real recursive laws on the finite tournament fibre,
not merely pictures. They still do not reconstruct which metric component,
owner assignment, or token stalk an LRC continuation occupies. The atlas is
therefore a strong fibre coordinate and rigidity detector, but remains
diagnostic until a pullback theorem attaches that metric LRC stalk.

The new black-flow normalization audit makes the analogy quantitative. At
`n=7`, raw boundary current is `2798` pure-black-to-mixed versus `1254` in the
opposite direction, but division by source black-mask mass gives rates
`10.99%` versus `18.45%` and reverses the arrow. Line multiplicity, node
support, and source-fibre density are different gauges. A flow statement must
name which measure it preserves, just as an exit statement must distinguish
raw wall count from metric extent and active interaction.

The former namespace collision is repaired: the later Hamiltonian-path
companion is now `THM-791-H-companion-laws-to-the-transitivity-flow.md`, while
the earlier reserved THM-790 remains the generally proved blue-parity theorem.
Their scopes remain separate.

The common strategic rule is now precise:

> quotient free refinement first; then recurse on globally competing,
> owner-labelled metric components.
