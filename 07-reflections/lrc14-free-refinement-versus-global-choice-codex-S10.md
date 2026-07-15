# LRC(14): free refinement versus global choice

The long history in this repository has not merely accumulated lemmas. It has
repeatedly corrected what we thought the underlying object was. The newest
prime-seven and two-sheet results make one common distinction visible:

> some dynamics only refine the description for free; progress requires
> quotienting those dynamics out and then making a global choice among the
> genuinely interacting pieces.

This is why two initially plausible opening targets now have different honest
statuses. A uniform good denominator `q<=25` is false: scale changes move the
first witness beyond every fixed raw ladder. Uniform emptiness of the n=12
sporadic branch remains open: the branch is finite in principle and sharply
sheet-stratified, but its remaining deep covers are not captured by a bounded
denominator, raw height, scalar measure, or one local phase packet.

## What each historical viewpoint really retained

| viewpoint | exact gain | information it loses |
|---|---|---|
| rational good periods | explicit witness and signed residue blocker deck | witness denominators scale with dilation |
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

The recurring error pattern is now unmistakable. Raw height, denominator,
component count, wall count, and phase-cell depth can all be inflated by a
free operation. Dilation inflates denominators and components. A very fast
owner inflates wall count. Refining a phase cell adds coordinates while staying
inside the same bad diamond. These are resolution variables, not compactness
coordinates.

## The prime-seven correction

THM-784 makes the failure exact. On `(5/16,7/20)`, the seven slow owners
`{1,2,3,4,5,8,10}` already form the complete token rainbow. The fastest owner
`560N+1` contributes `21N` consecutive covered walls without changing the
slow state. The event tournament simply becomes a longer transitive order.

THM-788 supplies the right quotient. Write a run as

```text
E_0, V_1, E_1, ..., V_A, E_A.
```

Each `E_j` is a fastest-owner block absorbing empty periods; each `V_j` is an
ordered zero-sum visitor packet. The proof cost is `A`, not the number of
events. Metric scale is retained by `f,g`, and exact inequalities convert an
active-period bound into wall-count and extent bounds.

This also clarifies the role of Tournament Analysis. Chronological comparison
gives a transitive tournament with one Hamiltonian path, zero cycles, and
singleton SCCs. Those fingerprints are correct but nearly content-free. The
proof carrier is the decorated path: owner labels, absorbed block lengths,
absolute slow scale, and zero-sum visitor hyperedges. Changing vertices from
runners to wall events was not enough; maximal free blocks must be contracted.

The remaining prime-seven theorem is therefore:

1. bound the number or metric density of active zero-sum packets;
2. do so with varying Beatty indices and possible order flips—MISTAKE-148
   rules out the fixed-index de-phase shortcut;
3. intersect the resulting blocked interval with the actual core-safe
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

1. **Shallow n=12:** construct a decreasing height/owner invariant extending
   THM-770 beyond lift height twelve; another raw box is not recursion.
2. **Deep two-sheet:** prove global erosion noncontainment by selecting a deep
   component, or force the dyadic seam tree into the existing finite bases.
3. **Prime-seven r=8:** bound active visitor packets with varying-index
   arithmetic and then use core incidence; do not count raw walls.
4. **Higher sheets / scale-normal families:** retain colour ownership and
   ramification through descent, rather than collapsing to residue counts.

The metagraph transitivity-flow work remains valuable as an organizational
atlas, especially for locating which finite fibres are balanced or rigid. Its
proper role is diagnostic until a reconstruction theorem attaches the metric
LRC stalk. The common strategic rule is now precise:

> quotient free refinement first; then recurse on globally competing,
> owner-labelled metric components.

