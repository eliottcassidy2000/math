# LRC(14): the missing object is one marked tooth word

**Session:** codex-2026-07-19-S78 coherent-word continuation

**Status:** frontier synthesis after THM-1252--1256.  This memo records a
strictly sharper common carrier and the precise remaining implication.  It
does **not** claim six-comb noncoverage, LRC(14), or uniform emptiness of the
`n=12` sporadic branch.

## 0. Outcome

Several routes that formerly looked only analogous now become literal
projections of one object.  Start with a hypothetical cover of a complete
`c`-safe gap by six faster danger combs.  Choose one deletion-minimal cover by
individual strict teeth, order it chronologically, and choose every centered
spoke blocker from this same cover.  The state is

```text
(G; ordered teeth I_a=(alpha_a,beta_a,s_a,n_a);
 every raw handoff W_a;
 six centered phases t_i and their marked teeth T_i;
 blocker-cycle incidence; relative digits and reflection).             (1)
```

Four losses disappear simultaneously.

1. Deletion minimality gives `beta_a<=alpha_(a+2)`, so **all** raw handoffs
   are pairwise disjoint.  The full word, not one spanning tree or a
   one-third Cayley average, is paid:

   ```text
   H_fast>=1/c+(7/12)sum_a 1/lcm(s_a,s_(a+1)),
   F_6>= (c/16)sum_a 1/lcm(s_a,s_(a+1)).              (2)
   ```

2. Every marked blocker tooth lies wholly inside the carrier gap and has a
   selected provider at each wall.  Those providers are its immediate word
   neighbours.  Thus every blocker edge carries two opposite, disjoint,
   lcm-quantized wall seams.  A common provider label gives the exact detuned
   toothpick rung

   ```text
   omega_-+omega_+
    =[h-(7r-1)j]/(7jh),
   1<=r<=335,
   h>(7r-1)j.                                         (3)
   ```

3. A blocker cycle has a speed-descent edge.  Its THM-1248 digit is binary.
   According to the order of its two marked teeth, the original or reflected
   mixed circuit has nonnegative affine factor.  The chronological/centered
   drift invoice is therefore unconditional; the former `c<=1171` branch is
   only a guardrail for deliberately arbitrary digit choices.

4. The invoice has now been interpreted honestly.  If the chronological
   endpoint teeth are `(s_0,n_0),(s_r,n_r)`, then

   ```text
   n_0 n_r Delta=n_r s_0-n_0 s_r,
   R-(n_r s_0-n_0 s_r)=s_0(P-n_r)=s_0(k+delta).       (4)
   ```

   Thus it is an exact address/position invoice, not a new independent
   overlap mass.  The binary phase order either agrees with marked-tooth
   order, in which case the intervening subword literally covers the segment
   between the centered phases, or disagrees, in which case the two marked
   teeth must be consecutive and overlap.  This is the first precise version
   of “transport the germ or charge its first failure.”

These results close the free-choice and sign obstructions.  The remaining
problem is no longer to find a blocker cycle, a second wall, a located tree,
or a bounded relative word.  It is to consume the resulting **typed closed
cell**: two wall seams, one binary phase landing, and the exact address
surplus in (4).

## 1. Why the common word was missing

The earlier proofs selected their useful structures independently:

- the centered blocker map chose one dangerous label at each sampled phase;
- the handoff theorem chose a minimal cover or a nerve tree;
- Fano and `chi_7` chose quotient incidences;
- the functional `H`-drift chose a density and a ratio chamber;
- the toothpick ladder chose one boundary crack and its next blocker.

Each selection was individually valid, but there was no reason for its edge,
tooth, or address to be the one used by another selection.  “A cycle plus a
tree” did not imply a cycle edge lived on the tree.  “A wall seam plus positive
holonomy” did not imply their orientations were compatible.

Choosing all centered blockers from one already fixed irredundant subcover
repairs this without strengthening any local hypothesis.  Every sampled
blocker is a literal marked tooth in a single linear order.  The functional
cycle and the continuous overlap chain are now two typed relations on the
same vertices, with addresses retained.

## 2. The exact interval laws

Let the selected teeth have increasing left and right endpoints.  If
`I_a` and `I_(a+2)` overlapped, their union would contain `I_(a+1)`, making
the middle tooth redundant.  Hence

```text
beta_a<=alpha_(a+2).                                  (5)
```

The consecutive overlaps

```text
W_a=(alpha_(a+1),beta_a)                              (6)
```

are positive, lie in `int(G)`, and are pairwise disjoint.  Integrating the
literal coverage multiplicity, rather than applying a graph inequality,
gives (2).  This is why the coefficient improves from `7/36` to `7/12`.

There is a second word law.  An immediate backtrack

```text
h,j,h                                                   (7)
```

places two distinct `h`-teeth at the walls of the middle `j`-tooth and
forces `h>6j`.  Therefore

```text
h,j,h,j                                                (8)
```

is impossible: its two middle positions would require both `h>6j` and
`j>6h`.  Every run of one unordered edge symbol has length at most two.  If
the word has `N` teeth, its number of nonbacktracking turns is at least

```text
ceil((N-1)/2)-1,                                      (9)
```

as well as the independent six-label floor `4`.  Toothpick
self-similarity is therefore real for one rung but cannot repeat on the same
two-label edge.  Every attempted second rung must branch to a third owner.

## 3. How the historical viewpoints land on (1)

### Functional `H`-drift

The twelve-piece envelope is not merely a scalar side condition.  Its density
lower bound charges the same disjoint handoffs as the harmonic discrepancy,
giving the second inequality in (2).  The useful functional coordinate is
therefore the lcm-decorated owner word, not the six ordered ratios alone.

### Toothpick ladders and Kakeya needles

Equation (3) is the literal finite-rung functional form.  The skeleton
`7r-1` is the self-similar shape; the positive detuning
`h-(7r-1)j` is the metric residue that prevents it from becoming a purely
combinatorial recursion.  Equation (8) proves that a two-label ladder cannot
iterate.  The Kakeya “direction” is `j/h`; its indispensable offset is the
address pair `(N,M_-,M_+)`.

### Fano and `chi_7`

Fano lines still organize seven obligations, and `chi_7` still records the
sign of zero seams.  But (3) shows what their quotient discards: the same
negative mod-seven skeleton can carry any positive detuning and any rung
`1..335`.  A colour becomes useful only after attachment to an actual wall
vertex of (1).  The still-open Fano/`chi_7` probe should therefore be read as
a possible organizer of several already placed forks, not a source of
placement.

THM-1260 now makes that warning exact.  Even after placing both walls of one
sharp fork, its binary phase side, and the next marked blocker, every one of
the `2^7` speed-colour words survives at every rung `1..335`.  The positive
seam digit satisfies a shifted congruence with no quadratic-character
content.  The first Fano-shaped object not refuted by this local surjectivity
is a triangle of *handoff occurrences* with their seam digits and separated
owner-reuse locations; runner point colours alone are exhausted.

### The `j=4` flood tail

The closed continuum four-comb tail supplies the positive global reservoir,
while the old flood atlas taught the operation-congruence lesson: a quotient
must survive rerooting and the next insertion.  Marked tooth position is the
slow-gap analogue of the missing rerooting coordinate.  The word updates by
local insertion; a bare Fano flag or edge current does not.

### Blocker tournaments

The speed tournament is transitive and loses everything decisive.  The
blocker relation has a directed cycle but is not a tournament.  The marked
position relation is transitive.  Overlaying the two typed relations forces a
backward edge; binary speed descent and reflection then select the valid
mixed circuit.  Collapsing them into one untyped orientation erases the proof.

## 4. The honest frontier after the sign collapse

The strongest closed local cell now consists of

```text
binary descent i->j;
marked target tooth J=(j,N) and its two adjacent wall handoffs;
original/reflected chronological subword;
phase-order alignment bit;
address surplus P-n_r=k+delta;
exact lcm/detuning sheets.                            (10)
```

If phase and word orders mismatch, (10) already creates a direct additional
handoff occurrence in (2).  If they align, the subword covers the actual
segment between the two centered phases.  What remains is one of the
following genuinely stronger statements.

1. **Alternate-gap descent.**  Show that the aligned segment reaches a new
   complete safe gap with a smaller positive carrier/address potential while
   retaining one wall fork.
2. **Metric conversion.**  Turn the address surplus in (4) into raw overlap
   or weighted functional mass beyond what (2) has already counted.  Merely
   renaming (4) as “drift” is not enough.
3. **Finite typed-word exclusion.**  Use the compact ratio box, rung alphabet,
   ABAB exclusion, binary phase landing, and exact kernel sheets to exclude
   every cyclic typed word without first bounding absolute tooth addresses.

The main guardrail is equally clear: `R>=n_0n_r Delta` simplifies to
`P>=n_r`.  It cannot be counted again as an independent Hunter credit.  Any
closing argument must use where the address lands, not just its sign.

THM-1262 now settles the first recursive base of that operation.  In a
blocker two-cycle the ascent target is protected from the reverse owner's
danger tooth.  The two marks are nonconsecutive, so the binary descent cannot
be a one-step inversion; it is aligned and its actual phase corridor exports
a protected third-owner seam.  The shortest cycle therefore does not form a
closed two-label stalk.  What remains is to iterate this forced third-owner
bridge with a well-founded address potential, or close several such bridges
into THM-1260's still-unresolved seam-digit incidence circuit.

## 5. Finite terminal bases and the separate `n=12` branch

The bounded-variation density method now gives exact terminal noncoverage
through carrier `44` (THM-1255, THM-1257, THM-1259, and THM-1261).  The
carrier-41 densities do not transport unchanged to carrier 42, and carrier
44 first needs cutoff depths 60 and 76.  The self-similar operation is density
*selection* with controlled variation on the changing cyclic phase orbit,
not reuse of a fixed stalk or finite atlas.
These computations are useful terminal bases for an eventual alternate-gap
descent, not substitutes for it.

The `n=12` sporadic branch remains separate.  THM-1249 and THM-1258 close the
AP-centred common-scale-35 and common-scale-36 Hamming-six faces; scale 37 is
prime-excluded and scale 38 is next.  The independent single-far absorption
atlas also proves the whole `N=12` first-gap single-far stratum empty.  The
remaining smooth H5/H6 banks, multi-defect/Freiman packets, and deep
reverse-content component cover remain.  No
coherent-word theorem above supplies AP extraction or arbitrary-height
ballot rigidity, so uniform sporadic emptiness must not be claimed here.

## 6. Target adjustment

The next proof attempt should not ask for another global mass inequality,
another sampled blocker, a stronger Fano count, or a literal repetition of a
carrier-41 density.  It should take the cell (10) as input and implement the
next legal operation on it.  The most falsifiable version is:

> An aligned binary phase segment with both target-wall seams either reaches
> a strictly smaller covered carrier/address cell, or contains a first
> nonaligned marked tooth whose raw overlap is a new occurrence in the full
> chronological invoice.

THM-1256 proves this statement for the first interval-order comparison.  The
remaining work is to make the iteration preserve a complete slow gap and a
well-founded potential.  That is now the narrowest known mathematical
frontier of the six-comb route to LRC(14).
