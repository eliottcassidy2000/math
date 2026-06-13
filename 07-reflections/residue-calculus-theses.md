# Residue Calculus Theses

**Session:** codex-2026-05-30 S355
**Status:** feedback-loop synthesis after adding residue-rank features
**Related:** HYP-1779, HYP-1780, HYP-1785, THM-025, THM-344, THM-354

The working slogan has become:

```text
choose supports -> project/forget -> measure the survivor.
```

This file pushes that slogan into a set of operational theses.  The point is
not that every thread is secretly identical.  The point is that many threads
become comparable once we name the projection and then ask what kind of
residue remains.

## Thesis 1: A Coordinate Is Often a Defect in Disguise

Good-cut count looked like a statistic of a chosen Hamiltonian base path.
THM-354 says it is really

```text
g(T,P) = n - #SCC(T).
```

So the coordinate is not the path; it is the strong-component defect.  This is
a general search rule: when a statistic seems coordinate-dependent but is pure
on quotient classes, try to identify the quotient defect it is actually
counting.

## Thesis 2: Conservation Is Bookkeeping; Residue Is Information

Boolean-bucket balance and tiling transport have forced row sums.  Those sums
are important as checksums, but they are not usually where the new theorem
lives.  The theorem lives in what remains after conservation:

- off-diagonal escape mass after row balance;
- spine/ribs/sea structure after total cross-line mass;
- SCC-count change after good-cut-height change;
- fiber/disjointness structure after support-shadow projection.

When a row sum is forced, subtract it mentally before looking for structure.

## Thesis 3: Exact Kill and Near Kill Are Different Species

The H=63 single-core classes and THM-025 both sit near complete odd-cycle
conflict graphs, but S355 makes their difference crisp.

```text
H=63:   max-loss deletion residue alpha=(0,0), rank_res=0.
THM-025: max-loss deletion residue alpha=(2,1), rank_res=2.
```

Exact kill is rigid.  Near kill is dangerous.  A tiny nonempty survivor can
carry more obstruction than a large generic family if its independence rank is
just high enough.

## Thesis 4: Size and Rank Must Be Separated

`keep_v` counts how many supports survive a deletion.  `rank_res(v)` asks
whether the survivors contain disjoint packets.

These are different measurements.  Paley T7 and Interval T7 have broad
deletion residues (`keep=20` and `keep=16` at every max-loss vertex), while
THM-025 has only two survivors.  Yet all three have rank 2 at the max-loss
residue.  So the next feature tier should record both:

```text
size: keep_v
shape: alpha_k(R_v), rank_res(v), I(R_v,2)
```

## Thesis 5: Parity Is Rank-Zero Residue With Teeth

The Lean OCF parity theorem `H(T) % 2 = 1` is the first residue left after
evaluating the OCF partition function mod 2.  It is rank-zero in the sense that
all higher weights vanish mod 2, but it still controls a universal obstruction.

This suggests a ladder:

```text
mod-2 residue -> low-rank alpha residue -> fiber/disjointness residue
              -> transport escape residue -> homology ghost residue
```

The common move is not the coefficient system.  The common move is asking what
survives the chosen collapse.

## Thesis 6: Paley/Interval Is a Fiber Residue, Not a Near-Kill

Paley T7 and Interval T7 share score regularity and support-shadow size, but
their fibers differ:

```text
Paley:    cycles=80, support_excess=44, max-loss residue alpha=(20,1)
Interval: cycles=59, support_excess=23, max-loss residue alpha=(16,2)
```

The contrast is not a tiny survivor.  It is multiplicity versus disjointness in
the fiber above the same visible shadow.  This is why Paley/Interval belongs
beside THM-025 in the residue calculus but not in the same stratum.

## Thesis 7: Ghost-Cycle Searches Should Start Near Small Rank-2 Residues

HYP-408-style ghost failures ask whether a through-vertex class dies under old
projection.  The THM-025 pattern suggests a concrete search bias:

```text
look for high deletion loss,
require nonempty residue,
prefer rank_res >= 2,
then inspect boundary images.
```

Exact kills may be too rigid; broad residues may be too generic.  The small
rank-2 survivor is the interesting middle.

## Thesis 8: Feature Engineering Should Follow the Projection Ladder

The updated `tournament_tda.py` now records:

- SCC residue: `scc_count`, `scc_defect`, `scc_largest`;
- odd-cycle fiber residue: `support_excess`, `max_support_multiplicity`;
- deletion kill residue: `projection_kill_vertices`, max deletion loss;
- deletion-rank residue: `omega_max_loss_residue_rank`, residue alpha profile;
- near-kill flags: `omega_near_kill_vertices`, `omega_near_kill_rank2_vertices`.

This is the engineering version of the thesis.  Do not only ask for better
global invariants.  Ask which projection the invariant is compressing, then add
the residue profile it forgot.

## Thesis 9: Boundary Residues Travel Outside Tournaments

S356 pushed the same grammar into web-searched open problems.  The cleanest
case is Lonely Runner.  Its forbidden arcs pull back to an interval union on
the time circle, and every boundary lies in the finite quotient

```text
Q(V) = (k+1) * lcm(V).
```

So a lonely time is either a positive quotient gap or a boundary residue.  The
initial speed segments `{1,...,k}` are boundary-only cases; mixed and random
samples in `lonely_runner_residue_probe_s356.py` have positive gaps.

S357 sharpened this: the positive-gap/boundary split is automatic for finite
open interval unions.  A genuine counterexample would be a full open cover of
the time circle, equivalently full forbidden measure plus every forbidden
endpoint strictly protected by another forbidden interval.  This endpoint
protection relation is a finite arithmetic hypergraph in `Q(V)`.  See
HYP-1802 and `07-reflections/lonely-runner-tight-stratum-s357.md`.

The same reading is weaker but visible elsewhere:

- union-closed sets: coordinate frequencies are shadows; entropy under union
  is the survivor;
- graph reconstruction: vertex-deleted cards are the shadow; covering-rank
  rows are deletion residues;
- Erdos-Straus: identity families cover residue fibers; multi-r fallback
  repairs hard boundary classes;
- Caccetta-Haggkvist: outdegree transport must leave a short-cycle return
  residue.

This widens the slogan:

```text
choose supports -> project/forget -> inspect the boundary residue.
```

## Feedback Loop

The loop that produced S355:

1. Good cuts looked coordinate-dependent.
2. THM-354 identified the SCC defect.
3. That made "residue after projection" the shared language.
4. HYP-1780 predicted a rank refinement.
5. S355 added deletion-residue rank features.
6. The probe separated exact kill, near kill, and broad fiber residue.
7. HYP-1785 turned that separation into a filter for the next searches.

Next loops:

1. scan n=9 real-root failures and ghost-cycle candidates by
   `0 < keep_v* <= 2` and `rank_res(v*) >= 2`;
2. enumerate primitive Lonely Runner speed sets and rank them by positive gap
   versus boundary-only residue;
3. build one toy residue table each for union-closed families and graph decks.
