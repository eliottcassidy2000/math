# Invert the whole fixed-h state and rule out closing the original fold at h=15

The earlier [fold obstruction](../chanpak_tail_return_20261115/FOLD_CLOSURE.md) fixes an ordinal boundary. Here the inverse allows **every** spine with a specified ordinary Hamilton count h, regardless of order, contact mask, or nonordinal edits used to make it. This separates a failure of one seam operation from a failure of the proposed target interface.

## Backward algorithm for any fixed ordinary count

An ordinary tournament spine decomposes uniquely as an ordinal sum of strong components. Its Hamilton count is the product of their counts. If a non-singleton factor has count a, [Busch's bound](https://doi.org/10.37236/1141) gives order at most 1+3 log_5 a, so finitely many strong last-component types need be checked for each divisor a of the prescribed h. Enumerate each apex mask on each literal component. Retain the finite local bank (a,Q,n,s,e,t,f) of ordinary count, two-slot count, marked count, marked start total and bias, and unmarked end total and bias. Include the two painted singleton types; these are the only unbounded size resource via repetitions.

Given a final state with ordinary h and weighted start S, try each last local profile with a|h and put u=h/a. The [exact cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md) forces the preceding start and the final bias:
```
S_old=(S−u s)/Q,     d_final=u e.
```
Reject nonintegral or negative starts, and reject any incorrect demanded bias. Recurse over all possible old (N,d) at (u,S_old), then use the actual physical update
```
N_final=a N_old+(n−a)u+t S_old+f d_old.
```
When S>0 each preceding S_old≤S/2, because Q≥2. When S=0 the entire marking is empty, as every vertex has a positive singleton start weight. There is only (N,d)=(h,0), attained if the fixed h bank can produce h by multiplication (with the empty spine allowed for h=1). In a witness search one stores any such representative. This terminates and is exact at all orders: every branch is an ordinal assembly of actual SCCs; conversely the last SCC of any physical witness shows up in the bank. Allowing the whole representative to be nonstrong enlarges the inverse, which makes its negative conclusions stronger.

For h=15, nontrivial factors are 3,5,15, all of order at most six under that bound. The attached program enumerates **every labeled** strong tournament through order six, uses endpoint-subset path DP for the local data, and separately counts one literal cap for each of the 168 distinct nonsingleton signatures. The bank sizes are 2,4,16,148 at a=1,3,5,15.

## The first eight-count fold, at either state or any state between them

The [first eight-count overlay](../chanpak_tail_return_20261115/TAIL_FOLDS.md) is all odds in 150401–150415. Its two original right starts are 272018 and 273938 with d=0 and h=15. At those exact original starts, the inverse finds:
```
S=272018:  150409,150411,150413,150415
S=273938:  150401,150403,150405,150407,150411
```
These lists are intersections with the requested eight values, not the full N fibers. A full inverse sweep over each integer S between them finds no shared run longer than four. The sweep also checks the neighboring width-five and width-six phases and their next dyadic occurrences; their receipts retain every scoped maximum. The inverse is independent of which arcs a recontact is allowed to edit.

## Stronger: no common h=15 continuation at any start, even after neutral preparation

A common continuation state would give a common oddpart(S+h). Indeed appending a painted 1 doubles S+h and preserves N. Thus it suffices to show no common normalized start among even a subset of the targets. For this use the five targets coprime to 15:
```
150401,150403,150407,150409,150413.
```
There is a finite, complete forward search for their normalized starts, not an imposed vertex cutoff.

Let X be the largest of these targets. At any partial SCC sequence, prune if (15/h)N>X: every future SCC count multiplication a increases N by at least the factor a. This follows from n≥a and tS+fd≥0 (both parity totals are nonnegative). In the finite signature bank let B=max over factors of (15/a)|e|, so |d|≤B; here B=15. Once S>X+B, a singleton painted 0 would add at least S−B>X. Any nonsingleton SCC with at least one unmarked end vertex adds at least 2S: an unmarked vertex has both its odd singleton path and an even two-vertex path into it, by strongness, so t−|f|≥2. Starts cannot decrease under extension. Thus from beyond the threshold the only possible future blocks are fully marked. A fully marked nonsingleton multiplies N by its factor a>1 dividing 15 and cannot contribute to any of the five unit targets. Fully marked singletons do not change the normalized start or N, so they can be stripped after recording the current state. At S=0 all histories for a given h coincide. These are exhaustive stop rules.

The complete signature search gives normalized-start counts 28,44,27,51,36. The running intersection has sizes 28,3,1,0,0. Therefore no h=15 spine bank supports the original whole eight-count fold under a common interface, **at any order or S**, even allowing independent count-neutral transitive preparations. The argument does not exclude a different h, or another range of eight targets. In particular a compensated return that changes h is not a cosmetic change: retaining h=15 is impossible for closing these specific counts.

## Random reading and what it changed

A seeded sample of September notes surfaced the [source-linear all-degree boundary obstruction](../continuing10_20260908_dg_linear_all_m.md), its [quadratic escape on the same surface](../continuing11_20260908_quadratic_unit_two.md), and the [adjacent-factor precision audit](../overnight14_20260906_jets_p31_adjacent_even_audit.md). The usable idea was to distinguish the *entire physical fiber at fixed response* from any particular visible presentation, and to test a stronger source change after the lower-order obstruction. No theorem about Jacobian maps or Smith ideals is used here.

Replay: compile `replay_inverse_fibres.cpp` with C++17 and compare with its text receipt. The finite stop rules certify the stated all-order negative result; the sampled further heights are not an all-height nonclosure theorem.
