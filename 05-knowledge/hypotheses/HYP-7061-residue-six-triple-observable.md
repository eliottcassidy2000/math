# HYP-7061 — the residue-six triple observable

**Status:** STANDALONE OPEN / BYPASSED FOR THE LIMITING SIGN (codex-2026-07-16-S18).
`THM-904`'s explicit pointwise triple reduction is finite-exact, but its universal
three-speed bound remains open.  `THM-909` now closes the limiting negative residue-six
sign directly at `32/343`, so this inequality is no longer required for that closure.
This file was renumbered from `HYP-7052` after pulling boxeph's earlier claim of that ID.

For the rational weight `beta` in `THM-904`, set

```text
q(a,b,c)=integral beta(sec(a x),sec(b x),sec(c x)) dx.
```

The target is

```text
q(a,b,c) <= 47/100
```

for every three distinct positive integers.  This would independently give
`-F_6<=47/490<0.097`; `THM-909` proves the stronger required limiting conclusion by a
different signed residue quotient.

Exact evidence currently covers all `28,876` primitive triples with largest speed at
most 60.  The unique maximum is `q(1,4,7)=81/175`.  The observed support stabilizes on
compact additive-relation shapes rather than scale-free generic triples.

The correct proof quotient is not runners, arcs, or the residue tournament.  It is the
84-state unordered triple-sector observable together with its additive-relation sidecar.
This quotient preserves the integrated negative `K_6` certificate but destroys the
identity of the fourth/fifth movers and finite-`t` wall chronology.  The alternative
four-runner box carrier in `THM-898-fourrunner-boxhit-relation-stratified` retains more
geometry but encounters the same nondecaying additive-relation plateaus.

The proposed `x -> 1-x` / negation transfer is refuted by `THM-903`: time reflection
acts within a fixed residue class.  In the pinned stationary-point frame it sends
`s -> 6-s`, so `{1,5},{2,4}` are exactly the reflection-fixed inner pairs.  The hard
support is therefore reflection-even, explaining structurally why the triple observable
must retain a cubic channel.

Next proof step: perform the exact ANOVA/Fourier split of `beta`; bound its one- and
two-coordinate projections using the 21 pair rays; then enumerate short primitive
three-term relations and bound the absolutely summable zero-marginal cubic tail.
