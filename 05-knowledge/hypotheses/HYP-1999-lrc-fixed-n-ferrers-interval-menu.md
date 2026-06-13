---
id: HYP-1999
status: OPEN
source: codex-2026-06-01-S525
related:
  - HYP-1981
  - HYP-1982
  - HYP-1987
  - HYP-1995
  - HYP-1996
  - HYP-1997
  - HYP-1998
  - THM-381
  - THM-382
  - THM-383
  - THM-384
---

# HYP-1999: Fixed-n LRC source targets are Ferrers interval tournament menus

## Statement

For total LRC denominator `n`, the open source target in the observer-marked
tournament model is a finite menu of tournament isomorphism classes on
`m=n-1` vertices.

This refines HYP-1996's circular-tournament menu and HYP-1998's round-necklace
body by imposing the fixed LRC safe-arc length `1-2/n` and by separating open
classes from compactified wall classes.  It is complementary to HYP-1997: the
metagraph walk describes the arithmetic motion, while the Ferrers interval
menu describes the fixed-`n` target set that the walk must hit.

After deleting the source observer, the `m` moving runners lie in the safe arc

```text
[1/n, 1-1/n],  length L = 1 - 2/n.
```

Ordering those moving runners in the arc, the runner sub-tournament is

```text
i -> j  iff  x_j - x_i < 1/2   for i < j.
```

Equivalently, the backward edges are exactly the pairs whose arc separation is
larger than `1/2`.  Monotonicity of separation implies that the backward-edge
set is a Ferrers filter in the triangular pair poset.  Thus LRC at fixed `n`
can be modeled as:

```text
Every primitive arithmetic clock must hit the Ferrers interval source menu
on m=n-1 tournament vertices, or its compactified wall boundary, inside the
A000568(m) isomorphism-class universe.
```

## Evidence

`lrc_fixed_n_tournament_class_atlas_s525.py` enumerates the open Ferrers
interval menus for `n=4..9` and compares them to open source classes reached
by bounded speed clocks for `n=4..8`.

Open menu counts:

```text
n=4: m=3, A000568=2,    open_iso=1
n=5: m=4, A000568=4,    open_iso=2
n=6: m=5, A000568=12,   open_iso=4
n=7: m=6, A000568=56,   open_iso=6
n=8: m=7, A000568=456,  open_iso=10
n=9: m=8, A000568=6880, open_iso=16
```

The interval-feasible Ferrers signatures are tiny compared with the Catalan
many monotone filters:

```text
n=5..9 feasible signatures = 8,16,32,64,128
```

Bounded speed-clock cross-checks match the geometric open menu exactly:

```text
n=4: 1/1 open classes match
n=5: 2/2
n=6: 4/4
n=7: 6/6
n=8: 10/10
```

The S512/S520 closed source menus are larger in some cases
(`2,2,6,6,>=12` for `n=4..8`) because equality-wall witnesses add compactified
boundary/tie-path classes.  Hence the proof object is not raw A000568 alone;
it is an open Ferrers tournament menu plus THM-383 boundary data.

Tournament Analysis over the fixed-`n` profiles is transitive:

```text
score_hist=((0,1),(1,1),(2,1),(3,1),(4,1),(5,1))
c3=0
SCCs=(1,1,1,1,1,1)
Hamiltonian_paths=1
```

## Predictions

1. The feasible-signature count `2^(n-2)` for `n>=5` should have a direct
   proof from the fixed length `L=(n-2)/n`.
2. The open isomorphism-class count sequence
   `1,2,4,6,10,16,...` agrees with HYP-1996's circular `2*Fib` menu from
   `n>=5`, while `n=4` exposes the safe-arc boundary cutoff.
3. The extra closed classes from S512/S520 should be exactly describable by a
   compactified tie-path boundary of the Ferrers interval menu.
4. A proof of LRC at a particular `n` should either force an open hit in this
   Ferrers menu or force one of the compactified wall classes from THM-383.
5. For hard frontiers such as `n=14` and `n=18`, the menu should remain small
   relative to `A000568(n-1)`, but its boundary fibers should carry the tight
   initial-segment and gate-ladder witnesses.

## Sources

- `04-computation/lrc_fixed_n_tournament_class_atlas_s525.py`
- `05-knowledge/results/lrc_fixed_n_tournament_class_atlas_s525.out`
- `04-computation/lrc_source_reachability_deep_s512.py`
- `04-computation/lrc_source_reachability_n8_s520.py`
- `05-knowledge/results/lrc_source_menu_saturation_s520.out`
- HYP-1981
- HYP-1982
- HYP-1987
- HYP-1996
- HYP-1997
- HYP-1998
- THM-381
- THM-383
