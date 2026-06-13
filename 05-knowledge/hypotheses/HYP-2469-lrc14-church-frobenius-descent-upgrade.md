# HYP-2469 - Church-style partial-Frobenius descent is the right LRC14 proof skeleton

**Status:** OPEN synthesis with exact finite anchors.

**Source:** codex-2026-06-13.  Upgrades HYP-2445 using HYP-2463, HYP-2464, and HYP-2465.

**External source:** Benjamin Church, arXiv:2508.14876,
`https://arxiv.org/abs/2508.14876`.

**Computation:** `04-computation/lrc14_church_frobenius_descent_codex.py`; stored output `05-knowledge/results/lrc14_church_frobenius_descent_codex.out`.

## Claim

The useful transfer from Church's product-quotient paper is not the `1092/91/78` numerology by itself.  The useful transfer is the proof grammar:

```text
lossy scalar passes
+ retained side channel on every twist/fiber
+ finite exceptions or strict descent
=> no global bad object
```

For Church's surfaces, the scalar shadow is Shioda supersingularity.  The retained channel is diagonal symmetric forms on every asymmetric partial Frobenius twist.  The descent sends non-exceptional curves through an Ekedahl/partial-Frobenius quotient and drops a projection degree.

For LRC14, the scalar shadow is raw denominator blocking, especially "all plain q<=27 shells are blocked."  The retained channel is the Q27 obligation/support ledger:

```text
safe twist obligations,
13-clock debt,
deleted-core address,
shell-27 antipodal class,
divisor fiber,
owner/Bprime target,
support-load data.
```

The proposed theorem is:

> Any primitive LRC14 row with no Q27 witness either enters one of the certified finite atlases, opens a named exception, or descends by a strict resource drop.

## Exact Anchors Now Available

The new script composes the current certified blocks:

```text
HYP-2444/HYP-2443:
  one-stranger family: 936/936 Q27 covered and 936/936 Bprime(any) covered.

HYP-2463:
  hard packet replacement hull: 77,520/77,520 Q27 covered.

HYP-2464:
  one-delete/two-add plain residuals: 877/877 Q27 covered.
  Every residual has a 13-clock speed and a deleted-core address;
  no residual deletes 7,21,49.

HYP-2465:
  near-core primitive set-cover: for |D|<=3, all 299 deletion cases
  are infeasible with add budget |D|+1.

HYP-2470:
  eight-core boundary: for |D|=4, Q27-only has exactly two feasible
  deletion addresses among 495, but both open at plain q=31/33 with
  Bprime and positive measure; Q27 union shell<=41 is infeasible for
  both addresses.
```

Thus, inside the carry window `1..1092`, any no-Q27 row must delete at least four of the twelve speeds in `CORE=7*{1,...,12}` unless it avoids the near-core normalization entirely.  More sharply, any row with no Q27 witness and no plain witness `q<=41` must delete at least five core speeds.

## Descent Skeleton

Assume `S` is a primitive 13-speed row with no Q27 witness.

The Church-style target is a decision tree:

```text
1. If S normalizes into 1..1092 and retains >=9 core speeds,
   HYP-2465 contradicts no-Q27.

2. If S compresses to shell-27/13-clock hard packets,
   HYP-2463 contradicts no-Q27.

3. If S is a one-delete/two-add plain residual in the carry window,
   HYP-2464 gives Q27 and the resource coordinates.

4. If support concentrates on a runner,
   owner-private deletion or Bprime(any runner) should open.

5. If S is AP, Vstar, or nonprimitive 2AP,
   it is a named wall atom, not a new blocker.

6. If S deletes exactly 4 core speeds,
   HYP-2470 gives Q27 or plain q<=41.

7. If S deletes >=5 core speeds,
   prove a below-eight-core support-load descent.

8. If S leaves the carry window,
   prove a large-speed Bprime/divisor/carry normalization.
```

The remaining live theorem is therefore not "find a bigger q."  It is:

```text
no-Q27 and no-plain-q<=41 blocker
=> below-eight-core or outside-window
=> strict descent / named exception.
```

## Proposed Resource Rank

The script proposes the live descent measure:

```text
R(S) =
  (outside_window?,
   deleted_core_count,
   Q27_setcover_defect,
   max_support_load,
   owner_private_width,
   low_clock_depth)
```

Generic moves should lower `R` lexicographically.  If a move does not lower `R`, it should land in a finite exception list: AP, Vstar, nonprimitive 2AP, `q=91`, `q=161`, owner-private/Bprime, or a low-clock witness.

The first computational target is the `|D|=4` analogue of HYP-2465.  A raw `e+1=5` set-cover MILP may be too large or too permissive; the better target is a typed budget that separates 13-clock, divisor-fiber, shell-27, owner, and low-clock resources.

## Paper Bridge, Carefully

The exact beacons remain worth keeping:

```text
|PSL2(F_13)| = 1092 = 84*(14-1) = 13*84.
D6/A4 index 91 = C(14,2), matching the LRC q=91 fiber.
D7 index 78 = C(13,2), matching the code72 lambda_5 beacon.
```

But HYP-2469 downgrades these to search beacons.  They do not prove anything unless the retained side channel is present.  This is exactly Church's lesson: scalar supersingularity can pass while unirationality fails.

For LRC14, the retained side channel is now concrete: Q27 set-cover obligations and marked support/owner data.

## Tournament Analysis

The computation uses proof routes as vertices.  The observable is:

```text
(retained channel, exactness, LRC leverage + descent drop, computability, -risk)
```

Stored outcome:

```text
score histogram: {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed 3-cycles: 0
edge flips versus scalar-only ranking: 4/28
leader: lrc14_Q27_obligation_setcover
second: church_diagonal_forms_all_twists
```

The ranking is a useful correction to HYP-2445: Church supplies the mature proof template, but HYP-2465 supplies the current strongest LRC-specific object.

Assumption challenge: candidate vertices considered were runners, individual denominators, unit twists, deleted core speeds, added speeds, subgroup indices, Frobenius twists, curve components, wall events, and proof obligations.  The selected quotient preserves descent/exception proof structure.  It destroys raw time geometry, explicit surface geometry, and arbitrary multi-stranger interactions.

## Next Moves

1. Build a below-eight-core typed scout for `|D|>=5`, with typed budgets rather than raw cardinality.
2. Prove or test an outside-window normalizer: a speed `v>1092` must open Bprime(any), dominate an existing core speed, or reduce into a divisor/carry fiber without losing blockedness.
3. Build support-transport curvature across `q -> 2q`, `q -> 7q`, and `27 -> 9 -> 3`.  Nonzero defect should be a witness; zero defect should force descent.
4. Make the exception catalogue explicit: AP, Vstar, nonprimitive 2AP, `q=91`, `q=161`, `q=31/33`, owner-private/Bprime, and low-clock exits.
