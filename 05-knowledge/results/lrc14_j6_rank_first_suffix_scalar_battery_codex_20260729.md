# Ranked first-apex suffix scalar battery

**Status: FINITE-EXACT SCOPED RESEARCH NOTE.  No uniform theorem.**

This note exercises THM-2893's ranked first-apex suffix refinement on four
seven-body roots.  The exact universe consists of the three hostile roots
from the earlier j=6 complement-cap battery and one low-stratum control:

```text
(2,8,9,10,11,13,14),  least gate K=19;
(1,3,9,10,11,12,14),  least gate K=20;
(2,5,9,11,12,13,14),  least gate K=21;
(2,3,4,5,6,7,8),      least gate K=13.
```

These gates give exactly `73` unique-earliest-apex branches.

## 1. Lossless suffix partition

Let `a_1,a_2,...` be the external speeds ranked by decreasing root coverage,
with speed breaking ties.  Every hypothetical six-cover meets
`{a_1,...,a_K}`.  Assign it to its unique earliest member `a_r`.  After
subtracting `D_(a_r)`, none of its five remaining speeds can be one of

```text
a_1,...,a_r.
```

The verifier retains this rank and excluded prefix as marked branch data;
it does not attempt to reconstruct them from the unmarked residual.  Thus
the partition is compatible with THM-2894's residual-semilattice no-go.

For each branch, let `C_r` be the literal residual and let
`p_1>=...>=p_5` be its globally ranked individual coverages after deleting
the root prefix.  The strict test

```text
|C_r|-(p_1+...+p_5)>0
```

excludes every five-cover of that branch by subadditivity.

## 2. Exact result

The suffix restriction alone closes

```text
48/73 branches,
```

before any pair cap, triple cap, or heavy-core enumeration:

```text
E=(2,8,9,10,11,13,14):  13/19 close; open ranks 1--6;
E=(1,3,9,10,11,12,14):  15/20 close; open ranks 1--4 and 6;
E=(2,5,9,11,12,13,14):  13/21 close; open ranks 1--8;
E=(2,3,4,5,6,7,8):       7/13 close; open ranks 1--6.
```

Only `25/73` early-rank branches survive.  On this deliberately hostile
battery, the suffix sidecar therefore removes `65.75%` of the post-gate
apex workload before the more expensive complement-cap recursion.

The smallest positive scalar margin is

```text
42457/149557408
```

at rank `9`, apex `78`, for `E=(1,3,9,10,11,12,14)`.  The closest open
margin is `-101/105105` at rank `4`, apex `52`, for the low control.  The
most negative margin is

```text
-2021405909/38818159380
```

at rank `1`, apex `19`, for `E=(2,8,9,10,11,13,14)`.

## 3. Nonmonotonicity

Branch closure is not monotone in the root rank because changing the apex
changes the literal residual.  For `E=(1,3,9,10,11,12,14)`:

```text
rank 5, apex 16:
  margin 5683263/4164400240 > 0,
  residual top-five speeds (15,34,46,26,38);

rank 6, apex 46:
  margin -2298437/841260420 < 0,
  residual top-five speeds (26,32,58,120,78).
```

Thus a uniform rank cutoff would be unsound; the full atlas must profile
each retained suffix branch exactly.

## 4. Tail and reconstruction controls

Every root top thirty and every residual top five is globally sealed using
the strict discrepancy estimate

```text
c_C(w)<|C|/7+(99/70)r/(7w).
```

The largest root threshold is `228324096/168781`.  The largest residual
threshold is `1300178880/469673`, at rank `13`, apex `168`, of
`E=(2,8,9,10,11,13,14)`; its exact scan ends at `2768` and the sealed tail
starts at `2769`.

The locked verifier performs:

```text
16 root vector/scalar controls,
292 residual vector/scalar controls,
73 literal-versus-direct residual reconstructions.
```

Normal and optimized runs are byte-identical.  The complete root/branch
ledger digest is

```text
6510bc4bc977f4c6815fe5024e6f032c0935cf2ce6ec6f9bfdbd17ad6ed91ac7.
```

## 5. Reproduction and scope

```text
04-computation/lrc14_j6_rank_first_suffix_scalar_battery_codex_20260729.py
SHA-256 6434f020c5aa4000ac81fa081881d93ac0b4190516f854fbd9d8493475baf539

05-knowledge/results/lrc14_j6_rank_first_suffix_scalar_battery_codex_20260729.out
SHA-256 d03f4be7ead138135447b4d720e91d215b1668c2182a099da92129289e605ee9
```

The exact result covers only four of the `3432` seven-body roots.  It does
not close the remaining `25` branches, prove the seven-body/six-slot rung,
or prove LRC(14).  Its consequence is architectural: the full successor
should compute scalar and complement caps on the strict root suffix, not on
the conservative universe of every speed other than the selected apex.
