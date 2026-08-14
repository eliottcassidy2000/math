---
id: THM-3373
title: "AMM 12592: the R=8 slack-kernel displacement has exact conformal locality width five"
status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED
source: kps-s180
depends_on:
  - THM-3371
related:
  - THM-3365
companion: 04-computation/amm12592_r8_slack_kernel_locality_kps_s180.py
output: 05-knowledge/results/amm12592_r8_slack_kernel_locality_kps_s180.out
audit: independent exact carrier, Farkas, allocation-equivalence, atom, and all-subset reconstruction
---

# THM-3373 -- R=8 slack-kernel conformal locality width five

## Carrier and definition

Write all `106` causal/terminal inequalities and all `8` initial-coordinate
bound inequalities of the THM-3371 `R=8,D0=0` polytope as

```text
A y <= b,                       A in Z^(114 x 42).
```

Adjoin the slack vector `s=b-Ay`.  Let `F` be the denominator-`103` vertex
and `I` the integer entry certificate of THM-3371.  Then

```text
Delta y = 103(I-F) in Z^42,
Delta s = 103(s_I-s_F) = -A Delta y in Z^114,
w=(Delta y,Delta s) in ker_Z [A I_114].                 (1)
```

The `42` `y`-coordinates are labelled by the seven causal rows `0,...,6`.
For `1<=L<=7`, an **L-window conformal decomposition** of `w` is

```text
w=sum_J (g_J,-A g_J),                                  (2)
```

where each `g_J` is supported on the `y`-coordinates of one interval `J` of
at most `L` consecutive causal rows, every summand lies componentwise in the
same orthant as `w`, and the summands reconstruct `w` exactly.  Slack support
is not discarded: all `114` labelled slack coordinates participate in the
conformality test.

## Theorem

The minimum `L` for which `(2)` exists over the reals, rationals, or integers
is

```text
                              L=5.                      (3)
```

More explicitly:

1. For `L=1,2,3,4`, no real conformal window decomposition exists.
2. At `L=5`, there is a decomposition of `w` into `12` nonzero **integer**
   kernel moves.  Their row windows are

   ```text
   0-1, 0-2, 0-3, 0-4, 1-5, 2-5, 2-6,
   3-5, 3-6, 4-6, 5-5, 5-6.
   ```

   The complete labelled atom vector has semantic hash

   ```text
   36de54cfff43b71efee10d33e6a4df06aa0fa7fea13819571176d5174cf32eb2.
   ```

Thus four consecutive causal rows are provably insufficient even after
relaxing integrality, while five rows suffice without dilation beyond the
already intrinsic factor `103` in `(1)`.

## Proof

For every interval `J` and every `y`-coordinate `j` in that interval, write

```text
g_(J,j)=sign(Delta y_j) z_(J,j),        z_(J,j)>=0.     (4)
```

Coordinate reconstruction gives one equality

```text
sum_(J contains row(j)) z_(J,j)=|Delta y_j|.           (5)
```

For each labelled slack row `q`, conformality of `-A g_J` with `Delta s_q`
is one homogeneous inequality; if `Delta s_q=0` it is an equality.  These
conditions are necessary and sufficient.  Indeed, moves with the same window
can be aggregated because a fixed orthant is a convex cone, while any
decomposition supplies `(4)--(5)`.

The resulting exact allocation systems have respectively

```text
L          1    2    3    4    5
windows    7   13   18   22   25
variables 34   92  166  246  320.
```

For `L=1,...,4`, the companion discovers and then verifies primitive integer
Farkas multipliers `(lambda,mu)` satisfying

```text
lambda E + mu U >= 0,       mu>=0,       lambda d<0.    (6)
```

The numbers of nonzero equality/inequality multipliers, largest absolute
multiplier, and contradictory right-hand side are

```text
L=1: (1, 3,  1, -2150)
L=2: (3, 6,  2,  -605)
L=3: (5,12,  4,  -111)
L=4: (6,27, 12,  -333).
```

All column signs and right-hand sides in `(6)` are checked with Python
integers after the numerical ray has been discarded.  Their certificate
hashes are pinned in the stored output.

For `L=5`, an integer allocation is discovered by MILP, rounded only after
an integrality tolerance check, and then rechecked entirely with integers:
all `42` `y` coordinates and all `114` slack coordinates sum to `(1)`, every
atom has the stated row support, and every nonzero atom coordinate has the
sign of the corresponding coordinate of `w`.  This proves both directions
of `(3)`.  Widths six and seven provide positive hostile controls.

Because each atom is conformal, any partial sum keeps every slack coordinate
between its two nonnegative endpoint values.  Hence the atoms give a monotone
atom-step lattice path inside the `103`-dilate of the augmented polytope.  This
path statement is about the fixed displacement `(1)`, not arbitrary vertices.

An independent audit rebuilt the `114 x 156` augmented carrier, checked every
Farkas multiplier columnwise, reconstructed all `42` `y` and `114` slack
coordinates of the twelve integer atoms, and exhaustively checked all `2^12`
atom subsets.  Every subset has nonnegative slack, so the monotone path claim
holds for every atom ordering.  Three normal and three optimized replays were
byte-identical to the stored output, and all pinned hashes were reproduced.

Normal and optimized runs are byte-identical.  The semantic hash of the
displacement, four Farkas certificates, and integer width-five atoms is

```text
75b989bd93a9f96f521809926fc44972a5e525504d111dc97c4af33d8689e8a5.
```

The LF-normalized source/output hashes are

```text
745fd7bbe1b1acbac915e9b068f9ac7931e6705466711005f42375c946663e9c
b75817fe24a1f587aea298863b12b90d5880ea94186d788f3dffcd6b42923c2f.
```

## Interpretation and boundary

This makes the toric transfer from THM-3365 precise but sharply limited.
The augmented kernel preserves all causal equations and inequality slacks;
the row labels and slack orthant are the mandatory sidecar.  In contrast, the
bare toric kernel forgets the right-hand side, capacities, golden profile,
feed chronology, and THM-3332 entry predicate.

- The twelve atoms are not claimed primitive Graver elements.
- Width five is exact only for the fixed `F`-to-`I` displacement at `R=8`.
- Translation of these atom shapes to `R=512` must separately preserve the
  two degree-step types, capacities, integrality, and entry-cone margins.
- No `R=512` point, general rounding theorem, or AMM deadline bound follows.
