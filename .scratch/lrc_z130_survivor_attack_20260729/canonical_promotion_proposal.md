# Canonical promotion proposal: THM-2941 k=4 ordered typed apex

Do not apply until the concurrent whole-`k=4` closure commit has landed and
its theorem scope is known.

## Exact tracked paths

Copy the audited scratch verifier byte-for-byte to

```text
04-computation/lrc14_j7_four_aligned_three_drift_typed_apex_thm2941.py
```

and its final ordinary output byte-for-byte to

```text
05-knowledge/results/lrc14_j7_four_aligned_three_drift_typed_apex_thm2941.out
```

Then add those two paths to THM-2941's YAML `verification` list.  The
verifier already discovers the repository root correctly from either its
scratch or canonical location, and its default output is the canonical path
above.

Reproduction commands after promotion:

```bash
env PYTHONHASHSEED=0 python3 \
  04-computation/lrc14_j7_four_aligned_three_drift_typed_apex_thm2941.py \
  --workers 16 --output /tmp/k4-ordinary.out
env PYTHONHASHSEED=271828 python3 -O \
  04-computation/lrc14_j7_four_aligned_three_drift_typed_apex_thm2941.py \
  --workers 16 --output /tmp/k4-opt.out
env PYTHONHASHSEED=314159 python3 \
  04-computation/lrc14_j7_four_aligned_three_drift_typed_apex_thm2941.py \
  --workers 1 --output /tmp/k4-one-worker.out
cmp /tmp/k4-ordinary.out /tmp/k4-opt.out
cmp /tmp/k4-ordinary.out /tmp/k4-one-worker.out
```

Final frozen controls before promotion:

```text
semantic digest:
  5303930eea27eb9713a81245efa0a975f7e2b21869edcd028e4b2f2bb56f8fb3
cap-witness digest:
  9c3469a3600a928aebcdae07d51a1c3f16fc221740505502250e59ca716efefe
cached-excess hostile-audit digest:
  95c91236dfb3e99ac2a28b4e030e1397dfcf046e5debda329a801f1a6300f544
source LF SHA-256:
  b5bd41551f12f9c781afb8444a7ac64c0eca21a8f64426e0dece7909605bcd4c
output SHA-256:
  26234efe820a2a045795895bcae33ddd4417d532a0c384c4a6b516d173bd7d99
```

## Proposed THM-2941 subsection

### Four aligned tails: an ordered physical finite reduction

The `k=4` projected-suffix bank admits a sharper physical-order reduction
than the general finite-tree statement.  Suppose a literal cover has four
aligned tails and three nonaligned drifts

```text
z_1<z_2<z_3.
```

The projected suffix leaves exactly `87,975` `(E,z_1)` rows on all `3,003`
roots, with `z_1<=182`.  For a positive-mass residual `R` of mass `h` and
`r` components, if `p<=6` integer labels cover `R`, the non-strict component
bound

```text
c_R(w)<=h/7+(6r/49)/w
```

forces an exact integer apex

```text
A_p(R)=floor(6pr/[7(7-p)h]).                           (K4-1)
```

Indeed, put `T=6pr/[7(7-p)h]`.  If every integer label exceeded
`floor(T)`, then every label would be strictly greater than `T`, including
when `T` is integral.  Each coverage would be strictly below
`h/7+(6r/49)/T`, so the `p` coverages would sum to strictly less than `h`,
a contradiction.  This direct non-strict lemma is sharper than enlarging
`6r/49` in order to invoke a strict tail statement.

Delete `D_(z_1)` and apply `(K4-1)` to the four aligned labels and
`z_2,z_3`.  If `z_2<=A_6`, record that bounded branch.  Otherwise the least
remaining actual aligned multiplier is at most `floor(A_6/L)`; branch on
that multiplier, record the forbidden smaller-multiplier prefix, delete its
danger set, decrement `p`, and repeat.  Selected multipliers are strictly
increasing, so this is a disjoint recursion, not a permutation quotient.
After four aligned deletions only `z_2,z_3` remain, and the complement of
the bounded `p=2` branch is impossible.

The exact recursion has `5,368` states:

```text
selected depth          1      2      3      4
states                758   1,052  1,427  2,131
bounded z_2 branches   758     990  1,293    543
closed complements       0      54    684  2,131.
```

For every one of the `87,975` suffix rows, the union of all necessary
bounded branches merges to one physical prefix interval beginning at
`z_1+1`, with aligned multiples of `L` omitted.  It contains `50,285,016`
row-labelled physical pairs and has the raw uniform cap `z_2<=6,515`.

The projected wall can be assigned to the third drift itself.  Rootwise

```text
alpha_4 L>max(E),             alpha_4=2366/21875,
```

with minimum `56784/3125>12` at `E=(1,2,3,4,6,12)`.  Therefore

```text
z_3>=H_E:=floor(alpha_4 L)+1.                          (K4-2)
```

The multiplicity surplus also requires

```text
delta(z_1)+delta(z_2)+delta(z_3)>=(51/1183)h.
```

For a fixed prefix put

```text
g=(51/1183)h-delta(z_1)-delta(z_2),
Z_0=max(z_2+1,H_E).
```

If `g>0`, the component bound forces the exact finite interval

```text
Z_0<=z_3<=floor((6r_E/49)/g).                          (K4-3)
```

If its upper endpoint lies below `Z_0`, the prefix is impossible.  If
`g<=0`, the scalar inequality supplies no positive `z_3` gap and the prefix
is conservatively retained; this does not assert that `delta(z_3)` is
nonnegative or that a cover exists.

The raw pairs split exactly into

```text
41,770,842  scalar-killed,
 2,042,669  with a finite nonempty z_3 interval,
 6,471,505  with no positive scalar gap.
```

Hence every literal four-aligned/three-drift cover has its ordered physical
prefix in an `8,514,174`-pair exact bank and satisfies the sharper uniform
cap

```text
z_2<=2,163.                                             (K4-4)
```

The raw cap has the unique witness

```text
E=(2,3,4,8,9,12), z_1=156, z_2=6,515, L=1,008,
```

which is scalar-killed.  The screened cap has the unique witness

```text
E=(2,3,4,6,8,12), z_1=44, z_2=2,163, L=336,
2,164<=z_3<=23,477.
```

Merging aligned-prefix branches forgets which aligned set generated an
interval.  Thus `(K4-4)` is a necessary physical-order finite reduction,
not a realized-cover census.  The denominator sidecar
`d(z)=L/gcd(z,L)` is independently unordered and cannot be identified with
the physical order `z_1<z_2<z_3`.

If the later Lorenz/activity/needle theorem closes the whole mixed
three-drift sector, cite that theorem immediately after this subsection:
it supersedes this route for emptiness, but not as an independent
physical-order reduction or as a cross-check of the denominator quotient.

The verifier checks all frozen censuses and both unique cap witnesses.  It
also makes `36,036` direct hostile comparisons—twelve boundary/interior
labels on every root—between the singleton-derived cached excess and an
independent exact numerator formula.  Ordinary, optimized, and one-worker
replays are byte-identical.

## Proposed theorem status/YAML wording

Append to the existing status paragraph:

```text
An ordered non-strict integer-apex recursion reduces the four-aligned/
three-drift face to an exact 8,514,174-prefix physical bank with
z_2<=2,163; this is an independent finite reduction, not by itself an
emptiness theorem.
```

If the whole `k=4` closure is already canonical, instead use:

```text
The now-empty four-aligned/three-drift face also has an independent ordered
integer-apex reduction to an exact 8,514,174-prefix physical bank with
z_2<=2,163; the later closure supersedes it only as an emptiness route.
```

## Proposed `05-knowledge/results/INDEX.md` entry

```text
| `lrc14_j7_four_aligned_three_drift_typed_apex_thm2941.out`
| [PROVED + FINITE-EXACT + VERIFIED; THM-2941]
| Ordered non-strict integer-apex recursion on the exact 87,975-row k=4
projected suffix.  Its 5,368 aligned-prefix states merge rowwise to one
physical z2-prefix; 50,285,016 raw pairs split into 41,770,842 scalar kills,
2,042,669 finite-z3 pairs, and 6,471,505 no-positive-gap pairs.  Therefore
every hypothetical four-aligned/three-drift cover lies in an exact
8,514,174-pair bank with z2<=2,163.  Unique cap witness
E=(2,3,4,6,8,12), z1=44, L=336, z3 in [2,164,23,477].
Necessary physical finite reduction, not a realized-cover census; preserve
as an independent cross-check if a later theorem closes the whole sector.
|
```

## Navigation wording after whole-sector closure

In `00-navigation/CURRENT-FRONTIER.md`, do not leave the old statement that
the mixed `k=4` census is open.  Route the actual emptiness theorem first,
then add a short independent-mechanism clause:

```text
THM-2941 also gives an ordered physical finite reduction of this sector:
8,514,174 necessary (E,z1,z2) prefixes and z2<=2,163, independently of the
denominator quotient.
```
