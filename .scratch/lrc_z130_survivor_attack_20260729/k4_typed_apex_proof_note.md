# k=4 ordered typed-apex reduction: proof note

Status: proof draft for audit.  The exact census is being replayed under
ordinary, optimized, and one-worker execution before promotion.

If the incoming Lorenz/activity/needle argument closes the whole mixed
three-drift sector, that later theorem supersedes this result only as an
emptiness route.  The ordered typed-apex theorem should still be preserved:
it is an independent physical-order finite reduction, exposes the exact
integer-apex mechanism, and supplies a cross-check that does not pass through
the denominator quotient.

## Statement

Fix a six-body root `E`, let `L=14*lcm(E)`, and suppose a literal seven-tail
cover has four aligned tails `La_1,...,La_4` and three nonaligned drifts

```text
15 <= z1 < z2 < z3,             L does not divide z_i.
```

The projected-suffix theorem restricts `(E,z1)` to its exact `87,975`-row
bank, on which `z1<=182`.  For every such row, the ordered recursion below
forces `z2<=6515`.  The multiplicity-excess scalar condition then leaves
exactly

```text
2,042,669 finite-z3-cap pairs
6,471,505 pairs with no positive scalar gap
--------------------------------------------
8,514,174 necessary (E,z1,z2) pairs,
```

with the sharper uniform cap `z2<=2163`.  These are necessary physical pair
prefixes, not realized covers.

## Exact non-strict integer-apex lemma

Let `R` be a positive-mass union of `r` intervals, of mass `h`, and suppose
`p<=6` integer labels cover `R`.  The component estimate is

```text
c_R(w) <= h/7 + gamma/w,             gamma=6r/49.
```

Put

```text
T = 7p*gamma/((7-p)h)
  = 6pr/[7(7-p)h],
A_p(R)=floor(T).
```

If all `p` integer labels were greater than `A_p(R)`, each would be strictly
greater than `T`, including when `T` itself is integral.  Hence every
coverage would be strictly less than `h/7+gamma/T`, and their sum would be
strictly less than `h`, contradicting coverage.  Therefore some label is at
most `A_p(R)`.

This direct non-strict argument is load-bearing.  Replacing `gamma` by
`gamma+1` to invoke a strict tail lemma proves finiteness but does not prove
the exact cap used here.

## Disjoint ordered recursion

Delete `D_z1` and call the literal residual `R_0`.  Six tail labels remain:
four aligned labels and `z2,z3`.  Apply the lemma with `p=6`.

1. If `z2<=A_6(R_0)`, record the bounded physical branch.
2. Otherwise both drifts exceed the apex.  The least remaining actual
   aligned multiplier must be at most `floor(A_6(R_0)/L)`.  Branch on that
   multiplier and record that all smaller multipliers are absent.
3. Delete its aligned danger set, decrement `p`, and repeat.  On a complement
   branch update the physical lower bound to
   `z2>=max(previous lower bound,A_p(R)+1)`.
4. Selected multipliers are strictly increasing.  Thus each hypothetical
   aligned set follows exactly one path: the recursion does not quotient
   distinct physical orders or count permutations.
5. After four aligned deletions only `z2,z3` remain.  The bounded `p=2`
   branch records `z2<=A_2(R)`; its complement is impossible because no
   aligned label remains to pay the apex.

The exact tree has `5,368` states at selected depths

```text
depth                 1      2      3      4
states              758   1052   1427   2131
bounded branches     758    990   1293    543
closed complements     0     54    684   2131.
```

The merged necessary `z2` set for each of the `87,975` suffix rows is,
unexpectedly, one prefix interval beginning at `z1+1` (aligned multiples of
`L` are omitted from the physical candidate count).  The raw bank contains
`50,285,016` ordered physical pairs and has cap `6515`.  Its unique cap row is

```text
E=(2,3,4,8,9,12), z1=156, z2=6515, L=1008.
```

## Projected wall and scalar screen

For `k=4`,

```text
alpha_4=2366/21875.
```

Rootwise `alpha_4 L>max(E)`; its minimum is

```text
alpha_4*168=56784/3125>12
at E=(1,2,3,4,6,12).
```

Consequently the projected-safe-arc wall bounds the largest drift itself,
not merely `max(E,z3)`:

```text
z3 >= high_floor=floor(alpha_4 L)+1.
```

The multiplicity surplus gives

```text
delta(z1)+delta(z2)+delta(z3) >= eta_4 h,
eta_4=51/1183.
```

For a fixed physical prefix put

```text
g=eta_4 h-delta(z1)-delta(z2),
Z0=max(z2+1,high_floor).
```

If `g>0`, the component bound forces

```text
z3 <= U=floor((6r/49)/g).
```

Thus `U<Z0` kills the prefix and otherwise gives an exact finite `z3`
interval.  If `g<=0`, this scalar inequality supplies no positive upper gap
for `z3`; the prefix is conservatively retained.  In particular this second
class does **not** assert that `delta(z3)>=0` or that a cover exists.

The `50,285,016` raw pairs split as

```text
41,770,842 scalar-killed
 2,042,669 finite z3 cap
 6,471,505 no positive scalar gap.
```

The surviving cap `z2=2163` is unique:

```text
E=(2,3,4,6,8,12), z1=44, L=336,
z3>=2164, U=23477.
```

## Information boundary and next composition

Merging branches deliberately forgets which aligned prefix produced a
bounded interval.  Therefore the result is a necessary physical pair bank,
not a census of realized aligned completions.  The next lawful screen joins
the physical order to the denominator sidecar

```text
d(z)=L/gcd(z,L).
```

Sorted denominator triples must be treated as multisets independent of the
physical order `z1<z2<z3`.  On finite-gap pairs, denominator existence alone
is insufficient: one must also exhibit an integer `z3` in `[Z0,U]` having
the remaining denominator.
