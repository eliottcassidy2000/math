# HYP-2453 - Triangular tower moment bridges expose fractional address carriers

**Status:** OPEN synthesis; exact tower identities, overlap scan, and
higher-moment leading address verified.
**Source:** codex-2026-06-12.
**Companions:** HYP-2452, HYP-2451, HYP-2444, HYP-2443, HYP-2438,
HYP-2436, THM-492, THM-491, THM-406.
**Computation:** `04-computation/triangular_tower_moment_bridge_codex.py`;
stored output `05-knowledge/results/triangular_tower_moment_bridge_codex.out`.
**Addendum computation:** `04-computation/triangular_tower_overlap_families_codex.py`;
stored output `05-knowledge/results/triangular_tower_overlap_families_codex.out`.
**Beatty/Pell addendum:** HYP-2456;
`04-computation/triangular_tower_beatty_pell_decomposition_codex.py`.
**External anchors:** OEIS A059270 and A059255.

## Statement

The user's two triangular towers are best read as adjacent unequal block
splits of one shell, decorated by the lowest moment that equalizes the split.

The first tower is

```text
n^2 + ... + (n^2+n)
  = (n^2+n+1) + ... + (n^2+2n).
```

It splits the square shell

```text
[n^2, (n+1)^2-1]
```

into `n+1` left terms and `n` right terms with equal first moment.

The second tower is

```text
T_{2n}^2 + ... + (T_{2n}+n)^2
  = (T_{2n}+n+1)^2 + ... + (T_{2n}+2n)^2.
```

It splits the triangular shell

```text
[T_{2n}, T_{2n+1}-1]
```

into the same `n+1` versus `n` counts, but with equal second moment.

Thus the addition-versus-multiplication bridge is not only symbolic.  The
first tower has a hidden square-sum shadow:

```text
A_common(n) = n(n+1)(2n+1)/2 = 3 * sum_{j<=n} j^2.
```

The second tower has an unsquared first-moment defect

```text
B_raw_left(n) - B_raw_right(n) = n(n+1) = 2T_n,
```

and squaring is exactly the moment lift that cancels that defect.

## Exact Evidence

The stored computation verifies the first rows:

```text
n=2 A: (4,6)=(7,8), sum=15;
    B^2: (10,12)=(13,14), sumsq=365, raw_defect=6

n=3 A: (9,12)=(13,15), sum=42;
    B^2: (21,24)=(25,27), sumsq=2030, raw_defect=12

n=4 A: (16,20)=(21,24), sum=90;
    B^2: (36,40)=(41,44), sumsq=7230, raw_defect=20
```

Closed forms:

```text
A_common(n) = n(n+1)(2n+1)/2
B_common_square(n) = n(n+1)(2n+1)(12n^2+12n+1)/6
B_raw_left(n) = n(n+1)(4n+3)/2
B_raw_right(n) = n(n+1)(4n+1)/2
```

The overlap scan through `n<=80` finds exactly one whole-side equality:

```text
B_3.L = A_4.R = [21,24].
```

This is not just a sampled coincidence.  Whole-side equality `B_L(n)=A_R(m)`
forces `m=n+1` and

```text
n^2 - 2n - 3 = 0,
```

so `n=3`; the other length/start equations have no positive solution.  The
user's `21+22+23+24` hinge is therefore unique among exact side intervals.

## Overlap Families

The additive tower covers all positive integers because

```text
A_n = [n^2,(n+1)^2-1]
```

partitions `[1,infty)`.  The square tower does not.  It covers only alternating
triangular shells

```text
B_m = [T_{2m},T_{2m+1}-1],
```

skipping

```text
[T_{2m+1},T_{2m+2}-1]
```

of size `2m+2` between them.

The follow-up overlap scout separates three useful families.

### Whole Equation Side-Aligned

The user's `10+11+12 = 13+14` example is the first nontrivial member of an
infinite Pell family:

```text
B_m.L subset A_n.L and B_m.R subset A_n.R
  iff T_n = 2T_m
  iff x^2 - 2y^2 = -1, x=2n+1, y=2m+1.
```

Starting from the trivial `(n,m)=(0,0)`, the recurrence is

```text
n' = 3n + 4m + 3
m' = 2n + 3m + 2.
```

The first nontrivial pairs are

```text
(m,n) = (2,3), (14,20), (84,119), (492,696), ...
```

For each pair, `B_m` is the middle of `A_n`; the outside padding on both ends
is `n-m`, and the B size is `2m+1` inside A size `2n+1`.

### Exact Whole-Side Equality

The less restrictive `21+22+23+24` pattern is an exact whole-side equality,
not a whole-equation containment:

```text
B_3.L = A_4.R = [21,24].
```

It is unique.  Length and start equations show `B_L(m)=A_R(n)` forces
`n=m+1` and `m^2-2m-3=0`, hence `m=3`; the other side-pair equations have no
positive solution.

### All Looser Single-Side Matches

All single-side containment instances have a one-line classifier.  For a B
side interval `I=[u,v]`, let

```text
n = floor(sqrt(u)).
```

Then:

```text
if v > (n+1)^2-1:  I crosses an A square-shell boundary;
elif v <= n^2+n:   I is inside A_n.L;
elif u >= n^2+n+1: I is inside A_n.R;
else:              I crosses the A left/right midpoint.
```

This predicts every looser match and reports its part of A plus padding.  The
side-containment word is Sturmian/Beatty-like because the start `T_{2m}` drifts
through square shells with slope `sqrt(2)`.

## Higher-Moment Address

For power `p`, ask for a start `A` satisfying

```text
sum_{i=0}^n (A+i)^p = sum_{i=0}^{n-1} (A+n+1+i)^p.
```

The first two powers give exact integer polynomial starts:

```text
p=1: A = n^2
p=2: A = 2n^2+n = T_{2n}
```

For `p>=3`, the balancing start is no longer integral.  It has the leading
form

```text
A_p(n) = p n^2 + (p-1)n + (p-1)(p-2)/(12p) + lower-order terms.
```

The computation verifies the leading cancellation for `p=3..8`; at `n=30`,
the positive roots have offsets

```text
p=3: 0.055553...  -> 1/18
p=4: 0.124991...  -> 1/8
p=5: 0.199983...  -> 1/5
p=6: 0.277751...  -> 5/18
p=7: 0.357106...  -> 5/14
p=8: 0.437453...  -> 7/16
```

This suggests the exact towers are the first two integral faces of a
Faulhaber/Bernoulli tower.  After that, a missing fractional address must be
attached before the scalar shell equality can be lifted.

## Transfer To LRC14

This does not prove LRC14.  It does sharpen a carrier principle.

HYP-2443/HYP-2444/HYP-2452 all say that scalar ledgers are too coarse:
`q` blocked, a coefficient row, or a weight enumerator is only a boundary
total.  The hard object is the hidden lift that remembers owners, carries,
fibers, supports, or factor-grid cells.

The moment-tower version says the same thing in one dimension:

```text
shell interval        -> scalar support
moment equality       -> observable boundary total
fractional address    -> missing lift coordinate
```

The `n=3` second-moment row is especially suggestive because its full
triangular shell is `[21,27]`, with unique exact hinge `[21,24]` against the
first-moment square shell and right cap `[25,27]`.  The numbers `21` and `27`
already recur in the repo as the unit-distance/LRC14 difficulty bridge and the
shell `2*14-1=27`.  The safe claim is only that this gives a new bookkeeping
language: LRC14 blockers should be tested not merely by scalar denominator
coverage, but by moment/resource defects and the fractional or fiber address
needed to lift a blocked shell into an actual support proof.

## Tournament Analysis

The computation uses proof routes as vertices:

```text
moment1_square_shell
moment2_triangular_shell
bernoulli_fraction_address
unique_hinge_overlap
convolution_factor_lift_bridge
lrc_resource_ledger
raw_oeis_sequence_match
```

Pairwise observable: majority win over exactness, computability, hidden-lift
value, LRC transfer, fractional-address value, and sequence anchoring.

Switch/gauge: orient toward the route retaining more proof-bearing
side-channel data.  Tie Hamiltonian path: listed order, then total score.

Fingerprint:

```text
score_hist = {0:1, 1:1, 3:3, 5:1, 6:1}
directed_3cycles = 1
scc_sizes = [3,1,1,1,1]
hamiltonian_paths = 3
leader = convolution_factor_lift_bridge
```

The nontransitivity is useful: the square-shell and triangular-shell exact
identities are cleaner, but the convolution/fractional-address routes retain
more hidden-lift data.

## Assumption Challenge

Alternate vertex sets considered: tower equations, block sides, shell gaps,
moment powers, fractional starts, OEIS sequences, LRC denominators, runner
resources, hidden convolution grid cells, and proof-route obligations.

Chosen quotient: adjacent unequal block splits of one shell, decorated by the
first moment that equalizes.

Preserved: exact shell interval geometry, first/second moment equality,
overlap/hinge data, and the fractional carry needed at higher powers.

Destroyed: runner positions, actual LRC circle arcs, Galois/root data, and
large-scale asymptotics beyond the leading Bernoulli address.  These must be
reattached before any proof claim.

Challenged assumption: the useful LRC clock need not be time `t` or runner
speed.  Here the clock is the moment order `p` and the hidden address that
makes two adjacent unequal blocks balance.

## Next Moves

1. Prove the higher-moment expansion to several terms and identify when
   Bernoulli denominators can be killed by changing shells or adding address
   coordinates.
2. Extend the HYP-2456 Beatty/Pell address system into a Q27-style LRC14
   ledger: shell address, owner support, carry residue, divisor fiber, and
   boundary atom.
3. Add moment/resource defects to the LRC14 Q27 blocker ledger: blocked twists,
   owners, divisor fiber, and the fractional address required to lift a scalar
   obstruction.
4. Compare with HYP-2452: coefficient convolution lifts and moment towers are
   both boundary-total shadows of hidden filtered support.
