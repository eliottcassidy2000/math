# HYP-2456 - The triangular crossover word has an exact Beatty-Pell normal form

**Status:** OPEN synthesis with exact finite classifier/proof scaffold verified.
**Source:** codex-2026-06-13.
**Companions:** HYP-2455, HYP-2454, HYP-2453, HYP-2444, HYP-2443,
HYP-2452, HYP-2241, HYP-2239.
**Computation:** `04-computation/triangular_tower_beatty_pell_decomposition_codex.py`;
stored output `05-knowledge/results/triangular_tower_beatty_pell_decomposition_codex.out`.

## Statement

HYP-2453 left the side-containment word between the square-shell tower `A`
and triangular tower `B` as "Sturmian/Beatty-like."  HYP-2455 then reframed
these phenomena as scalar-shadow/hidden-lift problems.  The sharper statement
for this exact toy is:

```text
the shell address is Beatty,
the address increment word is Sturmian,
the visible two-side crossover word is a carry-decorated six-interval
rotation with two zero-density Pell boundary atoms.
```

Let

```text
A_n   = [n^2, (n+1)^2-1],
A_n.L = [n^2, n^2+n],
A_n.R = [n^2+n+1, n^2+2n],

B_m.L = [T_{2m}, T_{2m}+m],
B_m.R = [T_{2m}+m+1, T_{2m}+2m].
```

The shell index of `B_m.L` is exactly

```text
n_m = floor(sqrt(2m^2+m))
    = floor(sqrt(2)*(m+1/4)).
```

Set

```text
d_m = n_m - m
    = floor((sqrt(2)-1)m + sqrt(2)/4),

r_m = m^2 - 2m*d_m - d_m^2.
```

Then the whole two-side state of row `m` is determined by the integer pair
`(d_m,r_m)`; no floating rotation approximation is needed.

## Exact State Formula

For `B_m.L`:

```text
L iff r <= d-m
M iff d-m < r <= d
R iff d < r <= 2d
S iff r > 2d
```

where `L/R/M/S` mean, respectively:

```text
L = contained in A_n.L,
R = contained in A_n.R,
M = crosses the A_n left/right midpoint,
S = crosses the square-shell boundary after A_n.
```

For `B_m.R`, let

```text
epsilon = 1 iff r >= 2d.
```

This is exactly the condition that `B_m.R` starts in `A_{n+1}` rather than
`A_n`.

If `epsilon=0`, then

```text
L iff r <= d-2m
M iff d-2m < r < d-m
R iff d-m <= r <= 2d-m
S iff 2d-m < r < 2d.
```

If `epsilon=1`, then

```text
L iff r <= 3d+2
M iff r > 3d+2.
```

The computation verifies this normal form against the direct floor-sqrt
classifier through `m=250000`.

## Why The Address Is Beatty

Let

```text
X = sqrt(2)*(m+1/4).
```

Then

```text
X^2 = 2m^2+m+1/8 = T_{2m}+1/8.
```

If `n=floor(X)`, then `n^2` is an integer at most `T_{2m}+1/8`, hence
`n^2 <= T_{2m}`.  Also `(n+1)^2 > X^2 > T_{2m}`.  Therefore

```text
n = floor(sqrt(T_{2m})).
```

Subtracting `m` gives the inhomogeneous Beatty sequence for `d_m`.

Consequently the first-difference word `d_m-d_{m-1}` is the genuine Sturmian
clock.  The stored factor-complexity sample gives

```text
p(k)=k+1 for k=1..20
```

for this binary increment word.

## The Visible Word Is Not Sturmian

The two-side token word is richer than the underlying Beatty increment.
Let

```text
eta_m = {sqrt(2)*(m+1/4)}.
```

Away from zero-density equality walls, the exact carry inequalities collapse
to a six-interval rotation coding:

```text
eta in [0, (2-sqrt(2))/4)                         -> LM
eta in [(2-sqrt(2))/4, (2-sqrt(2))/2)             -> MR
eta in [(2-sqrt(2))/2, 1/2)                       -> MS
eta in [1/2, 1-1/(2sqrt(2)))                      -> RS
eta in [1-1/(2sqrt(2)), (3sqrt(2)-2)/(2sqrt(2)))  -> SL
eta in [(3sqrt(2)-2)/(2sqrt(2)), 1)               -> SM.
```

The limiting lengths are

```text
LM, MR, RS, SL: (2-sqrt(2))/4
MS, SM:         (sqrt(2)-1)/2
LR, RL:         0
```

and the finite run through `m=500000` matches these lengths within about
`2e-5`.  However the token word has much higher factor complexity:

```text
token p(k), k=1..20 =
  [8,17,29,40,51,62,71,80,89,98,107,116,124,132,140,148,156,164,172,180].
```

Thus the correct slogan is not "the crossover word is Sturmian."  It is:

```text
the address increment is Sturmian, while the observer-visible token word is
a finite-interval/carry decoration of that Sturmian clock.
```

## Pell Boundary Atoms

Two rare token types are exact equality walls rather than positive-density
rotation intervals.

The `LR` atom is the HYP-2453 whole-equation side-alignment family:

```text
r = d-m
iff (2n+1)^2 - 2(2m+1)^2 = -1.
```

The first `m` values are

```text
2, 14, 84, 492, 2870, 16730, 97512, 568344, ...
```

with recurrence

```text
m_{k+2} = 6m_{k+1} - m_k + 2.
```

The `RL` atom is the adjacent endpoint wall where `B_m.R` starts exactly at
the next square:

```text
r = 2d
iff 2(n+1)^2 - (2m+1)^2 = 1
iff (2(n+1))^2 - 2(2m+1)^2 = 2.
```

The first `m` values are

```text
3, 20, 119, 696, 4059, 23660, 137903, 803760, ...
```

with the same recurrence.  The unique whole-side equality
`B_3.L=A_4.R=[21,24]` is the first member of this `RL` wall; later members are
endpoint alignments with padding, not whole-side equalities.

## Transfer To LRC14

This is a compact exact model for HYP-2455's repo-wide pattern:

```text
scalar/product collision -> attach address/carry coordinate -> unroll into a word.
```

For LRC14, this suggests replacing a scalar row such as "q is blocked" by a
live address tuple, for example:

```text
(denominator shell, unit quotient class, divisor fiber,
 owner support, carry residue, boundary atom type, opening/deletion target).
```

The point is not that LRC14 literally has slope `sqrt(2)`.  The point is that
the triangular toy shows how an apparently vague Sturmian/Beatty shadow can be
made exact only after keeping the Pell/carry remainder.  That is the same
information type HYP-2241 and HYP-2443 keep rediscovering for Q27: visible
scalar shadows mix strict and wall rows until owner/carry/private-deletion
addresses are restored.

## Tournament Analysis

The computation uses proof-carrier vertices, not interval entries:

```text
beatty_shell_address
sturmian_increment_clock
pell_wall_atoms
six_interval_rotation
exact_carry_normal_form
token_complexity_warning
lrc14_address_ledger
higher_moment_fractional_address
```

Pairwise observable: majority comparison over exactness, clock clarity,
hidden-address value, computability, LRC transfer, and proof potential.

Switch/gauge: orient toward the carrier retaining more hidden-address data;
tie Hamiltonian path is the listed order.

Stored fingerprint:

```text
score_hist = {0:1, 2:3, 4:1, 5:1, 6:1, 7:1}
directed_3cycles = 1
scc_sizes = [3,1,1,1,1,1]
hamiltonian_paths = 3
leader = exact_carry_normal_form
```

## Assumption Challenge

Alternate vertex sets considered: B rows, A shells, side intervals, square
boundaries, midpoint crossings, Beatty increments, Pell equality walls, carry
states, token factors, and LRC support obligations.

Chosen quotient: `(d_m,r_m,epsilon)` as the hidden address under the visible
two-side token.

Preserved: exact interval containment state, square-shell address, midpoint
crossing, boundary crossing, and endpoint equality atoms.

Destroyed: full moment sums, higher-power balance roots, actual LRC runner
geometry, and support-owner data.  Those must be reattached in the LRC transfer.

Challenged assumption: seeing a low-complexity Beatty/Sturmian trace does not
mean the observable word itself is Sturmian.  The proof-bearing coordinate may
be one layer lower than the visible token.

## Next Moves

1. Prove the six-interval model rigorously from the exact `(d,r)` inequalities
   and equidistribution of `eta_m`.
2. Add a `Q27` address-ledger analog: for each LRC14 blocker row, store the
   shell address, owner support, carry class, divisor fiber, and boundary atom
   instead of only the blocked denominator.
3. Compare HYP-2456's `LR/RL` zero-density atoms with AP/Vstar/2AP and
   HYP-2241's owner-private deletion rows: are LRC wall atoms also Pell-like
   equality events inside a higher-dimensional address word?
4. Use this as a toy regression test for any future "Sturmian" claim: identify
   which projected coordinate has `p(k)=k+1`, and which visible word is merely
   a decoration.
