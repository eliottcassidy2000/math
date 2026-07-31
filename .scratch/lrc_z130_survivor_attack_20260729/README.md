# LRC(14) low-drift phase-arc obstruction

**Status: scratch mathematics.  This is not a proof of LRC(14).**

## 1. The conservative full-cell wall

Let `E` be six body speeds, let

```text
L=14 lcm(E),
```

and suppose a seven-tail countercover has `k` distinct aligned speeds
`La_1,...,La_k` and `d=7-k` drift speeds.  Put

```text
M=max(E union {the d drift speeds}).
```

The partial body-plus-drift family has `6+d=13-k` speeds.  Settled
`LRC(13-k)` supplies a time `t_0` at clearance at least `1/(14-k)`.
The distance-to-integer functions are `M`-Lipschitz, so the closed arc

```text
I={t: |t-t_0| <= k/[14(14-k)M]}
```

is safe for every body and drift speed at level `1/14`.

If

```text
M <= kL/[14(14-k)],
```

then the radius of `I` is at least `1/L`.  Every closed arc of length at
least `2/L` contains a complete closed `1/L` grid cell.  On that cell all
body and drift combs are safe.  After `t=(j+u)/L`, the aligned combs become
the `k` normalized combs `D_(a_i)`.  Their union has measure at most
`k/7<1`, so they cannot cover the complete cell.  Contradiction.

Thus every countercover on this face obeys

```text
M > kL/[14(14-k)].                                    (FC)
```

Since direct-boundary body speeds are at most `14` and external drifts are
at least `15`, the exact integer largest-drift wall is

```text
z_d >= max(15, floor(kL/[14(14-k)])+1).                (FC-int)
```

The equality case in `(FC)` is excluded: the safe arc is closed and danger
uses the strict predicate `<1/14`.

## 2. The stronger projected-arc capacity wall

A full cell is more than is needed.  Let

```text
pi_L(t)=Lt mod 1,
P=pi_L(I).
```

The arc `I` has length

```text
ell=k/[7(14-k)M],
```

so

```text
mu(P)=min(1,L ell).
```

Every point of `I` is safe for the body and drift speeds.  A hypothetical
countercover therefore forces

```text
P subset U_A:=union_i D_(a_i).
```

The set `P` is compact, `U_A` is open, and the proved aligned safe floor
`u_A>0` makes `U_A` a proper subset of the circle.  Hence compact
containment is strict in measure:

```text
mu(P)<mu(U_A)<=1-u_A.
```

Using any uniform floor `u_A>=u_k` yields the stronger necessary wall

```text
M > kL/[7(14-k)(1-u_k)].                              (PA)
```

For the currently proved floors this gives

```text
k   u_k          coefficient alpha_k in M>alpha_k L
3   55/91        13/132
4   558/1183     2366/21875
5   478/1365     2275/18627.
```

This argument already includes the wrap case: if `L ell>=1`, then
`P` is the whole circle, impossible because `u_A>0`.

## 3. Uniform elimination of the `z_1=130` frontier

The exact excess suffix census left four `k=5` frontier rows:

```text
E=(1,4,8,10,12,13),    z_2=156, L=21840;
E=(1,6,8,10,11,14),    z_2=182, L=129360;
E=(1,6,8,10,12,13),    z_2=156, L=21840;
E=(2,6,8,10,11,14),    z_2=182, L=129360.
```

Here `M=z_2`.  Already the conservative condition gives

```text
5L>126M
```

in all four rows, contradicting `(FC)`.  The elimination is independent of
the aligned multiplier set, its phases, and whether the excess-census gap
is zero or positive.

As hostile exact controls, complete body-grid cells safe from both drifts
occur at

```text
L=21840,  z_2=156: j=1970;
L=129360, z_2=182, E=(1,6,8,10,11,14): j=10022;
L=129360, z_2=182, E=(2,6,8,10,11,14): j=4620.
```

Some left endpoints have clearance exactly `1/14`, which is safe because
the danger comb is open.  The conceptual proof above does not depend on
these witnesses.

## 4. Exact constrained suffix scan

`fullcell_constrained_suffix_scan.py` composes `(FC-int)` with the exact
singleton-excess suffix bank for `k=3,4,5`.  It integrates every label
through `7000` exactly and bounds omitted high labels by

```text
delta(w)<=6r/(49w)
```

starting at the larger of `7001` and the full-cell wall.  Its baseline
counts and maxima are locked against the prior suffix audit before any
constrained result is printed.

The stronger projected wall leaves the exact first-drift census

```text
k                              2          3        4      5     6
projected surviving rows 2,239,853    376,020   87,975  4,702     0
projected maximum z_1        2,142        380      182     66  EMPTY.
```

These are necessary-filter rows `(E,z_1)`, not realized covers.

## 5. The lossless projected residual

For a fixed drift packet `Z`, retain the whole safe remainder rather than
one arc:

```text
S_(E,Z)=C_E minus union_(z in Z)D_z,
P_(E,Z)=pi_L(S_(E,Z)).
```

If the aligned multipliers are `A`, then

```text
the drift and aligned tails cover C_E
iff P_(E,Z) subset union_(a in A)D_a.
```

This quotient is lossless for the completion predicate.  In particular, a
`k`-aligned completion forces

```text
mu(P_(E,Z))<1-u_k.
```

The inequality is strict even when the two measures are equal: `P_(E,Z)` is
compact, the aligned union is open and proper, and a proper open subset of
the connected circle cannot equal a nonempty compact set.

On a body-safe cell `j`, use normalized coordinate `u` and put

```text
E_z(j)={u: ||z(j+u)/L||<1/14}.
```

Then De Morgan gives the finite exact formula

```text
T minus P_(E,Z)
 = intersection_(body-safe cells j) union_(z in Z)E_z(j).
```

Intersecting only a prefix of cells gives an upper bound for the final
common-danger set and hence a lower bound for `mu(P)`.  Once this lower
bound reaches `1-u_k`, the pair is dead.  Every endpoint and comparison is
rational; early stopping never promotes a lower bound to an equality.

## 6. Exact closure of five aligned tails and two drifts

For `k=5`, write

```text
eta_5=88/1365,       1-u_5=887/1365.
```

The inclusive component-discrepancy cap produces `626,787` first-drift
rows.  The non-strict endpoint matters: there are five integral cap-equality
rows, all of which fail the later high-excess predicate.

The `4,702` projected-suffix rows split into two banks.

1. On the `4,084` rows with `delta(z_1)>=eta_5 h`, delete `D_(z_1)` and
   apply the six-tail first-apex gate.  It closes `3,827` rows and leaves
   `257` rows, hence `42,912` exact integral `z_2` candidates.  The
   projected-prefix certificate kills every candidate; the longest prefix
   uses `22` body cells.
2. On the remaining `618` suffix rows, put
   `g=eta_5 h-delta(z_1)>0`.  The component discrepancy bound gives
   `z_2<=floor((6r/49)/g)`.  The resulting finite intervals contain
   `7,218,110` labels.  Exact singleton integration retains `194,073`
   pairs that pay `g`, and the projected-prefix certificate kills all of
   them.  The longest prefix uses `871` cells.

Both banks have the same closest rational margin:

```text
180/277 - 887/1365 = 1/378105.
```

It occurs after one cell at

```text
high: E=(1,2,3,5,9,10), z_1=24, z_2=277, L=1260;
sub:  E=(2,3,4,5,6,9),  z_1=22, z_2=554, L=2520.
```

Thus the literal five-aligned/two-drift branch is empty.  This projected
proof is a short lossless-quotient route independent of the divisor/address
mask atlas in THM-2928.

## 7. Four aligned tails: ordered typed-apex reduction

The exact `87,975`-row `k=4` projected suffix now has an independent ordered
physical reduction.  The direct non-strict integer-apex lemma

```text
A_p(R)=floor(6pr/[7(7-p)h])
```

is applied recursively after deleting `z_1` and then the least remaining
actual aligned multiplier on every complementary branch.  Strictly
increasing selected multipliers retain the forbidden-prefix sidecar and
make the recursion disjoint.  Its `5,368` states merge, on every suffix row,
to one necessary `z_2` prefix interval beginning at `z_1+1`.

The raw physical bank has `50,285,016` pairs and `z_2<=6515`.  Combining the
projected wall `z_3>=floor((2366/21875)L)+1` with the exact multiplicity
surplus splits it into

```text
41,770,842 scalar kills,
 2,042,669 finite-z_3 prefixes,
 6,471,505 no-positive-gap prefixes.
```

Thus only `8,514,174` necessary physical `(E,z_1,z_2)` prefixes survive,
uniformly with `z_2<=2163`.  This is not a realized-cover census.  See
`k4_typed_apex_proof_note.md` for the proof and
`canonical_promotion_proposal.md` for exact canonical wording.  If the
incoming Lorenz/activity/needle theorem closes the entire mixed `k=4`
sector, it supersedes this route only for emptiness; the typed-apex result
remains an independent physical-order finite reduction.
