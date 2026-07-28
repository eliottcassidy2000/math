# The odometer is a skew product, not a tail controller

**Status: VERIFIED-EXACT FINITE SCOUT + STRUCTURAL REDUCTION.**  This note is
not a row exclusion and not an LRC(14) proof.  It records what the affine
integer lift can and cannot change, and it gives exact positive and hostile
controls for one four-letter near-carry alphabet.

The closest proved mechanisms are THM-2657's nonsplit `C_(13^6)` odometer,
THM-2682/2684's three-event failure for the scalar `D(x)={13x}` carrier, and
THM-2689's support/clock tradeoff for affine handoffs.  The canonical near
miss is the attractive alternating two-cycle: it has the right intrinsic
clock interfaces, but symmetry forces speed-14 resonance.  The least-used
sidecar is the terminal base-13 tail itself.

## 1. The exact state has a base and a fibre

Let

```text
R=13^6,                    T_k(x)={13x+k/R},
Rx=N+y,                    N in Z/RZ,  0<=y<1.
```

Writing `a=floor(13y)` gives the exact recurrence

```text
y'={13y},                  N'=13N+k+a  (mod R).          (1)
```

Thus the integer handoff label `k` never steers the terminal tail.  It acts
only in the odometer fibre.  For a tail of period `n`, the fibre-closing
equation is

```text
(13^n-1)N_0 + sum_i 13^(n-1-i)(k_i+a_i) = 0  (mod R).  (2)
```

Because `gcd(13^n-1,13^6)=1`, (2) has exactly one solution `N_0` for every
periodic tail and every carry word.  This is the useful reversal of viewpoint:
first find a recurrent allowed tail; only then solve its uniquely determined
odometer lift.  Searching carry words before checking tail recurrence spends
most of the computation decorating a base orbit that the delayed language
will later reject.

## 2. The symmetric cycle has an arithmetic resonance

Suppose a two-cycle is symmetric about `1/2`:

```text
x_+=1/2+u,   x_-=1/2-u,
x_+ --T_(-K)--> x_- --T_K--> x_+.
```

Then `(13+1)u=K/R`.  Hence

```text
{Rx_+}={1/2+K/14},         {Rx_-}={1/2-K/14},
```

and both tails lie exactly on the speed-14 resonance.  The failure is not a
bad choice of `K`; it is forced by the cyclotomic factor `13+1`.  Asymmetric
periods, rather than larger symmetric carries, are the first honest escape.

## 3. Exact asymmetric escape before the delayed word

The companion script exhausts words of lengths `3,...,7` in

```text
{-13^5-2,-13^5-1,13^5+1,13^5+2}.
```

There are no intrinsic-clock cycles of odd length.  There are `32` at length
four and `128` at length six, including respectively `24` and `120` of exact
period.  Apart from four lower-period presentations in each even census, the
cycles are non-speed-14-resonant.

One exact period-four witness uses

```text
(-371295,371294,-371295,371295).
```

Its intrinsic edges alternate `4->3,3->4`; its carries are `7,5,7,5`; its
future digits are `7,6,7,5`; and its root steps are `9,2,9,4`.  A common
source-one choice gives rail digits `0,12,0,12`, with every rail and present
inequality strict.  So variable affine handoffs genuinely evade scalar
three-event clock nilpotence at the rail/present level.

The hostile is equally sharp.  Every exact period-four and period-six state
in this alphabet fails the current delayed word already at its target-a
danger factor `D_(13^3)`.  The minimum positive failure margins are
`1433/4080` and `847565/2413404`.  The affine lift can repair the fibre, but
it cannot repair the tail language in (1).

## 4. The next proof object

The carrier should be treated as a finite-state skew product:

```text
allowed delayed-tail subshift  ----> recurrent SCC / periodic base orbit
             |                                      |
             v                                      v
unique odometer fibre lift     ----> carry, root, unit, endpoint sidecars.
```

The source is the delayed guard language; the target is a recurrent component
with a unique `Z/13^6` lift; the map is (1).  It preserves exact chronology and
all terminal digits.  Projecting to the tail destroys carry/root labels, while
projecting to the fibre destroys the delayed-word obstruction.  Both
coordinates are therefore necessary, but they should be solved in that order.

The companion depth-four work sharpens the hostile further: for the current
full delayed carrier, the base graph has no recurrent SCC at all.  That result
is uniform in every lift label, whereas the present scout is deliberately only
a finite four-letter positive control.  Consequently the next productive
search is not a larger carry alphabet on the same delayed language.  It is a
different delayed carrier or a lawful phase-local handoff whose base subshift
actually has recurrence; only then should the odometer, unit, and endpoint
fibres be attached.

## 5. Reproduction and scope

```bash
python3 04-computation/lrc14_affine_odometer_periodic_tail_scout_20260728.py
python3 -O 04-computation/lrc14_affine_odometer_periodic_tail_scout_20260728.py
```

Both runs byte-match
`05-knowledge/results/lrc14_affine_odometer_periodic_tail_scout_20260728.out`.
The script/output SHA256 values are respectively
`63550fcb0724888ad53c38a39cd33178cdea4d36dc4ee4daa00c72fde4064854` and
`a5f5226e36f09b1f62c7a4b5cbe7aa3623c2b7e5f4b158d9692ce071437c2903`.
No full packet transition, unit transport, endpoint gluing, row exclusion, or
LRC(14) conclusion is asserted.
