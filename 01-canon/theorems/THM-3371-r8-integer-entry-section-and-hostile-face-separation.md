---
id: THM-3371
title: "AMM 12592: the R=8 integer entry section and hostile fractional-face separation"
status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED
source: kps-s178
depends_on:
  - THM-3329
  - THM-3332
related:
  - THM-3340
  - THM-3342
companion: 04-computation/amm12592_r8_integer_entry_section_kps_s178.py
output: 05-knowledge/results/amm12592_r8_integer_entry_section_kps_s178.out
audit: independent exact block/capacity/parity and active-lattice reconstruction
---

# THM-3371 -- R=8 integer entry section and hostile-face separation

## Statement

Use the halved junk coordinates `y=j/2` of THM-3329 at the exact golden-floor
epoch `R=8, D0=0`.  The degree word is

```text
(4,5,5,6,7,7,8,8),
```

and the displayed sparse polytope has `42` variables and `106` causal or
terminal inequalities.

### A. Exact integer point and one-row capture

The following row vector is an integer point of that polytope:

```text
y0 = ( 6,-3, 2, 1)
y1 = (-1, 0, 0, 0, 0)
y2 = (-1, 0, 0, 0, 0)
y3 = ( 0, 0, 0, 0, 0, 0)
y4 = ( 0, 0, 0, 0, 0, 0, 0)
y5 = ( 0, 0, 0, 0, 0, 0, 0)
y6 = ( 0, 0, 0, 0, 0, 0, 0, 0).
```

The first feed-free transition is `i_pf=3`, of ambient degree `d=6`.  The
original junk state entering that row is

```text
j=2y2=(-2,0,0,0,0).
```

It satisfies THM-3332's one-row hypotheses:

1. `F1`: its nonzero support is negative, `m=0`, and `m+2=2<6`;
2. `F2`: `a0=-j0=2<=d-1=5`;
3. `F3`: the only nontrivial check is at `t=2`, where
   `2a1+a0=2 <= 2*C(5,2)=20`, with margin `18`.
4. The finite-scale `F4` budget, which is not automatic at `R=8`, is
   `i_pf+ceil(a0/2)=3+1=4 <= R-2=6`.

Pascal transport into row `3` gives the load

```text
(-2,-4,-2,0,0,0,0),
```

which lies completely inside the degree-six clamp box.  Choosing that load as
the clamp leaves zero junk, so capture is immediate.  This is an exact
instance of the THM-3332 entry compiler, not merely terminal LP feasibility.

### B. The integer certificate leaves the denominator-103 active face

The integer point has `25` active displayed constraints of row rank `22`.
The denominator-`103` hostile point has `45` active constraints of rank `42`,
so it is a vertex.  Exactly `8` named constraints are active at both points,
and their rank is `7`.  Thus `37` hostile active constraints become strict at
the integer point, while `17` new constraints become active.

Consequently the integer certificate is not obtained while remaining in the
hostile vertex's minimal face.  The positive object selected by the proof is
the feed-free **entry section** and its one-row cone, not the active face of a
generic LP vertex.  No claim is made that the entry section itself is a face.

### C. Exact `R=512` prefix compiler

At `R=512,D0=0`, the last feed row is `129` and the first feed-free row is
`i_pf=130`; both adjacent degrees are `383`.  Hence a search aimed at a
first-feed-free THM-3332 certificate needs only the rows `y0,...,y129`:

```text
prefix variables = 44,750,
full sparse variables = 234,117,
prefix causal inequalities = 89,146.
```

In the halved integer coordinates, the sufficient entry cone is the single
linear system

```text
y_t <= 0                                      (0 <= t <= 382),
y_381 = y_382 = 0,
-2 y_0 <= 382,
-2 y_(t-1) - y_(t-2) <= C(382,t)             (2 <= t <= 382).
```

Therefore an integer-feasible causal prefix landing in this cone extends,
by THM-3332, to an exact-floor `R=512` witness.  This reduces a concrete
sufficient search to about `19.1%` of the full sparse variables.  It is not a
necessary condition for every possible `R=512` witness.

## Proof

The companion independently rebuilds the golden profile with exact
Fibonacci--Lucas comparisons, the initial `T4` row, the two feed channels,
every width-two/three Pascal transition, coordinate bounds, and the terminal
row.  Direct substitution proves feasibility of the displayed integer point.
Exact rational elimination gives the three active-set ranks and the
intersection count.  It then reconstructs original junk by `j=2y`, checks
`F1--F3`, and performs the next clamp cell by cell.  The `R=512` counts follow
from the same exact profile and the first row at which both feed channels are
absent.  THM-3332 supplies the tail extension from the stated cone.

Normal and optimized replays are byte-identical.  The LF-normalized source
hash is

```text
2c3e86a23e137d1349763e9cb631fde6ba01caaef1cc11e661ee3c53d31fc9a9
```

and the LF-normalized stored-output hash is

```text
f492e5cbf60a82a61be8d11f0c3bc72cf391d54b6ba1405f4f42f00fe7e5e0cd
```

An independent reconstruction additionally checked every correction block's
capacity and parity, the exact identity

```text
sum_(i=0)^7 x^i Delta_i(x)=(1-x)^7,
```

the full seven-cell capture load, the finite `F4` clock, both active sets,
their ranks, and the determinant `1648=16*103`.  It found no mathematical
failure after the two display omissions above were repaired.

## Boundary and non-consequences

- The denominator-`103` vertex still disproves total unimodularity and whole-
  polytope integrality.
- No `R=512` integer point is produced here.
- The prefix cone is sufficient, not necessary, for arbitrary exact-floor
  witnesses.
- No new AMM deadline bound, lower bound, or value of `C*` follows.
- The earlier proposal to find a circuit inside the hostile active face is
  superseded: every route to this integer point must leave that minimal face.
