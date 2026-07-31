# LRC14 k=3 physical denominator bridge

**Status:** `FINITE-EXACT / SCRATCH / NECESSARY RELAXATION`.  This package does
not close the k=3 sector and does not realize a denominator profile by later
drift labels.

## Result

The bridge reconstructs the exact old THM-2941 projected physical ledger of
`376,020` pairs `(E,z1)` and assigns the distinguished reduced denominator

```text
L=14 lcm(E),                 d1=L/gcd(L,z1).
```

After separately deleting the canonically closed `z1=380` row, the special
`(E,z1)=((1,4,8,10,12,14),250)` row, and the nine canonically closed `z1=378`
rows, the live inherited ledger has `376,009` rows.  Every one admits at least
one support-transfer resolving denominator and a conditional denominator
profile.  The expected-spike screen leaves `247,615`; intersecting it with the
row-specific support-status screen leaves `247,614`.  Thus the bridge removes

```text
376,009 - 247,614 = 128,395
```

physical rows, or about `34.15%`, without enumerating later numerators.  At the
profile-occurrence level the same intersection reduces

```text
75,429,056,976  ->  18,778,875,909.
```

Support status alone leaves at least one profile on every physical row, though
it removes more than thirteen billion profile occurrences.  Its independent
contribution after expected spike is exactly one physical row:

```text
E=(1,2,3,5,6,10),  z1=105,  L=420,  d1=4.
```

Only `D=420` contributes to either surviving screen.  There
`|S_D|=198`, the `q=60` target counts are

```text
(N_0,N_1,N_2,N_3,N_4)=(60,60,48,14,10),
```

and the exact small-denominator top loads are `(99,66,51)`.  Expected spike
retains `262` profiles, all with `c=3`, but every one fails support status.
Support status retains `468` profiles (`276` with `c=1`, `192` with `c=2`),
but all fail expected spike.  Thus the two necessary screens have disjoint
conditional-profile support on this row.  This is a useful structural warning:
rowwise survival under separately chosen denominator profiles cannot be
intersected; the distinguished-denominator feature has to be retained.

The current post-closure physical ledger and the joint survivor ledger both
have maximum `z1=364`.  This is **not** a new cap improvement: after the nine
`z1=378` closures, the inherited ledger itself has no row with
`365<=z1<=377`.

## Conditional exact-lcm formula

Fix a physical row and a candidate resolving denominator `D`, put `q=D/7`,
and require `d1|D`.  For `E|D`, define

```text
U_D(E)=#{d|E : d>1 and d does not divide q},
V_D(E)=tau(E)-1-U_D(E),
epsilon_1=1_(d1 does not divide q),
Mult(a,j)=binom(a+j-1,j),  Mult(a,0)=1.
```

Then the number of later denominator three-multisets whose union with the
distinguished `d1` has exact lcm `D` and total uniform-one count `c` is

```text
A_(D,d1,c)
 = sum_(E|D, d1|E) mu(D/E)
     Mult(U_D(E),c-epsilon_1)
     Mult(V_D(E),3-c+epsilon_1).                 (1)
```

The implementation refines `(1)` by the additive feature

```text
(multiplicity of d=2,3,4; total d>4 ceiling capacity; c).
```

This is enough to apply expected spike and support status to the **same**
conditional profile.  A symbol-by-symbol current-lcm recurrence independently
reconstructs the complete refined coefficient table in `110` controls; a
literal three-multiset enumerator checks every `(D,d1)` with septimal
`D<=140` (`109` cases).

## Why expected spike is valid

For `q=D/7`, every denominator `d|D` has one of two traces.

- If `d` does not divide `q`, its enlarged mask contributes at most one point
  to every `q`-fibre.  These are the `c` uniform-one masks.
- If `d|q`, its literal common-phase mask fills `Y_d(u)` complete `q`-fibres,
  and exact phase averaging gives

  ```text
  integral_T Y_d(u) du = q/7.
  ```

Put `m=4-c` and let `N_c` be the number of target fibres whose projected body
support load exceeds `c`.  Coverage forces the spike masks to hit at least
`N_c` fibres throughout the compact three-aligned safe carrier.  Markov and
the compact/proper-open equality exclusion therefore give the strict necessary
condition

```text
N_c=0  or  55 N_c < 13 m q.                       (2)
```

The zero-demand clause in `(2)` is essential when `m=0`.

## Support-status intersection

For each row and `D`, the script computes the exact largest projected-support
residue-class load for `d=2,3,4`.  A small denominator contributes that load
when its open activity bit is on; its exact marginal is `d/7`.  Denominators
larger than four receive full ceiling capacity.  THM-2928's upward-event
fractional-cover theorem then gives the exact maximum one-threshold allowance
for this real marginal relaxation, with equality excluded by compactness and
proper openness.  Both `(2)` and this status screen are necessary for an actual
profile, so intersecting them coefficientwise is sound.  Their common-phase
bits are not optimized jointly, so the intersection is still an upper
relaxation.

## What is forgotten

The bridge preserves the physical first label and its denominator, but forgets
the later numerators, their strict order and distinctness, literal phase
locations, compatibility of one status table across all load thresholds, and
the projected high-wall condition beyond the already inherited physical row.
Consequently a surviving profile is not a surviving speed packet.

The most promising next refinement is the exact two-level floor/exception law
of each spike denominator.  The mean in `(2)` forgets whether the spike count
is the floor or ceiling of `d/7`; retaining these optional types should be
especially effective on the one-spike (`c=3`) residual.  Later labels and the
fixed first numerator should be imposed only after this scalable quotient is
exhausted.

## Reproduction

```bash
python3 .scratch/lrc14_k3_physical_denominator_bridge_20260730/bridge.py --workers 8
python3 -O .scratch/lrc14_k3_physical_denominator_bridge_20260730/bridge.py --workers 8
python3 .scratch/lrc14_k3_physical_denominator_bridge_20260730/independent_audit.py
python3 -O .scratch/lrc14_k3_physical_denominator_bridge_20260730/independent_audit.py
```

Ordinary and optimized transcripts must be byte-identical.  The bridge pins
the canonical THM-2941 physical source and THM-2928 body-projection source by
SHA-256, so a future canonical septimal promotion can retarget those two paths
without changing the conditional-GF interface.
