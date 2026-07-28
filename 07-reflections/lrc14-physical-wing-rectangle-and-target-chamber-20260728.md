# The physical wing is a 7-by-7 rectangle feeding a unit-valued 9-by-9 chamber

> **STATUS: FINITE-EXACT + VERIFIED + INDEPENDENTLY HOSTILE-AUDITED.**
> On all 81 lawful labels of THM-2754's canonical rail-8 root-zero bank, the
> forgotten-clock wing has a simple physical-present normal form.  Its source
> is the monomial `12z^3` on one Cartesian `7 x 7` rectangle and is zero on
> the complementary cross.  Its target is a `Phi_7` unit on all 81 labels and
> factors into three one-dimensional label chambers.  The 49 unit-to-unit
> cells carry four exact multipliers rather than one uniform gain.  This is a
> coefficient-algebra chamber theorem, not a physical wing action, endpoint
> current, row exclusion, or proof of LRC(14).

## 1. The two clocks are retained before taking the wing

Use THM-2754's canonical rail `8`, common label universe

```text
S=(0,1,2,3,8,9,10,11,12),
T=(3,4,5,6,7,8,9,10,11).                              (1)
```

For physical-present clock `e`, let `N^-_e(s,t)` and `N^+_e(s,t)` be the
source and target one-sided carriers obtained by imposing the lawful
`E3`-restricted full-target section `F_(e,s,t)` at that endpoint.  Let
`C^-_e(s,t)` and `C^+_e(s,t)` be the corresponding two-sided common carriers,
with both `F(x)` and `F(x+7/13^6)` imposed before integration.  The weighted
wing objects are the literal step-function differences

```text
L_e=N^-_e minus C^-_e,             R_e=N^+_e minus C^+_e. (2)
```

The implementation writes the old diagonal `e=j`, where `j` is the delayed
coefficient clock.  That diagonal is lawful here only because three separate
facts are checked:

1. after `E3`, all `3969` relevant carrier objects are independent of `j`;
2. the terminal-Q prefix is empty at `j=0` and identical at `j=1,...,6`;
3. the physical `e=0` wing amplitude vanishes independently.

Thus fixing `j=1` recovers exactly the same physical-`e` row.  The independent
audit used that fixed-`j` construction rather than trusting the diagonal.

Divide every wing coefficient by the inherited content `26`, normalize the
source at private root `12` and the target at private root `1`, and quotient
the uniform clock line.  This gives profiles

```text
ell_(s,t), r_(s,t) in A=F_13[z]/(Phi_7(z)).             (3)
```

All `567` literal object subtractions in `(2)` agree with the differences of
the independently integrated coefficients.

## 2. The source is one Cartesian rectangle

Put

```text
S*=(0,2,8,9,10,11,12),
T*=(3,4,5,6,7,10,11).                                 (4)
```

The complete source law is

```text
ell_(s,t)=12 z^3 1_(s in S*) 1_(t in T*).              (5)
```

Hence exactly `49=7*7` source profiles are units and the other `32` vanish.
The zero set is not irregular: it is the union of the two vertical rows
`s=1,3` and the two horizontal columns `t=8,9`, with their four intersections
counted once,

```text
2*9+2*9-2*2=32.                                       (6)
```

This is the first exact object-level explanation of the rank drop seen at the
chosen label.  The apparent sporadic loss is a product-mask boundary.

## 3. The target is a three-tooth chamber and is always a unit

For `s in S`, define `(a_s,b_s,c_s)` by

| `s` | `(a_s,b_s,c_s)` |
|---|---|
| `0,12` | `(9,2,2)` |
| `1` | `(9,2,0)` |
| `2` | `(7,0,4)` |
| `3` | `(7,2,3)` |
| `8,9,10,11` | `(7,4,4)` |

Then every target profile is exactly

```text
r_(s,t)
 =a_s z
  +b_s 1_(t<=9) z^2
  +c_s 1_(t notin {8,9}) z^3.                         (7)
```

Thus the target depends on `t` only through two overlapping tooth masks:

```text
{3,...,9},                    {3,...,7,10,11}.         (8)
```

Direct multiplication determinants in `A` show that all `81/81` target rows
are units.  The determinant-pair census `(det ell,det r)` is

```text
(0,1):2, (0,3):11, (0,4):5, (0,8):8, (0,10):2, (0,12):4,
(1,3):4, (1,5):20, (1,8):15, (1,11):10.              (9)
```

The target therefore fills the whole `9 x 9` bank even where the source
coefficient functional kills a physically nonempty wing.

## 4. Four multipliers, not one global gain

On the `49` cells of `S* x T*`, both profiles are units.  The exact ratio
`g_(s,t)=r_(s,t)/ell_(s,t)` takes four values:

| symbol | multiplier in the power basis | count | determinant |
|---|---|---:|---:|
| `A` | `(0,4,4,4,4,10)` | 20 | 5 |
| `B` | `(9,0,0,0,0,6)` | 15 | 8 |
| `C` | `(0,2,2,2,2,6)` | 10 | 11 |
| `D` | `(11,0,0,0,0,4)` | 4 | 3 |

With columns ordered by `t=3,...,11`, the complete chamber is

```text
s= 0: C C C C C X X D D
s= 1: X X X X X X X X X
s= 2: B B B B B X X B B
s= 3: X X X X X X X X X
s= 8: A A A A A X X B B
s= 9: A A A A A X X B B
s=10: A A A A A X X B B
s=11: A A A A A X X B B
s=12: C C C C C X X D D,                              (10)
```

where `X` means source-null/target-unit.  The multiplier of determinant `11`
from THM-2754 is exactly chamber `C`; it is not uniform.

Augmentation sharpens the earlier hostile.  On the 49 source-unit cells, the
target/source scalar gains have census

```text
gain  0: 10,             gain 2: 19,             gain 11: 20. (11)
```

The displayed THM-2754 cell `(s,t)=(0,3)` lies in the ten-cell gain-zero
chamber.  The old provisional gain `2` does survive lawfully, but on nineteen
different labels and only after choosing the matching physical-clock universe.
There is still no label-independent augmented scalar gain.

## 5. No nontrivial label mode is lost

Zero-extend `(5)` and `(7)` from `S x T` to `F_13^2`.  For each nonzero
label frequency `(f,g)`, bucket the coefficient at the single residue
`fs+gt mod 13`.  A vanishing characteristic-zero `13`th-root coefficient
would force all thirteen rational bucket sums to be equal.  Their reductions
modulo `13` are nonconstant in every one of the `168` nontrivial modes.

More precisely, the source certificate is nonconstant only in physical-clock
coordinate `z^3`, while the target certificate is nonconstant in each of
`z,z^2,z^3`:

```text
source: (0,0,0,1,0,0) on all 168 modes,
target: (0,1,1,1,0,0) on all 168 modes.               (12)
```

This proves nonvanishing of every nontrivial normalized label mode.  It is
not a THM-2334 exact relation-address current: `(12)` retains the two target
labels and physical clock quotient, but still forgets Fourier allocation,
exact relation address, endpoint determinant sector, and a physical packet
map.

## 6. Structural reading and next test

The `81`-label wing is best viewed as

```text
one 7x7 source rectangle
  -> four unit-ratio subchambers
  -> one everywhere-unit 9x9 target chamber,                         (13)
```

not as one scalar shear.  This is a toothpick-like self-similarity at the
coefficient level: two one-dimensional forbidden strips generate the entire
source-null cross, while two one-dimensional target masks generate the target
profile.  The irregularity is confined to the chamber boundaries and then
propagates multiplicatively.

The cheapest next fixed-triangle transplant should therefore avoid the
augmentation-zero `C` chamber and retain one of:

- the twenty-cell `A` chamber, with scalar gain `11`; or
- the nineteen cells whose augmentation gain is `2`.

One must insert that chamber selector before collapsing THM-2334's factor
harmonics and retain one THM-2625 determinant sector.  Equations `(5)--(13)`
do not supply the missing source-to-target carrier operation, exact address,
semantic endpoint current, scalar-row exclusion, or LRC(14) decrement.

## 7. Reproduction and independent audit

Run

```bash
python3 04-computation/lrc14_root_zero_all81_physical_wing_chamber_20260728.py
python3 -O 04-computation/lrc14_root_zero_all81_physical_wing_chamber_20260728.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_root_zero_all81_physical_wing_chamber_20260728.out.
```

The companion uses exact integer interval arithmetic, explicit `require`
checks, literal weighted subtraction, finite-field quotient arithmetic, and
no truth-bearing `assert`.  The independent audit held delayed clock `j=1`
fixed, reconstructed the physical-clock family, and used a separate event-
sweep weighted subtraction.  It recovered all `567` object differences,
all `3969` delayed-clock collapse checks, formulas `(5)` and `(7)`, all unit
determinants, and the four ratios in `(10)`.

LF-normalized SHA-256 after promotion:

```text
script  f860a75450dd5549db85dd79c670ce21f2ad5035447ce1860f7065c7d4682f4a
output  9a596ef60282941a9c15b4571b68bc93a5b08597fe9011804c990c328e731798
```
