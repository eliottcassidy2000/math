---
id: THM-2980
title: "Projected k2 band 1600--1679: finite packets, nonpositive rays, and cyclic-carrier closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The exact projected k=2 scalar atlas on 1600<=z_1<=1679 has 68 rows.
  Exact punctured torsion closes 69,602 packets on 65 finite rows and
  1,098,353 zero-high packets on the three remaining rows. Their 6,934
  exactly-one-high projective classes split into 6,884 small-torsion
  classes and 50 full-projection classes; a unit-dependent high-order
  two-cell chord closes all 754 relevant unit instances. The omitted
  zero/negative-ray lane is repaired by 27 fixed positive-low triples: a sharp
  translated-band capacity closes all 3,861 carrier--denominator pairs, while
  5,159,799 unit instances independently regress the same first-residue
  universe, and every later height has projected-safe mass at least 3/7. The
  band is empty and the projected k=2 cap descends to z_1<=1599.
  This is not LRC(14) and removes no profile from the 165-row ledger.
source: codex-lrc14-k2-band-1600-1679-closure-2026-07-30
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2972-all-four-high-punctured-packet-torsion-closure
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
related:
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
  - LEM-012-gap-splitting-good-period-bound
  - MISTAKE-129
  - MISTAKE-334
  - MISTAKE-338
script:
  - 04-computation/lrc14_j7_k2_scalar_band_1600_1679_thm2980.py
  - 04-computation/lrc14_j7_k2_1600_1679_zero_high_and_high_order_phase_closure_thm2980.py
  - 04-computation/lrc14_j7_k2_1600_1679_translated_carrier_capacity_thm2980.py
  - 04-computation/lrc14_j7_k2_1600_1679_integrated_promotion_audit_thm2980.py
output:
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1600_1679_thm2980.out
  - 05-knowledge/results/lrc14_j7_k2_1600_1679_zero_high_and_high_order_phase_closure_thm2980.out
  - 05-knowledge/results/lrc14_j7_k2_1600_1679_translated_carrier_capacity_thm2980.out
  - 05-knowledge/results/lrc14_j7_k2_1600_1679_integrated_promotion_audit_thm2980.out
script_sha256:
  - ae6610e52bb18057e7c15f8599141c0a738419dc7d364d56b6688d652e466922
  - 7cc93dfe2e9365aa1e313d97ab19e76bce0943fd9f6a243e8c1b33c80d439ef4
  - 528d2b5c15e9fedaea37995201319244fd7caa27d3a2f3864a716e49c63cb44d
  - 56a074db49eb30a5579522d34bf21fd434680ae2212843d49bf959931ae45bcd
output_sha256:
  - 86967e15677757e71d4962f5b369afcd40fc95ea09ff1c5a8585e3c6dbe55964
  - 84253223ee1ef0042ad6e5585cd48f7d16fa201a7041e351c71cca8a76b593eb
  - b0c3a6fb818f82ff21ab1ef8788bfb6736059c81f57d55e2e01712ec2e2cdc7a
  - 433145c039ee083ae1784ca2e89f12e6179520367f92e0a18b2f5931dd769fa8
hash_basis: LF-normalized bytes
---

# THM-2980 -- projected k2 band 1600--1679: finite packets, nonpositive rays, and cyclic-carrier closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## Statement

The exact projected `k=2` scalar atlas on

```text
1600<=z_1<=1679                                            (1)
```

has 68 survivor rows: 40 exact-suffix rows and 28 analytic `HIGH-TAIL`
rows. Every admissible packet on every row is incompatible with the
complete-body safe carrier. Consequently the entire projected band `(1)` is
empty. Combined with THM-2972, which empties the contiguous band above it,
the projected cap becomes

```text
z_1<=1599.                                                (2)
```

The conclusion is a projected scalar-atlas descent. It is not LRC(14), does
not alter the 165-profile ledger, and does not assert a uniform theorem below
height 1600.

The proof also gives a reusable cyclic-carrier lemma. Let `S` be the residue
set left after puncturing one factor of exact denominator `d`, and let `u` be
that factor's unit. The puncture closes whenever `uS` is not contained in one
strict-open circular danger arc of length `1/7`. Equivalently, if `G(uS)` is
the largest cyclic gap in integer residue units, then

```text
7(d-G(uS))>=d.                                            (3)
```

The equality case closes because the danger arc is open. Small torsion of
order at most seven and full cyclic projection are two different sufficient
cases of `(3)`.

## Proof

### 1. Exact scalar inheritance

The first companion imports the pinned exact THM-2941 body engine and runs it
on every six-subset body and every first label in `(1)`. It checks the direct
engine hash before evaluation. The exact census is

```text
candidate rows before the scalar wall: 240,240
survivor rows after the scalar wall:         68
exact-suffix rows:                           40
HIGH-TAIL rows:                              28.           (4)
```

The profile and survivor semantic digests are respectively

```text
8ef2eead4eb68eeff33e38152c09648032e5652886594c71214c54bd55e39d83
b94748a736bc0069efb3b091b112527ad42a146b7eb549bd15f101658d2249d6. (5)
```

No floating comparison enters the atlas.

### 2. Cyclic-carrier diameter lemma

Fix a modulus `L`, punctured label

```text
z=(L/d)u+mL,                  gcd(u,d)=1,                 (6)
```

and a nonempty set `S` of residues modulo `d` represented by complete cells
on which the body and every unpunctured label are safe. On the cell over
`s in S`, the punctured label's phase differs from a common phase by `us/d`.
Thus one common phase can make the punctured factor dangerous on every cell
exactly when the finite set `uS` fits in one open circular arc of length
`1/7`.

Delete a largest cyclic gap of `uS`. The complementary lift is the shortest
closed containing arc, of length

```text
(d-G(uS))/d.                                              (7)
```

An open arc of length `1/7` contains this set iff `(7)` is strictly less
than `1/7`. This proves `(3)`, including its endpoint convention.

If `S` contains two residues differing by `d/q` for a divisor `2<=q<=7`
of `d`, multiplication by `u` leaves a nonzero unit modulo `q`; the two
images have circular separation at least `1/q>=1/7`. This recovers the
punctured small-torsion lemma.

If instead `S=Z/dZ`, then `G=1` and the full-projection arc margin is

```text
7(d-1)-d=6d-7.                                           (8)
```

There is also a two-cell witness. Put `t=ceil(d/7)` and choose residues

```text
0,  u^{-1}t mod d.                                       (9)
```

Their punctured phases differ by `t/d>=1/7`. The selected-pair margin is

```text
7 ceil(d/7)-d,                                           (10)
```

not the larger full-projection margin `(8)`. The distinction is important.
For the denominators in this theorem, the exact data are

| `d` | `t` | additive order of the chord | pair margin `(10)` | full margin `(8)` |
|---:|---:|---:|---:|---:|
| 11 | 2 | 11 | 3 | 59 |
| 13 | 2 | 13 | 1 | 71 |
| 143 | 21 | 143 | 4 | 851 |

The rule `(9)` works for every unit, but its second residue depends on `u`.
No fixed nonzero pair works for all units at `d=13`: the unit orbit of any
fixed displacement contains `1`, whose separation `1/13` is below `1/7`.

This is the precise holotopy content. The `q<=7` incidence bank is a lossy
bounded-order shadow of the cyclic carrier: it erases the high-order chords
of orders 11, 13, and 143. Once the full gain is retained, the obstruction is
already a missing edge between two lifted danger arcs, not a higher nerve
cycle. THM-2658 supplies the broader gain-nerve language but is not a
load-bearing dependency here.

THM-2984 now supplies the reusable ambient mechanism. Its primitive-unit gate
makes phase data height-free on one signed ray; its affine maximum-gap formula
specializes to `(3)`, and its translated complete-cell projection corollary is
load-bearing in Section 5. The present theorem supplies the exact row carriers
and shows that every one of them clears the sharp capacity. THM-2984 alone
empties no atlas row, while the present finite band proves no uniform signed-
ray theorem outside these carriers.

### 3. The 65 finite rows

For a row, let `T=lower-delta(z_1)`. On every positive residue ray, the exact
increment is `a_r/(7Rz)` and decreases strictly under `z -> z+L`. Let
`B_1>=B_2>=B_3` be the three largest positive later-ray values and put

```text
eta=T-(B_1+B_2+B_3).                                    (11)
```

On 65 rows, `eta>0`. Any fourth label below `eta` remains inadmissible even
after the other three slots are replaced by the permissive global maxima.
Enumerating each ray only while its value is at least `eta` is exhaustive.

This produces exactly

```text
40 exact rows:             16,809 packets
25 finite HIGH-TAIL rows:  52,793 packets
total:                     69,602 packets.              (12)
```

For each packet the verifier tries all five punctures in deterministic order,
folds the complete fixed-safe cell bitmap modulo the punctured denominator,
and tests every divisor `2<=q<=7`. Every packet has a witness. The frozen
packet and witness digests are

```text
packets  fead04db235721418dc93d36e7e1894a88a40fec36e584037ea089040e64e028
witness  e1ecfb3b254e93dc044d28e6baf399f74a3f7a7c2ea4a679c90d225a8a93abc0. (13)
```

### 4. Positive-ray compactification of the three large rows

The direct finite engine leaves out exactly

```text
(1612,(1,8,10,12,13,14)),
(1612,(1,10,11,12,13,14)),
(1650,(1,10,11,12,13,14)).                              (14)
```

On the first two rows `eta` in `(11)` is nonpositive. On the third it is the
positive rational `2291/1482506025`, but the exact finite candidate bank has
`112,954` entries, so the same projective representation is far smaller.
The word *large* is therefore implementation-neutral; it does not assert
that all three positive cutoffs have the same sign.

Call a suffix low if it lies below the row's exact high floor. Let `H` be the
largest first representative over every positive high ray and let
`L_1>=L_2` be the two largest low values. The verifier first checks the
load-bearing ordering

```text
H<=L_2.                                                   (15)
```

It then checks

```text
T-(four best lows)<0,
T-(2H+L_1+L_2)>0.                                        (16)
```

Because of `(15)`, the second inequality excludes packets with three or four
high entries as well as exactly two: each extra high can be permissively
replaced by one of `L_1,L_2`. Thus every packet has either zero or exactly
one high suffix.

The exact row census is

| `z_1`, body | `L` | lows | `zero_gap` | `two_gap` | zero-high packets | one-high classes / small-torsion failures |
|---|---:|---:|---:|---:|---:|---:|
| 1612, `(1,8,10,12,13,14)` | 152880 | 5840 | `-58783805/61219412436` | `51816343/183658237308` | 428600 | 4488 / 21 |
| 1612, `(1,10,11,12,13,14)` | 840840 | 35683 | `-5575571/5930024100` | `836456293/1297489273080` | 614180 | 2279 / 28 |
| 1650, `(1,10,11,12,13,14)` | 840840 | 35668 | `-197690873/274757783300` | `160311421/209272463400` | 55573 | 167 / 1 |

All `1,098,353` zero-high packets close by punctured torsion.

For an exactly-one-high packet, fix its denominator `d`. The three low labels
determine one projective class; every height and every unit ray with enough
remaining scalar mass is then checked. There are

```text
projective denominator/low-triple classes: 6,934
closed by q<=7 torsion:                    6,884
residual classes:                              50.       (17)
```

Every residual fixed-safe carrier projects onto all of `Z/dZ`. The split is

```text
d=11:   9 classes
d=13:  30 classes
d=143: 11 classes.                                      (18)
```

The verifier checks all 754 relevant units, forms the pair `(9)`, evaluates
its cyclic span directly, and closes every positive ray.

### 5. Zero and negative rays: the sign repair

The first candidate stopped at the preceding positive-ray partition. That
was false. On the first row,

```text
(1612,1836,2004,2340,20384)                              (19)
```

is scalar-admissible and the last scalar numerator is zero. On the second,

```text
(1612,1736,1800,2340,210210)                             (20)
```

is another zero-ray witness. More generally a negative ray `A/z` increases
toward zero with height, so a positive-ray upper cutoff has the wrong
monotonicity. MISTAKE-338 records the omitted universe and the repair.

The repair is finite after quotienting height. First, the exact best-two and
mixed-shape inequalities are

```text
row A: T-(L1+L2)     = 6587261/10803425724 >0,
       T-(H+L1+L2)   = 40949945/91829118654 >0;
row B: T-(L1+L2)     = 1664387/2372009640 >0,
       T-(H+L1+L2)   = 291145997/432496424360 >0;
row C: T-(H+L1+L2)   = 166276211/209272463400 >0.        (21)
```

Hence a packet containing a nonpositive suffix cannot contain two
nonpositive suffixes, nor one positive high plus only two positive lows. It
must have exactly three positive companions, all below the high floor. The
top-three gaps are negative on rows A and B and positive on C:

```text
-11173775/61219412436,
-142111/1186004820,
 2291/1482506025.                                        (22)
```

Thus rows A and B have exactly `18` and `9` admissible fixed positive-low
triples, while row C has none.

Fix one of these 27 triples and puncture the arbitrary nonpositive label.
The complete carrier safe from the body, `z_1`, and the triple is nonempty;
the exact minimum over all 27 carriers is `24,348` complete cells. For every
first nonzero residue representative `r<L`, the unit phase is height-free.
The full unit-by-unit regression census is

```text
                         row A        row B          total
fixed triples                18            9             27
nonpositive residue rays  76,445      420,421        496,866
unit instances         1,376,010    3,783,789      5,159,799
small-torsion closes   1,375,902    3,783,150      5,159,052
full-projection closes       108          639            747
failures                        0            0              0. (23)
```

The 108 first-row full projections all have `d=13`. The second-row 639
split as `d=11:45`, `d=13:54`, and `d=143:540`. An independent targeted
residue-fold implementation reproduced these counts and the two fixed-cell
minima before the full companion was frozen.

The conceptual proof is THM-2984's translation-robust complete-cell
projection corollary. For

```text
kappa(d)=ceil(d/7),
```

it closes every primitive unit and every common phase translation
as soon as the fixed carrier has more than `kappa(d)` residue classes modulo
`d`.  The third exact companion folds all `27` carriers over all of their
nonpositive denominators and checks
`18*119+9*191=3,861` carrier--denominator pairs.  Every pair passes, with
minimum slack `|S|-kappa(d)=1`, attained at `d=2`, `|S|=2`.  The coarser
complete-cell count implication passes only `1,926/2,142` and `1,601/1,719`
pairs, so the actual residue cardinality is essential on small denominators.
This `3,861`-case quotient is the load-bearing proof for every first
representative.  MISTAKE-334 explains why the smaller centered capacity
`beta(d)=2 floor((d-1)/14)+1` would be invalid here: the cell carrier retains
an arbitrary common phase.  Conversely, the frozen `5,159,799` unit instances are
an independent exhaustive regression of the stronger conceptual quotient;
they are not needed to infer universal phase escape.

It remains to justify all later representatives and residue zero without an
infinite height scan. Let `a=z/L>=1` and choose one fixed complete cell. In
the normalized cell coordinate the punctured phase traverses an interval of
length `a`, hence contains `floor(a)` full turns. Each full turn has safe
mass `6/7`, so the projected safe image has mass at least

```text
(6/7) floor(a)/a >= 3/7 > 25/91.                         (24)
```

The last fraction is the inherited `k=2` aligned completion cap. For `z=L`
the exact safe mass is `6/7`. Thus `(24)` closes residue zero and every later
height on every nonpositive ray. The finite, positive-projective, and
nonpositive-projective lanes now exhaust all packets on all 68 rows.

### 6. Sharp strict-open and small-order boundaries

At order seven, two phases may have norm exactly `1/14`; the closed danger
arcs touch, but the strict-open dangers are disjoint. At order eight there is
an exact hostile: two phases both have norm `1/16`, and the common danger
overlap has mass

```text
1/7-1/8=1/56.                                            (25)
```

Thus neither closed-arc typing nor a universal `q=8` extension is lawful.
The selected-pair ceiling in `(9)` is sharp for the same endpoint reason.

LEM-012/MISTAKE-129 is a useful dual control, not a dependency. At `d=143`
the full arithmetic progression has no large empty gap; the corrected
good-period invoice is unavailable. That same equidistribution is favorable
here because a full residue ruler cannot fit inside one danger arc.

## Exact verification and integer audit

Run

```text
python 04-computation/lrc14_j7_k2_scalar_band_1600_1679_thm2980.py --output .scratch/thm2980.scalar.normal.out
python -O 04-computation/lrc14_j7_k2_scalar_band_1600_1679_thm2980.py --output .scratch/thm2980.scalar.opt.out
python 04-computation/lrc14_j7_k2_1600_1679_zero_high_and_high_order_phase_closure_thm2980.py --output .scratch/thm2980.closure.normal.out
python -O 04-computation/lrc14_j7_k2_1600_1679_zero_high_and_high_order_phase_closure_thm2980.py --output .scratch/thm2980.closure.opt.out
python 04-computation/lrc14_j7_k2_1600_1679_translated_carrier_capacity_thm2980.py --output .scratch/thm2980.capacity.normal.out
python -O 04-computation/lrc14_j7_k2_1600_1679_translated_carrier_capacity_thm2980.py --output .scratch/thm2980.capacity.opt.out
```

Normal, optimized, and stored outputs are byte-identical for all three
companions. The translated-capacity transcript has `1,052` LF bytes. Its
frozen hashes are

```text
capacity script   528d2b5c15e9fedaea37995201319244fd7caa27d3a2f3864a716e49c63cb44d
capacity output   b0c3a6fb818f82ff21ab1ef8788bfb6736059c81f57d55e2e01712ec2e2cdc7a
```

All scalar comparisons use `Fraction`; all bulk arrays have integer dtype;
no floating computation is used. Before every NumPy product the verifier
checks the exact maximum primitive product, amplitude-times-modulus product,
label-times-cell product, and safe linear form against `2^63`. On each of the
68 real atlas rows it independently compares the vectorized full-cell safe
mask with scalar Boolean evaluation and compares block-OR residue folding
with direct scalar folding. Enumeration order is fixed for rows, candidates,
packets, punctures, orders, and units.

The final type maxima are

```text
max primitive input             3,787,559,275,500
max amplitude times modulus    79,009,564,233,600
max label times cell               54,707,079,630
max safe linear form                12,791,926
safe/fold full-cell checks             68 / 68
int64 limit                    9,223,372,036,854,775,807   (26)
```

The LF-normalized hashes are

```text
scalar script   ae6610e52bb18057e7c15f8599141c0a738419dc7d364d56b6688d652e466922
scalar output   86967e15677757e71d4962f5b369afcd40fc95ea09ff1c5a8585e3c6dbe55964
closure script  7cc93dfe2e9365aa1e313d97ab19e76bce0943fd9f6a243e8c1b33c80d439ef4
closure output  84253223ee1ef0042ad6e5585cd48f7d16fa201a7041e351c71cca8a76b593eb. (27)
```

## Independent hostile audit

The independent audit rederived the three exhaustion steps rather than using
the displayed packet totals as a premise.  On the 65 finite rows, the global
top-three bound proves that every omitted positive candidate is too small to
occur in an admissible four-label suffix.  On each exceptional row, `(21)`
rules out two nonpositive entries and also rules out one positive-high entry
beside a nonpositive entry, leaving exactly the 27 displayed positive-low
triples.  For a negative numerator the audit explicitly reverses the usual
ray monotonicity; this reproduces the two MISTAKE-338 zero-ray witnesses and
confirms that the later-height argument `(24)`, not a positive cutoff, is the
required compactification.

The cyclic geometry was checked independently at both sharp boundaries.  A
translated open interval of length `d/7` contains at most `ceil(d/7)` residue
points, so the strict inequality certified in all 3,861 carrier--denominator
pairs is sufficient for every primitive unit and every common phase.  At
equality, the order-seven endpoint is safe because danger is strict-open;
the order-eight control has genuine overlap `1/56`.  The audit also reproduced
the `d=28` translated-capacity hostile to the smaller centered threshold and
the absence of a fixed unit-independent two-cell chord at `d=13`.

Fresh scalar and translated-capacity replays are byte-identical to the stored
outputs.  The completed normal and optimized closure replays are mutually
and byte-for-byte identical to the stored 21,678-byte transcript.  All six
LF-normalized hashes equal `(27)` and the frontmatter, the live THM-2972
utility has the pinned hash `b92f7c6c...9d97`, and the integer and direct-
Boolean controls in `(26)` pass.  No scope or consequence beyond the projected
`k=2` band was used in the audit.

### Independent monolithic set-join cross-audit

A fourth verifier independently reruns one `3,003`-body by `80`-height atlas,
so its input universe has exactly `240,240` body-height rows.  It recovers the
same `68` survivors and height distribution, then parses only certified typed
terminal records from the twelve earlier component referees.  The exact join
is

```text
monolithic atlas rows       68
component terminal rows     68
duplicates                   0
missing                      0
extra                        0
upper/lower overlap          0.
```

Its terminal types are `40 EMPTY`, `22 SCALAR-EMPTY`, one
`LITERAL-CLOSURE-BELOW`, one `LITERAL-TORSION-EMPTY`, and four
`SECTION-CLOSURE-BELOW`; the terminal handoff cap is `1599`.  This is an
independent reconstruction of the atlas and atlas-to-component set join, not
an independent implementation of every component proof.  Ordinary and
optimized runs are byte-identical to the stored transcript; the fourth
source/output hashes are recorded in the frontmatter.

**QED.**
