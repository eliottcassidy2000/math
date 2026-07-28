# The source-clocked natural root-zero sheet has a uniform endpoint gain but is not the common clutch

> **STATUS: FINITE-EXACT + VERIFIED + INDEPENDENTLY HOSTILE-AUDITED.**
> After the MISTAKE-313 repair, the natural sheet is built with the literal
> `source_present_section(E3, source_clock=1, s, t, ...)` constructor in each
> physical chart.  On all 81 labels lawful at the fixed THM-2744 endpoint, the
> source and target remain private units, every rowwise primitive-`t` and
> mixed `f,g!=0` target character survives, and the normalized target/source gain is uniformly
> `11`, not the retracted clock-blind value `7`.  The proved two-sided
> THM-2749 common section has equal raw vectors and mirror gain `-1`; the
> difference is a genuine one-sided right wing.  This is a coefficient-side
> physical-chart naturality boundary, not an endpoint current or LRC(14).

## 1. The MISTAKE-313 factor audit

The legacy helper `restrict_e3_and_sheet` inserted `E3` and the four shifted
`q1/c2/q2/c3` safe factors but omitted the source-one clock comb.  Its exact
gain `7` is therefore a clock-blind hostile, not the fully marked natural
sheet.  The repaired computation calls the proved THM-2742 constructor
directly:

```text
source_present_section(module, E3, source_clock=1, s, t, clock_combs).
                                                                  (1)
```

Thus every section includes

```text
E3,
d_c1(1/7),
g_q1(-s/13) g_c2(s/13),
g_q2(-t/13) g_c3(t/13).                                (2)
```

The section `(2)` is intersected independently with the source and translated
target rail-8/relative-present/seam carriers, in their own chart coordinates,
before delayed integration against `Q_(3,{1,2})`.  It is a genuine physical
source-clock-one subcarrier on both sides.  It is neither clock-blind nor the
two-sided fibre product `A intersect T_tau^(-1)A` of THM-2749.  The fixed
owner/source clock in `(1)` is also distinct from the seven delayed
coefficient-clock entries: clock zero has empty terminal prefix, while clocks
one through six have one common prefix.

The script checks that, after `(2)`, the physical interval list is literally
independent of that delayed coefficient clock.  Hence every computed vector
has the exact form

```text
(0,a,a,a,a,a,a).                                       (3)
```

This reduction is proved by tuple equality before the single delayed-prefix
evaluation; it is not inferred from sampled coefficients.

## 2. Exact support laws

Put

```text
S_clock={0,1,2,3,6,7,8,9,10,11,12}.                   (4)
```

The complete clocked `13 x 13` bank satisfies

```text
source support = S_clock x {3,...,11},                 count 99,
target support = S_clock x {3,...,12},                 count 110,
common support = S_clock x {3,...,11},                 count 99. (5)
```

The source-clock factor kills exactly `s=4,5`; this is the first visible
change from the retracted clock-blind census.  The fixed THM-2744 endpoint has
common lawful label set

```text
S_end={0,1,2,3,8,9,10,11,12},
T_end={3,4,5,6,7,8,9,10,11}.                           (6)
```

All `9*9=81` labels in `(6)` survive on both sides.  The extra global
coefficient sections `s=6,7` in `(5)` are not silently typed as labels at the
fixed endpoint.

## 3. Exact gain and wing anatomy

At `(s,t)=(0,3)`, the independently recomputed vectors are

```text
A=(0,C,C,C,C,C,C),     C=339633525654239542165440,
B=(0,D,D,D,D,D,D),     D=345341652135823400016960.     (7)
```

After content-26 division and root normalization at source root `12` and
target root `1`, their scalar profiles are `9` and `8`.  Therefore

```text
target/source = 8/9 = 11 mod13.                         (8)
```

The complete common-support gain histogram is

```text
gain 11: 81 labels,       gain 6: 18 labels.            (9)
```

The `81` gain-`11` labels are exactly the endpoint-typed window `(6)`; the
gain-`6` labels are the two extra `s=6,7` rows.  Thus the corrected natural
gain is uniform on the canonical endpoint window but not on the broader
coefficient bank.  The former value `7` is not a special subcase: it belonged
to the omitted-clock carrier and is fully retracted by MISTAKE-313.

THM-2749's two-sided common section has the same raw source coefficient `C`
on the endpoint window.  Hence its left wing is coefficient-null.  On the
target side, subtracting the common vector from `(7)` leaves

```text
D-C=5708126481583857851520,                             (10)
```

the nonzero right wing.  Root-normalized target profiles split as

```text
natural target 8 = common target 4 + right wing 4.      (11)
```

This gives the exact coefficient-vector anatomy:

```text
natural source coefficient vector = common source coefficient vector,
natural target = common target + right wing.            (12)
```

It explains why the natural sheet is physical and nonzero without being an
intertwiner: the one-sided target has an additional chart-supported wing.

## 4. Every rowwise primitive-`t` and mixed character survives

For each fixed active `s`, extend the source amplitude `a_(s,t)`, target
amplitude `b_(s,t)`, and cross amplitude
`c_(s,t)=a_(s,t)b_(s,t)` by zero to all of `F13`.  Exact reduction modulo
`Phi_13` gives every primitive `t`-character on every active row:

```text
full clocked bank: source 132/156, target 132/156, cross 132/156,
                    i.e. 12/12 on each of 11 active rows;
endpoint window:    source 108/108, target 108/108, cross 108/108. (13)
```

The two-coordinate character uses one common cyclotomic field:

```text
chi_(f,g)(s,t)=zeta_13^(f s+g t).                       (14)
```

Accordingly the companion buckets `c_(s,t)` by the single residue
`fs+gt mod13`.  A coefficient vanishes exactly when all thirteen rational
bucket sums agree.  All `144/144` mixed characters survive on the full bank
and again on the endpoint window.

The discarded first draft reduced in `s` and `t` separately, thereby working
in `Q(zeta_13) tensor Q(zeta_13)` rather than the character field of
`C13 x C13`.  Its count happens also to be `144/144`, both before and after
endpoint restriction.  This numerical coincidence is retained as a hostile:
matching a correct census does not repair a mistyped representation.

The raw cross coefficients

```text
K_(f,g)=sum_(s,t) a_(s,t)b_(s,t) zeta_13^(f s+g t)     (15)
```

are therefore genuinely nonzero in every primitive rowwise and mixed mode.
They pair two independently integrated coefficient tables; they do not
identify their physical interval carriers or turn either into the global
THM-2334 endpoint current.

## 5. Comparison and honest boundary

The exact comparison is

```text
clock-blind legacy sheet     -> gain 7 (retracted physical typing),
clocked natural sheet        -> gain 11 on all 81 endpoint labels,
two-sided common section     -> equal raw vectors, mirror gain -1. (16)
```

Thus retaining the missing source clock repairs the old typing but does not
make the natural sheet equal to the common clutch.  Character completeness
is strictly weaker than naturality, and naturality is strictly weaker than a
global endpoint current.  No row is excluded; the LRC(14) ledger remains 165.

## 6. Exact reproduction

Run

```bash
python 04-computation/lrc14_natural_single_sheet_root_zero_target_bank_20260728.py
python -O 04-computation/lrc14_natural_single_sheet_root_zero_target_bank_20260728.py
```

against
`05-knowledge/results/lrc14_natural_single_sheet_root_zero_target_bank_20260728.out`.
Both modes byte-match the stored transcript.  The script uses exact interval,
integer, and cyclotomic arithmetic, explicit `require` checks, and no
truth-bearing Python `assert`.

LF-normalized SHA-256:

```text
script  591241e912c452b1985c7fc700884183a6e50440a579e1d853123d276cd416a8
output  42d0013e15de975070a24e805dabfcd2b125dd1b1eb857931a632b51b6ba6b70
```
