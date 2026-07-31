---
id: THM-2972
title: "All-four-high punctured-packet torsion closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The 39
  exact-suffix rows left by THM-2970 contain exactly 19,640 literal
  all-four-high packets meeting the projected k=2 scalar wall. Every packet
  has a selectable puncture with a complete-cell order-two torsion pair.
  Together with THM-2970 this empties all 58 rows on 1680<=z_1<=1742 and
  gives the intermediate projected k=2 cap z_1<=1679. The later canonical
  THM-2941 descent reaches 1656. This is not LRC(14).
source: codex-lrc14-all-four-high-punctured-packet-torsion-2026-07-30
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2970-located-low-cell-torsion-tail-closure
related:
  - MISTAKE-331
script: 04-computation/lrc14_j7_k2_all_four_high_punctured_packet_torsion_closure_thm2972.py
output: 05-knowledge/results/lrc14_j7_k2_all_four_high_punctured_packet_torsion_closure_thm2972.out
script_sha256: b92f7c6cf295f29a3ec2c64730f6bd4486073b9197f322202df98df497369d97
output_sha256: f6b34a08c30964f2c6ebf71b5e4f580fb6c75a4cf62762cb530ef892f9d4af25
hash_basis: LF-normalized bytes
---

# THM-2972 -- all-four-high punctured-packet torsion closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## Statement

Consider the thirty-nine `EXACT`-suffix survivor rows in the proved projected
`k=2` scalar atlas on

```text
1680<=z_1<=1742.                                           (1)
```

Thus a row fixes a six-speed body `E`, its period `L`, a first later label
`z_1`, and four further distinct later labels. Every one of the five later
labels is at or above the projected high floor.

There are exactly

```text
19,640                                                     (2)
```

literal five-label packets meeting the exact scalar necessary inequality in
these rows. For every such packet `P` there is a distinguished label `z in P`
such that the other four labels are safe on two complete body-safe `L`-cells
whose difference has order two modulo the exact denominator of `z`. The
located torsion lemma gives projected safe residual of full mass one.
Consequently all thirty-nine exact-suffix rows are empty.

THM-2970 independently empties the complementary nineteen `HIGH-TAIL` rows.
The two theorems therefore empty all fifty-eight atlas rows in `(1)`. Composed
with the proved prior cap `z_1<=1742`, this gives

```text
projected k=2 cap: z_1<=1679.                             (3)
```

The canonical `z_1=1736` exact-descent replay closes fifteen of the fifty-eight
rows and is an overlap control. Thus forty-three row closures are new relative
to that shell theorem. The earlier hybrid replay is superseded evidence under
MISTAKE-331/332. Equation `(3)` is a projected scalar-atlas descent, not a
decrement of the 165-profile LRC ledger and not a proof of LRC(14).

## Proof

### 1. The finite cutoff is exhaustive

Fix one row. On a positive unit ray with residue `r mod L`, the exact scalar
increment has the form

```text
delta_r(z)=a_r/(7 R z),          z=r mod L,   a_r>0,       (4)
```

where `R=14 lcm(1,...,14)`. It decreases strictly under `z -> z+L`.
Nonpositive rays cannot help the lower scalar inequality.

Let

```text
T=lower-delta(z_1)>0,                                    (5)
```

and let `B_1>=B_2>=B_3` be the three largest positive later-ray values.
They occur among the first four representatives of each positive residue:
a fifth representative on one ray is dominated by four earlier distinct
representatives on that same ray and cannot enter the global top three. Put

```text
eta=T-(B_1+B_2+B_3).                                    (6)
```

The exact sweep proves `eta>0` on all thirty-nine rows. Any label with
`delta(z)<eta` is inadmissible even after the other three slots are replaced
permissively by `B_1,B_2,B_3`. Enumerating a positive ray only while
`delta(z)>=eta` is therefore exhaustive.

This produces

```text
1,401 finite candidates in total,
at most 456 candidates on one row.                       (7)
```

Exact branch-and-bound enumeration of distinct four-subsets retains precisely the
packets in `(2)`. No label horizon, floating comparison, analytic
`HIGH-TAIL` placeholder, or solver-selected certificate enters this step.
The largest retained label is `11580<L=11760`; separately, the whole-cell
inequality makes the fixed-safe mask empty for any label `z>=L`.

### 2. Punctured-packet torsion

Let `C_E` be the complete `L`-cells contained in the body-safe carrier. For a
five-label packet `P` and `z in P`, define

```text
S(P\{z})={j in C_E:
           every label in P\{z} is safe on the whole j-cell}.           (8)
```

Write the exact denominator of `z mod L` as `d`, so uniquely

```text
z=(L/d)u+mL,                    gcd(u,d)=1.               (9)
```

Suppose `j,k in S(P\{z})` and for some divisor `q` of `d`,

```text
2<=q<=7,                        k-j=d/q mod d.            (10)
```

The two `z`-phase centres differ by

```text
z(k-j)/L=u/q mod 1.                                      (11)
```

Since `u` is a unit modulo `q`, their circular separation is at least
`1/q>=1/7`. Two strict-open danger arcs of radius `1/14` are disjoint. The
body and the other four labels are safe on both whole cells by `(8)`, so at
every projected phase one of the two lifts is safe from the entire packet.
The projected safe residual has full mass one, contradicting the inherited
completion cap.

For every packet in `(2)`, the verifier tries all five choices of `z`, folds
the bitmap `(8)` onto `Z/d`, and tests each order `2<=q<=7` dividing `d`.
At least one puncture works in every packet. Moreover, its deterministic
selection rule always finds an order-two witness. This is existential in the
puncture: it does not say every working puncture has order two.

### 3. Exact census and sharp puncture boundary

The row and packet height histograms are

```text
row counts:
  1694:7, 1702:3, 1708:9, 1722:7, 1724:2, 1732:2, 1736:9;

packet counts:
  1694:1585, 1702:199, 1708:12023, 1722:516,
  1724:1653, 1732:306, 1736:3358.                       (12)
```

All five punctures work on `19,630` packets. Exactly ten partial packets
remain, all at

```text
z_1=1708,             E=(1,4,8,10,12,14).               (13)
```

The smallest all-five hostile is

```text
P=(1708,1836,2060,2172,10236),                           (14)
```

for which only puncture index four works. Thus any claim that every puncture
works is false.

The independent audit also found `3,075` valid punctures with no order-two
edge. The first occurs at

```text
z_1=1694, E=(1,4,8,10,12,14),
P=(1694,1708,1836,2060,2128), puncture index 4,
d=105, u=19, working orders={3,5,7}.                     (15)
```

This disproves the stronger statement that every valid puncture has `q=2`
while preserving the deterministic existential conclusion used in the proof.

The frozen semantic digests are

```text
rows     c5bbb617c3f1cd15fce322751a0d4857adf731545a0408cd07fb7f5b97a7cf55
packets  774e915486b0284718808af4dfb8ef59ddfdd3a7080aea38b2e1ef92e2a90236
witness  4d8e5ba05dd1f8effb1e61b544f054df13be89a1aed730d402b6b1f2843fb7b2
partial  ef3d7d97893b939c95dc61a9d1fc94fa257f122b1748c0f92fb3054e0a083fe6
failures 01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b. (16)
```

### 4. Strict-open seam and first failed order

The endpoint order seven is legal only for strict-open danger. Take

```text
L=14, d=7, q=7, u=1, m=1, z=16, j=0, k=1, theta=13/16.
```

Both high arguments have circular norm exactly `1/14`; the closed arcs meet,
but both strict-open dangers omit the seam.

Order eight is the first universal failure. For

```text
L=56, d=8, q=8, u=1, m=1, z=63, j=0, k=1, theta=5/6,
```

the two arguments both have norm `1/16<1/14`, and the danger arcs overlap in
mass `1/7-1/8=1/56`. Thus neither closed-arc typing nor `q=8` may be imported
into the torsion lemma.

## Exact verification

Run

```text
python 04-computation/lrc14_j7_k2_all_four_high_punctured_packet_torsion_closure_thm2972.py --output .scratch/thm2972.normal.out
python -O 04-computation/lrc14_j7_k2_all_four_high_punctured_packet_torsion_closure_thm2972.py --output .scratch/thm2972.opt.out
```

Normal, optimized, and stored output are byte-identical. The canonical files
are already LF; their hashes are

```text
script  b92f7c6cf295f29a3ec2c64730f6bd4486073b9197f322202df98df497369d97
output  f6b34a08c30964f2c6ebf71b5e4f580fb6c75a4cf62762cb530ef892f9d4af25.
```

The independent audit used a separate pure-Python implementation: heap-merged
rays, recursive packet enumeration, and direct set-based cell inequalities.
It reproduced all counts, histograms, ten partial hostiles, the empty failure
set, and a `q=2` witness in every packet. Its audit-record digest is

```text
b43e05d9f754a082739b7a1d9f5234d39c274ab7c30dae3d04b965036ec66488.
```

This theorem closes only the pinned exact-suffix atlas plus the explicit
composition with THM-2970. It does not assert a uniform theorem below height
`1680`, alter the strict-open convention, or bypass the deterministic-evidence
repair recorded in MISTAKE-331.

**QED.**
