# Typed T8--Prony fibre gate

Status: **FINITE-EXACT SCRATCH; NO CANON CLAIM OR ROW DISCHARGE.**

This companion tests whether the exact THM-2859 physical edge

```text
step 2 -> step 68 = +66h
```

can carry the THM-2863 `q=3,11` endpoint probes as one fully typed current.
It uses the twelve long-horn labels

```text
clock=1,  s in {0,3,12},  t in {5,6,9,10},
```

which collapse to one physical interval edge.  Its probe universe is the
`12 labels x 2 q-values = 24` edge rows, independently of multiplier.  The
four multipliers `Y=38,64,90,116` therefore see the same type gate.

The verdict is negative, but it localizes two distinct notions of fibre
product and reveals a stronger positive just before the failure.

## Pinned inheritance

LF-normalized SHA-256 pins:

```text
lrc14_horn_collar_endpoint_carry_thm2859.py
  6e062f3cc57c80fcff372c272bc138e280205bb953e484f1cc267340774260f0
lrc14_q3_q11_transverse_endpoint_horn_thm2847.py
  258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72
lrc14_literal_fixed_sheet_allocation_thm2806.py
  311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e
```

The promoted THM-2863 script/output are read as immutable Git blobs at
commit

```text
ab30502ddcd54fa9b4d597dd101184f2d3916b52
```

with hashes

```text
script  3d09e2f06cd17c38a34b00fa133f67963287173e5a69092e4e2d2b190f4459b9
output  2cf1902b19fd6ae9c6780521da0e431c8c7e1887235ea94ef4a8a5ca7cf9c7e4.
```

The calculation also replays THM-2859's exact base-root ancestry certificate

```text
counts=(966606,28534),
digest=15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd,
symmetric differences=(0,0).
```

## 1. Strict same-sheet/root fibre: `0/24`

Write `T=742586h`, so `T/13=57122h`.  Relative to the physical collar root,
the two q-sheets are displaced by

```text
q=3:   +171366h,
q=11:  -114244h       (centered modulo T).
```

Exact enumeration of every rooted path in the relevant labelled cell shows
that neither pulled q-interval is a common/right collar root at step 2 or
step 68.  All `48=12*2*2` state rows fail this first strict coordinate, hence
all 24 edge rows fail.  The `q=0` state is a positive control and is present
at both ends.

The pulled ancestry computation agrees with the root diagnosis:

```text
q=3:  (|U|,|V|)=(0,28535),
       digest=bdf95bb641f8d0ab0d8ebe764fecda782b9a7fe4995031c7f23e2f02868de931;
q=11: (|U|,|V|)=(0,28534),
       digest=89a69bdfdedb80eeabcf8eb3401012714a33cc2cb9690221f835c923d4cf9dd2.
```

These label sets are unchanged across `2 -> 68` and agree in the two
endpoint frames, but the `U` factor is empty.  This statement concerns the
pulled q-sheet.  It does not erase the nonempty, preserved ancestry on the
base root.

## 2. Weaker base-root current: a positive `24/24` pre-endpoint bank

There is a useful weaker construction: keep the physical step-2/68 root and
attach its q-shifted source and target currents.  On this interpretation,
every one of the 24 edge rows survives before endpoint transport.

At both edge ends, for both q-values and every label:

```text
root, native target, q-pulled source, q-pulled target
```

all retain the six factors

```text
(E3, clock, q1, q2, c2, c3)=111111,
```

and the shifted source and target carriers each contain the full atom with
weight

```text
27581135604.
```

The ancestry labels are constant across the edge.  Thus the relaxed failure
is not caused by `E3`, an ordinary factor, carrier support, or an ancestry
wall.  This `24/24` bank is the main new positive localization.

## 3. First relaxed loss: endpoint masks leave their translation orbits

The q0 endpoint masks give the hostile positive control: in both source and
target frames the unique translation is `(0,8)`.

For q3 and q11, by contrast, the masks at steps 2 and 68 are Cartesian
rectangles of equal mass but are not related by **any** translation of
`F_13^2`.  The source and target endpoint frames give identical results.

For q3 the rectangles change as

```text
A={0,1,2,3,4,5,6,7,8,9},
B={0,3,4,5,8,9,10,11,12}

                    ->

A={0,1,2,3,4,5,6,7,8,9},
B={0,5,6,7,8,9,10,11,12}.
```

Their cyclic second-axis gap necklaces are

```text
(1,1,1,1,1,3,1,1,3)
(1,1,1,1,1,1,1,1,5),
```

so translation equivalence is impossible.  Applying the physical `(0,8)`
map gives Hamming defect `80`, split as `40` lost and `40` gained addresses.

For q11 the rectangles change as

```text
A={0,1,2,3,4,5,8,9,12},
B={0,1,2,3,4,5,8,11,12}

                    ->

A={0,1,2,3,4,5,8,9,12},
B={0,1,2,5,6,7,8,11,12}.
```

The second-axis gap necklaces are

```text
(1,1,1,1,1,1,1,3,3)
(1,1,1,1,3,1,1,1,3).
```

Again there is no translation.  The `(0,8)` Hamming defect is `72`, split
as `36` lost and `36` gained addresses.

This is a structural obstruction, not a numerical mismatch: q3 changes the
gap multiset itself, while q11 preserves the gap multiset but changes its
cyclic ordering.

## 4. Why the scalar probe shadow looks healthy

The two named THM-2863 origins see only a lossy two-column restriction.  In
the right endpoint frame their occupancies are unchanged across the edge:

```text
q3:  (1,0) -> (1,0),
q11: (1,1) -> (1,1).
```

Consequently a two-origin or four-sample scalar calculation can remain
nonzero while the complete endpoint-address object has already left its
translation orbit.  Following the requested validity gate, the companion
does not evaluate endpoint sums, split Prony nodes, or form
`Pi_3(delta_8-delta_0)` after this failure.

For every `Y=38,64,90,116`, the census is therefore

```text
strict typed edges          0/24
relaxed pre-endpoint edges 24/24
relaxed fully typed edges   0/24.
```

## 5. Alignment with reserved THM-2868

The reserved stub

```text
THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
```

was read from origin commit
`3551e46666ae38708b02874fc49f0349f65c3846` (stub SHA-256
`b5830672d4ee0c0ecee4a7ce1c969bd3fa9853de73680bc652e68954eb6a16b6`).
It proposes a signed two-origin/two-multiplier projective projector on the
common 42-cell horn and explicitly does not claim one positive physical
packet or the q7 lift.

This audit does not duplicate or refute that target.  The unchanged
two-origin occupancy shadow leaves room for its signed coefficient-level
projector.  The immediate obstruction is narrower:

> A successful THM-2868 projective ratio cannot, by itself, serve as the
> missing physical vertical arrow from the THM-2859 `T8` edge, because on
> the exact twelve long labels the full q3/q11 endpoint masks admit no
> translation transport at all.

A future composition would need an additional rechart acting on the full
endpoint rectangle, with its action on signed coefficients proved, rather
than inferred from two origin columns.  It would also still need the q7
`E3`/complement lift.

## Reproduction

```bash
PYTHONHASHSEED=0 python3 \
  04-computation/lrc14_horn_collar_prony_typed_descent_gate_thm2859.py
PYTHONHASHSEED=0 python3 -O \
  04-computation/lrc14_horn_collar_prony_typed_descent_gate_thm2859.py
```

The integrated companion checks its own AST for executable assertion
statements.  An independent normal and optimized replay produced byte-identical
stdout, equal to the stored output.  LF-normalized SHA-256:

```text
script  c1563067fe0a6011fd512808fd377f1c1baeb64ad6cee598e277e9ed3b05fd87
output  5352e4b30b8822df4ffcabc5a35125eb9fc91c53444c7b68a8c1da6baafbfa20
```
