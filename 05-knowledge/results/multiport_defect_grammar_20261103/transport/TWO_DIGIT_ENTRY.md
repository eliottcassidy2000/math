# Both uniform digits from one actual infinite entry family

**PROVED construction family**, with literal checks. This supplies both
`3H+2` and `3H+4` from the same infinitely many named sources. It does **not**
prove that either output can repeat both constructions, or establish global
spectrum completeness. The sources are the inherited explicit family below,
not replacements of the170807 lineage by equal-count graphs.

Use the inherited seven-coordinate convention and ordered-sum law from
[`session182_physical_returns_20261104/induction/PROOF.md`](../../session182_physical_returns_20261104/induction/PROOF.md).
For `M(TT_t;B,W)`, t≥1, its exact count and outer-root defect are

```
H_t=(h-j)+2^t f,             b_r=j,
```

and the full tournament is strong exactly when W meets the terminal SCC of B.
Physical names and the actual port set remain part of the input.

## The three actual suffixes

Start with `B0=TT3`, names `0→1→2`, and `W0={0,1,2}`. For digit2, replace
the terminal singleton component2 by the triangle `2→3→4→2`, with only2
marked inside that triangle. The preceding0→1 remains fully marked. For
digit4 use the inherited `TT5` suffix `0→1→2→3→4`, with ports `{0,2,4}`.

| Suffix | Ports | Exact signature (h,p,f,g,j,e,s) |
|---|---|---|
| B0=TT3 |0,1,2|(1,8,7,7,6,1,1)|
| B2=TT2→C3 |0,1,2|(3,48,21,39,16,1,3)|
| B4=TT5 |0,2,4|(1,32,21,21,12,1,1)|

The B2 endpoint coordinate is **e=1**, not3: the terminal C3 has exactly one
Hamilton path ending at its single marked vertex2. This correction does not
change f, j or the return identity.

All terminal SCCs meet their actual ports. Therefore, for every t≥1,

```
H0(t) = -5+7·2^t,
H2(t) = -13+21·2^t = 3H0(t)+2,
H4(t) = -11+21·2^t = 3H0(t)+4.
```

The first four triples are `(9,29,31)`, `(23,71,73)`, `(51,155,157)` and
`(107,323,325)`. Root defects are transported as `6→16` and `6→12`.
This is the positive sharpness example for the purely transitive-suffix
obstruction: the digit2 output deliberately leaves that class.

## Literal name-preserving constructions

Keep the actual transitive prefix A, root r and old suffix names0,1,2.
For digit2:

1. Add4 with old in-neighbors `A∪{0,1,r}`. Thus4's sole old out-neighbor is2.
2. Add3 with in-neighbors `A∪{0,1,2,r}`, and with `3→4`.

Every old edge is retained. After the first birth the suffix order is
`0→1→4→2`, whose terminal vertex2 is still a port. After the second birth
the terminal component is `C3[2,3,4]`, still containing port2. Both completed
states are strong. The forest carrier can implement these two births with
the prefix anonymous and the displayed old names active; both newborns lose
uniformly to the prefix.

These births are the converse of the inherited adaptive two-birth reset.
The converse of [PROOF.md, T12](PROOF.md) gives `H2=H0+2S_2(2)` and a
start-zero new port. Comparing the displayed all-t count identity proves
that source port2 satisfies `S_2(2)-H0=1` for every t, and the digit2
output has `S_4(2)-H2=0`. This start-zero gauge at physical4 is different from
the outer-root end-defect16 at physical r; neither is silently reset.
Repeating this reset at port4 therefore gives `3H2`, not `3H2+2`: the macro
does not regenerate its own digit2 guard.

For digit4, the inherited all-strong schedule is: add4 as a terminal suffix
port; add nonport3 between old2 and4; reverse the old edge `1→r`. Only that
old edge changes. Terminal port4 certifies strongness throughout. These are
direct adjacency constructions; no Hamilton-path-complement guard certificate
is claimed for the reversal.

The two outputs also have a literal three-reversal sibling link. Starting
at B2's full attachment, reverse in order

```
r→4,       1→r,       4→2.
```

The first adds4 to the in-port set, the second removes1, and the third turns
the terminal triangle into the transitive order2→3→4. Every completed graph
is strong. The result is exactly the B4 graph on the same physical names,
not just an equal-count graph. The endpoint suffixes both have f=21;
H increases2, `E_r(2)` decreases2, and the root defect falls16→12. This sibling
link is not a recursively reusable ternary production by itself.

## Why a partially marked triangle removes the divisibility obstruction

For every actual port subset of a directed triangle, with k marked vertices,

```
sigma(C3,W)=(3,12,3k,3k,k(k-1),k,k),       0≤k≤3.
```

A homogeneous cyclic triangle with k=0 or3 retains its order3 rotation in
the full attachment, forcing H divisible by3. For k=1 or2 the new root's
incidence distinguishes vertices of the triangle, so that rotation no longer
extends while fixing the exterior. The obstruction's symmetry hypothesis
has genuinely been removed.

As the smallest control within the nonempty-prefix/terminal-singleton scheme,
take `B=C3→{z}`, keep z marked and mark k triangle vertices. Its signature is
`(3,24,12+3k,3+6k,k(k+3),3,k)`. A one-vertex outer prefix gives strong graphs
of order6 with counts `27,29,29,27` for k=0,1,2,3. More generally their counts
are `3-k(k+3)+(12+3k)2^t`, nonzero modulo3 exactly for the two partial port sets.

## Exact finite controls and remaining boundary

`check_two_digit_entry.py` checks all eight marked subsets of the literal
triangle, the three named suffix signatures, and both full construction
schedules at every t=1,...,8. It checks all completed strongness conditions,
old-name/edge retention, both physical defect gauges, and the literal sibling
endpoint equality. The root name in its receipts is99; prefix names start
at10. Normal and optimized receipts are `two_digit_entry.json` and
`two_digit_entry_optimized.json`.

The all-t statement follows from the proved ordered-sum and attachment
identities, not extrapolation from those checks. The missing next theorem is
closure of a specified actual component-word class containing both outputs.
The two arithmetic arrows above do not establish that closure.
