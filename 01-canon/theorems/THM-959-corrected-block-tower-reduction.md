---
id: THM-959
title: SHARPENED CORRECTED BLOCK-TOWER REDUCTION — any prescribed consecutive partition into blocks of size at most six is lonely when each junction satisfies the direct-THM-955 G0 source/target table (uniform worst case M/m=1, maximum G0(6,6)=1105); the proof is a direct nested-window induction and needs no within-block ratio cap
status: PROVED PAPER + FINITE-EXACT ARITHMETIC. This replaces the collided and overstated THM-956 block page: the original ratio monotonicity was reversed (2030 was too small), window_tail_glue did not compose arbitrary blocks, and the claimed exact description of the complement was unproved. The S62 correction gave the valid but conservative maximum 2198 after weakening both numerator and denominator of THM-955. The sharpened proof below uses THM-955 directly and lowers the maximum to 1105. The direct induction and G0 table are exact; 60/60 generated towers are an implementation referee, not the proof.
source: opus-2026-07-17-S339; corrected and renumbered codex-2026-07-17-S62 after the earlier THM-956 adaptive-q claim; direct-THM-955 sharpening codex-2026-07-17-S62 follow-up
depends_on: [THM-955]
scripts: 04-computation/block_reduction_G0_opus_S339.py -> 05-knowledge/results/block_reduction_G0_opus_S339.out
---

# THM-959 — the corrected block-tower reduction

Let the distinct positive integer speeds be partitioned, in sorted order, into
nonempty consecutive blocks

```text
B_1 | B_2 | ... | B_r,       1 <= k_i:=|B_i| <= 6,
m_i=min B_i,                 M_i=max B_i.              (1)
```

For `1<=k<=6`, choose the following rational entry lengths and define

```text
ell=(3/4,5/4,2,3,6,13),
mu_k=((1-k/7)ell_k-k/7)/(1+k+k ell_k).                 (2)
```

Their exact values are

```text
k       1       2       3       4       5       6
ell_k  3/4     5/4      2       3       6      13
mu_k   2/11   17/154   1/14   5/119   1/36   1/85.     (3)
```

Put

```text
G0(s,t)=ell_t/mu_s.                                    (4)
```

The full exact table is

```text
src\tgt     1        2        3        4        5        6
1          33/8     55/8      11      33/2      33      143/2
2         231/34   385/34   308/17   462/17   924/17   2002/17
3          21/2     35/2      28       42       84       182
4         357/20   119/4    238/5    357/5    714/5    1547/5
5           27       45       72      108      216       468
6          255/4    425/4    170      255      510      1105.   (5)
```

Thus the largest table entry is `G0(6,6)=1105`.  If every junction in
(1) satisfies

```text
m_(i+1)/M_i >= G0(k_i,k_(i+1)),                        (6)
```

then the full speed packet is lonely.

## Proof

For a block of `k` speeds in a safe input window of length `L`, with minimum
`m` and maximum `M`, assume the entry invariant

```text
mL >= ell_k.                                           (7)
```

THM-955 directly gives a nested safe output window of length `W` satisfying

```text
W >= ((1-k/7)L-k/(7m))/(1+k+L sum_(w in B) w).         (8)
```

The numerator is positive by (2)--(3) and (7).  Write `l=mL` and
`rho=M/m`.  Since `sum B<=kM`, multiplying (8) by `M` gives

```text
MW >= rho((1-k/7)l-k/7)/(1+k+k rho l).                (9)
```

For `l>=ell_k` and `rho>=1`, the right side is increasing separately in
`l` and `rho`: its two derivatives, after clearing their positive squared
denominators and the common positive factor `rho`, have numerators

```text
(1-k/7)(1+k)+k rho (k/7),
((1-k/7)l-k/7)(1+k),                                 (10)
```

respectively.  Hence (9) is at least its value at `(ell_k,1)`, namely
`mu_k`.

These rational entry lengths are chosen near the minima of the self-junction
cost `ell/mu_k(ell)`.  In the bottleneck `k=6`, its real minimizer is
`ell=6+sqrt(43)`, where the cost is approximately `1103.825`; the rational
choice `ell_6=13` gives `1105`.  Thus the displayed uniform maximum is within
`1.2` of the optimum available from this one-parameter direct-THM-955 scheme.

The first block needs no entry hypothesis.  On the full unit circle each
danger comb has exact mass `1/7`, so the common safe set has mass at least
`1-k/7`.  It has at most `1+sum B<=1+kM` interval components.  Its widest
component therefore has length `W` with

```text
MW >= (7-k)M/(7(1+kM)) >= (7-k)k/(7(1+k^2)) >= mu_k, (11)
```

where `M>=k` because the block has `k` distinct positive integers; the six
last comparisons are exact rational checks from (3).

Now induct on the blocks.  Suppose the nested window after `B_i` has length
`W_i` and `M_i W_i>=mu_(k_i)`.  From (6),

```text
m_(i+1)W_i >= G0(k_i,k_(i+1)) mu_(k_i)=ell_(k_(i+1)), (12)
```

so (7)--(10) apply to `B_(i+1)`.  Every output window is contained in its
input window, hence remains safe for all earlier blocks.  A point of the final
nonempty window is lonely for every speed. ∎

## Floor transfer

An interval of width `delta` contains at least

```text
ceil(delta D)-1                                        (13)
```

points of every `1/D` grid.  Indeed the number in `(a,a+delta]` is
`floor(D(a+delta))-floor(Da)>=floor(D delta)>=ceil(D delta)-1`.

This is the exact continuous-to-discrete bridge needed by live-multiplier
censuses.

## Scope guardrail

The theorem closes every packet for which a partition satisfying (1) and (6)
is supplied.  It does **not** prove that every packet outside this class has a
single block of seven comparable speeds.  Characterizing the complement of
the admissible-partition language—and then closing it by a seven-wall overlap
floor—remains open.  The natural vertices are blocks and junction obligations;
a runner tournament retains order but destroys block size, entry length, and
the quantitative `G0` labels.

The exact script recomputes (2)--(6), performs 60 generated end-to-end nested
window referees, and checks (13).  Those tests audit the implementation; the
proof is the induction above.  Lean formalization still requires the corrected
positive-numerator THM-955 interval lemma and this finite block induction.
