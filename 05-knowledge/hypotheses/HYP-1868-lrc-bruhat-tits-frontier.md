---
id: HYP-1868
status: OPEN
source: codex-2026-05-31-S400
related:
  - THM-357
  - THM-365
  - THM-366
  - THM-367
  - HYP-1859
  - HYP-1865
  - HYP-1866
  - HYP-1867
  - HYP-1891
---

# HYP-1868: LRC repair branches have positive Bruhat-Tits frontier mass

## Statement

For each prime `p`, view rational LRC endpoints as boundary points of the
Bruhat-Tits tree of `PGL_2(Q_p)`.  If an endpoint is `x=a/b`, define its
relative denominator-cusp depth above the base denominator `n` by

```text
h_p(x;n) = max(0, v_p(b) - v_p(n)).
```

The conjectural p-adic endpoint debt norm is the normalized frontier mass

```text
BT_p(E;n) = sum_{x in E} p^(-h_p(x;n)),
```

where `E` is the set of exposed or debt endpoints.  A genuine repaired LRC
branch should have positive frontier mass in at least one `p`-tree, or in the
adelic product of the relevant trees.

Equivalently, moving endpoints deeper into a denominator cusp can hide the
real gap only by increasing horosphere population.  Thus the HYP-1866 product
law should factor as

```text
ArchGap * raw_debt = (ArchGap * p^h) * BT_p
```

on single-tree branches, with a product-building version for mixed-prime
branches.

## Evidence

The S400 audit computes p-adic frontier profiles for the known debt-export
rows.

At `n=16`, the dyadic rows are a clean single Bruhat-Tits tree translation:

```text
n16 d=2:   h_2=1, raw=34,  BT_2=17,   Gap_A*2^h=2/33
n16 d=4:   h_2=2, raw=68,  BT_2=17,   Gap_A*2^h=2/33
n16 d=8:   h_2=3, raw=140, BT_2=35/2, Gap_A*2^h=2/33
n16 d=16:  h_2=4, raw=280, BT_2=35/2, Gap_A*2^h=2/33
```

So the `35/34` product jump from HYP-1866 is not a failure of tree
translation.  The tree translation factor `Gap_A*2^h` stays fixed; the jump is
a normalized frontier-mass jump from `17` to `35/2`.

The dyadic infinity-cusp residues modulo `16` also stabilize after the jump:

```text
d=8:  residues 1,7,9,15 have count 14; residues 3,5,11,13 have count 21
d=16: the same pattern doubled to 28 and 42.
```

This suggests a possible residue-branch certificate for the phase tax.

At `n=14` and `n=18`, the export rows are not one-tree phenomena.  They already
live in a product of p-adic trees:

```text
n14 d=7:  p=2 depths 0,2; p=7 depth 1
n18 d=9:  p=2 depths 0,4; p=3 depths 2,3
```

After doubling the scale, new p=2 depths appear while the odd payload tree
still carries mass.  This points toward an adelic product building rather than
a single Bruhat-Tits tree for composite denominators.

## Proof Route

The Bruhat-Tits proof target is:

```text
Every all-repaired endpoint core has positive normalized frontier mass
in at least one prime tree, after quotienting base denominator levels.
```

For `n=16`, this should reduce to a dyadic statement:

```text
an all-protected 16-gate branch either has positive real gap,
or exposes a nonzero dyadic frontier in the infinity cusp,
or pays a residue-branch phase tax visible in unit residues mod 16.
```

For mixed denominators, the likely object is a finite support in the product
of the relevant trees.  A counterexample would have to send the real gap to
zero while simultaneously erasing all product-building frontier mass.  That is
the Bruhat-Tits version of "gap=0 and debt=0 simultaneously."

HYP-1891 identifies the source of the mixed-prime complication: `n=14` and
`n=18` are adjacent first-even children in the odd-root `x+2` chain, while
`n=16` belongs to the pure vertical `x*2` chain over odd root `1`.  The former
need odd-payload transfer matrices; the latter needs dyadic row-mode frontier
mass.

## Sources

- `04-computation/lrc_bruhat_tits_frontier_s400.py`
- `05-knowledge/results/lrc_bruhat_tits_frontier_s400.out`
- `07-reflections/lrc-bruhat-tits-frontier-s400.md`
- HYP-1866.
- HYP-1859.
- THM-357.
- THM-365.
- THM-366.
- THM-367.
- HYP-1891.
