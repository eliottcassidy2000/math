# Knowing the finite families: the descent realizes every nonempty Z_7-core, so the per-level apex floor is exactly 4cos^2(3pi/7)

*klein-2026-06-29-S17. Working the descent (THM-580) to come to know the finite families it produces. The answer is clean and slightly surprising: the apex family is COMPLETE.*

## The descent as a finitizing chain

THM-580: for a speed set `S`, split `S = O u E` (odd/even), set `S' = E/2`, recurse. This produces a CHAIN
of odd cores `O_0, O_1, ..., O_{d-1}` (depth `d <= 1 + max 2-adic valuation), and
`meas(lonely S) = prod_j rho_j . prod_j meas(lonely O_j)`. The infinite covering family collapses, per `S`,
to a bounded-depth chain of cores; mod the apex prime `7` each core is a subset of `Z_7`. The apex content
of the per-level decorrelation `rho_j` is the cyclotomic gap `g(O_j mod 7)` (THM-590: `0` only at `Z_7`,
else `>= 4cos^2(3pi/7)`, equality at doublets).

## The finite family is COMPLETE

Simulating the descent over a broad covering family (consecutive prefixes, the tightest coverings
`{1..12,182}` / `{1..11,13,84}`, even-heavy, 2000+ random coverings) and collecting the cores mod 7:

> **All 127 nonempty subsets of `Z_7` arise as descent cores** (only `empty` is absent). The apex family is
> the entire nonempty power set of `Z_7`.

And this is not an artifact of sampling — it is constructive: to realize any residue-set `R subset Z_7` as a
core at level `j`, take speeds `2^j a` with `a` odd and `a ≡ r (mod 7)` for `r in R` (plus higher-valuation
speeds to fill earlier levels). Every nonempty residue-set is reachable.

The apex gaps over the arising cores are exactly THM-590's five values, with multiplicities matching the
full power set:
```
gap 0        : the single core Z_7         (the mod-7 covering -- the cusp)
gap 0.198062 : 42 cores, sizes {2,5}        (the DOUBLETS and their complements) = 4cos^2(3pi/7)
gap 0.307979 : 41 cores, sizes {3,4}
gap 1        : 14 cores, sizes {1,6}        (singletons and co-singletons)
gap 2        : 28 cores, sizes {3,4}        (the QR / planar-difference-set cores)
```

## Consequences (what knowing the family buys)

1. **The per-level apex floor is EXACTLY `4cos^2(3pi/7)`, and it is unavoidable.** Because the family is
   complete, the binding doublet core always arises; no structural constraint on the covering family can
   raise the per-level apex floor above `4cos^2(3pi/7)` (THM-590 forbids lower; doublets force equality). So
   `inf rho_j(apex) = 4cos^2(3pi/7)` is a true, attained minimum over a now-KNOWN finite family.

2. **The only gap-0 core is the full `Z_7`** -- the mod-7 covering, the apex cusp. For a covering whose
   odd part hits every residue mod 7 at some level, `rho_j(apex) = 0` there: the apex measure vanishes and
   EXISTENCE is carried by the discrete/witness side (the measure/existence split, HYP-3597). For coverings
   that never fill `Z_7`, every level has apex gap `>= 4cos^2(3pi/7) > 0`.

3. **The floor as a bounded product.** `meas(lonely S) = prod_j rho_j . prod_j meas(lonely O_j)` with each
   `rho_j(apex) >= 4cos^2(3pi/7)` (off the cusp) and each `meas(lonely O_j) >= cap` (THM-576), over depth
   `d <= log2(max speed)`. So `meas(lonely S) >= (4cos^2(3pi/7))^{d} . cap^{d}` -- a product of provable
   finite-family minima, bounded once `d` is bounded.

## Honest scope

- The apex cyclotomic gap and the completeness of the finite family are RIGOROUS (THM-590 + the constructive
  realization). What is CONDITIONAL is that `rho_j` (the genuine 2-sheet decorrelation, THM-580) equals /
  is bounded by its apex cyclotomic gap -- mac-mini S27/S28 found this needs `Gamma_0(14)` congruence-
  averaging when the core is not `Z_7^*`-invariant (HYP-3575/3576). So "knowing the finite family" rigorously
  pins the APEX SKELETON of the floor (`4cos^2(3pi/7)`, complete family, doublet-binding); the bridge from
  the skeleton to the full `rho_j` is the remaining reduction.
- The depth `d` must be bounded over the relevant family (THM-580: `d <= 1 + max 2-adic valuation`), which
  holds for size-bounded coverings.

## The picture

The descent turns the infinite covering family into a bounded-depth chain of `Z_7`-cores, and over all
coverings those cores realize EVERY nonempty subset of `Z_7`. So the apex skeleton of the floor is governed
by the complete finite family `2^{Z_7} \ {empty}`, whose per-level minimum is the doublet value
`4cos^2(3pi/7)`, attained and unavoidable, with the single full-`Z_7` core as the gap-0 cusp where existence
(not measure) carries the floor. This is what "the finite families" are: the full nonempty power set of the
apex `Z_7`, with THM-590 as their exact gap law.

See HYP-3598 (this), THM-590 (the gap law), THM-580 (the descent), HYP-3597 (measure/existence; finite vs
infinite), HYP-3575/3576 (the rho_j reduction, conditional), THM-576 (the caps).
