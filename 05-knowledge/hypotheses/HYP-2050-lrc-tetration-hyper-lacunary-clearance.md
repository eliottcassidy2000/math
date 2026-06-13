---
id: HYP-2050
status: SUPPORTED
source: codex-2026-06-01-S551
related:
  - HYP-2039
  - HYP-2041
  - HYP-2042
  - HYP-2047
  - HYP-2048
  - HYP-2049
  - HYP-2051
---

# HYP-2050: tetration speeds are the hyper-lacunary clearance extreme

**Claim.** On the LRC hyperoperation ladder, level-4 tetration speed families
are the hyper-lacunary extreme.  They should be read as stress tests for the
level-3 cascade product:

```text
runner clearance factors c_k       = level-2 local conditions,
prod_k c_k                         = level-3 repeated-x closure,
tetration speeds a^^k              = level-4 decorrelation stress.
```

The product is not a metaphor: after ordering the runners,

```text
F_k = F_{k-1} cap {||v_k t|| >= 1/n},
c_k = |F_k| / |F_{k-1}|,
|SAFE| = prod_k c_k.
```

Thus the cascade is the repeated multiplication of the runner clearances.  The
tetration experiment asks what happens when the speeds are generated far above
the additive/AP rung and even above the geometric rung.

## Evidence

`04-computation/lrc_tetration_hyper_lacunary_s551.py` compares six-runner
families on a midpoint grid with threshold `1/7`:

```text
L2_AP          (1,2,3,4,5,6)
L2_Fibonacci   (1,2,3,5,8,13)
L3_geom2       (1,2,4,8,16,32)
L3_geom3       (1,3,9,27,81,243)
L4_tetr2       (1,2,4,16,65536,2^^5)
L4_tetr3       (1,3,27,3^27,3^^4,3^^5)
```

The sampled clearance ledger gives:

```text
family          product    product/credit    mean H entropy
L2_AP           0.0000     0.0000            4.2321
L2_Fibonacci    0.0488     0.3672            4.3702
L3_geom2        0.2319     1.7458            3.4193
L3_geom3        0.1557     1.1724            3.9376
L4_tetr2        0.1671     1.2582            3.7503
L4_tetr3        0.1873     1.4102            3.1556
```

The AP has the expected zero cascade factor.  Fibonacci remains additive and
resonance-rich.  Geometric families open the cascade.  The base-2 geometric
family has the largest sampled product, while the base-3 tetration family has
the lowest mean half-turn tournament entropy.

Tournament Analysis on the family graph is transitive in this finite probe:

```text
L3_geom2 > L4_tetr3 > L4_tetr2 > L3_geom3 > L2_Fibonacci > L2_AP
```

The fingerprint is `score_hist={0:1,1:1,2:1,3:1,4:1,5:1}`,
`directed_3_cycles=0`, singleton SCCs, and one Hamiltonian path.

## Interpretation

Tetration is not merely "large speeds."  It is the regime where every new
speed is built by iterating the previous multiplicative closure so aggressively
that ordinary additive resonances are mostly skipped.  In LRC terms, later
tower runners tend to add no AP-like final zero factor.  In the S551 sample,
the final `L4_tetr3` tower runner has sampled conditional clearance `1.0000`
after the previous cuts.

The caution is just as important: finite moduli fold towers.  `L4_tetr3`
lowers entropy without automatically maximizing the product; it still has
visible danger-overlap (`1.1825`).  The right proof target is not "huge speeds
are independent"; it is:

```text
finite tower aliasing cannot create an AP-like zero factor outside
the additive/tight regime.
```

## Assumption challenge

Candidate vertices included runners, clearance factors, cascade prefixes,
speed families, gap cells, danger arcs, Fourier modes, residues, and proof
obligations.  The script uses speed-family vertices for the summary
tournament and clearance factors inside each family.

This quotient preserves sampled lonely measure, product-vs-credit behavior,
and half-turn tournament entropy.  It destroys exact endpoints, exact
set-vs-measure boundary witnesses, and residues beyond the chosen moduli.
The challenged assumption is that "larger speeds automatically make LRC
harder."  The finite probe suggests the opposite for hyper-lacunary towers:
they are easier by product and entropy, with the remaining debt concentrated
in modular aliasing.

## Predictions

1. Level-4 families should have lower mean tournament entropy than additive
   families at comparable runner count.
2. Any tetration obstruction should appear as a finite-modulus aliasing event,
   not as an additive AP-style resonance chain.
3. Cascade proofs should separate two debts: ordinary additive returns that
   can drive a zero factor, and tower aliasing that may raise overlap without
   closing the lonely set.
4. The useful proof quotient is likely a prefix-clearance or aliasing-event
   tournament, not a raw runner tournament.

## Files

`04-computation/lrc_tetration_hyper_lacunary_s551.py`  
`05-knowledge/results/lrc_tetration_hyper_lacunary_s551.out`  
`07-reflections/lrc-tetration-hyper-lacunary-s551.md`
