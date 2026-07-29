# All-hard exact global pair-cap census

**Status: THM-2901 PROVED + FINITE-EXACT + VERIFIED.**

On all `14,806` scalar-hard marked suffixes from THM-2899, the global
allowed pair-union maximum is attained and sealed by a second discrepancy
cutoff.  Every branch has positive self-seal gap

```text
H_0-q_1-h/7>0.
```

The exact branch partition is

```text
q_5+2beta_2<h direct closure                  1,835
beta_2<4h/7 finite H3 pair-child route       12,919
pair-cap exception                               52
                                                ------
                                                14,806.
```

The direct branch test closes five whole roots, bringing the current
THM-2895/2898/2899/2901 union to `15` and the official residual to
`3417`.  The `H_3` line is a finite child workload, not a closure.

The exact engine evaluates `212,869` literal pair unions against a raw
finite-head universe of `1,967,632,698`, and independently rechecks every
winner by direct residual subtraction.  The `52` exceptions are a sharp
hostile bank: `51` occur at rank one and one at rank two.

The next scoped target is the `61` non-direct members of the `65`-root
one-hard bank.  Pair flags should be generated from actual membership in
`H_3`, not from cutoff-universe binomial counts.  For a pair `L` and child
`R=C minus D_L`, the heaviness of every parent triple `L union {z}`
forces all three child labels into the residual link core

```text
J_L={z: |R intersect D_z| >= |R|-beta_2}.
```

This restricted link, with the ordered-prefix sidecar retained, is the
proper child object.

Artifacts:

```text
script SHA-256
7ba8244d8fc78ebc0d9381e05d69ca53c849d6008ff9cfb43f0efcbb4b394f81

full branch ledger SHA-256
5dea0eaa45dd52fbf1bef7cfcc328899a4789bc277b6e1e8ac2f4bdf192b85e4

semantic row digest
f3f63ac1953c0e2292488d91f59f831e0f7273b9c1eaeda32f74932d974e4ee4
```
