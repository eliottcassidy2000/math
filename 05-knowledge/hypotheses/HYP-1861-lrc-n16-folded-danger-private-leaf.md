---
id: HYP-1861
status: OPEN
source: codex-2026-05-31-S393
related:
  - HYP-1841
  - HYP-1842
  - HYP-1856
  - HYP-1857
  - HYP-1858
  - HYP-1859
  - HYP-1860
---

# HYP-1861: n=16 folded danger forces a private dyadic endpoint leaf

## Statement

At LRC denominator `n=16`, fold the time circle by the antipodal map
`t -> t+1/2` and write `s=2t`.

Even speeds descend to the quotient: speed `2w` is unsafe for both antipodal
preimages exactly when

```text
||w*s|| < 1/16.
```

Odd speeds are one-sided.  If `d=||u*s/2||`, then an odd speed `u` kills the
`t=s/2` side only when `d<1/16`, and kills the `t=s/2+1/2` side only when
`d>7/16`.

The folded antipodal test fails only on the **folded-danger** set:

```text
even_safe(s) and some odd_low(s) and some odd_high(s).
```

The hypothesis is:

```text
For primitive 15-speed n=16 systems satisfying the small-denominator sieve,
folded-danger can occur only with a private dyadic endpoint leaf.
```

Equivalently, the folded parity obstruction and the endpoint-core obstruction
should not both be defeated.  If folded parity alone does not produce a lonely
antipodal witness, the endpoint-protection graph should peel.

## Evidence

`lrc_n16_proof_gauntlet_s393.py` tried five further attacks beyond S389/S390/S391/S392.

The folded antipodal quotient is useful but not sufficient.  The initial
segment and the single `16`-gate row have zero folded witness mass, so the fold
can lose boundary-only information.  However, hard structured rows often retain
positive folded witness intervals:

```text
best 8-ladder: folded_witness=101/8580
min-cover seed: folded_witness=749/22880
S390 random near: folded_witness>0
```

The normalized `n=16` half-turn cube was exhausted exactly.  In the finite
scalar quotient with `1152` alpha patterns and `18432` candidate cells, the
best pure two-torsion defects still miss many cells:

```text
support 1: best missed 128, coordinate 15
support 2: best missed 160, coordinates 10 and 15
support 3: best missed 256
```

The local endpoint cover for a `16`-gate is rigid.  Covering all `32` endpoints
of `v=16` using lower residues requires exactly nine protectors:

```text
(1,3,5,7,8,9,11,13,15).
```

Treating that exact local cover plus the `16`-gate as an adversarial seed,
S393 gave it five more speeds from a structured pool.  The best completions
were sieve-complete but still positive-gap and core-empty:

```text
best min-cover completion:
  speeds=(1,2,3,5,7,8,9,10,11,12,13,15,16,42,120)
  fixed gap=1/384
  exact gap/th=0.041667
  unprotected=36
  coreE=0
```

A second sieve-aware beam from seed `{1,16}` found stronger high-coverage
sets.  The closest row was:

```text
speeds=(1,3,4,10,13,14,16,17,18,19,22,23,31,41,60)
fixed_gap=43/39360
exact gap/th=0.017480
unprotected=38
coreE=0
folded_witness=0
```

This is the important failure mode: some rows defeat the folded antipodal
witness test, but they remain positive-gap and endpoint-core-empty.  That is
exactly the pattern predicted by folded-danger forcing a private endpoint leaf.

## Predictions

1. A folded-danger `n=16` system should have a dyadic endpoint leaf in the
   first few peel layers, even when the small-denominator sieve is complete.
2. The nine-residue local cover of a `16`-gate should be a rigidity pattern:
   every five-speed completion has positive fixed-threshold gap or a private
   endpoint leaf.
3. The exact half-turn cube moat should be promotable to a compression lemma:
   rounding a near-dangerous residue vector to dyadic half-turn data does not
   increase coverage enough to erase the `128`-cell moat.
4. A proof of `n=16` should combine folded parity, local-cover rigidity, and
   endpoint-core peeling.  None of the three appears sufficient alone.

## Sources

- `04-computation/lrc_n16_proof_gauntlet_s393.py`
- `05-knowledge/results/lrc_n16_proof_gauntlet_s393.out`
- `07-reflections/lrc-n16-proof-gauntlet-s393.md`
