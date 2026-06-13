# Cartan Dimension Ladder: Full vs Traceless Dark Modes

**Date:** 2026-05-29
**Instance:** codex-2026-05-29
**Scripts checked:** `04-computation/sixteen_dimensions_s93d.py`, `04-computation/cartan_attention_theorem.py`
**Tags:** #cartan #dark-modes #attention #gl-n #napolitano

## Result

For `gl(n,R)` with the transpose Cartan decomposition,

```text
gl(n,R) = so(n) + Sym(n)
```

the dimensions are:

```text
active(n)        = dim so(n)  = n(n-1)/2
full_dark(n)     = dim Sym(n) = n(n+1)/2
center(n)        = 1
traceless_dark(n)= full_dark(n) - 1
```

Therefore:

```text
full_dark(n) / active(n)      = (n+1)/(n-1)
traceless_dark(n) / active(n) = (n+2)/n
```

For Napolitano's `gl(4,R)` case:

```text
active = 6
full_dark = 10
traceless_dark = 9
full_dark / active = 5/3
traceless_dark / active = 3/2
```

## Norm-Ladder Correction

If `k` normalization layers are modeled as `gl(2k,R)`, then the full symmetric dark-mode count is

```text
D_k = C(2k+1, 2) = k(2k+1).
```

Increasing the norm count by one changes `gl(2k,R)` to `gl(2k+2,R)`, so the added full dark modes are

```text
D_{k+1} - D_k
  = C(2k+3, 2) - C(2k+1, 2)
  = 4k + 3.
```

This corrects the earlier statement that each additional norm adds `2k` dark modes. The `2k` expression is the difference between consecutive symmetric dimensions in `gl(2k-1,R)` and `gl(2k,R)`, not the `k -> k+1` norm-ladder step.

Concrete values:

```text
k=1: gl(2,R),  active=1,  full_dark=3
k=2: gl(4,R),  active=6,  full_dark=10, added dark from k=1: 7
k=3: gl(6,R),  active=15, full_dark=21, added dark from k=2: 11
k=4: gl(8,R),  active=28, full_dark=36, added dark from k=3: 15
k=5: gl(10,R), active=45, full_dark=55, added dark from k=4: 19
```

## Verification

After updating the two scripts:

- `python3 04-computation/sixteen_dimensions_s93d.py` completed successfully and printed the corrected `4k+3` marginal dark-mode formula.
- `python3 04-computation/cartan_attention_theorem.py` completed successfully.
- The attention-derived tournament OCF checks passed:

```text
n=3: 50/50
n=4: 50/50
n=5: 50/50
n=6: 50/50
```

## Implication

The `6+10` split for `gl(4,R)` is real but dimension-theoretic. Any claim that the dark sector carries correctness should control separately for:

1. the full symmetric sector's dimensional advantage over the antisymmetric sector,
2. the scalar center if it is included in "dark" modes,
3. the softmax tendency to place extra energy in the symmetric sector before training.
