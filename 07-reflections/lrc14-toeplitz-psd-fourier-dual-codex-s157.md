# LRC14 Toeplitz-PSD Fourier Dual

Date: 2026-06-24
Agent: codex-2026-06-24-S157
Related: HYP-2974, HYP-2973, HYP-2972, HYP-2971, HYP-2970, HYP-2969, THM-534

The user pushed the right refinement: from the scalar multiplicity law to the
Fourier positivity of `C_S-1`.  If the danger arcs cover, `C_S-1` is
nonnegative, so its Toeplitz moment matrices must be PSD.

The useful certificate is curried:

```text
S -> fhat_S -> T_d(S) -> min eigenvector c
```

or equivalently

```text
S -> (p -> int (C_S-1)|p|^2).
```

A negative value forces a positive-measure safe interval.  This is cleaner than
another scalar mass estimate because it produces a trigonometric-square
multiplier that may localize near the actual witness interval.

The experiment was small but sharp.  AP and GW remain PSD through degree 160;
all positive hard rows tested fail.  The first failure is sometimes low
(`13->26` at 37, covering rows at 51) and sometimes delayed (`12->36` at 101,
`P10+GW` at 160).  That delay is useful information: it measures Fourier
certificate complexity and can be compared with HYP-2971 barrier degree,
HYP-2973 count-dual degree, and HYP-2972 first twist denominator.

Next work should decode eigenvectors into localized witness packets.
