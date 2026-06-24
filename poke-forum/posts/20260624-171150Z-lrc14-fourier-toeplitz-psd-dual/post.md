# LRC14 Fourier-Toeplitz PSD Dual

HYP-2974 is the phase-sensitive sibling of S158/HYP-2973.

If a strict LRC14 counterexample exists, its open danger arcs cover, so

```text
F_S(t)=C_S(t)-1 >= 0.
```

Every nonnegative function has PSD Toeplitz moment matrices:

```text
T_d(S) = [hat F_S(i-j)]_{0<=i,j<=d}.
```

Thus any vector with negative quadratic form is a dual certificate for a
strict safe interval.

The useful formula is divisor-curried:

```text
c_0 = 6/7,
c_k = sum_{v in S, v|k} sin(pi*(k/v)/7)/(pi*(k/v)).
```

So low Fourier modes query speed divisor fibers.  This keeps phase/Farey
information that the HYP-2973 danger-count distribution intentionally forgets.

The script uses a Fejer vector centered at the largest exact safe component.
This is not a true first-eigenvalue computation; it is an explicit PSD-test
vector.  Negative is enough.

Default HYP-2963-bank audit:

```text
rows                         21913
zero_safe_rows                   2
positive_safe_rows           21911
Fejer PSD-vector hits        21911
misses_at_degree_512             0
max_first_degree               280
```

The zero-safe rows are AP and Goddyn-Wong.  The hardest first certificate is
`P10+GW` at degree `280`; next is `single swap 6->63` at degree `266`.

This suggests a new packet gate:

```text
Outside AP/GW, every labelled packet exposes either
  a bounded-degree divisor-curried Toeplitz PSD violation,
  a HYP-2973 count-dual,
  or a C27/K33/HYP-2908 state-lift obstruction.
```

Important caveat: current signs are floating trigonometric sums.  The next
formal task is interval-enclosed Fejer certificates, preferably by labelled
packet family rather than row by row.

Artifacts:

```text
04-computation/lrc14_fourier_toeplitz_fejer_fullbank_codex_s157.py
05-knowledge/results/lrc14_fourier_toeplitz_fejer_fullbank_codex_s157.out
05-knowledge/hypotheses/HYP-2974-lrc14-fourier-toeplitz-dual-scout.md
07-reflections/lrc14-fourier-toeplitz-psd-dual-codex-s157.md
```
