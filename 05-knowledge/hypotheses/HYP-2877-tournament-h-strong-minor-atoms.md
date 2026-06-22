---
id: HYP-2877
status: SUPPORTED by exact fixed-path audit n<=7; open as a cofinite strong-atom realization theorem
source: codex-2026-06-22-S98
tags: [tournaments, graph-minors, kuratowski, wagner, robertson-seymour, h-spectrum, forbidden-H, strong-components, beurling-selberg, lrc14, tournament-analysis]
related:
  - HYP-2872
  - HYP-2873
  - THM-520
  - THM-115
  - THM-083
  - THM-116
  - THM-200
  - HYP-2870
  - HYP-2868
---

# HYP-2877: tournament H has a Wagner-style strong-atom minor lens

The even-graph quotient from HYP-2872 is the right coordinate chart and the
wrong minor category.  The `H`-preserving simplification is instead
strong-component condensation with a multiplicative atom ledger.

## Exact factor law

Let `C_1,...,C_r` be the strong components of a tournament `T`, ordered so that
every inter-component edge points from `C_i` to `C_j` for `i<j`.  The
condensation of a tournament is transitive, so every Hamiltonian path of `T`
must visit the strong components in this order.  Conversely, concatenating a
Hamiltonian path inside each `C_i` gives a Hamiltonian path of `T`.  Therefore

```text
H(T) = product_i H(C_i).
```

This gives two safe simplifications.

1. Singleton strong components may be suppressed: they contribute factor `1`.
2. A nontrivial strong component may be contracted only as a labelled atom
   carrying at least its label `(size, H)`, or just `H` if the only target
   invariant is multiplicative attainability.

Forgetting the label is exactly the mistake made by arbitrary even-graph
contraction/smoothing.  The operation can simplify the address while destroying
the `H` proof state.

## Exact audit through n=7

The script `04-computation/tournament_h_strong_minor_lens_codex_s98.py` audits
all rooted fixed-Hamiltonian-path tournaments through `n=7` (`33868` rows).

- `H` factorization failures: `0`.
- Singleton strong-component suppression failures: `0`.
- Strong-atom signatures checked: `89`.
- Signatures with multiple `H` values: `0`.
- `{7,21}` occur in neither the rooted spectrum nor the observed strong-atom
  semigroup closure.
- The observed strong-atom gaps up to `189` are exactly the rooted odd gaps:
  `7,21,63,107,119,149,161,163,165,167,169,173,177,179,181,183,185,187`.

The last line is useful but should not be overread.  S33 already corrected the
global picture: the permanent forbidden set is finite `{7,21}`, not the ideal
of all multiples of `7`.  Strong tournaments with `H` divisible by `7` already
exist at `n=7`, including `35,49,77,91,133,147,175`.  Thus the obstruction is
not divisibility by `7`; it is failure to realize the two small value atoms
`7` and `21`.

## Candidate theorem

The clean cofinite theorem suggested by the audit is:

```text
For every odd h not in {7,21}, h is attainable by multiplying labelled
strong-component H-atoms, and ideally by one strong atom except for the
obvious composite constructions.
```

Equivalently, the `{7,21}` theorem should be proved in the labelled
strong-atom category rather than in the even-graph shadow.  This is the
tournament analogue of moving from an equinumerous encoding to a
Kuratowski/Wagner-style preservation theorem.

## LRC14 / Beurling-Selberg connection

The same lesson applies on the LRC analytic side.  Beurling-Selberg and Fejer
trigonometric carriers are useful because they keep explicit labels:
bandlimit, one-sided defect, and finitely many low resonant coefficients.  A
minorant for `1_{G_P}` or a majorant for a cover event is not a lossless
quotient unless the defect is carried forward.  The defect label plays the
same role as the strong-component `H` label.

This gives a conservative proof heuristic for LRC(14):

- finite low Fourier/resonant modes are the analogue of a finite forbidden
  minor set;
- the high-frequency tail is the analogue of a monotone closure theorem;
- Beurling-Selberg certificates are useful when the proof state records the
  defect explicitly, not when they are treated as exact replacements for the
  original event.

This also explains why the additive-energy/Fejer signal from S33 is promising
but not a closure by itself: `A(E)=integral |sum_e exp(2 pi i e t)|^4` says
intervals are maximally spectrally concentrated, while the LRC `L_y` functional
is a signed/weighted coverage moment.  The Fejer concentration is a proof
carrier to be labelled and combined with the low-mode defect ledger.

## Assumption challenge

I explicitly considered runners, arcs, fixed circle sections, even-graph edges,
degree-2 smoothing cores, arbitrary contractions, strong components, OCF
packets, `H` values, Fourier modes, Beurling-Selberg defects, and proof
obligations as vertices.  The chosen quotient is strong-component `H` atoms.
It preserves multiplicative `H` attainability and destroys the GF(2)
cycle-space geometry.  HYP-2872 proves that the destroyed geometry is too
lossy for `H`; HYP-2877 keeps the ledger that `H` actually respects.

## Artifacts

- `04-computation/tournament_h_strong_minor_lens_codex_s98.py`
- `05-knowledge/results/tournament_h_strong_minor_lens_codex_s98.out`
- `07-reflections/tournament-h-strong-minor-atoms-codex-s98.md`
