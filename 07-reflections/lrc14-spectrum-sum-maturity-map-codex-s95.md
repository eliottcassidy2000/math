# LRC14 Spectrum Sum Maturity Map -- Codex S95

The current proof stack has a visible direction of travel.  Each false shortcut
failed because it scalarized the object before the cancellation had a chance to
act.

The bounded-denominator route failed first in principle: covering AP-core rows
can force the witness denominator up an unbounded ladder.  That rules out the
dream of a fixed rational atlas.  The Paley-7 story then failed in scale: the
QR/NQR structure is real and useful, but the floor covariance is not dominated
by the apex-7 shell.  The absolute Fourier tail also failed; only signed L2
keeps the right size.  Earlier sector-route work found the same lesson in
another language: per-band, per-cell, and monotone compression proofs all lose
the aggregate compensation.

The object left after these cuts is a signed low-frequency covariance:

```text
SPEC = sum_{n != 0} chat(n) conj(ghat(n)),
R' = 1 + SPEC / baseline.
```

The S95 consecutive-channel atlas makes this more concrete.  Across all
`2380` consecutive bounded-core rows, the worst exact floor is

```text
k=9, P={1,3,4,5}, R'=416640/779291.
```

The Paley ledger is positive there.  The damage comes from the joint residue-0
trunk and the nonzero-shell mean.  This repeats the S94 warning at atlas scale:
QR/NQR balance is a diagnostic, not the obstruction.

The abstract reframe I think is worth pursuing is:

```text
complement-even low-frequency covariance + complement-odd L2 fluctuation.
```

The just-pulled HYP-+2869 proof-assembly note is the optimistic version of the
same picture: universal Farey floor, `R' >= c`, LRC(<=13) for the small part,
then THM-565.  The S95 atlas should be read as supporting infrastructure for
the `R' >= c` line in that assembly, especially the finite bounded-core
low-channel side.

That is the nearest common language for the recent signals:

- KPS S35: binding pairs sum to `14`, so the extremal locus is controlled by a
  complement-like involution.
- S94/S95: the dangerous AP rows have reflection-symmetric bounded cores and
  high additive energy.
- HYP-2635/HYP-2676: coherent signed error repeatedly means Freiman/GAP
  structure.
- HYP-2860/HYP-2862: the generic tail is L2/Parseval, the same backbone as the
  cap route.
- HYP-2690's half-tiling quotient: mirror folds are useful for complement-even
  invariants and dangerous for orientation-sensitive data, exactly the warning
  needed here.

So the target should be stated as an involutive low-packet theorem, not another
numeric hunt:

```text
SPEC_low = SPEC_even + SPEC_odd,
SPEC_odd is mean-zero/L2-small,
SPEC_even either has bounded negative mass or forces an AP/GAP finite atlas.
```

The main caution is that `v -> 14-v` must be formalized in the finite-ruler or
apex-grid coordinate before treating it as an actual symmetry of the functions
on `[0,1)`.  Used carefully, though, it gives a plausible name for the hidden
object: the LRC floor is the positive mass left after quotienting a signed
low-frequency covariance by complement reversal.
