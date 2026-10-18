# Pre-registration: is there a largest prime where the natural Paley order survives every rotation order? (HYP-9287 revised)

2026-10-17, opus-2026-10-17-S40. Written and pushed **before** the run.

## Background

The Farey sweep (`paley_rotation_farey_sweep_20261017.py`) computes the best of all rotation orders `x -> {alpha x}` of
the Paley tournament `P_p` exactly, in `O(p^2)`. It rests on the jump formula
`Delta_b = 2 chi(b)[S(q) + (b - s)A(q) + s A(q + 1)]`, checked against direct excesses on every Farey interval for
`p <= 47`, and on the identity `sum phi(b) Delta_b = -2 p h(-p)`, checked for `p <= 400`. It reproduces every audited
rotation certificate.

Up to 479, the natural order survives the rotation family exactly at 3, 7, 23, 47, 71, 191, 239, 311, 359, 431, 479 (with
11, 19, 31, 43 beaten by other orders). The audit's 599 certificate shows that the survivor boundary rises quickly with
`p`.

## Prediction

Run: `python paley_rotation_farey_sweep_20261017.py range 480 20000 2.4` (all 245 primes `p ≡ 3 (mod 4)` in
`(479, 20000]` with `L(1, chi_(-p)) >= 2.4`; a prime with smaller `L` is beaten by the local search at every tested size,
and its natural excess is below the generic level).

* **F (finite-list form).** The natural order is beaten by some rotation order at **every** one of the 245 primes. In
  particular 479 would be the largest survivor up to 20000 (among `L >= 2.4`). Passes if there are no survivors.
* **Alternative (thin infinite set).** If survivors appear, they are expected at the largest `L` (the 87 primes with
  `L >= 3.6`, up to `L = 4.54`), and their list will be recorded as the revised conjecture's evidence.

Either outcome is reported. Surviving the rotation family is necessary, not sufficient, for optimality; being beaten is
a certificate (an explicit `alpha` is checked by a direct excess computation).

## Outcome (appended after the run)

Run as pre-registered (`paley_rotation_farey_sweep_480_20000_20261017.out` / `.json`): **F FAILED.** 46 of the 245
primes keep the natural order as their best rotation order (set R+, listed in HYP-9287), all with `L >= 3.58`; 199 are
beaten, each with a certificate `alpha` checked by a direct excess computation. The survivor boundary rises slowly with
`p` (minimum survivor `L` 3.58 below 1000, 4.35 near 20000) and overlaps the beaten range, so the rotation-survivor set
looks infinite and thin rather than finite. The thin-infinite alternative of the pre-registration is the one recorded.
