# Independent review of the Farey/native-conductance packet

23 October 2026. Reviewer: the original-guard/attachment lane.

**PASS, with the scope stated in the reviewed packet.** Read
[`FAREY_MASSES_AND_NATIVE_CONDUCTANCE.md`](../farey/FAREY_MASSES_AND_NATIVE_CONDUCTANCE.md)
and its complete checker. Replayed the frozen receipt in normal and optimized
Python; both passed **81,136** checks.

The all-size proof was checked independently of those computations:

- The incidence coupling is `D_r^-1 W D_c^-1 W^t`, and target elimination
  gives `L=D_r-W D_c^-1 W^t`. The native identical-row lump discards a
  zero eigenmode, so its remaining two-state gap really is the full gap.
- Both actual order-four reductions have raw conductance `4/5`; their
  masses `(3,4)` and `(2,6)` give gaps `7/15` and `8/15`. Exclusive leaves
  retain mass while adding no off-diagonal conductance. The factor seven
  between raw and probability-normalized conductance is correctly retained.
- Removing actual presentation code `2` deletes the entire bridge. The
  two-disjoint-copy hostile correctly shows why strong lumpability alone
  does not authenticate the full spectral gap.
- For reduced `g=n/(m+n)` with odd denominator and `0<g<1/2`, the three
  inverse ratio intervals select a unique positive parent. Coprimality,
  opposite parity and strict decrease of `m+n` exclude the only boundary
  obstruction and terminate at `(2,1)`. The rational branch formulas are
  correct. This proves completeness and uniqueness, beyond the finite box.
- A two-state reversible PSD chain is determined by `(gap, stationary
  masses)`, not by the gap alone. The cited hitting-time comparisons and
  the retained rejection-count clock make this information loss explicit.

The result does **not** claim that Berggren branches are native tournament
moves, that a conductance quotient authenticates the original mod-eight
guard, or that an arithmetic decoder supplies a deterministic Collatz
clock. The document retains these boundaries. Its new native statement
is the exact shared conductor and mass explanation, while the unbounded
ternary statement concerns rational spectral coordinates.

Reproduction:

```text
python -B 05-knowledge/results/native_guard_moves_farey_20261023/farey/check.py
python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/farey/check.py
```
