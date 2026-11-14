# Inverting the selective cut from an absolute target

This follows the [changing-source and reader session](
../chanpak_changing_sources_readers_20261114/README.md).
The output here is a constructive inverse within physical charts
plus an all-order sparsity obstruction. **Not a Chan–Pak proof.**

* [Dyadic target decoder](DYADIC_DECODER.md): the target's binary
  valuation picks a constant-size list of possible height pairs.
  No factorization, height scan, or post hoc choice of mask is needed.
  Two proved physical charts give the pair
  (120·2^{k+m}−57·2^k−78·2^m+{71,73}) and the triple
  (432·2^{k+m}−303·2^k−174·2^m+{225,227,229}).
  They start at 281,283 and 999,1001,1003 respectively. For
  example 13897 decodes to the second C5 mask at heights (4,3).
* [Independent physical switches](SWITCH.md): the three marked
  six-vertex cores have identical cover and endpoint aggregate
  weights. Their outputs stay exactly two apart after arbitrary
  painted transitive tapes on both sides. In a general ordinal
  word a local switch is scaled by the Hamilton count of the
  other blocks, so repeating one core does not give unit digits.
* [Four-coordinate state](CUT_STATE.md): a packing
  recurrence for all painted ordinal continuations. Four
  scalar contexts, each adding at most two vertices, recover
  the full continuation state. The earlier cut-vertex balance
  transfer gives a related way to track full path distributions.
* [Four-block inverse and quantitative sparsity](
  SPARSITY_AND_INVERSE.md): every two-contact base whose two
  contacts are singleton components has the indicated normal
  form, up to an initial multiplicative prefix. One can invert
  it once a physical terminal block is specified. But allowing
  source H up to 2^(0.5 sqrt(log₂ X)) and at most 0.1 log₂ X arcs
  on the smaller contact side still yields at most
  X^(0.916827…+o(1)) outputs up to X.
* [Singleton-contact mod-three skeleton](SINGLETON_CONTACTS.md):
  replace every nonempty losing gap by one or two transitive
  vertices according to parity. The scalar modulo three then
  factors as H(base) times the output of the reduced word. This
  gives an exact rejection test for all such marked sources,
  with any number of contacts, not an absolute-target guarantee.

## What remains

The decoder chooses cut, reader, and dimensions from the target
**within** its families. For a general target one still needs a
physical source (and possibly enough selective contacts) with
the requisite endpoint-profile coefficients. The small-complexity
theorem limits a suggested way to prove that supply; it does not
exclude witnesses using a larger source count or more contacts.
A next step is to vary the marked six-vertex profile using larger
nonordinal cores without compromising the cheap inverse; merely
stacking identical ordinal switches multiplies the step.

Run:

    python replay_dyadic.py
    python replay_state.py
    python replay_switch.py
    python replay_four_blocks.py

The scripts verify actual oriented graphs by subset path counting
and exercise exact inverse decisions. The all-height claims and
asymptotic bound are proved in the linked notes, not extrapolated.
