# LRC14 Multidenominator rho_D / Attacker Floor Reflection

HYP-3530 is a useful reframing rather than a closure.

The `rho_D` experiment says the global-witness threshold `1/7` has a rich
finite denominator reservoir on the hard THM-530 rows.  This is the right
successor to the old single-grid/multi-denominator survivor obstruction: the
object is no longer "does one chosen grid survive?" but "which denominator
families carry exact rational points in the phase-flow carrier?"

The caution is just as important.  On the known via-max refutation family,
continuous `rho*_{2/7}=0`, but `D=14` still has isolated grid witnesses.  So a
positive `max_D rho_D` cannot be read as a positive-measure statement.  It is a
certificate bank.  A theorem using it must retain at least one of:

```text
open interval witness
boundary witness with exact sidecar alignment
reconstructed route sidecar
```

This matches HYP-3529's lesson: the live proof should not chase a sharp scalar
constant when a finite sidecar or bounded remainder is enough.

On the union-bound side, the cleanest new packet is the attacker vocabulary.
For k=8..11 the threshold `thr_k=1-gp_min(k)` is so low that the exact
bounded-core bank has no attackers and not even near-attackers.  Consecutive is
the minimizer throughout the bank:

```text
k=8  span<=13  min=691/735
k=9  span<=14  min=247/294
k=10 span<=15  min=38/49
k=11 span<=16  min=1381/2205
```

The next proof target should be a large-span gentle lemma, not another random
scan.  A promising shape is:

```text
if span(E)>k+5, then either
  E has a tail gap that opens a 1/7 phase gap on a controlled denominator,
or
  E compresses to a bounded-core shadow without decreasing mu_1/7.
```

That would turn the bounded exact bank into an unconditional k<12 union-bound
floor.  It would also connect naturally to the finite `rho_D` bank: the same
large gap that makes a row gentle should produce a small-denominator or 14m
phase-flow witness.
