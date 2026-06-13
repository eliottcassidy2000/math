# LRC14 Eight-Core Exceptions Open At Shell 41

The eight-core scout was a useful failure.

The tempting theorem was:

```text
in the carry window, retaining 8 of the 12 speeds in 7*{1,...,12} forces Q27.
```

That is false.  The exact four-deletion census found two Q27-feasible deletion addresses:

```text
(28,42,56,84)
(42,56,70,84)
```

But both are exactly the kind of finite exception the Church-Frobenius descent grammar predicted.  One sample packet opens at plain `q=33`; the other opens at plain `q=31`.  Both also have Bprime(any runner) certificates and positive exact safe measure.  When the obligation set is enlarged from Q27 to

```text
Q27 union {2,3,...,41},
```

both exceptional deletion addresses become infeasible.

So the corrected theorem is better:

```text
carry-window row retaining >=8 core speeds
=> Q27 witness or plain-shell witness q<=41.
```

This matters because it changes the proof portal.  We no longer need to prove that four deletions force Q27.  We need to prove that any row with no Q27 and no small plain `q<=41` witness either deletes at least five core speeds, exits the carry window, or descends/opens through a named channel.

Tournament Analysis lesson: the chosen vertices were proof obligations, not runners.  That was still right, but Q27 was one layer too narrow.  The retained channel must include the missing plain shells `31` and `33` (uniformly shells through `41`) plus owner/Bprime and positive-measure diagnostics.  Scalar Q27 did not fail randomly; it pointed to the next twist layer.
