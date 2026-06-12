# LRC14, the Pisano quotient, and the fibered band ladder

Source: codex-2026-06-11-P6.  Companions: THM-491, THM-492, HYP-2438,
HYP-2443, HYP-2444, `04-computation/lrc14_pisano_band_ladder_codex.py`.

## What Changed

The previous S7 result said the single-stranger family

```text
S(r) = 7*{1,...,12} union {r}
```

has five old-toolkit evaders. They block all plain band-1 shells `q <= 27`
and the B' dodge restricted to multiples of 14, then appear at the next
plain band rung `q in {40,41}`.

The new run changes the emphasis. The plain band-2 witnesses are real, but
they are not the minimal proof object. The fibered lattice

```text
Q27 = {d*m : d | 14, m <= 27}
```

already covers all 936 primitive rows. Even the two rows whose first plain
witness is `q=41` are caught at `q=91 = 7*13`.

So the corrected sentence is:

```text
plain shells need band 2, but the fibered band-1 lattice already closes the
one-stranger family.
```

That is a better shape for a proof of HYP-2438. It says the divisor fibers are
not an optional enhancement to the shell ladder; they are the ladder.

## The Pisano Class Is the Residual's Name

Modulo 27, the unit quotient `(Z/27)^*/+-` has 9 classes. This is the same
number as `pi(27)/pi(3) = 9`, and multiplication by 2 runs through a single
9-cycle:

```text
1, 2, 4, 8, 11, 5, 10, 7, 13.
```

The core `7*{1,...,12}` hits every class except `10`. Its q=27 witness is
the inverse class `8`. Thus a stranger can kill the last shell-27 unit twist
only by being `0 mod 27` or by occupying the missing class `+-10`.

The exact finite family adds one more coordinate: every plain-shell blocker
also has `r mod 13 = 0`. The residual is therefore a product obstruction:

```text
missing Pisano class at 27  x  killed 13-clock.
```

This is what was hiding inside the old evader list
`{611, 702, 793, 962, 1053}`.

## Why q=91 Matters

The old mental picture was:

```text
27 fails, so try 40 or 41.
```

But the sharper picture is:

```text
27 fails on the plain shell; the 7-fiber over the 13-clock gives 91.
```

That explains why two rows with first plain witness `q=41` are still Q27 rows.
The witness is not "higher band" so much as "same m-horizon, different fiber."

This is good news for HYP-2438. A proof can try to show that blockers must
spend independent resources across:

- the shell-27 antipodal unit quotient;
- the 13-clock;
- the divisor fibers `d in {1,2,7,14}`;
- the B' width gaps.

The one-stranger family spends exactly two resources. Multi-stranger rows are
the live gap because they can spend resources in parallel, but the resource
ledger is now visible.

## Tournament Analysis Choice

I used proof gates as vertices, not runners or arcs:

```text
Q41, band2, Q27, B'(any), band1, B'(mult14), shell27-residue, 13-clock.
```

This quotient is intentionally lossy. It preserves exact S(r) coverage and the
residue mechanism. It destroys endpoint-owner geometry and all nonlocal
multi-stranger pressure, so it should not be mistaken for a proof of C'(14).
It is a map of the proof obligations.

The resulting tournament is transitive. That is not surprising: one observable
was deliberately chosen to rank proof gates by coverage. A future nontransitive
tournament should compare two incompatible observables, for example:

```text
coverage strength  vs  formalizability  vs  multi-stranger robustness.
```

That would reveal tradeoffs rather than just rank the current gates.

## Next Concrete Move

The next useful computation is no longer "does Q41 help the single-stranger
family?" It does not add coverage there. The next computation should target
two-stranger rows, for example

```text
7*{1,...,11} union {r,s},
```

and record the first actual simultaneous blockers of Q27 if they exist. The
right output is a resource vector:

```text
(blocked shell-27 classes, blocked low clocks, blocked divisor fibers,
 B' target type, first surviving q).
```

If no two-stranger row can block Q27 and B'(any), HYP-2438's height-bound
lemma becomes much more plausible. If one exists, its resource vector is the
next theorem-shaped obstruction.

Mini-scout result: the first obvious two-stranger attempt is not hard. Pairing
the eight one-stranger plain-shell blockers

```text
260, 351, 442, 611, 702, 793, 962, 1053
```

over `7*{1,...,11}` gives 28 primitive pairs, and all 28 have first Q27
witness `q=12`. The reason is instructive: deleting the core runner `84`
opens a low divisor clock. So the next search should not merely combine
hard residues; it must maintain low-clock coverage while spending shell-27
classes.
