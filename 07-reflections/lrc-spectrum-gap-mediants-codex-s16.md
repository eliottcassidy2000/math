# LRC Spectrum Gap Mediants - Codex S16

The user's doubled-top family is real and useful:

```text
{1,...,k-1,2k}  ->  M = 2/(2k+1).
```

That makes the reciprocal gap exactly `(k+1)(2k+1)`, not a heuristic lift-depth.
The exact envelope script verifies it cleanly for `k=2..30`.

The more interesting find is that the next constant is hiding in an AP defect,
not in the top-replacement family.  The family

```text
{1,...,k}\{k-1} union {3(k-1)}
```

splits by residue:

```text
k = 1 mod 30          -> AP-tight again
k = 7,13,19,25 mod30 -> M = 3/(3k+2).
```

So the gap can be

```text
1 / ((k+1)(3k+2)),
```

which beats the doubled-top constant but stays at `Theta(1/k^2)`.  The `k=31`
row was the sanity slap: it looked like the same `1 mod 6` family, but it fell
all the way back to the AP floor.  The mod-`30` address is not decorative.

The immediate next route is a multiplier-ladder classifier:

```text
A_{k,r} = {1,...,k}\{k-1} union {r(k-1)}.
```

In the stored probe, `r=4` is no longer isolated: on the `k=1 mod30` class where
`r=3` goes tight, `r=4` gives

```text
M = 4/(4k+3)
```

through `k=181`.  That suggests a ladder of constants may exist, but it is
sparse and residue-governed.  To threaten the lower bound order, we would need
constants `r` growing with `k`, not just fixed `2,3,4`.  The next computation
should scan by residue classes in `(k,r)` and try to recognize formulas of the
shape `r/(rk+c)`.

This is also the cleanest Tournament Analysis lesson from the session.  Runner
vertices are too local; pair-crossing denominators are better; but for the
lower-bound question the right tournament vertices are proof routes/families.
They preserve exactly the scalar we care about, the gap above `1/(k+1)`, while
keeping just enough address data to guess the next recurrence.
