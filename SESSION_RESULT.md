# Session Result

## Task Chosen

I chose one small Hamiltonian-path sanity check: re-verify the
near-transitive single-arc flip formula for the transitive tournament
`T_n`.

The checked claim is THM-250:

```text
H(T_n^(a,b)) = 2^(b-a-1) + 1,  if b-a >= 2
H(T_n^(a,b)) = 1,              if b-a = 1
```

where `T_n^(a,b)` is obtained from the transitive tournament by reversing the
single arc `a -> b`.

## What I Did

I ran a transient exact Held-Karp dynamic program for Hamiltonian path counts.
For every `n = 3..11` and every single flipped arc `0 <= a < b < n`, I computed
`H(T_n^(a,b))` and compared it with the formula above.  No floating-point
arithmetic was used, and no retained script was added.

## Concrete Result

There were zero mismatches across all checked flips for `n = 3..11`.

The values depend only on the distance `d = b-a`:

```text
d=1  -> H=1
d=2  -> H=3
d=3  -> H=5
d=4  -> H=9
d=5  -> H=17
d=6  -> H=33
d=7  -> H=65
d=8  -> H=129
d=9  -> H=257
d=10 -> H=513
```

This extends the computational verification range stated in
`01-canon/theorems/THM-250-flip-formula.md` from `n <= 9` to `n <= 11` for
this independent DP check.

One minor documentation issue surfaced: the theorem statement is consistent
with the computation, but its corollary display includes `2` in the listed set
of single-flip values.  Under the theorem's own distance convention, no checked
single flip gives `H=2`; the set through `n` should be
`{1, 3, 5, 9, 17, ..., 2^(n-2)+1}`.

## Confidence Note

Confidence is high for this narrow verification.  The computation used exact
integer counts and independently checked every single flip through `n=11`.
This is not a new proof of THM-250; it is a finite audit plus a small corollary
typo detection.
