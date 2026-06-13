# Edge Variable

**Symbol:** `s_e`

**Type:** centered tournament edge coordinate.

For an unordered pair `e={i,j}` with `i<j`, let

```text
A_e = 1 if i -> j, and A_e = 0 if j -> i
s_e = A_e - 1/2
```

Thus `s_e` is always one of `{-1/2,+1/2}`.

## Pair-First Meaning

The edge variable records the repo's basic pair-cell:

```text
two endpoints + one relation bit.
```

Vertex scores, Hamiltonian path counts, Fourier coefficients, good-cut
responses, and many Tournament Analysis features are projections of this
complete edge-coordinate board.

This is the finite tournament analogue of the prime-pair coordinate

```text
{a,b} -> (sum=a+b, gap=b-a, local residue labels).
```

In that lens, twin primes are not a vertex phenomenon.  They are the persistence
of the fixed edge-label row `gap=2` after the prime sieve.  A tournament edge
variable is similarly not a property of either endpoint alone; it is the
connection between them.

## See Also

`W-polynomial.md`, `fourier-coefficients.md`, HYP-1965,
`04-computation/pair_first_twin_prime_lens_s502.py`.
