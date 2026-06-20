# LRC14 Cube-Root Order Filter

codex-2026-06-20-S52

## Summary

The cube-root hint gives an exact diagnostic, not a standalone proof.  For a
far triple, label the simultaneous-peel residual packets

```text
A,B,C = one-far
D,E,F = two-far
G     = three-far.
```

The actual cap correction is `A+B+C+D+E+F+G`.  The recursion
`A+B+C-D-E-F+G` is the pair-tax shadow `H(1)-2(D+E+F)`.

The cube-root refinement keeps the Eisenstein cyclic modes

```text
S_omega = A + omega B + omega^2 C
P_omega = D + omega E + omega^2 F.
```

All of this is exact over `Q(omega)`, represented by `a+b*omega` with norm
`a^2-ab+b^2`.

## Finding

The lenses separate.  In the `(15,16,17)` all-core bank, actual residual is
positive in `2821/3003` rows, but the pair-tax shadow is negative in
`1753/3003` rows.  Actual residual, pair-tax shadow, Eisenstein imbalance, and
direct `p0` choose different leaders.

That means the user's `A+B+C-D-E-F+G` recursion is useful because it exposes the
two-far layer, but it should not replace the actual residual.  The proof should
carry a three-coordinate packet:

```text
actual residual,
pair-tax shadow,
Eisenstein imbalance.
```

## Next Move

Use this packet as a finite-atlas coordinate for the rank-one triple relation
`u-2v+w=0`.  The scalar rank does not determine phase; S51 already showed
scaled AP-relation triples flip sign words with the shift.  S52 adds that the
cube-root imbalance and pair-tax shadow isolate different parts of that phase.

No LRC(14) proof is claimed.
