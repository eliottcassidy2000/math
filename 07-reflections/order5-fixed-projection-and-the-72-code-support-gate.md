# Order-5 fixed projection and the 72-code support gate

The earlier `[72,36,16]` work kept saying "the scalar gate passes; the support
gate is the problem." This session found a small support gate that is sharp
enough to be useful.

For an order-5 automorphism of the known possible type `5-(14,2)`, the fixed
subcode projects to length `16`: fourteen orbit coordinates plus two fixed
coordinates. Since `5 == 1 (mod 4)`, doubly-evenness survives the projection.
So the fixed projection is one of the two Type II `[16,8]` codes, `e8+e8` or
`d16+`.

The two actual fixed coordinates become two marked projected coordinates. A
projected tetrad containing both marks lifts to original weight `12`, because
it uses two 5-cycle orbit coordinates and both fixed coordinates. That is the
whole obstruction.

The exact computation says:

- `d16+` is dead: every coordinate pair is covered by a tetrad.
- `e8+e8` survives exactly for split marked pairs, one mark in each `e8` block.

So the order-5 branch, if it exists, is not just "some symmetric case." It is
forced into a very specific shape:

```text
one fixed point + seven 5-cycles + Fano-plane tetrads,
one fixed point + seven 5-cycles + Fano-plane tetrads.
```

The fixed minimum words are then exactly `14`, leaving `49967` moving orbits of
minimum words. This is satisfying because `249849 == 4 (mod 5)` and `14 == 4
(mod 5)`, but it is also constraining: the residual `5-(72,16,78)` design has
to be carried almost entirely by nonfixed orbits.

The next real problem is the `F_16` component. The group algebra factor
`x^4+x^3+x^2+x+1` asks for a Hermitian self-dual length-14 code over `F_16`,
but now with an `e8+e8` split-heptad fixed boundary and exact residual design
counts. That is a finite support-completion problem, not a modular-form problem.

## Assumption challenge

The Tournament Analysis vertices were not runners or arcs. I considered:
automorphism cases, projected code types, marked pairs, tetrads, Fano heptads,
residue classes, nonfixed Fourier modes, and proof obligations. The useful
quotient here was marked-pair branch classes. It preserves the fixed-code
low-weight leak and the `A_16 mod 5` ledger, but it destroys the nonfixed
eigenspace and coset-glue information. The resulting tournament is transitive;
that is a warning that the next quotient needs at least one nonfixed-glue
observable.

## Source touchpoints

- Feulner-Nebe: automorphism restrictions leave `Aut(C)` either order `5` or
  dividing `24`: https://arxiv.org/pdf/1110.6012
- Bouyuklieva/O'Brien/Willems and subsequent summaries give the prime-order
  cycle types, including `5-(14,2)`: https://www.researchgate.net/publication/3084801_On_the_Automorphism_Group_of_a_Doubly-Even_723616_Code
- Borello-Dalla Volta-Nebe exclude larger subgroup patterns and sharpen the
  "near-rigid" context: https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/S3A4D8.pdf

## Files

- `04-computation/order5_fixed_projection_72_codex.py`
- `05-knowledge/results/order5_fixed_projection_72_codex.out`
- `05-knowledge/hypotheses/HYP-2439-order5-fixed-projection-marked-type-II-16-gate.md`
- `05-knowledge/hypotheses/HYP-2440-d16-plus-is-excluded-from-order5-fixed-projection.md`
- `05-knowledge/hypotheses/HYP-2441-order5-heptad-split-and-minimum-design-ledger.md`
