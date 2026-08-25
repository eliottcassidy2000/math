# A fixed-gauge obstruction can be the derivative of a mixed carrier

**Date:** 2026-08-25
**Status:** research reflection; THM-4058 and THM-4060 are the truth sources
for mathematical claims

## Portfolio and inheritance

- **Anchor:** decide whether THM-4058's all-`m` fixed-`a` obstruction survives
  when the first target output is allowed to move.
- **Niche:** replace a sequence of larger coupled matrices by a circulation
  identity on the exceptional three-branch triangle.
- **Wildcard:** factor a linearly constructed closed target two-form only
  after solving the branch restriction problem, instead of integrating the
  nonlinear pair equation term by term.

The closest proved mechanism was THM-4054's rank-`59` mixed tangent escape.
The hostile example was THM-4058 itself: the frozen-output cokernel has one
genuine class in every degree, so a claim of mixed survival needed an
all-degree carrier rather than another cutoff-five computation.  The corrected
near miss was MISTAKE-493: parity of the fold polynomial cannot be imported for
the exceptional non-even quartic.  The least-used sidecar was the fact that a
nonvanishing closed two-form in three formal variables can be factored after
its branch restrictions have been solved.

## Concept board

The live board had five objects:

1. the normalized fixed-output period `Lambda`;
2. the signed triangle circulation `Pi` before division by `delta A^5`;
3. first-output carriers `f` and their `A` derivatives on the branches;
4. the ambient closed-form complex in `(A,c,u)`;
5. formal presymplectic Darboux factorization.

The decisive comparison was between items 2 and 3.  If `f` is one common
target germ, differentiation of its circulation has no endpoint residue:
the value at each moving vertex occurs once with each sign.  Thus the
first-output part is not a new unrelated matrix block.  It is the derivative
of the same triangle period that obstructs the fixed-output primitive.

## What changed after the pull

The triangle has edge scale `A^5` and area scale `A^10`.  Therefore every
common target germ has normalized circulation in `A^5 K[[A]]`.  Applying the
first-output derivative changes this to `A^4 K[[A]]`.  This explains the
otherwise surprising rank table exactly:

```text
fixed output: one cokernel class in every retained degree;
mixed closed forms: only degrees 0,1,2,3 survive.
```

The powers `u^s` are triangular carriers.  Their leading responses never
vanish because they are nonzero rational multiples of the already audited
field element `rho`.  Hence every fixed-output monomial obstruction lies in
the mixed transgression range.

This also reorganized the nonlinear question.  Directly solving
`dF wedge dG` would introduce the quadratic term `df wedge dg` at every
stage.  The cleaner order is:

```text
branch density
  -> one ambient closed two-form
  -> formal Darboux factorization
  -> target reparametrization u=H(w).
```

The restriction problem is linear, and decomposability is deferred until a
nonvanishing closed form is already available.  That separation is the real
reason the dual-number escape can be integrated formally.

## Hostile controls and boundary

The finite probes were still useful, but only as hostiles:

- the complete mixed ranks have cokernel `min(N+1,4)` from the smallest
  cutoffs onward;
- the cutoff-five checkpoint is exactly `59/63`, agreeing with THM-4054;
- all four roots modulo `137` must show the same rank and response behavior;
- a target-degree cap can fail at the next source cutoff even though the
  complete jet bank succeeds, so no bounded polynomial ansatz may be confused
  with the formal theorem.

The result remains local and formal.  Darboux factorization supplies formal
target germs, not convergent or polynomial ones.  The rational coordinate
`a=e/(b+4)` still has a global denominator, other fibres are uncontrolled,
and no injectivity, properness, or behavior-at-infinity statement follows.

## Research-move assessment

The move “differentiate the obstruction period with respect to a newly freed
coordinate” is a plausible meta-pattern candidate.  It should not yet enter
`META-PATTERNS.md`: there is evidence from only this one thread.  A second
independent example should verify both its trigger (a moving closed cycle with
endpoint cancellation) and its counterindication (a boundary value that is not
the restriction of one common target germ).
