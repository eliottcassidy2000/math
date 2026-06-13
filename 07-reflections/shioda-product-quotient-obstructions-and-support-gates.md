# Shioda Product-Quotient Obstructions and Support Gates

Source: codex-2026-06-12.  Companions: arXiv:2508.14876
(`https://arxiv.org/abs/2508.14876`), HYP-2445, OPEN-Q-069,
`04-computation/product_quotient_support_gate_atlas_codex.py`.

## The Paper's Lesson

Church's construction is a clean warning against scalar optimism.  The
surfaces are arranged so that the scalar cohomological shadow is as friendly as
one could want: simply-connected and Shioda supersingular at infinitely many
good reductions.  Yet they are not unirational.

The obstruction comes from a channel that the scalar shadow forgets:
diagonal symmetric differential forms on all asymmetric partial Frobenius
twists.  Rational and elliptic curves cannot wander freely.  They are forced
into finite exceptional families, or they descend along partial Frobenius while
one projection degree drops.  The proof is not "supersingular implies X"; it is
"supersingular passes, but a diagonal side channel still sees the curves."

That is exactly the support-gate grammar the repo has been finding from the
other side.

## The 1092 Hinge

The genus-14 Hurwitz curve in the paper has automorphism group
`PSL2(F_13)` of order

```text
1092 = 84*(14-1) = 13*84 = 2^2*3*7*13.
```

This is the same `1092` used in the LRC14 one-stranger cutoff.  The subgroup
indices sharpen the echo:

```text
D6, A4 index 91 = C(14,2) = LRC14 q=91 fibered rescue.
D7     index 78 = C(13,2) = [72,36,16] lambda_5.
```

I do not want to oversell that.  It is numerology until it pays rent.  But it
is unusually organized numerology: genus 14, the Hurwitz bound, the `13*84`
LRC cutoff, the `91` fiber, and the `78` design parameter all sit on the same
small ledger.

## LRC14 Transfer

The partial-Frobenius descent proof suggests a way to phrase the missing LRC14
resource theorem.

For Church's surfaces:

```text
curve not exceptional -> tangent to a standard foliation
-> pull back from a partial Frobenius twist
-> projection degree drops.
```

For HYP-2444, the hoped-for analogue is:

```text
configuration blocking Q27 and Bprime(any)
-> a runner/support set pays for shell-27 class, 13-clock, and divisor fibers
-> either some resource strictly drops or a named Bprime/owner-deletion
   exceptional route opens.
```

This gives a better proof target than "try bigger q."  It asks for a monotone
resource under the same `Q27={d*m:d|14,m<=27}` horizon.

## Code Transfer

For `[72,36,16]`, the scalar lesson is already established.  The Type II
Gleason enumerator is healthy; the problem is support realization.  The
surface paper gives another mature example where scalar positivity/cohomology
does not decide the geometric question.

The `78` echo is the one to test carefully:

```text
index_{PSL2(F13)}(D7) = 78
lambda_5 of the length-72 minimum design = 78.
```

The right experiment is not to force `PSL2(F13)` as a code automorphism.  That
would fight the near-rigidity expectations for an extremal length-72 code.  The
right experiment is incidence arithmetic: can a coset/design ledger with a
`78`-sized local layer generate useful obstruction tests for the minimum-word
support or for the order-5 `F_16` completion problem?

## Working Rule

When a future route has a tempting scalar invariant, ask immediately:

```text
what side channel survives the quotient?
what twists/deletions/fibers preserve it?
what descent measure strictly drops?
what exceptional types are finite and named?
```

If the answers are absent, the invariant is probably only a shadow.  Church's
paper is a valuable new example because it has all four answers in a serious
geometric setting.
