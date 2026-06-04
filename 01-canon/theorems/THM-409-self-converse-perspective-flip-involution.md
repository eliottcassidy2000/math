---
id: THM-409
title: Self-converse tournaments have a canonical perspective-flip involution
status: PROVED
source: codex-2026-06-04-S629
depends_on:
  - HYP-2121
  - THM-408
related:
  - HYP-2204
  - HYP-2205
  - HYP-2133
---

# THM-409 - Self-Converse Perspective Flip Involution

## Statement

Let `T` be a tournament, let `Aut(T)` be its automorphism group, and let
`Anti(T)` be the set of anti-automorphisms, meaning permutations `sigma` of the
vertex set satisfying

```text
x -> y in T  iff  sigma(y) -> sigma(x) in T.
```

Equivalently, `sigma` is an isomorphism from `T` to its converse `T^op`.

If `T` is self-converse, so that `Anti(T)` is nonempty, then:

1. `Anti(T)` is a left and right coset of `Aut(T)`.
2. Every `sigma in Anti(T)` has `sigma^2 in Aut(T)`.
3. `sigma` induces a well-defined map on rooted perspectives, i.e. on the
   `Aut(T)`-orbits of vertices:

   ```text
   [v] -> [sigma(v)].
   ```

4. This induced map is independent of the chosen `sigma in Anti(T)`.
5. The induced map on rooted perspectives is an involution.

Thus edge reversal need not canonically swap individual vertices back and
forth.  What is canonical is the back-and-forth involution on rooted
perspectives.

## Proof

Choose `sigma in Anti(T)`.  If `a in Aut(T)`, then `a sigma` and `sigma a` are
anti-automorphisms, because composing an arc-preserving permutation with an
arc-reversing permutation reverses arcs.

Conversely, if `sigma,tau in Anti(T)`, then `tau^{-1} sigma` preserves arcs:
both `sigma` and `tau` reverse arcs, so their composition reverses twice.  Thus
`tau^{-1} sigma in Aut(T)`, and `sigma` and `tau` lie in the same `Aut(T)`
coset.  This proves the coset statement.

The square `sigma^2` is a composition of two anti-automorphisms, hence an
automorphism.

Now let `[v]` denote the `Aut(T)`-orbit of `v`.  If `w` lies in the same orbit
as `v`, then `w=a(v)` for some `a in Aut(T)`.  Since anti-conjugation preserves
automorphisms,

```text
sigma a sigma^{-1} in Aut(T),
```

and therefore

```text
sigma(w) = sigma a(v) = (sigma a sigma^{-1})(sigma(v))
```

lies in the same `Aut(T)`-orbit as `sigma(v)`.  So `[v] -> [sigma(v)]` is
well-defined.

If `tau` is another anti-automorphism, then `tau=a sigma` for some
`a in Aut(T)`.  Hence `tau(v)=a(sigma(v))`, which lies in the same rooted
perspective as `sigma(v)`.  So the map on perspectives is independent of
which anti-automorphism was chosen.

Finally,

```text
[v] -> [sigma(v)] -> [sigma^2(v)] = [v],
```

because `sigma^2 in Aut(T)`.  Hence the perspective flip is an involution.
QED.

## Programmatic Form

To compute the object:

1. Enumerate `Aut(T)`.
2. Enumerate any one `sigma in Anti(T)`.
3. Compute the vertex orbits of `Aut(T)`.
4. Map each orbit `O` to the orbit containing `sigma(v)` for any `v in O`.

The theorem proves that step 4 is independent of both the representative
`v in O` and the chosen anti-automorphism `sigma`.

## Verification

`04-computation/sc_perspective_flip_cyclotomic_s629.py` verifies this theorem
exhaustively for tournament isomorphism classes through `n=6`.  It also shows
why the perspective formulation is necessary: at `n=6`, some
anti-automorphisms have vertex cycle type `(6,)`, so the individual vertex
motion is not always an involution.  The induced rooted-perspective map is an
involution in every self-converse class checked.

**Artifacts:** `04-computation/sc_perspective_flip_cyclotomic_s629.py`,
`05-knowledge/results/sc_perspective_flip_cyclotomic_s629.out`,
`07-reflections/sc-perspective-flip-cyclotomic-carriers-s629.md`.
