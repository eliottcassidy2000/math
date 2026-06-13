# Unordered Bucket Balance as an Orbit Theorem

**Session:** kind-pasteur-2026-05-30-S2
**Status:** reflection after Lean orbit-parity formalization
**Related:** THM-346, THM-348, THM-350, HYP-1775, HYP-1778

The S1 half-line formalization made THM-346 look like two lemmas:

```text
finite conservation + fixed-point-free involution = unordered balance.
```

The S2 Lean layer sharpened that further, and the new orbit lemma closes the
abstract gap. Lean now proves that a finite fixed-point-free involution has
even cardinality, then applies it to internal half-lines. The strengthened
theorem gives

```text
2*internalLineCount + crossHalf = |fiber|*|moves|.
```

So the genuinely missing object is no longer a bucket theorem or an orbit
theorem. It is the Boolean-mask action theorem for the tiling cube.

For each internal half-line `(x,u)`, the partner is

```text
(x,u) -> (step(u,x),u).
```

If `step(u)` is an involution, this partner map comes back in two steps. If
`step(u)` has no fixed points, no internal half-line can be its own partner.
Lean now proves both local facts and the global finite-cardinality statement:

```text
fixed-point-free involution on a finite set => even cardinality.
```

That lemma is the conceptual place where unordered lines enter. The unordered
edge is not primitive; it is the two-element orbit of the oriented partner
action.

This gives a useful way to think about quotient transport. A bucket boundary
is measured by half-lines because half-lines are what the source bucket can
see. Internal geometry is measured by orbits because both endpoints live in
the same bucket. Cross geometry is one-sided from a fixed bucket, so it is not
halved. The formula is exactly the combination of those two viewpoints.

Mathematically, this suggests a hierarchy worth reusing:

1. A finite bucket map gives a row conservation law.
2. A move semigroup with reversible generators gives partner maps.
3. Fixed-point-free reversible generators make internal mass even.
4. Quotient-specific geometry then explains where the escaping mass lands.

The merged tiling case occupies level 4 with spine/ribs/sea, projection
defect, and good-cut height changes. But levels 1-3 are generic. That means
future exact computations should be organized around the same split: prove
the row law once, prove the orbit parity once, then let the quotient geometry
carry the interesting combinatorics.

The next Lean move should therefore be cube-specific:

```text
step(u,x)=x xor u is a fixed-point-free involution for u != 0
```

Once that exists, THM-350 should collapse the remaining formal gap to full
THM-346.
