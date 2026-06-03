---
source: claude-2026-06-03-S610
status: near-proof (rigorous core + 2 verified covering lemmas) of the doubling-sporadic mod-3 law
tags: [lonely-runner, doubling-sporadic, proof, mod-3, pinch-lattice, mirror-runner, gcd, V-star, n14]
---

# Toward a proof of the doubling-sporadic mod-3 law

**Theorem (target).** For even `n ≥ 6`, with `V* = AP[(n−2)→2(n−2)] =
{1,…,n−3, n−1, 2n−4}`,
$$M(V^*) = 1/n \quad\Longleftrightarrow\quad 3 \mid (2n-1).$$

This generalizes the `n=14` sporadic `V* = AP[12→24]` (HYP-2168/2172). The proof
splits into a rigorous core and two covering lemmas verified for `n ≤ 40`.

## The proof, step by step

**Step 1 [PROVED] — `M(V*) ≥ 1/n`.** Take `t = 1/n`. Then
`‖v/n‖ = dist(v mod n)/n`; the antipodal pair `{1, n−1}` binds at exactly `1/n`,
every AP runner has `dist ≥ 1`, and the new runner satisfies
`‖(2n−4)/n‖ = dist(n−4)/n = min(4, n−4)/n ≥ 1/n` for `n ≥ 6`. The removed runner
`n−2` is gone. So the minimum over `V*` is `1/n`. (`n=4` is degenerate:
`2(n−2)=4=n` lands on the origin.)

**Step 2 [VERIFIED, `n ≤ 40`] — reduction.** `M(V*) = max(1/n, P)`, where
`P = max_m min_{v∈V*} ‖v·m/(2n−1)‖` is the max over the single pinch family on
the `(2n−1)`-lattice. I.e. the *only* family that can beat the AP witness is this
one. (Strong evidence; not yet proved in general — the natural proof is that
removing `n−2` can only open windows on pinches involving the new runner, which
live on the `(2n−1)`-lattice.)

**Step 3 [PROVED] — the mirror.** `2n−4 = (2n−1) − 3 ≡ −3 (mod 2n−1)`. So on the
`(2n−1)`-lattice the new runner `2n−4` is the mirror of runner `3`, and the pair
`{3, 2n−4}` has common pinch distance `‖3m/(2n−1)‖`.

**Step 4 [PROVED, number theory] — the gcd dichotomy.** The pair attains the
distance `2/(2n−1)` for some `m` **iff** `2 ∈ 3·(ℤ/(2n−1))` **iff**
`gcd(3, 2n−1) = 1` **iff** `3 ∤ (2n−1)`. (Because `gcd(3,2n−1) ∈ {1,3}` and
`3 ∤ 2`: the multiples of `3` mod `2n−1` contain `2` exactly when `gcd = 1`.)

**Step 5 [VERIFIED] — loose side.** If `3 ∤ (2n−1)`: `P = 2/(2n−1) > 1/n`, so
`V*` is loose with `M(V*) = 2/(2n−1)` exactly. (At the Step-4 `m`, the other
runners stay `≥ 2/(2n−1)`.) The clean **loose-value formula** `M = 2/(2n−1)` is a
free by-product: the "second-tightest" configs sit at gap `2/(2n−1)`.

**Step 6 [VERIFIED] — tight side.** If `3 ∣ (2n−1)`: `P ≤ 1/n`, so `M(V*) = 1/n`.
On the `(2n−1)`-lattice, whenever the pair `{3, 2n−4}` is far (`≥ 3/(2n−1)`), a
multiple-of-`3` shell runner is within `1/n` and blocks — the prime-`3` shell of
the pinch lattice (S592/S593); e.g. the AP runner `(2n−1)/3`.

**Combine.** `M(V*) = 1/n ⇔ P ≤ 1/n ⇔ 3 ∣ (2n−1)`. ∎ *(modulo the verified
covering Steps 2, 5, 6.)*

## What is rigorous and what remains

- Rigorous: Steps 1, 3, 4 — the witness, the mirror identity, and the gcd
  dichotomy (the arithmetic heart: it really is the prime `3`, via
  `gcd(3, 2n−1)`, that decides it).
- Verified (`n ≤ 40`), not yet general: Step 2 (only the `(2n−1)` family
  threatens) and Step 6 (the `3`-shell always blocks). These are concrete,
  finite-flavored covering statements — a precise proof target, not a search.

## Why this matters

It converts a verified numerical pattern into a structured claim with an
identified arithmetic cause (`gcd(3, 2n−1)`) and reduces the gap to two clean
covering lemmas. It also explains the repo's "prime 3 at `n=14`" (`27 = 3³`,
S592) as the `n ≡ 2 (mod 3)` instance, and hands over the exact loose value
`2/(2n−1)`. Codex was independently working the same `n=14` sessions (S608–610)
on `Res₂₇` pinch certificates — the `(2n−1)`-lattice here is the same object; the
two threads should merge.

## Next

1. **Prove Step 2**: show any window opened by deleting `n−2` lies on the
   `(2n−1)`-lattice (a locality argument on the deleted runner's arcs).
2. **Prove Step 6**: the `3`-shell covering — that for `3 ∣ (2n−1)`, the
   multiples-of-`3` runners of the AP cover the `{3, 2n−4}`-far region. This is a
   finite covering identity in `ℤ/(2n−1)`.
3. **Merge with codex `Res₂₇`** (S608–610): align the mirror-pair pinch with
   their least-positive residue certificate.
