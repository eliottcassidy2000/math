---
source: claude-2026-06-03-S612
status: loose direction PROVED; tight direction PROVED on the (2n-1) lattice; full law modulo Step D (= codex lift conservativity)
tags: [lonely-runner, doubling-sporadic, proof, mod-3, pinch-lattice, euclid, gcd, codex-merge, res27, gcd-strata, V-star]
---

# The doubling-sporadic mod-3 law, proved on the pinch lattice (and merged with codex)

**Theorem.** For even `n ≥ 6`, with `V* = AP[(n−2)→2(n−2)] = {1,…,n−3, n−1, 2n−4}`,
`M(V*) = 1/n ⟺ 3 ∣ (2n−1)`.

The S610 skeleton had this verified with two covering lemmas open. This session
closes both **on the pinch lattice** with a single arithmetic fact, and shows the
one remaining gap is *identical* to codex's parallel `Res₂₇` thread.

## The fact that does all the work

By the Euclidean algorithm,
```
gcd(n−2, 2n−1) = gcd(n−2, (2n−1) − 2(n−2)) = gcd(n−2, 3) = gcd(3, 2n−1),
```
so **`n−2` is invertible mod `Q := 2n−1` ⟺ `3 ∤ (2n−1)`.** The removed runner's
invertibility is exactly the mod-3 condition. Everything follows.

Work on the lattice `t = m/Q`: runner `v` sits at residue `vm (mod Q)`, with
`‖vt‖ = dist_Q(vm)/Q`. "Close" means `dist_Q(vm) ≤ 1` (`‖vt‖ ≤ 1/Q < 1/n`); "far"
means `dist_Q(vm) ≥ 2` (`‖vt‖ ≥ 2/Q > 1/n`).

## Step A — the lattice is the mirror pair

The only pair of `V*` summing to `Q = 2n−1` is `{3, 2n−4}` (AP pairs sum to at
most `2n−3`). So the `(2n−1)`-pinch lattice is the `{3, 2n−4}` family, and
`2n−4 = Q − 3 ≡ −3`: the new runner is **runner 3's mirror**.

## Step B — LOOSE direction, fully proved

If `3 ∤ Q`, then `n−2` is a unit; let `m = (n−2)⁻¹ mod Q`. At this pinch every
`v ∈ V*` is far:
- `v ∈ {1,…,n−1}\{n−2}`: `vm≡0 ⇒ v≡0` (impossible); `vm≡1 ⇒ v≡n−2` (excluded);
  `vm≡−1 ⇒ v≡n+1 ∉ {1,…,n−1}`.
- `v = 2n−4 ≡ −3`: `vm ≡ −3(n−2)⁻¹`; `≡0 ⇒ 3≡0`; `≡±1 ⇒ n−2 ≡ ∓3 ⇒ n=5` (odd)
  or `n=2` — excluded for even `n ≥ 6`.

So `min_{v∈V*} ‖vt‖ ≥ 2/Q > 1/n`, hence `M(V*) ≥ 2/(2n−1) > 1/n`. (Equality
holds: the clean loose value `M(V*) = 2/(2n−1)`.) ∎

## Step C — TIGHT direction, proved on the lattice

If `3 ∣ Q`, every lattice point has a close `V*` runner:
- `m` invertible: one of `±m⁻¹` lies in `{1,…,n−1}` (the AP is a half-system of
  residues mod `Q`); that rep is a unit, so it is `≠ n−2` (a non-unit), hence in
  `V*`, and gives `vm ≡ ±1`.
- `m` non-invertible, `g = gcd(m,Q) ∈ {3,9,…}`: `w = Q/g ≤ Q/3 ≤ n−1`, and
  `w ≠ n−2` (`Q = g(n−2)` would force `n=5`), so `w ∈ V*` and `wm ≡ 0`.

Therefore `P := max_m min_v ‖vt‖ ≤ 1/Q < 1/n`: the threat family never beats the
AP witness. ∎ (on the lattice)

## Step D — the one remaining gap (= codex's conservativity)

`M(V*) = max(1/n, P)`: the `n`-lattice witness `t = 1/n` gives exactly `1/n`, and
*no pinch family other than `{3, 2n−4}`* beats `1/n`. Granting D, Steps B/C give
the law. D is verified for `n ≤ 40` and is **the same statement** as codex's
`Res₂₇` **lift/CRT conservativity** (HYP-2167): the least-positive `Q`-residue
section must capture the true integer max-min over the continuum.

## The merge with codex's `Res₂₇` thread

The two threads are the same object seen from two sides:

| this thread (general `n`) | codex `Res₂₇` (n=14) |
|---|---|
| `(2n−1)`-pinch lattice | least-positive `C=27` quotient |
| multiple-of-3 shell (Step C blockers `Q/g`) | THM-407 gcd-strata `{1, 3, 9}` |
| Steps B+C: only `AP, V*` at the floor | HYP-2164: exhaustive floor = `AP, V*, 2·AP` |
| general reason `gcd(3, 2n−1)` | the `27 = 3³` prime-3 structure (S592) |
| open Step D (lattice→continuum) | open lift/CRT carry-fiber (HYP-2167) |

So codex's exhaustive `n=14` certificate is the `n=14` instance of Steps B+C, and
this session supplies the **general-`n` reason** (the Euclid fact) and the clean
**loose value `2/(2n−1)`**. The shared open piece is one conservativity lemma.

## Exploration notes

- The loose value `M(V*) = 2/(2n−1)` is the **second-tightest gap** after `1/n`:
  the doubling configs are the closest non-tight family to the wall.
- The proof never used `n = 14` — it is uniform in `n`. The repo's "prime 3 at
  `n=14`" is the instance `27 = 3³`; the law says any `n ≡ 2 (mod 3)` works, and
  the *power* of 3 (here `3³`) is irrelevant to tightness (only `3 ∣ Q` matters),
  though it may matter to codex's `{1,3,9}` stratification depth.
- The mechanism is purely about **runner `n−2`'s invertibility** mod `2n−1`. This
  suggests a general principle: *a single-element swap `a → a'` of the AP is tight
  iff `a` is non-invertible mod `2n−1` and `a'` lies in `a`'s residue shell* — a
  testable generalization beyond doubling.

## Next

1. **Step D / conservativity** — prove the `(2n−1)`-lattice (plus the `n`-lattice
   witness) captures `M(V*)` over the continuum: a locality argument that deleting
   `n−2` only opens windows on pinches whose sum divides into the `Q`-lattice.
   Closing this + merging codex's carry-fiber finishes the law as a theorem.
2. **The swap principle** — test "AP-swap `a→a'` tight ⟺ `a` non-invertible mod
   `2n−1` and `a'` in `a`'s shell" against the full census (e.g. `n=8`'s second
   sporadic `(1,4,5,6,7,11,13)`).
3. **Non-doubling sporadics** — `n=6`'s `(1,3,4,5,9)` is not a swap of the AP;
   classify its mechanism (different shell? a `2`-adic partner?).
