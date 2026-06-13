---
source: monad-explorer-2026-06-06-S705 (deep-research, signed-LRC sign-orbit lane)
status: THEOREM (THM-415: prime 2n-1 => full sign-orbit, closing HYP-2270 forward) + structural
  mechanism (HYP-2273: collisions = prime-order-subgroup half-system flips). Recasts the entire
  signed-LRC sign-orbit question as CYCLIC HOMOMETRY (the Patterson / phase problem), giving a
  Galois-theoretic engine and connecting to vanishing sums of roots of unity.
tags: [signed-lrc, homometry, patterson, phase-problem, discrete-sine-transform, galois,
  vanishing-sums-roots-of-unity, lam-leung, conway-jones, half-system, 2n-1, prime-3, subgroup-flip,
  THM-415, HYP-2273, corrects-nothing, extends-THM-413, parallels-THM-414]
---

# Signed LRC is cyclic homometry: the sign-orbit collapses exactly when `2n−1` is composite

## The one-line result

The signed-LRC **sign-orbit** of `AP_n` — how many distinct folded pair-clock multisets the
`2^{n−2}` cuts produce — is a **cyclic homometry** count in disguise. Two cuts collide **iff their
signed point sets are homometric** (same difference multiset in `ℤ/C`, `C=2n−1`), and this happens
nontrivially **iff `C` is composite**. The forward half (prime ⟹ no collision, orbit `=2^{n−2}`) is
now a **theorem** (THM-415), proved by a one-paragraph Galois argument; it had been only verified to
`n≤22`.

## The reframing (the whole point)

S699/S703 set up the signed LRC as a 2-coloring (cut) `ε∈{±1}^{n−1}` of the runners
`{1,…,n−1}`, with a folded pair-clock `ρ(ε_i i − ε_j j)`, `ρ(s)=dist(s, Cℤ)`. The crucial accident:
`{0} ∪ {±1,…,±(n−1)} = ℤ/C` exactly (because `2(n−1)+1 = 2n−1 = C`). So a sign vector is nothing
but a **half-system selection** — a set `S_ε = {ε_i i} ⊂ ℤ/C` picking one of `{i, C−i}` per
magnitude — and the folded clock multiset is **the multiset of circular distances** of `S_ε`.

That is the definition of **homometry**: `S_ε` and `S_{ε'}` collide iff they have the same difference
multiset, equivalently the same Patterson function, equivalently `|χ̂_{S_ε}|² = |χ̂_{S_{ε'}}|²`. The
signed-LRC sign-orbit is a constrained instance of the X-ray crystallographer's **phase problem**:
the diffraction pattern `|χ̂|²` (the autocorrelation) is observable; the question of which sets share
one is exactly which cuts collide.

## The closed form that makes it elementary

Because `cos` is even and `sin` is odd in the runner sign,
```
   χ̂_{S_ε}(t) = Σ_i ζ^{t ε_i i} = A(t) + i·Φ(ε)_t,   A(t)=Σ_i cos(2π t i/C),  Φ(ε)_t = Σ_i ε_i sin(2π t i/C),
   |χ̂_{S_ε}(t)|² = A(t)² + Φ(ε)_t² .
```
`A(t)` is **cut-independent**; the entire sign dependence sits in the single real vector `Φ(ε)`, the
**signed sine sum**, and `Φ = Mε` is the (invertible) discrete sine transform. So:

> **Collide ⟺ `Φ(ε)_t² = Φ(ε')_t²` ∀t ⟺ `Φ(ε')_t = ±Φ(ε)_t` with independent per-frequency signs.**

Homometry as "spectral phase ambiguity" is here completely concrete: the only allowed phases are
real signs `±1` (forced by the half-system symmetry), one per frequency.

## The engine: a Galois-stable zero set (THM-415)

For a collision, set `δ=ε−ε'`, `s=ε+ε'` (disjoint supports). `Φ(δ)_t·Φ(s)_t = 0` for all `t`. Writing
`h(x)=Σ_{i∈Diff}ε_i x^i∈ℤ[x]`, `Φ(δ)_t = 0 ⟺ h(ζ^t)∈ℝ`. Galois `σ_a` (`a∈(ℤ/C)^*`) sends this zero
condition at `t` to the same at `at` — so the **zero set is closed under `t↦at`**.

- **Prime `C`:** `(ℤ/C)^*` is transitive on nonzero residues ⟹ the zero set is `∅` or everything;
  `M` invertible kills "everything"; so `Φ(δ)` is nowhere zero, forcing `Φ(s)≡0`, i.e. `ε'=−ε`. No
  nontrivial collision. ∎ (orbit `=2^{n−2}`).
- **Composite `C`:** `(ℤ/C)^*` preserves `gcd(t,C)`, so the zero set can be a proper union of
  divisor-classes — exactly the room a collision needs.

This is the precise reason the prime/composite dichotomy is the **vanishing-sums-of-roots-of-unity**
dichotomy (Lam–Leung, Conway–Jones): a sine sum over a half-system vanishes nontrivially iff the
modulus is composite.

## The shape of every collision: prime-order-subgroup flips (HYP-2273)

Computing every collision for composite `C≤39` shows the flipped runner set is **always the
half-system `H_q` of a prime-order subgroup** `(C/q)ℤ/C`: flip the multiples of `C/q`. This move is
**frequency-localized** — it perturbs `Φ` only off the `q`-divisible frequencies (because the order-`q`
subgroup's sine sum collapses there). The `q=3` case `H_3={C/3}` is THM-413's order-3 silent single
flip. For squarefree `C=qr` the deficiency is `N_{(q)}+N_{(r)}`; for prime powers (e.g. `C=27=3³`)
the subgroup **lattice** generates a small group of moves (`{9}`, `{3,6,12}`, `{3,6,9,12}`), which is
why `27` is the most degenerate small modulus (69 mergers). The clean count
```
   deficiency(3p) = 2^{(p−1)/2}     (p>3 prime; verified p=5,7,11,13: 4,8,32,64),
```
splits as `2^{(p−3)/2}` order-3 flips `+ 2^{(p−3)/2}` "multiples-of-3" flips.

## Hidden connections this exposes

1. **To THM-413 (order-3 silent flip).** THM-413's Eulerian-graph criterion and "blind spot at
   `x=C/3`" is the `q=3` instance of the subgroup-flip mechanism. THM-413 proved `3∣C ⟹` degenerate;
   THM-415 proves the *converse-free* half (prime ⟹ non-degenerate) that THM-413 explicitly left
   open. Together they pin HYP-2270 down to a single missing lemma (composite ⟹ some `H_q` realized).

2. **To THM-414 (S704, additive face = multiplicative energy in the CM field = degree-2 Krawtchouk).**
   The peer's S704 face and this one are **two readings of the same `|DFT|²`**: S704 reads it as
   additive energy `E_+` / a Krawtchouk cut in the CM field `ℚ(ζ)`; S705 reads its imaginary part as
   the signed sine sum `Φ(ε)` and asks when the **squared** spectrum coincides. The matching-cap and
   "popular pair-sum mirror resolves negatively" of THM-414 sit on the same modulus arithmetic that
   makes the homometry zero set Galois-stable here. Bridging the two is an open, promising lane.

3. **To "everything is the triangle" / the prime-3 rosette.** The most degenerate sign-orbit is at
   `C=27=3³` (`n=14`) — the same `π/3` / Eisenstein / prime-3 locus that governs the worry-set
   shells (THM-401), the unit-distance density-6 optimum (THM-412), the Harborth `+3` frontier, and
   the strong-component scissors volume. The signed invariant's **blind spot** is the order-3 point;
   its **discriminating power** (the `V*` shell-partner at `n=14`, S699) is the same `3∣C` structure
   read the other way. Loneliness (cyclic order `ℤ/(2n−1)`) and kissing/unit-distance (Euclidean
   order, roots of unity `2,4,6`) keep turning out to be the same theorem in two metrics; homometry is
   the bridge that makes "same `|DFT|²`" literal on both sides.

## Honest status

- **PROVED:** homometry reframing; `|χ̂|² = A² + Φ²`; `M` invertible; **prime `C` ⟹ full sign-orbit**
  (THM-415); frequency-localization of subgroup flips.
- **VERIFIED:** clock = diff-multiset = `|DFT|²` partition equality (`n≤13`); the closed form (all
  cuts `n≤11`); prime ⟹ orbit `=2^{n−2}` (`n≤22`); subgroup-flip anatomy (all composite `C≤39`);
  `deficiency(3p)=2^{(p−1)/2}` (p=5,7,11,13).
- **CONJECTURE (HYP-2273):** the converse of HYP-2270 (composite ⟹ collision) via "some `H_q` is
  realized," and a uniform count law beyond `C=3p`.
- **Untouched (by design):** the observer gap `M` (sign-invariant, T1); the real (unsigned) LRC.

**Artifacts:** `04-computation/signed_lrc_homometry_s705.py`, `…sine_product_s705b.py`,
`…collision_anatomy_s705c.py`, `…c39_check_s705d.py` (+`.out`s). THM-415, HYP-2273, T760. Extends
THM-413/THM-401; parallels THM-414 (S704).
