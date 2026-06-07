# HYP-2318 — The general n=2p CRT fiber bundle: the p-clock base section (proven LRC(p)) + the mult-of-p fiber, whose doubling automorphism is the cube root (= Paley-7 at n=14)

**Session:** S640
**Status:** CONFIRMED (base section formalized + verified; the fiber-fit kernel is open, as in HYP-2346)
**Provenance forward:** math-lean `Math/LonelyRunner/SevenClockFiber.lean` (sorry-free)
**Extends:** the external "fiber bundle over the 7-runner base" prompt, refined by HYP-2346 (S643, the
mod-7 obstruction) and HYP-2317 (S639, the mod-2 half-turn is orthogonal); ties HYP-2225, S637/S638.

---

## 0. The point

The external prompt's instinct — fiber `n = 14` over the `7`-runner base via `14 = 2·7` — is correct
once you put the obstruction in the right summand (HYP-2346: it's the **mod-7** fiber; HYP-2317: the
half-turn is a benign **mod-2** tool, orthogonal). This session **generalizes** it to all `n = 2p`
(`p` an odd prime), proves the **base section** rigorously, identifies the **fiber** as a recursively
smaller LRC, and shows the **fiber's automorphism is the cube root** — which at `p = 7` is exactly the
Paley/`μ₃` structure of S638. So the fiber bundle is real and clean; the only open kernel is the same
"does the fiber dodge fit the window" of HYP-2346.

---

## 1. The base section (formalized, general p)

CRT: `ℤ/2p ≅ ℤ/2 × ℤ/7`. At the **p-clock** `t = b/p` (`gcd(b,p) = 1`), a runner of speed `v` has clock
distance
```
  ‖v·b/p‖ = (1/p)·min(s, p−s),    s = (v·b) mod p.
```
If `p ∤ v` then (with `gcd(b,p)=1`) `s ≠ 0`, so `min(s, p−s) ≥ 1` and `‖v·b/p‖ ≥ 1/p`. For `n = 2p` the
threshold is `δ = 1/(2p)`, and `1/p = 2δ > δ`:

> **Every non-multiple-of-`p` runner is automatically clear at the p-clock, with a factor-2 margin.**
> Hence **LRC(p)** (the base, *proven* for `p ≤ 7`, far inside the `k ≤ 6` frontier) is a genuine
> section of the bundle over the proven base — it disposes of all non-mult-of-`p` runners at once.

Verified for `n ∈ {6, 10, 14, 22, 26}`: the worst-case non-mult-of-`p` margin is exactly `1/p = 2δ`
every time (`n2p_fiber_bundle_doubling_qr_s640.py`).

**Formalized (math-lean, sorry-free):** `pclock_margin {p} [NeZero p] (s : ZMod p) (hs : s ≠ 0) :
1 ≤ min s.val (p − s.val)` — the general base-section margin, for *all* `p`.

---

## 2. The fiber is a recursively smaller LRC (the protection chain `2p → p`)

The runners *not* handled by the base are the **mult-of-`p`** speeds `v = p·w`. Near `t = b/p + ε`:
```
  v·t = w·b + p·w·ε  ≡  p·w·ε  (mod 1),   so  ‖v·t‖ = ‖p·w·ε‖,
```
i.e. in the perturbation window the speeds `{p·wᵢ}` at time `p·ε` behave **exactly as an LRC on the
reduced speeds `{wᵢ} = {vᵢ/p}`**. So:

> **LRC(2p) ⟸ [ LRC(p) clears the base ] ∧ [ the `{v/p}`-LRC fits the perturbation window ].**

This is the genuine "fiber bundle over the `p`-runner base": base = LRC(p), fiber = LRC on the
`mult-of-p`/`p` speeds, depth-1 protection chain `2p → p` (recursive). For speeds `1..13` the fiber is a
*single* runner (`v = 7 → w = 1`), trivially dodgeable; the hardness lives only in speed-**sets** with
several mult-of-`p` (HYP-2346's distribution `1:319, 2:689, 3:417, 4:75`). **Open kernel** (= HYP-2346):
prove the `{v/p}`-LRC lonely time always lands inside the window where the base stays clear.

---

## 3. The fiber automorphism is the cube root (= Paley-7 at n=14)

The `2` of `14 = 2·7` acts on the 7-clock by **doubling** `t ↦ 2t`, i.e. `×2` on `ℤ/7`, of order
`ord₇(2) = 3` — a **cube root / 3-cycle**. Its two orbits on `(ℤ/7)\{0}` are
```
  1 → 2 → 4 → 1      and      3 → 6 → 5 → 3,
```
which are **exactly the quadratic-residue / non-residue cosets** `{1,2,4}` / `{3,5,6}` = the **Paley-7
connection set** (`paleySet`, S638) = the **cube roots of unity `μ₃`** and their complement. So:

> The 7-fiber's doubling symmetry **is** the Paley / `μ₃` / `ω`-resonance the whole arc converges on
> (`7 = Φ₃(2) = N(3+ω)`). The "2-and-3 seam" of `14 = 2·7` (HYP-2225) is literally: the **half-turn `7`**
> carries the `ℤ/2` summand (a mod-2 detector, blind to the fiber — HYP-2317), while **doubling `2`**
> carries the `ℤ/7` summand as the **cube root**. The two prime structures of `14` are the two halves of
> the perspective key — `σ` (order 2, the apex) and `ω` (order 3, the resonance).

**Formalized:** `seven_doubling_cube` (`2³ = 1` in `ℤ/7`), `doubling_threecycle_QR`/`_nonQR` (the two
3-cycles), `doubling_orbit_eq_paley` (`{1,2,4} = paleySet`, linking to PaleyRado.lean / S638).

---

## 4. The coincidence sweep: where doubling = Paley (new)

When is the doubling orbit equal to the QR coset (`⟨2⟩ = QR`, so doubling fills a whole cube-root
coset)? Iff `2` is a QR (`p ≡ ±1 mod 8`) **and** `ord_p(2) = (p−1)/2` (2 generates the QR subgroup). When
additionally the **Paley tournament exists** (`p ≡ 3 mod 4`, so `−1` is a non-residue, S638), both hold
iff `p ≡ 7 (mod 8)`:
```
  ⟨2⟩=QR:        p = 7, 17, 23, 41, 47, …
  Paley (p≡3/4): p = 3, 7, 11, 19, 23, 43, 47, …
  BOTH (p≡7/8):  p = 7, 23, 47, 71, …     ← n = 2p = 14, 46, 94, …
```

> **`p = 7` (n = 14) is the smallest prime where the LRC(2p) fiber's doubling-dynamics fill the entire
> cube-root coset *and* the Paley tournament exists.** The fiber automorphism and the Paley/CM/cube-root
> structure (S638) coincide exactly at the `n = 14` frontier — not a coincidence the arc imposed, but one
> the arithmetic forces. The next such frontier is `n = 46 = 2·23`.

---

## 5. Connections through the repo

- **HYP-2346 (S643):** this is the general form of S643's `n=14` reduction; §1 proves the base section
  for all `p`, §2 makes the fiber recursion explicit, §3/§4 add the doubling=cube-root layer.
- **HYP-2317 (S639):** the `ℤ/2` summand (half-turn) is the orthogonal benign part; here the `ℤ/7`
  summand (doubling) is the active cube-root part. Together = the full CRT Künneth split.
- **HYP-2225 (the 2-and-3 seam):** `ord₇(2) = 3`, doubling = two 3-cycles with σ-pairs — now grounded as
  the fiber automorphism of the LRC(14) bundle.
- **S637/S638:** `7 = Φ₃(2) = N(3+ω) =` Paley-7; the fiber's QR/μ₃ structure is the *same* 7.
- **HYP-2185 (perspective key):** `σ` (order 2, apex, mod-2) and `ω` (order 3, resonance, mod-7) are the
  two CRT halves of `14`; the perspective key splits along the divisors of `n`.

## 6. Open / handoffs
- The open kernel (shared with HYP-2346): prove the `{v/p}`-fiber-LRC lonely time fits the base's clear
  window — a quantitative "window vs reduced-LRC-depth" inequality. This is the real remaining content of
  LRC(2p); everything else (base section, fiber identification, automorphism) is now settled/formalized.
- Whether the `p ≡ 7 mod 8` frontier (n = 14, 46, 94, …) is genuinely *harder* (maximal cube-root
  alignment) or *easier* (most symmetric) than generic `n = 2p` — a sharp, testable question.
