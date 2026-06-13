---
source: opus-2026-06-03-S580 (remote-control)
status: STRUCTURE — translated APs as the binary IFS {x→2x, x→2x+1}; the doubling map x→2x mod n is the phase dynamics; clean cyclic iff n odd; degenerate at n=14=2·7 (CRT: 2-collapse × prime-7 cycle)
tags: [LRC, translated-AP, binary-IFS, doubling-map, primitive-root, CRT, 2-adic, even-n, n14, fractal]
---

# Translated APs, the binary IFS, and the doubling degeneracy at n = 2·7

**Prompt (user):** consider the recursive fractal concept of translated APs.

Two affine maps generate everything: `D(x)=2x` (double) and `T(x)=2x+1`
(double-and-**translate**). They are the binary IFS whose attractor is `ℕ`; the AP is
its depth-truncation, and the **translate** branch `T` produces exactly the odd binders
of the 2-adic tower (S579). Pushed to the phase circle, `D` becomes the **doubling map**
`x↦2x mod 1` — the canonical fractal/chaotic map — and its arithmetic
`x↦2x mod n` is a clean cycle iff `n` is odd, degenerating *precisely* at the even
frontier, with `n=14` the first composite-`2q` failure.

## 1. The binary IFS: translated APs build the integers

```
D(x) = 2x        (even branch — the scaled AP)
T(x) = 2x + 1    (odd branch — the AP translated by 1, then scaled)
```
From the root `1`, `{D,T}^*` generates every positive integer (its binary expansion).
`AP_{n-1}={1,…,n-1}` is the breadth-truncation. The 2-adic strata of S579 are the
images: even stratum `= D(·)`, odd stratum `= T(·)` — so the **odd binders are exactly
the `T`-translated branch**. "Translated AP" is the `+1` in `T`.

## 2. Literal translates `{c+1,…,c+n-1}`: the base AP is the unique tight one

Sliding the window (a literal translation of the AP):

| n | c=0 | c=1 | c=2 | c=3 | … |
|---|---|---|---|---|---|
| 14 | `M=1/14` **(tight)** | .125 | .167 | .200 | ↗ |

For `c=0` the translate is the tight AP; for every `c≥1` the length-`(n-1)` window
**slides a multiple of `n` inside** and `M` rises *monotonically* (margin `+0.05,
+0.10, +0.14,…`). So **the base AP is the unique tight translate**, and translation has
a measurable cost — a clean instance of C′ (multiple ⟹ loose) along the translation
orbit. (`lrc_translated_AP_fractal_s580.py`.)

## 3. The doubling map `x↦2x mod n` is the phase dynamics

At `t=1/n` the AP phases are `{v/n}`, and `D` acts on them as `v↦2v mod n` — i.e. the
**binary shift** on the phase circle. Its orbit structure:

| n | orbit of the binder `1` under `x↦2x` | nature |
|---|---|---|
| 7 (prime) | `1→2→4→1` | `ord_2=3`, clean cycle |
| 13 (prime) | `1→2→4→8→3→6→12→11→9→5→10→7→1` | **2 a primitive root**, one full 12-cycle |
| **14 = 2·7** | `1→2→4→8→2…` | **degenerate**: `1` is a *transient source*, never returns |
| 16 | `1→2→4→8→0` | collapses to `0` |

> **`x↦2x mod n` is a permutation of `{1,…,n-1}` iff `n` is odd (2 invertible).** For
> `n` prime it is rich (often a single long cycle; primitive root at `13`). At every
> **even** `n` it degenerates — and `n=14` is the first composite `2q` where the binder
> `1` falls out of every cycle.

This is why **odd/prime `n` is "fractal-regular"** (the doubling cycles organise the
phases — exactly the multiplicative structure the 2025 polynomial/sieve method
exploits) and **even `n` is the hard frontier** (no doubling cycle to ride).

## 4. CRT: `n=14` splits into a clean prime cycle and a 2-collapse

By CRT (`14 = 2·7`), `v mod 14 ↔ (v mod 2, v mod 7)`, and the doubling map factors:
```
x↦2x mod 14  =  (x↦2x mod 2)  ×  (x↦2x mod 7)
             =  ( 1 → 0,  collapse )   ×   ( 1→2→4→1,  clean order-3 cycle ).
```
The **mod-7 factor is the clean prime doubling** (the solved case); the **mod-2 factor
is the irreversible `1→0` collapse** — the entire obstruction. This is *exactly* the
2-adic tower's seam (S579: even layer = `D(`AP_6`)` = scaled prime-7; odd layer = the
binder) and the apex-`2q` (HYP-2063) in dynamical dress: **the "2" in `2·7` is the
non-invertible coordinate that breaks the fractal regularity.**

## 5. Synthesis — the recursive fractal in one picture

- **Generators** `{D=×2, T=×2+1}` build `ℕ` (binary tree); the AP is a truncation; the
  worry-set's translated copies are the `T`-branch (odd binders) ⊔ flip-translations
  (the `2^((n-2)/2)` family, S570).
- **Phases** evolve by the doubling map `x↦2x mod 1`; on the `n`-clock it is `x↦2x mod
  n`, a clean permutation ⟺ `n` odd.
- **The frontier** is even `n` = degenerate doubling; `n=14=2·7` = (clean prime-7 cycle)
  × (mod-2 collapse) by CRT. The literal translates confirm the base AP is the unique
  tight point, with translation cost = the C′ margin.

## 6. Extension (the actionable handle)

> **Proof split by CRT along the doubling factorisation.** For `n=2^a·q` (`q` odd), ride
> the clean doubling cycle on the `q`-coordinate (the prime/odd part — where the
> multiplicative method works) and treat the `2^a`-coordinate's collapse separately.
> For `n=14`: the mod-7 part is the solved prime cycle; the residual is the single
> mod-2 collapse — the *same* odd-layer / C′ / `Φ>0` residual reached from the
> 2-adic tower (S579) and the multiple split (THM-398). The "translated AP" `T`-branch
> is exactly where that residual lives.

## 7. Honest status

- **Binary IFS, literal-translate uniqueness, doubling-map orbits, CRT factorisation:**
  **verified** (`n=7,13,14,16,…`). Base AP unique tight translate; `2|n ⟺` degenerate
  doubling.
- The dynamical re-description unifies S579 (2-adic tower), HYP-2063 (apex `2q`), and
  C′; it explains *why* odd/prime `n` is tractable (doubling cycles) and `14` is the
  first hard `2q`. The CRT proof-split (Extension §6) is a structural direction, not a
  theorem.

**Artifacts:** `04-computation/lrc_translated_AP_fractal_s580.py` (+`.out`),
`lrc_translated_AP_crt_s580.out`. Builds on S579 (2-adic tower), HYP-2063 (2q apex),
HYP-2097 (flips), THM-398 (C′), THM-369. New: **HYP-2117**.
