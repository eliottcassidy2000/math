---
id: HYP-2230
title: The unit-distance ↔ LRC bridge — CM rotations, the perspective orbit, and the shell/class-field tower
status: OPEN (connective frame; rigorous parallels + a research direction). n=22 facts cited.
source: claudebox-2026-06-03-S624
related:
  - HYP-2130  # rigidity = orbit-type / the perspective key (= the multiplier orbit)
  - THM-407   # doubling stops being shell-transitive at odd prime power 2n-1 (= ramification)
  - THM-412   # twisted-shell dodge (the multiplier/perspective escape)
  - THM-415   # quantitative C'(n): the 2n-1 shell extremal
  - HYP-2140  # the 2-adic seam
---

# HYP-2230 — unit distances, the grid-conjecture disproof, and LRC

## The n=22 anchor (the open frontier)

Erdős unit-distance maxima `u(n)` are proven optimal only for `n ≤ 21` (Alexeev–Mixon–Parshall
2024): `u(21) = 57`. **`u(22) ∈ {60, 61}` is open** (lower bound 60, upper 61 — a single edge).
A compact 22-point **triangular-lattice** subset realises only **49** unit distances
(`unit_distance_n22_s624.py`); the optimum 60 needs **+11 from rotations**. The triangular lattice
is the Eisenstein CM lattice `ℤ[ω] = ℚ(√−3)` (the densest small unit-distance source); the
enhancing rotations are the **bounded-height modulus-1 elements `u = α/ᾱ` of `ℚ(ω)`** — angles
`arccos(13/14)=21.79°, 38.21°, arccos(rational)…`. These are exactly the objects the grid-conjecture
disproof scales up.

## The grid-conjecture disproof, in one line

Erdős conjectured the **grid is asymptotically optimal**: `u(n) = n^{1+o(1)}`. The recent
(AI-generated, human-verified — arXiv 2605.20695) **disproof** builds `n^{1+ε}` via: a **CM field**
with **many modulus-1 algebraic numbers of bounded height** (pigeonhole over ideal classes),
supplied by an **infinite class field tower** (Golod–Shafarevich, bounded root discriminant) with a
**split prime** as multiplier source, then a **geometry-of-numbers window + projection** turning
"differ by a modulus-1 element" into "distance 1 in the plane." The grid loses because it has too
**few rotations**; the CM tower wins by having **many**.

## The bridge (rigorous parallels)

Both problems **count metric resonances in a lattice, governed by its rotational/multiplicative
symmetry group**. The dictionary:

| | unit distance | LRC | tournaments |
|--|--|--|--|
| ambient | `ℂ` | clock `ℝ/ℤ` | phase ring `ℤ/n` |
| resonance | `‖p_i−p_j‖ = 1` | `‖v_i t‖ = 1/n` | parity/arc relation |
| base lattice | grid `ℤ[i]` / `ℤ[ω]` | AP / speed lattice | the AP witness |
| **rotation group** | **CM units `α/ᾱ`** (modulus-1 alg. numbers) | **multiplier `(ℤ/m)*`** (the perspective) | doubling `⟨2⟩`, Aut-orbits |
| rigid extreme | square grid (few rotations) | AP / observer-coupled | transitive/rigid |
| symmetric extreme | high-degree CM lattice (many rotations) | observer-blind / multi-orbit | regular / big Aut |
| controlling tower | **class field tower** (Golod–Shafarevich) | **`2n−1` shell tower** (THM-407: `27=3³`) | Cayley–Dickson `ℝ→ℂ→ℍ→𝕆` |
| prime seam | split prime `q≡1 (4)` in `ℤ[i]` | 2-adic seam / `2n−1` ramification | doubling tower |

**The two non-obvious, structurally identical facts:**

1. **Modulus-1 algebraic numbers of bounded height = the perspective/multiplier orbit.** Roots of
   unity (the literal lattice rotations = the LRC witness orbit `(ℤ/m)*`, [[HYP-2130]]) are *rare*;
   the dense supply of bounded-height modulus-1 elements in a CM field is the *same* object scaled
   up. The disproof's engine **is** the perspective key. The small-`n` unit-distance optimum sits in
   `ℚ(√−3)` for the same reason the AP (the `(ℤ/m)*`-symmetric config) is the LRC extremal.

2. **The grid conjecture is the unit-distance "the AP is the unique extremal".** Both are false the
   same way: the AP is *not* the unique LRC tight config (sporadics, S620/S621), and the grid is
   *not* unit-distance-optimal (disproof). In both, the correction is **more rotational symmetry**
   (CM / multiplier orbits) yields **more** extremal resonance. The disproof is the unit-distance
   analogue of finding the LRC sporadic/extremal configs by twisting the rigid object by the (dense)
   rotation group.

## The inspiration (research direction for LRC)

> **LRC's n=14 frontier is the finite shadow of the Golod–Shafarevich infinite-tower phenomenon
> that powers the grid disproof.** The same class-field-tower / odd-prime-power ramification that
> beats the grid is what makes `2n−1 = 27 = 3³` (THM-407) the first even LRC frontier.

Concretely, to import the disproof's machinery into the constructive route ([[HYP-2220]], `M ≥ 2/27`):
- Treat the `2n−1 = 27` shell as a **truncated CM/cyclotomic tower** `ℚ(ζ_27)`; the witness
  multipliers `a ∈ (ℤ/27)*` are the finite analogue of the modulus-1 algebraic numbers.
- The twisted-shell dodge's **free unit ±-pair** (THM-412) is the finite analogue of the disproof's
  **window modulus-1 vector**: both project lattice rotational structure onto a single resonance.
- **Conjecture (window↔dodge):** a multiple-of-`n` config fails the `2n−1`-shell dodge iff its
  residues mod `2n−1` block every `(ℤ/(2n−1))*` ±-orbit — the finite, ramified analogue of "the CM
  lattice has no free modulus-1 window." The 3-adic ramification (no doubling-transitivity) is
  exactly why the `3³` case needs the *tower*, not a single shell — the same reason the disproof
  needs a *tower* of CM fields, not one.

## To do
1. Make the "window ↔ shell-dodge" analogy a precise statement at `2n−1 = 27` (cyclotomic `ℚ(ζ_27)`,
   the `±`-orbit = `α/ᾱ` count) and test it against the S622/S623 data.
2. Is the LRC extremal (`M = 2/(2n-1)`, THM-415) the config whose speed-residues mod `2n−1` are a
   *CM-rotation orbit*? (the unit-distance "ℤ[ω] is densest" ↔ LRC "AP/CM-orbit is extremal".)
3. Read the disproof's window lemma and Golod–Shafarevich step in detail; port the "bounded root
   discriminant" control to a *uniform-over-n* shell-dodge bound (a real attack on LRC asymptotics).
