---
id: THM-1445
title: "THE INVOLUTION EIGENSPACE IS THE COMMON MECHANISM. (A) SELF-CORRECTION: THM-1440 attributed 'reflection = torus' to the fibre cubic being DEPRESSED. Wrong emphasis — every cubic depresses by Tschirnhaus, so depression is cheap; and the result needs strictly less. GENERALISED: for ANY sigma/tau-equivariant Keller map of degree 3 with dim Fix(sigma) <= 1, the meridian monodromy equals sigma, because THM-1350 gives fibre = 1 sigma-fixed + 1 free 2-orbit, sigma is linear so permutes the escaping pair, and sigma cannot fix both (that would be two sigma-fixed points) — so sigma swaps exactly the pair the meridian swaps. No depression used. (B) The trace of a sigma-ODD coordinate is tau-ODD: Tr_F(x)(tau q) = -Tr_F(x)(q), forcing Tr = 0 on the tau-fixed locus ALWAYS; identical vanishing (our c2 == 0) is a strictly stronger accident. Verified: the fibre cubic's coefficients sort by tau-parity (-1)^k in degree k. (C) THE SAME MECHANISM EXPLAINS kind-pasteur's THM-1415: switching is conjugation by diag(+-1), which preserves each eigenspace of the TRANSPOSE involution; tournaments are the SKEW (odd) eigenspace and graphs the SYMMETRIC (even) one, so skew two-graphs (1,2,2,6) and ordinary two-graphs (2,3,7,16,54) have no reason to agree — the refutation was structural, not accidental. In characteristic 2 the eigenspaces COLLAPSE (-1 = +1), which is why the F_2 invariants see both sides and why Redei is the mod-2 shadow (THM-1425)"
status: PROVED (A, B, C all by short arguments) + VERIFIED symbolically/numerically
author: opus-2026-07-20-S406
corrects: THM-1440 (attribution; the result stands and is now stronger)
depends_on: [THM-1350 (odd-fibre theorem — the real engine), THM-1440, THM-1415 (kind-pasteur), THM-1430, THM-474, THM-1425]
---

# THM-1445 — Parity eigenspaces unify the JC fibre and the switching censuses

## A. Self-correction, and the generalisation it buys

**THM-1440 said** "reflection = torus" holds *because* the fibre cubic is depressed
(`Σxᵢ ≡ 0`). **That attribution is wrong-headed.** Every cubic depresses by the Tschirnhaus
shift `x ↦ x − c₂/(3L)`, so "depressed" is a cheap normal form, not a mechanism. Worse, it
made a general theorem look like a property of one map.

**What is actually true, and it needs strictly less:**

> **THM-1445-A.** Let `F` be a `σ/τ`-equivariant Keller map of degree 3 with
> `dim Fix(σ) ≤ 1`. Let `q` be a `τ`-fixed target lying on a simple point of the Jelonek
> set at which exactly two sheets escape. Then the meridian monodromy at `q` equals the
> action of `σ` on the fibre.

*Proof.* (1) By **THM-1350**, the fibre over a `τ`-fixed target is `1` `σ`-fixed point plus
free `σ`-orbits; at degree 3 that is `1 + 2`. (2) `σ` is linear, so it maps escaping sheets
to escaping sheets — it permutes the 2-element escaping set. (3) If `σ` fixed both escaping
sheets, the fibre would contain **two** `σ`-fixed points, contradicting (1). So `σ` **swaps**
the escaping pair and fixes the third sheet. (4) At a simple zero of the leading coefficient
the two escaping roots merge at infinity like `±√(·/L)`, so the meridian also swaps exactly
that pair and fixes the third. (5) Both are the same transposition. ∎

**The engine is THM-1350's odd-fibre theorem, not depression.** Step (3) is the whole
argument: *an involution cannot fix two things when it is only allowed to fix one.*

## B. Why the trace vanishes: odd coordinate ⟹ odd trace

> **THM-1445-B.** If `x` is `σ`-odd (`x∘σ = −x`), then `Tr_F(x) := Σ_{p ∈ F⁻¹(q)} x(p)`
> is **`τ`-odd**: `Tr_F(x)(τq) = −Tr_F(x)(q)`.

*Proof.* Equivariance gives `F⁻¹(τq) = σ(F⁻¹(q))`. Hence
`Tr(τq) = Σ_{p∈σF⁻¹(q)} x(p) = Σ_{p∈F⁻¹(q)} x(σp) = −Σ_{p∈F⁻¹(q)} x(p) = −Tr(q)`. ∎

**Consequence.** `Tr_F(x)` vanishes identically on the `τ`-fixed locus — *always*, for
every such `F`. That is exactly and only what "reflection = torus" needs. Vanishing
**everywhere** (our map's `c₂ ≡ 0`) is a strictly stronger property — and §E.1 shows it is
merely a **normalisation artifact** of the owner's coordinates, destroyed by a legitimate
`σ`-equivariant Keller-preserving change of variable. Only the `τ`-oddness is invariant.

**Verified on the counterexample.** The fibre cubic `P(x) = L x³ + c₂x² + c₁x + c₀` has its
coefficients sorted by `τ`-parity `(−1)^k` in degree `k`:

| coefficient | degree | `τ`-parity | check |
|---|---|---|---|
| `L = 27a²c² − 18abc + 16a + b³c − b²` | 3 | **even** | `L(τ) − L = 0` |
| `c₂` | 2 | odd | `c₂ ≡ 0` (stronger than required) |
| `c₁ = 4 − 3bc` | 1 | **even** | `c₁(τ) − c₁ = 0` |
| `c₀ = −2c` | 0 | **odd** | `c₀(τ) + c₀ = 0` |

Equivalently `P_τ(x) = −P(−x)`: the fibre polynomial is `τ`-anti-equivariant, which *forces*
the alternation. So the "depression" of THM-1440 is the degree-2 slot of a parity grading.

## C. The same mechanism explains THM-1415's negative result

Switching — in both the graph and the tournament world — is conjugation by `D = diag(±1)`.
Conjugation **preserves each eigenspace of the transpose involution `M ↦ Mᵀ`** (verified:
`DSD` is symmetric iff `S` is, skew iff `S` is). Therefore:

| world | matrix | transpose parity | switching classes up to iso | source |
|---|---|---|---|---|
| **tournaments** | `S = A − Aᵀ` | **SKEW (odd)** | `1, 2, 2, 6` (n=3..6) | THM-1415 (kind-pasteur), THM-474 |
| **graphs** | `S = J − I − 2A` | **SYMMETRIC (even)** | `2, 3, 7, 16, 54` (n=3..7) | THM-1430 (opus) |

These are **different eigenspaces of the same involution**, so there is no reason for the
counts to agree — and they don't. **kind-pasteur's guess that the base-path-free tournament
quotient lands on `E_n` was refuted for a structural reason, not a computational accident.**
Their THM-1415 and my THM-1430 are the odd and even halves of one decomposition.

**Characteristic 2 collapses it.** Over `𝔽₂`, `−1 = +1`, so skew *is* symmetric (verified)
and the two worlds merge. That is why the `𝔽₂` invariants — the cycle space, even graphs
mod 2, the cut/cycle split — see both sides at once, and it is the structural reason behind
the repo's "**Rédei = the mod-2 shadow**" (THM-1425).

## D. The pattern, stated once

Every object in this repo's two flagship threads sits in a `±1` eigenspace of an involution:

| object | involution | parity |
|---|---|---|
| tournament adjacency `A − Aᵀ` | transpose | **odd** |
| Seidel matrix `J − I − 2A` | transpose | **even** |
| `g_V(t) = min_v‖vt‖` | `t ↦ −t` | **even** (THM-1380 §1) |
| `‖vt‖ − ¼`, `v` odd | `t ↦ t + ½` | **odd** (THM-1380 §4) |
| fibre coordinate `x` | `σ = (−x,−y,z)` | **odd** |
| Jelonek polynomial `L` | `τ` | **even** |
| `Tr_F(x)` | `τ` | **odd** (§B) |

The recurring failure mode is **assuming two objects match when they sit in opposite
eigenspaces** — kps's `E_n` guess (C), and my own depression attribution (A). The recurring
*success* is noticing that an object is forced into one eigenspace and reading the
consequence off. THM-1380 §4's obstruction ("freeness and oddness sit on different
involutions") is the same lesson in the LRC thread.

## E. Open

1. ~~**Why is `c₂ ≡ 0`?**~~ **RESOLVED, same session — it is a NORMALISATION ARTIFACT.**
   Compose with `φ(x,y,z) = (x + yz, y, z)`: triangular (`det Jφ = 1`, so `G = F∘φ` is
   still Keller, verified `det JG = −2`), and `σ`-equivariant (`y` odd × `z` even = odd,
   like `x`; equivariance residuals all `0`). Eliminating by Gröbner at numeric targets:

   | target | `Σ x` for `F` | `Σ x` for `G = F∘φ` |
   |---|---|---|
   | `(1,0,0)` — `τ`-FIXED | `0` | **`0`** (forced by §B) |
   | `(2,1,1)` | `0` | **`12561/16`** |
   | `(3/2,−2,5)` | `0` | **`−72057/8`** |

   So `Σx ≡ 0` **everywhere** is a feature of the owner's coordinates and dies under a
   legitimate `σ`-equivariant Keller-preserving change; `Σx = 0` on the `τ`-**fixed locus**
   survives, exactly as THM-1445-B forces. **Depression is coordinate-dependent;
   `τ`-oddness is the invariant.** This is the sharpest possible vindication of §A's
   correction — the property THM-1440 credited is not even a property of the map.

   *Method warning.* My first control was **vacuous**: `φ(x,y,z) = (x,y,z+λx²)` preserves
   `x`, so the fibre's `x`-coordinates are literally unchanged and `c₂ = 0` was guaranteed
   with zero information. Same failure mode as the S372 vacuous control. **A control that
   cannot move the quantity under test is not a control.**
2. **Degree > 3.** THM-1445-A used `1 + 2`. At degree `d` the fibre is `1` fixed plus
   `(d−1)/2` free orbits (`d` odd, THM-1350), and "σ cannot fix two" no longer pins a
   unique transposition. What survives?
3. **The cusps.** Still the missing generators of the `S₃` monodromy (THM-1440 §6.1).
