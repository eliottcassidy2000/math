---
id: THM-1605
title: "THE OUTSIDE INFINITE FAMILY: THEIR m=2 IS OUR MAP, AND OUR THM-1350 EXPLAINS THEIR HEADLINE COUNT. An outside party exhibits E_m : C^3 -> C^3, Keller and non-injective for every m >= 2, with deg E_2 = 7 and #E_m^{-1}(c,0,0) = 2m-1, concluding 'every odd fibre cardinality >= 3 occurs'. VERIFIED IDENTIFICATION: our THM-1300 counterexample has deg 7, det JF = -2, fibre 3 -- so their m=2 IS our case and m >= 3 is genuinely new to us. HONEST VERDICT ON EXTENT: their family is STRICTLY BROADER than anything we hold (infinitely many degrees 7,13,26,43,64,...; every odd cardinality), and we should not pretend otherwise. HONEST VERDICT ON MECHANISM: they EXHIBIT odd cardinalities; THM-1350 PROVES oddness is FORCED -- for a sigma/tau-equivariant Keller map with dim Fix(sigma) <= 1, F restricted to Fix(sigma) is injective by JC_1, so every fibre over a tau-fixed target contains EXACTLY ONE sigma-fixed preimage and the rest fall into free 2-orbits, giving |fibre| = 1 + 2k: a DOUBLE is impossible and 3 is the minimum. Their 2m-1 factors as 1 + 2(m-1) on the nose. So their 'every odd cardinality >= 3 occurs' is the EXISTENCE half of a statement whose NECESSITY half we proved independently, and their write-up shows no sign of the necessity direction"
status: VERIFIED (the identification: degree, determinant, fibre size, and the 1+2(m-1) split checked exactly on our map)
author: opus-2026-07-20-S415
depends_on: [THM-1300 (our counterexample), THM-1350 (the odd-fibre forcing theorem), THM-1330, THM-1440, THM-1445]
---

# THM-1605 — Extent versus mechanism: the outside infinite family

## 1. The identification (verified)

| | outside family at `m = 2` | our THM-1300 map |
|---|---|---|
| degree | **7** | **7** (`u³z` with `u = 1+xy`) |
| fibre over `(c,0,0)` | `2m−1 = 3` | **3** (`F⁻¹(1,0,0)` = `{(0,0,1)} ∪ {(±i/2,±3i,−26)}`) |
| `det J` | constant | `−2` |

**Their `m = 2` is our case.** Their `m ≥ 3` — degrees `13, 26, 43, 64, …`, fibres
`5, 7, 9, …` — is genuinely new to us.

## 2. Where they are ahead, stated plainly

**Their family is strictly broader than anything this repo holds.** We have *one* Keller
counterexample; they have an infinite family realising infinitely many degrees and *every*
odd fibre cardinality `≥ 3`. Nothing in our canon produces a second counterexample, let
alone a family. On extent, they win, and the repo should record that rather than hedge.

## 3. Where we are ahead: necessity, not existence

Their headline is *"#F⁻¹(c,0,0) = 2m−1: every odd fibre cardinality ≥ 3 occurs."* That is an
**existence** statement — the construction realises those values. **THM-1350 proves the
matching necessity statement**, which their write-up shows no sign of:

> **THM-1350 (opus-S399, repaired S400).** Let `F` be `σ/τ`-equivariant with
> `det JF ∈ ℂ*` and `dim Fix(σ) ≤ 1`. Then `F` maps `Fix(σ) → Fix(τ)` and the restriction is
> again a constant-Jacobian polynomial map in dimension `≤ 1` — i.e. `JC₁`, which is **true**
> — so `F|_Fix` is injective **unconditionally**. Hence every fibre over a `τ`-fixed target
> contains **exactly one** `σ`-fixed preimage, and the remaining preimages fall into **free
> `σ`-orbits of size 2**. Therefore
> ```
> |fibre| = 1 + 2k :  ODD.  A collision can never be a DOUBLE.  The minimum is a TRIPLE.
> ```

**Their count factors on the nose:**

```
2m − 1  =  1 + 2(m−1)   =   one σ-fixed preimage  +  (m−1) free σ-orbits
```

So the *shape* of their result is predicted by ours. They show every odd value **is
attained**; we show **no even value can be**. Those are the two halves of one theorem, and
only one of them is a construction.

**Verified on `m = 2`.** Over the `τ`-fixed target `(a,0,0)` the fibre cubic (THM-1440) is
`16a·x³ + 4x = 4x(4ax²+1)`, with roots

```
x = 0          <- the σ-fixed sheet (σ = (−x,−y,z)); it is the point (0,0,a), since F(0,0,z) = (z,0,0)
x = ± i/(2√a)  <- a single free σ-orbit
```

`1 + 2·1 = 3 = 2m−1` at `m = 2`. Exactly the predicted split.

## 4. The specific example that carries the extra structure

Take our `m = 2` map and ask what can be said about it that the family write-up does not say
about any of its members:

1. **The fibre is the root set of an explicit cubic** whose leading coefficient *is* the
   Jelonek polynomial:
   ```
   P(x) = L·x³ + (4 − 3bc)·x − 2c ,   L = 27a²c² − 18abc + 16a + b³c − b²
   ```
   `L = 0` is exactly where sheets escape to infinity (THM-1440 §1, cross-validated against
   THM-1330's independently derived `−L|_{c=0} = b² − 16a`).
2. **The Jelonek set is Zariski's 1929 three-cuspidal quartic**, zero nodes — the worst
   rational quartic, not merely non-nodal — with monodromy `S₃` and self-normalising
   stabiliser (THM-1330/1375). "Non-injective" is the coarse fact; *which* degeneration
   occurs is this.
3. **Two sheets are lost over `{L = 0}`, not one**: the fibre drops `3 → 1`, surviving root
   `x = 2c/(4−3bc)` (THM-1440 §2b).
4. **Reflection = torus.** The meridian monodromy around the Jelonek quartic *equals* the
   `λ = −1` torus involution — and by THM-1445-A this holds for **any** `σ/τ`-equivariant
   degree-3 Keller map, because `σ` is linear so it permutes the two escaping sheets, and it
   cannot fix both (that would be two `σ`-fixed points, contradicting THM-1350). *An
   involution cannot fix two things when it is only allowed to fix one.*

Item 4 is the sharpest contrast: it is a statement about **every** member of their family at
`m = 2`-type degree, derived with no construction at all.

## 5. A concrete prediction for their family

If each `E_m` is `σ/τ`-equivariant with `dim Fix(σ) ≤ 1` — which its `(t = 1+xy, p = x²h)`
coordinates strongly suggest — then THM-1350 forces, for every `m`:

- the `2m−1` preimages split as **1 `σ`-fixed + (m−1) free `σ`-orbits**;
- `σ` acts on the fibre as a product of **exactly `m−1` disjoint transpositions**, fixing one
  sheet;
- consequently, at a `τ`-fixed simple Jelonek point where exactly two sheets escape, the
  meridian monodromy is a **transposition inside one of those orbits** (THM-1445-A).

None of this requires their construction; all of it is testable against it. **That is the
form our advantage takes: we cannot build their family, but we can predict its fibre
combinatorics before looking.**

## 6. Verdict

- **Extent:** theirs, decisively. Infinite family vs our single map.
- **Mechanism:** ours. Necessity of oddness, the `1 + 2k` split, the Jelonek singularity
  type, the monodromy identification — none of which follows from having a construction.
- **Not in competition:** their `m = 2` *is* our map. The right reading is that they extended
  the object and we explained it; the two results compose.

## Verification

`04-computation/family_comparison_opus_S415.py` — degree, determinant, and fibre of our map;
the `1 + 2(m−1)` factorisation; the `τ`-fixed fibre roots `{0, ±i/(2√a)}` exhibiting the
one-fixed-plus-one-free-orbit split. Output in `05-knowledge/results/`.
