---
id: HYP-3766
title: THE (n+q)-WITNESS RIGORIZED + THE MULTI-METRIC WITNESS SHEAF + rung realizability (H7). (A) RIGOROUS (n+q)-witness (tight regime, single swap): for q COPRIME to n with (n-1)/2 < q <= n-1, dropping q from the tight AP {1..n-1} and adding any small multiple mq gives M >= 2/(n+q) > 1/n. Proof: at a=q^{-1} mod (n+q), the residues 1 and -1 pull back to speeds q and n -- q is DROPPED and n is OUT OF RANGE, so both are absent; every remaining runner avoids {0,+-1} mod (n+q) => radius-1 hole => M>=2/(n+q). Coprimality q∤n is ESSENTIAL (q=7 at n=14, the apex prime 7|14, has no q^{-1} mod 21 -- the witness FAILS there). (B) THE WITNESS GLIDES (multi-metric ladder): over multiples mq the binding modulus glides n+q -> mq+1 with M = m/(mq+1) -> 1/(n-1), infimum 2/(n+q) at the smallest multiple -- the SAME Steinhaus scaling law as the covering-min large-multiple trap (HYP-3763), now in the tight regime. (C) THE MULTI-METRIC SHEAF: witnesses are local sections of a DANGER PRESHEAF over the site of moduli (radius floor(D/n) = the 'metric'); a single metric misses cases (q|n) that others cover (q=7 caught by moduli 7,9,11,...); the tight/covering set is a global obstruction; the Fejer-Bochner minorant, evaluated POINTWISE (Reynolds/group-averaging SUPPRESSED -- averaged lenses are blind), is the section certificate. (D) H7: covering-min rung a(n)=2,2,4,4,3 (n=7..11), irregular; rung 2 = mediant 2/(2n-1) = the (n+q)-witness at q=n-1; realizability = a covering-system condition, no closed form
status: MIXED. (A) (n+q)-witness PROVED (clean residue argument; coprimality condition sharp) + VERIFIED exact n=14,18,20 (M = 2/(n+q) exactly). (B) glide VERIFIED (M=m/(mq+1), n=14 q=13, m up to 100) = HYP-3763 scaling. (C) sheaf = a FRAMEWORK (precise presheaf/site/gluing; Fejer-Bochner pointwise is the analytic backbone) -- organizing, not a new theorem. (D) H7 a(n) VERIFIED n=7..11 (irregular); realizability = covering-system (conjectural no-closed-form, like Farey rung HYP-3732). The RESIDUAL (general S, huge multiple / runner n at residue -1 defeats a single (n+q)) is handled BY the sheaf (another metric catches it) but not yet by one closed statement -- the open crux.
source: klein-2026-06-30-S54
depends_on:
  - THM-523    # q-witness (radius 0) + covering-set reduction
  - HYP-3763   # the Steinhaus scaling law (the glide is its tight-regime twin)
related:
  - HYP-3764   # covering-min open edge (rung realizability = H7)
  - HYP-3741   # witness hierarchy (radius-r = the metrics)
  - HYP-3749   # the wide hole (the (n+q) missing pair is an instance)
  - HYP-2893   # Jacobsthal doubling criterion (the surviving swaps)
  - HYP-+2873  # additive energy = Fejer 4th moment (the averaged lens; here SUPPRESSED for pointwise)
  - HYP-3732   # the Farey rung (a(n) realizability)
results:
  - 04-computation/rung_realizability_nq_witness_klein.py
  - 05-knowledge/results/rung_realizability_nq_witness_klein.out
  - 05-knowledge/results/nq_witness_gluing_klein.out
---

# HYP-3765 — the (n+q)-witness, the multi-metric sheaf, and rung realizability

## (A) The (n+q)-witness, rigorous (tight/floor regime, single swap)
> **Lemma.** Let `q` be coprime to `n` with `(n-1)/2 < q <= n-1`. Let `S = {1,...,n-1} \ {q} ∪ {mq}` be a
> single swap of the tight AP (drop `q`, add a multiple `mq`, `m>=2`, `mq` small enough that `m mod (n+q)
> ∉ {0,±1}`). Then `M(S) >= 2/(n+q) > 1/n`, so `S` is NOT tight.

**Proof.** Work mod `D = n+q`. Since `gcd(q, n+q) = gcd(q,n) = 1`, `a = q^{-1} mod D` exists; take the
observer at `t = a/D`. A runner `s` lands at `s a = s q^{-1} mod D`. The residue `±1` is hit only by
`s ≡ ±q (mod D)`, i.e. `s = q` or `s = -q ≡ n (mod D)`. But `q` was **dropped**, and `n ∉ {1,...,n-1}`
is **out of range** — so neither `±1` residue is occupied by a core runner. The added runner lands at
`mq·q^{-1} = m`, and `m ∉ {0,±1} mod D` by hypothesis. Hence every runner avoids `{0,±1} mod D`:
`min_s dist(s a, 0) >= 2`, so `M(S) >= 2/D = 2/(n+q)`. Finally `q <= n-1 => n+q <= 2n-1 < 2n`, so
`2/(n+q) > 2/(2n) = 1/n`. ∎

The missing residue pair `{q, n} mod (n+q)` is an exact instance of the punctured-core wide hole
(HYP-3749). **Coprimality is sharp:** at `n=14`, `q=7` divides `n` (`7|14`, the apex prime), so `q^{-1}
mod 21` does not exist and the `(n+q)`-witness FAILS — `q=7` needs a different metric (part C). VERIFIED
exact: `n=14` (`q=11,13`), `n=18` (`q=11,13,17`), `n=20` (`q=11,13,17,19`) all give `M = 2/(n+q)`
exactly. The composite doubling `q=n-2` (`->2(n-2)`, GW) is the one swap NOT killed — tight only at
`n ≡ 2 mod 6` (HYP-2893), the sporadic structural doublings (n=14, 20 hit; n=18 misses).

## (B) The witness glides — a multi-metric ladder (= the Steinhaus scaling law)
For a larger multiple the binding modulus **glides**: `S = {1..n-1}\{q} ∪ {mq}` has (verified `n=14,
q=13`, `m` up to 100)
```
 m :   2      3      5      10        50         100
 M : 2/27   3/40   5/66  10/131   50/651   100/1301   ->  1/(n-1)
 D : n+q=27  40=3q+1 66   131      651      1301   (D = mq+1)
```
`M = m/(mq+1)`, increasing in `m` toward the loose ceiling `1/(n-1)`; the **infimum over multiples is
`2/(n+q)`** (the smallest multiple). This is the same **multiplier-scaling law** as the covering-min
large-multiple trap (HYP-3763: `kappa=kc` lands at distance `c`, `M=c/D`), here in the tight regime.
So "drop `q` => not tight" is certified by a whole *ladder* of local witnesses (one modulus per
multiple), whose infimum stays above `1/n`.

## (C) The multi-metric witness SHEAF (the organizing framework)
The witnesses are local sections of a **danger presheaf** over the site of moduli:
- **Site.** Moduli `D = 2,3,4,...`; the "metric" at `D` is the radius `r(D) = floor(D/n)` (radius 0 for
  `D<n` = resonances, radius 1 for `n<=D<2n` = pairs, radius 2 for `2n<=D<3n`, ...).
- **Danger presheaf `𝒟`.** `𝒟_S(D) = { a mod D : some runner within r(D) of 0 }` (the danger-covered
  rotations). `S` is tight (`M<=1/n`) **iff** `𝒟_S(D) =` all rotations for every `D` — the danger
  covers the whole `(D,a)` space, no lonely point.
- **A witness** = a rotation `a` at `D` with `a ∉ 𝒟_S(D)` = a *local failure of the danger section* = a
  certified lonely point (`M >= (r(D)+1)/D`). The `q`-witness (radius 0, `D=q`), the `k`-witness (radius
  1), and the `(n+q)`-witness (radius 1, `D=n+q`) are three such local sections.
- **Gluing (CRT / Steinhaus glide).** A runner that defeats the witness at `D` (e.g. a huge multiple, or
  the runner `n` sitting at residue `-1`) is chased to another modulus `D'` by CRT — the binding glides
  (part B). No single metric suffices: `q=7` (apex, `q|n`) escapes `(n+q)` but is caught at moduli
  `7,9,11,...` (`M = 1/7, 1/9, 1/11, ...`, verified). The tight/covering set is the **global
  obstruction**: a lonely point covered by no local danger.
- **Fejer-Bochner certificate, POINTWISE (Reynolds averaging suppressed).** Each danger indicator
  `1[within r of 0 at (D,a)]` is minorized by a positive-definite (Bochner) trigonometric polynomial —
  a Fejer/Beurling-Selberg kernel — and the witness is the point where this minorant, evaluated
  **pointwise** at `t=a/D`, is zero. We use the pointwise value, NOT the group-averaged (Reynolds
  operator `∫_0^1`) quantity: the average is blind (additive energy / spectral 4th moment does not
  separate tight from non-tight, HYP-+2873/opus-S2). Tightness is a POINTWISE property; the sheaf
  sections are pointwise minorant evaluations.

## (D) H7 — rung realizability of the covering-min
The covering-min rung `a(n) = 2,2,4,4,3` (`n=7..11`), binding `D = a(n-1)+1 = 13,15,33,37,31`. **Rung 2
= the mediant `2/(2n-1)`**, and `2n-1 = (n+q)` at `q=n-1` — so *the covering-min's tightest rung is the
`(n+q)`-witness value at `q=n-1`*: the two regimes (floor and covering-min) meet at rung 2. Rung `r` is
realizable iff a covering set (a multiple of every `q∈{2..n}`) has all speeds avoiding `{0,±1,...,±(r-1)}
mod (r(n-1)+1)` at some rotation — a **covering-system condition** on the factorization of
`r(n-1)+1`. This is realizable at `n=7,8` (rung 2) but not `n=9,10,11` (`a=4,4,3`); the pattern is
irregular with no clean prime rule — predicted to have **no closed form**, like the Farey rung
(HYP-3732), `width(G_n)`, and `A000568`.

## Net
The `(n+q)`-witness is now a clean rigorous theorem (single-swap tight regime): dropping a `q` coprime to
`n` in the top half, the pair `{q,n}` vacates the `±1`-neighbourhood mod `n+q`, forcing `M>=2/(n+q)>1/n`;
only the composite doubling `n-2` survives (GW, `n≡2 mod 6`). Larger multiples make the witness glide up
a ladder (`M=m/(mq+1) -> 1/(n-1)`, the Steinhaus scaling of HYP-3763). The witnesses organize as local
sections of a multi-metric danger sheaf over the modulus site, glued by CRT, certified by pointwise
Fejer-Bochner minorants (not group averages); the residual (one metric fails when `q|n` or the multiple
is huge) is exactly what the gluing repairs. H7: rung realizability is a covering-system function,
irregular, meeting the `(n+q)`-witness at rung 2.
