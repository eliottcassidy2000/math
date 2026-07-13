---
id: THM-735
title: The far-peel disc_v certificate, proved via Farey / three-gap — for the interval-core single-killer covering family S = {1,…,12, 182m} (m ≥ 1), in particular the DEEP WELL (m=1, the covering-min), the loneliness L(S) > 0, hence M(S) ≥ 1/14 (LRC(14) holds for S). This is the FIRST time klein's disc_v certificate (THM-731) is turned into a THEOREM (not "verified"): the missing analytic upper bound on disc_v is supplied at the FAR peel by the exact Farey-arc (three-gap) structure of SafeSet({1..12}) — N = φ(13) = 12 arcs, |G'| = H₁₂/91 — and the whole certificate collapses to H₁₂ > √2.
status: PROVED (rigorous; elementary modulo klein's THM-731 inequality, which is itself rigorous). Scope: the interval-core single-killer family (deep well + far-dilates) — the covering-min extremal, opus-S270's "certificate-poorest / binding" peel. NOT all covering families: general (non-consecutive) cores need the general Farey/three-gap count and are the remaining work. Verified end-to-end numerically (mac-mini-S88): |G'|=H₁₂/91 exact, N=12 exact, disc₁₈₂=4.13e-4 ≤ bound 1.45e-3, L_cert=0.0159>0.
source: mac-mini-2026-07-13-S88 (owner: "prove klein's disc_v bound at the far peel via three-gap")
depends_on:
  - THM-731   # klein's peeling certificate L=(6/7)|G'_~v|-ε_v, |ε_v|²≤(6/49)disc_v (RIGOROUS); this THM supplies its open disc_v bound at the far peel
  - THM-732   # opus's disc_v ≤ r²/(3v²) jump-bound form (r = #intervals); this THM supplies the rigorous r=φ(13)=12 via Farey
  - THM-724   # single-killer covering-min rigidity (independent balance proof that the deep well is 1/14-lonely; this is the disc-route proof)
related:
  - HYP-6530  # the #far reduction + far-peel three-distance pinpoint that set this up
  - HYP-2566  # the covering-min = LRC(14); this proves its extremal family, not the whole
external: LRC(≤13) SETTLED. Steinhaus three-gap / Farey dissection (classical).
---

# THM-735 — The far-peel disc_v certificate, proved via Farey / three-gap

**One line.** klein's certificate (THM-731) says `L ≥ (6/7)|G'_{~v}| − √((6/49)·disc_v)` for every
peel `v`, rigorously, **except** for an analytic upper bound on `disc_v` (the one piece THM-731
lists as "verified-not-proved"). At the **far peel** of the interval-core single-killer family
`S = {1,…,12, 182m}`, that bound is supplied **exactly** by the Farey/three-gap structure of the
leftover good set `SafeSet({1..12})`, and the certificate reduces to the trivial inequality
`H₁₂ > √2`. Hence `L(S) > 0`, i.e. **`M(S) ≥ 1/14`** — in particular for the **deep well**
`{1,…,12,182}`, the covering-min extremal.

## Setup (klein's peel, THM-731)

`c = 1/14`; danger `D_w = {t : ‖wt‖ < c}`, `|D_w| = 1/7`. Peel the far element `v = 182m`; the
leave-one-out good set is `G' := G'_{~v} = ⋂_{w∈{1..12}} D_w^c = SafeSet({1,…,12})` (already
inside the "middle" because `D_1` removes `[0,c]∪[1−c,1]`). THM-731 gives, **rigorously**
(Cauchy–Schwarz + Wiener–Khinchin + Poisson `v`-grid sampling):
$$L \;=\; \tfrac67|G'| - \varepsilon_v,\qquad |\varepsilon_v|^2 \le \tfrac{6}{49}\,disc_v,\qquad
disc_v = \textstyle\sum_{m\ne0}|\hat 1_{G'}(mv)|^2 \ge 0,$$
so `L ≥ (6/7)|G'| − √((6/49)·disc_v)`. It remains to **bound `disc_v` from above** and show the
right side is positive. That is the whole content below.

## Step 1 — The Farey-arc (three-gap) structure of `SafeSet({1..12})`

**Lemma (Farey dissection of the consecutive safe set).** For `c = 1/14`,
`SafeSet({1,…,k}) = {t : ‖it‖ ≥ c, i=1,…,k}` is the disjoint union, over pairs of **adjacent**
Farey fractions `a/b < c'/d` in `F_k` (so `bc'−ad = 1`, `c'/d − a/b = 1/(bd)`), of the arcs
$$\Big(\tfrac ab + \tfrac{c}{b},\ \tfrac{c'}{d} - \tfrac{c}{d}\Big),\qquad\text{of length }
\ \tfrac1{bd}-\tfrac cb-\tfrac cd=\frac{1-c(b+d)}{bd}=\frac{14-(b+d)}{14\,bd},$$
**which is present iff `b+d ≤ 13`.**

*Proof.* The danger set `⋃_{i≤k} D_i` has teeth centred exactly at the fractions `{n/i : i≤k}`,
i.e. at the points of `F_k`. The widest tooth at a reduced fraction `a/b` comes from `D_b`
(half-width `c/b`); the teeth from `D_{2b}, D_{3b}, …` sit inside it. Between Farey **neighbours**
`a/b, c'/d` there is, by definition of adjacency, **no** point of `F_k`, hence no tooth centre; and
a computation (using `b,d ≤ 13`) shows neither neighbour's own tooth, nor any tooth centred outside
`[a/b, c'/d]`, reaches into the open gap beyond `a/b+c/b` on the left and `c'/d−c/d` on the right.
So the gap is safe exactly on the stated arc, of the stated length, positive iff `c(b+d)<1`. ∎

**Specialisation `k = 12`.** Adjacent fractions in `F_k` always satisfy `b+d ≥ k+1` (else their
mediant `(a+c')/(b+d) ∈ F_k` would lie between them). For `k = 12` this is `b+d ≥ 13`, so combined
with `b+d ≤ 13` the safe arcs occur **exactly at `b+d = 13`**. Each fraction of denominator `13`
(there are `φ(13)=12` of them, all reduced since 13 is prime) is the mediant of exactly one such
adjacent pair, so:
$$\boxed{N = \varphi(13) = 12\ \text{arcs},\qquad
|G'| = \sum_{\substack{b+d=13}} \frac{1}{14\,bd}
= \frac{1}{14}\sum_{b=1}^{12}\frac{1}{b(13-b)} = \frac{2H_{12}}{14\cdot 13} = \frac{H_{12}}{91}.}$$
(The last identity uses `1/(b(13−b)) = (1/13)(1/b+1/(13−b))` and `Σ_{b=1}^{12}(1/b+1/(13−b)) =
2H_{12}`.) Numerically `|G'| = 6617/194040 = 0.0341012`, and both `N = 12` and this measure are
confirmed by direct grid computation (mac-mini-S88).

## Step 2 — The Fourier jump bound on `disc_v`

`1_{G'}` is a sum of `N` interval indicators, so it has `2N` jumps `±1` at the arc endpoints
`x_p`, and for `ℓ ≠ 0`
$$\hat 1_{G'}(\ell) = \frac{1}{2\pi i \ell}\sum_{p} s_p e^{-2\pi i \ell x_p},\qquad
|\hat 1_{G'}(\ell)| \le \frac{2N}{2\pi|\ell|}=\frac{N}{\pi|\ell|}.$$
Hence (this is opus's THM-732 form, `r = N` intervals, now with `N` known exactly):
$$disc_v = \sum_{m\ne0}|\hat 1_{G'}(mv)|^2 \le \frac{N^2}{\pi^2 v^2}\sum_{m\ne0}\frac1{m^2}
= \frac{N^2}{\pi^2 v^2}\cdot\frac{\pi^2}{3} = \frac{N^2}{3v^2}
\ \le\ \frac{12^2}{3\cdot 182^2}=\frac{12}{8281}\quad(v=182m,\ m\ge1).$$
(The true value is `disc_{182} = 2629220219/6363107150400 = 4.13\times10^{-4}`, well under the
bound `1.45\times10^{-3}`; klein/opus, HYP-6510.)

## Step 3 — The certificate collapses to `H₁₂ > √2`

`L_cert > 0 ⇔ (6/7)^2|G'|^2 > (6/49)disc_v ⇔ 6|G'|^2 > disc_v`. Using Steps 1–2:
$$6\Big(\frac{H_{12}}{91}\Big)^2 \;>\; \frac{12}{8281}.$$
Since `182 = 2·91`, `182^2 = 4·91^2`, so `12/8281 = 12/91^2`, and the inequality is
$$\frac{6H_{12}^2}{91^2} > \frac{12}{91^2}\iff 6H_{12}^2 > 12 \iff H_{12}^2 > 2 \iff H_{12} > \sqrt2.$$
`H_{12} = 3.10321 > 1.41421 = √2` (indeed `H_{12}^2 = 9.63 > 2`, margin `4.8×`). Therefore
$$L(S) \ \ge\ \tfrac67\cdot\tfrac{H_{12}}{91} - \sqrt{\tfrac{6}{49}\cdot\tfrac{12}{8281}}
\ =\ 0.02923 - 0.01332 \ =\ 0.01591 \ >\ 0.$$
So `L(S) > 0` and **`M(S) ≥ 1/14`** for every `m ≥ 1`. ∎

## General interval-core form (why 12 is the tight case)

For a consecutive core `{1,…,k}` the same three steps give `disc ≤ k²/(3v²)` with `v` a multiple of
`lcm(k+1,…,14) ≥ 14(k+1)` (the far element must carry all missing moduli), and — when `k+1` is
prime — `|G'| = H_k/(7(k+1))`, so the certificate reduces to
$$H_k \;>\; \frac{k}{6\sqrt2}\qquad(\text{equivalently } H_k^2 > k^2/72).$$
This holds for **all `k ≤ 12`** (margins `3.5× → 2.2×` for `k=6→12`; mac-mini-S88), so every
interval-core single-killer family is 1/14-lonely by the far peel. Since single-killer **forces**
`k = 12` (one outlier ⇒ `|core| = 12`), the theorem's family `{1,…,12,182m}` is the complete
interval-core single-killer class, and `k=12` is the extreme (thinnest `|G'|`, tightest margin) —
consistent with the deep well being the covering-min.

## Honest scope

- **Proved:** `L > 0` (hence LRC(14)) for the interval-core single-killer family — the covering-min
  extremal and its far-dilates — via klein's certificate with a **rigorous** analytic `disc_v`
  bound at the far peel. This is the first theorem-grade instance of THM-731's open piece.
- **What is new vs THM-724:** THM-724 already proves the deep well is 1/14-lonely (via balance).
  The novelty here is the **route**: it validates klein's disc_v certificate program by proving its
  one "verified-not-proved" ingredient in the cleanest case, with exact Farey closed forms
  (`N = φ(13) = 12`, `|G'| = H_{12}/91`) and a one-line reduction (`H_{12} > √2`).
- **Not proved:** the general covering family (non-consecutive core). There the good set is still a
  Farey/three-gap arc union, but the count `N` and measure `|G'|` are governed by the general Farey
  dissection at level `1/14`; bounding `disc_v ≤ N²/(3v²) < 6|G'|²` there is the remaining step
  toward the full covering `L > 0` = LRC(14). The template is exactly Steps 1–3.

*Artifacts:* `04-computation/lrc14_farpeel_farey_proof_macmini_S88.py` (+`.out`) — exact Farey arcs,
closed-form `|G'| = H_{12}/91`, `N = 12`, `disc` bound, certificate `L_cert = 0.0159 > 0`, and the
`H_k > k/(6√2)` table. Credits: klein-S287 (THM-731 certificate + the rigorous `|ε|²≤(6/49)disc`),
opus-S270/271 (THM-732 jump-bound form; "only the far peel certifies"), mac-mini-S87 (HYP-6530:
the far peel = three-distance consecutive core pinpoint), THM-724 (independent balance proof).
