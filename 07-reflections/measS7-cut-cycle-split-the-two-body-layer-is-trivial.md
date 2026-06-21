# The cut/cycle split of measS7: the 2-body layer is trivial, extremality is many-body

**Author:** mac-mini-2026-06-20 (opus, THREAD E)
**Status:** PROVED (the split is an exact identity; the 2-body structure is exact) + VERIFIED (extremality scans, exact Fractions, k=8,9,10)
**Companion scripts:**
- `04-computation/lrc14_routeE_cutcycle_split_measS7_opus_s2.py`
- `04-computation/lrc14_routeE_2body_layer_structure_opus_s2.py`
- `04-computation/lrc14_routeE_2body_extremality_scan_opus_s2.py`
**Results:** `05-knowledge/results/lrc14_routeE_*_opus_s2.out`
**Depends on / twins:** THM-559 (c3 = 2-body cut-space Ising), THM-560 (OCF odd degree ladder),
THM-555 (cut/cycle wall), HYP-2708 (death-chain transfer operator), the Fourier-on-Z/7 identity
(mac-mini-0620 `lrc14_fourier_z7_resonance`).

---

## The question (THREAD E)

THM-559 split the tournament datum c3 into a **cut-space 2-body Ising energy** (regular tournament
= ground state, provably extremal) plus a **cycle-space many-body** remainder. The LRC(14) cover
`measS7(E)` is the offset-side analog (a Z/7 object, not a Z/2 object). Does measS7 admit the same
cut/cycle split, with the **cut part already consec-extremal by a 2-body argument** and the cycle
correction bounded?

## The split (EXACT IDENTITY)

Using the proved Fourier-on-Z/7 identity `measS7(E) = sum_a Khat(a) J(E,a)` with
`J(E,a)=int_0^1 prod_e omega^{a_e c_e(x)}dx`, `c_e(x)=floor(7 frac(ex))`, define the **single-clock
product (decorrelated) moment** `J_cut(E,a)=prod_e m(a_e)` where `m(a)=int_0^1 omega^{a c_1(x)}dx`.

> **Apex-prime vanishing (verified exact):** each clock is *individually* uniform on Z/7, so
> `m(a)=[a=0]`. Hence `J_cut(E,a)=[a=0]` and
> $$\boxed{\ \text{measS7}(E) \;=\; \underbrace{\widehat K(0)}_{\text{CUT}} \;+\; \underbrace{\sum_{a\ne 0}\widehat K(a)J(E,a)}_{\text{CYCLE}=\,\mathrm{corr}(E)} \;=\; \mathrm{iid}_k + \mathrm{corr}(E).\ }$$
> The CUT part is `Khat(0)=7!\,S(k,7)/7^k = iid_k`, the iid surjection probability — **independent of E**.

So the brief's "wide-bound mechanism" (`measS7 = iid_k + corr`) IS the cut/cycle split. The cut part
is the single-particle / decorrelated piece (the leading mode of HYP-2708's death-chain operator).

## The 2-body cycle layer is an EXACT pairwise energy — but trivial

Split the cycle part by Fourier weight `t = #nonzero coords of a`: `corr = sum_{t>=1} corr_t`.

> **(EXACT, verified additivity)** `corr_2(E) = sum_{pairs {i,j}} P_k(e_i,e_j)` — a genuine **2-body
> energy on the offsets, zero external field**, the LRC twin of THM-559's line-graph Ising.
>
> **(EXACT, verified)** The per-pair coupling `P_k(e_i,e_j) = 0` unless `e_i,e_j` are both nonzero
> and `e_i \equiv e_j \pmod 7` (a **resonant** pair), and `P_k < 0` there. Non-resonant pairs are
> **exactly pairwise-independent**: `P(c_i = c_j) = 1/7` exactly. Resonant pairs over-agree:
> for a pair `(e, e+7)` with smaller speed `e\in\{1..7\}`,
> $$P(c_i=c_j) = \tfrac{2}{e+7}\quad(\text{closed form, exact}),\qquad \text{excess}=\tfrac{7-e}{7(e+7)}>0.$$
> The negative coupling is the antiferromagnetic penalty for this excess agreement (two clocks that
> agree waste coverage of Z/7). `P_k(\text{resonant})` is `-7.6e-4, -2.3e-3, -4.2e-3, -5.9e-3` for
> `k=7,8,9,10`.

This is the precise, rigorous THM-559 analog: a 2-body antiferromagnetic Ising energy whose
couplings are zero on the "decorrelated" (non-resonant) pairs and negative on the resonant pairs.

## The honest verdict: extremality is NOT 2-body — it is many-body

Here is where the analogy **breaks**, instructively. For c3 the 2-body cut energy carries *all* the
extremality (regular = ground state). For measS7 the 2-body layer is **trivial**:

- **k <= 8:** `corr_2`-max `= 0`, attained by every **7-gap-free** shape (<=1 offset per residue
  class mod 7). consec_8 is one of 64 such shapes — so consec **ties** for the 2-body max
  (vacuously extremal: `|corr_2| <= 0.0052`, max value `0`). PROVABLE by the per-pair ground-state
  argument (all couplings `<=0`; consec realizes the all-zero state).
- **k >= 9:** within span `<=13`, the nonzero offsets `{1..13}` occupy residues
  `{1,2,3,4,5,6,0,1,2,3,4,5,6}` mod 7; a 7-gap-free set has `<=7` nonzero offsets, so **`k<=8`**.
  For `k>=9` **EVERY** span-`<=13` shape has a resonant pair, so `corr_2 < 0` for all, and **consec
  is NOT the 2-body max** (598/1287 shapes beat it at k=9; the least-penalized 2-body shape is
  `[0..7,13]`, not consec). consec carries resonant penalties `(1,8)`, `(2,9)`,... that the 2-body
  layer alone *dislikes*.

So the deliverable is a **clean cut/cycle split with the 2-body cut/decorrelated piece exactly
characterized and proved consec-extremal only in the degenerate (constant cut + zero-max 2-body)
sense — and an honest refutation that the 2-body layer drives consec's extremality.** The genuine
advantage of consec lives in the **odd / many-body cycle layers**: over spread8, consec's `+0.139`
lead is carried by `corr_3 (+0.007 rel)`, `corr_5 (+0.025 rel)`, and the `t>=6` residual
`(+0.092 rel)`; `corr_2` contributes essentially nothing.

## Why this matters (the cross-theorem resonance)

- **The cut/cycle seam is universal but the WEIGHT shifts.** THM-559: c3's action is 2-body (cut).
  THM-560: the OCF's odd-cycle degree ladder pushes the action to high (odd) body order (cycle).
  measS7 sits on the **OCF/cycle side**: its 2-body layer is inert and its extremality is irreducibly
  many-body and **odd-dominated** (corr_3, corr_5 carry the signal; corr_4 is tiny). This is the
  exact LRC echo of THM-560's "oddness IS the degree drop / magic onset".
- **Confirms THM-555's verdict in a new currency.** "Cut cheap, cycle dear" becomes "the 2-body
  (cut) layer of measS7 is `O(0.005)` and flat; everything that decides consec lives in the cycle
  space." The cut/cycle split does NOT give a 2-body handle on the cover's extremality — and now we
  know *precisely why*: the resonance that creates 2-body couplings (`7 | gap`) cannot even be
  avoided once `k>=9` within the operative span, so the 2-body layer turns *against* consec.
- **It re-derives, from the Ising side, the project's standing negative results** (consec-max is
  irreducibly aggregate; per-block/monotone/convexity all dead): the only *clean* 2-body extremal
  statement that survives is the vacuous `k<=8` one.

## One genuinely new, exact object

The closed form `P(c_i=c_j)=2/(e+7)` for a resonant pair `(e,e+7)` and the exact pairwise
independence `P(c_i=c_j)=1/7` for every non-resonant pair are new exact facts. They give the
**pairwise** layer of the apex-prime decorrelation (HYP-2708 is the *marginal/1-body* layer; this is
the *2-body* layer), with an explicit, signed, antiferromagnetic coupling table — usable as the
2-body term in any future variational/Ising attack, even though that term alone does not close LRC.
