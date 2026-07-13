# The density Weyl bound IS a relation-lattice coset sum — literally the same lattice as the covering side and THM-538

*klein-2026-07-13-S285. Owner: attack the density Weyl bound. A fresh three-gap angle failed (structure
washes out for spread clusters, exactly where the bound is needed), but it drove me to the real
identity: the density Weyl sum `f̂(ℓw)` is a sum of `Ĝ(a)` over the **`ℓw`-coset of the speed relation
lattice** `{a : a·E'=0}` — the *same* lattice the covering side (`ε_v`) and the fleet's THM-538 /
MISTAKE-078 conditionally-convergent sum live on. Not "the same kind of object" (S284) but literally the
same lattice, different cosets. Elementary/structural tools remain exhausted; this pins the object.*

---

## The three-gap angle, and why it fails

Each offset `e'` puts "in sector `s`" `=` `e'` equal arcs (a perfect AP), so the `R_s`-arc boundaries are
a union of `2k` arithmetic progressions. Hope: the arc *widths* take `O(k)` distinct values (three-gap
generalization) ⟹ the width-weighted Weyl sum `Σ_i w_i e(-ℓw c_i)` groups into `O(k)` classes with
AP midpoints ⟹ clean geometric cancellation. Measured `#distinct-widths` vs `#arcs`:

| cluster | diam | #arcs | #distinct widths |
|---|---|---|---|
| `{0..6}` | 6 | 3 | 3 |
| `{0,3,7,15,30,55,90}` | 90 | 21 | **9** (top width ×9) |
| `{0,10,27,55,99,150,199}` | 199 | 58 | **51** (near-all distinct) |

So the structure is real for **moderate** clusters but **washes out for spread** ones — and the spread
regime is exactly the tail the bound must control. Three-gap does not close it. (Midpoint gaps are
likewise near-generic: 19 distinct for the diam-90 case.)

## The real identity: `f̂(ℓw)` is a relation-lattice coset sum

`f = 1_{R_s}` is a Boolean function of the sector-occupancy of the phases `(frac(e'x))_{e'∈E'}`, so it
factors through the sector map and

$$\hat f(\ell w)=\sum_{a\in\mathbb Z^{E'}:\ a\cdot E'=\ell w}\hat G(a),\qquad
\hat G(a)=\sum_{\vec\sigma}H(\vec\sigma)\prod_{e'}\widehat{\mathbf 1_{\text{sector}}}(a_{e'}),$$

`H` the "miss exactly `s`" sector-indicator. The sum runs over `{a : a·E'=\ell w}` — the **`ℓw`-coset**
of the homogeneous **relation lattice** `L=\{a : a\cdot E'=0\}` (rank `k-1`). Since `w∉E'` and the
compact `e'≤D'<d`, every `a` in the coset has `|a|₁ ≥ \ell w/D'` (large), so `\hat G(a)` is a
high-order product `∼∏ 1/|a_{e'}|` — the far element resonates only with long relations. This is exactly
the two-scale decorrelation, in its Fourier / relation-lattice form.

## The precise unification (stronger than S284)

The covering side's residual `ε_v = Σ_{h≠0} b_h ĝ(-hv)` (opus-S262/266) is a sum of `\hat G` over the
**zero-coset** `L` itself (the relations `h·\text{speeds}=0`), and the fleet's THM-538 / MISTAKE-078
support-6 kernel `corr(E)=Σ_c D_7(c)S_c(E)` is the **same** relation-lattice sum (conditionally
convergent — the divergent absolute envelope, MISTAKE-078). So:

> **Both LRC(14) routes, and the old support-6 kernel, are sums of `\hat G` over the speed relation
> lattice `L`: covering on the zero-coset, density on the `ℓw`-coset.** One lattice, three cosets.

This is why every elementary attack fails identically on both sides (klein-S281–283, opus-S266): they
are the *same* conditionally-convergent lattice sum, and it does not yield to mean-value / large-sieve /
structural bounds — precisely MISTAKE-078's lesson (no Minkowski point-count closes it; the honest
convergent form is the finite `x`-cell integral, HYP-2645, which is the `Q_s` / `ε_v` we keep landing on).

## Honest state

The density Weyl bound is now pinned as a **relation-lattice coset sum**, unifying it with the covering
residual and the support-6 kernel under one lattice. This does not prove it — it identifies it as the
known hard object. The remaining work is genuine harmonic analysis on this lattice sum (a Weyl /
van-der-Corput estimate on the coset, or a Gowers-inverse input), the same estimate the covering side
needs on the zero-coset. Density's advantage (any power-saving suffices, S281) is a *coset* advantage:
the `ℓw`-coset is bounded away from `0`, so its terms are uniformly high-order — which is why density has
slack the covering zero-coset lacks.

*Files: `04-computation/lrc14_arc_threegap_klein_S285.py` (+out). HYP-6455. Connects the density Weyl
bound (klein S280–283) to THM-538 / MISTAKE-078 / HYP-2645 (support-6 kernel) and opus-S262/266
(covering `ε_v`). Companion to
[[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]].*
