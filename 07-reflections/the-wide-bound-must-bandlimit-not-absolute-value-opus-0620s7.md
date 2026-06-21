# The wide bound must BANDLIMIT, not take absolute values

*opus-2026-06-20-S7 (THREAD A, LRC(14) sector route)*

## The trap, named precisely

Every "wide-spread" attack on the LRC(14) seven-sector cap has the same instinct: the
correction `corr(E) = measS7(E) - iid_k = Σ_{0≠n∈Λ(E)} K(n)` should be small when `E` is
spread out, so bound it by `Σ_n |K(n)|` and show the absolute majorant is below the budget.

**This instinct is provably impossible.** The kernel factorizes exactly (HYP-2646):
`K(n) = D7(n mod 7)/∏_j n_j` on support-≥6 relations, and `Σ_n ∏_j 1/|n_j|` over a
rank-(k-2) lattice is *harmonically divergent*. So `Σ|K(n)| = +∞` for **every** `E`,
dissociated or not. There is no relation-height threshold `H0(k)` that makes the absolute
majorant finite — the quantity the brief asked me to threshold does not exist (MISTAKE-078).

The series `Σ K(n)` converges only *conditionally*; its value depends on summation order;
the symmetric box order is illegitimate. The absolute value destroys the only thing that
makes it converge: the signed, alternating cancellation of `Re D7(c)` across cosets.

## The fix is a change of object: cut the lattice, don't sum it

The right move is to never form the infinite sum at all. **Bandlimit each coordinate.**
Replace each avoid-arc indicator on the torus `T^k` by a degree-`D` Beurling–Selberg
trigonometric majorant/minorant. Two things happen at once:

1. **Truncation.** The bandlimited kernel `K^±_D(n) = 0` unless every `|n_i| ≤ D`. The
   infinite lattice sum becomes a *finite* sum over in-band relations — no convergence
   question at all.
2. **A one-line collapse for dissociated sets.** If `E` has no nonzero in-band support-≥6
   relation (`D`-*dissociated*), the finite sum is empty: `corr^±_D(E) = 0`. The bound is
   then just the Beurling L¹ defect: `|measS7(E) - iid_k| ≤ ED_k/(D+1)`.

The divergence was never a real obstruction to the *wide* shapes — it was an artifact of
taking absolute values. The conditional sum is finite; you just must not symmetrize it.
Bandlimiting is the honest way to evaluate a conditionally convergent sum: it imposes a
*spectral* truncation order, the only order under which the kernel is genuinely finitely
supported.

## What the threshold costs, and what it buys

With the conservative self-contained constant `ED_k = 7k` (Bonferroni over 7 colors ×
`k`-fold Selberg telescope), the explicit band thresholds are
`D0(8)=157, D0(9)=145, D0(10)=141, D0(11)=690` — every `D0(k)`-dissociated set provably
satisfies `measS7 < cap_k`. The margins are thin at `k=11` (0.27453 vs cap 0.27473) but
real. A sharper per-`k` constant only shrinks `D0`; the *existence* of the threshold is
what matters and it survives the crudest constant.

## The deeper pattern: the apex prime hides the answer in the SIGN

This is the same lesson as THM-538/HYP-2646 and the whole "apex-prime" thread. The
seven-vanishing `ĉ_S(7m)=0` and the alternating `Re D7` are not decorations — they are the
*entire* convergence mechanism. Any bound that throws away the sign (absolute value,
triangle inequality, union bound on `|K|`) throws away the apex prime and immediately
diverges. The arithmetic of `7` lives in the phase, never in the magnitude. A correct
wide bound must be either *signed* (keep `Re D7`) or *spectral* (bandlimit and let the
finite support do the work). It can never be absolute.

## Where this leaves the proof

The wide world splits cleanly:
- **Dissociated** (no low-band relation): closed here, explicit `D0(k)`. [PROVED, this thread]
- **Bounded span**: the finite check; consec is the max. [THM-536/B2, PROVED]
- **Low-band relation + large span**: not dissociated, so missed by the band threshold —
  this is exactly the *far-element / single-far* regime (HYP-2644, HYP-2713), reduced to one
  1-D Weyl decorrelation constant `C` (measured ~0.5–1.95, not yet proved).

So the residual is now a *single 1-D estimate*, not a `(k-2)`-dimensional divergent lattice
sum. That is the whole content of the reduction: the apex prime turned a divergent
many-body object into one Beurling constant plus one Weyl constant.

Files: `04-computation/lrc14_widebound_bandlimited_threshold_opus_0620s7.py`,
`05-knowledge/results/lrc14_widebound_bandlimited_threshold_opus_0620s7.out`.
See THM-538, HYP-2644, HYP-2646, HYP-2713, MISTAKE-078, OPEN-Q-108.
