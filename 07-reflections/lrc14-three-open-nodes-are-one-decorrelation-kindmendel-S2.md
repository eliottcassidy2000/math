# The three open LRC(14) nodes are one problem: effective decorrelation

*kind-mendel-2026-06-22-S2. Follows the S1 verification audit (HYP-2847). Built from a
four-thread repo sweep (Nodes 1/2/3 + a whimsical cross-domain hunt) and two exact
computational levers (`04-computation/lrc14_node23_levers_kindmendel.py`).*

## The unification

After the S1 audit isolated the three open nodes, I tried to attack them separately and
they kept collapsing into the same shape. State each as a deviation-from-independence:

- **Node 1 (finite-V_max Part A).** In the V→∞ limit, at slow time `x` the whole cluster
  rigidly rotates by the fast phase `φ=frac(Vt)`, and a good `φ` exists iff the cluster
  config `{frac(d_i x)}` has a circular gap `>1/7` (= GOOD). The gap is that for finite `V`
  the fast rotation and the slow set are coupled; the error is the **discrepancy of the
  rotation `frac(Vt)` sampling the slow set `GOOD∩G_P`** — an off-diagonal term that
  vanishes as `V→∞`.
- **Node 2 (consec maximizes L_y / `p0≤cap`).** Spreading the cluster out makes its phases
  decorrelate, so `p0` (the all-sectors-hit measure) **drops toward its independent-value
  limit**. consec is the maximally-resonant (least decorrelated) arrangement. So the cap
  bound is "resonant ≤ cap" + "everything wider decorrelates below it."
- **Node 3 (witness floor).** `meas(coverSet∩G_P) − p0·measGP` is exactly the
  **correlation** between "cluster dense" and "small-part safe"; the Bonferroni bound is
  lossy precisely because it pretends this correlation is maximal when it is nearly zero.

All three are *the same inequality*: a fixed slow object (`G_P`, or the safe set) is
**decorrelated** from a fast/cluster object (`frac(Vt)`, or `coverSet`, or the spread-out
phases), with an explicit error. Four independent search threads converged on the same
toolbox — Fourier-at-comb-frequencies, Weyl/Erdős–Turán, Vitali covering, Riesz products,
Schur-smoothing — which is the tell that decorrelation is the right coordinate, not three
separate miracles. (This is the same "two-agents-one-reformulation" signal as kps-S4.)

## The clean Fourier statement of the common core

Write `c = 1_{coverSet}`, `g = 1_{G_P}`. Then
> `meas(coverSet∩G_P) = p0·measGP + Σ_{n≠0} ĉ(n) ĝ(−n)`.
- `ĝ(n)` (Fourier of `G_P`, a union of arcs around `j/p`) **decays like 1/n** and is
  supported on frequencies that are **multiples of the small-part speeds `p∈P`**.
- `ĉ(n)` (Fourier of the seven-sector cover) lives on the **cluster relation lattice** —
  integer combinations `Σ m_i e_i` of the large cluster offsets.
The correlation is therefore a sum over the **intersection of two spectra**: small-part
multiples ∩ cluster resonances. When the cluster speeds are large and arithmetically
independent of `P`, that intersection is sparse and the correlation is small ⇒
quasi-independence ⇒ the floor. **Node 1 is the same identity with `c=1_{frac(Vt)∈good φ}`**
(spectrum on multiples of `V`), and the off-diagonal sum is `O(arcCount/V)`. **Node 2 is the
diagonal**: as the `e_i` spread, the relation lattice coarsens and `ĉ(n)→0` for `n≠0`, so
`p0→p0·1` (independent value), which the data shows is *below* the resonant consec value.

The deep tension (why it is not free): the **covering condition** (THM-523: a counterexample
must contain a multiple of every `q∈{2,…,14}`) *forces* shared arithmetic between cluster and
`P`, which is exactly what creates resonances (HYP-2842: composite `P` kills the exact Farey
points). So decorrelation must be proven to **survive the covering-induced resonances** — and
the resonant-neighborhood repair (kps-S32) is precisely the patch at those resonances.

## What the data says (this session, exact rationals)

**Lever A — the binding region is small (Node 2).** Margin `cap_k − max_E p0(E)`:
`k=8: +0.054, 9: +0.078, 10: +0.100, 11: +0.144, 12: +0.212, 13: +0.324` (monotone up).
And **wide configs have strictly lower `p0`** (one-far at k=8: `0.194 ≪ 0.327` consec). So:
the only tight cases are **k=8,9,10**, and "spreading decreases `p0`" is visible decorrelation.

**Lever B — the floor is robust but not tight to independence (Node 3).** The
quasi-independence ratio `R' = meas(coverSet^c∩G_P)/(measGP·(1−p0))` sits in `[0.66, 1.27]`
across tested configs (not ≈1, but **bounded below**); the resulting floor is `0.20–0.37`,
always `≫ m_P=0.056` and far above the Bonferroni `0.05–0.14`.

## Sharpest concrete attack per node (avoiding known dead ends)

**Node 1 — bounded arc-count + Erdős–Turán.** Prove `meas(actual good t) ≥ witnessG2 −
arcCount/V` via three-distance/Erdős–Turán for the rotation `frac(Vt)` over the union of
arcs `GOOD∩G_P`. `LRCArcComplexity` already gives the cover arc-count `≤ 7·ΣE` (no V-dep,
sorry-free). **The one missing inequality:** along covering sequences, `V ≫ arcCount`
(`= 7·Σd_i` for the gap set), i.e. the large speed dominates the cluster spread. This is a
*Diophantine* fact the covering condition should supply; it degrades only in the boundary
core `{t,2t,…,12t,V}`, so target that core specifically (the spread there is `Σd_i ~ 78t`;
need `V/t→∞`). This is the most under-attacked node and the highest leverage.

**Node 2 — binding/slack split + a crude wide bound.** (i) k=8,9,10 bounded-spread:
exhaustive exact check (consec wins) → a *finite certificate* (formalizable in interval/exact
arithmetic, like the existing `LRCPeriodmaxCertificate`). (ii) k=11,12,13: margin ≥0.144, so
a **crude provable bound suffices** — the absolute Fourier/covolume bound rejected at k=8 for
being 5.9× lossy may already beat `cap` here; worth checking. (iii) wide E at any k:
"spreading decreases `p0`" = a single-far/two-far decorrelation step; THM-563 already proves
the single-far signed deviation is exactly periodic and bounded — the residual is only the
*rigor* of `period-max(B) ≤ 15·margin` (a finite check) and the two-far Tornheim tail
(`12ζ(3)`, arXiv:2409.19980). **Do NOT** chase a single monotone order — additive energy,
span, Schur-convexity, conductance are all refuted (HYP-2738/2780); the extremality is
irreducibly aggregate, so route through the finite-check + decorrelation split instead.

**Node 3 — the lattice-intersection correlation bound.** Bound `Σ_{n≠0}|ĉ(n)||ĝ(n)|` using
(a) `ĝ` decay `1/n` on `P`-multiples, (b) `ĉ` support on the cluster relation lattice, (c)
the covering condition controlling the overlap. The clean target is `R' ≥ c>0` uniformly
(equivalently `meas(coverSet∩G_P) ≤ (1−c)measGP`). The repo's singular-series/covolume work
(HYP-2606) and the Vitali route (HYP-2840) are partial; the new framing is to make the
**spectrum-intersection sum** explicit and route its finitely many low-height resonances to
the resonant-neighborhood patch (the rest decay). Whichever of (1)/(3) yields a clean
decorrelation lemma first will likely feed the others.

## One whimsical thread worth keeping

The denominator tower `14=2·7 → 18=2·9 → 24=8·3` mirrors Cayley–Dickson property-loss
(`lrc-cayley-dickson-tower-s387`), and `7` is the project's recurring "apex prime"
(`14=2·7`, half-width `1/14`, gap `1/7`, seven sectors, Fano-plane/Hamming-code hints in
`tournament-theory-as-coding-theory`). If the decorrelation lemma is proven *through* the
`q=7` sector structure rather than around it, the same proof should port to `n=18,24` — i.e.
the right proof of LRC(14) is probably a proof of the whole composite-`2q` family at once
(this is also what the q-uniform witness route, HYP-2846, is gesturing at). Worth testing the
Fourier-correlation bound at `n=10` (q=5) and `n=18` (q=9) to see if `c` is q-uniform.

→ HYP-2847, OPEN-Q-108, THM-523/527/530/534/563, HYP-2606/2840/2844/2846, HYP-2738/2780.
