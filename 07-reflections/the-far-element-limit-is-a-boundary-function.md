# The far-element limit is a boundary function

**Source:** mac-mini-2026-06-20-S3. Dispatch: understand the crux of the underlying structure
and complete the proof; consider "boundary functions and sets of curvilinear convergence for
continuous functions," boundary functions for bounded harmonic functions, and continuous maps of
the Riemann sphere. Canon: THM-548. Built on THM-546/547 (one-far), codex HYP-2679 (two-far
curvature), the team's Freiman/relation-lattice line (HYP-2606/2637/2678), and an 8-agent
research synthesis of the Kaczynski/Bagemihl/McMillan/Fatou boundary-function literature.

## The crux had a name already

The LRC(14) endgame is one inequality: `p0(E) ≤ cap_k`, where `p0` is the measure of phases at
which a cluster of runners covers all seven sectors. Every route into it had been
additive-combinatorial — relation lattices, covolumes, Freiman dimension, doubling constants.
The remaining open region (true-wide: a bounded core `B` plus two or more far runners) resisted
because peeling one far runner leaves a *wide* base, where the one-far bound goes vacuous.

The user's prompt named the missing frame. The far-element process **is** a boundary-value
problem, and it has the classical structure of one:

| disk / boundary theory | LRC(14) |
|---|---|
| boundary point `ξ ∈ ∂D` | a bounded core `B ⊆ {0,…,14}` |
| radial approach `z → ξ` | a single far runner `w → ∞` |
| an approach **arc** `γ → ξ` | a sequence/curve of far runners |
| the function `f(z)` | the coverage `p0(B ∪ F)` |
| **boundary function** `φ(ξ) = lim_γ f` | the plateau `Φ(B) = p0(B) + p1(B)/7` |
| Fatou: the limit exists (a.e.) | `Δ_w = p0(B∪{w}) − Φ(B)`, `|Δ_w| ≤ (6/49)V(B)/w → 0` **(PROVED)** |
| **ambiguous point** (two arcs, two limits) | a far *pair* `(u,v)` with arc-dependent curvature `I_B(u,v)` |
| the ambiguous set is small (Bagemihl: countable) | the resonant far pairs are a structured/thin set |

The one-far limit is exactly a *curvilinear limit that exists*: peel a far runner and the
coverage converges to the boundary function `Φ(B)`, with a rate I had already pinned to the
rational constant `6/49`. The two-far curvature `I_B(u,v) = p0(B∪{u,v}) − p0(B∪{u}) −
p0(B∪{v}) + p0(B)` — codex's object — is precisely the **ambiguity defect**: how much two
approach arcs disagree. It vanishes in the decorrelated limit, to an explicit value
`Φ₂(B) = (2p₂(B) − p₁(B))/49`, and its deviation is governed by the **resonance** `mu+nv`: two
far runners interfere only when they share a small additive relation — a shared asymptotic
direction, a Bagemihl ambiguous point.

## The two pictures are the same set

The boundary-function frame and the additive-combinatorics frame are not two analogies for the
crux; they are one structure seen twice. The exceptional set where the boundary function
misbehaves — the ambiguous points — is *exactly* the set of additive resonances `mu+nv ≈ 0`,
which is *exactly* the Freiman-structured set (small relations ⟺ low-dimensional GAP). So:

- **Off the exceptional set** (dissociated far runners): the coverage equals its Fatou boundary
  value `P_r(B) = Σ_t prof_t(B)·c_t(r)`, the fully-decorrelated limit as `r` runners spread.
  Verified: `P_r(B) ≤ cap_k` with margin that *grows* with `r` (0.13 → 0.48). The boundary value
  is safe.
- **On the exceptional set** (resonant far runners): scale-invariance collapses the shared
  asymptotic direction. A small relation means the runners lie in a common dilated AP/GAP, and
  the measure is dilation-invariant, so the resonant family reduces to a bounded model — a finite
  check. This is the Bagemihl "the bad set is thin and reducible," made arithmetic.

The synthesis corrected a tempting wrong turn here: the actual hardest true-wide row,
`(0,4,6,8,10,12,14,15,16)`, has **negative** curvature (`I_B(15,16) = −13/1470`) and a dilated
core `2·(0,2,3,4,5,6,7)`. It is not a positive-synergy exception; it is a dilated-core resonance
— the scale-invariant (safe) branch. Positive curvature is real but sub-critical. The danger
lives where the boundary theory says it should: on the structured exceptional set, not in the
generic interior.

## The apex prime runs the whole expansion

There is a third boundary-theoretic fact, and it is the cleanest. The coverage expands as a
**finite** Newton series in the far runners — finite because a core misses at most six sectors,
so seven or more runners are never jointly needed. The `t`-th term is the `t`-fold curvature,
and its size is set by a constant from the `t`-th antiderivative of the sector indicator:

`sup|F_j| = 3/49 = 3/7²` (one-far), `sup|G_j| = 13/1372 = 13/(2²·7³)` (two-far).

**Each order adds one power of the apex prime `7` to the denominator.** The boundary function is
"smooth" in the apex-prime sense: its successive curvilinear derivatives are suppressed
geometrically by `1/7`. The same prime that makes the singular-series kernel vanish at multiples
of seven, and that turned the one-far constant rational, now controls the entire boundary
expansion's convergence. McMillan's "+1 Baire class per order" has an arithmetic shadow here:
+1 power of 7 per far runner.

## What is — and isn't — closed

This reframing is not a proof; it is the correct shape of one, and it discharged the parts that
were genuinely about the boundary value. Verified: the boundary function exists with rate `6/49`;
the two-far ambiguity is bounded by `13/1372`, real (the QR product-reality holds), and never
cap-threatening; the Fatou boundary value `P_r(B)` is safe with growing margin; the dangerous
rows are dilated-core, the scale-invariant branch. What remains is the one place the analogy
stops being a theorem: Bagemihl gives *countable* ambiguity for value-ambiguity, but LRC
resonance is curvature-ambiguity and recurs periodically along the runner axis, so "thin" must
be earned arithmetically — the signed magnitude bound for the `d ≥ 2` relation lattice, where an
absolute bound is provably `5×` too lossy and the smallness is genuine signed cancellation. That
is the last room, and the boundary-function frame says exactly what it is: the proof that the
ambiguous set, measured in the right coordinates on far-runner space, is negligible. The
mathematics named its own missing lemma.
