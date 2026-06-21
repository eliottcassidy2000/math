# The multiplicative group Z_7* is the shared spine of both extremalities

**Session:** mac-mini-2026-06-21-S19 (high-throughput trawl, concurrent with kps/codex). Builds on
the unification (HYP-2747), kps's dilation symmetry (HYP-2754) and Paley=QR result (HYP-2755), and my
resonance localization (HYP-2753) + LAYERS 1-2 (HYP-2756).

## The two extremizers are the same kind of object

The project has two open extremalities, long tracked separately and recently shown to share one
Lovász-θ′ ceiling (HYP-2747):

- **LRC(14):** the consecutive block (the AP) maximizes the Z/7 sector cover `measS7`.
- **Tournaments:** the Paley/QR tournament maximizes `H`.

This session pinned *why they are the same object*: **both are the maximally MULTIPLICATIVELY-symmetric
object on Z_7.** The shared symmetry is the multiplicative group `Z_7* ≅ Z_6`, which acts on **both**
sides — as **dilation** `E ↦ cE` on offset sets (LRC), and as the **circulant/QR structure** on
connection sets (tournaments).

## LAYER 2 is forced by multiplicative symmetry (the new rigorous step)

The consec-max wall localizes (HYP-2753) to three layers; LAYER 2 said "consec doubles the *identity*
residue 0," verified but unexplained. It is now a symmetry consequence, **proved**:

> A full-residue `k=8` shape has residue multiset `{0,1,…,6}` plus one repeat `r*`. That multiset is
> `Z_7*`-(dilation)-invariant **iff `r* = 0`** — because dilation fixes residue 0 and freely permutes
> `{1,…,6}`, so a repeat anywhere but the fixed point is moved by some `c`. (Verified: only `r*=0`.)

So "double the identity residue" = "be dilation-symmetric," and **consec is the minimal-magnitude
dilation-symmetric full-residue shape.** The dilation also acts on the resonance decomposition itself:
`W_a(cE) = W_{ca}(E)` (HYP-2757), so `Z_7*` permutes the six resonance phases and `measS7 = Σ_a W_a`
is the orbit-sum — manifestly dilation-invariant (= THM-531, now seen per-phase).

On the tournament side the mirror is exact (kps HYP-2755): Paley = the QR cyclic *linear* code is the
unique `Z_7*`-structured (its connection set is the index-2 subgroup of `Z_7*`) circulant H-maximizer.

## Additive vs multiplicative — the honest asymmetry

The duality is not a symmetry: the two sides sit on the two group structures of `Z_7`.

- **LRC consec = AP** is *additively* closed (a full coset of the additive group), but its extremal
  symmetry is the *multiplicative* dilation `Z_7*`. It is **non-linear** as a code — which is exactly
  why the CJJ/Delsarte hierarchy cannot certify it until the dilation (affine group) is factored out
  (kps HYP-2754); after factoring, consec's orbit is pinned.
- **Tournament Paley = QR** is *multiplicatively* closed (a subgroup of `Z_7*`), and **is** a linear
  code — so its extremality is the CJJ-linear, potentially certifiable case (kps HYP-2755).

The bridge between additive and multiplicative characters on `Z_7` is the **Gauss sum** (|g| = √7),
and that is the precise sense in which "whatever proves one extremality is owed to the other": a single
argument must intertwine the additive picture (Walsh/Krawtchouk, the AP) with the multiplicative one
(QR, Gauss sums). The shared object is `Z_7*`; the asymmetry (one linear, one not) is the Gauss-sum
twist between the two character groups.

## What this buys, honestly

It does **not** close LAYER 3 — the genuinely-aggregate step that the *most* multiplicatively-symmetric
point maximizes the dilation-invariant `Σ_a W_a`. Pure symmetrization cannot, because `measS7` is
*constant* on dilation orbits, not Schur-concave under them: the symmetry **narrows** the maximizer to
the dilation-symmetric family (LAYERS 1-2, now both rigorous) but the residual cross-phase trade-off
survives.

## CORRECTION (S20, harvested from the creative-lead trawl) — the unifier is the ADDITIVE AP, not Paley

The thesis above leaned on "Paley = the tournament extremizer," and that is **wrong for large p** — a
small-`p` coincidence I should not have generalized. The canon already records it (HYP-479 / THM-135,
the **Paley Crossover**; HYP-455): Paley/QR maximizes circulant `H` only at `p = 3, 7, 11`; at `p = 19`
the **cyclic Interval** (the *additive* AP / consecutive block) strictly beats it
(`H_Interval = 1.184e12 > H_Paley = 1.173e12`), and the Interval is the unique H-maximizer thereafter.
The trawl's ANGLE 6 re-derived this independently (Schur-concavity of `H` in `|λ|²` holds at p=7,11 and
**fails at p=19**, exactly where Interval majorizes Paley but has higher `H`).

So the honest unification is **stronger and simpler than the multiplicative one**: the extremizer on
*both* sides is the **additive arithmetic progression** — `consec = {0,…,k-1}` for the LRC cover and the
cyclic **Interval** `{⌈p/2⌉,…,p-1}` for tournament `H`. "Whatever proves one proves the other"
(HYP-2747) becomes literal: both are *"the AP Schur-maximizes."* The `Z_7*` multiplicative/dilation
symmetry is a genuine **feature** of `consec` (it forces LAYER 2 — double the identity residue —
verified and Lean-formalized) but it is **not the source of extremality**: the pure-multiplicative
object (Paley/QR) is *not* extremal once `p ≥ 19`, while the pure-additive object (the AP) is. The
"Gauss-sum twist" framing is therefore demoted: it explains why the two *coincide* at small `p`, not why
either is extremal. The driver is additive (AP / minimal-conductance), and the proof tools that match
the data are the conductance vector `c_r = Σ_{e≡r} 1/|e|` and the windows/binding-speed harmonic sum
(S20 trawl ANGLEs 1–2), not multiplicative symmetry.
