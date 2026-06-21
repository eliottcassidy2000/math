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
survives. What is new is that the wall is now a clean statement in invariant-theory language —
*the `Z_7*`-fixed configuration extremizes a `Z_7*`-invariant sum of cycle-graph survival widths* —
and that it is literally the same `Z_7*` story as the tournament Paley extremality. The remaining
content on both sides is the single Gauss-sum-twisted aggregate gap.
