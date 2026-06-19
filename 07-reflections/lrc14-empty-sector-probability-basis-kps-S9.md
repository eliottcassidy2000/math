# The empty-sector functional in the probability basis: separability restored, resonance dissolved

**Source:** kind-pasteur-2026-06-19-S9 (angle: moment-convex-order). Attacks HYP-2607 — "consec
maximizes the empty-sector moment functional `L_y(E)=E[g(N_E)]`" — the single scalar lemma that
six independent angles converged on as the LRC(14) finish line.

## The reframe: probability basis, not moment basis

`N_E(x) = #{ sectors j∈{1..6} empty of the orbit {frac(e_i x)} } ∈ {0,…,6}`, and
`meas(S7(E)) = P(N_E=0) = p_0` is the quantity the whole reduction needs `≤ cap_k`.

The canon (mac-mini-2026-06-19-S1) recorded a **non-separability obstruction**: in the binomial
moment basis `S_r = E[C(N,r)]`, consec is *not* separately extremal — over the k=8 bank one shape
beats `S_1`, fourteen beat `S_2`; only `S_4` is consec-extremal. The alternating dual `L_y = Σ y_r
S_r` trades these off, so HYP-2607 looked like an irreducibly *joint* statement on the law of `N`,
calling for a convex-order/coupling argument.

The convex-order angle's finding is that **this non-separability is a basis artifact.** Re-express
everything in the *probability* basis `p_t = meas{N_E=t}` (the natural basis for a coupling). The
k=8 dual is `g(0)=g(6)=1, g(3)=1/10, else 0`, so

> `L_y(E) = p_0 + (1/10) p_3 + p_6`.

In this basis the two unit-weight terms separate cleanly, and BOTH are consec-extremal:

- **`p_0 = meas(S7)` is consec-maximal among same-k sets** (Lemma B below);
- **`p_6 = P(N=6)` is consec-maximal among same-k sets**, equal to `1/(7(k−1))` (Lemma A).

Only the middle `(1/10)p_3` term is not separately consec-extremal — but it carries weight `1/10`
and is dominated by the `p_0,p_6` surplus. The "joint" character is confined to one small term.

## Lemma B is stronger than HYP-2607 — and it is the target itself

`p_0 = meas(S7)` IS the LRC quantity. If consec maximizes `p_0` directly among same-k clusters, we
never need `L_y` at all: `meas(S7)(E) ≤ meas(S7)(consec_k) ≤ cap_k` (the last step is THM-535's
proved finite check at k=8,9,10). The moment dual `L_y` was a *relaxation* introduced because the
direct extremality looked intractable; the probability-basis view says the direct extremality is
the honest target.

VERIFIED-EXACT this session: consec is the **unique** primitive same-k maximizer of `meas(S7)` at
all three binding rows — k=8 exhaustive over all 6434 primitive clusters with `max(E)≤15` (plus
random wide clusters to spread 200), k=9 exhaustive `max(E)≤12`, k=10 exhaustive `max(E)≤13`. Zero
beaters anywhere. The tightest non-consec primitive competitor at k=8 is `E={0,2,3,4,5,6,7,8}` with
`meas(S7)=0.2736`, a full `0.054` below consec and `0.108` below cap. The finite check is
comfortably non-tight off the AP; only the AP itself approaches cap.

**Honest scope caveat (a real failure, harmlessly placed).** Direct `meas(S7)` extremality of
consec is NOT a same-k fact for *all* k: at **k=12** the shape `{0,1,…,10,12}` *beats* consec
(`0.645 > 0.624`), consistent with HYP-2604's "AP non-extremality at large k". This is harmless
because THM-535 proves k≥9 closes via the subadditive bound `meas(S7)(consec_k) < (k−6)/7 ≤ cap_k`,
and the genuinely tight rows are exactly **k=8,9,10** — where Lemma B holds. So Lemma B is a
*binding-row* statement: consec is the same-k `meas(S7)`-maximizer precisely where it must be. (The
moment functional `L_y` has the same harmless non-extremality pattern at k≥11; the dual degree
drops to 2 there and consec is no longer the `L_y`-max either, per THM-534.)

This does not *prove* the same-k extremality — that is the irreducible analytic core, the "aggregate"
phenomenon THM-536 and HYP-2608 identified. But it relocates the open statement to its cleanest
form: a direct same-k extremality of the target measure, with the relaxation stripped away.

## Resonance dissolved by primitivity

HYP-2608's live route (a) — a wide-spread bound — carried one specific fear: the "resonant" wide
configs `w ≡ 0 mod 7` (apex-prime-7, THM-503), where the wide-spread collapse might fail. The
probability-basis stress test shows **these resonant configs are exactly the non-primitive ones.**
`{0,7,14,…,49}` has `gcd=7`; by scale-invariance (THM-531) its `meas(S7)` *equals* that of its
primitive representative `{0,1,…,7}` = consec — so it sits at consec, not above cap. Every
`{0..6}+7m` tail for `m≥2` is *primitive but strictly below consec*; only the `m=1` "tail" coincides
with consec because it IS consec.

Since the LRC reduction operates on primitive clusters (we divide out the gcd), the resonance is a
non-issue. The wide-spread bound only ever has to beat *primitive* wide configs, and those collapse
cleanly: primitive `meas(S7)` is `≈0.03–0.08` at spread 40–200, an order of magnitude under cap.
The thing that looked like an obstruction was a scale-copy of the extremiser in disguise.

## Lemma A: a genuinely rigorous building block

`P(N=6)(E) = meas(G_E)`, `G_E = {x : frac(ex)∈[0,1/7) ∀e∈E}` — the set where the whole cluster
orbit is confined to the first 1/7-arc. Two facts:

1. **Per-component length bound (PROVED, elementary).** On any connected component of `G_E` the
   winding integers `⌊ex⌋` are constant, so the constraints read `x ∈ ⋂_e [n_e/e, n_e/e + 1/(7e))`
   — an interval of length `≤ min_e 1/(7e) = 1/(7 max(E))`. (0 violations / 3176 primitive E.)
2. **Component-count bound (VERIFIED, the remaining gap).** `#components(G_E) ≤ max(E)/(k−1)` for
   primitive E. (0 violations.)

Together they give `meas(G_E) ≤ 1/(7(k−1)) = P(N=6)(consec_k)`, extending THM-535's
`Φ(c,k)=c/(k−1)` from the consecutive cluster to the same-k *maximizer*. The first factor is
airtight; the open piece is a clean lattice-point count, not an analytic estimate.

## What this points to

The honest distance to LRC(14) is unchanged in length but clearer in shape. The single open lemma
is best stated as **"consec maximizes `meas(S7)` among primitive same-k clusters"** — the target
measure directly, the moment relaxation discarded, the resonance worry dissolved into a primitivity
remark, and one of its two extreme-event terms (`p_6`) reduced to a proved per-component bound times
a lattice count. The probability basis is the right coordinate system for the coupling that finishes
this; the moment basis hid the structure behind a false non-separability.

Files: `04-computation/lrc14_wsb_{moment-convex-order,convexorder_stage2,p6formula_stage3,
p0p6_extremal_stage4,pN6_proof_stage5,primitive_stage6,resonance_widespread_stage7,
lemmaA_components_stage8,FINAL_summary}_kps-S9-wf.py` (+ `.out` in `05-knowledge/results/`).
