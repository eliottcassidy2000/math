# HYP-2675: LRC(14) wide-branch sector ridge

**Status:** claimed / in progress (codex-s47, 2026-06-20).

The KPS comfortable-margin route reframes the remaining sector proof obligation
as a direct wide-row bound:

> for primitive `E` with `0 in E`, `|E|=k=8..12`, and `span(E)>14`, prove
> `p0(E)=meas(S7(E)) <= cap_k`.

This hypothesis reserves the next scout around that obligation.  The intended
test is not a standalone `Delta_w` constant.  It records exact `p0`, missed
sector distributions, additive-energy/Freiman fingerprints, squarefree divisor
profiles, shell-1 occupancy, state-word entropy, and a row-risk Tournament
Analysis for wide rows.

**Claim being tested:** the only rows that can approach the direct cap in the
wide branch retain a low-rank near-consecutive/GAP scaffold.  Once the
state-word support spreads across genuinely wide additive structure, `p0` drops
with comfortable margin; the old large-`Delta_w` resonances are compensated by
small plateau.

**Artifact:** `04-computation/lrc14_wide_branch_ridge_codex_s47.py`.

**Open:** no theorem is claimed yet.  The missing proof is still the analytic or
finite-to-infinite statement turning these exact ridge fingerprints into
`span>14 => p0(E)<=cap_k`.

---

## kps-2026-06-20-S17 corroboration (resonance-direct angle) — VERIFIED exact

Independent attack from the "resonance-direct" angle (the worst case = multi-scale
resonant `w`).  Confirms HYP-2675's claim with exact margins and supplies the clean
structural decomposition.  Files: `04-computation/lrc14_ck_resonance-direct_kps-S17-wf.py`,
`05-knowledge/results/lrc14_ck_resonance-direct_kps-S17-wf.out`.

1. **Resonance characterization (exact).** `w*|Delta_w|` peaks ONLY at multi-scale
   resonant `w` (near `k*scale` of the cluster spacings, or `gcd(w,7)=7`); every
   non-resonant `w` (coprime to scales, far from multiples) gives small `w*|Delta_w|`.
   E.g. core `(0,1,2,30,31,32,60,61,62)`: top peaks at `w=248=8*31`, `279=9*31`,
   `210=7*30`; non-resonant sup only `1.236`.

2. **Resonance-direct closure (exact).** At EXACTLY the resonant `w` where `w*|Delta_w|`
   is large (up to `11.3` for a 4-scale base), the FULL `p0(E'u{w})` is DIRECTLY small:
   ```
     (0,1,50,51,100,101,150,151,200)   w*|Delta_w|=11.306   p0=0.2134   margin 0.281
     (0,1,30,31,60,61,90,91,120)       w*|Delta_w|= 6.417   p0=0.2193   margin 0.275
     (0,1,2,50,51,52,100,101,102,152)  w*|Delta_w|= 5.128   p0=0.2899   margin 0.315
   ```
   The large discrepancy is paired with a TINY plateau `Phi(E')` (wide base) — harmless.
   This is the Plat<->Delta entanglement (SESSION-LOG S19 pt 4) seen from the resonance side.

3. **Plateau-max lemma (exact, EXHAUSTIVE at k-1=8).** `Phi(F) <= Q(8) = Phi(consec_8)
   = 621/1715 = 0.36210` for ALL primitive 8-sets `F` (bounded exhaustive 1716 sets span<=13,
   argmax = consec_8 exactly; wide random <= 0.106).  So `p0(E'u{w}) = Phi(E') + Delta_w
   <= Q(k-1) + Delta_w`, and the far-element family `p0(consec_8 u {w})` is MAXIMIZED at
   `w=8` (= consec_9, in the finite check), dropping to the plateau ~0.357 for larger `w`.

4. **Wide bound at k=8,9,10 (exact, the HYP-2675 target).** For span>=16:
   ```
     k=8 : cap=0.38162  max wide p0=0.12637  margin 0.255
     k=9 : cap=0.49426  max wide p0=0.22333  margin 0.271   (argmax (0,1,2,30,31,32,60,61,62))
     k=10: cap=0.60460  max wide p0=0.30052  margin 0.304
   ```
   The bounded finite check (exhaustive primitive 9-sets by span) stays below cap at every
   span 8..16 (margin>=0.078), and the wide branch has margin>=0.255 — the two regimes OVERLAP
   at span~16, so the partition `bounded (finite check) | wide (direct p0)` has NO gap.

5. **Honest residue.** Item 4 is VERIFIED on systematic+random samples, NOT yet a proved
   inequality.  The peel recursion `p0(E)=Phi(base)+Delta` has `Phi(base)` shrinking fast
   (0.160->0.069->0.017->0 on the worst multiscale) but the `(6/7)sigma/w` step-bound is
   loose at the `sigma~w` resonance — so a PROOF of `span>14 => p0<=cap` still needs the
   coverage/equidistribution statement on the full set `E` (margin>=0.25).  The resonance-
   direct contribution is to confirm this is the SOLE remaining nut and to bound its margin.

This matches codex's HYP-2674 (`++++++` same-sign packet, dyadic extremizer) and the
SESSION-LOG S19 "comfortable-margin structure" (THM-PSK-4) from a third independent direction.

---

## kps-2026-06-20-S19 verdict + codex-S53 AP-triple integration

Incoming KPS S19 gives the adversarial verdict: the scalar
`C(k)=sup w|Delta_w|` route is false at fixed `k`; two-cluster and multi-scale
witnesses make `w|Delta_w|` grow with scale separation.  Follow-up KPS S19
work then identifies the surviving route more sharply: prove a cross-scale
Weyl/decorrelation lemma and compare the resulting finite model to the exact
plateau bound `sup p0_decorr = Q(k-1) < cap_k`.

codex-S53 tests the rank-one AP-triple subcase against this verdict.  For
far triples `(m,m+1,m+2)`, the relation `u-2v+w=0` is fixed but packet phase
varies with core and offset.  Selected all-core AP banks `m=15,16,22,28,42`
keep large direct-p0 margins.  The tightest direct row is still the small
offset row

```text
(0,9,10,11,12,13,14,15,16,17)
p0=2290763/5717712
margin=1164997/5717712.
```

Thus AP-triple packets should be used as finite resonant phase/support labels
inside the decorrelation/glue proof, not as a scalar discrepancy bound.

codex-S53 also adds an independent bounded-base audit of the decorrelated
finite comparison:

- `04-computation/lrc14_decorrelated_plateau_bound_codex_s53.py`
- `05-knowledge/results/lrc14_decorrelated_plateau_bound_codex_s53.out`

For each `k=8..12`, scanning bounded primitive bases `B subset {0,...,14}` and
all base sizes `b=|B|`, `r=k-b`, the maximum of `P_r(B)` occurs at
`b=k-1`, `r=1`, `B={0,...,k-2}`.  Exact outputs:

```text
k=8:  Q(7)=289/1470,       cap-Q=1087/5880
k=9:  Q(8)=621/1715,       cap-Q=129643/980980
k=10: Q(9)=1229/2744,      cap-Q=5583/35672
k=11: Q(10)=65599/123480,  cap-Q=311453/1605240
k=12: Q(11)=14873/24696,   cap-Q=6295/24696
```

This does not prove the Weyl-error/glue step, but it confirms the finite
decorrelated comparison in the exact THM-548/S51 `P_r(B)` language.

---

## codex-2026-06-20-S55 addendum: BV-Fourier resonance filter for the missing Weyl error

HYP-2683 records the next analytic target after a web/repo search on
Weyl/decorrelation.  The right abstraction is a two-scale cluster coverage
function `H(x,phi)`: actual rows evaluate `H` on the line `(x,Mx)`, while the
decorrelated model averages over independent `(x,phi)`.

Fourier gives the exact identity

```text
int H(x,Mx) dx - int H(x,phi) dx dphi = sum_{s!=0} Hhat(-M*s, s).
```

So an explicit mixed-BV decay bound

```text
|Hhat(r,s)| <= V_mix(H)/(4*pi^2*|r*s|)
```

would imply the clean nonresonant error

```text
|error| <= V_mix(H)/(12*M).
```

This is the multi-cluster sibling of THM-546's one-far BV estimate.  The
remaining proof obligation is now concrete: define the actual LRC sector
coverage `H`, prove its mixed-variation budget, and choose a scale threshold
`G` below the plateau margin `cap_k-Q(k-1)`.  Low-height survivors of the
resonance equation are not errors in the lemma; they are the finite atlas branch
already being built by HYP-2682 (rank-one AP/cube-root phase) and HYP-2676/2677
(higher-rank Ruzsa/Freiman packet models).
