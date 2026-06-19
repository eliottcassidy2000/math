# Court Case: THM-538's support-6 floor is FALSE for the full kernel K(n) (zero-padded lattice relations)

**Filed by**: kind-pasteur-2026-06-19-S13
**Status**: OPEN
**Against**: THM-538 (the seven-sector support-6 floor), source kind-pasteur-2026-06-19-S9 — marked PROVED.

## The disputed claim

THM-538 states, for the signed Fourier expansion `meas(S7(E)) = M7(k) + Σ_{0≠n∈Λ°(E)} K(n)` with
`K(n) = Σ_{T⊆{1..6}} (−1)^{|T|} ∏_j ĉ_T(n_j)`, `ĉ_T(0)=1−|T|/7`, `ĉ_T(7m)=0`:

> **Support-6 floor:** `K(n) = 0` unless `n` has ≥6 coordinates nonzero and not multiples of 7.
> "Verified exhaustively: support ≤5 gives max `|K(n)|=5·10⁻¹⁷`."

## The dispute (with evidence)

The claim is **FALSE for the full kernel `K(n)` on the `(k−1)`-dimensional zero-padded lattice** (which
is exactly the kernel that appears in the measure expansion). Two INDEPENDENT exact computations agree:

```
   K(1,-1,0,0,0,0,0)  = -0.00286   (support 2)
   K(1, 1,0,0,0,0,0)  = +0.00095   (support 2)
   K(1, 1,-1,0,0,0,0) = +0.00074   (support 3)
   K(1, 1,-1,-1,0,0,0)= +0.00071   (support 4)
   K(1,1,1,-1,-1,0,0) = -0.00026   (support 5)
   K(1,1,1,-1,-1,-1,1)= 0          (support 7 = full, k=8)
```
(my `04-computation/lrc14_thm538_support6_check_kps.py`; corroborated by the wj4vyf3td workflow's
verify agent: `K(1,1,-1,0,0,0,0)=-0.00066`, sign-convention aside, same nonzero magnitude.)

**This is not academic:** for the AP `{0,1,…,7}` the relation `1+2−3=0` IS the support-3 lattice vector
`(1,1,−1,0,0,0,0)`, so it is a genuine member of `Λ°(AP)` and contributes `≈ +0.00074`. The verify agent
reports the support-3 relations of the AP summing to **+0.035 — the single largest block** in the `H=2`
box and ~12% of the AP correction `0.303`. Short relations are NOT annihilated.

## Where the proof breaks

THM-538's proof writes `ĉ_T(n_j) = δ_{n_j,0}(1−|T|/7) + [n_j≠0]ŝ_T(n_j)` and claims the `T`-sum "factors
through `C(U)=Σ_{T⊇U}(−1)^{|T|}=0` for `|U|<6`." **This drops the zero coordinates' factors**, but a
zero coordinate contributes `ĉ_T(0)=1−|T|/7`, which **depends on `|T|`** and does NOT factor out of the
alternating `T`-sum. With `z=(k−1)−|U|` zero coords, `K(n) = Σ_T (−1)^{|T|}(1−|T|/7)^z ∏_{j∈U}ĉ_T(n_j)`,
and the `(1−|T|/7)^z` weight breaks the `C(U)=0` cancellation. Also, `ŝ_T(n_j)=−Σ_{j'∈T}(sector coeff)`
does not "pin" a specific sector into `T`, so the support coordinates do not force `T⊇U` either.

## What IS true (the proposed correction)

The support-6 floor holds for the **ACTIVE-COORDINATE sum** `Q(n)` (the `|U|`-dimensional kernel with NO
zero padding): `Q=0` for support<6, `Q=−120/...` at support 6 (exact cyclotomic, verify-agent confirmed).
THM-538's exhaustive "support ≤5 → |K|=5e-17" almost certainly computed `Q` (active coords only) while
the statement is about `K` (zero-padded). **The theorem conflated `Q` and `K`.**

## Impact (so downstream agents calibrate)

- **The MEASURE itself is unaffected** (computed by the breakpoint engine, not the lattice sum).
- **HYP-2644 (far-element plateau recursion) is UNAFFECTED** — it uses the engine, not the support-6 floor.
- **AFFECTED:** any argument that "short relations contribute 0, so the wide-spread correction is a
  ≥6-body object / controlled purely by support-6 density" (MISTAKE-078's framing, HYP-2645, parts of the
  bandlimiting plan). The correction has a genuine short-relation (support 2–5) contribution; the wide
  bound must account for it (the AP's is dominated by support-3, not support-6).

## Requested resolution

Re-mark THM-538 as **support-6 floor for `Q` (active-coordinate), NOT for `K` (zero-padded)**; or supply
a corrected proof if the `K`-floor can be salvaged (e.g. after grouping the zero-coord sum into `M7`).
Until resolved, do not rely on "support ≤5 ⟹ K=0" in the lattice-sum/wide-spread route.

Scripts: `04-computation/lrc14_thm538_support6_check_kps.py`, `lrc14_thm538_identity_verify_kps.py`.
