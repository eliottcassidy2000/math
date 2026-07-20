---
id: THM-1390
title: "LEVERAGE AUDIT OF FIVE UNFAMILIAR THREADS — what the repo's novel machinery can and cannot reach, with four decisive computations. (0) THE BIG ONE, filed as a court case: the path-homology corpus uses a NON-GLMY regularity convention (distinct vertices, where GLMY requires only consecutive distinctness), so every degree->=3 result is computed in the wrong complex — beta_4(T_7) = 0 under GLMY, not 6, and chi(T_p) = p is not merely wrong but UNDEFINED since dim Omega_p does not terminate. Independently verified from the definitions with a calibrated control. THM-103 (beta_1 <= 1) and THM-108 (beta_2 = 0) are CONVENTION-SAFE and survive, by a subspace argument. (I) ERDOS-HAJNAL IS UNREACHABLE from the repo's central invariant: hp (Hamiltonian-path count) does NOT control tr (largest transitive subtournament) — explicit n=6 violations in BOTH directions. The invariant that does control it is the feedback arc set, via the elementary tight bound tr >= n - fas, but that is vacuous exactly in the EH regime since fas ~ n^2/4. (II) THE hp-SPECTRUM STABILISATION LEMMA, which makes the {7,21} gap question well-posed: a dominant vertex preserves hp, so the spectra are NESTED and 'the spectrum' is a single increasing union rather than an n-indexed family. Exhaustive n<=6 plus fas<=4 search at n<=9 reproduces 7 and 21 as the only missing odd values <= 60, with 35 and 39 filling in at n=7. (III) THM-488's rigidity certificate stopped at |S| <= 6 because that is where its INSTRUMENT stopped, not the mathematics: |z_min(S_j)| leaves the truncation-safe zone at exactly j = 7. Replacing the truncated series by the Euler PRODUCT removes the artifact and extends the certificate to j = 14, where |z_min| SATURATES at ~0.9723 rather than approaching 1"
status: >
  (0) VERIFIED-EXACT and filed as 02-court/active/CASE-path-homology-regularity-convention.md.
  Written from the GLMY definitions independently of repo code, exact linear algebra
  over a prime field, with the directed 3-cycle as a calibrated control returning the
  textbook beta_1 = 1.  Reproduces an agent report's numbers exactly.  This CONTRADICTS
  canon (THM-129, THM-130, degree->=3 tables of THM-096/099/154/226), hence a court case
  rather than an override.
  (I) PROVED (the lemma tr >= n - fas, three lines) + VERIFIED-EXHAUSTIVE (n <= 6, all
  2^C(n,2) labelled tournaments; 0 violations; tight 100/74/37% at n = 4/5/6).  The
  non-control of hp is established by EXPLICIT WITNESSES, not by failure to find a bound.
  (II) The stabilisation lemma is PROVED (elementary).  The spectrum claim is
  VERIFIED-EXHAUSTIVE at n <= 6 and a SEARCH at n = 7..9 — searches establish attainment,
  never non-attainment, and the detection floor is stated.  {7,21} non-attainment itself
  is canon's (THM-343, THM-115), not re-proved here.
  (III) VERIFIED, with a stated numerical caveat: the bare Euler product's winding is
  unreliable above r = 0.95 (underflow to ~1e-72), though F_S is perturbation-dominated
  there.  The saturation is the substantive finding and it is robust across j = 12,13,14.
  NOTHING here proves a famous open problem.  Two lines are CLOSED as unreachable, which
  is the point.
source: kind-pasteur-2026-07-20-S128c106 (owner: explore the repo's novel math in unfamiliar areas and see what it can leverage toward open problems)
depends_on:
  - THM-103     # beta_1 <= 1 -- shown convention-safe here
  - THM-108     # beta_2 = 0  -- shown convention-safe here
  - THM-488     # the rigidity theorem whose certificate is extended
  - THM-343     # H != 7
  - THM-115     # H != 21
related: [THM-129, THM-130, THM-096, THM-099, THM-154, THM-226, THM-455, MISTAKE-021, MISTAKE-196]
script: 04-computation/glmy_convention_audit_kps_S128c106.py, hp_vs_trans_stability_kps_S128c106.py, fas_controls_transitivity_kps_S128c106.py, hp_spectrum_gaps_kps_S128c106.py, rigidity_Sj_stress_kps_S128c106.py, rigidity_product_form_kps_S128c106.py (+ .out)
---

# THM-1390 — leverage audit: what reaches an open problem, and what does not

Five threads I had never worked were surveyed for leverage toward *named* open
problems. Most of the value is negative — two plausible bridges are closed — and
the largest single item is a convention error in canon.

## 0. The path-homology convention (filed as a court case)

Every path-homology implementation here takes an elementary `p`-path to be a
sequence of **distinct** vertices. GLMY require only **regularity**, `i_k ≠
i_{k+1}`; non-consecutive repeats are allowed. In a tournament `a→b→a` is
impossible, but `a→b→c→a` is allowed whenever `abc` is a 3-cycle — so the
conventions must first diverge at `p = 3`.

Computed from the definitions, with the directed 3-cycle as a control (`β₁ = 1`,
textbook, both conventions):

| Paley `T₇` | `dim Ω₀..Ω₅` | `β₀..β₄` |
|---|---|---|
| **GLMY** | 7, 21, 42, **70, 105, 147** | 1,0,0,0,**0** |
| **repo** | 7, 21, 42, **63, 63, 42** | 1,0,0,0,**6** |

`Ω₀,Ω₁,Ω₂` agree, so **THM-103 (`β₁ ≤ 1`) and THM-108 (`β₂ = 0`) survive** — and
`β₂ = 0` survives *a fortiori*: the repo's `Ω₃` is a subspace of the true one, so
`im ∂₃` only grows while `ker ∂₂` is unchanged. What falls is everything in degree
≥ 3, including `β₄(T₇) = 6` and `χ(T_p) = p` — the latter not merely false but
**undefined**, since `dim Ω_p` does not terminate under GLMY.

Full charge, relief sought, and the calibration protocol for anyone reversing this:
`02-court/active/CASE-path-homology-regularity-convention.md`.

## I. Erdős–Hajnal is not reachable from `hp`

The repo's central invariant is `hp(T)`, the Hamiltonian-path count (Rédei: odd;
`hp = 1` iff transitive). The Erdős–Hajnal quantity is `tr(T)`, the largest
transitive subtournament. Both measure distance from transitivity from opposite
ends, and the extremes agree, so a stability bound `tr ≥ φ(n, hp)` looks plausible.

**It is false.** Exhaustively at `n = 6`, `min tr` and `max tr` are both
non-monotone in `hp`, with explicit witnesses:

> `min`: `(hp,tr) = (9,4) → (11,5)` and `(17,4) → (19,5)`  ·  `max`: `(13,4) → (15,5)`

i.e. **more Hamiltonian paths with more transitivity**. The controlling invariant
is instead the feedback arc set:

> **Lemma.** `tr(T) ≥ n − fas(T)`.
> *Proof.* Reversing a minimum feedback arc set `A` makes `T` transitive. Delete one
> endpoint of each arc of `A` — at most `fas(T)` vertices. The survivors meet no arc
> of `A`, so they induce the same tournament as the reversed one: transitive. ∎

Exhaustive `n ≤ 6`: 0 violations; tight for 100% / 74% / 37% of tournaments at
`n = 4/5/6`. But `fas ~ n²/4` for a random tournament, so `n − fas < 0` exactly in
the regime EH cares about. **Both invariants are dead ends for EH**, one by
counterexample and one by vacuity. Recording this so the bridge is not rebuilt:
`hp` is a *global count*, `tr` a *local extremal* quantity, and the flip experiment
shows `hp` moving by an order of magnitude while `tr` does not move at all.

## II. The `hp`-spectrum stabilisation lemma

The claim "the attained set is `{odd} \ {7,21}`" is not well-posed as stated — the
attained set depends on `n`. It becomes well-posed:

> **Lemma.** If `v` is a dominant vertex of `T'` (out-degree `n`, in-degree 0) then
> `hp(T') = hp(T' − v)`.
> *Proof.* Every Hamiltonian path contains `v`; `v` has no in-arc so it is the
> initial vertex; the remainder is a Hamiltonian path of `T'−v`, and conversely each
> extends uniquely. ∎

So the spectra are **nested** and "the spectrum" is their increasing union — a
single set, and "is 7 attained?" is one question, not one per `n`. Verified nested
at `n = 3..6`.

Independent reproduction: exhaustive `n ≤ 6` plus every `fas ≤ 4` perturbation at
`n = 7,8,9` plus sampling gives, among odd values `≤ 60`, exactly `{7,21}` missing —
with `35` and `39`, absent at `n = 6`, filling in at `n = 7`. Non-attainment of
`{7,21}` is canon's (THM-343, THM-115); the searches here can only confirm
attainment, and the detection floor is `fas ≤ 4`, `n ≤ 9`.

## III. THM-488's certificate stopped where its instrument stopped

THM-488 certifies the hard half of its rigidity theorem on `|S| ≤ 6`, and names its
own near-boundary sets — all initial segments of `S_j = {1,3,…,2j−1}`. Measuring
`|z_min(S_j)|` with its truncated-series winding:

> `j` = 1…6 → `0.7815, 0.8801, 0.9101, 0.9262, 0.9411, 0.9497`, then **nothing
> detectable from `j = 7`**.

`0.9497` is just inside the truncation-safe ceiling `r ≤ 0.95`. **The certification
range coincides with the reach of the instrument, not with a change in the
mathematics** — the same failure mode as MISTAKE-196, one level up.

The artifact is avoidable. The perturbation is *finite*, so writing
`F_S = Π(1−xⁿ) + 2Σ_{k∈S}(−1)^{k+1}(x^{g_k}+x^{ḡ_k})` and evaluating the first term
as a **partial product** — zero-free in the disk, no spurious zeros — removes the
0.95 ceiling. Extending to `j = 14`:

> `|z_min|` = 0.9552, 0.9611, 0.9650, 0.9678, 0.9709, 0.9721, 0.9723, 0.9723
> for `j = 7…14` — **saturating at ≈ 0.9723, not approaching 1.**

So Theorem A's hard half survives well past `|S| ≤ 6` on its own worst family, and
the "pushes the zero toward the boundary" diagnostic is a small-`|S|` effect that
stalls. Structural reason: near `|z| ≈ 0.97` the Euler product is `~10⁻²⁴`, so `F_S`
is dominated by the finite polynomial `P_{S_j}`, whose least zero modulus converges.

**Numerical caveat, stated rather than buried.** The control — the *bare* Euler
product, which must have winding 0 everywhere — returns 0 at `r ≤ 0.95` but 12 at
`r = 0.99`, because the product underflows to `~10⁻⁷²` and the winding is then pure
noise. `F_S` is perturbation-dominated there and well conditioned, so its windings
are trustworthy where the bare control's are not; but the control as designed does
not certify that, and a sharper run should use interval or extended precision.

## Named next

- **Redo the homology corpus in true GLMY** before any external write-up, then check
  overlap with arXiv:2602.04140 (Paley `T_p` are circulant digraphs) — canon's own
  novelty explainer flags this and it has never been done.
- **Close THM-108's embedded bridge**, "key equivalence verified n ≤ 12", which is an
  unproved step inside a theorem marked PROVED.
- The `S_j` saturation at `≈ 0.9723` deserves a closed form: it should be the least
  zero modulus of `lim_j P_{S_j}`, a lacunary series in the odd pentagonal lags.
