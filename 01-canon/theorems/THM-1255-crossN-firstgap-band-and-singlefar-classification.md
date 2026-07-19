---
id: THM-1255
title: The cross-N first-gap band — the single-far stratum of the first gap (1/(N+1), 2/(2N+1)) is COMPLETELY CLASSIFIED for N = 6..13 (exactly two members, the canonical mediant attainers at N = 7 and N = 13), the N = 8..11 first gaps are EMPTY across a validated multi-species census (the July-6 O-depth-monotone middle, finally executed), and the populated N in [6,13] are exactly {6} ∪ {N ≡ 1 mod 6}. Corollary at N=13 — no single-replacement family enters opus's interval (1/14, 3/41).
status: >
  PROVED — the single-far classification at N ∈ {6,...,13} (absorption lemma + exact finite
  check; base floors cited from SETTLED LRC(≤13)); the N=12 case is an independent
  replication of THM-1004 (klein), discovered post hoc — see the convergence note.
  VERIFIED (evidence-grade, gates passed) — first-gap emptiness at N = 8, 9, 10, 11 over
  ~0.6M–1.7M validated-species families per N. CONJECTURE — the cross-N gate law and the
  N=31 prediction.
source: death-star-2026-07-19-S59 (HYP-7885)
depends_on:
  - THM-1002  # pair-sum denominator lemma (exact evaluator; general-M form)
  - THM-1004  # Hamming-1 rigidity at N=12 (klein) — the N=12 row; my run replicates it
  - THM-1005  # Hamming-2 rigidity at N=12 (klein) — replicated by the nested two-far run
  - THM-633   # the i=12 ladder m/(12m+1) at N=12 (mac-mini, Lean)
external: LRC(≤13) SETTLED (owner directive 2026-07-02) — supplies M(B) ≥ 1/N > 2/(2N+1) for every base B with ≤ 12 speeds.
related:
  - HYP-4506  # opus-S118 cross-N map (N=8..11 unverified middle — closed here)
  - HYP-4516  # opus-S119 mod-30 binder gate (canonical family) — its solutions proved UNIQUE in-stratum here
  - HYP-4542  # mac-mini-S25 depth 2->1->0 (the conjecture this theorem's data refines)
  - HYP-7840  # opus-S398 bound-D = bound-speeds (the N=13 corollary feeds this)
  - THM-1240  # opus (1/14,3/41) not settled — single-far stratum now closed
scripts:
  - 04-computation/lrc_firstgap_crossN_census_deathstar_S59.py -> 05-knowledge/results/lrc_firstgap_crossN_census_deathstar_S59.out
  - 04-computation/lrc_singlefar_absorption_atlas_deathstar_S59.py -> 05-knowledge/results/lrc_singlefar_absorption_atlas_deathstar_S59.out (+ lrc_singlefar_dip_atlas_N12_deathstar_S59.out)
  - 04-computation/lrc_twofar_nested_absorption_N12_deathstar_S59.py -> 05-knowledge/results/lrc_twofar_nested_absorption_N12_deathstar_S59.out
---

# THM-1255 — the cross-N first-gap band

Setting: N distinct positive integer speeds; `M(A) = max_t min_{v∈A} ‖vt‖`; the **first
gap** at N speeds is the open window `W_N = (1/(N+1), 2/(2N+1))` between the conjectured
floor and its Farey neighbour. Known before this session (opus-S118, 2026-07-06): `W_N`
nonempty at N = 6 (5/33), N = 7 (3/23), N = 13 (3/41); empty at N = 12 (fleet censuses).
**N = 8, 9, 10, 11 were explicitly unverified** (mac-mini-S25: "depth 2→1→0, n=9..12
unverified"; the bordered-AP enumeration was named as the next step and never run).

## 1. The far-element absorption lemma (elementary)

> **Lemma.** Let `B` be a finite set of positive integers, `θ ∈ (0, 1/2)`, and let
> `ℓ(B,θ)` be the length of the longest closed interval contained in
> `S(B,θ) = {t ∈ [0,1] : ‖ut‖ ≥ θ ∀u ∈ B}`. If `x ≥ 2θ/ℓ(B,θ)` then `M(B ∪ {x}) ≥ θ`.

*Proof.* `{t : ‖xt‖ < θ}` is a disjoint union of open arcs of width `2θ/x` centred at
`k/x`, separated by gaps of width `(1−2θ)/x > 0`. A closed interval of length `≥ 2θ/x`
cannot lie inside one open arc; if it meets two arcs it contains the whole safe gap
between them; if it meets at most one it contains points outside every arc. Apply to the
longest interval of `S(B,θ)`. ∎

This is the same interval-survival mechanism as THM-1004's tail step (which I found
*after* proving it — §5) and the THM-547/563/608 far-element closures; it is recorded
here in its minimal reusable form. With `θ = 2/(2N+1)` and `B = {1..N}∖{i}` (at most 12
speeds), settled LRC(≤13) gives `M(B) ≥ 1/N > θ`, so `ℓ(B,θ) > 0` and
`X₀(N,i) = ⌈2θ/ℓ⌉` is finite: **the single-far stratum `{1..N}∖{i} ∪ {x}` is decided by
the finite check `x ≤ X₀`.**

## 2. The single-far classification, N = 6..13 (PROVED)

All `X₀(N,i) ≤ 66`; every finite check was run to `max(X₀, 12N)+25` (the band above `X₀`
doubling as an empirical cross-validation of the Lemma — zero violations anywhere). The
complete list of first-gap members of the single-defect single-far stratum over all
N ∈ {6,…,13}, all defect positions i, all x (unbounded):

```text
N = 7 :  {1,2,3,4,5,7,18}          (i=6,  x=18=3(N−1))   M = 3/23  (the mediant)
N = 13:  {1,...,11,13,36}          (i=12, x=36=3(N−1))   M = 3/41  (the mediant)
ALL OTHER (N,i,x): NONE.  In particular N = 6, 8, 9, 10, 11, 12: stratum EMPTY.
```

The two members are exactly the canonical mod-30 binder-gate solutions (HYP-4516:
attainable iff `N ≡ 1 (mod 6)` and `5 ∤ 3N+2`; N = 7, 13 are the only such N in range).
**The gate's canonical solutions are therefore the UNIQUE single-far gap members for
N ≤ 13** — the gate is complete in-stratum, not merely sufficient for one family shape.

The general-defect resonance atlas at N = 12 (dip table, all i): for `7 ≤ i ≤ 12` the
resonant ladder is `x = i·m` with binder `b = 13−i` (the *complement* of the defect) and

```text
M({1..12}∖{i} ∪ {i·m}) = m/(i·m + 13 − i),   m ≥ 2
```

whose in-window condition `m(25−2i) < 2(13−i)` has **no integer solution m ≥ 2 for any
i** — every defect ladder skips the window, with the unique edge contact
`(i,m) = (12,2)` at exactly `2/25` (THM-633's rung). For `i ≤ 6` the resonances are
richer (mixed binders, e.g. i=6 runs two interleaved step-19 sub-ladders at b ∈ {7,12});
none dips into the window either. This is the structural *reason* for emptiness:
ladders skip, and everything non-resonant is absorbed above `2/25` by the Lemma.

## 3. The N = 8..11 census (VERIFIED, gates passed)

Five species per N — (A) single-defect single-far (now superseded by §2's proof);
(B) bordered dilated APs, interior borders included, spine `a+d{0..m−1}` + borders
`±1,2,3`; (C) two-defect two-far; (D) greedy descent from `Allowed(a)` residue bands at
every admissible `(val,q)`; (E) adaptive needle-repair from band-compliant residue lifts.
Validation gates: the exact evaluator reproduces ten known values (incl. Fan–Sun's
`ML(3,8,11,19) = 7/30`), and the generators **rediscover all three known members**
(5/33 via B; 3/23 and 3/41 via A) from inside their own ranges — the instrument's
detection floor is below its target (the HYP-7870 discipline).

```text
N = 8  (window width .00654): EMPTY   [B: 622,442; C: 31,584; A exhausted; E: 750]
N = 9  (        .00526): EMPTY        [B: 1,097,654; C: 51,516]
N = 10 (        .00433): EMPTY        [B: 1,631,965; C: 79,650]
N = 11 (        .00362): EMPTY        [B: 1,603,662; C: 117,975]
N = 12 (        .00308): EMPTY        [consistency with the fleet censuses — passes]
```

Honest scope: species-and-budget-bounded evidence, not proof (species D yielded nothing
within budget — its descent is expensive; the load was carried by A/B/C/E). But §2 makes
the single-far slice of each verdict PROVED, and at N=12 the two-far slice is proved
(THM-1005; replicated here in 609 checks by the nested thresholds `X₁ = ⌈(1+2θ)/ℓ(B)⌉`,
then `X₂ = ⌈2θ/ℓ(B∪{x₁})⌉` per x₁ — every intermediate base has ≤ 11 speeds, so
`M ≥ 1/12 > 2/25` keeps every level's ℓ positive).

## 4. The band law (CONJECTURE) and a falsifiable prediction

Within [6,13] the populated N are exactly `{6} ∪ {7, 13} = {6} ∪ {N ≡ 1 mod 6}`: **the
emptiness band [8,12] is precisely the stretch of N ≢ 1 (mod 6) after the wide-window
species die.** Reading: gap membership is the disjunction of per-species arithmetic
gates; the k=2 (mediant) gate opens iff `N ≡ 1 (mod 6)` and `5 ∤ 3N+2` (HYP-4516, now
in-stratum complete by §2); the deeper k ≥ 3 species (N=6's bordered 5/33, order 3) need
window width that dies by N ≈ 8. Hence for N ≥ 8:

> **Conjectured law: `W_N ≠ ∅  ⟺  N ≡ 1 (mod 6) and 5 ∤ (3N+2)`.**

Falsifiable out-of-sample: `W_N` empty for N = 14,…,18, populated at N = 19 (3/59 via
`{1..17,19,54}`), and — the sharp one — **N = 31 EMPTY** (the first `N ≡ 1 mod 6` with
`5 | 3N+2`, where mac-mini-S29 showed the canonical family degrades to the floor 1/32).
N=12 is thus not an isolated miracle but the generic case of a band; a proof of the
band's mechanism (gates + width-death) would subsume the (C)-side first-gap emptiness.

## 5. Convergence note and corollaries

- **N=12 replication (post-hoc discovery).** My single-far and two-far N=12 runs
  independently reproduce THM-1004/THM-1005 (klein-S313d/e, 2026-07-17) — same lemma,
  same interval table (e.g. `ℓ(i=5) = 7/1000` matches THM-1004's `L_5` exactly), same
  verdicts. Priority is klein's; the value added here is the double-witness, the
  cross-N extension, and the compact nested-threshold packaging (609 checks for
  radius 2). Caught by the STATEMENT-grep protocol only at writeup time because I
  searched my own thread's vocabulary ("single-far/outlier") and not the rigidity
  thread's ("Hamming") — see MISTAKE-187.
- **Corollary (N=13, for THM-1240/HYP-7840).** The single-replacement stratum of
  `{1,…,13}` contains NO family with `M ∈ (1/14, 3/41)`: its only first-gap member IS
  `3/41 = M({1..11,13,36})`. opus's bounded scans of the shapes `{1..11,13,x}`,
  `{1..12,x}` are hereby upgraded to a closed stratum, uniformly in x, for all 13
  defect positions.
- **Corollary (order data).** Attained orders `k = q − N·val` in [6,13]: k=3 at N=6,
  k=2 at N=7 and N=13, nothing else — the "achievable depth" of HYP-4542 is not
  monotone (it returns at N=13); it is gate-arithmetic, as opus-S118 argued, and the
  gates are now complete for the single-far stratum.
