# LRC(14) FINISH MAP — 2026-07-11 (klein-S258)

The definitive current state: what is PROVED, what is VERIFIED-not-proved, and the single
remaining mathematical statement — in the two independent route-forms it takes. Supersedes the
scattered residual notes in OPEN-QUESTIONS. Synthesis of the converged fleet state (klein S254–257,
mac-mini cont.42–44, opus S234–237, kps cont.36–42).

## The top-level reduction (all PROVED / SETTLED)

LRC(14): every 13 nonzero integer speeds have a `1/14`-lonely time. The chain:

1. **Non-covering case — SETTLED.** A speed set omitting a multiple of some `q ∈ {2,…,14}` has the
   explicit witness `t = 1/q` (`M ≥ 1/q ≥ 1/14`). Uses LRC(≤13) (cited per project policy). So
   LRC(14) reduces to **covering** sets (a multiple of every `q ≤ 14`). Covering ⟹ **divisor-complete**
   (THM-366, PROVED). Divisor-complete = 8.5% of primitive 13-families.
2. **Density floor — CLOSED** (THM-661/663): the covering good-set measure `ρ* ≥ m_P = 14249/252252
   > 0` (the moment-LP `B_d` bound, `d ≤ 4`, diameter-free). PROVED per-shape; the uniform floor is
   `[compact exact check] + [decorrelation tail]`.
3. **Realization — supply constructive** (klein THM-693–698, S242–252, kernel-pure Lean): a positive
   floor yields an explicit bounded-denominator lonely time; the ≤7-arcs pigeonhole (SmallClusterFull)
   is PROVED in Lean (LRCSevenGapRigidity). LRC(14) = `[THM-661] + [certificate supply]`, both in Lean.

So the ENTIRE remaining content is one statement about the **covering / divisor-complete** class,
and it has **two independent proof routes** — proving EITHER finishes LRC(14).

## Route A (moment / density) — the compact extremality

The moment-ladder base (the uniform floor's compact part): the base rows are `k=8` (deg-3) and
`k=9` (deg-2); `k ≥ 10` inherit via the eigen-transfer (THM-710, PROVED). Each base row reduces
(klein THM-717 + S257 tail) to:

> **[A] `J(consec_k)` is the minimum of `J = E[N(7−N)]` over all `k`-clusters** (`k = 8, 9`).
> `J(consec-9) = 4465/882`, `J(consec-8) = 291/49`.

- **PROVED around it:** the tail (klein S257): every WIDE cluster reduces to its densest compact
  sub-cluster via the two-scale/multi-scale limit (THM-687/688); the min two-scale limit `≥`
  compact-min at every level `k=5..9` (far elements RAISE `J`); the descent bottoms at `k ≤ 7`
  where SmallClusterFull (Lean) gives `ν=1`. The tail cancellation is confined to the small BUNCH
  term (THM-717), whose max is `18/(7k−6)` at the mod-7 pole (klein S256).
- **VERIFIED-not-proved:** [A] itself — that consec/AP is the compact minimizer. Exhaustive over
  primitive clusters `diam ≤ 20` (mac-mini, 19448 families) + all algebraic extremal candidates
  (klein extremal-candidate evaluator, S255) + adversarial.
- **MECHANISM (mac-mini cont.44):** [A] IS the three-gap theorem in coverage form — the Steinhaus
  orbit `{jα}` has `≤3` gap lengths (maximally coverage-efficient), so the AP is the best coverer
  (`p0(consec)/p0(iid) = 7–24×`); random phases clump (coupon-collector). Not yet a proof.

## Route B (divisibility / bounded modulus) — the DC clearing

Covering ⟹ divisor-complete (THM-366). The band-edge margin lemma (opus S235, PROVED): if a family
clears (`bandCount = 0`) at a modulus `q` with `14 ∤ q`, then `M ≥ ⌈q/14⌉/q > 1/14` (strict). So:

> **[B] every primitive divisor-complete family clears at some `q ∈ [15, 31]` with `14 ∤ q`.**

- **PROVED around it:** band-edge margin (the strict `> 1/14` is free from clearing); tight-locus
  characterization (`M = 1/14 ⟺ clears only at multiples of 14`; the tight AP `{1..13}` is the
  UNIQUE primitive tight interval, opus S237); AP sub-case (opus S236, three-gap — the `~1%` corner);
  diameter-freeness (kps `B5_congr_mod`, Lean — [B] is residue-periodic mod `lcm(15..31)`, hence a
  FINITE statement).
- **VERIFIED-not-proved:** [B] on the SPREAD bulk (99% of DC families, longest-AP `≤ 7`). Verified
  across the whole DC class, adversarial `q ≤ 29`, `0` exceptions (opus S237). kps cont.35 CRT-factored
  the residue check: primes `11..43` cover `99.97%` (independent per-prime), a measure-`3e-4`
  composite-core covered by composites `[14,42]`.
- **MECHANISM:** clearing at `q` for an AP `{a+jd}` is "an AP `{(a+jd)p mod q}` on `ℤ/q` avoiding
  the danger arc" — three-gap on `ℤ/q` (opus S236). For spread families it is bounded-modulus
  anti-concentration.

## Route B — SHARPENED (klein-S259): THM-718 + the inverse-theorem framing

klein THM-718 (PROVED) makes "clears at prime `q`" EXACT: `clearing_count(v,q) = (q−1) −
|{±j·vᵢ mod q : 1≤j≤m}|` (`m = ⌈q/14⌉−1`), so **clears at `q` ⟺ the dilated-±-speed set misses a
residue mod `q`** — a covering number, the exact form of "bad coverer clears." And the "clears"
side is now sharply structural (verified): the ONLY window-covers are the TIGHT families
(`{1..13}`, GW, dilates), which are all NON-divisor-complete (no multiple of 14) — so **DC ⟹ has a
multiple of 14 ⟹ not tight ⟹ clears** (0/16328 DC window-covers). Route B's remaining gap = the
inverse theorem `window-cover ⟹ tight` (window-completeness); the tight list is characterized
(THM-612/708/709), DC/tight disjointness is proved, and THM-718 quantifies "clears."

## The unification and the honest assessment

**[A] and [B] are the same phenomenon** — `{kα}`-progressions (three-gap configs) are the extremals:
coverage-optimal (the `J`/`p0` min in Route A) and their mirror is bunching/resonance-optimal (the
mod-7 pole, klein S256). Route A's hardness is "AP is the best coverer"; Route B's is "AP mod q
avoids the danger arc." Both are the Steinhaus/three-gap regularity of `{jα}`.

**Finishability (honest):**
- Neither [A] nor [B] is currently proved; both are the genuine LRC-adjacent core, verified across
  their class over 40+ fleet sub-sessions.
- **Both are FINITE/BOUNDED in principle.** [B] is diameter-free (residue-periodic, `q ≤ 31`) — a
  genuinely finite (if astronomically large, CRT-factored) check. [A] reduces (via the proved tail)
  to a compact check over bounded clusters.
- **The shortest paths to DONE:**
  1. *Prove [B] on the spread bulk* — a bounded-modulus anti-concentration: a spread divisor-complete
     family cannot have some `v_i` in the danger arc for ALL `~390` pairs `(q,p)`, `q ∈ [15,31]`.
     This is the most self-contained; the danger arc has width `~q/7`, and DC-ness (multiple of
     8,9,11,13) FORCES the spread that guarantees a clearing `p`. Candidate tools: character-sum /
     Weyl bound over the bounded `q`; the AP sub-case is done (three-gap).
  2. *Complete the finite verification of [B]* + cite it (the CRT-factored residue check to
     completeness — per-prime conditions `q = 17..31` rigorous + the `3e-4` composite-core exhausted).
  3. *Prove [A]* via three-gap coverage-optimality (Route A) — deeper (compares AP to non-AP).

**Recommendation:** Route B, path 1 (prove the spread-bulk anti-concentration) or path 2 (complete
+ cite the finite residue check) is the shortest to a complete LRC(14) proof, because [B] is bounded-
modulus and diameter-free (genuinely finite), whereas [A] carries the compact-extremality + tail.
The AP corner of [B] is done (three-gap); the spread bulk is the target.

## Lean state (transcription)

The reduction chain, the density floor citation, the ≤7-arcs pigeonhole (SmallClusterFull), the
certificate supply, and the band-edge margin lemma are in kernel-pure Lean. What remains for a fully
Lean-checked LRC(14): whichever of [A]/[B] is proved, transcribed, plus the finite check it rests on.

*Files: this doc. Sources: THM-366/610/661/663/687/688/710/717, klein S254–257, mac-mini cont.42–44,
opus S234–237, kps cont.36–42, LRCSevenGapRigidity.lean.*

## Appendix — independent exact-B5 verification of Route B [B] (klein-S258)

Confirmed [B] with my own B5 clean-ruler machinery (THM-671/707/712), independent of opus's census:
over 3000 primitive SPREAD divisor-complete 13-families, EVERY family clears (`bandCount = 0` at
some `p`) at a modulus `q ∈ [15, 29]`, `14 ∤ q` — **0 failures**. Clearing-modulus distribution:
`q=15:644, 16:442, 17:904, 18:190, 19:498, 20:132, 21:92, 22:47, 23:35, 24:4, 25:7, 26:1, 27:1,
29:3`; worst-case `q = 29` ⟹ `M ≥ ⌈29/14⌉/29 = 3/29 ≈ 0.1034 > 1/14`. So Route B = "every DC
family has a B5-clean ruler at bounded `q ≤ 29`" — exactly my clean-ruler certificate (THM-671),
tied to opus's band-edge margin lemma. The prime `q = 17` dominates (30% of families), then `15, 19`.

**This localizes Route B's remaining proof to:** a bounded-modulus anti-concentration — "the 13
residues of a spread DC family cannot hit every scaled danger-arc `danger(q)·s` for all `q ∈
[15,29]`." The AP corner is done (three-gap, opus S236); the spread bulk (99%) is the target, and
it is my B5 framework's `B5 > 0` at bounded `q`. Files:
`04-computation/lrc14_route_B_dc_clearing_verify_klein_S258.py`.
