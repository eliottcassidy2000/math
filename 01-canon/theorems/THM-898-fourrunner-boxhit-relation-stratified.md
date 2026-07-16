---
id: THM-898
title: THE FOUR-RUNNER BOX-HIT LAW IS RELATION-STRATIFIED — the owner-directed "uniform four-runner box-hit bound" for the residue-6 kernel EXISTS but is stratified by the quadruple's ADDITIVE RELATIONS, not scale-decaying: (i) exact strata machinery (hit = product + pair stratum [codex's D-law, |D| ≤ 12/7·(1/AB)] + remainder R, all in exact ℚ); (ii) REFUTATION: the naive uniform R = O(1/(v-triple)) is FALSE — R·minTriple grows unboundedly (8 → 492) along the kernel's (1,5,2,4)-residue quadruples, because each carries the additive resonance v₁+v₂ = v₃+v₄ whose (1,1,−1,−1) lattice vector contributes an O(1) PLATEAU (~5·10⁻³, scale-independent; certified across v = 8..54 and by relation-typed controls: fixed relation ⟹ R → const, taller relations ⟹ smaller R); (iii) THE REPLACEMENT: |R| ≤ Σ_{primitive k ⊥ v, no sub-pair} c_B(k) + O(1/v), c_B(k) = the box-Fourier product over the relation — uniform PER RELATION-STRATUM; (iv) the kernel reading: codex-THM-891's A15+A24 target = [pair D-laws] + [per-relation O(1) constants on the E₃-degenerate quadruples] + [O(1/v) off them] — the residue-6 crux is ADDITIVE-ENERGY-stratified, and THM-730 (AP uniquely maximizes additive triples) is the counting tool for the degenerate stratum. The week's frame-heuristic strikes a seventh time: the naive uniform constant privileged the relation-free frame
status: (i) exact machinery + battery banked; (ii) REFUTATION certified (exact ℚ, growth table + relation attribution); (iii) stated with the lattice/Fourier proof shape, constants named as the next computation (c_B(k) for the height-4 relations); (iv) the kernel consequence composed with codex's framework honestly (their normalization to be reconciled by them or next session)
source: boxeph-2026-07-16-S33 (owner: turn the residue-six crux into a uniform four-runner box-hit bound)
depends_on: [codex-S17 THM-891 (the kernel + pair D-law), THM-730 (the E3 counting tool)]
script: 04-computation/lrc14_fourrunner_boxhit_boxeph_S33.py -> 05-knowledge/results/lrc14_fourrunner_boxhit_boxeph_S33.out
---

# THM-898 — the relation-stratified four-runner law

hit(v, B) − ∏|Bᵢ|/7 − Σ_{pairs}[D/(vᵢvⱼ)-terms] = R(v, B), computed exactly. The
kernel battery (residues (1,5,2,4), near-AP): R·minTriple = 8.3, 30.1, 63.9, 132.5,
264.0, 491.8 — UNBOUNDED: no 1/v³ (nor any power) uniform law. Attribution: every
plateau quadruple satisfies v₁+v₂ = v₃+v₄; for a FIXED primitive relation k ⊥ v the
remainder tends to a scale-independent constant c_B(k) (measured ≈ 5·10⁻³ for
(1,1,−1,−1) at singleton boxes; taller relations measurably smaller: 7·10⁻⁴–1.6·10⁻³
at height-3/4 with a coefficient 2). The uniform law is per-stratum:
> |R| ≤ Σ_{primitive k ⊥ v, no sub-pair relation} c_B(k) + O(1/min v),
with c_B(k) the Fourier mass of the boxes along the relation (explicit; evaluation =
named next step). For the residue-6 kernel: A15+A24 inherits this stratification —
its four-runner content is pairwise D-data plus O(1) contributions EXACTLY on the
additively-degenerate quadruples, which THM-730 counts (AP-extremal). The crux is not
"a constant" but "a constant per additive relation" — the seventh instance of the
week's frame-correction pattern.

## Evidence log
- [x] exact strata + battery; refutation growth table; relation attribution table
- [ ] evaluate c_B(k) exactly for height ≤ 4 relations (closed form via box Fourier)
- [ ] compose with codex's normalization to close K6({1,5}) = K6({2,4}) = −12
