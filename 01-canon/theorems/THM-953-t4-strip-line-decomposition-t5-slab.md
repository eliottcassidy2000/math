---
id: THM-953
title: THE STRIP-LINE DECOMPOSITION — the T₄ resonance strip priced, the T₅ slab schemed; THM-946(IV)'s two remaining stones worked at coarse-constant rigor. SETTING (support 4, speeds v₁<v₂<v₃<v₄): fixing outer coordinates (u,t), the inner two-pole separation is Δ ∝ |k|/(v₁′v₂′) with k = v₃u + v₄t; the RESONANCE STRIP is {|k| ≤ K₀ ≈ v₁′v₂′} where the inner sum does not decay. (I) THE STRUCTURE THEOREM: the strip is a UNION OF ≤ 2K₀+1 CONGRUENCE COSET LINES — for each value k, the outer solutions (u,t) form a single 1-D progression (u,t) = (u_k,t_k) + m(v₄″,−v₃″) — so the strip is not a 2-D region but a k-indexed line family, the exact analog one rank up of THM-952's congruence orbits (the same reason purely analytic strip bounds fail: the strip's population is arithmetic); (II) THE STRIP PRICE: on each k-line the outer weight decays as 1/(m²v₃″v₄″) away from the origin, and the DISSOCIATION FLOOR covers the origin (small |u|,|t| with small inner near-solutions would constitute a small support-≤4 relation — forbidden; surviving points carry a coordinate > H): strip mass ≤ C·L²/H; (III) THE COMPLEMENT: |k| > K₀ gives inner Δ ≥ |k|/(v₁′v₂′) and THM-946(I)'s corrected atom decays; summing over the k-line family: ≤ C·L³/H. TOTAL: **T₄(H) ≤ C₄L³/H at coarse-constant rigor** (case tree at THM-952's level; audit explicitly invited — the constants are NOT sharp and the gcd/zero-coordinate bookkeeping follows THM-952's pattern). REFEREED (exact enumeration, box 250, three dissociated quadruples, H ∈ {10,40,160}): the strip is ≤ 22% of the mass and decays FASTER than 1/H (the floor eats its origin); the total envelope T₄·H/L³ is BOUNDED (4×10⁻⁷–8×10⁻⁶, no growth); (IV) THE T₅ SLAB SCHEME: the slab {|v₃u+v₄t+v₅r| ≤ K₀} is a k-indexed family of COSET PLANES (rank-2); each plane prices by the strip result one rank down (the plane's own strip/complement recursion) plus the same dissociation floor: T₅ ≤ C₅L⁴/H — SCHEME with the recursion precise, the plane-level bookkeeping named as the remaining audit surface (NOT claimed proved); (V) STATUS OF THE EXHAUSTION ROUTE (THM-946(V)): T₃ unconditional (THM-952 + codex's atom), T₄ at coarse-constant rigor (this file, audit invited), T₅ schemed — the route is one audit and one plane-bookkeeping pass from closing
status: (I) proved (elementary: k fibers the outer lattice into progressions); (II)-(III) proved at coarse-constant rigor with the case tree following THM-952 (constants not tracked sharply; codex-style audit INVITED per the MISTAKE-157/THM-944 lessons — this file's claims are calibrated to survive one); (IV) scheme + referee; referee: t4_strip_referee_kps_S128c41.py (envelopes bounded)
source: kind-pasteur-2026-07-17-S128 (cont.41; owner: work the T4 resonance strip and the T5 slab)
depends_on:
  - THM-946 (codex: the corrected atom + the strip/slab problem statement)
  - THM-952 (the congruence-orbit + dissociation-floor mechanism this lifts one rank)
related:
  - THM-935/948 (the E_s frame and exact packet audits)
  - the exhaustion route THM-946(V)
script: 04-computation/t4_strip_referee_kps_S128c41.py -> 05-knowledge/results/t4_strip_referee_kps_S128c41.out
---

# THM-953 — the strip-line decomposition; the slab scheme

## (I) The strip is a line family

k = v₃u + v₄t fibers Z² into progressions: for each k (divisible by g₃₄), the solutions
form one coset line with primitive direction (v₄″,−v₃″). The strip |k| ≤ K₀ is therefore
≤ 2K₀+1 lines — its population is arithmetic, not planar. (This is why analytic bounds
treating the strip as a 2-D region overpay, and why no purely analytic estimate closed it.)

## (II)–(III) Pricing

Per line: outer weight 1/(|u||t|) ≤ C/(m²v₃″v₄″) off the origin (ζ(2)-summable); at the
origin the dissociation floor applies — a strip point with all four coordinates ≤ H IS a
small relation. Surviving origin-adjacent points carry a coordinate > H, contributing the
1/H factor. Complement: the corrected atom with Δ ≥ |k|/(v₁′v₂′); Σ_k 1/(1+Δ_k) costs one
L. Assembly: T₄ ≤ C₄L³/H. Referee: envelope bounded across all tested quadruples/horizons;
the strip portion decays faster than 1/H as the floor predicts.

## (IV) The slab scheme (T₅ — NOT claimed proved)

The slab is a k-indexed family of rank-2 coset planes; each plane recursively splits into
its own strip (lines, by the second linear form) and complement (atom-decaying). The
recursion terminates in two-pole atoms with all separations explicit. The remaining work
is the plane-level gcd/zero-coordinate bookkeeping and constant assembly — one focused
audit-grade pass.


## The hardened composition (cont.44 — audit-grade structure)

**The composed-atom theorem.** Both levels of T₄ are two-pole sums: the INNER slice at
separation Δ^in_k = |k|/(v₁′v₂′) (THM-946(I) verbatim), and — the hardening's key — the
OUTER k-line itself: u(m) = u_k + mv₄″, t(m) = t_k − mv₃″ are two linear forms with poles
separated by Δ^out_k = |k|/(v₃″v₄″). Hence

> T₄ ≤ Σ_{g₃₄|k, k≠0} (1/π²)·Atom(v₂′,v₁′; Δ^in_k) · (1/π²)·Atom(v₄″,v₃″; Δ^out_k)

with Atom = THM-946's corrected lemma (1). REFEREED (t4_composed_atom_referee, box 200,
3 quadruples): true ≤ composed bound everywhere (the H-free structural check; the floored
version supplies the 1/H exactly as in T₃ — cont.41's envelope measurement).

**The floor cases** mirror T₃'s hardened tree on each level independently: which of the
four coordinates exceeds H splits into inner-pair cases (THM-952's subcases A/B/C verbatim
on (h₁,h₂)) and outer-pair cases (the same subcases on (u,t), with the congruence-averaging
lemma applied to the OUTER orbit — r_out(k) = the least-abs residue of the k-line's
near-pole in u). The averaging lemma is rank-free and covers both levels; no new mechanism
appears at support 4. The k-sum converges against both atoms' (1+Δ)⁻¹ decays with the two
L factors, giving the L³ of the assembled bound.

With this, T₄'s tree matches T₃'s audit grade structurally; the remaining delta to a
fully-audited T₄ is transcribing the eight (inner × outer) floor combinations with their
constants — mechanical, each a copy of a THM-952 subcase.

## Named next
- The codex-style audit of (II)-(III)'s case tree (invited; claims calibrated for it).
- The T₅ plane bookkeeping pass.
- The cluster-gap leaf: `cite_cluster7_lonely` (LRCSevenGap.lean) is a CITATION Prop —
  "frozen arcs can't perfectly pack" — an open math leaf of the covering program's pages,
  distinct from the cluster-gcd modules (which are kernel-pure; their grep-"sorry" hits
  are doc comments). Sharpened as the next covering-page target.

## Evidence log
- [x] strip/complement referee: 3 quadruples × 3 horizons, envelopes bounded
- [ ] audit of the case tree; T₅ plane pass; cluster7 citation discharge
