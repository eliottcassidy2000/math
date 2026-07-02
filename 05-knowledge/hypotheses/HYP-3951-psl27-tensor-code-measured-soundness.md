# HYP-3951: The explicit PSL_2(F_7) tensor code — MEASURED soundness s vs the 1/36 family

**Status:** VERIFIED (code built exactly; soundness measured; 1/36 identified as the query-grid quantum,
NOT the extremal soundness) — kind-pasteur-2026-07-01-S28.

## The explicit code (definition)
G = PSL_2(F_7) (|G| = 168). A = symmetric set of ORDER-7 elements (inverse pairs), B = symmetric set of
ORDER-3/4 elements. Order-disjointness => zero degenerate squares: **faces = triples (g,a,b) under the
4-way identification (g,a,b)~(ag,a^-1,b)~(gb,a,b^-1)~(agb,a^-1,b^-1)**, so F = |G||A||B|/4 EXACTLY
(verified all configs) and every vertex link is the full |A| x |B| grid of distinct faces (fixes
opus-S31's corner-set construction, which conflated/dropped squares). Inner codes C_A (columns),
C_B (rows); **C = {f : F -> F_2 : every vertex's grid view lies in C_A ⊗ C_B}**. Tester: read a uniform
vertex's |A||B|-bit view, accept iff in the tensor code. Soundness of f: rho(f) = [viol(f)/V]/[dist(f,C)/F].

## Measured results (psl27_tensor_code_measured_soundness_kps.out)

| config | comp | F | dim C | d_ub | s_meas | note |
|---|---|---|---|---|---|---|
| SURFACE rep2⊗rep2 (deg 2) | **8** | 168 | 8 | 21 | 0.26 | disconnected (mac-mini regime); s decays with support |
| par4⊗par4 (deg 4) | 1 | 672 | 172 | 12 | 1.23 | connected, λ₂=3.414 < Ramanujan 3.464 |
| par6⊗par6 (deg 6 = φ(7)) | 1 | 1512 | 678 | **6-8** | 2.14 | DISTANCE COLLAPSE via torsion tubes |
| rep6⊗rep6 (deg 6) | 1 | 1512 | 1 | 1512 | 2.05 | trivial code (all-ones) |
| ham8⊗ham8 (deg 8, [8,4,4]) | 1 | 2688 | **16** | **936** | 3.13 | the genuinely good one |

(weight-1/2 rows exhaustive = true coset leaders; hill-climb rows use wt as dist proxy — conservative-low.)

## Findings

1. **The disconnection and the upgrade, quantified in one table.** SURFACE (single-gen pairs) has 8
   components (= double cosets <a>g<b>, each the Frobenius-21) and the worst measured s (0.26 at
   support 69, decaying — the poly-soundness leak). Connected tensor configs hold s = 1.2-3.1, O(1),
   INCREASING with degree/inner strength. Both mac-mini's correction (HYP-3824) and opus's mechanism
   (HYP-3823) are now measured, not just certified.

2. **TORSION TUBES = the distance obstruction at small q (the sharp finite-q mechanism).** The
   par⊗par code (inner distance 2) has weight-6/8 codewords: closed tubes over SHORT RELATIONS in the
   generators — b^3 = 1 gives a 3-cycle x inverse-pair tube (6 faces); a commuting A-pair gives a
   4-cycle tube (8 faces; inspected explicitly: 8 faces on 8 vertices, each view a 2x2 rectangle).
   These are klein-S86's b1-cosystoles materialized as codewords. **Girth(Cayley) <= element order**
   at small q, so inner distance 2 => code distance <= 2·(min torsion) — DELLM's "one q cannot exhibit
   growing distance" caveat, with the exact mechanism. **Inner distance >= 3 (their actual requirement)
   kills precisely the torsion tubes**: ham8⊗ham8 jumps to d_ub = 936 (relative distance ~0.35!).

3. **Explaining opus's 16.** The surface complex's dim H^1 = 16 (opus-S30) decomposes as 8 components
   x H^1(torus) = 2: each Frobenius-21 component has V=21, E=42, F=21, chi = 0 — a TORUS. The
   ham8⊗ham8 tensor code ALSO has dim exactly 16 (= |A| + |B|?) on the connected deg-8 complex — a
   numerical echo worth a look (lead, not a claim).

4. **THE 1/36 VERDICT (the session's question).** The measured EXTREMAL soundness is O(1) (0.26-3.1)
   — it does NOT recover 1/36 = 0.0278 as its value, and it shouldn't: an LTC's soundness is designed
   to be a large constant. **Where 1/36 does appear, exactly: the query-grid quantum.** At the natural
   degree phi(7) = 6 config, the single-defect ratio is s(w=1) = |A|·|B| = 36.0000 exactly, i.e. one
   flipped face carries 1/36 of the tester's 6x6 grid resolution. The analytic floor 1/36 = (1/6)^2 has
   the same anatomy: 6 = inner sectors (safe arcs per runner at band 1/14), and the r=2 moment
   relaxation resolves the 2-far configuration on a 6x6 sector-pair grid. **So the soundness/measure
   mirror is a RESOLUTION identity (both floors = one cell of the phi(7) x phi(7) grid), not an
   extremal-value identity.** The mirror of the census constant (0.0323 = 1.16 x 1/36) is the
   TIGHTNESS of the grid bound, which on the code side corresponds to how close the true minimal
   detection sits above the single-cell quantum — structural analogy, now with the exact seam located.

## Honest scope
Single q = 7 instance; hill-climb s values are lower estimates (wt >= dist beyond coset-leader range);
d_ub values are randomized-bank upper bounds (weight-1/2 exhaustive rows are exact). The DELLM
asymptotic theorem is not (and cannot be) reproduced at one q; what IS established: the explicit code,
its exact complex (F = |G||A||B|/4), the measured soundness ordering, the torsion-tube mechanism, and
the location of 1/36 in the construction.

## Artifacts
- 04-computation/psl27_tensor_code_measured_soundness_kps.py (+ .out)

## Depends on / relates to
HYP-3823 (opus), HYP-3824 (mac-mini disconnection — confirmed + used), HYP-3832 (klein cosystoles =
our torsion tubes), kps-S25/S27 (sqrt-p faces, census), HYP-3950 (the analytic side), DELLM
(Annals 2026 203-2 / STOC 2022).
