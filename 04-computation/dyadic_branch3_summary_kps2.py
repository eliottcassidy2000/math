"""
kind-pasteur-2026-06-09-S2 : BRANCH III consolidated summary writer.
Collates the session's frontier tables and verdicts into one results file.
(All numbers produced and independently verified by the sibling _kps2 scripts.)

Output -> 05-knowledge/results/dyadic_branch3_summary_kps2.out
"""
SUMMARY = r"""
====================================================================================================
BRANCH III -- THE DYADIC-GAP HUNT (HYP-2359, Erdos-Gyarfas / Erdos 64)  kind-pasteur-2026-06-09-S2
====================================================================================================

LITERATURE ANCHORS (verified against primary sources this session):
  * Markstrom 2004 (Congr. Numer. 171; PDF from author site, parsed): ALL cubic graphs on
    <29 vertices contain a cycle of length 4, 8 or 16 => cubic counterexample needs >=30
    vertices. Smallest cubic graphs with no C4 and no C8: n=24, EXACTLY 4 classes
    (one planar); counts 4 @ n=24, 23 @ n=26, 251 @ n=28 (his Table 3).
  * Royle: all C4-free min-deg-3-structure graphs on <16 vertices checked, no counterexample.
  * Liu-Montgomery [LiMo20]: min degree >= some absolute constant => power-of-2 cycle.
  * Erdos-Gyarfas 1995, $100/$50. Status: OPEN (erdosproblems.com #64).

CHECKER VALIDATION (dyadic_cycle_checker_kps2): exact spectra match Markstrom's psi values:
  Petersen 57 (c5=12 c6=10 c8=15 c9=20), Heawood 213 (c8=21), McGee 5608, Tutte-Coxeter 41400.
  All counts independently reproduced by networkx.simple_cycles on every key specimen.

CORRECTION OF S710: McGee (n=24, girth 7) HAS 34 eight-cycles (c7=32, c8=34, |Aut|=32
verified). S710's 'McGee -> C16' read the first power-of-2 cycle FOUND in enumeration
order, not the smallest. Petersen has 15 C8s (not 10).

FRONTIER TABLE 1 -- girth >= 5 cubic, min #C8 found (SA, 8 restarts x 60K + 16 x 250K deep):
   n=20: 6   (16/16 deep restarts hit exactly 6)
   n=22: 3   (16/16)
   n=24: 3   (16/16)  [PROVEN >0: all 4 C4+C8-free classes at n=24 have girth 3]
   n=26: 1   [PROVEN >0: all 23 C4+C8-free classes at n=26 have girth 3 -- census complete]
   n=28: 0   <-- FIRST girth-5 C8-free cubic graphs; >=2 iso classes, 3-connected,
                 c16=614/731; independently verified
   n=30..40: 0 (increasingly easy)
  => THE MINIMUM ORDER OF A GIRTH->=5 CUBIC GRAPH WITH NO 8-CYCLE IS EXACTLY 28
     (conditional only on Markstrom's exhaustive Table-3 counts).

FRONTIER TABLE 2 -- girth >= 5, C8-free, min #C16 found (stage-2 descent):
   n=28: 614   n=30: 736   n=32: 780   n=34: 855   n=36: 879   n=38: 961   n=40: 970
  => NO {8,16}-avoider; C16 mass GROWS along the frontier. The {girth>=5, C8-free}
     stratum is nearly frozen under edge switches (>=4000 consecutive rejected moves typ.).

FRONTIER TABLE 3 -- C4-free (triangles allowed; the true counterexample regime):
   n=22: min c8 = 2 (matches Markstrom: impossible to reach 0 below 24)
   n=24: 0; ALL 4 classes rediscovered (girths 3,3,3,3; planar one = Markstrom graph,
         |Aut|=3, c16=228; only power-of-2 cycles are C16 in all four: 207/228/315/330)
   n=26: 0; ALL 23 classes collected (SA restarts + in-stratum walk), girth histogram:
         23 x girth-3, c16 range [161, 454], none planar
   n=30: min #C16 with c4=c8=0: 271   (n=30 is the first OPEN order beyond Markstrom)
   n=32: min #C16 with c4=c8=0: 379
  => no {4,8,16}-avoider found at n=30/32; quantitative gap to a counterexample remains huge.

THE GIRTH LADDER OF C8-AVOIDANCE (min order of a cubic C8-free graph with girth >= g):
   g=3 (C4-free): 24  EXACT (Markstrom exhaustive; 4 classes, all girth 3)
   g=5:           28  EXACT (this session; see Table 1)
   g=6:           <= 32 (specimen c16=925 c32=87; n=28,30,34 SA floors at 1)
   g=7:           > 46 by SA (floors 18/13/7/2/1 at n=30/34/38/42/46); <= 58 trivially
   g=8:           IMPOSSIBLE at girth exactly 8 (girth-8 means an 8-cycle exists);
                  girth>=8 C8-free <=> girth>=9 => n=58 EXACT ((3,9)-cage bound)

HIGH-GIRTH CAGES (data: Brouwer aeb.win.tue.nl/graphs/cages, 24 graphs): every (3,9)-cage
(n=58, all 18) has c16 in [1923,2193]; (3,10)-cages (n=70): 3855/3855/3900; (3,11)-cage
(n=112): 1956; (3,12)-cage (n=126): 3780. NO high-girth cage avoids C16: the first dyadic
window above the girth is realized in EVERY known cubic cage (g5,6,7 -> C8; g8 -> C8 = girth;
g9..12 -> C16).

GENERALIZED PETERSEN: all 380 GP(n,k), 3<=n<=40: the ONLY C8-free one is GP(3,1) (6
vertices -- too small for an 8-cycle). EVERY GP(n,k) with 4<=n<=40 contains an 8-cycle.

FLOWER SNARKS: J5/J7/J9 have c8 = 15/21/27 (all carry C8).

CIRCULANTS Cay(Z_n,{+-a,n/2}), n<=80: 779 checked; every CONNECTED one on >=8 vertices
has an 8-cycle (C8-free hits are all disconnected or on 6 vertices).

INTERPRETATION (HYP-2359): within n<=40 the dyadic gaps {8,16} cannot be opened
simultaneously even in the most favorable (girth-5..7) strata; C16 abundance grows as C8
is suppressed. Structured families (GP, circulant, dihedral Cayley, snarks, cages) are
uniformly WORSE than SA-optimized graphs at avoiding C8. Everything found is consistent
with Erdos-Gyarfas; the quantitative 'rigidity' of the C8-free stratum is the new signal.

====================================================================================================
PART 2 -- COMPLETION SUB-SESSION (same day; the original run was interrupted during the
dihedral Cayley scan). New scripts: dyadic_cayley_b2b4_kps2, dyadic_exoo_h7_k4_kps2,
dyadic_exoo_g78_kps2, dyadic_dihedral_c32_kps2.
====================================================================================================

LITERATURE LOCKED (exact, fetched this sub-session):
  * Markstrom, K., 'Extremal graphs for some problems on cycles in graphs',
    Congressus Numerantium 171 (2004) 179-192 (Wikipedia pagination; Exoo's reference
    list misprints it as 2179-2192). Author PDF (abel.math.umu.se/~klasm/Uppsatser/
    cycex.pdf) was parsed by Part 1 but the host now refuses connections (ECONNREFUSED)
    -- Table 1/3 numbers stand on Part 1's parse plus our recomputations.
  * Wikipedia 'Erdos-Gyarfas conjecture' (citing Royle + Markstrom searches): any
    counterexample needs >= 17 vertices; any CUBIC counterexample needs >= 30 vertices;
    Markstrom found four graphs on 24 vertices whose only power-of-two cycles are
    16-cycles, one of them planar. (CORRECTS Part 1's 'Royle <16' phrasing.)
  * Exoo, G., 'Three Graphs and the Erdos-Gyarfas Conjecture', arXiv:1403.5636 (2014):
    define f(k) = min order of a cubic graph with no 2^m-cycle for all m <= k. Then
    f(2)=10 (3 graphs, incl. Petersen -- matches our mindeg3 census exactly), f(3)=24
    (Markstrom's 4), 54 <= f(4) <= 78 (lower bound = UNPUBLISHED Markstrom computation;
    upper = Exoo's G78), f(5) <= 450 (Tutte-Coxeter blown up by H15); plus G420 =
    3-connected cubic PLANAR {4,8,16}-free graph (buckyball with H7 gadgets), which
    forces m >= 5 in the Heckman-Krakovski planar bound.
  * Heckman-Krakovski, Electron. J. Combin. 20(2) (2013) P7: every 3-connected cubic
    planar graph has a 2^m-cycle with m <= 7.
  * Liu-Montgomery, arXiv:2010.15802 (J. AMS 2023): there is d such that average degree
    >= d forces a power-of-2 cycle (dense regime of Erdos-Gyarfas SOLVED; min-deg-3 open).
  * Carr, A., arXiv:2605.22844 (May 2026): >= 4/7 of the vertices of any minimal
    counterexample have degree exactly 3; a regular minimal counterexample must be cubic.
  => CONTEXT FOR TABLE 3: f(4) >= 54 means NO {4,8,16}-avoider exists at n <= 40; our SA
     floors (c16=271 @ n=30 etc.) probe exactly the regime the bound says is empty, and
     the first open {4,8,16} territory is 54 <= n <= 76 (Exoo holds 78).

EXOO H7 / MARKSTROM PLANAR GRAPH TRIPLE-LOCK (dyadic_exoo_h7_k4_kps2):
  H7 verified (9 edges, deg seq 2,2,2,3,3,3,3; cycles c3=2,c5=1,c6=2,c7=1 -- no 2^m-cycle;
  d(u,v)=d(u,w)=2, d(v,w)=3). K4 with three H7s + one K3: of all 216 attachment wirings,
  EXACTLY 16 are {C4,C8}-free and they form ONE iso class = planar, |Aut|=3, c3=7,
  c16=228 = ISOMORPHIC to the planar n=24 class our SA hunt rediscovered blind.
  Literature construction, SA search, and explicit gadget assembly all agree.

McGEE RECHECK (dyadic_cayley_b2b4_kps2 B0): networkx.simple_cycles independently gives
  c7=32, c8=34, c9=16, girth 7, |Aut|=32 -- the S710 correction stands twice-verified.

CAYLEY COMPLETION (dyadic_cayley_b2b4_kps2 B2-B4) -- REVISES Part 1's 'uniformly worse':
  * B2 dihedral D_m, three reflections {s, s r^j, s r^k} (1841 gap-classes, 2m <= 80):
    672 C8-free, 638 connected specimens on >= 8 vertices. TWO clean strata:
    - girth-6 stratum: 577 connected C4+C8-FREE cubic Cayley graphs, orders 38-80
      (smallest: D_19 refl(0,1,8), n=38), c16 in [1739, 3159] -- {4,8}-avoidance is
      GENERIC in this bipartite family, but the C16 mass is enormous and C32 exists
      (spot checks D_19(0,1,8), D_22(0,1,5)).
    - even-m stratum D_m refl(0,j,m/2): 54 specimens with c8 = c16 = 0 (girth 4),
      orders 36..80: dyadic spectrum EXACTLY {4, 32} (C32 exists in every one checked,
      dyadic_dihedral_c32_kps2; n < 64 excludes C64). A min-deg-3 family whose only
      dyadic contact below 32 is C4.
  * B3 dihedral {s, r^a, r^-a} (380 checked): every connected one on >= 8 vertices has C8.
  * B4 Z2 x Zm (190 checked): same -- all C8-free hits are disconnected.
  => the bipartite 3-reflection family is the ONLY structured family found that opens
     the {4,8} gap; it pays with c16 ~ 2000-3000, mirroring the SA frontier's C16 wall.

EXOO G78 RECONSTRUCTED + NEW DYADIC DATA (dyadic_exoo_g78_kps2):
  Petersen -> replace one vertex by K3 (G12: c3=1, c5=6, c6=10, matching Exoo's counts)
  -> replace 11 of 12 vertices by H7; SA over (bare vertex) x (6^11 wirings).
  * {4,8,16}-free wirings found ONLY when the bare vertex is one of the three K3
    vertices (matches Exoo's Figure 5 solid vertex); non-triangle bares floor at
    c16 = 38-48. Found 3 wirings, forming 2 iso classes; both n=78 cubic connected,
    girth 3, c3=22, c4=c8=c16=0. This independently CONFIRMS f(4) <= 78.
  * NEW (not in Exoo's paper): both classes CONTAIN 32-cycles (found instantly), and
    class 1 also contains 64-cycles; so G78-type graphs die at the NEXT dyadic gate --
    an Erdos-Gyarfas counterexample at n=78 < 128 must also avoid {32, 64}. The
    {4,8,16}-free stratum at n=78 sits hard against the C32 wall.

REVISED INTERPRETATION (HYP-2359): the dyadic gaps close in a LADDER: opening {8} costs
c16 ~ 600-1000 (girth 5, n=28-40) or ~2000-3000 (Cayley girth 6); opening {8,16}
(possible only at n >= 54 by Markstrom's unpublished bound, realized at n = 78 by Exoo)
costs immediate C32 (verified on both reconstructed G78 classes). Every stratum we can
reach pays the next power of two. Erdos-Gyarfas unchallenged; the new quantitative signal
is the geometric growth of the next-gate cycle mass as each gate is closed.
"""

def main():
    out = "05-knowledge/results/dyadic_branch3_summary_kps2.out"
    with open(out, "w", encoding="utf-8") as f:
        f.write(SUMMARY.strip() + "\n")
    print(SUMMARY)
    print("saved ->", out)

if __name__ == "__main__":
    main()
