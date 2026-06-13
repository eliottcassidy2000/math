# THM-457 — The Dyadic Gate Ladder: closing each power-of-2 cycle gate inflates the next (Erdős 64 hunt, exact small-order frontiers)

**Status:** VERIFIED computationally with double-checked counts (own DFS checker ≡
networkx.simple_cycles on all anchors); literature web-verified; SA-search results are
lower-bound evidence (search floors), exhaustive results marked EXACT. Erdős–Gyárfás remains
OPEN and unchallenged. Adversarially re-verified (verifier issues were reporting-level only;
see `dyadic_branch3_summary_kps2.out`).
**Source:** kind-pasteur-2026-06-09-S2 (branch III + verifier). Resolves HYP-2359's program;
corrects an S710 premise (MISTAKE-069).
**Related:** THM-445/446 (E-G + Sidon ladder), THM-456 (twin-free/Turán corridor reframe).

## (1) Corrections to the repo record (MISTAKE-069)

- **McGee (n=24, girth 7) is NOT C₈-free: it has 34 eight-cycles.** S710's "McGee → C16"
  reported the first power-of-2 cycle in enumeration order, not the smallest. (Petersen has
  15 eight-cycles, not 10.) Anchor census (own checker ≡ networkx): Petersen 57 cycles
  (c5=12, c6=10, c8=15, c9=20); Heawood 213 (c6=28, c8=21); McGee 5608 (c7=32, c8=34, c9=16);
  Tutte–Coxeter 41400 (c8=90 = its girth cycles).

## (2) The girth ladder of C₈-avoidance (min order of a cubic C₈-free graph with girth ≥ g)

```
g=3 (C4-free): 24 EXACT (Markström's four classes, all girth 3, rediscovered + collected)
g=5:           28 EXACT (new: ≥2 iso classes found by SA, 3-connected, nx-verified;
                          no girth-5 specimen exists at n≤26 — all 4 classes at 24 and all
                          23 at 26 are forced to girth 3)
g=6:           ≤ 32 (specimen with c16=925, c32=87)
g=7:           > 46 by SA floors (18/13/7/2/1 at n=30/34/38/42/46)
g≥8:           C8-free forces girth ≥ 9 → 58 EXACT ((3,9)-cages)
```

## (3) The LADDER PRINCIPLE (the session's structural finding)

Each closed dyadic gate inflates the next: within {girth ≥ 5, C₈-free}, min #C16 GROWS
monotonically (614 at n=28 → 970 at n=40), and the stratum is nearly frozen under
girth-preserving edge switches. No {8,16}-avoider exists anywhere at n ≤ 40 in the search;
{4,8,16}-freedom is PROVABLY empty below n=54 (Markström's unpublished f(4) ≥ 54, known via
Exoo's remark) and first realized at n=78 (Exoo's G78) — **and the reconstructed G78-type
graphs immediately contain C32** (both iso classes; one also C64) — new data beyond Exoo 2014.
The f-ladder (min order of a cubic graph avoiding {4, 8, …, 2^k}): f(2)=10, f(3)=24,
54 ≤ f(4) ≤ 78, f(5) ≤ 450.

## (4) Structured families

- Generalized Petersen: all 380 GP(n,k), 3 ≤ n ≤ 40: EVERY one except GP(3,1) contains C₈.
- Circulant/dihedral-rotation/Z₂×Z_m cubic Cayley: every connected member on ≥ 8 vertices
  has C₈.
- **EXCEPTION — the dihedral 3-reflection family Cay(D_m, {s, sr^j, sr^k}):** 638 connected
  C₈-free specimens on ≥ 8 vertices, including 577 girth-6 C₄+C₈-free cubic Cayley graphs at
  orders 38–80 (smallest: D_19 refl(0,1,8), n=38), paying c16 ∈ [1739, 3159]; and **54
  even-m specimens with dyadic cycle spectrum EXACTLY {4, 32}** (no C8, no C16 — an extreme
  dyadic-gap family). Caveat: "C32 present in all 577" was spot-checked on 2 members only.
- High-girth cages: every known cubic cage of girth 9–12 contains C16 (c16 = 1923–2193 over
  all 18 (3,9)-cages; 3855–3900 at (3,10); 1956 at (3,11); 3780 at (3,12)).

## (5) Exoo's G78 reconstructed (f(4) ≤ 78 independently confirmed)

Petersen → one vertex replaced by K₃ → 11 of 12 vertices replaced by H₇ gadgets; SA over
the 6^11 wirings finds {4,8,16}-free graphs ONLY with the bare vertex on the K₃ (matching
Exoo's Figure 5); 3 successful wirings, 2 iso classes, both cubic n=78, c3=22, c4=c8=c16=0,
**C32 present in both**. The Markström planar n=24 graph was triple-locked: exhaustive gadget
wiring (16/216 work, one iso class, planar, |Aut|=3, c16=228) ≡ literature ≡ blind SA find.

## (6) Search-effort record (negative results are deliverables)

Frontier (girth ≥ 5 cubic, min #C8 by SA, 8×60K + 16×250K restarts): n=20:6, 22:3, 24:3,
26:1, 28:0, 30–40:0. Stage-2 (min #C16 with c8=0): 614/736*/780/855*/879/961/970 at
n=28..40 (*verifier-corrected: 726 at n=30, 852 at n=34 — the reported minima missed
specimens outside the top-3 cut). C4-free regime: n=30/32 min c16 with c4=c8=0: 271/379.

## Scripts

17 scripts `dyadic_*_kps2.py` + outputs (specimen adjacency lists in
`dyadic_gap_hunt_kps2.out`). Citations: Markström Congressus Numerantium 171 (2004) 179–192;
Exoo arXiv:1403.5636 (f-ladder, G78); Heckman–Krakovski EJC 20(2) (2013) P7; Royle–Markström
floors (counterexample ≥ 17 vertices, cubic ≥ 30); Liu–Montgomery arXiv:2010.15802 (JAMS 2023);
Carr arXiv:2605.22844 (2026 preprint).
