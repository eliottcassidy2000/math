# THM-925 — The complete n=8 invisible census: invisibility at n=8 is EXACTLY the cone stratum

**Status:** VERIFIED (exhaustive machine census of all 6880 iso classes; exact integer panels; canonical forms, not hashes)
**Instance:** klein-2026-07-16-S316
**Files:** `04-computation/n8_invisible_census_klein_S316.py`, `05-knowledge/results/n8_invisible_census_klein_S316.out`
**Companions:** THM-918 (coning tower, the construction), THM-924 (walk reciprocity, why cpK is omitted).

---

## Method

Every n=8 tournament contains an n=7 subtournament, so extending the 456 n=7 class representatives by all 2⁷ arc patterns of a new vertex (58,368 tournaments) covers all iso classes. For each, the **extended panel** was computed exactly: (cpA, cpL, H, sorted τ_in, sorted τ_out) — cpK omitted because THM-924 proves it is a function of cpA. Panels via batched integer Faddeev–LeVerrier (divisibility asserted), DP path counts, and batched matrix-tree minors (float dets validated against exact Bareiss on a sample). Panel buckets were then deduplicated by **exact canonical forms** (refined-partition minimization — an isomorphism decision, not an invariant hash). Validation: total classes = **6880 = A000568(8)** on the nose.

## Results

1. **One-eyed panel (cpA, cpL, H, τ_in): exactly 27 invisible pairs** — precisely the 27 cones THM-918 manufactured from the equal-H deep-tied n=7 pairs. The manufactured list is COMPLETE: no other pair of n=8 classes ties on the census panel.
2. **Two-eyed panel (adding τ_out): exactly 4 invisible pairs** — precisely the double-blind four (H = 23, 29, 31, 43). **Zero non-cone invisible pairs exist at n=8.** The laundering principle (collapse the last witness) is not just *a* source of spectral invisibility — at n=8 it is the *only* source.
3. Every invisible class contains **both a sink and a source** (score sequences (0,2,3,3,4,4,5,7) and (0,3,3,3,4,4,4,7)) — the bases carried a source (score 6 at n=7), so the cones are suspension-type. On all four pairs **τ_in = τ_out as vectors** (both collapsed to (0,…,0,D) with the *same* D ∈ {7560, 7616, 7644, 8568}) — the two tree collapses agree numerically, consistent with the L5 reversal-closure of the panel.
4. **Spectral-tie growth:** cpA-tie groups at n=8: 1460, covering 6817/6880 classes (**99.1%**, up from 89% at n=7). Distinct cpA: 1523.
5. **The symmetrization-collapse sequence** (# distinct skew charpolys, by THM-924 = # distinct symmetrized cpA): **1, 2, 2, 6, 11, 50** for n = 3..8 — vs 2, 4, 12, 56, 456, 6880 classes and (at n=8) 1523 distinct cpA. New sequence, OEIS check queued.

## Reading

Combined with THM-918: the invisible stratum at n=8 is **exactly** the image of the cone construction — spectral blindness is manufactured, never accidental, at this order. The first "wild" (non-cone) invisible pair, if it exists, lives at n ≥ 9. Conjecture-shaped follow-up: does the cone stratum stay complete at n=9 (where cones of the 4 + cones of new n=8 equal-H deep ties + possibly wild pairs compete), or does the sea produce its first native tie?
