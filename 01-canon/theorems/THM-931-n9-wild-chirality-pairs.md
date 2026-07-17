# THM-931 — The n=9 invisible census: the first WILD pairs exist, and they are CHIRALITY pairs (T, T^op)

**Status:** VERIFIED (exhaustive machine census of all 191,536 n=9 iso classes; exact canonical forms; every wild pair re-verified with exact rational arithmetic)
**Instance:** klein-2026-07-16-S317
**Files:** `04-computation/n9_wild_hunt_klein_S317.py`, `04-computation/wild_eight_analysis_klein_S317.py`, results `n9_wild_hunt_klein_S317.out`, `wild_eight_analysis_klein_S317.out`
**Companions:** THM-918 (coning tower), THM-924 (walk reciprocity — cpK omitted as cpA-determined), THM-925 (n=8 completeness: invisibility = cone stratum exactly)

---

## Census

Method of THM-925 one rung up: 6880 × 2⁸ = 1,761,280 extension tournaments cover all n=9 classes; exact extended panels (cpA, cpL, H, sorted τ_in, sorted τ_out) via batched integer Faddeev–LeVerrier / DP / matrix-tree minors (spot-validated against exact rational recomputation); exact individualize-and-refine canonical forms inside panel buckets. **Validation: 191,536 classes = A000568(9) exactly.**

- **One-eyed panel (cpA, cpL, H, τ_in): 264 invisible groups, 277 pairs** (n=8 predictions: 265 cone pairs + 27 cocone pairs, overlapping).
- **Extended panel (+τ_out): 35 groups, 35 pairs**, splitting into:
  - **27 manufactured pairs**, all with sink AND source. The two predicted families — cones of the 27 τ_out-tied n=8 pairs and cocones of the 27 τ_in-tied n=8 pairs — **coincide**, because cone and cocone commute: cocone(cone(X)) = cone(cocone(X)) = the suspension Σ(X). The manufactured stratum at n=9 is exactly the suspension images of the n=7-rooted tower (the 4 double-blind cones included).
  - **8 WILD pairs** — no sink, no source (near-regular scores (3,3,4,4,4,4,4,5,5)-type), members from distinct n=8 base classes, each verified non-isomorphic by canonical form and panel-tied by exact arithmetic.

## The wild eight are chirality pairs

**All 8 wild pairs have T₂ ≅ T₁^op** (verified 8/8): each is a non-self-complementary tournament panel-tied with its own reversal. The mechanism, exactly:

- cpA(T^op) = cpA(T) and H(T^op) = H(T) — automatic;
- cpL(T^op) = cp of L_in(T); τ_in(T^op) = τ_out(T) — so the pair ties **iff cpL is self-dual (cpL_in = cpL_out as polynomials) and the tree vectors balance (τ_in = τ_out as multisets)**. Verified 8/8.

Consequences, all verified:
- The **cpA-deck ties automatically** (cards of T^op = op-cards; cpA is op-invariant per card) — the deck heuristic that separated every invisible pair at n ≤ 8 is *structurally blind* to chirality. Likewise c₃, c₄, c₅ (139–165, equal per pair), tr A^k, |Aut| (1 for seven pairs; 3 for one).
- The **full-panel deck** (each card's (cpA, cpL, H, τ_in, τ_out), as a multiset) **splits all 8** — per-card tree vectors detect the handedness.
- Every pair is op-closed as a set; no member is self-complementary; all groups have size exactly 2.

## Reading

Spectral blindness bifurcates into exactly two species at n=9:

1. **Laundered** (THM-918's mechanism): collapse the last witness with a sink/source; lives on the suspension tower; complete at n=8 (THM-925), 27 pairs at n=9.
2. **Chiral** (new): the panel — though not literally reversal-invariant — cannot certify orientation for tournaments whose L-spectrum is self-dual and whose tree vectors balance. First occurrences: 8 pairs at n=9, none at n ≤ 8.

The census question "what can the spectrum not see" has the sharp answer: **at n ≤ 9 it cannot see exactly (i) what a sink/source erases and (ii) handedness on the τ-balanced locus.** Open (shaped): is every wild pair at n=10 still a chirality pair, or does a genuinely non-chiral native tie appear? (The n=10 census needs ~9.7M extensions of 191,536 reps — feasible with the same machinery, ~2-3 h.)
