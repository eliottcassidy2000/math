# THM-453 (renumbered from THM-449; that number was concurrently reserved by mac-mini-2026-06-09-S1) — H(T[K₂]): strong-component product law, twin-lift cycle laws, and the failure of odd-cycle determination

**Status:** Product law PROVED (sketch below) + VERIFIED 74/74 iso classes n=3..6. Lift laws
VERIFIED 74/74 (c7' coefficient 128 analytic-only — c7=0 at n≤6, flagged). IP-determination
REFUTED with exact counterexample. Adversarially re-verified with three independent Hamiltonian
counters including a full 10! permutation scan (`verify_A_H_formula_kps1.out`).
**Source:** kind-pasteur-2026-06-09-S1 (branch A + verifier). Directly answers OPEN-Q-045 Q1.
**Related:** THM-447, OPEN-Q-045, OCF (THM-002), HYP-2334 (resolved into this), HYP-2347.

## (1) REFUTED: I(Ω(T), x) does NOT determine H(T[K₂])

Counterexample (n=5): idx4 (bits 0001000000, scores (1,1,2,3,3)) vs idx6 (bits 0101000000,
scores (2,1,1,3,3)). Both have IP = 1 + 4x AND the same length-typed IP (three 3-cycles, one
5-cycle), yet H(T[K₂]) = 3225 vs 2785. **Cause: the odd-cycle conflict graph omits EVEN cycles
(c4 = 3 vs 2), and even cycles of T become odd cycles of T[K₂] via twin insertion.** The
doubling is a parity mixer: OCF data of T is insufficient for OCF data of T[K₂].

## (2) PROVED: the strong-component product law

```
H(T[K₂]) = ∏_{C nontrivial strong component of T} H(C[K₂])
```

Proof sketch: the strong components of T[K₂] are exactly C[K₂] for |C| ≥ 2 plus the twin pairs
of vertices in singleton components; the condensation is transitive; Hamiltonian paths of a
tournament factor over the transitive chain of strong components (Moon). Verified 74/74
(e.g. n=6 idx26 = disjoint union of two 3-cycles: H = 45² = 2025). **This — not the IP — explains
the cross-n collisions** H(T[K₂]) ∈ {1, 45, 393, …}: every tournament whose only nontrivial
strong component is a single 3-cycle gets H(T[K₂]) = 45 at EVERY n.
NOTE: the product law FAILS for D_skew and SC-blowup (they hold only trivially in the 43
all-strong classes; both doublings strongly connect across components).

## (3) VERIFIED 74/74: twin-lift laws for the odd cycles of T[K₂]

Odd cycles of T[K₂] = twin-decorated lifts of closed walks of T with vertex multiplicity ≤ 2:
```
c3' = 8·c3
c5' = 32·c5 + 32·c4 + 6·c3
c7' = 128·c7 + 192·c6 + 80·c5 + 8·c4 + 64·p331 + 48·p332 + 64·p341 + 32·p342
```
where p(a,b,s) = number of cycle pairs of lengths (a,b) sharing s vertices. Single-cycle
coefficient: a k-cycle with t twin insertions contributes C(k,t)·2^(k−t) (the 128 = C(7,0)·2^7
term is untestable at n ≤ 6 where c7 = 0 — analytic extrapolation, flagged for an n=7 test).

**New combinatorial identity discovered via rank deficiency of the feature matrix (74/74):**
```
p(3,4,3) = 2·p(3,3,2)
```
(each triangle-pair sharing an arc yields exactly two (triangle, 4-cycle) all-3-shared
incidences — mechanism asserted, not yet bijectively proved).

## (4) VERIFIED 74/74: congruences

```
H(T[K₂]) ≡ 2·H(T) − 1  (mod 8)     [NOT mod 16: census {0:39, 8:35}]
H(SCblow(T)) ≡ 1 (mod 4)
```
A Rédei-parity refinement transported through the blowup.

## (5) VERIFIED (n ≤ 6): the full cycle spectrum (c3,c4,c5,c6) determines H(T[K₂]) — and more

Across all 74 iso classes n=3..6 (32 spectrum groups, 22 with ≥ 2 members, 0 broken): equal
cycle spectrum ⟹ equal H(T[K₂]), and in fact equal ENTIRE independence polynomial of Ω(T[K₂])
(13/13 same-spectrum strong pairs at n=6 have identical c3'..c11' and (i1..i4) at 12 vertices).
CAVEAT (verifier): at n ≤ 6 spectrum-equal classes also share pair-intersection features, so
"spectrum alone" vs "spectrum + pair statistics" is NOT yet separated — n=7 stress test needed
(HYP-2347). No global polynomial functional exists (disjoint-union multiplicativity forbids it).

## (6) Contrast: H(D_skew(T)) is NOT determined by any cycle statistic

Counterexample (n=4): idx1 vs idx3 (C3 + dominated vs dominating apex): identical cycle
spectrum (1,0,0,0), identical pair features, H(D) = 189 vs 333. D is op-ASYMMETRIC:
H(D(T)) ≠ H(D(T^op)) in 50/74 classes — the doubling's chirality (cf. THM-451). Odd-n locking:
H(D(T)) ≡ 1 (mod 4) at n odd (≡ 1 mod 8 at n=5, 12/12; ≡ 5 mod 8 at n=3).

## (7) New integer sequences (NOT in OEIS, live-checked)

```
H(D(transitive_n))      n=3..7:  13, 95, 1033, 15611, 313285
H(SCblow(transitive_n)) n=3..7:  41, 629, 14937, 513669, 24104937
```

## Scripts

`ip_doubling_hunt_kps1.py`, `k2_spectrum_functional_kps1.py`, `c7_lift_law_verify_kps1.py`,
`verify_A_H_formula_kps1.py` (+ .out in 05-knowledge/results/).
