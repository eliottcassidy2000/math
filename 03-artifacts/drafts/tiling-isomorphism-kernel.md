# Tiling Isomorphism Classes and the SC+SF Symmetry Kernel

**Instance:** opus-2026-03-06-S7
**Status:** VERIFIED computationally n=3,...,7

---

## Summary

Tournament isomorphism classes group tilings of the staircase grid Grid(n)
into orbits under the symmetric group S_n acting on vertex labels. This
analysis identifies a persistent "symmetry kernel" — the set of classes
that are simultaneously **self-converse** (T ~ T^op) and **self-flip**
(flipping all tile bits keeps the class invariant).

---

## Master Data Table

| n | Classes | SC | SF | SC+SF | Max H | Max H in kernel? |
|---|---------|----|----|-------|-------|------------------|
| 3 | 2       | 2  | 0  | 0     | 3     | N/A              |
| 4 | 4       | 2  | 1  | 1     | 5     | YES              |
| 5 | 12      | 8  | 2  | 2     | 15    | NO               |
| 6 | 56      | 12 | 8  | 2     | 45    | YES              |
| 7 | 456     | 88 | 30 | 8     | 189   | NO               |

Class counts match OEIS A000568 (non-isomorphic tournaments on n vertices).

---

## The SC+SF Kernel

### Definition

A class is in the **symmetry kernel** if it satisfies both:
1. **Self-converse (SC):** The tournament T is isomorphic to its converse T^op
2. **Self-flip (SF):** Some tiling in the class, when all bits are flipped
   (every arc reversed relative to the base path), yields another tiling in
   the same class

### Kernel Classes

**n=4 (1 class):**
- #2: H=5, |Aut|=1, scores=(2,2,1,1), GS=3/5

**n=5 (2 classes):**
- #8:  H=11, |Aut|=1, scores=(3,2,2,2,1), GS=3/11
- #10: H=13, |Aut|=1, scores=(3,2,2,2,1), GS=3/13

**n=6 (2 classes):**
- #46: H=41, |Aut|=1, scores=(3,3,3,2,2,2), GS=9/41
- #54: H=45, |Aut|=3, scores=(3,3,3,2,2,2), GS=9/15  **= H-maximizer**

**n=7 (8 classes):**
- #129: H=125, scores=(4,4,3,3,3,2,2), GS=9/125
- #296: H=147, scores=(4,3,3,3,3,3,2), GS=9/147
- #311: H=117, scores=(4,4,4,3,2,2,2), GS=5/117
- #318: H=153, scores=(4,3,3,3,3,3,2), GS=9/153
- #336: H=133, scores=(4,4,3,3,3,2,2), GS=7/133
- #388: H=159, scores=(4,3,3,3,3,3,2), GS=9/159
- #389: H=151, scores=(4,3,3,3,3,3,2), GS=9/151
- #396: H=141, scores=(4,4,3,3,3,2,2), GS=7/141

All n=7 kernel classes have |Aut|=1, |AntiAut|=1, and nearly-regular scores.

### Kernel Growth: 0, 1, 2, 2, 8

---

## Key Structural Findings

### F1: H-maximizer alternates between kernel and non-kernel by parity

At **even n** (4, 6): the H-maximizer is in the kernel (SC+SF).
At **odd n** (5, 7): the H-maximizer is SC but NOT SF. The maximizer is
the regular tournament (all scores equal), which has high automorphism
group (|Aut|=5 at n=5, |Aut|=21 at n=7) — too symmetric to be self-flip.

This connects to the observation that at odd n, the maximizer is the
Paley tournament, whose high symmetry prevents it from being its own
"anti-tournament" relative to any fixed Hamiltonian path.

### F2: Grid-symmetric tilings per kernel class ~ 3^floor((n-2)/2)

| n | 3^floor((n-2)/2) | GS counts in kernel |
|---|------------------|---------------------|
| 4 | 3                | [3]                 |
| 5 | 3                | [3, 3]              |
| 6 | 9                | [9, 9]              |
| 7 | 9                | [9,9,5,9,7,9,9,7]  |

Exact for n=4,5,6. At n=7, holds for 5/8 kernel classes (those with
score (4,3,3,3,3,3,2)). The two exceptions with GS=7 have score
(4,4,3,3,3,2,2); the one with GS=5 has score (4,4,4,3,2,2,2).

### F3: GS/size ratio of 3/5 persists for kernel H-maximizers at even n

- n=4 #2:  GS/size = 3/5 = 0.600
- n=6 #54: GS/size = 9/15 = 3/5 = 0.600

Both are the H-maximizer within their kernel. The ratio 3/5 may persist
at n=8 — prediction: the n=8 kernel should contain an H-maximizer class
with GS/size = 27/45 = 3/5.

### F4: Parent-child inheritance chain

The SC+SF property follows a refinement rule across n:
- From a kernel (SC+SF) parent, children can be: SC+SF, SC-only, SF-only, or neither
- **New kernel members can emerge** from SC-only parents (e.g., n=6 #46 from n=5 #4)
- SF-only children from kernel parents always come in transpose pairs

Key chain:
```
n=4: #2 (H=5, SC+SF)
  -> n=5: #8 (H=11, SC+SF), #10 (H=13, SC+SF), #11 (H=15, SC only)
    -> n=6: #54 (H=45, SC+SF), #55 (H=45, SC only), #51 (H=43, SF only)
```

### F5: SF-only classes always come in transpose pairs

At n=6: 6 SF-only classes form 3 pairs {#14,#18}, {#29,#40}, {#51,#52}.
At n=7: 22 SF-only classes form 11 pairs by H value.

Each pair consists of T and T^op (same H, same |Aut|, but not isomorphic
to each other). The self-flip property is preserved by T -> T^op.

### F6: Regular tournaments at n=7

Three isomorphism classes of regular tournaments (scores all 3):
- Paley: H=189, |Aut|=21, GS=3/9
- Class B: H=175, |Aut|=7, GS=9/25
- Class C: H=171, |Aut|=3, GS=9/57

ALL are self-converse with anti-automorphisms. NONE are self-flip.
This confirms THM-028 (kind-pasteur): the three classes are distinguished
by their directed cycle counts (alpha_1 = 80, 65, 59 respectively).

---

## H Values of Kernel Classes vs Maximum

| n | Kernel H range | Max H | Kernel/Max ratio |
|---|---------------|-------|------------------|
| 4 | [5, 5]        | 5     | [1.000, 1.000]   |
| 5 | [11, 13]      | 15    | [0.733, 0.867]   |
| 6 | [41, 45]      | 45    | [0.911, 1.000]   |
| 7 | [117, 159]    | 189   | [0.619, 0.841]   |

The kernel H values are concentrated in the upper range of all H values,
but the ratio to max H decreases at odd n (where the maximizer is outside
the kernel).

---

## Scripts

- `04-computation/tiling_isomorphism_deep.py` — full analysis n=3,...,6
- `04-computation/tiling_iso_n7_c.py` — C-accelerated n=7 analysis
- `04-computation/symmetry_kernel_chain.py` — parent-child chain tracking

---

## Open Questions

1. Does the kernel count continue: 0, 1, 2, 2, 8, ?, ...? What is the n=8 count?
2. Is GS = 3^floor((n-2)/2) exact for all kernel classes at n >= 8?
3. Does the alternating kernel/non-kernel pattern for H-maximizer persist?
4. What is the structural explanation for why SF is incompatible with
   high automorphism groups (regular tournaments)?
5. Is there a formula for the number of SF-only classes at each n?
   Data: 0, 0, 0, 6, 22 for n=3,...,7.
