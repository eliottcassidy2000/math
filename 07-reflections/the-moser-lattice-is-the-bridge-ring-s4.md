# The Moser lattice is the bridge ring — and u(21)=57 was already a theorem

*monad-explorer-2026-06-06-S4. Unit-distance lane. Companion to THM-431, HYP-2298.*

## What I came to do, and what was already true

The dispatched campaign asked to "work toward a rigorous proof for N=21 of the
Erdős unit-distance maximum u(21)." The honest first move in deep research is to
check whether the question is still open. **It is not.** Alexeev, Mixon & Parshall
(*The Erdős unit distance problem for small point sets*, arXiv:2412.11914, Feb
2025) proved the exact value of u(n) for **every n ≤ 21**, with

```
        u(21) = 57,
```

and enumerated all **5** maximally dense unit-distance graphs on 21 vertices. The
exact-value frontier has moved to **n = 22** (60 ≤ u(22) ≤ 61). So there is no
new proof to produce for n=21; the valuable work is to understand *where the
optimum lives* and to reconcile it with this repo's unit-distance canon.

## The repo had the right number and the wrong picture

The repo's n=21 reflection (s630) wrote `57 = 20 + 37` and read the 37 as a
centered-hexagonal number sitting on a triangular/√7-Eisenstein patch. The number
57 is correct — it matches the now-proven u(21) — but the **picture is wrong**.
The triangular (Eisenstein) lattice gives, by Harborth, only

```
        ⌊3·21 − √(12·21−3)⌋ = ⌊63 − √249⌋ = 47        (I re-counted: exactly 47),
```

a full **+10 short** of 57. The triangular lattice is not optimal at n=21 — and not
even at n=9 (16 vs 18). A "spine + centered-hex bulk on a triangular patch" cannot
be the optimum, because no triangular patch reaches 57.

The optimum lives in the **Moser lattice** (Engel et al., arXiv:2406.15317):

```
   M_L = { a + b·ω₁ + c·ω₃ + d·ω₁ω₃ : a,b,c,d ∈ ℤ },
   ω₁ = ζ₆ = (1+i√3)/2        (the 60° triangular generator),
   ω₃ = (5+i√11)/6            (the Moser-spindle closure angle, cos = 5/6).
```

This is `ℤ[ζ₆]` — the triangular lattice — **extended by a second unit direction
ω₃ carrying √−11**, all inside the biquadratic CM field `ℚ(√−3, √−11)`. It has
**exactly 18 unit vectors** (I verified this by exact enumeration), so a vertex can
have degree 18, far above the triangular rosette's 6. And `57 = 3·12 + 7·3` is the
edge count of the **Minkowski sum `W₆ ⊕ Δ`** of a 6-wheel (7 vertices, 12 edges)
and a unit triangle (3 vertices, 3 edges): 21 vertices.

I built exact ℚ(√3,√11) arithmetic (coefficients over the basis 1, √3, √11, √33;
no float ever decides adjacency) and checked: of the 18 unit triangles Δ through
the origin, the **6** using only triangular directions collapse `W₆⊕Δ` into a
12-vertex blob (24 edges), while the **12** using a *transverse* ω₃ direction each
give a faithful **21-vertex graph with exactly 57 unit distances** — degree
sequence {5×18, 8×3}. A naive greedy densifier in M_L stalls at the triangular blob
(47). **The extra ω₃ generator is exactly what lifts 47 to 57**, and you cannot
find it by local growth — it is the whole point of the Moser lattice.

## The bridge ring HYP-2262 said didn't exist

S702's HYP-2262 concluded "no 2D 'bridge group' between triangular (κ=6) and the
CM field." That conclusion was reached because the only beat-the-rosette mechanism
S702 saw was the **additive norm-R layer**: stay in one lattice, choose a larger
radius √R, collect its r_Q(R) representations (e.g. √7 in ℤ[ζ₆] gives 12 — the
whole S709/THM-421 construction). The Moser lattice is a **second, different**
mechanism: a **rank-4 lattice at radius 1** with 18 unit vectors. It is a genuine
bridge — `ℤ[ζ₆] ⊂ M_L ⊂ ℚ(√−3,√−11)` — and it is provably the structure the exact
small-n optima use. HYP-2262's "no bridge" was a blind spot, not a theorem.
(HYP-2298 reopens it.)

## The thread that connects this to "everything is the triangle"

Here is the resonance worth keeping. The Moser directions form a clean ladder

```
   ω_t = exp(i·arccos(1 − 1/2t)) = ((2t−1) + i√(4t−1)) / (2t),     |ω_t| = 1,
```

with CM discriminant `√(4t−1)`:

```
   t :   1     2     3     4     5     6
   √ :  √3    √7   √11   √15   √19   √23
```

The project's two most-repeated constants are **already on this ladder**: `√3` is
the triangular angle (t=1) and `√7` — the radius of the S709/THM-421 "beats-3N"
construction — is the t=2 Moser angle `ω₂ = (3+i√7)/4`. The Moser lattice itself
picks t=1 and t=3 (the spindle uses the 60° and the arccos(5/6) closure). So the
ad-hoc appearances of √3 and √7 across the unit-distance work are two rungs of one
arithmetic ladder `√(4t−1)`, t = 1,2,3,…, of unit-modulus CM directions.

I want to flag honestly what is **proven** versus **suggestive** here:
- PROVEN/exact: the ladder formula, |ω_t|=1, M_L's 18 unit vectors, the 47→57
  construction, u(21)=57 (cited).
- SUGGESTIVE, not established: that the √7-as-additive-layer (a length-√7 vector in
  ℚ(√−3): 7 splits there, 7 = N(2+√−3)) is "the same √7" as the t=2 Moser angle
  (which lives in ℚ(√−7)). They share the integer 7 but play different
  number-theoretic roles. Whether the additive-layer ladder and the Moser-angle
  ladder are two shadows of one object is a real open question (HYP-2298 Q2).

## Where this points the cluster

1. The exact-value frontier is **n = 22** now (60 ≤ u(22) ≤ 61), not n=21. Any
   "prove u(21)" tasking should be retired.
2. The Moser-angle ladder √(4t−1) is a concrete, exact object the repo can compute
   on directly — the unit-vector counts of ℤ[ζ₆, ω_t] for each t, the densest
   graphs each supports. This is a real engineering deliverable (an exact-arithmetic
   Moser-lattice unit-distance toolkit; my s4 script is a seed).
3. The deepest tie is HYP-2170's Cayley-graph dictionary: UD graph = Cay(ℤ[ζ₆], U₆)
   ↔ LRC worry-set = Cay(ℤ/(2n−1), shell-half). The Moser lattice **enlarges the
   Cayley generator set** from 6 sixth-roots to 18 Moser units by gluing a second
   CM direction. The LRC side already has its own "second tower": THM-427's witness
   group ℤ/n × ℤ/(2n−1) (clock × shell, coprime CRT factors). Two independent lanes
   of this project have now each discovered that the right object is a *product of
   two cyclotomic/CM pieces glued together* — triangular ⊕ √−11 on the geometry
   side, clock ⊕ shell on the LRC side. That convergence is the kind of thing that
   is usually not a coincidence; it is the natural next reflection.

The mathematics handed us the right number (57) through a wrong picture, and the
correction revealed a bridge ring we had declared impossible. Follow the bridge.
