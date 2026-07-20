# The shift-machine: triangular continuations, shifts, and the relational thesis

**Instance:** klein-2026-07-20-S331 (owner: the Proth grid as shifted triangles; find
all meaningful continuations of the triangular numbers; compare shifts 0/1/2/3 with
sums/products/metrics; "the triangular numbers arise from relation itself"). Frozen
survey: `shift_machine_survey_klein_S331.out`. Fleet context (pulled): death-star
THM-1355 (the Proth n·2^x+1 unification: observer 2n+1 × hypotenuse 2^x+1 axes;
Fermat/observer/Cullen slices), kps HYP-8165 (the owner's triangle is OEIS-new in
all readings — the "Rosetta triangle").

## 1. The machine

For a one-parameter family S_j(n) (j = family parameter) and shift d, form the
triangle whose row m collects S_{j₀+t}(m − d·t); take row sums, alternating sums,
products. Calibration on Pascal (simplex family): shift-1 sums = 2^n − 1 ✓,
alternating sums ≡ 1 ✓, and **shift-d sums = partial sums of the d-bonacci family
a(n) = a(n−1) + a(n−d)** (Fibonacci at d = 2, Narayana's cows at d = 3, verified by
exact linear-recurrence identification) — the owner's "shift-2 makes it Fibonacci-
analogous" is the general law.

## 2. The two clean discoveries

- **The Proth shift-1 law:** row sums of the shifted n·2^x+1 grid are EXACTLY
  2^{m+1} − 2 (doubled Mersenne) — closed-form: Σ_{x<m}(m−x)2^x + m = 2^{m+1} − 2.
  Complements THM-1355's slice picture: the antidiagonal mass of the Proth family is
  Mersenne. (Shift-2/3 Proth sums 2,3,7,10,18,25,41,56,88,119 / 2,3,4,8,11,14,22,…
  are engine-unidentified — filed as leads; likely interleaved Jacobsthal-flavored.)
- **The owner's triangle deconstructed:** the PURE power-sum shift-1 triangle
  (family PS_p, shift 1) reproduces the owner's triangle EXCEPT at exactly the two
  break entries — pure row 7 = {7,21,55,100,98,33,1} vs owner {7,21,55,101,99,33,1}.
  So the Rosetta triangle = shift-machine output + the two-point Moser break, and
  the PURE version's row sums are 1, 3, 7, 16, 39, 105, 315, 1048, 3829 (note
  315 = 3·105 — a clean ×3 step; sequence itself unidentified order ≤ 4, lead).

## 3. The catalog of triangular continuations (compared)

Surveyed as families through the machine: **polygonal** P_k(n) (2D figurate; k-axis),
**simplex** C(n+d−1,d) (dimension axis — Pascal's columns; the RELATIONAL reading:
d = arity of the relation, T = C(n,2) = binary relations = edges), **pyramidal/
polyhedral** (3D figurate), **centered polygonal** (k·T(n−1)+1 — note: centered
polygonals are LITERALLY k copies of the triangle around a center — the triangle as
universal part), **power sums** PS_p (the Faulhaber axis, S330), **Proth** n·2^x+1
(the arithmetic axis, THM-1355). Sums/alternating/products tabulated for shifts
1/2/3 in the frozen out; identified hits above; ~14 unidentified sequences filed as
a lead-bank (products are the wildest — factorial-Proth hybrids with 2-adic
structure worth mining).

## 4. The relational thesis, grounded

"T_n is n in the relational sense": T_{n−1} = C(n,2) = edges of K_n = the number of
PAIRS — the repo's whole architecture agrees: the staircase has C(n−1,2) tiles (the
tiling coordinate), tournament score sums are C(n,2) (every game distributes one
point), the odd sector's Faulhaber coordinate is T (S330), and the JC master cubic's
trace law lives on the pair-invariant w = 1 + xy. Every "continuation of the
triangular numbers" in the catalog is a different completion of the same primitive —
relation — along a different axis: arity (simplex), gon (polygonal), dimension
(pyramidal), centering (centered), moment (power sums), arithmetic (Proth).

## 5. Leads filed

(i) identify the shift-2/3 Proth sums and the pure power-sum row sums (OEIS pass);
(ii) the product rows' 2-adic valuations (Proth products = 2, 9, 100, 2835, 202878 —
v₂ pattern); (iii) polygonal shift-1 sums 1,4,11,25,50,91 (degree-4 polynomial fit —
extend engine order); (iv) the centered family's near-miss with the Rosetta columns;
(v) cross with kps HYP-8165's three readings and death-star THM-1355's slices.
