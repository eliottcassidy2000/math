# The shear spectrum: every T-continuation family carries a ladder of growth constants, and the parameter's seat decides everything

**kind-pasteur-2026-07-20-S128c103** (HYP-8170) · owner: shear the n·2^x+1 grid
down 1 (triangle) or 2 (Fibonacci-analogous); sum or product columns; find all
meaningful continuations of the triangular numbers and compare at shears 0,1,2,3+;
"the triangular numbers are key, they arise from relation itself, the edges of
complete graphs."

## 1. The machinery and the controls

A family is a grid G(c, n); the shear-s row sum is R_s(m) = Σ_c G(c, m−sc).
Pascal/simplex is the paradigm and validates the instrument exactly: s=0 gives
central binomials C(2m,m), s=1 gives 2^{m−1}, s=2 gives Fibonacci, s=3 gives
Narayana's cows — the **Pascal ladder** of dominant roots x^s = x^{s−1}+1:
2, φ, 1.4656, … ↓ 1.

## 2. The Proth theorems (exact, GF-verified)

For the owner's grid G(x, n) = n·2^x + 1:

    GF of R_s  =  t/((1−t)²(1−2t^s)) + t/((1−t)(1−t^s))     (verified s = 1,2,3)

- **The Proth shear spectrum is 2^{1/s}**: 2, √2, 2^{1/3}, … — the owner's
  "Fibonacci-analogous" shear-2 grows at **√2, the repo's hypotenuse/
  pseudo-doubling constant**. The analogy to Fibonacci is constructional; the
  constant is better: it is the triangle foundation's own.
- **Pascal dominates Proth at every shear s ≥ 2** (φ > √2, Narayana > 2^{1/3});
  they tie at s = 1 (both 2). Relations beat doubling under shearing.
- **Self-similarity**: the square-array row sums are again Proth numbers:
  R₀(m) = m·2^{m+1} + 1. The grid is closed under its own row-summation.
- Shear-1 row sums = 2(2^m − 1) — twice Mersenne. The Proth triangle's rows sum
  to Mersenne doubles while its rims are 2n+1 and 2^x+1: Mersenne inside,
  Fermat-form on the rim.

## 3. The Faulhaber shear closes the c102 loop

The Rosetta triangle IS the Faulhaber grid (Σ_{j≤n} j^c) sheared down 1 — so
its NE-diagonals are Faulhaber-shear-2, and the pure values continue 1, 2, 4,
7, 12, 21, 37, 68, **129, 254, 520, 1099** (OEIS-NEW). Against c102's
deviation-inclusive 130, 255: **the owner's three +1 deviations are exactly
what tunes the diagonal residual from (…, 34, 74, 165) to (…, 34, 75, 166) =
A045648** — the chiral-polyomino law and the deviation rule are one hypothesis,
discriminated by the owner's T(8,4) (218 vs 225/226/227). Likewise the c102 row
sums 106, 317 sit exactly 1, 2 above pure Faulhaber-shear-1 (105, 315) — one
per deviation. Bookkeeping closes perfectly.

## 4. The dichotomy (sharpened from "relational vs positional")

Measured across simplex, Proth, Faulhaber, polygonal, centered, pyramidal,
pronic: shear sums are **exponential iff the family parameter sits in an
EXPONENT** (binomial upper index, 2^x, j^c) and **polynomial iff it sits in a
COEFFICIENT** (gonality k in ((k−2)n²−(k−4)n)/2, centered k, pronic c). The
owner's relational reading survives refined: T_n = C(n,2) = |E(K_n)| lives in
the exponent-seat families (sub-relation counting), and those are exactly the
ones with Fibonacci-like ladders. The three exponential ladders found:
Pascal's (roots of x^s = x^{s−1}+1), Proth's (2^{1/s}), Faulhaber's
(super-exponential drift ~2.1–2.3, constant unnamed — open).

## 5. Products, parity, and the honest ledger

- Simplex shear-2 products = **A073617** (Pascal +45° slope products) —
  construction validated against a known object.
- Proth shear-2 products 2, 3, 12, 25, 210, 567, 10296, … : OEIS-NEW (210 = the
  primorial makes a cameo).
- Fibonacci-diagonal alternating sums cycle with period 6 (1,1,0,−1,−1,0) — the
  ζ₆ shadow; Proth alternating sums are wild (sign changes, no small period).
- **Correction (my own probe)**: the klein-T1532 sequence-D "MATCH: simplex s=2"
  in the .out is a FALSE POSITIVE — an 8-term window caught Fibonacci's prefix;
  D diverges at term 9 (33 vs 34). D remains its own object, matching no
  catalog family/shear at depth 14. The 8-term-window lesson: match windows
  must exceed the family's known divergence depth.
- OEIS-NEW harvest this session: Proth s=2 sums (2,3,7,10,18,25,41,56,…),
  pure Faulhaber diagonals (…,129,254,520), Proth s=2 products, Faulhaber s=3
  (1,2,3,5,8,12,18,… — Fibonacci prefix, Moser-flavored break at 12).
  Candidates for the owner's OEIS-submission pipeline alongside klein's T1532
  batch.
