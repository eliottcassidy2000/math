# The second moments: the LRC atom pair-correlation is squared-Ramanujan, and the SC-spine's blue spectrum is metallic (golden/silver), bipartite at odd n

*opus-2026-07-01-S24. The owner asked to work the LRC atom pair-correlation (Ramanujan sums = the atoms' 2nd
moment) and the SC-spine's blue-degree vs blue-graph spectrum. Both are the **second-moment** layer above the
Lefschetz/trace (1st-moment) work of S23 — and both carry the same odd/even parity, with strikingly different
arithmetic characters.*

## Task 1 — the atoms' 2nd moment is |Ramanujan|² (the pair-correlation)
The AP lonely atoms are the units `(ℤ/N)*` (for N=14: the 6 primitive residues {1,3,5,9,11,13} = φ(14)). Their
**1st moment** is the trigonometric moment `μ̂(k) = c_N(k)` — the **Ramanujan sum** (klein HYP-3793), which S23
identified as a Lefschetz/Frobenius **trace**. Their **2nd moment** is the pair-correlation
`A(m) = #{unit pairs with difference m}`, and by Wiener–Khinchin its power spectrum is **exactly `|c_N(k)|²`**
(verified N=6,10,14). So:
```
    atoms' 1st moment  =  c_N(k)          (Ramanujan sum = Lefschetz trace)
    atoms' 2nd moment  =  |c_N(k)|²       (pair-correlation power spectrum = SQUARED Ramanujan)
```
For N=14 the pair-correlation is `A(0)=6=φ`, `A(m)=5` on **even** m, `A(m)=0` on **odd** m — because every unit
is odd, so all differences are even. The power spectrum `|c₁₄(k)|²` is a flat background `1` with **spikes of 36
at k=0 and k=7=N/2** — the `k=N/2` spike is the **antipode/ι frequency** (the complement involution). So the 2nd
moment sees exactly the ι-symmetry: DC (k=0) and Nyquist (k=N/2) carry the mass, the rest is flat — the
signature of an equidistributed, ι-symmetric atom cloud (the far-element/three-distance regime).

## Task 2 — the SC spine: T-join degrees vs a metallic, bipartite spectrum
The blue subgraph (= the SC spine = the half-tiling metagraph, HYP-3810) on the SC nodes, with the blue
flip-lines as edges:
- **Degrees = the T-join** (all **odd**, T = SC nodes) with **even multiplicities** (twin-pairing, mac-mini/kps):
  `[1,1]`, `[1,1,1,1,3,3,3,3]`, `[1,1,3,3,5,5,7×6]`, `[1×10,3×12,5×16,7×32,9×18]` for n=4,5,6,7.
- **Spectrum is bipartite (symmetric about 0) for ODD n** and self-looped for EVEN n — because **blue self-loops
  are even-n-only** (n=4:1, n=5:0, n=6:2, n=7:0). At odd n the blue graph is a genuine bipartite (SC↔SC across
  the T-join) with `+λ ↔ −λ` symmetry; at even n the self-loops break it.
- **The eigenvalues are quadratic irrationals — metallic ratios.** n=4: the **golden ratio** `(1±√5)/2` (the
  matrix `[[1,1],[1,0]]` = one edge + one self-loop). n=5: the **silver ratio** family `±{1, √2−1, 1+√2}`. So
  the SC-spine adjacency spectrum is **algebraic of degree 2 (metallic)**, *not* roots of unity.
- It is a **multigraph**: `Σλ² = Σ B_ij²` exceeds `Σ deg` at n=7 (552 > 512), so some SC pairs carry multiple
  blue lines. The 2nd spectral moment `Σλ² = trace(B²)` is the graph's own **pair-correlation** (closed
  2-walks) = degrees + twice the multi-edge excess.

## The parallel — two 2nd moments, one parity, two arithmetics
The two computations are the same object on the two pillars:

| | LRC atoms (circle) | SC spine (blue graph) |
|---|---|---|
| 1st moment | `c_N(k)` = Ramanujan trace | `Tr(B)` = #self-loops (0 at odd n) |
| **2nd moment** | `\|c_N(k)\|²` = pair-correlation | `Σλ² = Tr(B²)` = degree/2-walk pair-correlation |
| parity | units all **odd** ⇒ even differences; ι-spikes at k=0, N/2 | self-loops **even-n only** ⇒ **bipartite at odd n** |
| arithmetic | roots of unity / Ramanujan / **√p Gauss** (S22/23) | **metallic** (golden, silver — degree-2 irrationals) |

Two clean facts fall out:
1. **The odd/even parity is the same on both sides.** On the circle it is "the atoms are odd residues, so the
   pair-correlation lives on the even lattice with an ι (k=N/2) spike." On the spine it is "blue self-loops are
   even-n only, so the spine is bipartite exactly at odd n." Both are the complement/ι-symmetry (S22: is-1-in-
   the-spectrum) read at the 2nd-moment level.
2. **The arithmetic characters differ and are informative.** The *Cayley* spectra are cyclotomic/quadratic-Gauss
   (roots of unity, √p — S22/S23); the *adjacency* spectrum of the SC spine is **metallic** (golden/silver,
   degree-2 units). The spine is a recursive fold (half-tiling), and metallic ratios are the fixed points of
   `x ↦ k + 1/x` — the continued-fraction/self-similar signature — which is exactly the flavor of the LRC deep
   well `t* = [0; n−1, n]` (a metallic-type continued fraction). So the SC spine's spectrum and the LRC covering
   modulus share the **continued-fraction (metallic) arithmetic**, while the Paley obstruction shares the
   **Gauss-sum (√p) arithmetic** — the two hard-side flavors (S23) show up as two spectral flavors here.

## Status
- **Verified (Task 1):** atoms' pair-correlation `A`, `Â(k)=|c_N(k)|²` (Wiener–Khinchin, N=6,10,14); the ι-spike
  at k=N/2; the odd-residue ⇒ even-difference structure.
- **Verified (Task 2):** blue-degree = T-join (odd, even-multiplicity); blue spectrum bipartite at odd n /
  self-looped at even n; golden ratio (n=4), silver family (n=5); multigraph `Σλ²>Σdeg` (n=7).
- **Synthesis:** both are the 2nd-moment (pair-correlation) layer above the S23 traces; the complement/ι parity
  is shared; the arithmetic splits into **metallic (continued-fraction / spine / covering-modulus)** vs
  **Gauss-sum (√p / Paley)** — the two hard-side flavors.
- **Honest:** the pair-correlation = |Ramanujan|² is exact (Wiener–Khinchin); the metallic-spectrum reading is
  verified at n=4,5 (clean golden/silver) and qualitative at n=6,7 (radius 6.49, 7.35 — degree-2 but not a named
  metal); the "spine ≈ covering-modulus continued-fraction arithmetic" is a structural resonance, not a proof.

Related: HYP-3815 (Lefschetz/trace = the 1st moment; this is the 2nd), HYP-3793/klein (moments = Ramanujan
sums), HYP-3810 (blue subgraph = half-tiling = the SC spine), HYP-3811 (odd/even parity), S22 (Cayley roots of
unity), the eigenvalues-of-the-merged-metagraph reflection. HYP-3816 (this). Scripts:
04-computation/{lrc_atom_paircorrelation_ramanujan, sc_spine_bluedegree_vs_spectrum}_opus_20260701.py.
