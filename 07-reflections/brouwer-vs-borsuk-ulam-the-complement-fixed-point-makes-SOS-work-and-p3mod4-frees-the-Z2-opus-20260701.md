# Brouwer vs Borsuk–Ulam: the complement fixed point makes SOS work, and p≡3 mod 4 frees the ℤ₂ (so the three pillars are the hard-regime tools)

*opus-2026-07-01-S22b. The owner proposed a topological organizing principle: the complement reflection is an
automorphism, so a Brouwer fixed point exists and SOS works; the hard p≡3 mod 4 regime is where the ℤ₂ becomes
free (Borsuk–Ulam), which is exactly what the three pillars are for. I pushed the open extensions and verified
the pivot. It holds — and it explains *which* LRC/tournament problems are SOS-tractable and which are not.*

## The extensions, pushed
- **c₃ is a spectral moment (dictionary third face, verified n≤5).** The number of cyclic 3-cycles equals
  `c₃ = Tr(A³)/3`. So the same quantity is: algebra (`C(n,3) − ΣC(sᵢ,2)`, Kendall), geometry (**grid triangles**
  with cyclic orientation on the staircase), and **spectrum** (`Tr(A³)/3`, a moment of the tournament
  adjacency). Three faces, as with complement.
- **The "i" is intrinsic to the Cayley map, not to ρ=R∘V.** A tournament's skew matrix `S` has *imaginary*
  eigenvalues `iλ` — the order-4 element `i` (`i²=−1`) is already inside the Cayley transform (`S` skew ⇒
  spectrum imaginary ⇒ the Cayley image on the unit circle). The geometric `ρ=R∘V` (a 90° tile-grid rotation,
  `ρ⁴=id`) is **not** a "√complement" at the algebra level, because `R=complement` is a *reflection* and no plane
  isometry squares to a reflection; and `V` is iso-blind (S21). So the honest home of the Cayley–Dickson `i` is
  the skew-imaginary spectrum, not `ρ`.

## The pivot: complement is an involution — Brouwer where it has a fixed point, Borsuk–Ulam where it's free
The complement (`Tᵒᵖ`; antipode `ι:t↦1−t`; conjugation `z↦z̄` on the Cayley circle; negation `x↦−x`) is a **ℤ₂
involution** — a continuous automorphism of order 2 of the configuration space. Everything turns on whether it
acts with a fixed point or freely.

### Brouwer side (fixed point ⇒ SOS works)
On the convex, compact body of moment sequences (the spectrahedron of the Lasserre/SOS relaxation), the
involution `ι` is an affine self-map, so it has a **fixed point** — the `ι`-invariant configuration. When the
*optimizer* is such a fixed point, the **`ι`-invariant (symmetry-reduced) SOS relaxation is exact at the same
degree** (Gatermann–Parrilo symmetry-adapted SOS). Concretely:
- The AP lonely measure is `ι`-invariant (`G(t)=G(1−t)`, klein) — a fixed point — so the symmetric SOS captures
  it, and the `M(S)≥1/n` bound is SOS-certifiable there. **"SOS works."**
- `SC` tournaments are the complement-fixed points on the staircase (grid-symmetric tilings) and on the circle
  (conjugation-symmetric spectra) — the Brouwer fixed set is exactly the self-complementary world.

So the chain the owner gave — *complement is an automorphism → Brouwer fixed point exists → SOS works* — is the
**easy regime**: a symmetric optimum, certified by the invariant SOS.

### Borsuk–Ulam side (free ℤ₂, p≡3 mod 4 ⇒ the hard regime)
The pivot is arithmetic and I verified it exactly: **`p ≡ 3 (mod 4) ⟺ −1 is a quadratic non-residue ⟺ negation
`x↦−x` swaps QR↔QNR *freely*** (verified p=3,7,11,19,23 free; p=5,13,17 fixed). This is precisely where:
- the **Paley tournament exists** (antisymmetry needs `−1 ∈ QNR`; at `p≡1` you get a *graph*, not a tournament),
- the Paley tournament is the **flip-rank obstruction** (HYP-3805, the |Aut|=21 heptagon at p=7) and the
  covering-min hardness lives (klein HYP-3812, metric-irreducible at Φ₆=3·61),
- the complement is realized by negation, a **free ℤ₂** on the QR-structure (one residual fixed vertex `0` = the
  single real-axis Cayley eigenvalue of the odd-`n` anti-automorphism, S22/HYP-3811).

A free ℤ₂ action is the domain of **Borsuk–Ulam**: there is *no* ℤ₂-equivariant map to a lower-dimensional
sphere. Translated to certification: **there is no symmetric SOS certificate for a free-orbit optimum** — the
invariant relaxation cannot "choose" one of the two complement-swapped extremals continuously. So when the
extremal is a free ℤ₂ orbit (Paley / QR / p≡3 mod 4), the SOS/moment machinery does *not* close by symmetry
reduction — exactly the obstruction my MOSEK campaign hit (every naive `p₀/L_y` relaxation returns the trivial
bound; the coupling is the singular series, HYP-3791) and exactly the metric-irreducibility klein found.

## Why the three pillars are the free-regime tools
`⟨fold-thinking diagnoses obstruction⟩` (klein HYP-3812) plus this: the **three pillars** (mac-mini's
Kaczmarz+Christoffel+Blaschke merge, HYP-3796) are precisely the *constructive* methods that survive a free ℤ₂,
where symmetric SOS cannot:
- **Kaczmarz / POCS (alternating projections):** a *constructive* witness search that converges to the extremal
  intersection **without** needing a fixed point — the standard replacement for a Brouwer certificate when the
  action is free.
- **Christoffel / CD kernel (flat extension):** the atomic (Curto–Fialkow) structure — in the free regime the
  atoms come in **complement-pairs** with no central atom (the 2-atom `{t*, 1−t*}` covering-min), i.e. a free
  orbit, handled atom-by-atom rather than by a symmetric moment matrix.
- **Blaschke dynamics (degree-1 circle map):** the runner dynamics — a free ℤ₂ ⇔ **no fixed point** of the
  circle map ⇔ Diophantine/Herman rigidity (the deep-well isolation `t*=[0;n−1,n]`, mac-mini HYP-3796). The
  Borsuk–Ulam obstruction *is* the absence of a fixed rotation number.

So the owner's chain completes: **Brouwer (fixed point) → symmetric SOS → the easy/AP/p≡1 regime; free ℤ₂
(p≡3 mod 4, Paley) → Borsuk–Ulam → no symmetric certificate → the three pillars (POCS / flat-extension /
Blaschke) are the required constructive tools.** The topology sorts the problem by whether the complement fixes
the extremal or moves it freely, and the arithmetic switch is `−1 ∈ QR?` = `p mod 4`.

## The unified dictionary (now four columns)
| Algebra | Staircase | Circle (Cayley) | Topology / method |
|---|---|---|---|
| complement `Tᵒᵖ` | anti-diagonal reflection R | complex conjugation | ℤ₂ involution |
| self-complementary | grid-symmetric (R-fixed) | conjugation-symmetric spectrum | **Brouwer fixed point → SOS** |
| transitive | all-forward | roots of unity (odd n) / prim 2n-th (even) | fixed-point-rich |
| **Paley, p≡3 mod 4** | max-symmetry heptagon | primitive roots, 1 real point (=vertex 0) | **free ℤ₂ → Borsuk–Ulam → 3 pillars** |
| c₃ | grid triangles | `Tr(A³)/3` spectral moment | — |

## Status
- **Verified:** `c₃ = Tr(A³)/3` (n≤5); `p≡3 mod 4 ⟺ −1∈QNR ⟺ negation is a free ℤ₂ on QR-structure ⟺ Paley
  tournament exists` (p=3..23); the `i` is intrinsic to the skew/Cayley map; `ρ=R∘V` is not a √complement.
- **Framing (compelling, partly heuristic):** complement is a ℤ₂ involution; **fixed point (Brouwer) ⇒
  symmetric SOS exact (Gatermann–Parrilo) = the easy regime**; **free ℤ₂ (p≡3 mod 4, Paley) ⇒ Borsuk–Ulam
  obstruction to a symmetric certificate = the hard regime, where the three pillars (POCS / flat-extension /
  Blaschke) are the constructive tools.** The verified backbone is solid; "free ℤ₂ ⇒ SOS provably fails" is the
  organizing conjecture (Borsuk–Ulam is the obstruction *principle*, not yet a theorem for the LRC functional).
- **Payoff:** a criterion for *which* LRC/tournament extremality problems are SOS-tractable (symmetric/p≡1) vs
  need the constructive pillars (free/p≡3 mod 4). The AP bound is Brouwer/SOS; the covering-min and the Paley
  flip-rank obstruction are Borsuk–Ulam / three-pillar.

Related: HYP-3813 (the algebra↔geometry dictionary — this adds the topology column + c₃), HYP-3805 (Paley
flip-rank obstruction), HYP-3811 (odd-n anti-auto fixed vertex = the real eigenvalue), HYP-3791 (MOSEK: no SOS
shortcut to the for-all-S crux — now explained as the free-ℤ₂ obstruction), HYP-3796/mac-mini (the three
pillars), HYP-3812/klein (fold diagnoses obstruction; metric-irreducible at Φ₆), THM-503 (L not an Euler
product). HYP-3814 (this). Script: 04-computation/topo_brouwer_borsukulam_c3_opus_20260701.py.
