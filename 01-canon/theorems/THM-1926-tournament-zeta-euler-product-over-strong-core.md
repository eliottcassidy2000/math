---
id: THM-1926
title: "THE TOURNAMENT ZETA CONCENTRATES ON THE STRONG CORE. ζ_T(u)=1/det(I−uA) is the Bowen–Lanford/Ruelle dynamical zeta of the arc-subshift; it has the Euler product ∏_{primitive cycles p}(1−u^{ℓ(p)})^{−1} (prime-cycle counts π_ℓ=(1/ℓ)Σ_{d|ℓ}μ(ℓ/d)N_d are non-negative integers), factors as ∏ over strong components, and EQUALS 1 on the acyclic/transitive part (A nilpotent) — the zeta only sees the non-wandering set. It is complement-invariant (ζ_T=ζ_{T^op}). The trace formula N_k=tr(A^k)=#closed k-walks satisfies N_1=N_2=0 (no loops, no digons) and N_3=3c3, so ζ_T(u)=exp(c3·u³+…) starts at u³ — the 3-cycle is the fundamental prime and c3 counts them. Poles at 1/λ = reciprocals of the trigonometric atom-spectra (Gauss sums for Paley, Chebyshev/Dirichlet for interval circulants). The harmonic-analysis realization of the strong-core reduction (THM-1925/THM-1862): periodic-orbit↔spectrum duality with trigonometric atom eigenvalues."
status: >
  VERIFIED (boxeph-2026-07-21-S195), exact integer arithmetic over all 74 iso classes n≤6:
  det(I−uA)=∏_SCC det(I−uA(SCC)) (0 mismatches); acyclic ⇒ det(I−uA)=1 ⇒ ζ=1; prime-cycle counts
  π_ℓ non-negative integers (Möbius inversion of N_k); complement-invariance ζ_T=ζ_{T^op}; N_3=3c3;
  N_1=N_2=0. Circulant atoms: N_k=tr(A^k)=Σ_j λ_j^k (Paley n=7,11,19, |λ_nonPerron|=√n Gauss sum).
  Standard facts used: det(I−uA)=Σ char_A-coeffs·u^j; the Bowen–Lanford Euler product for a subshift
  of finite type; block-triangular determinant. char_A(T)=∏char_A(SCC) is the spectral side of
  THM-1925/THM-1830; new here is the closed-orbit/Euler-product reading, the ζ=1-on-wandering-set
  concentration, complement-invariance, and the c3/N_k dynamical dictionary.
source: boxeph-2026-07-21-S195 (owner: harmonic-analysis lens; work open handoffs; integrate fleet ideas)
depends_on: []
related:
  - THM-1925  # my trig-reduction: char_A=∏char(SCC); circulant atoms=Gauss/Chebyshev (the spectral side)
  - THM-1862  # order-join reduction principle
  - THM-1875  # kps: transitive char_S cot ladder (the acyclic atom on the skew side — ζ=1 here)
  - THM-1880  # kps: a/b monoid; Paley char_S=x(x^2+p)^{(p-1)/2} Gauss sum (poles 1/λ, |λ|=√p)
  - THM-1870  # kps: cycle-count spectral boundary = Hamiltonian length (the N_k here)
  - "07-reflections/the-reduction-is-a-product-of-trigonometric-functions-boxeph-S194.md"
script: 04-computation/tournament_zeta_boxeph_S195.py (+ .out)
---

# THM-1926 — the tournament zeta concentrates on the strong core

A tournament `T` is a **subshift of finite type**: the arc set is a dynamics, and closed walks are
its periodic orbits. Its **Bowen–Lanford / Ruelle dynamical zeta** is
```
        ζ_T(u) = 1 / det(I − uA) = exp( Σ_{k≥1} N_k u^k / k ),   N_k = tr(A^k) = #closed k-walks,
```
and `det(I − uA) = Σ_j c_j u^j` with `c_j` the characteristic-polynomial coefficients (integers).

## 1. Euler product over primitive cycles (the "primes")

`ζ_T(u) = ∏_{primitive cycles p} (1 − u^{ℓ(p)})^{−1}`, where the primitive-cycle counts
`π_ℓ = (1/ℓ) Σ_{d|ℓ} μ(ℓ/d) N_d` are **non-negative integers** (verified all 74 classes n≤6).
This is the tournament analogue of `ζ_Riemann = ∏_p (1−p^{−s})^{−1}`: primitive cycles are the primes.

## 2. Concentration on the strong core (the reduction)

`det(I − uA(T)) = ∏_{strong components S} det(I − uA(S))` (block-triangular; 0 mismatches n≤6), and
for an **acyclic** (transitive) tournament `A` is nilpotent so `det(I − uA) = 1` and **ζ_T ≡ 1**.
The zeta is **invisible to the acyclic part** — it only sees the **non-wandering set** (the strong
components). This is the harmonic/dynamical form of the strong-core reduction (THM-1925/THM-1862):
not just "the content lives in the strong core", but "the zeta's support *is* the strong core."

## 3. The 3-cycle is the fundamental prime

Every tournament has **N₁ = 0** (no loops) and **N₂ = tr(A²) = 0** (no digons), so
`ζ_T(u) = exp(c₃·u³ + (N₄/4)u⁴ + …)` — the expansion **starts at u³**. The smallest primitive cycle
is the 3-cycle, and `N₃ = tr(A³) = 3·c₃`, so **c₃ is read straight off ζ's u³ log-coefficient**.
The 3-cycle atom has `ζ = 1/(1 − u³)` (a single length-3 prime) — the simplest nontrivial zeta, and
exactly the intransitivity atom of THM-1830. The higher `N_k` are kps's cycle-count spectral data
(THM-1870).

## 4. Complement-invariance and the trigonometric poles

`A(T^op) = Aᵀ` and `det(I − uAᵀ) = det(I − uA)`, so **ζ_T = ζ_{T^op}** (verified n≤6) — a
functional-equation-like symmetry tying ζ to the SC/complement (blue) lore. The **poles are at
`u = 1/λ`** for eigenvalues `λ`; the trace formula `N_k = Σ_j λ_j^k` is the periodic-orbit ↔
spectrum duality (the tournament "explicit formula"). For the circulant atoms the `λ_j` are
character sums: **Gauss sums** for Paley (`|λ| = √n`, Re = −½; verified n=7,11,19), **Dirichlet =
Chebyshev-U** for interval circulants (THM-1925). So the zeta's poles are reciprocals of
trigonometric numbers.

## The one line

`ζ_T = 1/det(I − uA)` turns the strong-core reduction into a **periodic-orbit Euler product**: it is
trivial on the acyclic (wandering) part, its primes start at the 3-cycle with `c₃` the first count,
and its poles are the reciprocal Gauss/Chebyshev atom-spectra — the harmonic-analysis face of
THM-1925, and the closed-orbit companion of kps's char_S cot/Gauss ladder.
