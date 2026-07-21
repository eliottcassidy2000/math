---
id: THM-1925
title: "THE REDUCTION IS A PRODUCT OF TRIGONOMETRIC FUNCTIONS — a trigonometric/harmonic lens unifying the repo's reduction principles. (a) SPECTRAL REDUCTION: char_A is multiplicative under order-join (block-triangular), so char_A(T)=∏ char_A(strong component) and the adjacency spectrum is the disjoint union of the strong components' spectra (EXACT, all 74 iso classes n≤6). (b) SINE-PRODUCT: the signed tournament partition function Σ_T(−1)^{back}x^{score}=∏_{i<j}(x_j−x_i); on the unit circle x_k=e^{iθ_k} it is ∏ 2i·sin((θ_j−θ_i)/2) — a product of sines, so mac-mini's transitive-core involution IS a trigonometric factorization; at the n-th roots of unity |∏(ωʲ−ωⁱ)|=n^{n/2}=√|disc(xⁿ−1)|. (c) TRIGONOMETRIC ATOMS: circulant strong tournaments have eigenvalues that are explicit trig sums — Gauss sums (Re=−1/2) for Paley, and Dirichlet-kernel = Chebyshev-U values (Dₘ(2πj/n)−1)/2 = (U₂ₘ(cos πj/n)−1)/2 for interval connection sets. UNIFYING FRAME: reduction principles are decompositions along the atoms' group characters, and characters are trigonometric — the tournament mirror of the LRC covering sum ∏sinc / Fejér-certificate / Chebyshev-equioscillation side."
status: >
  VERIFIED (boxeph-2026-07-21-S194). (a) EXACT via Faddeev–LeVerrier integer characteristic
  polynomials over all 74 iso classes n≤6 (0 mismatches); the order-join case is exact to 0.0
  (block-triangular). (b) signed-partition = Vandermonde verified to 1e-13 (n=3,4,5), the circle
  sine-product to 4e-14, and |∏(ωʲ−ωⁱ)|=n^{n/2} exactly (n=5,7). (c) Paley-7/11/19 non-Perron
  Re=−1/2 (Gauss sum) and interval-C={1..m} Re(λ_j)=(Dirichlet_m(2πj/n)−1)/2 verified exactly
  (n=7,9,11). The synthesis (a)–(c) as one "product over trigonometric atoms" principle, and its
  bridge to the LRC side, is a UNIFYING LENS (reflection), resting on these verified pieces plus
  standard facts (char-poly of block-triangular matrices; Vandermonde; Dirichlet=Chebyshev-U;
  Gauss sums). char_A multiplicativity under reducibility is already implicit in THM-1830; stated
  here as the spectral companion of THM-1862 with the explicit trigonometric atoms.
source: boxeph-2026-07-21-S194 (owner: another archeology session; pursue more reduction principles; think trigonometric)
depends_on: []
related:
  - THM-1862  # order-join reduction principle (this is its spectral/trigonometric companion)
  - THM-1830  # reducible => char_A = product of char(SCC); SCCs are the atoms
  - THM-456   # blow-up (cycle-length) spectrum law — a different blow-up reduction
  - "07-reflections/the-reduction-is-a-product-of-trigonometric-functions-boxeph-S194.md"
  - "07-reflections/the-covering-min-is-a-chebyshev-equioscillation-and-why-greedy-has-no-shortcut.md"
  - "07-reflections/the-cyclotomic-magic-function-is-the-fejer-kernel-kps.md"
  - "07-reflections/the-sign-reversing-tournament-involution-as-a-repo-wide-engine-macmini-S159.md"
script: 04-computation/trig_reduction_tournament_boxeph_S194.py (+ .out)
---

# THM-1875 — the reduction is a product of trigonometric functions

`A[i][j]=1` means *i beats j*; scores are out-degrees; `back(T)` = # arcs from the smaller to the
larger index. The repo's reductions ("reduce to the strong core", "reduce to the covering core")
are the same move — **factor a global object over irreducible atoms** — and on the unit circle each
factorization is a **product of trigonometric functions**.

## (a) Spectral reduction — char_A is multiplicative over strong components

In condensation order a tournament's adjacency matrix is **block upper-triangular** with the
strong-component adjacency matrices on the diagonal. Hence
```
        char_A(T) = ∏_{S : strong component} char_A(S),   spec(T) = ⨆_S spec(S).
```
Verified EXACT (integer Faddeev–LeVerrier) on all 74 iso classes n≤6, and to 0.0 for every
order-join `T₁▷T₂` (n=3×3). The **spectral atoms are the strongly connected tournaments** — the
spectral companion of the order-join reduction THM-1862 (and implicit in THM-1830). *(numpy's float
`eigvals` is off by ~√3 on the nilpotent transitive blocks — eigenvector ill-conditioning, not a
real gap; the exact char-poly check is definitive.)*

## (b) The transitive-core reduction is a product of sines

mac-mini's signed involution engine is the Vandermonde:
```
        Σ_{T}  (−1)^{back(T)}  ∏_k x_k^{score_k(T)}  =  ∏_{i<j} (x_j − x_i).
```
On the unit circle `x_k = e^{iθ_k}`,
```
        ∏_{i<j}(x_j − x_i) = ∏_{i<j} 2i · sin((θ_j − θ_i)/2) · e^{i(θ_i+θ_j)/2},
```
a **product of sines**. So the sign-reversing involution that collapses all tournaments to the
transitive survivors *is* the trigonometric factorization of the Vandermonde. At the n-th roots of
unity the surviving product has modulus `|∏_{i<j}(ωʲ − ωⁱ)| = n^{n/2} = √|disc(xⁿ−1)|` (verified
n=5,7) — the transitive core reduces to **cyclotomic** arithmetic.

## (c) The atoms are trigonometric — Gauss sums and Chebyshev/Dirichlet kernels

A circulant (rotational) strong tournament on ℤ/n with connection set `C` (`i` beats `i+k`,
`k∈C`) has eigenvalues `λ_j = Σ_{k∈C} ωʲᵏ` (`ω=e^{2πi/n}`) — a **character sum**:

- **Paley** (`C` = quadratic residues, `n≡3 mod 4`): `λ_j` for `j≠0` is a **Gauss sum**, with
  `Re λ_j = −1/2` (verified n=7,11,19) — the identity `Σ_{k∈QR} cos(2πk/n) = −1/2`.
- **Interval** (`C={1,…,m}`, `m=(n−1)/2`): `Re λ_j = (D_m(2πj/n) − 1)/2` where `D_m(θ)=
  sin((m+½)θ)/sin(θ/2)` is the **Dirichlet kernel** `= U₂ₘ(cos(θ/2))`, a **Chebyshev polynomial of
  the second kind** (verified exactly n=7,9,11).

> **REFINED by THM-1955 (boxeph-S196):** in fact `Re λ_j = −1/2` for **every** ℤ/p circulant
> tournament (j≠0) — the pair `{k, p−k}` shares a cosine, so the real part is `½·Σcos = −1/2`
> regardless of `C`. The interval case gives `−1/2` because the Dirichlet kernel `D_m(2πj/n)=0` at
> the roots of unity. So Paley and interval both sit on `Re=−1/2`; they differ only in the
> **imaginary (sine) part** — Paley is *flat* (`|λ|=√p`, Gauss sum), interval is *spread* (Dirichlet).

Via (a) the spectrum of *any* reducible tournament is assembled from these trigonometric atoms.

## The unifying frame (why trig is everywhere)

All three are **decompositions along group characters**, and group characters are trigonometric:
- circulant tournament ↔ ℤ/n action → eigenvalues are `χ_j` on `C` = trig sums (a,c);
- the Vandermonde ↔ the sign character of `S_n` → sine products on the torus (b);
- **LRC** ↔ the torus ℝ/ℤ → the covering sum `Σ_{k·v=0} ∏_j sinc(k_j δ)` is a product of sincs, the
  covering-min is a **Chebyshev equioscillation**, and the certificate is the **Fejér kernel**
  `= |Dirichlet|²` (kps, `the-cyclotomic-magic-function-is-the-fejer-kernel`).

The **Dirichlet/Fejér kernel is the shared object**: it is the interval-circulant tournament's
eigenvalue (c) and the LRC certificate (Fejér). Reduction principles are character decompositions;
the hard analytic content lives where the trigonometric product **vanishes** (LRC covering nullcone)
or **equioscillates** (covering-min) — exactly where the tournament sine-product hits a zero.
