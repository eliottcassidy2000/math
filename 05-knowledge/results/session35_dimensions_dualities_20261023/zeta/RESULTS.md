# Session 35: tournament spectra, theta lifting, and the two different origins of 2

Status: PROVED elementary statements and standard theta/Mellin consequences below; finite checks explicitly separated. No RH claim. Prepared outside the repository; no Git mutation.

## Inheritance and correction ledger

Closest inherited mechanism: `01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md` defines the adjacency dynamical zeta `1/det(I-uA)`, with closed-walk Euler product and strong-component factorization. Its Paley adjacency modulus is wrong in both the YAML status and section 4: eigenvalues are `(−1 ± i sqrt(p))/2`, hence modulus `sqrt(p+1)/2`, not `sqrt(p)`. The latter belongs to skew eigenvalues `±i sqrt(p)`. This does not damage the determinant/Euler-product theorem. Exact characteristic-polynomial checks below cover p=3,7,11,19.

Less-used sidecar: `01-canon/theorems/THM-1760-the-tournament-pencil-three-faces-functional-equation.md` already requires an adjugate sidecar for its complement equation and explicitly makes no completion claim. Likewise `05-knowledge/hypotheses/HYP-2645-lrc14-poisson-theta-dual.md` warns that a valid Fourier/Poisson translation can merely recover the existing representation, without creating the desired analytic bound. The relevant conceptual lesson is that a named duality does not supply its missing analytic sidecar.

Concept board: (1) adjacency dynamical zeta / periodic walks; (2) positive skew-square spectrum / finite heat trace; (3) binary tower / discrete scale invariance; (4) integral lattice / Fourier inversion; (5) centered reflection / localization of zeros. Never substitute one object for another without the map below.

## 1. Finite-spectrum completion obstruction

Let M be any real skew-symmetric n-by-n matrix, n>=1, and Q=I−M²=I+MᵀM. Thus Q is positive definite, with eigenvalues q_j>=1. Define the finite energy spectral zeta

`Z_Q(s)=tr(Q^(−s))=Σ_j exp(−s log q_j)`.

It is entire. For every integer m>=0,

`Z_Q(−m)=tr(Q^m)>0`.

Consequently `Γ(s)Z_Q(s)` has a genuine simple pole at each nonpositive integer, with residue `(-1)^m tr(Q^m)/m!`. Equivalently the length spectral zeta `L_Q(s)=Z_Q(s/2)` has `L_Q(−2m)>0`, and `π^(−s/2)Γ(s/2)L_Q(s)` has genuine poles at all `s=−2m`. This is the exact opposite of Riemann's trivial-zero cancellation.

The finite heat trace `h_Q(t)=tr exp(−tQ)` has the entire Taylor expansion `Σ_m (-t)^m tr(Q^m)/m!`. Its nonzero coefficients are the Mellin pole residues. No finite positive Q,Q' and d>0 can satisfy `1+h_Q(t)=c t^(−d/2)(1+h_Q'(1/t))` for fixed c>0: the left tends to n+1 as t→0, the right diverges. Without the 1's the right tends to zero, again impossible. A finite heat trace cannot be silently identified with a lattice theta function.

The adjacency dynamical zeta is a third object: as a rational function `1/det(I-uA)` it has no finite zeros at all. Its interesting singularities are poles. Complement invariance A→Aᵀ is an identity at the same argument, not a functional equation reflecting an analytic variable.

## 2. Binary tower: exact scaling creates a pole line, not a zero line

The actual skew doubling operation is

`D(M)=[[M,M+I],[M−I,−M]]`, so `Q(DM)=2(Q⊕Q)`.

For the direct sum of all finite generations k>=0, take the infinite positive operator with q_j 2^k repeated 2^k times. Its spectral zeta converges absolutely for Re(s)>1 and is exactly

`Z_tower(s)=Z_Q(s)/(1−2^(1−s))`.

This gives a meromorphic continuation. Potential poles lie at `s_k=1+2πik/log 2`; a pole is canceled exactly when `Z_Q(s_k)=0`. At k=0 cancellation is impossible, and the residue is `Z_Q(1)/log 2>0`. At every nonpositive integer, `Z_tower(−m)=tr(Q^m)/(1−2^(m+1))<0`; again there are no trivial zeros. At length-variable p=2s the abscissa is p=2, so the tower has spectral dimension 2. This is a rigorous sense in which doubling multiplicity and doubling squared energy select a dimension, but it does not select Riemann's critical line.

A sharper time-domain statement exposes the missing inversion. Let

`H(t)=Σ_(k>=0) 2^k h_Q(2^k t)`,
`P(t)=Σ_(k in Z) 2^k t h_Q(2^k t)`.

The latter converges absolutely for t>0 and satisfies P(2t)=P(t). Splitting the negative k terms and Taylor expanding their exponentials gives the exact convergent formula

`H(t)=P(t)/t − Σ_(m>=0) [(-t)^m tr(Q^m)/(m! (2^(m+1)−1))]`.

In particular `H(t)=P(t)/t−n+O(t)`. The leading coefficient is log-periodic rather than a forced constant. Discrete dilation t→2t is not Fourier inversion t→1/t. For Q=[1], every vertical candidate pole is genuine because Z_Q≡1; this is the cheapest hostile control against a presumed reflection law.

## 3. Constructive bridge: a genuine lattice theta lift of every tournament

Now explicitly add the missing infinite lattice. For a skew tournament matrix M, define

`Θ_M(t)=Σ_(v in Z^n) exp(−πt vᵀQv)`, `Q=I−M²`.

Poisson summation for the Gaussian gives

`Θ_Q(t)=det(Q)^(−1/2)t^(−n/2) Θ_(Q⁻¹)(1/t)`.

This is a genuine, dimension-dependent inversion law. Its Fourier dual is Q⁻¹; edge reversal M→−M leaves Q unchanged and is not this duality. Q⁻¹ need not equal the Q-matrix of any tournament, so the target category must be enlarged to positive quadratic forms.

The *actual tournament doubling* becomes the exact theta identity

`Θ_(D(M))(t)=Θ_M(2t)^2`.

Proof: the quadratic form on Z^(2n) is `2vᵀQv+2wᵀQw`; independent summation factors. Thus doubling induces convolution of lattice representation counts, not merely doubling of the finite spectral trace. This is the precise operation-preserving bridge.

Define the length-normalized Epstein zeta

`E_Q(s)=Σ_(v≠0) (vᵀQv)^(−s/2)`, Re(s)>n,
`Λ_Q(s)=π^(−s/2)Γ(s/2)E_Q(s)`.

Termwise Mellin integration gives `Λ_Q(s)=∫_0^∞(Θ_Q(t)−1)t^(s/2−1)dt` in that half-plane. Splitting at 1 and applying Poisson on (0,1) yields

`Λ_Q(s)=−2/s+2 det(Q)^(−1/2)/(s−n)`
` + ∫_1^∞(Θ_Q(t)−1)t^(s/2−1)dt`
` + det(Q)^(−1/2)∫_1^∞(Θ_(Q⁻¹)(t)−1)t^((n−s)/2−1)dt`.

The two integrals are entire. Hence Λ_Q has only simple poles at 0,n and obeys

`Λ_Q(s)=det(Q)^(−1/2)Λ_(Q⁻¹)(n−s)`.

At s=−2m, m>=1, the right side is strictly positive and finite (positive Dirichlet series at n+2m). Since Γ(s/2) has a simple pole there, E_Q has a simple zero there. At s=0, Λ_Q itself has a pole and E_Q(0)=−1, so 0 is not a trivial zero. This establishes precisely where the negative-even progression comes from: the archimedean Gamma factor plus analytic continuation, not binary graph growth.

For a self-reflecting completion, normalize B=Q/det(Q)^(1/n), so det B=1, then set `C_B(s)=Λ_B(s)+Λ_(B⁻¹)(s)`. It obeys C_B(s)=C_B(n−s). Multiplying by s(s−n) removes the two poles and preserves reflection. This construction does not constrain all zeros to Re(s)=n/2.

Dimension and scaling: n is the lattice dimension, and its length-variable reflection center is n/2. Increasing tournament order by 2 changes that center by +1; doubling order changes it from n/2 to n. Reparameterizing w=s/n moves the center to 1/2 but moves trivial zeros to w=−2m/n. Except in dimension one, one cannot simultaneously preserve the raw Riemann labels 1/2 and −2m by this normalization. At n=1, Q=[1], E_Q(s)=2ζ(s), recovering exactly the classical case.

## 4. Reflection plus positivity is still insufficient

For any a,b>0, the entire real polynomial

`F(s)=((s−1/2−a)^2+b²)((s−1/2+a)^2+b²)`

obeys F(s)=F(1−s), is strictly positive for every real s, and has zeros `1/2±a±ib` off the symmetry axis. Choose a=1/4,b=1 to keep all roots inside the critical strip. This does not refute RH; it refutes an inference from reflection, reality, positivity on the real axis, or symmetry of paired factors to RH. A successful localization mechanism needs additional information not present in these properties.

## 5. Loss ledger and boundary

Source: tournament M with labeled integer vertex coordinates. Target: positive integral Q and its full lattice theta/Epstein pair. Map M→I−M²→(Z^n,Q). Preserved: exact skew doubling via Θ_D(t)=Θ_M(2t)²; reversal and signed-permutation switching invariance; positivity; dimension. Lost: orientation/sign of M and generally tournament reconstruction. Added sidecar: all integral lattice vectors plus Fourier-dual quadratic form Q⁻¹. Not preserved: the original closed-walk Euler product or its primitive cycles; no claim of Euler product for arbitrary Epstein zeta; no RH implication. In the scalar-Q families, determinant normalization sends every Q=cI to I, explicitly collapsing all scale information.

Cheapest tests: n=1 and its two-vertex double (theta squaring); Q=[1] tower (uncanceled vertical poles); three-vertex cycle and transitive tournament (different graph zeta, yet test finite positivity); a=1/4,b=1 centered polynomial (off-axis roots). General analytic proofs above are not inferred from finite enumeration.

## 6. Reproduction and sources

Run `python3 session35-zeta/check.py`; see check.out. Exact labeled tournament universe n=1..5 checks Q doubling and negative moments; exact Paley adjacency characteristic polynomials; exact polynomial hostile and first 21 theta coefficients. Scalar heat identity checks use 45-digit arithmetic and are labeled NUMERICAL, not a proof.

Primary/reference sources checked 2026-10-23:
- NIST DLMF §25.4, especially equations 25.4.2–25.4.4: https://dlmf.nist.gov/25.4 (Riemann reflection and completed xi).
- NIST DLMF §25.10: https://dlmf.nist.gov/25.10 (negative-even trivial zeros and RH distinction).
- NIST DLMF §5.2: https://dlmf.nist.gov/5.2 (Gamma poles, zero-free Gamma, reciprocal Gamma zeros).
- NIST DLMF §20.7(viii): https://dlmf.nist.gov/20.7#viii (theta modular transformations).
- Elkies, Harvard Math 229 course notes and linked zeta1.pdf: https://people.math.harvard.edu/~elkies/M229.18/index.html (Poisson/theta derivation). The generalized matrix argument is supplied in full above rather than relying on any secondary-source formula.

Open direction, explicitly not proved: determine when the lattice theta lift retains distinguishing tournament information beyond the finite Q spectrum, and classify tournament Q-forms whose determinant-normalized lattices are isodual. The scalar-Q family is a known collapse, not evidence for injectivity.
