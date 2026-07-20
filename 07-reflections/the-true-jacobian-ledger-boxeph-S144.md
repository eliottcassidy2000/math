# The true-Jacobian ledger: what survives the collapse (boxeph-2026-07-20-S144)

Owner brief: pull from the fleet, fully understand the counterexamples, or find
sharper versions of JC that are TRUE. The overnight swarm covered the anatomy
(klein: symplectic ℂ⁶ lift; death-star: THM-1300 Weyl endomorphism + THM-1305
equivariant anatomy with k=3 emptiness by mod-p/Newton — independently matching
my S143 D/W emptiness, two methods, one theorem; mac-mini: THM-1315 surjectivity,
the fiber cubic, S₃ pinned by syzygy, the caustic; kind-pasteur: verification +
Chebotarev; opus-S414: normal-form/moduli lane). The untaken lane was the
SALVAGE ledger — what remains TRUE. Assembled here, with two new results.

## A. Classical floor (cited, placed)
1. Keller + proper ⟹ automorphism (étale + proper over simply-connected = iso).
2. Keller + injective ⟹ automorphism (injective polynomial maps are automorphisms).
3. Keller of degree ≤ 2 ⟹ automorphism (Wang). The kernel has degree 7.
4. Filtration-preserving Weyl endomorphisms are automorphisms — and our φ_F
   visibly violates filtration preservation (degree growth): consistent.

## B. NEW: the dimension-2 mechanism no-go (theorem, verified 8/8 + 4-line proof)
In dim 2 the equivariant class is F = (yA(s), xh(s)), s = x^k y, and
det JF = −[Ah + s(kAh′ + A′h)] (identity, verified for k = 1..4, all degree
shapes). Keller forces the degree-(a+b) leading coefficient (1+a+kb)·A_a·h_b to
vanish — impossible — so A, h are constants and F is linear. **The machinery
that killed JC₃ provably cannot descend to ℂ²: the JC₂ island is safe from the
entire equivariant z-linear mechanism.**

## C. NEW: the mod-p decision theorem (salvage; Jordan + Chebotarev + Lang–Weil)
**A Keller map over a number field is an automorphism ⟺ it is bijective mod p
for all sufficiently large p ⟺ for infinitely many p.** Not-auto ⟹ non-injective
(A.2) ⟹ generic degree d ≥ 2 ⟹ the Galois closure's group acts transitively on
d ≥ 2 sheets ⟹ (Jordan) it has a fixed-point-free element ⟹ (Chebotarev +
Lang–Weil) deficiency ≥ δp^n + O(p^{n−1/2}) ⟹ non-bijective for all large p.
Conversely automorphisms are bijective mod every good prime. COROLLARY: one
bijective large prime CERTIFIES automorphy — the S140 census was not evidence
but a decision procedure; demo: the kernel's deficiency/p³ = .32325, .32480,
.32559 at p = 31, 37, 41, monotone to the S₃ value 1/3. JC is FALSE
geometrically but TRUE as an arithmetic test: automorphy is mod-p decidable.

## D. In-class classification corollaries (S142/S143 + THM-1305, convergent)
m=2 class: Keller ⟹ automorphism-type or the (1,1;6,4)-kernel orbit (unique).
m=3 class: Keller ⟹ automorphism-type (empty; DOUBLE independent proof).
Corollary: **Galois Keller in-class ⟹ automorphism** (kernels are S₃, non-Galois).

## E. Conjectures (labeled)
(a) kernels exist only at z-weight 2 (test m = 4, 5 — D/W calculus is routine);
(b) every Keller non-automorphism has FULL symmetric monodromy S_d;
(c) the decision constant in C is effectively computable from deg F alone;
(d) WILD: Keller + empty Jelonek set ⟺ proper (known) — and the kernel's
    asymptotic set is the caustic mac-mini located: escape IS ramification at ∞.

Files: jacobian_true_versions_boxeph_S144.py + .out (both new results frozen).
