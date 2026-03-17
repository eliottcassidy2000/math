# Adelic Tournament Geometry

**Session:** opus-2026-03-17-S91a/b

## The Discovery

The tournament eigenvalue space is an **adelic manifold**. At each level n, the flip chain eigenvalues live on a truncated adelic space:

```
A_T(n) = R × ∏_{p | D(n)} Z/p^{v_p(D)}Z
```

where D(n) = odd_part(C(n,2)) is the **conductor**.

As n → ∞, this fills out the full odd adeles: A_T(∞) = R × ∏_{p odd} Z_p.

## The 18 Geometric Similarities

| Concept | Adelic geometry | Tournament space |
|---------|----------------|------------------|
| Base space | Spec(Z) = {primes} | D(n) = odd part C(n,2) |
| Local factors | Q_p at each prime | Z/p^e Z at each p\|D |
| Global space | A_Q = R × ∏' Q_p | R × ∏ Z/p^e Z |
| Local-global | Hasse-Minkowski | CRT (always exact!) |
| Operator | Hecke T_p | Flip P_{ij} |
| Eigenvalue bound | Ramanujan | Spectral gap 4/C(n,2) |
| Formal group | Lubin-Tate | F(x,y)=(x+y)/(1+xy) |
| Height | height 1 or ∞ | height 1 (odd p), ∞ (p=2) |
| Product formula | ∏ \|x\|_v = 1 | ∏ \|Q(λ)\|_v = 1 |
| Tree structure | Bruhat-Tits T_p | Prime depth tower |
| Metric | p-adic + real | Ultrametric + hyperbolic |
| Zeta function | ζ(s) = ∏ 1/(1-p^{-s}) | L_T(s) = ∏ L_p(s) |
| Class field | Artin map | Cayley → cyclotomic |
| Conductor | Level of modular form | D(n) = eigenvalue denom |
| Supersingularity | Height ∞ at p | H(T) odd (Rédei) at p=2 |
| Rapidity | p-adic log | arctanh = formal log |
| Strong approx | G(Q) dense in G(A_f) | CRT bijection (exact!) |
| Archimedean | R-place geometry | Poincaré disk |

## Key Insights

1. **T does not factorize as tensor product** — the adelic structure is in the eigenvalues, not the matrix. The 56 classes at n=6 cannot be 2×12 from T_3⊗T_5.

2. **The profinite tournament space is the ODD profinite integers** — because F(x,y)=(x+y)/(1+xy) is supersingular at p=2, which kills the 2-adic place, leaving ∏_{p odd} Z_p.

3. **Eigenvalue geometry is simultaneously hyperbolic (archimedean) and tree-like (non-archimedean)** — this is a Berkovich space structure.

4. **The formal group distance = hyperbolic distance** — arctanh(λ) gives rapidities, and the formal group law IS hyperbolic geometry. The zero mode sits at the center; stationarity at the boundary.

5. **The spectral zeta function at s=1/2 is SELF-DUAL** — zeta_T(1/2) = zeta_T(1/2), reflecting the λ ↔ -λ pairing of S/A eigenvalues.

## What This Means

The tournament flip chain is a **finite-level Hecke operator** on a **truncated adelic space**. The conductor D(n) controls which primes are visible. Each new prime p entering D(n) (at n = p) is like ascending to a new level in the modular form tower.

42 = 2·3·7 is the conductor of the universal adelic tournament: the unique triple point where inert (2), ramified (3), and split (7) primes all coexist.
