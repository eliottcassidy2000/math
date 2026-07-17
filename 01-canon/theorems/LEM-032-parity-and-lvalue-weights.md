---
id: LEM-032
title: THE PARITY LAW, THE SUPPORT LAW, AND THE L-VALUE WEIGHTS — the structure theory of the frame-cross spectrum ĉ(χ). (A) PARITY: W(−a)=W(a), X(−m)=X(m) ⟹ cross(−w)=cross(w) ⟹ ĉ(χ)=χ(−1)ĉ(χ): every ODD character carries ZERO mass (and per-factor Ŵ_g(χ)=X̂_g(χ)=0 — double kill). DECODE OF S60: Legendre mod p is odd iff p≡3 (4); 7≡3≡3 (4) ⟹ the S60 vanishings at mod 7/mod 3 are PARITY zeros — NOT the seven-section structure (LEM-031's conjecture CORRECTED): the even mod-7 characters (the cubics) carry mass on all three referee clusters. (B) SUPPORT: the class-g term of ĉ(χ) vanishes unless cond(χ) | P/g (fiber orthogonality). (C) THE TWISTED JORDAN LEMMA: T_q(χ) = Σ_{u∈(Z/q)*}χ(u)csc²(πu/q) = 0 (χ odd) = J₂(q)/3 (χ₀: THM-892(C*) is the trivial row) = (2q²/π²)·L_q(2,χ) (χ even; Euler factors at p|q removed) = 2τ(χ)B_{2,χ̄} (χ even primitive) — proof: tent-kernel inversion (THM-892(K)) + Gauss sum + quadratic Bernoulli moment; the partial-fraction face gives the L-value; the classical evaluation L(2,χ)=π²τ(χ)B_{2,χ̄}/q² drops out self-contained. (D) THE CLOSED-FORM WEIGHT SIDE: Ŵ_g(χ) = (2/g²)·L_{P/g}(2,χ) — pure 1/g² times a Dirichlet L-value at s=2. The factorization law (LEM-031) becomes ĉ(χ) = Σ_{g|P, g<P, cond(χ)|P/g} (2/g²)·L_{P/g}(2,χ)·X̂_g(χ): the weight side is CLOSED FORM; the ONLY cluster-dependent factor left is X̂_g (twisted coincidence sums). REFEREED EXACT everywhere (worst 7·10⁻¹⁴)
status: PROVED ((A) 3 lines; (B) fiber orthogonality; (C) kernel inversion + Gauss + Bernoulli, both faces; (D) assembly) + REFEREED EXACT — Part C battery: 144 characters over q∈{5,7,9,12,13,15,20,36,63,180}, all four faces (odd max 3e-11, L-face 9e-15, Gauss–Bernoulli 8e-15); exact-rational faces T₅=8/√5 (B_{2,χ}=4/5), T₁₃=8√13 (B_{2,χ}=4); parity law at machine zero on 336+144+48 odd characters across three clusters; Parseval exact (masses sum = frame variance) on all three; closed-form spectrum reproduction worst 7e-14; S60 regression −8.3403/var 3563.1 reproduced (two-owner cluster)
source: boxeph-2026-07-17-S61 (owner directive: prove the mod-7 vanishing; evaluate the Ŵ_g character sums)
depends_on: [LEM-031 (the factorization law being closed-formed), THM-892 ((K) kernel identity, (C*) = the trivial row), LEM-030 (the χ₀/mean row)]
related: [THM-884/879 (the N(h) coincidence family — X̂_g is its χ-twist), klein's discrepancy lane (the character masses are the factorized discrepancy objects)]
script: 04-computation/lrc14_parity_lvalue_boxeph_S61.py -> 05-knowledge/results/lrc14_parity_lvalue_boxeph_S61.out
---

# LEM-032 — parity, support, and the L-value weights

Setting (LEM-031): cross(w) = Σ_{a≠0} W(a)X(aw) on w ∈ (Z/P)\*, with
W(a) = (π²/P²)csc²(πa/P), X(m) = |Σ_e S_e(m)|² − Σ_e |S_e(m)|²,
ĉ(χ) = φ(P)⁻¹Σ_w cross(w)χ̄(w), and per gcd-class g (q := P/g):
Ŵ_g(χ) = Σ_{u∈U_q} W(gu)χ(u) (unnormalized), X̂_g(χ) = φ(q)⁻¹Σ_u X(gu)χ̄(u).

## (A) The parity law

W(−a) = W(a) (csc² is even about P/2) and X(−m) = X(m) (S_e(−m) = conj S_e(m),
signs real, so every |·|² is even). Hence cross(−w) = cross(w): **the frame-law
factors through U_P/{±1}**. Substituting w ↦ −w in ĉ gives ĉ(χ) = χ(−1)ĉ(χ), so

> **every odd character (χ(−1) = −1) carries zero mass** — and both factors die
> separately: Ŵ_g(χ) = χ(−1)Ŵ_g(χ), X̂_g(χ) = χ(−1)X̂_g(χ). ∎

**Decode of the S60 vanishings.** The quadratic (Legendre) character mod p is odd
iff p ≡ 3 (mod 4). Since 7 ≡ 3 ≡ 3 (mod 4), the S60 zeros at mod 7 and mod 3 are
**parity zeros**; 5 ≡ 1 (mod 4), so mod 5 survives (−8.3403). The LEM-031
conjecture "the seven-section structure annihilates its own character" is
**CORRECTED**: the even conductor-7 characters — the two cubics — carry mass on
every referee cluster (balanced +1.974±1.287i; two-owner −2.640±1.889i; family
−5.042±0.344i). Parity is the whole story. Bonus symmetry: ĉ(χ̄) = conj ĉ(χ)
(cross is real), so the census needs only even characters up to conjugation.

## (B) The support law

The class-g term of cross(w) is a function of w mod q. Summing χ̄ over a fiber of
U_P → U_q (a coset vK of the kernel K) gives χ̄(v)Σ_{k∈K}χ̄(k) = 0 unless χ is
trivial on K, i.e. **cond(χ) | q = P/g**. ∎ (Referee: six spot instances at
10⁻¹⁶–10⁻¹⁸.) So ĉ(χ) = Σ_{g : cond(χ) | P/g} Ŵ_g(χ)X̂_g(χ) exactly.

## (C) The twisted Jordan lemma

For χ a character mod q, T_q(χ) := Σ_{u∈(Z/q)*} χ(u)csc²(πu/q):

| χ | T_q(χ) |
|---|---|
| odd | **0** |
| trivial | **J₂(q)/3** (THM-892(C\*) is the trivial row) |
| even | **(2q²/π²)·L_q(2,χ)** (L-series with Euler factors at p \| q removed) |
| even primitive | **2·τ(χ)·B_{2,χ̄}** (Gauss sum × generalized Bernoulli) |

*Proof.* (i) Kernel inversion (THM-892(K)): csc²(πn/q) = −2qΣ_j K_j e(−nj/q),
K_j = (j/q)(1−j/q). (ii) For primitive χ: Σ_u χ(u)e(−uj/q) = χ̄(−j)τ(χ).
(iii) Quadratic Bernoulli moment: Σ_j χ̄(j)K_j = −B_{2,χ̄}/q for χ nonprincipal
(B₂(x) = x²−x+1/6, B_{2,χ̄} = qΣ_j χ̄(j)B₂(j/q)). Assembling:
T_q(χ) = 2χ(−1)τ(χ)B_{2,χ̄}. (iv) Partial fractions csc²(πx) = π⁻²Σ_{n∈Z}(x+n)⁻²
give T_q(χ) = (1+χ(−1))(q²/π²)L_q(2,χ); imprimitive moduli inherit by Euler
restriction. The trivial row recovers J₂(q)/3 via ζ(2)Π(1−p⁻²). ∎
*Corollary (self-contained rederivation of the classical evaluation).* For even
primitive χ mod q: L(2,χ) = π²τ(χ)B_{2,χ̄}/q².

Exact-rational faces from the referee: T₅(quad) = 8/√5 with B_{2,χ₅} = 4/5;
T₁₃(quad) = 8√13 with B_{2,χ₁₃} = 4 exactly.

## (D) The closed-form weight side

Ŵ_g(χ) = Σ_{u∈U_q} W(gu)χ(u) = (π²/P²)T_q(χ), so by (C):

> **Ŵ_g(χ) = (2/g²) · L_{P/g}(2, χ)** — the g-dependence is pure 1/g²; the
> arithmetic content is a Dirichlet L-value at s = 2.

The factorization law (LEM-031) in final form:

> **ĉ(χ) = Σ_{g | P, g < P, cond(χ) | P/g} (2/g²) · L_{P/g}(2,χ) · X̂_g(χ).**

Referee: the full spectrum reproduced from pure L-values (no csc² sums) at worst
rel. err 7·10⁻¹⁴ over the top-mass + named characters of all three clusters,
including the S60 regression ĉ(Legendre mod 5) = −8.3403 on the two-owner
cluster. **The weight factor is now closed form for every (g, χ); the only
cluster-dependent object left in the frame spectrum is X̂_g(χ) — the χ-twisted
signed coincidence sums (the N(h)/Ĉ family of THM-892/THM-879).**

## (E) The census (measured)

All characters of (Z/P)\*, three clusters; Parseval Σ_{χ≠χ₀}|ĉ(χ)|² = Var_w(cross)
exact (3202.7 / 3563.1 / 4060.3). Masses are spread — top character carries only
2.4% / 3.2% / 10% of the variance — consistent with THM-892's generic-frame
picture. **Measured pattern (not proved): the co-resonant conductors contain the
full 7-part of P** — balanced (P = 2940, 7-part 49): top conductors 735, 245,
147, 980, 49 all ≡ 0 mod 49; two-owner (7-part 7): 315, 35, 1260; family: 420,
140, 35, 105. The seven-section skeleton surfaces by CONCENTRATING mass on
conductors containing its prime, not by killing its own conductor. Named open:
prove the 7-part concentration from the comb structure of X.

## Evidence log
- [x] Part C battery: 144 characters, 10 moduli, four faces, worst 9.3e-15
- [x] Parity at machine zero: 528 odd characters across three clusters
- [x] Support law: 6 instances ≤ 2e-16; Parseval exact ×3; conjugation symmetry ×3
- [x] Closed-form spectrum: worst 6.9e-14; S60 regression reproduced
- [x] the 7-part concentration law — PROVED as LEM-033 (S62): valuation–conductor pairing, per prime, exact
- [x] X̂_g(χ) closed form — LEM-033(3) (S62): cross-pair Gauss expansion over compatible grade cells
