# Factorization algebras on derived stacks, vertex algebras via Beilinson-Drinfeld, and supersymmetric QFT: definitions, and three computed chains from the repository's objects

2026-10-18, opus-2026-10-18-S43.

**Owner prompt.**
> Keep going on strongest open leads as they emerge. If you see an opportunity for a meaningful proof or hypothesis,
> run many experiments to investigate all around it, related concepts and such. Maybe define factorization algebras on
> derived stacks, connect them to vertex algebras via the Beilinson-Drinfeld approach, and eventually apply the whole
> machinery to supersymmetric QFTs.

**What this note does.**
- Section 1 gives the definitions, precisely and with sources. These are standard (CITED); nothing in section 1 is
  claimed as new.
- Sections 2-4 are three chains, each computed from scratch in this session. Each starts from an object already in
  the repository and ends in a supersymmetric QFT:
  - Section 2: a Landau-Ginzburg derived critical locus, the K3 elliptic genus, Mathieu moonshine, M24, and back to the
    Golay code built from the Paley tournament `P_23`.
  - Section 3: Mason's eta-products of M24, Apéry-like supercongruences and CHL string BPS counting. This chain
    produced the session's theorem, **THM-4669**.
  - Section 4: 4d `N = 2` SCFTs and their vertex algebras, and the Lee-Yang / Argyres-Douglas parametrisation of
    Apéry's ζ(2) numbers.

**Audit.** Independently audited in S43; the report is 04-computation/experiments/s43_audit_20261018/REPORT.md, and
the findings are logged as MISTAKE-629.
- MAJOR M1: THM-4670 is a corollary of the classical supercongruences.
- MAJOR M2: the K3 sigma-model claim about Paley automorphisms is retracted.
- MAJOR M3: the twining genera are VERIFIED against Cheng, not derived.
- MAJOR M4: a false battery sentence is corrected.
- 12 minor findings, all applied. Every numerical claim was reproduced by independent code.

**Typing.**
- PROVED: THM-4669 (novelty unverified); THM-4670 (an elementary proof of a classical equivalence); THM-4685 (the twin
  classes; see the separate note).
- OPEN: HYP-9310 (the CM harmonic log law).
- VERIFIED: numerical or exact checks of known results.
- FINITE-EXACT: the censuses and batteries.
- CITED: the literature, including section 1.
- OBSERVATION: a composition of known facts that was not found stated.
- ANALOGY: section 5.

## 1. Definitions (CITED)

### 1.1 Derived stacks and derived critical loci

Work over `C`. A **derived stack** is a functor `X : cdga_(<=0) -> ∞Grpd` from connective commutative dg algebras (or
simplicial commutative rings) to ∞-groupoids satisfying étale hyperdescent (Toën-Vezzosi, HAG II, Mem. AMS 902, 2008;
Lurie, DAG; Gaitsgory-Rozenblyum, *A study in derived algebraic geometry*, 2017).

The example used below is the **derived critical locus** of a function `W : X -> A^1` on a smooth scheme. It is the
derived zero locus of `dW`, a section of `T*X`. Its algebra of functions is the Koszul complex
```
O(Crit^h W) = ( Sym_(O_X)(T_X[1]),  ι_(dW) ),      ι_(dW)(∂_i) = ∂_i W.
```
For an isolated singularity the cohomology is concentrated in degree 0 and equals the Jacobian ring
`Jac(W) = O_X/(∂_1 W, ..., ∂_n W)`. For `W = x^n` this is `C[x]/(x^(n-1))`, of dimension `n - 1`. This dimension is
the Witten index of the Landau-Ginzburg model with superpotential `W` (section 2.1).

### 1.2 The Ran space

For a scheme `X`, `Ran(X)` is the prestack sending a test scheme `S` to the set of nonempty finite subsets of
`Hom(S, X)`. Equivalently, `Ran(X) = colim_(I, surjections) X^I` over nonempty finite sets.

Two subsets are *disjoint* when they have disjoint graphs, which defines the open `(Ran × Ran)_disj`. Union gives the
map `∪ : Ran × Ran -> Ran` (Beilinson-Drinfeld, *Chiral algebras*, AMS Colloq. Publ. 51, 2004, §3.4; Francis-
Gaitsgory, *Chiral Koszul duality*, Selecta 18, 2012).

### 1.3 Factorization algebras and chiral algebras

A **factorization algebra** on a curve `X` is a D-module (or ind-coherent sheaf) `A` on `Ran(X)` together with
isomorphisms
```
(A ⊠ A)|_disj  ≅  (∪^! A)|_disj,
```
a unit, and the associativity and commutativity coherences. Restricting to `X^I` gives the familiar compatible system
`A_(X^I)` (BD §3.4).

A **chiral algebra** is a Lie algebra in D-modules on `X` for the chiral tensor structure. It is a D-module map
`j_* j^*(A ⊠ A) -> Δ_! A` satisfying antisymmetry and the Jacobi identity, together with a unit.

**Theorem (BD 3.4.9; Francis-Gaitsgory).** On a smooth curve, factorization algebras and chiral algebras are
equivalent.

Commutative factorization algebras correspond to commutative `D_X`-algebras.
- **Higher dimensions.** For a higher-dimensional `X`, factorization algebras on `Ran(X)` are defined the same way.
  Francis-Gaitsgory prove the factorization = chiral equivalence in every dimension (arXiv:1911.01627, Theorem 2.15(3)
  and Remark 2.16(3), which confirm BD 3.4.9).
  - *Chiral Koszul duality* is a different statement. It relates chiral Lie algebras and factorization coalgebras.
  - See also factorization homology (Lurie, *Higher Algebra* §5.5; Ayala-Francis).
  - (Corrected after S43 audit m12.)
- **Derived stacks.** D-modules (crystals) do not see derived structure, so on a derived stack one must use
  ind-coherent sheaves (IndCoh) on `Ran(X)` (Gaitsgory-Rozenblyum).

### 1.4 Vertex algebras as translation-equivariant chiral algebras

- A `G_a`-equivariant chiral algebra on `A^1` is the same as a vertex algebra `(V, |0>, T, Y(-, z))`: the fibre at 0
  is `V`, and the polar part of the chiral bracket is the operator product expansion (Frenkel-Ben-Zvi, *Vertex
  algebras and algebraic curves*, 2nd ed., §19.2; Huang-Lepowsky; explicit functors in the appendix of
  Bakalov-De Sole-Heluani-Kac 2019). The earlier pointer "BD §0.15" could not be verified (audit m12).
- Conversely, a quasi-conformal vertex algebra `V` defines a chiral algebra `V_X` on every smooth curve `X`. The
  construction twists `V` by the `Aut O`-torsor of formal coordinates.

### 1.5 Costello-Gwilliam factorization algebras and holomorphic field theories

- A **prefactorization algebra** on a manifold `M` assigns a cochain complex `F(U)` to every open `U`. Disjoint
  inclusions `U_1 ⊔ ... ⊔ U_k ⊂ V` give multiplications, compatibly.
- BV quantization of a perturbative field theory produces one (Costello-Gwilliam, *Factorization algebras in quantum
  field theory*, vols. 1-2, CUP 2017/2021).
- On `M = C`, a holomorphically translation-invariant, `S^1`-equivariant prefactorization algebra satisfying finiteness
  conditions has a vertex algebra structure on its cohomology (CG vol. 1, Theorem 5.3.3; corrected from "5.2.3.1",
  audit m6).
- **Examples.**
  - The βγ system gives the βγ vertex algebra.
  - BV quantization of the curved βγ system with target a complex manifold `Y` gives the sheaf of chiral differential
    operators (Gorbounov-Gwilliam-Williams, Astérisque 419, 2020).
  - The chiral de Rham complex of `Y` is a sheaf of vertex algebras (Malikov-Schechtman-Vaintrob, CMP 204, 1999).
  - Costello's construction of the Witten genus uses exactly this: the derived mapping stack from an elliptic curve
    into `Y` (*A geometric construction of the Witten genus I/II*, 2010/2011).

### 1.6 Supersymmetric QFT: twists, chiral algebras and indices

- **Twists.** A nilpotent supercharge `Q` (`Q^2 = 0`) makes the theory topological or holomorphic in some directions
  (Costello; Elliott-Safronov-Williams classify the twists).
- **2d theories.** For a 2d `N = (0,2)` or `(2,2)` theory, the `Q`-cohomology of the half-twisted theory is a vertex
  algebra. Its graded character with `(-1)^F` is the **elliptic genus** (Witten 1987; Witten 2005; Kapustin 2005).
  - For a Landau-Ginzburg model with quasi-homogeneous `W`, the twisted theory is the βγbc system with the differential
    of `W`.
  - Its elliptic genus depends only on the R-charges, `Z = Π_i θ_1(τ, (1-w_i) z)/θ_1(τ, w_i z)` (Witten 1993; Kawai-
    Yamada-Yang, NPB 414, 1994, including orbifolds).
- **4d `N = 2` SCFTs** (Beem-Lemos-Liendo-Peelaers-Rastelli-van Rees, CMP 336, 2015).
  - The cohomology of a nilpotent combination of Poincaré and conformal supercharges, restricted to a plane, is a vertex
    algebra with `c_2d = -12 c_4d`.
  - The Schur index equals its vacuum character.
  - Characters satisfy modular linear differential equations (MLDEs) (Beem-Rastelli, JHEP 2018; Arakawa-Kawasetsu
    2018).

### 1.7 The chain, schematically

```
derived stack (Crit^h W, Map(E, Y), ...)  --BV quantization-->  holomorphic factorization algebra on Ran(Σ)
   --(BD 3.4.9)-->  chiral algebra on Σ  --(fibre, G_a-equivariance)-->  vertex algebra V
   --(trace of q^(L_0) (-1)^F y^(J_0))-->  elliptic genus / Schur index  -->  modular and Jacobi forms.
```

Sections 2-4 compute the last arrows explicitly in three cases. They also identify where repository objects enter.

## 2. Chain 1: Landau-Ginzburg, K3, Mathieu moonshine, M24, Golay, Paley `P_23` (VERIFIED, FINITE-EXACT)

Scripts: k3_lg_elliptic_genus_mathieu_20261018.py, m24_paley_golay_eta_20261018.py.

### 2.1 Derived critical locus and Witten index

For `W = x^n` the elliptic genus is `Z(τ, z) = θ_1(τ, (1 - 1/n) z)/θ_1(τ, z/n)`. Numerically `Z(τ, 0) = n - 1` for
`n = 3, 4, 5, 6`. This equals `dim Jac(W)`, the Euler characteristic of the derived critical locus (§1.1).

### 2.2 The quartic orbifold is the K3 elliptic genus

The Kawai-Yamada-Yang orbifold of `W = x_1^4 + x_2^4 + x_3^4 + x_4^4` by `Z_4` is
```
Z_orb = (1/4) Σ_(a,b=0)^3 q^(a^2) y^(2a) Π_i Z_(1/4)(τ, z + aτ + b)       (ĉ = 2).
```
It equals the K3 elliptic genus `2 φ_(0,1) = 8 Σ_(i=2,3,4) (θ_i(τ,z)/θ_i(τ,0))^2` to relative error
`1.8 × 10^(-15)` at six random points.

This is the Gepner `(k=2)^4` model: an `N = (2,2)`, in fact `(4,4)`, SCFT on the quartic K3. Its holomorphic twist is
the chain of §1.7 starting from the derived critical locus of `W`, orbifolded.

### 2.3 The `N = 4` decomposition: Mathieu multiplicities

Eguchi-Hikami's formula `H(τ) = (-2 E_2 + 48 Σ_(r>s>0, r-s odd) (-1)^r s q^(rs/2))/η^3` gives
```
H = 2 q^(-1/8) (-1 + 45 q + 231 q^2 + 770 q^3 + 2277 q^4 + 5796 q^5 + 13915 q^6 + 30843 q^7 + 65550 q^8
                + 132825 q^9 + 260568 q^10 + 494385 q^11 + 915124 q^12 + ...),
```
the Mathieu moonshine sequence (Eguchi-Ooguri-Tachikawa 2010).

### 2.4 Twining genera: moonshine congruences, checked against Cheng 2010 (VERIFIED)

For `g ∈ M23` of prime order `p` with `χ_g = 24/(p+1)` fixed points:
```
Z_g = (χ_g/12) φ_(0,1) + T_g φ_(-2,1),     H_g = (χ_g/24) H - T_g/η^3,     T_g(i∞) = 2 - χ_g/12.
```
The traces `Tr(g | K_n) = 2 A_n(g)` must be integers congruent to `dim K_n` modulo `p`.

| class | fixed points | weight-2 part `T_g` | `A_1(g), A_2(g), ...` | how determined |
|---|---|---|---|---|
| 2A | 8 | `(4/3) φ_2^(2)` | -3, 7, -14, 21, -28, 43, ... | forced: `M_2(Γ_0(2))` is 1-dimensional |
| 3A | 6 | `(3/2) φ_2^(3)` | 0, -3, 5, 0, -9, 10, ... | forced |
| 5A | 4 | `(5/3) φ_2^(5)` | 0, 1, 0, -3, 1, 0, ... | forced |
| 7A/B | 3 | `(7/4) φ_2^(7)` | -1/2, 0, 0, 2, 0, -1, 1, ... | forced; traces -1, 0, 0, 4, 0, -2, 2 equal Cheng's Table 1 |
| **11A** | 2 | `(11/6) φ_2^(11) - (22/5) η(τ)^2 η(11τ)^2` | 1, 0, 0, 0, -1, 0, -1, 1, ... | `c ≡ -22/5 (mod 11)` from the congruences; the minimal solution equals Cheng's formula |
| **23A/B** | 1 | `(23/12) φ_2^(23) - (23/11)[2 η(τ)^2 η(23τ)^2 + η(τ)η(23τ) θ_(x^2+xy+6y^2)]` | -1, 1, -1/2, 0, 0, ..., 1 | a solution of the congruences, not unique (see below); equals Cheng's combination coefficient by coefficient. Cheng uses the rational basis `f_(23,1)`, `f_(23,2)` of sums and differences of the two eigenforms, with `f_(23,2) = η(τ)^2 η(23τ)^2`, and `dim S_2(Γ_0(23)) = 2` |

Here `φ_2^(N) = (N E_2(Nτ) - E_2(τ))/(N - 1)`.

**Method and its limit** (corrected after S43 audit M3). Integrality and `Tr(g | K_n) ≡ dim K_n (mod p)` determine
the cusp-form coefficients only modulo `p`. The forms involved have integral `q`-coefficients, so the congruences
leave lattices of solutions:
- 11A: every `c ∈ -22/5 + 11Z`. The session's own output lists `-22/5` and `33/5`.
- 23A/B: every `(c_1, c_2) ∈ (-46/11, -23/11) + 23Z^2`. The earlier "unique solution in the scan" was an artefact of
  the scan window.

The table entries are the minimal representatives. They **coincide** with the literature (Cheng 2010;
Gaberdiel-Hohenegger-Volpato 2010; Eguchi-Hikami 2010), which is what identifies them: the table is VERIFIED against
those sources, not derived independently of them.
- The audit confirmed 2A, 3A and 7A against Cheng's Table 1 for `q^1..q^7`, and 5A except one misprint: Cheng's printed
  `q^7` coefficient `2 × 7` must be `2 × 3`, as forced by `M_2(Γ_0(5))`.
- It confirmed the 11A and 23A/B traces integral and `≡ dim (mod p)` for `n <= 60`.

### 2.5 From the Paley tournament `P_23` to M24

**The Golay code.** The quadratic residues mod 23 are the out-neighbourhood of 0 in the Paley tournament `P_23`. The
translates `a + ({0} ∪ QR)`, together with the all-ones word, span a 12-dimensional binary code on `P^1(F_23)`. Its
weight enumerator is `1 + 759 y^8 + 2576 y^12 + 759 y^16 + y^24`: it is the Golay code.

**M24.** Conway's generators `x+1`, `2x`, `-1/x` and `δ` preserve the code.
- Schreier-Sims (audit) gives `|<α, β, γ, δ>| = 244823040 = |M24|` exactly.
- In a census of 400,000 group elements (long random words), the cycle-shape frequencies match `1/|C(g)|` for all 21
  cycle shapes, and no other shape occurs. (There are 26 classes; five A/B pairs share a shape.)

The census also fixes the labels:
- `2^4 4^4` is class 4A (frequency 0.00264 against 1/384);
- `1^4 2^2 4^4` is class 4B;
- `4^6` is class 4C.

**OBSERVATION: the Paley automorphisms lie in the M24 classes 11A and 23A/B.** `Aut(P_23)` is the group
`{x -> ax + b : a ∈ QR}` of order 253. Inside M24 it consists of
- the identity;
- 22 elements of class 23A/B (the translations);
- 230 elements of class 11A.

The corresponding Mathieu-moonshine twining functions are those of §2.4.

**Retracted** (S43 audit M2). An earlier version said that these symmetries "act on the K3 sigma model". That is
false.
- Symmetry groups of K3 sigma models are subgroups of `Co_0` fixing a rank `>= 4` sublattice of the Leech lattice
  (Gaberdiel-Hohenegger-Volpato 2011).
- The twining genera of 12B, 21A/B and 23A/B cannot be realised by any symmetry of a K3 sigma model
  (Gaberdiel-Volpato, arXiv:1206.5143).
- `Aut(P_23)` contains elements of class 23A/B, so it acts on no K3 sigma model.
- 11A is realised, but only non-geometrically (Mukai: symplectic automorphisms have order `<= 8`).
- The 23A/B functions are McKay-Thompson series of Mathieu moonshine, not traces of a sigma-model symmetry.

These twining functions contain:
- the newform `η(τ)^2 η(11τ)^2` of the elliptic curve 11a;
- the level-23 cusp forms built from Ramanujan's weight-1 form `η(τ)η(23τ)`.

That weight-1 form encodes the Hilbert class field of `Q(sqrt(-23))`, with class number `h(-23) = 3`. By THM-4666 the
same class number gives `P_23`'s natural-order ranking excess, `23 · 3 = 69`.

## 3. Chain 2: Mason eta-products, Apéry-like congruences, CHL BPS counting

Scripts: m24_paley_golay_eta_20261018.py, apery_like_m24_supercongruences_20261018.py,
apery_like_halfperiod_clusters_20261018.py, s7_apery_*_20261018.py, swap_lemma_halfperiod_20261018.py.

### 3.1 The 21 eta-products of the 26 classes (FINITE-EXACT + CITED)

For each class with cycle shape `Π a^(n_a)`, `η_g = Π η(aτ)^(n_a)` has weight `#cycles/2` and level
`(max a)(min a)`.
- **Eigenforms.** All 21 distinct eta-products are Hecke eigenforms. This is CITED (Mason 1985;
  Dummit-Kisilevsky-McKay 1985).
- **Spot checks.** Multiplicativity for coprime `m, n < 40` with `mn <= 600`, and the Hecke relation at `p <= 23` with
  the character fitted as `±1`.
- (Corrected after audit m4: the earlier "checked to `q^600`" overstated the script.)
- **Elliptic curves.** The audit matched the weight-2 eta-products to 11a1, 14a1, 15a1, 20a1, 24a1 and 36a1 (point
  counts at all good `p < 200`).

| weight | classes | level and eigenform |
|---|---|---|
| 4 | 3B, 4A, 5A, 6A | levels 9, 8, 5, 6 |
| 3 | 4C, 7A/B, 8A | levels 16, 7, 8 |
| 2 | 6B, 10A, 11A, 12A, 14A/B, 15A/B | elliptic curves 36a, 20a, 11a, 24a, 14a, 15a |
| 1 | 12B, 21A/B, 23A/B | levels 144, 63, 23 |

### 3.2 Half-period supercongruences (FINITE-EXACT; THM-4669 PROVED)

Test: `u((p-1)/2) ≡ (d/p) a_g(p) (mod p^2)` for all primes `5 <= p <= 900`. The battery has 12 distinct Apéry-like
sequences plus 2 controls, the 21 distinct eta-products and 12 twists.

It was corrected after S43 audit m2:
- 'Zagier D' was Zagier's #E. Zagier's #D is Apéry's `b`.
- 'Cooper s18' was misdefined. It is now OSS 2016 eq. (7): 1, 6, 54, 564, 6390, ...
- 'Yang-Zudilin' (`Σ C(n,k)^2 C(2k,k)^2`, not Apéry-like) was removed. The Yang-Zudilin sequence is `Σ C(n,k)^4`,
  which is already in the battery.
- 'Zagier A' duplicated Franel.

The corrected battery was rerun at `p <= 300` and `p <= 900`, and the conclusions below survive.

| sequence | modulus | class and twist | status |
|---|---|---|---|
| Apéry ζ(3) numbers `A(n)` | `p^2` | 4A, `(1/p)` | Ahlgren-Ono 2000, recovered |
| Apéry ζ(2) numbers `b(n)` (= Zagier #D) | `p^2` | 4C | van Hamme 1986 / Ishikawa 1990 / Ahlgren 2001 (conjectured by Stienstra-Beukers 1985), recovered |
| **Cooper's `s_7`** | `p^2` | **4A, `(-1/p)`** | **THM-4669**, proved below |
| Domb; Franel-4 (`= s_10`) | `p` only | 4A, `(-1/p)`; 4A | |
| Franel (= Zagier #A); Zagier #E (`Σ C(n,k) C(2k,k) C(2n-2k,n-k)`, labelled "Zagier D" at filing) | `p` only | 8A, `(-1/p)` or `(2/p)`; 8A | |
| Almkvist-Zudilin; Cooper `s_18` (corrected); Zagier C, F | none | no M24 eta-product, even mod `p` | AZ and `s_18` agree with each other mod `p` at the half-period, twist `(-1/p)` |
| Zagier B | `p` only | 4C | |
| controls `C(2n,n)^2`, `C(2n,n)^3` | none | | |

Only the 2-power classes 4A, 4C and 8A ever occur.

**Classes.** Pairwise clustering of the sequences finds one class modulo `p^2` beyond the known mechanisms:
`{A, s_7 · (-1/p), the truncated 4F3 Σ C(2k,k)^4/256^k}`.

Modulo `p`, Apéry `A`, Domb, `s_10` and `s_7` fall into the 4A class. Almkvist-Zudilin, `Σ C(n,k)^2 C(2k,n)^2` and
`s_18` match no M24 eta-product. The earlier sentence "every ζ(3)-type sequence falls into the 4A class" was false
(audit M4).

**THM-4669.** For every odd prime `p`, `s_7((p-1)/2) ≡ (-1)^((p-1)/2) A((p-1)/2) (mod p^2)`.
- Hence `s_7((p-1)/2) ≡ (-1/p) a_(4A)(p)`.
- Equivalently, `s_7((p-1)/2) ≡ (-1)^((p-1)/2) [x^((p-1)/2)] (2K(sqrt x)/π)^2`.

Proof outline:
- exact: `C(n,k)C(2k,n) = C(2k,k)C(k,n-k)`;
- exact: Clausen's formula `3F2(1/2,1/2,1/2;1,1;4x(1-x)) = 2F1(1/2,1/2;1;x)^2`;
- a `p`-adic perturbation `-1/2 = n - p/2`, together with the exact antisymmetry of `H_k - H_(n-k)`.

Every step was checked for all odd `p < 400`, and the theorem for all 237 primes `5 <= p < 1500`. The audit
independently confirmed all 549 odd primes `p < 4000`.
- **Specificity.** At `r = 1` the congruence fails modulo `p^3`. At `r >= 2` it fails except at the non-ordinary prime
  `p = 11`. The swap `C(n+k,k) <-> C(2k,n)` fails for the nine other weights tested.
- **Frobenius points.** The congruence holds at every `n = (p^r - 1)/2`. HYP-9291 was filed after 13 numerical cases
  with `r = 2, 3` and resolved by the same proof with `p`-adic valuation bounds; THM-4669 now covers every `r >= 1`.

**THM-4670, the ζ(2)-class twin.** An exhaustive search over three- and four-factor binomial summands
(clausen_twin_search_20261018.py) found exactly one more half-period twin of the Apéry numbers. Every other hit is a
rewriting, or `s_7` in disguise.

**Audit M1.** This twin is an immediate corollary of classical results. Its Step 1 gives
`Σ C(n+k,k)^3 - b ≡ 4(Σ a_k^3 - b) (mod p^2)`, and both `b((p-1)/2) ≡ c(p)` and `Σ_(k<=(p-1)/2) a_k^3 ≡ c(p)` are the
classical guises of Beukers' ζ(2) supercongruence (van Hamme 1986, Ishikawa 1990, Ahlgren 2001; Mortenson 2008). What
THM-4670 adds is an elementary Dixon-derivative proof that the two guises agree.
- **Statement.** For every prime `p >= 5`, `Σ_(k<=(p-1)/2) C(n+k,k)^3 ≡ b((p-1)/2) (mod p^2)`, with `b` Apéry's ζ(2)
  numbers, governed by `η(4τ)^6`, class 4C.
- **Proof.** The congruence reduces to `Σ a_k^3 O_k ≡ 0 (mod p)` and is proved by differentiating Dixon's well-poised
  summation at `a = -n`:
  - odd `n`: `Σ (-1)^k C(n,k)^3 H_(n-k) = Γ(1-n/2)Γ(1+3n/2)/(3n Γ(1+n/2)^2)`, which is divisible by `p`;
  - `n = 2m`: the defect is `(D_n/2)(H_m - H_(3m))`, and `H_(3m) ≡ H_m (mod 4m+1)` by reflection.

Each of the two Apéry classes at the half-period (4A for ζ(3), 4C for ζ(2)) thus has a proven twin. The ζ(2) twin is
classical in substance. The name "Clausen twin" fits only THM-4669: THM-4670 uses Dixon, not Clausen (audit m11).

**All twin classes (THM-4685, HYP-9310).** A census of all 648 distinct three- and four-factor binomial sums finds
exactly seven half-period twin classes. They are governed by the M24 classes 4A, 4C and 7AB
([note](halfperiod_twin_classes_cm_log_law_20261018.md)).
- The three 4C classes are PROVED from the Dixon lemma and the reflection `k -> n-k`. Their values are `F`, `4^n F`
  and `16^n F`, with `F = Σ C(2k,k)^3/64^k`. Per audit M1, the Dixon lemma is equivalent to the classical guise
  agreement, so these are corollaries of the classical supercongruences plus elementary shadows.
- The 7AB classes reduce to `Σ_(k<=n) C(2k,k)^3 (H_(2k) - H_k) ≡ 0 (mod p)`. This is the case `m = 1` of a new
  pre-registered CM harmonic law: `6 Σ C(2k,k)^3 m^(-k) (H_(2k) - H_k) ≡ q_p(m) Σ C(2k,k)^3 m^(-k) (mod p)` at all
  eight of Sun's CM values.
  - The law is proved at the self-dual points `m = ±64` and is invariant under `m -> 4096/m`, the 2-isogeny.
  - It has no failures to `p < 20000` and fails at generic `m`.
- At `m = 16, 256` the governing form `η(2τ)^3 η(6τ)^3` is not an M24 cycle shape.

### 3.3 Supersymmetric meaning of the eta-products (CITED + OBSERVATION)

**Twisted Fock spaces.** `1/η_g(τ)` is the trace of `g q^(L_0 - 1)` on the Fock space of 24 chiral bosons: the rank-24
Heisenberg vertex algebra, permuted by `g` (Mason). `η_g` itself is a supertrace on a fermionic Fock space.

**CHL strings.** For `g` in M23, with cycle shapes `1^8 2^8, 1^6 3^6, 1^4 2^2 4^4, 1^4 5^4, 1^2 2^2 3^2 6^2, 1^3 7^3,
1^2 2 4 8^2`, `1/η_g` is the half-BPS (electric) counting function of the `Z_N` CHL orbifold of string theory on
`K3 × T^2` (Jatkar-Sen 2005; David-Jatkar-Sen 2006; Govindarajan-Krishna for `N = 4, 6, 8`).

Consequences:
- **OBSERVATION (Z_8 CHL).** The form governing the Franel numbers modulo `p`, `η_(8A)`, is the inverse half-BPS
  partition function of the `Z_8` CHL model.
- **No geometric CHL model for the Apéry forms.** The two Apéry forms `η_(4A)` and `η_(4C)` have no fixed points, so
  they are not in M23 and correspond to no geometric CHL orbifold. They are twining characters of the 24-boson Fock
  space, and of Mathieu moonshine for the non-geometric classes 4A and 4C.

## 4. Chain 3: 4d `N = 2` SCFTs and vertex algebras (VERIFIED)

Script: schur_index_vertex_algebra_mlde_20261018.py.

### 4.1 SU(2) with four flavours

From the Haar integral of the plethystic single-letter indices:
```
I(q) = 1 + 28 q + 329 q^2 + 2632 q^3 + 16380 q^4 + 85764 q^5 + 393589 q^6 + 1628548 q^7 + ...
```
The half-integer powers cancel. `f = q^(7/12) I` satisfies `D^(2) f = (35/144) E_4 f` at every computed order. This is
the vacuum character of `so(8)` at level `-2`, with `c = -14` (BLLPRR; Beem-Rastelli; Arakawa-Kawasetsu).

### 4.2 The `(A_1, A_2)` Argyres-Douglas theory and Lee-Yang

The Schur index is the Lee-Yang vacuum character `q^(11/60) H(q)`, with
`H = Π 1/((1-q^(5n+2))(1-q^(5n+3)))` (Cordova-Shao 2016). The `h = -1/5` character is `q^(-1/60) G(q)`. Both satisfy
the MLDE with `α = 11/3600` through `q^60` (Mathur-Mukhi-Sen).

### 4.3 OBSERVATION: Apéry's ζ(2) numbers in the language of the Argyres-Douglas chiral algebra

Beukers' Hauptmodul for the Apéry ζ(2) numbers is
```
t(τ) = q Π (1 - q^n)^(5 (n/5))  =  (χ_0 / χ_(-1/5))^5,
```
the fifth power of the ratio of the two Lee-Yang characters, i.e. of the Rogers-Ramanujan continued fraction. We
verified to `q^50` that
```
Σ_n b_n t^n = 1 + Σ_n Re[(3 - i) Σ_(d|n) χ(d)] q^n,
```
the weight-1 Eisenstein form on `Γ_1(5)` (`χ` the quartic character mod 5 with `χ(2) = i`). So Apéry's ζ(2)
generating function is a weight-1 modular form written in the ratio of the two characters of the `(A_1, A_2)` chiral
algebra.
- Each fact is classical: `t = r(τ)^5`, with `r` the Rogers-Ramanujan continued fraction, is in Beukers 1987 and in
  Zagier 2009, case D.
- The reading through the Argyres-Douglas characters was not found stated.

## 5. ANALOGY: the ROOT current as a Witten index

The S42 ROOT current is `R(N) = Σ_(N < k <= 4N) w(k) f(k)`, with `f(k) = (-1)^(stopping time)`.

| index side | ROOT side |
|---|---|
| `(-1)^F` grading | the clock parity `f` |
| `Q`-doublets with equal energy, cancelling in the index | sibling pairs `(s, 2s+1)`: opposite parity, the same transport image |
| the index localising on unpaired states | clock_duality's R5 removes `10N/21` in `l1`; what remains is singletons and boundary sources |
| localisation | the frontier lemma: Cauchy-Schwarz over first-entry basins |

No implication is claimed. The dictionary suggests asking for a "supercharge" on the free abelian group of sources
whose cohomology carries `R(N)`. That would make the square-root problem a question about the growth of a cohomology
group.

## 6. What is new, credited, open

**New.**
- **THM-4669** (PROVED at every Frobenius point `(p^r - 1)/2`), with the Clausen-identity mechanism. Novelty is
  unverified: no prior statement was found in the S43 or audit searches.
- **THM-4685** (PROVED): the census of seven twin classes, the three 4C classes, the reflection lemma, and the CM
  harmonic law at `m = ±64`.
- **HYP-9310** (OPEN, pre-registered and passed): the CM harmonic log law at the other six CM points.
- The classification of Apéry-like half-period congruences against the M24 eta-products (FINITE-EXACT, battery
  corrected after audit m2).
- The Paley `P_23` → Golay → M24 construction from scratch: the code, the generators, Schreier-Sims order and the
  census. Also `Aut(P_23) = 1 + 22 × 23A/B + 230 × 11A` inside M24 (group theory only; see the retraction in §2.5).
- The Apéry ζ(2) ↔ Lee-Yang / `(A_1, A_2)` composition.
- Independent numerical verifications of the LG/K3, Mathieu, SU(2) `N_f = 4` and Lee-Yang statements.

**Credited.**
- All of section 1, and the literature cited in sections 2-4.
- The twining genera of §2.4 (Cheng 2010; Gaberdiel-Hohenegger-Volpato 2010; Eguchi-Hikami 2010). The congruences fix
  them only modulo `p`; the minimal representatives coincide with these sources.
- **THM-4670** (PROVED): an elementary Dixon-derivative proof of the agreement of the two classical guises of Beukers'
  ζ(2) supercongruence (van Hamme / Ishikawa / Ahlgren). It is not a new congruence (audit M1).
- **Prior repository result on the same newform:** [THM-4652](../../01-canon/theorems/THM-4652-paley-pattern-sums-vanishing-and-series-parallel-polynomiality.md)
  / [HYP-9274](../hypotheses/HYP-9274-paley-rate-three-halves-wheel-w4-apery-form.md) (S34, the Paley "Apéry wheel"
  `Σ*(W_4) = -q(q-1)(a_q + 4)`). It is built on `η(2z)^4 η(4z)^4` with the same Ahlgren-Ono input (audit m13).
- **Concurrent work on the same owner prompt:** kind-pasteur-2026-10-18-S1 on the classical freeness of vertex
  algebras.
  - THM-4671: the first quantum correction of Vir_{p,q}.
  - THM-4672: the free hypermultiplet with a `Z_2` gauging.
  - THM-4673: the `c = 1`, `N = 2` minimal model.
  - HYP-9293..9297.
  - Their 4d chiral algebras overlap with §4 (SU(2) `N_f = 4`, Argyres-Douglas). The repository's earlier moonshine chain
(collatz_procgen_20260922_synthesis.md, "6 = 5+1 → S_5 → M_12 → M_24 → Golay → Leech") and Ogg's primes
(special_primes_curling_genus_zero_langlands_20261015.md) are the prior repository context.

**Open.**
- **Higher powers.** HYP-9291 (the Frobenius points) is resolved within THM-4669. Open: a corrected form modulo
  `p^3`; the plain congruence fails there.
- **A conceptual reason for THM-4669.** Both supercongruences point to the weight-4 motive 8.4. An explanation via the
  K3 families of the levels 6 and 7 is open.
- **A factorization-algebra meaning of the M24 twining.** The trace of a symmetry of a vertex algebra is a
  factorization-algebra trace. An earlier sentence read the twining of §2.5 as the trace of `Aut(P_23)` on "the chiral
  algebra of the K3 sigma model". That is retracted (audit M2): 23A/B is realised in no K3 sigma model. A
  factorization-algebra home for the Mathieu-moonshine module itself is open. Only 11A is realised in sigma models,
  and only non-geometrically.
