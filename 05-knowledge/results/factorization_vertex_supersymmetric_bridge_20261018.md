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

**Typing.**
- PROVED: THM-4669.
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

Commutative factorization algebras correspond to commutative `D_X`-algebras. For a higher-dimensional `X` (or a
derived stack), factorization algebras on `Ran(X)` are defined the same way. The chiral side becomes chiral Koszul
duality (Francis-Gaitsgory) and factorization homology (Lurie, *Higher Algebra* §5.5; Ayala-Francis).

### 1.4 Vertex algebras as translation-equivariant chiral algebras

- A `G_a`-equivariant chiral algebra on `A^1` is the same as a vertex algebra `(V, |0>, T, Y(-, z))`: the fibre at 0
  is `V`, and the polar part of the chiral bracket is the operator product expansion (BD §0.15; Frenkel-Ben-Zvi,
  *Vertex algebras and algebraic curves*, 2nd ed., ch. 19-20).
- Conversely, a quasi-conformal vertex algebra `V` defines a chiral algebra `V_X` on every smooth curve `X`. The
  construction twists `V` by the `Aut O`-torsor of formal coordinates.

### 1.5 Costello-Gwilliam factorization algebras and holomorphic field theories

- A **prefactorization algebra** on a manifold `M` assigns a cochain complex `F(U)` to every open `U`. Disjoint
  inclusions `U_1 ⊔ ... ⊔ U_k ⊂ V` give multiplications, compatibly.
- BV quantization of a perturbative field theory produces one (Costello-Gwilliam, *Factorization algebras in quantum
  field theory*, vols. 1-2, CUP 2017/2021).
- On `M = C`, a holomorphically translation-invariant, `S^1`-equivariant prefactorization algebra satisfying finiteness
  conditions has a vertex algebra structure on its cohomology (CG vol. 1, Thm 5.2.3.1).
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

### 2.4 Twining genera derived from moonshine congruences

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
| **23A/B** | 1 | `(23/12) φ_2^(23) - (23/11)[2 η(τ)^2 η(23τ)^2 + η(τ)η(23τ) θ_(x^2+xy+6y^2)]` | -1, 1, -1/2, 0, 0, ..., 1 | unique solution in the scan; equals Cheng's (her basis is two Hecke eigenforms; the `q` and `q^2` coefficients agree in the 2-dimensional `S_2(Γ_0(23))`) |

Here `φ_2^(N) = (N E_2(Nτ) - E_2(τ))/(N - 1)`.

**Method.** The congruence `Tr(g | K_n) ≡ dim K_n (mod p)` determines the cusp-form coefficients without any
character table. This reproduces the literature (Cheng 2010; Gaberdiel-Hohenegger-Volpato 2010; Eguchi-Hikami 2010)
from moonshine structure alone.

### 2.5 From the Paley tournament `P_23` to M24

**The Golay code.** The quadratic residues mod 23 are the out-neighbourhood of 0 in the Paley tournament `P_23`. The
translates `a + ({0} ∪ QR)`, together with the all-ones word, span a 12-dimensional binary code on `P^1(F_23)`. Its
weight enumerator is `1 + 759 y^8 + 2576 y^12 + 759 y^16 + y^24`: it is the Golay code.

**M24.** Conway's generators `x+1`, `2x`, `-1/x` and `δ` preserve the code. In a census of 400,000 group elements (long
random words), the cycle-shape frequencies match `1/|C(g)|` for all 26 classes, and no other shape occurs. The census
also fixes the labels:
- `2^4 4^4` is class 4A (frequency 0.00264 against 1/384);
- `1^4 2^2 4^4` is class 4B;
- `4^6` is class 4C.

**OBSERVATION: the Paley automorphisms twine the K3 elliptic genus by 11A and 23A/B.** `Aut(P_23)` is the group
`{x -> ax + b : a ∈ QR}` of order 253. Inside M24 it consists of
- the identity;
- 22 elements of class 23A/B (the translations);
- 230 elements of class 11A.

So the symmetries of the Paley tournament act on the K3 sigma model with exactly the twining genera of §2.4. These
contain:
- the newform `η(τ)^2 η(11τ)^2` of the elliptic curve 11a;
- the level-23 cusp forms built from Ramanujan's weight-1 form `η(τ)η(23τ)`.

That weight-1 form encodes the Hilbert class field of `Q(sqrt(-23))`, with class number `h(-23) = 3`. By THM-4666 the
same class number gives `P_23`'s natural-order ranking excess, `23 · 3 = 69`.

## 3. Chain 2: Mason eta-products, Apéry-like congruences, CHL BPS counting

Scripts: m24_paley_golay_eta_20261018.py, apery_like_m24_supercongruences_20261018.py,
apery_like_halfperiod_clusters_20261018.py, s7_apery_*_20261018.py, swap_lemma_halfperiod_20261018.py.

### 3.1 The 26 eta-products (FINITE-EXACT)

For each class with cycle shape `Π a^(n_a)`, `η_g = Π η(aτ)^(n_a)` has weight `#cycles/2` and level
`(max a)(min a)`. All 26 are Hecke eigenforms; this was checked to `q^600` (Mason 1985; Dummit-Kisilevsky-McKay 1985).

| weight | classes | level and eigenform |
|---|---|---|
| 4 | 3B, 4A, 5A, 6A | levels 9, 8, 5, 6 |
| 3 | 4C, 7A/B, 8A | levels 16, 7, 8 |
| 2 | 6B, 10A, 11A, 12A, 14A/B, 15A/B | elliptic curves 36a, 20a, 11a, 24a, 14a, 15a |
| 1 | 12B, 21A/B, 23A/B | levels 144, 63, 23 |

### 3.2 Half-period supercongruences (FINITE-EXACT; THM-4669 PROVED)

Test: `u((p-1)/2) ≡ (d/p) a_g(p) (mod p^2)` for all primes `5 <= p <= 900`, with 17 sequences, 26 classes and 12
twists.

| sequence | modulus | class and twist | status |
|---|---|---|---|
| Apéry ζ(3) numbers `A(n)` | `p^2` | 4A, `(1/p)` | Ahlgren-Ono 2000, recovered |
| Apéry ζ(2) numbers `b(n)` | `p^2` | 4C | Stienstra-Beukers, recovered |
| **Cooper's `s_7`** | `p^2` | **4A, `(-1/p)`** | **THM-4669**, proved below |
| Domb; Franel-4 (`= s_10`) | `p` only | 4A, `(-1/p)`; 4A | |
| Franel (`= Zagier A`); Zagier D | `p` only | 8A, `(-1/p)` or `(2/p)`; 8A | |
| Zagier B | `p` only | 4C | |
| controls `C(2n,n)^2`, `C(2n,n)^3` | none | | |

Only the 2-power classes 4A, 4C and 8A ever occur.

**Classes.** Pairwise clustering of the sequences finds one class modulo `p^2` beyond the known mechanisms:
`{A, s_7 · (-1/p), the truncated 4F3 Σ C(2k,k)^4/256^k}`. Modulo `p`, every ζ(3)-type sequence falls into the 4A
class.

**THM-4669.** For every odd prime `p`, `s_7((p-1)/2) ≡ (-1)^((p-1)/2) A((p-1)/2) (mod p^2)`.
- Hence `s_7((p-1)/2) ≡ (-1/p) a_(4A)(p)`.
- Equivalently, `s_7((p-1)/2) ≡ (-1)^((p-1)/2) [x^((p-1)/2)] (2K(sqrt x)/π)^2`.

Proof outline:
- exact: `C(n,k)C(2k,n) = C(2k,k)C(k,n-k)`;
- exact: Clausen's formula `3F2(1/2,1/2,1/2;1,1;4x(1-x)) = 2F1(1/2,1/2;1;x)^2`;
- a `p`-adic perturbation `-1/2 = n - p/2`, together with the exact antisymmetry of `H_k - H_(n-k)`.

Every step was checked for all odd `p < 400`, and the theorem for all 237 primes below 1500.
- **Specificity.** The congruence fails modulo `p^3`. The swap `C(n+k,k) <-> C(2k,n)` fails for nine other weights.
- **Frobenius points.** The congruence holds at every `n = (p^r - 1)/2`. HYP-9291 was filed after 13 numerical cases
  with `r = 2, 3` and resolved by the same proof with `p`-adic valuation bounds; THM-4669 now covers every `r >= 1`.

**THM-4670, the ζ(2)-class twin.** An exhaustive search over three- and four-factor binomial summands
(clausen_twin_search_20261018.py) found exactly one more half-period twin of the Apéry numbers. Every other hit is a
rewriting, or `s_7` in disguise.
- **Statement.** For every prime `p >= 5`, `Σ_(k<=(p-1)/2) C(n+k,k)^3 ≡ b((p-1)/2) (mod p^2)`, with `b` Apéry's ζ(2)
  numbers, governed by `η(4τ)^6`, class 4C.
- **Proof.** The congruence reduces to `Σ a_k^3 O_k ≡ 0 (mod p)` and is proved by differentiating Dixon's well-poised
  summation at `a = -n`:
  - odd `n`: `Σ (-1)^k C(n,k)^3 H_(n-k) = Γ(1-n/2)Γ(1+3n/2)/(3n Γ(1+n/2)^2)`, which is divisible by `p`;
  - `n = 2m`: the defect is `(D_n/2)(H_m - H_(3m))`, and `H_(3m) ≡ H_m (mod 4m+1)` by reflection.

Each of the two Apéry classes at the half-period (4A for ζ(3), 4C for ζ(2)) thus has a proven Clausen twin.

### 3.3 Supersymmetric meaning of the eta-products (CITED + OBSERVATION)

**Twisted Fock spaces.** `1/η_g(τ)` is the trace of `g q^(L_0 - 1)` on the Fock space of 24 chiral bosons: the rank-24
Heisenberg vertex algebra, permuted by `g` (Mason).

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
algebra. Each fact is classical; the composition was not found stated.

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
- **THM-4669** (PROVED at every Frobenius point `(p^r - 1)/2`), with the Clausen-identity mechanism.
- **THM-4670** (PROVED), the ζ(2)-class twin, via the derivative of Dixon's summation.
- The classification of Apéry-like half-period congruences against the M24 eta-products (FINITE-EXACT).
- The congruence method that derives twining genera.
- The Paley `P_23` → 11A/23A twining observation.
- The Apéry ζ(2) ↔ Lee-Yang / `(A_1, A_2)` composition.
- Independent numerical verifications of the LG/K3, Mathieu, SU(2) `N_f = 4` and Lee-Yang statements.

**Credited.** All of section 1, and the literature cited in sections 2-4. The repository's earlier moonshine chain
(collatz_procgen_20260922_synthesis.md, "6 = 5+1 → S_5 → M_12 → M_24 → Golay → Leech") and Ogg's primes
(special_primes_curling_genus_zero_langlands_20261015.md) are the prior repository context.

**Open.**
- **Higher powers.** HYP-9291 (the Frobenius points) is resolved within THM-4669. Open: a corrected form modulo
  `p^3`; the plain congruence fails there.
- **A conceptual reason for THM-4669.** Both supercongruences point to the weight-4 motive 8.4. An explanation via the
  K3 families of the levels 6 and 7 is open.
- **A factorization-algebra meaning of the M24 twining.** The trace of a symmetry of a vertex algebra is a
  factorization-algebra trace. The twining of §2.5 is the trace of `Aut(P_23)` acting on the chiral algebra of the K3
  sigma model, through the Golay code.
