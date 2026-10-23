# Positive realization: the A1 lattice VOA is the semion bridge

Status: CITED standard construction and representation-category equivalence; elementary specializations/identities DERIVED; coefficient checks FINITE-EXACT; floating-point modular tests NUMERICAL controls only. No novelty claim. This addendum leaves the audited main note unchanged.

## 1. Primary-source gate: passed

Use the standard conformal vector, not a shifted Feigin–Fuchs vector.

- Allen, Lentner, Schweigert, Wood, *Duality structures for representation categories of vertex operator algebras and the Feigin–Fuchs boson*, Selecta Math. 31 (2025), 36. Proposition 3.8, equations (3.28)–(3.31), gives the lattice VOA, central charge, and module weights. Theorem 3.12 gives a ribbon equivalence for conformal parameter `ξ=0`; equations (3.11), (3.13) identify braiding trace and twist from the lattice pairing. Theorem 3.4 states the quadratic-form classification of the abelian cocycle class. These exact statements were inspected, not merely the abstract. [Primary full text](https://orca.cardiff.ac.uk/id/eprint/178289/1/s00029-025-01032-x..pdf).
- Höhn–Scheithauer, *A new generalized Kac-Moody algebra of rank 14*, §2.1, Theorem 2.1 and Lemma 2.2, records the discriminant-labelled irreducible modules, fusion, cocycle class, and theta/eta character formula. The displayed fusion line has a repeated lambda on its left in the retrieved text; the adjacent quantified statement and group-ring identification make the intended two-input formula unambiguous. [Primary author-hosted text](https://www.math.ksu.edu/~gerald/papers/gkm12v16.pdf).
- Dong–Li–Mason, *Regularity of rational vertex operator algebras*, proves regularity for positive-definite even lattice VOAs. [Primary paper](https://arxiv.org/abs/q-alg/9508018). This supplies the semisimplicity setting; no assertion about arbitrary vertex algebras is being used.

Dong's original 1993 DOI landing page was inaccessible in this check. No claim rests on pretending that it was inspected: the explicit primary papers above supply the needed statements.

## 2. Lattice, conformal weights, and the actual category

Let `L=Zα`, with `(α,α)=2`, so `L=√2 Z`. It is positive definite, even, and full rank. Set `ξ=0`. Then

`c=rank L=1`, `L*=Z(α/2)`, `D=L*/L=Z/2`.

There are two irreducible ordinary modules:

`M_0=V_L`, `M_1=V_(L+α/2)`.

Their ground-state weights follow directly from the momentum norm:

`h_0=min_(n∈Z) n²=0`,
`h_1=min_(n∈Z)(n+1/2)²=1/4`.

The nonvacuum ground space has dimension TWO, from momenta `±α/2`; having one simple module does not mean its ground-state space is one-dimensional.

Fusion is `M_a ⊠ M_b ≅ M_(a+b mod 2)`. The discriminant quadratic form is

`q_D(a)=exp(πi(aα/2,aα/2))=i^(a²)`, for `a=0,1`.

Consequently the standard ribbon twists are `(1,i)`. The explicit cocycle in the main note,

`ω(a,b,c)=(-1)^(abc)`, `r(a,b)=i^(ab)`,

has exactly this braiding trace `r(a,a)=q_D(a)`. The cited quadratic-form classification and lattice ribbon equivalence therefore identify its category, with the standard twist, with the module category of `V_(√2Z)`. This is an equivalence up to coherent gauge, not a claim that every arbitrary normalization of lattice intertwining operators produces those identical scalar entries.

This is the positive semion category, not anti-semion: reversing the braiding conjugates `i` to `-i`. The discriminant group alone would not select between them; the quadratic form does.

## 3. Characters and why the exponent is 1/24

Fix `q=exp(2πiτ)` with `Im τ>0`, and

`η(τ)=q^(1/24) Π_(m>=1)(1-q^m)`.

Use the normalized VOA character `χ_a(τ)=Tr_(M_a) q^(L_0-c/24)`. The momentum-plus-oscillator decomposition gives

`χ_0(τ)= [Σ_(n∈Z) q^(n²)] / η(τ)`,
`χ_1(τ)= [Σ_(n∈Z) q^((n+1/2)²)] / η(τ)`.

The leading exponents are `h_a-c/24`, hence

`χ_0=q^(-1/24)(1+3q+4q²+7q³+13q⁴+19q⁵+29q⁶+...)`,
`χ_1=q^(5/24)(2+2q+6q²+8q³+14q⁴+20q⁵+34q⁶+...)`.

Thus the earlier `q^(1/24)` is concretely the eta normalization of one bosonic oscillator. The character has its inverse `q^(-1/24)` times momentum and oscillator counts. This is not an arbitrary numerical fit and not a new identity.

The unnormalized trace `Tr q^(L_0)=q^(1/24)χ` is different. Multiplying by this common factor destroys the constant-matrix weight-zero modular transformation in that normalization. Keep the distinction explicit.

## 4. Modular transformations, derived directly

Write `Θ_a(τ)=Σ_n exp(2πiτ(n+a)²)` for `a=0,1/2`.

Poisson summation of the Gaussian gives

`Θ_a(-1/τ)=sqrt(-iτ/2) Σ_(m∈Z) exp(πiτm²/2) exp(2πima)`.

Splitting `m` into even and odd values yields

`Θ_0(-1/τ)=sqrt(-iτ)/sqrt(2) (Θ_0(τ)+Θ_(1/2)(τ))`,
`Θ_(1/2)(-1/τ)=sqrt(-iτ)/sqrt(2) (Θ_0(τ)-Θ_(1/2)(τ))`.

Together with the standard eta identity `η(-1/τ)=sqrt(-iτ)η(τ)` (the S-specialization of [NIST DLMF §23.18.5–6](https://dlmf.nist.gov/23.18)), this proves

`χ(-1/τ)=Sχ(τ)`, `S=(1/sqrt(2)) [[1,1],[1,-1]]`.

Use the consistent analytic square root on the upper half-plane. Its factors cancel in the quotient.

Under `τ -> τ+1`, the theta multipliers are `1,i`, since `n²` is integral and `(n+1/2)²=n(n+1)+1/4`. The eta product gives `η(τ+1)=exp(πi/12)η(τ)`. Therefore

`χ(τ+1)=T_char χ(τ)`,
`T_char=exp(-πi/12) diag(1,i)=diag(exp(-πi/12),exp(5πi/12))`.

**Do not conflate this with the ribbon matrix:**

`T_ribbon=diag(1,i)`.

The difference is exactly the common central-charge anomaly factor `exp(-2πic/24)`.

With positive dimensions `d_0=d_1=1`, normalized categorical S is the same displayed matrix: the double-braiding pairing is `b(a,b)=q_D(a+b)/(q_D(a)q_D(b))`, giving the unnormalized matrix `[[1,1],[1,-1]]`. It is nondegenerate. Both labels are self-dual, and direct matrix multiplication gives

`S²=I`, `(S T_char)³=I`,
`(S T_ribbon)³=exp(πi/4)I`.

The normalized quadratic Gauss sum `(1+i)/sqrt(2)=exp(2πi/8)` agrees with `c=1 mod 8`. This is the arithmetic/topological information retained by the ribbon category, while the character normalization retains the actual `c/24` factor. Verlinde's formula with this S gives exactly `N_ab^c=1` when `c=a+b mod 2`, zero otherwise; the checker verifies all eight entries.

## 5. Certificate and boundary

Run `python3 session34-vertex/check_lattice_voa.py`.

`LATTICE_VOA_CERTIFICATE.json` contains:
- Exact coefficients through oscillator level 30, using integer partition dynamic programming independently cross-checked against the Euler pentagonal recurrence.
- Exact algebraic verification of the three matrix relations and all Verlinde coefficients.
- Three 80-digit numerical S/T transformation checks, using 121 momentum terms and 400 eta factors. Residuals are below `10^-65`. These are controls, not interval-certified analytic proofs; Poisson summation and the eta transformation supply the identities.

This is a genuine VOA realizing the finite coherence example. It does not realize the earlier tournament winner algebra, recover SCC data lost by its associative quotient, solve a Yang–Mills equation, or establish a general tournament-to-VOA construction. It closes only the specific positive arithmetic–VOA–modular-category bridge requested here.
