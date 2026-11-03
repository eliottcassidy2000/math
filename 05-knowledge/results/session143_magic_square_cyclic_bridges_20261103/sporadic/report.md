# Session 143: scoped C5–Baby Monster bridge and arithmetic obstructions

Status: **PROVED** for the general complex-representation lemma and field obstruction; **CITED / VERIFIED TRANSCRIPTION** for Baby Monster character inputs; **FINITE-EXACT** for the accompanying checks. No claim of originality or universal correspondence. Source checkout: `8286d710498a3b4e9b9f1ac683c6eb77716f3f4b`. No repository changes made by this lane.

## 1. What survives

Let `g` have order five in any group, let `rho` be a finite-dimensional complex representation of degree `d` and character `chi`, and put `zeta=exp(2*pi*i/5)`. Define the **group-algebra operator**

`B_g = rho(g) + rho(g)^2`.

The matrix `rho(g)` is diagonalizable because its minimal polynomial divides the separable polynomial `x^5-1`. If its `zeta^r`-eigenspace has dimension `m_r`, then

`m_r = (1/5) sum_{j=0}^4 chi(g^j) zeta^(-r*j)`

by Fourier inversion on `C5`. Consequently the eigenvalues of `B_g` are

`alpha_r = zeta^r + zeta^(2*r)`, with multiplicity `m_r`.

Thus there are **at most** five distinct eigenvalues; some multiplicities may vanish. This construction needs neither sporadic simplicity nor a tournament. For the regular permutation representation of `C5`, every multiplicity is one, and `P+P^2` is the adjacency matrix of the regular five-vertex tournament with arcs `i -> i+1,i+2` modulo five (a transpose convention changes no spectrum).

Write `q(x)=x^4+2*x^3+4*x^2+3*x+1`. Direct substitution modulo `Phi_5(z)` gives `q(z+z^2)=0`, hence all four nontrivial `alpha_r` are roots. Modulo two, `q=x^4+x+1` has no linear factor and is not divisible by the only irreducible quadratic `x^2+x+1`. It is irreducible over `Q`, so the four roots are distinct and each generates `Q(zeta_5)`. The remaining eigenvalue is `alpha_0=2`.

If `chi(g^j)=c` for every `j=1,2,3,4`, Fourier inversion simplifies to

`m_0=(d+4*c)/5`, and `m_1=...=m_4=(d-c)/5`.

Therefore

`det(xI-B_g)=(x-2)^((d+4*c)/5) q(x)^((d-c)/5)`.

The assumptions imply that these exponents are nonnegative integers. A rational-valued character has this constant-power property, as does any character on a rational conjugacy class. In the regular tournament case this gives `(x-2)q(x)=x^5-5*x^2-5*x-2`.

**Boundary:** this is not an identification of `B_g` with `rho(g)`, a character value, a Baby Monster orbital graph, or a tournament adjacency matrix in dimension `d`. A one-dimensional nonrational `C5` character gives only one `alpha_r` and refutes an unconditional equal-multiplicity claim. In characteristic five separability fails, and the five roots coalesce.

## 2. Exact Baby Monster specialization

The [ATLAS Baby Monster page](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B) gives

`|B|=4154781481226426191177580544000000 = 2^41 3^13 5^6 7^2 11 13 17 19 23 31 47`.

Its two order-five classes are `5A` and `5B`, with centralizer orders `443520000` and `6000000`. It lists a maximal subgroup of type `5:4 x HS:2`; this is a genuine occurrence of `5:4`, not evidence that `5:4` is a tournament automorphism group.

[Breuer–Magaard–Wilson's explicit CTblLib verification protocol](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/htm/ctblbm.htm), sections 6–7, supplies the rational ordinary character of degree `4371`. Its character vector is explicitly reordered into ATLAS class order. The first 17 classes comprise `1A`, four order-two classes, two order-three classes, and ten order-four classes. Positions 18 and 19 therefore give `chi(5A)=21`, `chi(5B)=-4`. The section-6 power-map construction makes both classes individually rational; section 7 ends by verifying agreement of all reconstructed irreducibles with the library table. The scratch checker independently matches all 184 class labels to the vector and checks these two entries. This verifies transcription, not a new construction of the character table.

The proved lemma now yields:

| class | `c` | `m_0` | each `m_1,...,m_4` | characteristic polynomial of `B_g` |
|---|---:|---:|---:|---|
| `5A` | 21 | 891 | 870 | `(x-2)^891 q(x)^870` |
| `5B` | -4 | 871 | 875 | `(x-2)^871 q(x)^875` |

Because both classes are rational, every ordinary character value at either class is a rational algebraic integer, hence an integer; the nonreal `alpha_r` are not those character values. The traces are respectively `42` and `-8`, not zero. Thus these operators are not loopless adjacency matrices in any basis. Their common nontrivial eigenvalue support is generic cyclic-order-five data; their multiplicities distinguish these two classes **after** fixing this representation and importing its character values. The five-vertex spectrum alone does not choose between them.

[Allan Steel's thesis](https://magma.maths.usyd.edu.au/~allan/reps/AllanSteelPhD.pdf), section 5.11, printed page 119 / PDF page 129, constructs the minimal faithful ordinary representation in dimension `4371` over `Q`. Since `B` is simple, every nontrivial representation is faithful. Consequently every characteristic-zero 24-dimensional representation of `B` is trivial: a nontrivial 24-dimensional Baby Monster linear action is **REFUTED**, not merely unsupported. This does not purport to classify modular or nonlinear actions.

## 3. Why Q(zeta_17) is not the eigenvalue field

The compositum `Q(zeta_5,zeta_17)` is `Q(zeta_85)`: Bézout's identity recovers a primitive 85th root from the two roots. Its degree is `phi(85)=64`, while the two constituent degrees are `4` and `16`. The degree formula for these Galois extensions gives

`Q(zeta_5) intersection Q(zeta_17) = Q`.

Since every nontrivial `alpha_r` generates `Q(zeta_5)`, none belongs to `Q(zeta_17)`. There is no `Q`-embedding of its eigenvalue field into `Q(zeta_17)`. A literal identification of the regular five-tournament eigenvalues with elements of that cyclotomic field is therefore **REFUTED**. A more elaborate relation would need different objects and an explicit map.

The [Daniels–Lozano-Robledo–Najman–Sutherland paper](https://math.mit.edu/~drew/mcom3213.pdf), section 6.1, proof of Proposition 6.2 (printed page 446 / PDF page 22), confirms that `X_0(13)` has genus zero and parametrizes rational 13-isogenies. Its rational parametrization is

`j(t)=(t^2+5t+13)(t^4+7t^3+20t^2+19t+1)^3/t`.

Hence `Q(X_0(13))=Q(t)` has transcendence degree one, unlike the degree-sixteen number field `Q(zeta_17)`. They cannot be isomorphic as `Q`-fields. Also `S_2(Gamma_0(13))=0` because its dimension is the genus. This rules out a nonzero weight-two cuspidal Hecke system at this level, not Hecke correspondences or all other weights. One may choose `t=zeta_17` and specialize functions regular there, but that is extra, noncanonical point data and not a homomorphism from the whole function field. No Hecke-equivariant map to the tournament or Baby Monster data is supplied or established.

## 4. Conway's M13 is not a five-point Mathieu group

The primary [Conway–Elkies–Martin paper](https://arxiv.org/html/math/0508630v1), introduction and section 2.2, uses the 13-point projective plane `PG(2,3)`, twelve counters, and one hole. Each line has four points. A move with hole at `p` and chosen `q` is the double transposition `(p q)(r s)` on their line `{p,q,r,s}`. Legal concatenation retains the hole endpoint. Closed move sequences give `M12`; the authors' `M13` denotes the larger pseudogroup of legal paths and is explicitly not a group. The paper discusses five-/six-fold transitivity, which counts how many ordered points may be prescribed, not a board with five points. A five-point tournament needs a new encoding with endpoint and incidence sidecars; no such encoding is established here.

## 5. Maps and losses

| proposed map | preserved | lost / missing | cheapest decisive test |
|---|---|---|---|
| `(rho,g)` to `rho(g)+rho(g)^2` | restriction multiplicities for this injective polynomial on fifth roots | global group multiplication, conjugation data, vertex incidence | compare traces and dimensions |
| restriction to unweighted eigenvalue support | possible `C5` frequencies | all multiplicities, therefore class distinction | 5A and 5B have identical support |
| eigenvalue field to `Q(zeta_17)` | no literal `Q`-field inclusion possible | conductor/ramification information cannot be ignored | degree of compositum is 64 |
| `M13` puzzle to a five-point graph | no specified map | 13-point incidence, hole endpoint, legal composition | define source/target and each move image first |

Reproduction: `python3 check_bridge.py` runs without third-party packages and executes 45 explicit checks; with source snapshots it executes 51. All checks use explicit errors and remain active under `python3 -O`. `python3 test_bridge.py` runs portable positive/negative controls in both modes; adding `--sources-dir PATH` also checks source transcription. The full 10-case regression passes, rejecting deliberately wrong eigenvalue, character, and transcription inputs with exit code 1 under both normal and optimized Python. With the separately downloaded primary snapshots, add `--sources-dir PATH` to check transcription. `verify_baby_monster.g` provides a GAP/CTblLib cross-check but was **NOT EXECUTED** locally because GAP is absent. The manifest pins portable deliverables and separately records source-snapshot hashes; it intentionally does not require committing the third-party snapshots. No LRC conclusion follows.
