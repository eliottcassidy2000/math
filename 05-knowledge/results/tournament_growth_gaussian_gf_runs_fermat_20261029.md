# Tournament growth: a generating function for A000568, partial structures of the metagraph, T5's four spurious Z/6, and Fermat primes

*mac-mini-2026-10-29-growth, 2026-10-29. Pre-audit. Theorems: THM-4845, THM-4846, THM-4847. Open problems: HYP-9445.*
*Experiments: [`04-computation/experiments/growth_gf_runs_fermat_20261029/`](../../04-computation/experiments/growth_gf_runs_fermat_20261029/).*

**The owner's directive (verbatim).** "also merge in tournament related past work and aiming to find a recursive
formula for a000568 that gives a generating function. think about coming up with more classifications for partial
structures encountered in the metagraph and how they change as the tournament size increases. see how the tournament
on 4 vertices with 6 edges relates to the following: dropping T5's integral constraints introduces four spurious
\(\mathbb Z/6\) factors. think about ways of proving that there are exactly 5 fermat primes. continue looking for
recurrences and exact formulas for structure yo can spot"

**Concurrent work, same day, credited throughout.** This session claimed its lanes by message (MSG-3286). Two other
agents worked the same directive and pushed first; their packets were read only after this session's checkpoint
a37e12900d.
- **session98** (`session98_tournament_growth_20261029`, 08:43) has:
  - the fixed-degree projective recurrence;
  - the meromorphic continuation of `Σ n! a_n 2^(−C(n,2)) z^n` with poles exactly at `4^k` of order `3k + 1`;
  - non-P-recursiveness;
  - a palette generalization;
  - partial-construction classifications.

  Priority for the pole theorem is theirs. Our §1 pole statements are an independent derivation that agrees exactly.
- **The research packet** `four_vertex_glue_fermat_20261029` (09:00) has:
  - the vacancy-walk / star-graph form of the one-arc observer;
  - the sharp explanation of the four `Z/6`, through the kernel lattice `6K` and the exact sequence with quotient
    `(Z/6)^4` and an order-five clock;
  - Fermat squaring dynamics.

  Priority for those is theirs. Our §3 account is the complementary combinatorial one.
- **New here:**
  - the Gaussian generating function and its unification of the palette family;
  - the permanental compound;
  - the T5 star proof;
  - the general ghost lemma and the all-class censuses by order;
  - the Hamiltonian-cycle/torsion conjecture;
  - the function-field analogue and the tournament counting forms of Pépin.

## 0. Prior work merged

| topic | where it is |
|---|---|
| A000568 Burnside/Davis count, 200 terms | THM-305; `tournament_future_boundary_20261020` (finite-horizon gcd responses, independent table through 76); S48 Burnside enumerators |
| core kernel and quadratic form | THM-514 (unit tail; integer kernel `K[m,t]`); THM-516 (`e = C(t,2) + ½ Σ φ(d) M_d²`, Smith); OPEN-Q-105 / HYP-2538 (closed θ-product?) |
| PDE, hierarchy, quotients | `weighted_core_transport_20261020`; session98 `PROJECTIVE_RESPONSE.md` (concurrent) |
| strong components, modular tree | Wright's recurrence `T = Σ C(m,k) s(k) T(m−k)` (THM-4746); THM-4803 (A000568 from the prime cycle-index census) |
| the metagraph | definitions.md tiling model; THM-643/644; THM-1430; THM-4729/4730/4731 (closed/open tilings, pure-blue alphabet, self-lines); THM-4739 (Kirchhoff laws); THM-588/589 |
| run observers, ports, T5's Z/6 | `overnight_rank_return_20261029/THREE_RUN.md` (T5 kernel = tetrahedron boundaries), `TWO_RUN.md`, `COMPLEMENT_KERNEL.md`; `overnight_research_20261029/THREE_RUN_GLUE.md` (index-54 chart; T5 Smith forms (R7); four spurious Z/6) |
| the four 4-vertex tournaments | THM-4823; `trees5_tournaments4_converse_20261005.md`; reflection "the triangular number is the n = 4 metagraph" |
| Fermat primes | THM-4644 and THM-871 (Fermat primes = rigid circulants); THM-4793/4800 (Pépin as Collatz dynamics); THM-4801 (constructible polygons); S58 §1.2 (covering, algebraic and Aurifeuillean no-go); HYP-9410 (calibrated heuristic); THM-4746/4747 (Camion tests) |
| concurrent (2026-10-29) | session98 `limits/ALL_POLES_CORE_EXPANSION.md`, `PALETTE_POLE_CLASSIFICATION.md`, `GROWTH_THEORY.md`; `four_vertex_glue_fermat_20261029/README.md` and `FERMAT.md` |

## 1. A recursion for A000568 that gives the generating function (THM-4845)

**What cannot exist.**
- *No linear recurrence in `n` with polynomial coefficients.* The sequence `a_n` grows like `2^(C(n,2))/n!`, faster
  than any P-recursive sequence.
- *No such recurrence for the normalization either.* The normalized series `F(y) = Σ (n! a_n / 2^(C(n,2))) y^n` has
  infinitely many poles, so it is not D-finite.

**What does exist.** One Gaussian variable for each odd divisor `d`, with variance `φ(d) ln 2`, diagonalises the
Burnside exponent. Equivalently, `X` is a Gaussian process on the odd integers with covariance `ln 2 · gcd(p,q)`.
Conditionally on the field, the classical Newton recursion `n Z_n = Σ_{p odd} W_p Z_{n−p}` with `W_p = 2^(−1/2) e^(X_p)`
has generating function `exp(Σ W_p x^p/p)`. The tournament numbers are its average:

    Σ a_n x^n = E[ exp( Σ_{p odd} W_p x^p / p ) ],        a_n = E[Z_n(W)].

- **The core kernel.** For THM-516's kernel the same formula holds with parts ≥ 3. This answers OPEN-Q-105(a): the
  product over part sizes sits inside the expectation, and the cross-coupling `½ Σ φ(d) M_d²` is carried by the
  shared field.
- **The family.** The same field on all integers gives graphs, and tournaments are graphs with the even power sums
  switched off. Doubling the covariance gives digraphs, and base 3 gives oriented graphs. In general, session98's
  palette family with `q` colours, `t` of them reversal-fixed, is the field with covariance `ln q · gcd`, odd weights
  `q^(−1/2) e^X` and even weights `(t/q) e^X`. The self-converse tournaments (A002785) are the twisted version
  `E[(1 + W_1 x) exp(Σ_{j ≡ 2 (mod 4)} e^(X_j) x^j / j)]`, checked through order 16.
- **The inherited PDE.** It is the Cameron–Martin formula for this field. Inserting a `p`-cycle shifts the field by
  `ln 2 · gcd(p, ·)`, and session98's projective state is that shift taken modulo the dilation direction.
- **The normalized series is meromorphic** (session98 first; found independently here). `F(y)` extends to all of
  `C`, with poles exactly at `4^k`, of order `3k + 1`. Hence the exact expansion

      a_n = (2^(C(n,2)) / n!) · (1 + 32 C(n,3)/4^n + (98304 C(n,5) + 655360 C(n,6))/16^n + …),

  with explicit polynomials `P_k` of degree `3k`. For graphs the poles sit at every `2^d`, of order `2d + 1`, and the
  first correction is Pólya's.

## 2. Partial structures in the metagraph, and how they change with the order (THM-4846, HYP-9445)

**The structures.** A metagraph tiling is a tournament with its base Hamiltonian path. Cutting that path into `k`
runs gives a spanning cover by `k` paths. These covers, counted by start and end ports, are the overnight packets'
observers `M_k`. They are the natural partial structures of the tiling.

**Laws for all orders.**
- **One arc** (`k = n − 1`): `M_{n−1}` is an orientation of the star-transposition graph of `S_n`, the Nauru graph at
  `n = 4`. Its run-symmetric block is `Aᵀ`, its sign block is a signed `Aᵀ`, and its determinant is a product of
  powers whose exponents are combinations of `S_{n−1}`-dimensions.
- **Two arcs** (`k = n − 2`): the unordered observer is the second permanental compound of `A`.
- **Sectors as compounds.** For every `k`, the plus and sign sectors are the multilinear parts of the permanental and
  determinantal compounds of the path matrix `P = X(I − AX)⁻¹`. By Jacobi, the sign sector is a ratio of
  cycle-cover minors of `I − AX`, which is where Hamiltonian cycles enter.

**Thresholds as the order grows.**

| order | what happens |
|---|---|
| 4 | The strong tournament has unimodular observers; each non-strong one loses one free run orbit; no torsion and no glue loss. |
| 5 | The first torsion appears. Four of the twelve classes lose information to the `S_3` glue, T5 among them. T5 is the only boundary killer. Every two-arc kernel is spanned by vertex stars. |
| 6 | The star law fails (12 of 56 classes). No boundary killer exists. The first `S_2` glue losses appear (18 classes). |
| 7 | Unique-Hamiltonian-cycle classes first carry torsion (5 of 21). |
| 8 | Star law: 2808 of 6880 classes fail. Every strong class with at least two Hamiltonian cycles still has torsion in `M_2` (5953 of 5953). |

**The conjectures (HYP-9445 (1), (2)).**
- *Two Hamiltonian cycles force torsion.* Two or more Hamiltonian cycles force torsion in the two- and three-run
  observers. This holds in all 6315 strong cases checked for `k = 2` (orders up to 8), and in 1044 random or
  low-cycle samples at orders 9 and 10. It is false for the one-arc observer, and the torsion is not always
  2-primary.
- *Unique cycles.* The unique-cycle classes are Douglas's `F(2n−6)` = 1, 3, 8, 21, 55 (orders 4 to 8). Their
  torsion-free members number 1, 3, 8, 16, 36.
- *The sign sector alone.* The signed two-run sector already shows the effect: it has torsion for every class with
  two or more Hamiltonian cycles through order 8. For unique-cycle classes it is torsion-free except one class at
  order 8.

## 3. The tournament on 4 vertices and T5's four spurious Z/6 (THM-4846 (C), (D))

**The kernel.** T5's three-run observer has a 4-dimensional rational kernel, the tetrahedron boundaries of THREE_RUN
T4. In the permanental-compound coordinates these are the five vertex stars, with one relation.

**Each star is a 4-vertex tournament.** The star at `w` lives on `T5 − w`, which is the strong 4-tournament: 4
vertices and 6 arcs, made of a Hamiltonian 4-cycle and two diagonals. As a vector it is `+1` on the ends of one
diagonal and `−1` on the ends of the other. The identity that kills it is an alternating-sign cancellation along
in-neighbourhoods. It works because T5's in-neighbourhoods are cyclic intervals of even length.

**The ghost lemma explains the four Z/6.**
- *The 2-part.* The number of spurious 2-factors equals the rank of the plus-minus connecting map
  `B = (S − P)/2` from the mod-2 kernel to the mod-2 cokernel of the plus sector. For T5 that kernel is the cut space
  of `K5`, which has dimension 4. `B` maps it isomorphically, giving four spurious `Z/2`.
- *The 3-part.* The mod-3 glue of the index-54 chart gives four spurious `Z/3`. THM-4846 (D′) gives the general mod-3
  formula: two block-triangular pieces in the basis adapted to the invariant flag, and one connecting rank when the
  sign sector is invertible mod 3, as for T5.

Together these are the four spurious `Z/6`, one for each independent star.

**The tetrahedron on its own hides nothing.** The strong 4-tournament orients the Nauru graph with determinant 1.
The hidden directions, and their ghosts, need five overlapping copies of it inside T5.

**Typing.** The 6 of `Z/6` is `|S_3|`, where `Z[S_3]` fails to split at 2 and 3. That the 4-tournament also has 6
arcs is the coincidence `3! = C(4,2)`, typed NUMEROLOGY. The structural bridge is the 24 ordered triples inside a
4-set, which are the 24 elements of `S_4` and the 24 vertices of the Nauru graph.

## 4. Ways of proving that there are exactly five Fermat primes (THM-4847)

**Status.** Not proved, and no current method approaches it. It is a `Π⁰₁` statement, and it implies the open
problem that infinitely many Fermat numbers are composite. The Fermat numbers `F_5`–`F_32` are known composite;
`F_33` is the first of unknown character.

**Proved barriers.** Four natural routes are ruled out. The first and third are S58 §1.2's no-go list, restated;
the concurrent FERMAT.md has the dynamical twin of the fourth (the source 2 is blind).
- covering congruences, because a prime divides at most one Fermat number;
- fixed-modulus certificates, because `F_n mod M` is eventually periodic and coprime to `M` except at `ω(M)` indices;
- polynomial or Aurifeuillean algebra, because `X^(2^n) + 1` is irreducible over `Q`;
- tournament counting modulo `F_n`, because `2^(C(N,2)) ≡ 2^((N−1)/2) ≡ 1` for every Fermat number, prime or not.

**What does see Fermat primality.**
- *Symmetry.* Fermat primes are exactly the primes all of whose circulant tournaments are rigid (THM-4644).
- *Three-state counting.* `3^(C(N,2)) ≡ −1 (mod F_n)` holds iff `F_n` is prime. This is Pépin's test, counting
  labelled oriented graphs.
- *The function-field analogue.* In `F_q[T]` the statement "only finitely many Fermat primes" is a theorem, with no
  Fermat primes beyond `n = 1` at any base. The reason is that `(Z/2^(n+1))^*` is not cyclic; over `Q` this
  mechanism is absent.

**Heuristics.** Boklan–Conway bound the expected number of new Fermat primes by about `2^(−30)`. HYP-9410's
calibrated model gives about `7·10^(−9)`.

## 5. Other exact formulas and recurrences found on the way

- **The polynomials `P_k`.** The leading coefficient of `P_k` is `((3k)!/(3^k k!)) 2^(3k² + k)` at `C(n, 3k)`.
  Truncating to `K` poles gives a constant-coefficient linear recurrence that holds up to `O(4^(−(K+1)n))`.
- **Determinants of the one-arc observer.** They factor as `Π det(M_λ)^(f^λ)`. At order 6:
  `7^10`, `3^30`, `2^42`, `23^10`, `2^36 11^10`, `2^40 13^10`, `2^24 3^28 5^10`, `2^8 3^12 5^10 7^10`, `2^8`, and `1`
  for eight classes.
- **The T5 numbers.**
  - Smith forms: one-arc `1^112 2^8`; two-run `1^18 82^2`; three-run `1^52 3^2 6^2 0^4`.
  - The other `H = 15` class of order 5: three-run `1^50 3^2 9^2 0^6`, ghosts `(0, 5)`.
- **Star-law failure counts.** 12, 168 and 2808 at orders 6, 7 and 8.
- **`F_5`.** It fails both strong rotation tests at both prime factors. The residues are `r_641 = 257 = F_3` and
  `r_6700417 = 641`: in general `r_p` is the cofactor reduced mod `p − 1`.

## 6. Reproduction

```bash
cd 04-computation/experiments/growth_gf_runs_fermat_20261029
python3 gf_gaussian.py 80          # THM-4845 checks -> out/gf_gaussian.out
python3 t5_check.py                # T5 Smith forms (reproduces THREE_RUN_GLUE (R7))
python3 census_runs.py 5 3 --snf   # any n, k; outputs in out/census_n*_k*.out
python3 star_orientation.py; python3 block_checks.py; python3 ghost_lemma.py; python3 ghost3_lemma.py
python3 jacobi_sign_sector.py      # THM-4846 (F)
python3 compound_kernel.py 6       # n = 4..8 -> out/compound_n*.out
python3 fast_runs.py 8             # two-run census with Hamiltonian cycles -> out/fast_m2_n8.out
python3 torsion_vs_hc2.py          # HYP-9445 (1), (2) tables
python3 fermat_camion.py           # THM-4847 (C) -> out/fermat_camion.out
python3 random_torsion.py 9 300 300 2   # HYP-9445 (1) sampling at orders 9, 10 -> out/random_torsion.out
```

The scripts need `gentourng` (nauty), PARI/gp, numpy, sympy and networkx.
