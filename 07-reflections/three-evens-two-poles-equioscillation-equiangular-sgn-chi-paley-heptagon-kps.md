# Three evens, two poles: equioscillation ↔ equiangular is the sgn ↔ χ axis, pinned at the Paley heptagon

*kind-pasteur-2026-07-01-S11. The owner flagged a conflation — "even degree" and "even automorphism parity" are not the same — and named a triad to chase: equiangular lines, equioscillation, the Paley heptagon. Pulling the threads together (grounded by computation and by the repo's own `determinant-lens`, `twentyeight`, `forbidden-seven`, `equioscillation-Chebyshev`, and `even-graph-equinumerosity` reflections): there are **three** distinct "evens," two of them orthogonal; and **equioscillation and equiangular are the two poles of one axis** — the sgn↔χ axis of odd-function space — which pinches at n=7, the Paley heptagon, where a fourth "28" (max equiangular lines in ℝ⁷) joins the octonion apex.*

## Three evens, carefully separated

The word "even" is doing three different jobs. Keep them apart:

- **(F) even FUNCTION.** A `±1` function on vertex pairs. **Graph = even (symmetric); tournament = odd (skew)**, `S = A − Aᵀ`. This is the `det(I+S) = Σ_{K even} Pf(S[K])²` lens (`determinant-lens`). The whole project lives on the *odd*-function side; the even-function side is graphs.
- **(D) even DEGREE.** All vertices even degree = an **even graph** = the cycle space `H₁(K_n; F₂)`. This is the `E_n` metagraph, the Ising/flow/two-graph side.
- **(P) even PARITY.** An automorphism that is an **even permutation** (`sgn = +1`, in `A_n`).

**(D) and (P) are orthogonal, and (F) links them.** Computed (`even_degree_vs_even_parity_kps.py`):
- For **tournaments**, `|Aut|` is always odd (a transposition reverses an arc, so it is an *anti*-automorphism), hence every automorphism has odd order = a product of odd cycles = an **even permutation**: **`Aut(T) ⊆ A_n` always**. Verified: transitive, carousel `R₇`, Paley `P₇` all have `Aut` signs `{+1}` and anti-automorphisms (complementation `ι`) signs `{−1}`. So *being an odd function forces even-parity symmetry*, and the odd coset is exactly complementation.
- For **even graphs**, this fails: `C₄` has all-even degrees but the odd-parity automorphism `(1 3)`; at `n=4` *no* even-degree graph has `Aut ⊆ A₄`. Even degree and even parity cross-classify graphs into all four cells — **orthogonal**.

So the "even" of `E_n` (degree) and the "even" of automorphism parity are independent coordinates, and the tournament world is special precisely because its odd-function nature *couples* them: `sgn` becomes the orientation-preserving/reversing `Z₂` = the complement/converse mirror = the **dihedral reflection of the heptagon's `D₇`** = the `σ`-involution that splits Rédei-odd (witness) from LRC-even (existence) (THM-581/582). **Odd versions:** odd-degree graphs (all degrees odd) exist only for even `n` (a T-join coset of the cycle space); odd-parity automorphisms are `Aut ∩ (S_n∖A_n) ≠ ∅` — for tournaments this set is empty and its "shadow" is the anti-automorphism coset.

## Two poles: equioscillation and equiangular are the ends of the sgn↔χ axis

The `determinant-lens` shows odd-function space runs between two canonical odd functions:
- **`sgn`** (the carousel, `J = {1,…,(n−1)/2}`, `= R₇`): the *smoothest* odd function, locally transitive, determinant **floor** `d=1`.
- **`χ`** (the quadratic character, `J = QR`, `q ≡ 3 mod 4`, `= Paley`): the *flattest* odd function (perfect Legendre autocorrelation), doubly regular, determinant **ceiling** `d=max`.

**These two poles are exactly equioscillation and equiangular.**

- **`sgn` pole = EQUIOSCILLATION.** The smooth/low-frequency extremal is the Chebyshev/Fejér–Bochner object: the LRC covering floor *is* a Chebyshev equioscillation extremal + Cohn–Elkies magic function (HYP-3132, `the-cap-is-a-chebyshev-equioscillation-extremal`). Concretely, at the AP extremal the lonely function `g(t)=min_v‖vt‖` **equioscillates** — it hits its max `1/14` at the `φ(14)=6` atoms with *equal height*. Equal heights at many points = equioscillation. The atoms are the heptagon's units; the count `6` is the Verblunsky termination and the `MFAS` depth (last session). Equioscillation is the *smooth/circle/carousel* pole.
- **`χ` pole = EQUIANGULAR.** Flat autocorrelation = conference matrix = **Seidel two-graph = equiangular lines**. `A002854` counts even graphs = two-graphs = switching classes; the *regular* two-graphs are conference matrices = Paley. Equal angles between lines = equiangularity. Equiangular is the *flat/Hadamard/χ* pole.

The two poles **pinch at n=7**: `QR₇` is the unique H-maximizer *and* sits where "the two determinant ceilings touch" (`argmax H ⊂ argmax d`, the 6-class switching family). The **Paley heptagon is the equioscillation↔equiangular pinch** — the one place the smooth extremal and the flat extremal are the same object.

## 28 = N(7): a new face of the octonion apex

The `twentyeight` reflection collects the apex identities `28 = C(8,2) = T(7) = ` 2nd perfect number `= dim so(8)`, tower `7,14,21,28 = Im(𝕆),G₂,so(7),so(8)`, with `7 = ` Mersenne ∩ Heegner ∩ `(3 mod 4)`. **Add one more, from the equiangular pole:** the maximum number of equiangular lines in `ℝ⁷` is **`N(7) = 28`**, and the *absolute bound* `N(d) ≤ C(d+1,2)` is **tight exactly at `d ∈ {2,3,7,23}`** (the "tight" dimensions of regular two-graphs / tight spherical designs). So:

> `χ(E₇) = 28 = N(7) = C(8,2) = T(7) = dim so(8)` — the even-graph metagraph's chromatic number, the equiangular-line maximum in `ℝ⁷`, and the octonion arc-count are the *same* 28.

And `d = 23` is the *other* tight equiangular dimension — the project's other distinguished prime (Paley `T₂₃`, band prime 23). The equiangular "tight dimensions" `{2,3,7,23}` overlap the project's special primes `{3,7,23}`. This is the sharpest form of the unifying bet: **the absolute equiangular bound is tight precisely where the project's apex lives.**

## The unifying bet, sharpened

One self-dual object, two faces, pinned at Paley:

- **EVEN face** (graphs / even functions): the `K_n` Tutte/Ising system mod `S_n`. Even *degree* (`E_n`, flow, Curie–Weiss Ising — verified identity) ⊥ even *parity* (`Aut ⊆ A_n`). Its extremal node = regular two-graph = **equiangular lines**.
- **ODD face** (tournaments / odd functions): the sgn↔χ determinant axis, running from **equioscillation** (`sgn`/carousel/atoms/floor) to **equiangular** (`χ`/Paley/conference/ceiling). Odd-function nature *forces* even parity; complementation `ι` is the odd coset = the heptagon's `D₇` reflection.
- **The pin:** Paley `P₇` is the H-maximizer (odd face), the regular two-graph (even face), the equioscillation↔equiangular pinch, and the `|Aut|=21` obstruction to iso-covering (last session). `28 = N(7)` is where the equiangular bound goes tight, `n=8` (`= 7+1`) is the octonion arc-count, and `7` supplies all three proof pillars (Mersenne/Heegner/`3 mod 4`).

**n=7/8 is a single master threshold** seen from many sides: `E_7` loses perfection (odd holes); iso-covering breaks (`k_min(6)=7`, `k_min(7)=12`, Paley obstruction); the equiangular absolute bound goes tight (`N(7)=28`); the octonion apex sits at `28 = C(8,2)`; `sgn` and `χ` (equioscillation and equiangular) coincide. The even/odd of `n` itself (7 odd prime vs 8 = 2³) is the last hinge — the parity descent `14 = 2·7` peels the one free `Z₂` onto the all-odd apex-7 face.

## Honest status

- **Computed/verified:** the three-evens separation and `(D)⊥(P)` orthogonality; `Aut(T) ⊆ A_n` with `ι` odd (transitive, `R₇`, `P₇`); the Ising identity (prior); `28 = C(8,2) = T(7) = dim so(8)` (repo); `N(7)=28`, absolute bound tight at `{2,3,7,23}` (classical Lemmens–Seidel).
- **Grounded but not proven here:** `χ(E₇)=28 = N(7)` is a numerical coincidence with a *plausible* mechanism (both count a regular-two-graph / `C(8,2)` design object) — the mechanism is the concrete next computation (Seidel spectrum of `E₇`'s extremal node; is it the conference two-graph?).
- **Conjectural reframe:** "equioscillation = sgn pole, equiangular = χ pole" is a lens, not a theorem — the precise statement to prove is that the LRC equioscillating floor certificate (HYP-3132) and the equiangular/Seidel ceiling are the two ends of one extremal problem on odd-function space. If true, OPEN-Q-108 becomes an **equiangular-lines / tight-frame extremal bound at `d=7`**.

— Related: `the-determinant-lens-sgn-vs-chi-and-the-three-geometries.md` (sgn↔χ, the two poles), `twentyeight-the-octonion-apex-and-the-three-pillars.md` (28, the apex), `forbidden-seven-in-all-senses.md` (|Aut| odd ⟹ Aut⊆A_n), `the-cap-is-a-chebyshev-equioscillation-extremal-...` (HYP-3132, equioscillation), `even-graph-equinumerosity-one-cube-four-faces-eight-wild-lines-kps.md` (two-graphs = even graphs, the Ising/Tutte face), `two-axes-of-the-tournament-metagraph-...` (Paley as shared extremum), THM-580/581/582, HYP-3802 (heptagon), OPEN-Q-108. Scripts: `even_degree_vs_even_parity_kps.py`, `even_graph_equinumerosity_probes_kps.py` (+ .out). Not a HYP reservation (synthesis); the `χ(E₇)=N(7)` mechanism and the equioscillation=sgn / equiangular=χ statement go to INVESTIGATION-BACKLOG.
