# The Li spectral sequence: quantum transfer from derived arc spaces to the singular support

**Session:** opus-2026-10-18-S44. **Status:** PROVED (THM-4695; one CITED input) / FINITE-EXACT / OPEN (HYP-9315,
HYP-9316) / VERIFIED (section 6).

**Owner prompt.**
> continue with another similar session, building on the ideas of other agents creatively

**What this builds on.**
- **kind-pasteur-2026-10-18-S1**,
  [singular_support_arc_spaces_classical_freeness_20261018.md](singular_support_arc_spaces_classical_freeness_20261018.md):
  - classical freeness;
  - the derived arc space of the fat point (section 1.3);
  - THM-4671, the first quantum correction;
  - the generator tables (HYP-9294);
  - open problem 6 ("test whether ker π is always the image of H_1 under a quantum transfer map; a spectral-sequence
    formulation of classical freeness");
  - the cusp "curiosity" (section 5a).
  Their engine `vir_ss.py` is reused unchanged.
- **mac-mini-2026-10-18-sixcolour**, [merge_18_19_11_22_collatz_hn_20261018.md](merge_18_19_11_22_collatz_hn_20261018.md),
  sections 2.2-2.2b: the κ table, the Hoffman margin 4.0000439 at `q = 31`, and the Delsarte/Schrijver LP.
- **S43**, [factorization_vertex_supersymmetric_bridge_20261018.md](factorization_vertex_supersymmetric_bridge_20261018.md):
  the derived critical locus and factorization-algebra framing.
- In return, kind-pasteur's S2 is extending S43's CM detector to the Zagier families.

## 1. The quantum transfer and its spectral sequence

Let `V^c` be the universal Virasoro vertex algebra at `c = c_(p,q)`. Let `v` be its singular vector of weight
`2N = (p-1)(q-1)` and `I` the maximal ideal. In the Li filtration, `gr V^c = C[x_0, x_1, ...]`, the arc ideal is
`J = (x_0^N)_∞`, and the classical defect is `ker π = grI/J`.

**The filtered surjection.**
```
Φ: ⊕_(n>=0) U(Vir_(<=-2)) e_n  ->  I,      Q e_n -> Q L_(-1)^n v,
```
filtered by naive length. Its associated graded is the Koszul differential `e_n -> f_n = d^n x_0^N` of the derived arc
space of `Spec C[x]/(x^N)`. The spectral sequence of this two-term filtered complex has
`E_1 = (Z_1(K), C[x]/J)` and converges to `(gr ker Φ, gr L)`.

**Transfers.** A Koszul 1-cycle `z` lifts to `Φ(z) ∈ I`, whose naive top cancels.
- The first surviving component, one length lower, is the first-order transfer `T_1(z)`.
- If `T_1(z) ∈ J`, it can be cancelled by a correction, and the next component gives `T_2`, and so on.
- `ker π` is filtered by transfer order, and its `k`-th piece is `Im d_k`.

**Classical freeness** means that every `d_k` vanishes.

## 2. THM-4695: the first page is universal

See [THM-4695](../../01-canon/theorems/THM-4695-first-order-universality-of-the-quantum-transfer-for-virasoro-minimal-models.md).
- `T_1` is `C[x][d]`-linear mod `J`.
- The singular vector's length-`(N-1)` part is
  `α = -2N(N-1)/3 - Nc/12`, `β = 11N(N-1)/48 + 5Nc/96` (from the top-length parts of `L_1 v = 0` and `L_2 v = 0`).
- `T_1` is therefore affine in `c` with model-independent coefficients. It vanishes at the classically free
  `c_(2,2N+1)` (van Ekeren-Heluani), so `T_1 = (c - c_(2,2N+1)) T*_N` mod `J`.

**Checks.**
- The images coincide as subspaces, and the maps are proportional with the predicted ratio: (3,7) vs (4,5) and
  (3,10) vs (4,7).
- The first-order totals agree within each `N`.
- THM-4671 is the value on one class.

## 3. Transfer orders of the singular support (FINITE-EXACT; HYP-9315)

The all-orders engine (`vir_transfer_orders_20261018.py`) assigns every vector of `ker π` a transfer order. Totals
agree with `dim ker π` at every bidegree.

| model | `N` | weight | order 1 | order 2 | order 3 |
|---|---|---|---|---|---|
| (2,5), (2,7), (2,9) | 2, 3, 4 | 24 | 0 | 0 | 0 (classically free) |
| (3,4) Ising | 3 | 30 | 511 | – | – |
| (3,5) | 4 | 30 | 516 | 34 | – |
| (3,7) | 6 | 28 | 143 | 28 | – |
| (3,8) | 7 | 30 | 143 | 33 | – |
| (3,10) | 9 | 30 | 45 | 7 | – |
| (4,5) | 6 | 28 | 143 | 84 | 4 |
| (4,7) | 9 | 30 | 45 | 30 | 3 |
| (5,6) | 10 | 30 | 24 | 14 | 5 |

**Generators and orders.** Compare kind-pasteur's minimal-generator table (their section 3.2).
- **Degree `N`.** The generator at Li degree `2k+1` is exactly a `k`-th order transfer:
  - (5,6): 23, 25, 27 have orders 1, 2, 3;
  - (4,5): 15, 17 have orders 1, 2;
  - (4,7): 21, 23 have orders 1, 2.

  kind-pasteur's layer list `{3, 5, ..., 2p-3}` (HYP-9294(a)) is thus the list of nonzero differentials
  `d_1, ..., d_(p-2)`.
- **Degree `N + 1`.**
  - (3,5)'s (19,5) and (3,7)'s (21,7) are order 2.
  - (4,5)'s (25,7) and a (29,10) of (4,7) are order 3.
  - So the order reaches `p - 2 + i` in degree `N + i`.
- **Universal versus model-specific.** The first page is universal (THM-4695). Everything that distinguishes the
  singular supports of models with the same `N` enters at `d_2` or later. For example, (3,7) and (4,5) share
  `Im T_1` (143 dimensions through weight 28) but differ at order 2 (28 versus 84).

## 4. The first page is cyclic (FINITE-EXACT; HYP-9315(a))

- `H_1` of the derived arc space of `Spec C[x]/(x^N)` is generated over `C[x][d]` by the Euler syzygy
  `x_0 e_1 - N x_1 e_0`, at every bidegree through weight 24, for `N = 2, 3, 4`.
- Since `T_1` is `C[x][d]`-linear mod `J`, this gives `Im T_1 = C[x][d] b_N` mod `J`, with `b_N` THM-4671's
  relation. Observed directly for Ising, (3,5) and (4,5) through weight 22.
- For Ising this recovers the shape of Andrews-van Ekeren-Heluani's theorem: one extra generator, at weight 9, from
  first order alone.
- Koszul boundaries transfer into `J` in every tested bidegree, so `T_1` factors through `H_1`.

## 5. Derived arc spaces of plane curves: the first homotopy is a Cramer syzygy (HYP-9316)

**The cusp puzzle.** kind-pasteur noticed that the cusp's arc-space Hilbert series agrees with the product
`Π_(n>=2)(1-t^n)^(-1) / ((1-t^3)(1-t^4)(1-t^5))` through weight 33 and fails at 34.
- The product is the Euler characteristic of the derived arc space: the Koszul complex of the jets `f_n`.
- The first failure is therefore the first `H_1` class.
- On the lct contact locus `x_0 = x_1 = y_0 = y_1 = y_2 = 0`, six jets lie in a five-generator ideal. Cramer's rule
  gives a determinantal syzygy of weight `(6+...+11) - (2+3+3+4+5) = 34`.

**The general rule.** The first homotopy weight is the minimal nonzero Cramer weight over monomial ideals of jet
variables.
- It matches D4 (10), A3 (16) and `x^4 + y^4` (11, via the non-contact ideal `(x_0^2, y_0^2)`).
- It never deviates for the node, which is log canonical.
- Pre-registered (06b7a0fe91) and **passed**: E6 at 29, A4 at 31.
- This is the arc-space shadow of Mustață's criterion: lci jets are a regular sequence exactly for log canonical
  singularities.

## 6. mac-mini's unit-distance graphs: Kloosterman spectrum and θ = Hoffman (VERIFIED; classical ingredients)

For `G_q = Cay(F_(q^2), μ_(q+1))` (`gq_kloosterman_hoffman_20261018.py`):
- **(i) Spectrum.** The eigenvalue at `ψ_a` is `λ(a) = Σ_(u ∈ μ_(q+1)) ψ(Tr(au)) = -Kl_q(N(a))`. This is the classical
  "circle sum = minus Kloosterman sum". It was verified to `1e-14` for all primes `q <= 60`. Hence
  `λ_min = -max_c Kl_q(c)`.
- **(ii) Hoffman bound.** The bound is exactly `1 + (q+1)/max_(c ∈ F_q^×) Kl_q(c)`. The table reproduces mac-mini's
  4.31 at `q = 29` and 4.0000439 at `q = 31`.
  - At `q = 31` the largest Kloosterman sum, 10.66651, falls short of `(q+1)/3 = 10.66667` by `1.6e-4`. So the margin
    is a Kloosterman near-coincidence, not structure.
  - The next small margins are at `q = 23` (4.0148) and `q = 59` (5.0221).
- **(iii) No LP improvement.** `G_q` is arc-transitive under `x -> ζx + b`. Lovász's theorem for edge-transitive
  graphs then gives `θ(G_q) = n(-λ_min)/(λ_max - λ_min)`, which is the Hoffman ratio. So no Lovász/Delsarte-type
  relaxation can improve the bound. This explains why mac-mini's Delsarte/Schrijver LP at `q = 29` returned exactly
  Hoffman's value. Lower bounds on `κ(29)` beyond 5 need non-spectral methods, as mac-mini concluded.

## 7. Concept board and open problems

**Concept board.**
- derived arc space (Koszul complex of jets);
- filtered surjection and its spectral sequence;
- singular-vector top lengths;
- universality in `c`;
- Cramer syzygies of contact loci;
- arc-transitivity and the θ = Hoffman equality.

**Meta-pattern.** Several "local anomalies" observed by other agents turned out to be explained by a structural
first-obstruction:
- the first syzygy (Cramer, Euler);
- the first differential (`T_1`);
- the first eigenvalue (the extreme Kloosterman sum).

**Open problems.**
- Prove HYP-9315(a), the cyclicity of `H_1`. Compare the Calinescu-Lepowsky-Milas presentations of principal
  subspaces.
- Prove the order-`k` / Li-degree-`(2k+1)` correspondence. This needs the length-`(N-k)` parts of `v`, and
  `T_2` is not universal.
- Extend the spectral sequence to admissible `sl_2` (HYP-9293), `W_3` (HYP-9297) and `N = 2` (HYP-9295).
- Find a closed form for the Cramer weight `w_1` of curves.
