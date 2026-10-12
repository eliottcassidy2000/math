# Four logic lenses on the planar Jacobian conjecture, Collatz and Hadwiger–Nelson

**Session:** mac-mini-2026-10-12-logic4 (isolated worktree `math-wt-logic4`, branch `codex/session-logic4`).
**Status:** PROVED (THM-4620 dictionary, THM-4622 choice chain, THM-4623 order/endomorphism/Koopman
dictionary, THM-4624 identifications and classification) + FINITE-EXACT (THM-4621's exact fibre
computations and controls; THM-4624's nauty census to `q = 81`) + CITED (Bavula; Gao; Nagata–Nowicki;
Beckman–Quarles; Kiermaier–Kurz) + OPEN (HYP-9255, HYP-9256, HYP-9257). **Audit:** pending at the time
of the first checkpoint; see section 7.

## 0. The four sources and the lens each supplies

| source | content | lens used here |
|---|---|---|
| Pinsky, *Frucht's theorem without choice*, arXiv:2305.11382 | every group is Aut of a graph holds in ZF but fails in permutation models; Babai's strengthening fails too | automorphism groups and rigidity as a predicate; *support* (a graph realizing a group must be finitely supported) |
| Xu–Ji, *Every countable meet-continuous lattice is Scott sober*, arXiv:2609.19612 | countability allows successive finite choices; the non-sober cardinal spectrum is upward closed | order topology of orbit spaces: dcpo / sober / Scott; countable = choice-free |
| Fernandes–Mezabarba–Rodrigues, *Existence of bases implies AC, foundation-free*, arXiv:2609.20140 | "every vector space has a basis" ⟹ AC without Foundation, also with atoms | where choice enters colouring and dimension statements; absoluteness of countable certificates |
| Chini, *Generic vector space endomorphisms V: the o-minimal case*, arXiv:2609.20184 | theories `T_θ ∪ {Σ∩Ker ρ[θ] = Σ∩Ker η[θ]}`; model companions; TP₂/SOP/NATP, distal non-o-minimal | the kernel lattice of an endomorphism as the invariant to compute: dual fields of Keller maps, the Collatz Koopman operator |

## 1. Inheritance pass

* **Closest proved mechanisms.** THM-4562 (slice lemma: a locally nilpotent dual field splits
  `A^n = X_j × A^1`; THM-1300 is fully non-slice), THM-4523 (the Collatz graph is rigid, by 3-adic
  separation of backward trees), THM-4558 (residue colourings of field planes, `χ(G(L)) ≤ κ(q)`).
* **Canonical hostile examples.** THM-1300's map (every claim about Keller modules must survive it);
  the padding `F' × id` (kills "surjective ⟹ invertible" in dimension `≥ 4`); the Clebsch graph
  (kills "every unit-circle automorphism is affine" at `q = 4`).
* **Corrected near misses.** MISTAKE-205 (the Alpöge attribution of THM-1300 was ruled a
  hallucination by the owner; now a human-authored arXiv paper, Gao 2608.00222, makes the same
  attribution — reported in THM-4621, not adopted); MISTAKE-586 (KNOWN theorems typed as new
  packaging — every KNOWN ingredient below is labelled).
* **Least-used sidecars.** The `κ(q)` table of THM-4558 (now read through automorphism groups);
  the dual fields' *cokernel* (THM-4562 used only kernels).

## 2. Concept board

1. `V_j = ∂/∂F_j`, the dual fields of a Keller map, as endomorphisms of `K[x]`.
2. The kernel lattice `ρ ↦ Ker ρ(V_j)` (Chini) and the cokernel `coker V_j = H^1_dR(K[x]/K[F_(≠j)])`.
3. Punctured-line fibres: a degree-`d` function on `P^1` with its critical points deleted.
4. Frames and bases in `F_(q^2)` over a subfield: when the unit circle is "too small".
5. The reachability preorder of a self-map and its Alexandrov/Scott topology.
6. Finite support: a certificate (5-chromatic graph, residue colouring, collision) lives in `R_alg`
   or `Q`, hence is absolute.
7. Endomorphisms versus automorphisms: rigidity (Frucht) is about the group, not the monoid.

How they moved each other: (2) applied to (1) gave the `t`-primary lattice and the surjectivity
criterion; (3) is what (2) sees on THM-1300 and is Gao's tangent sweep in fibre coordinates;
(4) is the finite version of Beckman–Quarles that (7) asked for; (5) was asked by Xu–Ji and answered
by "no divergence"; (6) settled which parts of Hadwiger–Nelson are choice-sensitive.

## 3. Connection ledgers

| | source | target | map | preserved | lost | sidecar | decisive test |
|---|---|---|---|---|---|---|---|
| A | Chini's kernel identities | Keller dual fields | `θ := V_j` on `V = K[x]` | `Ker ρ(θ)` for all `ρ`; surjectivity sentence | model companion, o-minimal base, dp-rank | THM-4562's dual fields; Nagata–Nowicki | eigenvector lemma (THM-4620(2)); residues on fibres (THM-4621(8)) |
| B | Xu–Ji sobriety | orbit space of `T` | reachability preorder, `T_0` quotient | directed sets ↔ orbits that merge | the lattice structure (orbit spaces are posets, not lattices) | THM-4523 (rigid); THM-4581 (coalescence) | divergent orbit ⟺ non-principal irreducible closed set (THM-4623(1)) |
| C | bases ⟹ AC | chromatic number of the plane | finite certificates ⟹ `R_alg` | lower bounds, residue colourings, `κ(q)` | the uncountable plane itself | de Bruijn–Erdős; Tarski transfer | `χ(R^2) = χ_alg` in ZFC, `≥` in ZF (THM-4622) |
| D | Frucht without choice | unit-distance graphs | `Aut` as a predicate | Beckman–Quarles over `R`; affine rigidity over `F_q` | the permutation-model failure (no analogue) | nauty; Kiermaier–Kurz | census to `q = 81`; three exceptions explained (THM-4624) |

## 4. Results

### 4.1 THM-4620 — Keller dual-field modules (PROVED; CITED parts labelled)

For a Keller map `F` and dual field `V_j`: kernel ladder `Ker V_j^m = ⊕_(i<m) F_j^i R_j`, locally
nilpotent part `R_j[F_j]`; `Ker ρ(V_j) = Ker V_j^(ord_0 ρ)` for every `ρ` (no eigenvector with nonzero
eigenvalue, even in `K̄(x)`, by a minimal-polynomial-plus-ODE argument); `coker V_j` is the relative de
Rham cohomology of `F_(≠j)`, holonomic, of generic rank `Σ(2γ_i + r_i − 1)` over the geometric
components of the generic fibre; a surjective `V_j` forces all components to be lines and
`deg F = [L_j : K(F_(≠j))]`; **in the plane, `(f, g)` is an automorphism iff `D_f = {f, ·}` is surjective
iff `H^1_dR(A^2/A^1; f) = 0`** (KNOWN: iff `D_f` locally nilpotent). Bavula's holonomicity/length/"generated
by 1" is CITED. The `n ≥ 4` failure of the surjectivity criterion (padding) is recorded; `n = 3` is
HYP-9255.

### 4.2 THM-4621 — THM-1300's coordinate-pair fibrations (PROVED by exact computation + elementary arguments)

All three generic fibres are punctured lines (`P^1` minus 3, 5, 3 points) on which the remaining
coordinate is a degree-three function whose critical points are exactly the punctures:
`F_3 = −(y − b)(2y^2 − by + 18a − b^2)/(27a^2)` on `A^1_y ∖ {y^2 − by + 3a = 0}`;
`F_2 = 2aX + 2/X − c/X^2` on `P^1_X ∖ {0, ∞, aX^3 − X + c = 0}`;
`F_1 = s(cs^2 − 2s + b)/2` on `P^1_s ∖ {∞, 3cs^2 − 4s + b = 0}`. The fibre cubic is Gao's tangency
polynomial under `w = 3ac/(b − y)`; the discriminant is `−4c_3/(729a^6)` with Gao's `c_3`
(independent re-derivation of his Theorem 3.4: fibre sizes 3/1/0, Jelonek set `{c_3 = 0}`, image
`A^3 ∖ C`); Galois group `S_3`; `R_j = K[F_(≠j)]`, nilpotent parts `K[F]`, `x, y, z ∉ K[F]`; cokernel
ranks `2, 4, 2` with explicit witnesses `x^2 ∉ V_3K[x]`, `x ∉ V_2K[x]`, `x^2 ∉ V_1K[x]`; `A_3`-length `≥ 3`.

### 4.3 THM-4622 — choice and the chromatic number of the plane (KNOWN in substance; PROVED as written)

Countable de Bruijn–Erdős holds in ZF (least-child König); `χ_alg = χ(UD(R_alg^2))` is arithmetical
and absolute with `5 ≤ χ_alg ≤ χ(R^2) ≤ 7` in ZF and `χ(R^2) = χ_alg` in ZFC; residue colourings and
countable field planes are choice-free; `Aut(UD(R^2)) = Isom(R^2)` (KNOWN) versus `Aut(UD(Q^2)) ⊇ Sym(ω)`
(Pythagorean denominators: the points `(3^(−k), 0)` lie in distinct components).

### 4.4 THM-4623 — Collatz under order and module lenses (PROVED, elementary)

For any self-map: `P` dcpo ⟺ `Alex(P)` sober ⟺ directed sets have greatest elements ⟺ no divergent
orbit; Collatz ⟺ `P` is a dcpo with a top. `|End(N, T)| = 2^(ℵ_0)` (collapse one ancestor tree onto
the other preimage's tree) while `Aut = 1` (THM-4523). The Koopman operator on locally constant
functions on `Z_2` has `Ker ρ(θ) = K·1` if `ρ(1) = 0` and `0` otherwise.

### 4.5 THM-4624 — finite Beckman–Quarles (FINITE-EXACT census + PROVED structure)

`Aut(G_q) = Aff(q)` (order `q^2(q+1)·2e`) for every prime power `q ≤ 81` except `q = 4, 8, 9`
(ratios 6, 24, 6). The unit circle is a union of subfield lines iff `p^f − 1 | q + 1`, possible only for
`(p, f) ∈ {(2,1), (2,2), (3,1)}`, and is a basis/frame exactly for `q ∈ {2, 3, 4, 8, 9}`:
`K_4`, `K_3 □ K_3`, Clebsch, `H(3,4)`, ternary five-frame (automorphism groups `2^4 ⋊ S_5`, `S_4 ≀ S_3`,
`F_3^4 ⋊ (S_5 × C_2)`). `κ(9) = 3` via a nowhere-zero functional; `κ(4) = κ(8) = 4`.

### 4.6 THM-4625 — the whole tangent-sweep family (PROVED)

For every Keller map of Gao's shape with `deg p = d ≥ 2` (THM-1300 is `d = 2`, Gao's `G` is `d = 3`,
Gallagher's family by Gao's account): the `(F_2, F_3)`-fibres are `P^1` minus `d+1` points carrying
`F_1` as a degree-`(d+1)` polynomial minus its `d` critical points (`F_1' = (cb − p)/(2c^2)`), the
`(F_1, F_3)`-fibres are `P^1` minus `d+3` points carrying `F_2` of degree `d+1`
(`F_2' = −2(c^2a − q)/(cw^2)`), and the `(F_1, F_2)`-fibres are birational to `t^2 = b^2 − 8aP(w)/w^2`
(`P = ∫p`; discriminant identity `disc = p^2(b^2w^2 − 8a(wp − 2q))`), irreducible hyperelliptic of
genus `⌊(d−2)/2⌋`. So `deg F = d+1`, all rings of constants are `K[F_(≠j)]`, no dual field is
surjective or locally nilpotent, no coordinate-pair fibre is a union of lines (HYP-9255 holds on the
family), and the Galois group is `S_(d+1)` whenever the sweep curve is parametrized birationally
(`S_3` for THM-1300, `S_4` for `G`).

### 4.7 Hypotheses

* HYP-9255 — surjective dual field ⟹ automorphism in dimension three (refuted for `n ≥ 4` by padding; holds in fibre form on the whole tangent-sweep family by THM-4625).
* HYP-9256 — affine rigidity of `G_q` for all `q ∉ {4, 8, 9}`.
* HYP-9257 — `length_(A_3) V_F = 4` for THM-1300.

## 5. Computations and reproduction

    cd 04-computation/experiments/logic4_20261012
    python3 thm1300_fibrations.py      # parametrizations, F_3|fibre, residues, discriminant = Gao's c_3, image complement, iterates
    python3 thm1300_module_checks.py   # Ker V_3 in degree <= 12 (dim 4), tangency polynomial, x,y,z outside K[F]
    python3 gq_automorphisms.py        # nauty census q <= 81 (needs dreadnaut on PATH)
    python3 gq_exceptional.py          # Clebsch / H(3,4) / ternary frame isomorphisms, chi by SAT, frame colouring
    python3 sweep_family_identities.py # derivative identities for generic p, d = 2..6
    python3 sweep_family_discriminant.py # discriminant identity, wp - 2q = int p, generic genus, S_4 for Gao's G
    python3 gao_degree4_fibrations.py  # Gao's degree-four map: three fibrations (its 'genus 2' line is a wrong sqf_part reading; see THM-4625)
    python3 gq_census_extension.py     # nauty census for 83 <= q <= 128 (HYP-9256 evidence; background)
    python3 kappa_sat.py 13 5 ; python3 kappa_sat.py 17 5   # CaDiCaL attempts at kappa(13), kappa(17) (background)

Outputs are committed beside the scripts (`*.out`). Hostile controls: (i) `dim Ker V_3` on polynomials
of degree `≤ 12` is `4 = #{1, F_2, F_2^2, F_1}` (a larger kernel would refute `R_3 = K[F_1, F_2]`);
(ii) the discriminant of the fibre cubic reproduces Gao's `c_3` exactly; (iii) the stabiliser orbit
counts in the census equal `1 + #(Frobenius orbits on F_q^×)` in every affine case; (iv) the
exceptional graphs' automorphism groups are bounded below by explicit subgroups and above by nauty,
and the two bounds agree; (v) `det JF = −2` and the integer triple collision were re-verified.

## 6. Literature placed

* Bavula (arXiv:2112.03177; C. R. Math. 2024): the twisted module is holonomic, of finite length with
  an explicit bound, cyclic; JC ⟺ generated by `1`. (THM-4620(6) is CITED to this.)
* van den Essen (Pacific J. Math. 151, 1991): bounded fibre cardinality of Keller maps via D-modules.
* Zhang (arXiv:2609.10180, Sept 2026): characteristic-`p` Keller endomorphisms of `D` preserve holonomicity.
* Gao (arXiv:2608.00222, July 2026; "Claude Fable 5 assisted"): tangent sweep; Theorem 3.4 is the
  fibre stratification re-derived in THM-4621. Gao attributes the dimension-three counterexample to
  Alpöge (July 19, 2026), the infinite family to Gallagher, the geometric explanation to Speyer;
  the repo's provenance ruling (THM-1300 owner correction, MISTAKE-205) is not changed here.
* Nguyen Van Chau (arXiv:1502.00328): Keller pairs over `Q` and integral points (not used).
* Kiermaier–Kurz (Discrete Math. 309, 2009): affine rigidity of the integral-distance graph on `F_q^2`.
* Shelah–Soifer (JCTA 103, 2003): axiom-dependence of the chromatic number of a distance graph on `R`.
* Beckman–Quarles (Proc. AMS 4, 1953).

## 7. Audit

Pending at the first checkpoint. Two independent adversarial audits were requested (A: THM-4620,
THM-4621, HYP-9255, HYP-9257; B: THM-4622, THM-4623, THM-4624, HYP-9256); their findings and the
fixes applied are recorded below when they return.

## 8. Next obligations

1. Test HYP-9255 on Gao's degree-four map `G` (THM-4621's elimination method; count components of
   the generic coordinate-pair fibres).
2. Extend the nauty census to `q ≤ 256` (HYP-9256) and attempt the common-neighbourhood proof route.
3. The cusp-curve local computation for HYP-9257.
4. For the plane: the surjectivity criterion says a counterexample `f` has `b_1(generic fibre) ≥ 1`;
   the Gauss–Manin system `coker D_f` is a holonomic `A_1`-module of that rank — its regular
   singularities at the atypical values are the next object to type (what does `D_g` do to it?).
5. Is `μ_j = 0` enough for `coker V_j = 0` (HYP-9255's sub-question)?

## 9. META-PATTERNS candidate

**Lens transfer from a logic paper.** Trigger: a paper classifies structures by a lattice of
first-order-definable subspaces (kernels, fixed points) or by an axiom's strength. Move: compute
that lattice for the repo's concrete endomorphism (dual field, Koopman operator) and ask which
first-order sentence separates the standard model from the hostile one; for axiom-strength papers,
ask which certificates have finite support. Counterindication: when the object is uncountable and
the predicate is not first-order (the Collatz conjecture itself; `χ(R^2)` in ZF), the lens only
reformulates. Evidence from distinct threads: THM-4620/4621 (Jacobian), THM-4623 (Collatz), THM-4622
(Hadwiger–Nelson).
