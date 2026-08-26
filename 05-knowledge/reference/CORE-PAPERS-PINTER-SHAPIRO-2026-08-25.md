# Pintér four-flows and Shapiro double-well instantons

> **Freshness:** primary arXiv records and full texts checked 2026-08-25.
> Both sources are fresh preprints. Imported claims are **CITED**, not
> independently proved by this source ledger.

## Pintér -- nowhere-zero four-flows below a proper Petersen minor

- **Primary/status:** József Pintér,
  [*Nowhere-zero 4-flows in graphs excluding a proper minor of the Petersen
  graph*, arXiv:2607.22267v2](https://arxiv.org/abs/2607.22267v2), submitted
  2026-07-24 and revised 2026-08-09. **PREPRINT v2.**
- **Exact imported theorem:** let `P` be the Petersen graph, `e` an edge, and
  `Q=P/e`. Every finite bridgeless `Q`-minor-free multigraph has a nowhere-zero
  `Z_2^2`-flow and hence a nowhere-zero integer four-flow (Theorem 1.1).
  Combining this with the cited Thomas--Thomson theorem for `P-e` gives:
  for every proper minor `R` of `P`, every finite bridgeless `R`-minor-free
  graph has a nowhere-zero four-flow (Corollary 1.2).
- **Proof architecture:** a minor-minimal no-flow graph is reduced, using
  cited obstruction lemmas, to a girth-five almost-four-connected core. The
  Thomas--Thomson structure theorem leaves Triplex, Petersen, Dodecahedron and
  Basket bases. Triplex, Petersen and Basket contain `Q` directly. For the
  exceptional Dodecahedron, the Norin--Thomas jump/cross theorem reduces the
  nonplanar extensions to three noncofacial-pair orbits and one facial-cross
  orbit; explicit deletion/contraction tables put `Q` in every representative.
  The planar branch is closed by four-colour duality.
- **Repo-derived exact bridge:** in a loopless cubic graph, a nowhere-zero
  `Z_2^2`-flow labels the three incident edges at every vertex by the three
  distinct nonzero colours. The matrix

  ```text
  C=[[0,-1],[1,1]],             C^3=-I,             C^6=I
  ```

  reduces modulo two to a three-cycle on those colours. The normalized
  rational quadratic cycle

  ```text
  -2 -> 1 -> -1/2 -> -2
  ```

  is the same abstract `C_3`-set. Therefore Pintér's cubic consequence can be
  restated as existence of a locally bijective edge labelling by that cycle
  palette. This restatement retains the all-distinct cubic-vertex condition;
  the multiplicative product-one shadow alone is not equivalent.
- **Typed loss:** the palette bijection destroys graph incidence and additive
  conservation on one side, and rational coordinates/quadratic interpolation
  on the other. Outside cubic degree, “all three colours once” no longer
  captures the full `Z_2^2` flow law.
- **Does not prove:** Tutte's Petersen-minor four-flow conjecture, a
  characterization of all no-flow graphs, sufficiency of containing `Q`, a
  planar-Jacobian statement, or a transfer from graph minors to decorated
  resolution dual graphs. In particular, `Q` itself has four-flows; it is not
  a no-flow obstruction merely because every no-flow graph must contain it.

## Shapiro -- Poisson-distributed double-well instanton counts

- **Primary/status:** Jacob Shapiro,
  [*Instantons in a Double-Well are Poisson Distributed*,
  arXiv:2608.23342v1](https://arxiv.org/abs/2608.23342), submitted
  2026-08-24. **PREPRINT v1.**
- **Setting/hypotheses:** `H_lambda=-Delta+lambda^2V` on `L^2(R^n)`, with a
  reflection-symmetric nonnegative `C^2` potential having exactly two
  nondegenerate global minima, positive lower bound at infinity, a global
  low-energy spectral gap, and the stated one-well Dirichlet ground-state,
  gap and normalization assumptions. The two lowest global levels are
  discrete below the threshold `lambda^2 c_infinity`, with `E_0` carried by a
  reflection-even state and `E_1` by a reflection-odd state, as explicitly
  assumed in the paper.
- **Exact imported result:** after decomposing the normalized localized
  Feynman--Kac trace by the number `k` of alternating passages between
  shrinking well neighbourhoods and defining the hopping coefficient
  `rho_lambda`, Theorem 2.1 proves, for each fixed `N` and `k`,

  ```text
  P_(lambda,N/|rho_lambda|){k} -> exp(-N)N^k/k!.
  ```

  It also supplies a summable uniform tail `C_N 2^(-k)`. Corollary 2.2 gives
  `E_1-E_0=2|rho_lambda|(1+o(1))`; the action estimate identifies the
  exponential scale of `rho_lambda`.
- **Proof architecture:** exact stopping-time sectorization and strong-Markov
  concatenation reduce separated passages to convolution powers of a
  normalized one-passage kernel. One-well spectral projection extracts a
  rank-one term. Three defect controls -- integrated mass, first moment and a
  pointwise convolution bound -- suppress excited and clustered passages.
  The all-good word produces the simplex volume `beta^k/k!`; a generating
  function supplies the uniform tail needed to sum over all `k`.
- **Imported methodological role:** fixed-complexity factorization and a
  theorem uniform over all complexities are separate obligations. A planar
  Jacobian program that closes each fixed weight, jet or monodromy-word length
  still needs a summable/global algebraic majorant before it can take their
  union.
- **Central-sign analogy:** passage parity is the exact `Z/2` character used
  to distinguish even and odd spectral sectors. For the repo-derived lift
  `B^3=-I`, projective exponent modulo three plus parity reconstructs exponent
  modulo six. The paper supplies a clean model of the missing sign sidecar,
  not a theorem about the deterministic algebraic lift.
- **Typed loss:** replacing the passage kernel by its rank-one term forgets
  incoming boundary point, residence time, excited-state content and clustered
  passage geometry. The three defect norms and uniform-in-`k` bound are the
  required sidecars. By contrast, THM-2230's planar target-shear quotient is
  exact; Shapiro's asymptotic quotient licenses no coarser Keller quotient.
- **Does not prove:** a Poisson point-process theorem for passage locations,
  a polynomial-Keller spectral gap, positivity for monodromy transfer
  operators, a statement about rational dynamics, `3:4:5`, `63`, or `JC(2)`.

## Current consumers and firewall

- [THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md)
  proves the arithmetic/order-six extension independently of either paper;
  its graph/semiclassical comparisons are explicitly analogy-only.
- The session synthesis proposes two experiments only: decorated minor-orbit
  certificates for finite Keller boundary passports, inspired by Pintér; and
  labelled monodromy transfer operators with an ordered-word/commutator
  sidecar, inspired by Shapiro. Neither is an imported consequence.
