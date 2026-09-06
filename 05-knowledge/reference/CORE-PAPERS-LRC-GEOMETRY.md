# LRC geometry inputs: sources, consumers, and boundaries

Current bibliographic sidecar, routed from CORE-PAPERS. The source entries
below are preserved verbatim from checkpoint0d6ff3d1c; their individual
freshness dates and scope restrictions remain in force. This relocation
keeps the startup surface bounded without deleting source detail.

<a id="vaaler"></a>

### Vaaler — *Some extremal functions in Fourier analysis*

- **Primary / freshness:** J. D. Vaaler, *Bulletin of the American Mathematical
  Society (N.S.)* **12** (1985), 183--216, Theorem 19. **PUBLISHED / stable.**
- **Imported role:** supplies the one-dimensional degree-`H` trigonometric
  majorant/minorant sandwich for an interval that THM-2085 tensors with signed
  coordinate defects.
- **Repo consumer:** [THM-2085](../../01-canon/theorems/THM-2085-explicit-height-57-rank-seven-selberg-gate.md).
- **Does not prove:** the relative-Hunter inequality, the signed tensor
  bookkeeping, `H=57`, optimality of that height, or LRC(14). Those are
  repo-derived arguments and constants.

<a id="ungar"></a>

### Ungar — *2N noncollinear points determine at least 2N directions*

- **Primary / freshness:** Peter Ungar, *Journal of Combinatorial Theory,
  Series A* **33** (1982), 343–347,
  [DOI 10.1016/0097-3165(82)90045-0](https://doi.org/10.1016/0097-3165(82)90045-0).
  **PUBLISHED / stable; bibliographic record checked 2026-07-21.**
- **Imported role:** historical alternate lens only.  Applying the
  even-cardinality direction bound to the symmetric signed-column
  configuration can supply a nonradial secant whose perpendicular projection
  has a repeated absolute speed.  The current proof of THM-2053 does **not**
  depend on Ungar: its adjacent-normalized-column construction produces the
  full-support repeat projection elementarily.  The standing lower-dimensional
  LRC citation, not this paper, is the remaining external input to the torus
  floor `M_T>=1/13`.
- **Repo consumers:**
  [THM-2053, rank-two geodesic terminal](../../01-canon/theorems/THM-2053-rank-two-parameter-plane-geodesic-terminal.md),
  [HYP-8846, finite tangent-disk completion](../hypotheses/HYP-8846-lrc14-pointed-plane-transport.md).
- **Does not prove:** LRC(14), the repeat-projection lemma now used in
  THM-2053, the determinant gate, necessity of that gate, or emptiness of any
  tangent disk. The anisotropic estimate and disk identity are separate
  in-repo arguments; disk membership remains only an uncertified case.

<a id="malikiosis"></a>

### Malikiosis--Santos--Schymura — *Linearly-exponential checking is enough...*

- **Primary / freshness:** [arXiv:2411.06903v2](https://arxiv.org/abs/2411.06903),
  published as *Forum of Mathematics, Sigma* **13** (2025), e164,
  [DOI 10.1017/fms.2025.10107](https://doi.org/10.1017/fms.2025.10107).
- **Imported role:** reduces LRC up to `n+1` runners to an explicit finite
  integer-velocity range of order
  `binom(n+1,2)^(n-1) <= n^(2n)`, dramatically improving Tao's earlier
  `n^(O(n^2))` range.  The zonotope formulation is the standard finite-check
  backend used in the repo.
- **Repo consumers:**
  [finite-check feasibility ledger](../../00-navigation/LRC14-FINITE-CHECK-FEASIBILITY-LEDGER-2026-07-19.md),
  [LRC14 frontier](../../00-navigation/LRC14-FRONTIER-2026-07-15.md),
  [THM-599 torus-band theorem](../../01-canon/theorems/THM-599-torus-band-theorem.md).
- **Does not prove:** LRC(14), a feasible small practical search, or the shifted
  theorem unconditionally in every dimension; the paper's shifted statement
  retains its stated Lonely Vector Problem dependency.
