# August 26, 2026 intake: p-adic-zeta density and critical percolation

## 1. Long -- density-one irrationality of p-adic zeta values

- **Source:** Christopher D. Long,
  [*Multiweight arithmetic holonomy and density-one irrationality of p-adic
  zeta values*](https://github.com/octonion/p-adic-zeta-density), commit
  [`4c87bcdf`](https://github.com/octonion/p-adic-zeta-density/commit/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57).
- **Status:** **PREPRINT CLAIM / UNDER SPECIALIST AUDIT**. No density or
  irrationality claim is promoted to repo canon.
- **Claim:** quantitative density-one irrationality among positive odd
  arguments for every fixed prime, with simultaneous finite-packet bounds.
- **Current audit:** the proposed `u-f` specialization counterexample is
  correct for a chosen-section map but no such map appears in the manuscript.
  Proposition 6.2 is written on a universal flat torsor chart; Proposition
  12.3 is written in a universal finite-free coefficient lattice. Their
  exact completed-ring/module diagrams should be stated explicitly. See the
  [full specialization audit](p-adic-zeta-density-specialization-kernel-audit-20260826.md)
  and
  [THM-4255](../../01-canon/theorems/THM-4255-specialization-kernel-and-transverse-hasse-jet-repair.md).
- **Reproducibility:** the one tracked certificate replays exactly, but seven
  cited artifacts are missing and the manuscript's printed hashes for the
  tracked pair are stale.

The exact new repo contribution is independent of the preprint claims:
single-arc evaluation has principal kernel, finite Hasse jets have exact
power kernels, a universal-slope pencil is injective, and finite-box
Kronecker substitution has a sharp degree boundary.

## 2. Cerf -- critical percolation and finite-cluster tails

- **Source:** Raphaël Cerf,
  [*theta(p_c,Z^d)=0?*](https://arxiv.org/abs/2608.23661v1), arXiv
  `2608.23661v1`, submitted August 24, 2026; manuscript dated August 22, 2026.
- **Frozen PDF SHA-256:**
  `37c5759165df93e14da20db1de9f0117d67e0d16a1d1dd385ce16d85c0422a97`.
- **Status:** **PREPRINT CLAIM / UNDER SPECIALIST AUDIT**. Only the theorem
  statements and proof guide were checked in this intake; the 708-page proof
  was not independently verified.
- **Identity check:** this is a percolation paper, not a p-adic-zeta or
  Cartier paper. Its Propositions 6.2/12.3 are unrelated to the GitHub
  manuscript.

For Bernoulli site percolation on `Z^d`, `d>=3`, let `C(0)` be the cluster of
the origin. Cerf's Theorem 1.2 states that at the critical point either

```text
theta(p_c,Z^d)=0,
```

or

```text
liminf_(n->infinity) (1/log n)
  log log(1/P_(p_c)(n<=|C(0)|<infinity)) = 0.            (1)
```

Thus, if critical percolation exists, the finite-cluster tail is not
stretched exponential. The stronger Theorem 1.3 states that, for
`p in (0,1)`,

```text
theta(p,Z^d)>0 and the liminf in (1) is positive  =>  p>p_c. (2)
```

The question mark in the title is load-bearing: the paper does **not** claim
to prove `theta(p_c,Z^d)=0` for all `d>=3`.

The proof program uses one-step renormalization. The extra tail hypothesis
rules out negative higher pivots which are far apart; multiple intertwined
explorations control local higher pivots. The author identifies control of
far-apart/nonlocal higher pivots as the missing step toward the full
conjecture.

## 3. Typed connections and firewalls

| source | target/map | preserved | lost / required sidecar |
|---|---|---|---|
| universal torsor identity | p-adic Cartier coefficient algebra | literal zero and pole grade before scalarization | completed-chart/module diagram |
| one chosen formal arc | `ev_g:A[u][[f]]->A[[f]]` | scalar value on that arc | kernel `(u-g)`; retain Hasse jets or another section family |
| planar-JC `DE` edge | `w->T` in THM-4120 | edge leading value and fixed weight | `(w-T)` normal class; THM-4120 retains order-one/two re-entry |
| Cerf local-pivot reduction | locality-filtered exploration | local pivotal geometry | far-apart pivot configurations; no p-adic theorem transfer |

The shared lesson is a proof-obligation split, not a mathematical equivalence:
local control becomes global only after the nonlocal fibre has been ruled out
or retained. Cerf's higher pivots are Boolean/probabilistic derivatives;
Long's pole grades are `ell`-adic associated-graded terms. The common word
`pivot` licenses no identification of their objects, measures, or theorems.

For the planar Jacobian frontier the connection is exact rather than
analogical: THM-4120 already computes the principal specialization kernel and
its normal re-entry. THM-4255 turns any future bounded-degree edge response
into a finite sharp observer problem using all normal Hasse jets or enough
lawful fibres. This does not construct a Keller pair or prove `JC(2)`.

## 4. Research targets exposed by the pair

1. Write the universal coefficient-module diagram needed by the general-curve
   Cartier step and test one explicit torsor/pole-grade example.
2. Determine the smallest transverse-jet family that retains the actual
   bounded fibre degree in each p-adic block.
3. On the planar-JC `DE` layers, derive the exact `w`-degree cap and replace a
   scalar edge value by the sharp Hasse-jet observer.
4. Treat Cerf's far-apart higher-pivot problem only as an independent
   percolation frontier; do not turn locality rhetoric into an LRC, Cartier,
   or Jacobian implication without a typed map and preserved predicate.
