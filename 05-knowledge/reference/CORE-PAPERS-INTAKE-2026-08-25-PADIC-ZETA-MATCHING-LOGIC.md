# 2026-08-25 intake: p-adic zeta draft and basic matching logic

**Audit date:** 2026-08-25.  This file separates source claims, published
baseline, replayed finite certificates, and new in-repo deductions.  A source
being recent or AI-assisted is provenance, not evidence for or against its
mathematics.

## 1. Long -- hybrid arithmetic holonomy for p-adic zeta values

### Immutable pin and status

- **Primary repository:**
  [`octonion/p-adic-zeta-irrationality`](https://github.com/octonion/p-adic-zeta-irrationality).
- **Audited commit:**
  [`b46a1770901551961710e155d775aae7c5ea39e7`](https://github.com/octonion/p-adic-zeta-irrationality/commit/b46a1770901551961710e155d775aae7c5ea39e7).
- **Status:** **AUTHOR-CLAIMED / UNREFEREED RESEARCH DRAFT.**  At audit time
  the repository was untagged, one day old, had no CI or formalization, and
  explicitly requested specialist review.
- **Numerical substatus:** the final real interval margins are
  **FINITE-EXACT / independently replayed**, conditional only in the sense
  that they are the final numerical implication of the manuscript's proposed
  geometric and adelic formulas.  Replaying them does not validate those
  formulas.

The draft uses Calegari's normalization of the Kubota--Leopoldt values and
claims the following `22` singleton irrationalities:

```text
zeta_2(s), s=3,5,...,29;                  14 cells
zeta_3(s), s=3,5,...,11;                   5 cells
zeta_5(3), zeta_5(5), zeta_7(3).           3 cells
```

The mechanism is a high-dimensional Hermite--Pade/Apéry analogue: raw
Eichler rows and Dwork relations, a BGG lift, small-prime Hasse-kernel
determinant saving, large-prime de Rham-torsor product digits and no-backflow,
Bost/CDT slopes, analytic continuation, and modular Jensen collision energy.
It does not supply a scalar Apéry recurrence, a transcendence theorem, an
algebraic-independence theorem, or an irrationality measure for the new cells.

### Published baseline and quantifier boundary

- [Calegari, *Irrationality of certain p-adic periods for small
  p*](https://arxiv.org/abs/math/0408214) gives the classical singleton
  controls `zeta_2(3)` and `zeta_3(3)` in the relevant normalization.
- [Lai--Lupu--Sprang, *On the irrationality of certain p-adic zeta
  values*](https://arxiv.org/abs/2505.23088) proves existential
  irrationality in finite packets.  An “at least one” conclusion does not
  imply any named singleton.
- [Lai--Sprang--Zudilin, *A note on the irrationality of
  zeta_2(5)*](https://arxiv.org/abs/2505.05005) supplies a later singleton
  control by an Apéry-type construction.
- [Calegari--Dimitrov--Tang, *Arithmetic holonomy bounds and effective
  Diophantine approximation*](https://arxiv.org/abs/2510.04156) is one of the
  slope/height interfaces used by the draft; its applicability to the new
  vector-bundle source is itself part of the specialist-audit boundary.

Accordingly, cells recorded as `OPEN` in the prior published-literature
ledger are now more precisely **OPEN in the audited published baseline;
AUTHOR-CLAIMED by the pinned unrefereed draft**.  They are not promoted to
`PROVED` here.

### Exact repo consumers

- [THM-4089](../../01-canon/theorems/THM-4089-hybrid-padic-zeta-margin-optimization-and-next-case-obstruction.md)
  treats the pinned margin formula as a self-contained real function.  It
  proves unique global optimizers in both continuous variables and gives
  strictly negative global upper bounds for the four adjacent cells
  `(2,31),(3,13),(5,7),(7,5)`.  This is an obstruction to retuning that one
  formula, not a no-go for other templates or methods.
- [THM-4093](../../01-canon/theorems/THM-4093-rational-edge-diagonal-gauge-and-padic-tournament-zeta-tangent.md)
  proves that vertex-ratio edge weights are diagonal gauge and identifies the
  directed-triangle `p`-adic tangent of adjacency zeta.  Its zeta is a graph
  determinant, not a Kubota--Leopoldt value.
- The literal denominator connection to THM-4056 is the clearing clock
  `L_N^e`: it preserves valuation depth but loses numerator cancellation,
  distinguished-prime deletion, and analytic decay.

### What still requires specialist audit

The first load-bearing interfaces are the BGG normalization and genuine
multiplier; descent and saturation of the global source; Hasse degree/kernel
and adelic lattice passage; fixed-bundle Bost/CDT transfer; torsor
product-digit no-backflow; continuation radii at `p=5,7`; and primitive-row
Jensen multiplicities.  The smallest displayed margin is about `0.1318`, so
constant and normalization errors are not automatically harmless.

## 2. Chen--Rosu -- completeness and incompleteness of basic matching logic

### Identity and status

- **Primary:** Xiaohong Chen and Grigore Roșu,
  [*Completeness and incompleteness of basic matching
  logic*](https://arxiv.org/abs/2608.13306v1), arXiv:2608.13306v1,
  submitted 2026-08-13.  **PREPRINT v1.**
- The supplied URL is a logic paper, not a number-theory paper.  Its source
  bundle contains one TeX file and no mechanization; the paper lists
  mechanization as future work.

### Exact imported results

In one-sorted, fixpoint-free, definedness-free basic matching logic over an
arbitrary finitary signature, with no set variables, the paper proves for
closed `Gamma,phi`

```text
Gamma |= phi    iff    Delta_Gamma |=_loc phi    iff    Gamma |- phi.
```

The model-theoretic mechanism takes the largest backward-closed core and, for
a nonempty proper core, replaces it by a double cover `C x {0,1}`.  The proof
imports the displayed system's soundness and strong local completeness.

With least fixpoints, the paper proves non-recursive-enumerability of validity
already over one unary and two binary symbols, by a Diophantine reduction.
This is a non-r.e. result, not an upper complexity classification and not a
barrier to calculi whose proof relation is itself non-r.e.

For many sorts the paper gives a satisfiable three-sort counterexample to its
standard Hilbert system.  This is failure of that calculus, not
nonaxiomatizability; fixpoint-free global consequence still translates to
first-order consequence.  Its nominal obstruction likewise applies only to
localization-respecting calculi.

### Exact repo sharpening

[THM-4090](../../01-canon/theorems/THM-4090-two-sort-matching-logic-global-completeness-obstruction.md)
shows that the third sort is unnecessary.  With `f:b->a`,

```text
Gamma={forall x:b forall y:b. f(x and y)},
phi=forall x:b. x,
```

one has satisfiable `Gamma`, semantic entailment, and nonderivability.  The
hypothesis forces `|M_b|=1`, while the feed graph has `b->a` and no path
`a->b`.  The paper's one-sort completeness makes two the exact minimal sort
count existentially for this fragment and Hilbert system.

The feed carrier is an intrinsic directed graph, not generally a tournament.
Forcing a total orientation can invent reachability and invalidate the
localization proof.

## 3. Shared research lesson

The two sources meet at a typed boundary, not at a common theorem.

- Matching logic separates **local reachability** from **global totality**;
  the missing sort-feed path blocks proof transport.
- The p-adic draft separates **finite numerical margin** from **global
  arithmetic interfaces**; a correct last inequality cannot backfill a
  missing BGG, lattice, continuation, or Jensen map.
- THM-4093 supplies the tournament version: a multiplicative vertex potential
  is diagonal gauge, so closed-walk data see none of its endpoint arithmetic.

In each case the decisive question is not whether two expressions look alike,
but which map preserves the target predicate and which coordinate was lost.

## 4. Does not prove

This intake proves no new Kubota--Leopoldt irrationality, transcendence, or
algebraic independence; no p-adic analogue of matching-logic completeness;
no tournament classification of rational, irrational, or transcendental
numbers; and no LRC(14), planar Jacobian, Mahler `3/2`, or Rule-30 prize.
