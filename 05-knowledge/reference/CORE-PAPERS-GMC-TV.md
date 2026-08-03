# Core papers: GMC counterexamples and TV homogenization (opened 2026-08-03)

Overflow entry for [`CORE-PAPERS.md`](CORE-PAPERS.md), which is at its bounded
startup-surface byte budget.  Same contract: what the repository imports, where
it is consumed, and what the source does **not** establish.

## Kontorovich — *TV homogenization inequalities*

- **Primary / freshness:** [arXiv:2601.04079](https://arxiv.org/abs/2601.04079),
  v1 2026-01-07, v2 2026-02-08, **v3 2026-02-25** (the version pinned here).
  math.PR.  **PREPRINT; not confirmed peer reviewed as of this check
  (2026-08-03).**
- **Imported role:** the block inequality
  `delta_N <= delta_I + delta_J - delta_I delta_J` for the *homogenization*
  map on inhomogeneous Bernoulli product measures, equivalently
  supermultiplicativity of the affinity `1-delta`.  The paper's point is that
  homogenization -- replacing each Bernoulli parameter by the cumulative mean
  -- is **not** a data-processing map, so this is not the usual coupling bound;
  it nevertheless contracts `TV` up to a universal constant, via two-sided
  control of `TV` between Poisson binomials.
- **Repo consumer:**
  [THM-3291](../../01-canon/theorems/THM-3291-two-block-tv-homogenization-rigidity.md),
  which proves the two-single-observation-block case from a box constraint plus
  AM-GM and, going beyond the cited statement, classifies the equality locus.
- **Does not prove (in repo terms):** anything about LRC(14), the Gaussian
  moment lane, or the AMM 12592 biased-coin lane.  THM-3291 records the
  relation to THM-3290 as an explicit **contrast, not a bridge**: `TV` is a
  positive functional and can only be contracted by an orbit average, whereas
  a signed moment functional can be annihilated outright.
- **Scope warning:** THM-3291 proves only the two-block case.  The general
  block statement remains CITED.

## Long — *Small Counterexamples to the Gaussian Moments Conjecture*

- **Primary / freshness:** [arXiv:2607.18186](https://arxiv.org/abs/2607.18186).
  Also carried in the main `CORE-PAPERS.md` GMC section; repeated here only for
  the size comparison that THM-3290 depends on.
- **Imported role / size fact:** GMC is false in every dimension `n>=3`; the
  smallest three-variable example reported there has **degree 4 with five
  terms**, plus a six-term cubic example in four variables.
- **Consequence for repo claims:** the `nu=1` object proved in
  [THM-3290](../../01-canon/theorems/THM-3290-archimedes-flatness-and-the-gmc3-gvc3-counterexample-family.md)
  has degree 12 and 23 terms and is therefore **not minimal**.  No minimality
  or priority claim may be attached to it.

## Unresolved provenance: the supplied GVC(3) object

The owner supplied, on 2026-08-03, the three-variable object

```text
rho=t^2+xy,  A=rho+x^2,  C=y rho^2-2x t^2 rho-x^3 t^2,  P=A C^2,
Delta=4 d_x d_y+d_t^2,  Q=x^2,
Delta^(6m)(P^m)=0,  Delta^(6m+1)(Q P^m)=2^(8m+1)(6m+1)!(2m)!(12m+3)!!/(4m+1)!!
```

attributed to `arXiv:2606.17854`.  **That identifier is wrong.**  It resolves
to Ajwani--Gajjala--Raman--Ray, *Counterexamples to Wegner's Conjecture for
Rectangles* (cs.CG; v1 2026-06-16, v2 2026-07-09), a rectangle piercing/packing
paper containing none of this material.  The object itself is correct --
THM-3290 proves every displayed claim, including the closed form, for all
`m>=1` -- but its source is **UNRESOLVED**.

Rules for later sessions:

1. Do not cite `arXiv:2606.17854` for anything in the Gaussian-moment lane.
2. Do not attach a priority or first-discovery claim to the `nu=1` object.
   THM-3290 claims only its own proof, closed form, threshold, and family.
3. If the true source is later identified, add it here and update THM-3290's
   `external:` block rather than editing the theorem body.

## Standing dictionary (used by THM-3290, worth reusing)

For `rho` a nondegenerate quadratic form in `n` variables with Laplacian
`Delta` and `L(f)=(exp(Delta/2)f)(0)`:

```text
f homogeneous of degree 2k  =>  L(f)=Delta^k f/(2^k k!)=(2k+1)!!_n * <f>,
```

where `<f>` is the mean over the unit sphere of `rho` and the radial factor for
`n=3` is `(2k+1)!!`.  Hence a *Generalized Vanishing Conjecture* statement
about `Delta^j` in `n` variables is exactly a *Gaussian Moments Conjecture*
statement in `n` variables; they are not separate lanes.  Zhao's equivalence
with the Jacobian Conjecture, however, is for `j=1` only.
