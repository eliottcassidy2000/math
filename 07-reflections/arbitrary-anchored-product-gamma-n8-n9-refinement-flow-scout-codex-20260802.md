# Arbitrary anchored product-Gamma: exact `N=8,9` refinement-flow scout

Status: **FINITE-EXACT / theorem-grade scout, not yet a proved dependency.**

This note independently extends the monomial-fibre refinement mechanism being
packaged after `THM-3112-cycle-weighted-young-subgroup-gap-and-uniform-octopus-boundary.md`.
It neither edits nor anticipates the status of `THM-3115`.

## 1. Inheritance and the live board

The closest proved mechanism is THM-3112's identity

```text
beta_N(S)=sum_(mu partition N) m_mu(S) Lbar_mu,
```

where `Lbar_mu` is the uniform average of `I-P_Hpi` over set partitions of
block-size type `mu`.  If `mu` refines `nu`, the biregular incidence coupling
between the two set-partition orbits and subgroup inclusion give

```text
Lbar_mu <= Lbar_nu.                                             (1)
```

The canonical near miss is THM-3112's negative degree-five coefficient at
`mu=(2,1,1,1)`: coefficientwise positivity in the monomial basis is false.
The least-used relevant sidecar is the exact positive-forest no-go in
`gmc-product-gamma-rank4-forest-holotopy-boundary-codex-20260802.md`: a
successful lift must retain signed intra-fibre cancellation rather than add
another same-sign forest resolution.

The resulting concept board is:

1. Young-subgroup gaps `Lbar_mu`;
2. the integer-partition refinement poset;
3. chamber-Newton coefficients in the two support variables;
4. exact positive transport on that poset;
5. the zero operator `Lbar_(1^N)=0`; and
6. restricted row/hook/two-row carriers as hostile sub-cones.

## 2. Exact normalization

For either THM-3110 invariant and either integer chamber, put

```text
G_mu = [Phi(h_N)m_mu(Q)-Phi(m_mu)h_N(Q)]/base.                  (2)
```

Then the row-normalized central operator satisfies

```text
h_N(Q) D_N = sum_mu G_mu Lbar_mu.                              (3)
```

The normalization in the first display is exact.  For a fixed set partition
`pi` of type `mu`, summing labelled colourings with fibre partition `pi`
produces

```text
[prod_i(mu_i!) prod_r(mult_r(mu)!)/N!] m_mu(S).
```

There are exactly the reciprocal number of such set partitions after
multiplication by `N!`; hence averaging over the type orbit cancels the
factor and leaves precisely `m_mu(S)Lbar_mu`.

At one chamber-Newton slot, write the coefficient vector of `(2)` as `g_mu`.
Whenever a positive type `nu` coarsens a negative type `mu`, an amount `t`
can be replaced by

```text
t(Lbar_nu-Lbar_mu) >= 0.                                       (4)
```

Thus an exact bipartite flow from positive coarser types to all negative finer
types proves that slot positive.  The type `(1^N)` is deleted on both sides
because its Young subgroup is trivial and `Lbar_(1^N)=0`.

## 3. Stable degree bound

The apparent `4N-9` interpolation bound has a sharper structural form.
Inside either chamber every residual alphabet is a signed union of integer
intervals whose endpoints are linear in the two chamber coordinates.  Its
power sum `p_k` therefore has degree at most `k+1`.  The augmented-monomial
inclusion-exclusion formula, or equivalently the triangular power-sum
recurrence used by the companion, gives

```text
deg m_mu <= N+length(mu).                                      (5)
```

Both products in the numerator of `(2)` consequently have degree at most
`3N+length(mu)`.  The inherited forced divisor has degree nine, so

```text
deg G_mu <= 3N+length(mu)-9.                                  (6)
```

Only `mu=(1^N)` can reach `4N-9`, and its operator is zero.  Every relevant
type obeys the sharper bound

```text
deg G_mu <= 4N-10.                                             (7)
```

This explains, rather than merely observes, the vanished top shell in the
exact runs.  Polynomiality follows from the inherited forced divisor for
homogeneous symmetric functions (equivalently, expand `m_mu` in the Schur
basis before applying `Phi`).

## 4. Exact `N=8` and `N=9` results

The companion reconstructs the full bivariate Newton triangle, checks the
next shell, and checks three remote off-grid points.  At every nonzero slot it
runs an exact rational max flow along `(1)`.

For `N=8`, each bank/chamber has `300` nominal slots:

```text
I1 tight: 276 flow + 24 zero; 13 sign patterns; 37 used edge types
I1 wide : 276 flow + 24 zero; 14 sign patterns; 47 used edge types
I2 tight: 276 flow + 24 zero; 13 sign patterns; 38 used edge types
I2 wide : 276 flow + 24 zero; 14 sign patterns; 46 used edge types.
```

For `N=9`, each bank/chamber has `406` nominal slots:

```text
I1 tight: 378 flow + 28 zero; 16 sign patterns; 66 used edge types
I1 wide : 378 flow + 28 zero; 17 sign patterns; 78 used edge types
I2 tight: 378 flow + 28 zero; 16 sign patterns; 70 used edge types
I2 wide : 378 flow + 28 zero; 18 sign patterns; 77 used edge types.
```

Thus the zero census continues the exact law `4(N-2)`, while every other
slot transports.  In all eight cases the singleton coefficient is negative
or zero, never positive.  This independently rules out the subtle invalid
certificate in which mass from the zero operator might be used as supply.

The restricted-carrier hostile is also stable.  At `N=8`, hooks alone pass
only `230--241` of the `276` live slots and two-row types pass only `219--229`.
At `N=9`, the corresponding ranges are `299--313` and `284--298`.  The union
of hooks and two-row types performs exactly like hooks in every case.  Hence
the full refinement poset is not cosmetic: genuinely higher-width positive
types are needed.

## 5. What is proved by the scout and what is not

Combining `(1)--(7)` with the exact flow files is a finite exact certificate
for positivity of the row-normalized Young-gap operators at `N=8,9`, for both
THM-3110 invariant banks and every integer support `0<a<b`.  The files remain
a reflection/scout until integrated into audited canon.

The source-to-target map is:

```text
source: chamber-Newton coefficient vector on occupancy types
target: positive Young-subgroup gap operator
map:    rational transport from negative fine types to positive coarsenings
keeps:  full block-size type and chamber coefficient
loses:  labelled fibre, roots, insertion path, and cross-degree compatibility
sidecar needed for all N: a symbolic flow or inductive insertion kernel.     (8)
```

No uniform row, hook, two-row, octopus, or same-sign forest reduction follows.
The computation does not prove the all-degree row extremum, arbitrary anchored
product-Gamma goodness, SFC in width three, GMC(2), NC2, LRC(14), JC(2), or
DC(2).  Its strongest forward signal is the stable `4(N-2)` zero law plus
complete refinement transport through `N=9`; the next useful object is an
insertion-compatible flow from degree `N` to `N+1`, not a larger static cone.

## 6. Reproduction

Run

```text
python 04-computation/gmc_monomial_fibre_refinement_n8_scout.py
python -O 04-computation/gmc_monomial_fibre_refinement_n8_scout.py
python 04-computation/gmc_monomial_fibre_refinement_n8_scout.py 9
python -O 04-computation/gmc_monomial_fibre_refinement_n8_scout.py 9
```

and compare the first pair with

```text
05-knowledge/results/gmc_monomial_fibre_refinement_n8_scout.out
```

and the second pair with

```text
05-knowledge/results/gmc_monomial_fibre_refinement_n9_scout.out.
```

The script uses exact integers and rational numbers only.  The certificate
digests in the stored outputs hash every nonzero transport edge and amount.
At this checkpoint the default `N=8` transcript has independently replayed
under both normal and optimized Python.  The `N=9` transcript records one
complete exact normal run; its optimized byte replay was deliberately stopped
to release the shared machine for the canonical `N=5,6,7` audit, so that last
evidence check remains pending even though no assertion failed.
