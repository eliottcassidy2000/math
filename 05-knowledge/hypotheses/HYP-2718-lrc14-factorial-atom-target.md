---
id: HYP-2718
title: LRC14 carrier-product gap is a factorial-origin atom target
status: OPEN; exact basis theorem proved by THM-561
source: codex-2026-06-21
depends_on:
  - THM-561
  - HYP-2717
  - HYP-2716
  - HYP-2715
  - HYP-2714
related:
  - THM-406
  - THM-534
  - HYP-2680
  - HYP-2698
  - THM-558
  - OPEN-Q-108
---

# HYP-2718 - Factorial-Origin Atom Target

## Claim

The remaining HYP-2714 multi-block carrier-product lemma should be proved in
the factorial atom basis.  Let `T` be the number of missed inner sectors in
`{1,...,6}`.  For a split row compare the actual line law with the
shared-slow-`x` independent-carrier product law, and write

```text
W_j = Delta E[binom(T,j)]
    = sum_{|R|=j} (z_prod(R)-z_actual(R)),

Q_t = Delta Pr(T=t).
```

By THM-561,

```text
Q_t = sum_{j=t}^6 (-1)^(j-t) binom(j,t) W_j.
```

Therefore the HYP-2716 top character is not only a Boolean-cube shadow; it is
the `T=0` atom after exact factorial-moment inversion:

```text
ProductCover - p0 = Q_0
                  = sum_j (-1)^j W_j
                  = M_6.
```

The next proof obligation can be stated as the atom budget

```text
|Q_0(E)| <= cap_k - ProductCover(E)
```

for every HYP-2714 moderate-span balanced split row.

## Why This Is Sharper

The raw miss-zeta coordinates and the raw missed-count law both contain larger
motions than the cover error.  The factorial transform isolates the boundary
atom that matters.

In the exact scout, the aggregate missed-count atom tournament is

```text
t=2 > t=1 > t=0 > t=4 > t=3 > t=5 > t=6.
```

So `Q_0` is not the largest atom discrepancy.  Bounding all atoms, or all
zeta coordinates, asks for more than the LRC14 cap needs.  The cover predicate
is the alternating finite difference at the origin, exactly like the old
finite-difference handoff where falling-factorial divisibility exposed hidden
roots.

This also clarifies how to read rising/falling factorials here:

```text
binom(T,j) = (T)_j/j!,
(T)_j = (-1)^j(-T)^{overline j}.
```

Thus the LRC14 cover gap is a signed boundary evaluation in the falling basis,
or equivalently a rising basis in the reflected variable `-T`.

## Exact Evidence

The updated scout

```text
04-computation/lrc14_multiblock_miss_zeta_layers_codex_20260621.py
05-knowledge/results/lrc14_multiblock_miss_zeta_layers_codex_20260621.out
```

now prints `Q_t=Pr_product(T=t)-Pr_actual(T=t)` for each tested split row and
asserts `Q_0=ProductCover-p0`.

Representative atom profiles:

```text
two 4-blocks, moderate gap:
  Q_0=-3495299/477209040
  Q_1=-1520515/63627872
  Q_2=33516137/954418080

3+3+2 split:
  Q_0=-51629953/57697542630
  Q_1=-56452381/15386011368
  Q_2=4048085941/461580341040

five 2-blocks:
  Q_0=-2447628624709/93106921650624
  Q_1=-43399251952823/1047452868569520
  Q_2=102517895003/1518047635608
```

Across the six-row bank,

```text
sum |Q_2| ~= 0.140231675,
sum |Q_1| ~= 0.095127546,
sum |Q_0| ~= 0.058978420.
```

The cover atom is smaller than the dominant `t=1,2` atom motions but still
larger than the high-tail atom motions.  That profile suggests a proof using
an origin-atom finite-difference estimate plus a finite resonance ledger, not
a stochastic-dominance statement for the whole missed-count distribution.

Incoming S68 evidence is compatible with this target.  The exact
moderate-span scout

```text
04-computation/lrc14_moderate_multiblock_budget_codex_s68.py
05-knowledge/results/lrc14_moderate_multiblock_budget_codex_s68.out
```

extends finite windows to `k=8,9,10` and finds no over-cap or near-cap rows.
The newly added moderate-balanced leaders have large margins:

```text
k=8:  margin=1/6
k=9:  margin=17499/140140
k=10: margin=30437/194922
```

Its carrier-product diagnostic has both signs (`product-exact` ranges from
negative examples such as `-27973/352800` to positive examples such as
`2423/44100` in the top-row bank), reinforcing that product is not an
envelope.  The proof needs a signed scalar error bound for the cover atom.

## Proof Route

Combine THM-561 with HYP-2717 as follows.

1. Expand the top/origin atom `Q_0=M_6` in carrier Fourier modes.
2. Split the modes by the carrier relation lattice `n.M=0` versus `n.M!=0`.
3. Use factorial-basis cancellation before absolute values are taken.
4. Route low-height exact relations and small-denominator nonrelations to the
   HYP-2714 finite ledger.
5. Prove the high-height/high-denominator tail fits under `cap_k-Product`.

The old THM-406 and THM-534 sector-moment identities say the same thing in
different coordinates: cover depth is controlled by all factorial moments, but
the useful certificate is a low-dimensional polynomial in `binom(T,j)`.  Here
the target polynomial is the delta atom `1[T=0]`, evaluated on the difference
between actual and product laws.

## Assumption Challenge and Tournament Analysis

This session considered vertices as runners, carrier blocks, residual masks,
residual-size layers, Boolean characters, carrier Fourier modes, missed-count
atoms, factorial moments, and proof obligations.  The productive vertex set
for HYP-2718 is the seven atom discrepancies `Q_0,...,Q_6`.

Pairwise observable: larger aggregate `sum |Q_t|` over the tested split rows.
Switch/gauge: invert miss-zeta radial weights from factorial moments to
missed-count atoms.  Tie Hamiltonian path:

```text
t=2 > t=1 > t=0 > t=4 > t=3 > t=5 > t=6.
```

Fingerprint: transitive tournament, score histogram
`{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, zero directed 3-cycles.

The quotient preserves the exact cover correction at `t=0` and destroys
sector labels and which block created each miss.  Those labels must re-enter
only through HYP-2717's finite low-relation ledger.

Challenged assumption: the carrier-product error should be dominated as a
positive cone or as a full missed-count law.  The factorial-origin atom target
requires only one boundary finite difference.

## Status

THM-561 proves the basis identity.  HYP-2718 is still open as an LRC14 analytic
bound: the missing work is a relation-filtered signed estimate for `Q_0`.
