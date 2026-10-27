# Named routes, integral carriers, and exact energy invoices

Session 75 final mathematical notes, 27 October 2026. These notes explain the
results and retain a small selection of standalone computations. They are not
an archive of execution records or a claim that the included examples reproduce
all earlier censuses. No claim of historical priority is made.

## What is proved, and what is finite evidence?

* [Black routes](BLACK_ROUTES.md): every prescribed vertex of a regular
  tournament of order `2m+1 >= 7` can be separated from a strong complementary
  core in exactly `m` black native moves, with either native source/sink role.
  A more general minimum-semidegree region and a separate explicit two-move
  family are proved with their precise dependencies.
* [Integral Paley carriers](INTEGRAL_CARRIERS.md): for field-Paley order `p^3`,
  the minimal maximal-order-stable enlargement has additive index
  `p^(((p+1)/2)^3-1)`. Stability is not metric integrality or integral similarity.
* [Clocked occupation](CLOCKED_ROOT.md): a finite clocked killed graph has a
  least admissible occupation weight and an exact norm-defect identity. This
  is not a contraction theorem in the original coherent-source norm.
* [Pairing and guard information](PAIRING_GUARDS.md): three four-child pairings
  have an exact average energy, integral reconstruction has a parity condition,
  and linear guards on the full bit cube have an exact extension cost.
* [Eisenstein invoices](EISENSTEIN_INVOICES.md): an integer-slack indicator lift
  yields global, leaf-corrected merge potentials. The special catalyst optimum
  at total 11 is separately identified as an audited finite certificate result.

The included programs check an explicit family route, finite-field and
clocked-weight identities, and elementary Eisenstein/pairing controls. They
use explicit exceptions rather than optimization-disabled assertions, do not
read another research directory, and do not write result archives. Running
with ordinary Python and `python -O` gives the same mathematical checks.

General black-core distance, all-scale native ROOT debt, universal integer
arrival, and a complete unrestricted tournament atlas are **not** established.

## Small presentation fibres: a credited, conditional classification

Write `tc(T)=H(T)/|Aut(T)|`, with ordinary directed isomorphism, not converse
identification. The independently checked extension of the repository's
[THM-4812](https://github.com/eliottcassidy2000/math/blob/3858d2afca72fb97a46a3601a6ae4b1e6e99a812/01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md)
has

`S_23(x) = 4x^6 + 6x^7 + 4x^8`.

Thus there are 14 strong ordinary classes with `tc=23`, and 45 strong classes
in the union of odd fibres through 23. This all-order conclusion uses named
prior orbit-product, vertex-transitive and Busch lower-bound inputs, plus a
finite orbital window. It is not inferred from a census through order eight
alone. The repaired vertex-transitive bound is strictly `tc>23` outside the
atoms `K1,C3,R5,P7`; a weak `>=23` bound would not suffice.

The key retained coordinate is the group index. For the actual automorphism
orbits, put `P=product Aut(T[O])` and `I=[P:Aut(T)]`. The inherited orbit-product
inequality implies

`tc(T) >= I product_O tc(T[O])`.

The index cannot be dropped when deciding which orbit configurations survive.
The infinite tail is closed by comparing Busch's strong lower bound `f(n)`
with the maximum compatible product-group order `B(n)`: beyond the finite
base range, `f(n+3)=5f(n)` but `B(n+3)=3B(n)`. The finite comparison and its
strict inequalities were checked separately; the three programs supplied here
are **not** a replacement for that census or those group inputs.

For all tournaments, the ordered strongly-connected-component word gives

`A_23(x) = S_23(x)/(1-x-x^3)^2`.

Indeed normalized path counts multiply along the uniquely ordered SCCs;
23 is prime, so exactly one component has value 23 and every other component
has value one. The value-one strong alphabet is `K1,C3`, giving an arbitrary
word on either side. This is a counting bijection, not a natural operation
that changes an arbitrary tournament's `tc` value. The coefficients through
orders 6–12 are `4,14,28,50,92,162,272`; the later values are consequences of
the series, not fresh full censuses at those orders.

## Earlier public stages

The earlier [checkpoint 3 research note](https://github.com/eliottcassidy2000/math/blob/318d10ff81b7e9f147e22162d32762059ef15931/05-knowledge/results/session75_overnight_black_root_spectral_20261027/README.checkpoint3.md)
contains the prior black-cut counterfamilies, native cutoff calculations and
normalized `tc` exclusions. The present notes preserve the distinctions
between raw `H`, normalized `tc`, path-reversal and native-complement gauges,
source-labelled and coherent norms, and node-positive and level-positive debt.
