# EI17 over finite fields: the last small four-chromatic triangle-free unit graph against Q(zeta55)

2026-10-17, opus-2026-10-17-S39. **PROVED (reduction lemma, exhaustive elimination) + FINITE-EXACT (counts).** Canon:
[THM-4665](../../01-canon/theorems/THM-4665-ei17-has-no-faithful-realization-over-the-split-and-conjugation-fixed-primes-of-q-zeta55-tested-so-every-q-zeta55-realization-has-bad-reduction-there.md).
Not an exclusion proof; HYP-9247 is updated, not settled.

## Why

[HYP-9247](../hypotheses/HYP-9247-every-finite-unit-distance-graph-in-q-zeta55-is-three-colourable.md) (every finite
unit graph in `L = Q(zeta55)` is three-colourable) is decisive for the six-colour programme (THM-4614, THM-4661). Of the
four smallest triangle-free four-chromatic unit graphs, EI19 (S24), the fish and EI21 (THM-4662) are excluded from `L`;
EI17 was excluded only numerically, to coordinate height `1e63`
([S24 note](ei17_placement_q_zeta55_20261012.md)), and the exact elimination of its six-angle system did not finish.

## The device

This is the repository's finite-field counting device (the `H(P_31)` and Kloosterman-moment computations), applied
to a rigid graph. In the gauge of S24, the realizations are six free unit directions with six norm-one conditions on
integer combinations of them. Over `F_p` the system is triangular: two directions free, four determined by quadratics
(the codegree lemma of THM-4663), two checks. That is `≈ 16 p^2` operations per prime, so primes up to `10^4` take
seconds to minutes. Two models:

* **split** (`F_p`, conjugate `x ↦ 1/x`): where realizations in `L` land at primes `p ≡ 1 (mod 55)`;
* **twisted** (`μ_(p+1) ⊂ F_(p^2)`): embeddings of EI17 in the finite unit-distance graph `G_p`; where realizations in
  `L` land at the conjugation-fixed primes `p ≡ 54 (mod 55)`.

## Results

| primes | faithful points |
|---|---|
| `p ≡ 1 (mod 55)`, below 12100 (32 primes) | none at 26; 8 or 16 at `2861, 3191, 4621, 4951, 7481, 10781` |
| `p ≡ 54 (mod 55)`, below 6000 (20 primes) | none at all 20 (EI17 is not a subgraph of `G_p`) |
| controls: all primes below 700 (split) / 600 (twisted) | split: 11 primes with points (`103, 139, 281, ...`); twisted: 14 (`127, 149, 227, ...`) |

Every faithful solution passes the 31-edge check; degenerate branches identify two vertices (proved) and are enumerated
or skipped accordingly. **Consequence (THM-4665 (3)):** a realization of EI17 with all directions in `L` has bad
reduction (a direction that is not a unit, or two vertices congruent modulo a prime and its conjugate) at every prime
of `L` above each of these 46 rational primes.

## Reading, typed

* The finite-field route answers the question S24 posed ("the elimination polynomial mod good primes gives the
  Frobenius data") without the elimination. It does not by itself give an exclusion: a realization of very large
  height could be degenerate at all these primes. Combined with the LLL bound, a realization in `L` would have to have
  height above `1e63` and degenerate reduction at the primes above 46 specific primes.
* EI17 does embed in many finite unit-distance graphs (`G_127, G_149, ...`): the obstruction is specific to the classes
  `±1 (mod 55)`. (A further pattern, that the split model never had points at the 38 primes `p ≡ 10 (mod 11)` below 2600,
  has Poisson probability of a few percent against the base rate and no Galois explanation consistent with the twisted
  data; it is recorded as an observation only.)
* The counts come in multiples of 4 or 8, suggesting a free symmetry group of order 8 on faithful realizations (mirror,
  the graph automorphism, a third involution not identified).

## Next

An unconditional exclusion would need a bound on the bad primes of the realization scheme over `Z`, for example from
the discriminant of an elimination polynomial computed modulo many primes and reconstructed. The residue-degree-ten
primes of `L` remain untested.

## Reproduction

```
python -B 04-computation/experiments/ei17_finite_field_realizations_20261017.py 700
python -B 04-computation/experiments/ei17_finite_field_realizations_20261017.py 0 331,661,...   (primes 1 mod 55)
python -B 04-computation/experiments/ei17_twisted_unit_graph_embeddings_20261017.py 109,439,...  (primes 54 mod 55)
```
