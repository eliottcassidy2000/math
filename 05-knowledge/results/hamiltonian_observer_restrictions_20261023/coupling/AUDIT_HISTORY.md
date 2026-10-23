# Independent audit of the recovered cycle and gain lift

Status: **PASS**, with the sign-precision note below. This audit inspected
[`VALUATION_CYCLE_GAIN_LIFT.md`](../history/VALUATION_CYCLE_GAIN_LIFT.md), its
checker, and the inherited rotation-system carrier at
[`geometry_567_flag_carriers_20261004.md`](../../geometry_567_flag_carriers_20261004.md).
It does not promote the supplied gain representation to an arithmetic
iteration or a tournament/Collatz theorem.

1. Identifying the endpoints of a five-edge path creates five quotient
   states. The added reciprocal reset is a sixth *new* edge and not a sixth
   original valuation-map iteration. The `N/2` normalization gives period
   two and the same endpoint ratio, as stated.
2. The quotient-cycle consistency equation is exactly `H=1`. Rational
   vertex gauge factors telescope on a cycle and concentrate all debt at
   one seam. The complex constant-gain gauge is legitimate because the
   recursively constructed frame satisfies `h_r=h_0`.
3. A change of logarithm branch cannot remove the nonzero real displacement
   `log(4/13)`. The document correctly distinguishes one full change of
   argument on an open spiral from a closed-curve winding number.
4. The rational valuation/unit homomorphism handles denominators divisible
   by 2 or 3 before reducing unit residues. Directly,
   `rho_(2,3)(4/13)=(2,5)` and `rho_(3,2)(4/13)=(0,1)`.
   The hostile `73=1+72` lies in both register kernels. Neither the full
   valuation-map quotient nor exact gain closure follows from this residue
   agreement.
5. The inherited geometry requires the actual dart rotation permutation.
   The formula `chi=d(1/2+1/3+1/q-1)` and its signs at `q=5,6,7` are correct
   for the specified cubic map. The document correctly refuses to infer
   that rotation system from the scalar gain cycle or the shared graph.
6. The all-permutation mod-9 sum and the Hamiltonian-subset sum have
   different universes. On the transitive six-vertex control, the former
   is `-189 omega`, while the latter is `-1-omega`, exactly as required by
   the hostile.

Minor scope precision: all prime valuations determine a **positive**
rational gain, while an arbitrary nonzero rational also requires its sign.
The exact rational itself includes the sign. The original construction
uses positive labels, so this clarification does not change its proofs.

Reproduction completed independently:

```text
python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/history/check_gain_lift.py
python -O -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/history/check_gain_lift.py
```

Both runs passed **66,786** checks: 19,525 positive label words, 4,096 ordered
fraction pairs tested in two registers, normalization and modular-hostile
controls, and the 720-permutation transitive order-six control. The script
uses explicit errors rather than optimization-sensitive assertions.
