# A joint Frobenius chart and its remaining name-transport obligation

**2026-10-29. FINITE-EXACT construction and obstruction for one specified
lift; PROVED consequences of certified matrix identities. Same-agent
audit. This does not exclude a different compatible lift.**

The source is the eighteen-dimensional ternary module in
[TERMINAL_SELECTOR S5](TERMINAL_SELECTOR.md), with its literal 35-state
boundary transfer and six terminal assignments. Incoming
[PERIODIC_FROBENIUS_LIFT F1--F2](../overnight_research_20261029/PERIODIC_FROBENIUS_LIFT.md)
supplies X->X^487 on the old cyclic module. Its exponent and CRT proof are
inherited. Here the new test is compatibility with the enlarged name action.

Put h=(X+1)g, with g=X^5+X^4-X^3-X^2-1 over F3. Let b0=b_012 and b1=b_021.
The vector

    v=b1-(2+B^2+B^3+2B^4+2B^5)b0

satisfies h(B)v=0. The columns

    K=(b0,Bb0,...,B^11 b0, v,Bv,...,B^5 v)

have rank 18 and give an explicit chart for F3[X]/h^2 + F3[X]/h.
Let T be multiplication by X in this chart, and let S substitute X^487
in each summand. Exact arithmetic verifies

    S^5=I,    S!=I,    ST=T^487 S.

This extends the compatible lift to the entire joint clock module. For
each native name permutation pi, let M_pi be its matrix in K. Although
T M_pi=M_pi T, the ranks over F3 of S M_pi-M_pi S are:

| pi | 012 | 021 | 102 | 120 | 201 | 210 |
|---|---:|---:|---:|---:|---:|---:|
| Rank | 0 | 5 | 10 | 10 | 10 | 10 |

For 021, entry (0,0) of the commutator is already one. These are exact
modular ranks. The receipt stores K,T,S, all six native name matrices,
selector coordinates and witnesses.

The physical source functional lambda=aK satisfies lambda S!=lambda.
The identity selector is fixed by S in this observable-vector chart.
This places the transport cost dually to the incoming chart, which fixes
its initial clock vector and moves its scalar observation. The exact
comparison transports all ingredients:

    z -> Sz,    lambda -> lambda S^-1,
    T -> T^487,  M_pi -> S M_pi S^-1.

Every evaluated response is then preserved. Keeping the old name matrices
while changing the clock fails the displayed commutators. The conjugated
matrices are not established as native re-markings or vertex operations.
The construction therefore does not supply a free field-phase compression.

A different splitting or lift might commute with the native M_pi. That
existence question, and which source or selector data it must move, remain
open. Failure of this specific S is not a universal obstruction theorem.

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/joint_phase.py
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/joint_phase.py
```

Both checks rebuild the common chart from native receipts and reduce exact
integer arithmetic modulo three after every product. No graph realization
or floating-point rank is used.
