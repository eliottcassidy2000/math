---
id: HYP-8802
title: "Direct A2 quantization of the rank-two Poisson witness"
status: >
  OPEN / first obstruction and one-pair repair computed. Weyl-symmetric
  quantization of THM-2044 has a nonzero cubic Moyal anomaly in five of the six
  canonical relations. For the (D,R) relation the anomaly lies in the
  R-centralizer and a finite two-step polynomial correction restores the exact
  Weyl relation. Simultaneous correction of all six relations, with termination
  rather than an infinite formal hbar series, is the DC(2) gate.
source: codex-2026-07-21-DC2-JC2
related:
  - THM-2044
  - THM-2045
  - THM-1345
script: 04-computation/poisson_rank2_symplectic_suspension_codex_20260721.py
output: 05-knowledge/results/poisson_rank2_symplectic_suspension_codex_20260721.out
---

# HYP-8802 -- direct quantum descent

Let `P^n(f,g)` denote the `n`-fold bidifferential power of the canonical
Poisson tensor. In Weyl symbols, the normalized Moyal commutator at `hbar=1`
is

```text
M(f,g)=P(f,g)+P^3(f,g)/24+P^5(f,g)/1920+... .       (1)
```

For the four THM-2044 symbols, `P^3` is nonzero for

```text
(D,R), (S,T), (R,S), (D,T), (D,S),
```

and vanishes only for `(R,T)`. Their respective nonzero term counts are
`42,42,3,165,273`. Thus naive symmetric ordering is not an `A_2`
endomorphism.

There is nevertheless a finite exact repair of the first Weyl relation. Put

```text
A=P^3(D,R),
C1=-A D0/24,
B=P^3(C1,R),
C2=-B D0/24,
Dq=D+C1+C2.                                         (2)
```

Here `{D0,R}=1`, `{A,R}=0`, and `R` is cubic, so (1) terminates at `P^3` when
its second argument is `R`. Exact calculation gives

```text
B=108x^12(3xq-2)(3xq-1)(27x^2q^2-2)                (3)
```

and

```text
M(Dq,R)=1.                                          (4)
```

The corrected symbol has 332 terms. Equation (4) is meaningful progress but
not DC(2): correcting `D` perturbs its three other relations, and corrections
to `T,S,R` may generate higher odd Moyal orders. The decisive target is a
simultaneous polynomial solution whose `hbar` recursion terminates. A merely
formal power-series quantization is insufficient at `hbar=1`.

For planar JC, THM-2045 closes the most direct de-stabilization route. A useful
positive program is to turn its Laurent-sector obstruction into a criterion on
the inner-polynomial/Newton-polygon candidates of Lee--Li: stabilization can
repair a unimodular gradient row, but the planar mate requires the unique
constant-producing diagonal sector to satisfy an integrability condition.
