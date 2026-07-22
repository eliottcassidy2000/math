---
id: HYP-8802
title: "Direct A2 quantization of the rank-two Poisson witness"
status: >
  OPEN / exact R-column repair and formal S,T boundary lift computed. Weyl-symmetric
  quantization of THM-2044 has a nonzero cubic Moyal anomaly in five of the six
  canonical relations. Finite polynomial corrections restore the entire
  R-column: M(Dq,R)=1 and M(T,R)=M(Sq,R)=0. THM-2049 proves that the Ore
  boundary complex for the S,T relation is associated-graded acyclic and gives
  a formal beta-adic lift. Polynomial termination and coupling to D remain.
  Simultaneous finite correction of all six relations is the DC(2) gate.
source: codex-2026-07-21-DC2-JC2
related:
  - THM-2044
  - THM-2045
  - THM-2046
  - THM-2049
  - THM-1345
  - HYP-8803
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

There is nevertheless a finite exact repair of the whole `R`-column. Put

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

The corrected symbol has 332 terms. For `S`, let

```text
A_S=P^3(S,R),             Sq=S-A_S D0/24.           (5)
```

Here `{A_S,R}=0`, while `A_S D0` has momentum degree below three, so
`P^3(A_S D0,R)=0`. Therefore

```text
M(Sq,R)=0.                                          (6)
```

The corrected `Sq` has 85 terms. Since `R` is cubic, `P^1` and `P^3` are the
only possible odd layers against it; moreover `P^3(T,R)=0`. Thus (4), (6), and
`M(T,R)=0` are exact, not truncated checks.

This is meaningful progress but not DC(2): the remaining residual
`M(Sq,T)-1` has 59 terms and momentum degree two (and, as Jacobi predicts,
star-commutes with `R`). Correcting it while preserving the exact `R`-column
must then be coupled to `M(Dq,T)=M(Dq,Sq)=0`. The decisive target is a
simultaneous polynomial solution whose recursion terminates. A merely formal
power-series quantization is insufficient at `hbar=1`.

## The first T-column descent and the six-weight cascade

The separate `--quantum-t` audit rewrites the residual in the adapted
centralizer coordinates `(x,q,D0,ell)`. It is independent of `D0`, has 31 terms,
and has `ell`-degree two. Its leading coefficient is

```text
-(3/2)x^10(3xq-2)(18(xq)^3-39xq+19).                (7)
```

Let the star denote the Weyl/Moyal product. The star-central correction

```text
C_ST,1=f star ell,
f=-(x^6/16)(108(xq)^3-72(xq)^2-93xq+58)            (8)
```

kills (7) exactly while preserving the entire `R`-column. The remaining
`ell`-coefficient splits into weights four and seven; both one-variable
equations

```text
(2-3s)K'(s)+3mK(s)=H_m(s),          m=4,7,          (9)
```

have polynomial solutions. The resulting multiplication correction kills the
whole `ell`-layer. What remains is the nonzero constant

```text
-(3/16)x^6(26163q^2x^8-31455qx^7+330qx+9402x^6-188). (10)
```

The homogeneous freedom in (8) is `f_m=xR^(m-1)`. After normalizing its
`ell`-layer by (9), exact audits give the constant-weight propagation

```text
m=6:  {6,12},       m=12: {12,18},       m=18: {18,24}. (11)
```

The modes `6,12,18` do not solve (10) simultaneously. This rules out only that
fixed-`T`, tangent-linear repair. THM-2049 supplies the required hostile
control: after `T` is allowed to move, the simultaneous associated-graded map
is surjective in every grade. Hence the `m -> m+6` cascade is a gauge artifact
of the fixed-`T` slice, not a general obstruction. The live question is
whether the resulting exact correction ladder terminates as a polynomial and
can be coupled to the `D` relations.

For planar JC, THM-2045 closes the most direct de-stabilization route. A useful
positive program is to turn its Laurent-sector obstruction into a criterion on
the inner-polynomial/Newton-polygon candidates of Lee--Li: stabilization can
repair a unimodular gradient row, but the planar mate requires the unique
constant-producing diagonal sector to satisfy an integrability condition.

## Ore localization update

HYP-8803 proves that the `x,q,ell` part of this calculation is the exact Ore
extension `Q[x,q][ell;delta]`, with `R` central, and that after inverting `x`
the slice `t=-1/(3x)` satisfies `[ell,t]=1`. It also proves that Weyl ordering
is the unique scalar PBW ordering that lowers the raw `[S,T]-1` residual from
`ell`-degree three to degree two. THM-2049 then proves the associated-graded
simultaneous correction complex is acyclic and explicitly advances the exact
residual from beta grades 6 through 14. Thus the cascade above should be
attacked as a finite-termination and `x`-adic boundary-integrality problem,
rather than by searching for another scalar ordering, local Darboux
coordinate, or first cohomology obstruction.
