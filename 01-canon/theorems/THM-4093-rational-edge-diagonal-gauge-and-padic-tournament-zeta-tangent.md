---
id: THM-4093
title: "Rational-edge diagonal gauge and p-adic tournament-zeta tangent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Vertex-ratio
  edge weights are diagonal gauge and leave every closed-walk spectral
  invariant unchanged. For a tournament, the first adjacency-zeta term is
  the directed-triangle count at degree three, giving a sharp p-adic tangent.
  Nontransitive tournament zeta sends transcendental inputs to transcendental
  outputs, but raw vertex labels and transitive tournaments remain
  arithmetic-type blind. This is not a p-adic L-value irrationality theorem.
source: codex-padic-zeta-tournament-20260825
audit: >
  PASS WITH SCOPE REPAIRS INCORPORATED. The primary companion checks exact
  similarities and every labelled tournament through order six. The
  independent standard-library verifier checks 6,148 determinant, 7,680
  power/trace, and 102,700 valuation gates, including zero-weight and
  p=2,c3=2 strict-depth hostiles. Normal and optimized outputs agree.
depends_on:
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
related:
  - THM-1926-tournament-zeta-euler-product-over-strong-core
  - THM-4059-stern-brocot-depth-packet-character-and-divisor-star-convolution
  - MISTAKE-409
script: 04-computation/rational_edge_diagonal_gauge_padic_tournament_zeta_thm4093.py
output: 05-knowledge/results/rational_edge_diagonal_gauge_padic_tournament_zeta_thm4093.out
independent_audit_script: .scratch/gauge_p13_referee_20260825/independent_audit.py
independent_audit_report: .scratch/gauge_p13_referee_20260825/REPORT.md
script_sha256: ef82e05ed85f30fef95c79c24c8b972d36cfaea990b930c381220b9c66153c0f
output_sha256: f3deed30e5cb68a869d80e88596a2268f5dabc894eaa6a284ad9bf2b1071c0ac
independent_script_sha256: a5b0861f1056d880dc8b7af341cf3b4493a8370ed67eb6b628c9cf6e73b9ba6e
hash_basis: raw LF bytes
---

# THM-4093 -- rational-edge gauge and the tournament-zeta tangent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Vertex-ratio weights are exact diagonal gauge

Let `A=(A_ij)` be an `n` by `n` matrix over a commutative ring `R`, and let
`d_1,...,d_n` be units of `R`.  Put

```text
D=diag(d_1,...,d_n),              W_ij=A_ij d_i/d_j.    (1)
```

Then

```text
W=D A D^(-1).                                             (2)
```

Consequently, whenever the displayed expressions are defined,

```text
det(lambda I-W)=det(lambda I-A),
tr(W^k)=tr(A^k),
det(I-uW)=det(I-uA).                                     (3)
```

For a directed walk `i_0,...,i_k`, its weighted product is

```text
product_(r<k) W_(i_r,i_(r+1))
  = (d_(i_0)/d_(i_k)) product_(r<k) A_(i_r,i_(r+1)).     (4)
```

Thus all vertex ratios cancel on a closed walk.  Over a valued field,
writing `h_i=v(d_i)` gives the exact coboundary

```text
v(W_ij)=v(A_ij)+h_i-h_j.                                (5)
```

Equations (2)--(5) prove the claim.  Nonzero is sufficient over a field; over
a general ring each `d_i` must be a unit.

This is the precise boundary for the proposed rational-edge decoration.  A
weight “numerator label divided by denominator label” can retain endpoint
data on an open path, but it is spectrally pure gauge.  Flipping an adjacency
orientation is a different operation and is not covered by (1).

## 2. The directed-triangle tangent

Let `T` be a finite tournament with zero-diagonal adjacency matrix `A`, and
let `c_3(T)` be its number of cyclic three-vertex subtournaments.  Define the
adjacency determinant and its formal reciprocal by

```text
P_T(u)=det(I-uA),                   Z_T(u)=1/P_T(u).     (6)
```

Then for an integer polynomial `Q_T`,

```text
P_T(u)=1-c_3(T)u^3+u^4 Q_T(u).                          (7)
```

Indeed, the linear coefficient vanishes because `A` has no loops.  The
quadratic coefficient vanishes because a tournament has no directed
two-cycle.  A principal `3` by `3` minor is `1` precisely on a cyclic triple
and `0` on a transitive triple; the determinant sign in (6) supplies the
minus sign.  All remaining terms have degree at least four.

Fix a prime `p` and nonzero `x in p Z_p`, with `v_p(x)=m>=1`.  Equation (7)
shows that `P_T(x)` is a `p`-adic unit and

```text
Z_T(x)-1
 = x^3 (c_3(T)-x Q_T(x)) / P_T(x).                     (8)
```

Therefore

```text
v_p(Z_T(x)-1) >= 3m,                                   (9)
v_p(Z_T(x)-1)  = 3m       iff       p does not divide c_3(T). (10)
```

If `p` divides `c_3(T)`, the valuation is strictly greater than `3m`, with
the value `+infinity` allowed when the numerator vanishes.  This proves the
necessity as well as the sufficiency in (10).  The exact hostile

```text
P_T(u)=1-2u^3-u^4,             p=x=2                  (11)
```

has `v_2(Z_T(2)-1)=5>3`, so the unit condition cannot be dropped.

## 3. What the zeta can and cannot see about arithmetic type

For rational `alpha` away from a pole, `Z_T(alpha)` is rational.  For
algebraic `alpha` away from a pole it is algebraic.  If `T` is nontransitive,
then it contains a directed triangle, so `c_3(T)>0` and `P_T` is nonconstant.
For every transcendental `alpha`, `P_T(alpha)` is nonzero and

```text
Z_T(alpha) is transcendental.                            (12)
```

For if `Z_T(alpha)=beta` were algebraic, then `beta` would be nonzero and
`alpha` would satisfy the nonzero polynomial
`P_T(X)-beta^(-1)` over the algebraic numbers, contradicting
transcendence.

This is a functional preservation statement, not a classification from the
bare tournament.  A transitive tournament has nilpotent adjacency matrix,
`P_T=1`, and constant `Z_T=1`, erasing every input type.  Even for a
nontransitive tournament, algebraic output alone need not distinguish a
rational input from an algebraic irrational input.  Arithmetic labels that
do not enter the intrinsic relation or a non-gauge decoration remain
invisible under relabelling.

## 4. Stern-tournament controls and quotient loss

Apply (6) to the finite Stern-depth tournaments of THM-4057.  The exact
companions give

```text
T_5:   apex imbalance B=0,    c_3=1,   P=1-u^3,
T_9:   apex imbalance B=0,    c_3=17,
T_13:  apex imbalance B=0,    c_3=52.                  (13)
```

Thus the one-vertex apex imbalance can vanish while the closed-walk zeta is
highly nontrivial.  This is a hostile example to recovering cycle data from
the THM-4057/4059 lower-star statistic.  A lawful arithmetic refinement needs
at least non-exact cycle holonomy, character-resolved rows, or a retained
height/error sidecar; a vertex potential alone cannot supply it.

The source, target, and loss of the rational-edge connection are now exact:

```text
ordered coprime arcs -> vertex-ratio weighted digraph -> closed-walk zeta
preserved: adjacency cycle products and spectrum
destroyed: endpoint scale, open-path ratio, arithmetic label type.        (14)
```

## 5. Verification and hostile audit

The primary companion checks `126` exact rational diagonal similarities,
including characteristic polynomials and powers.  It exhausts all `33,867`
labelled tournaments through six vertices and all `665,864` of their
three-vertex slots for (7), then checks `190` exact `p`-adic tangent rows,
including unit multiples of `p^m`, the three Stern controls, (11), and the
`B(T_9)=0` hostile.

The independent standard-library route does not import SymPy or the primary
code.  It exhausts all `512` looped/digoned three-vertex digraphs against
three rational gauge rows, all tournaments through six vertices, `102,700`
valuation gates, the sharp hostile (11), and the exact `p=13` shell sidecar
used only in the session reflection.  Both normal and optimized routes pass.

## 6. Scope

`Z_T` in this theorem is an adjacency-determinant zeta, not a Kubota--Leopoldt
`p`-adic zeta or `L`-value.  The p-adic tangent (9)--(10) proves no
irrationality or transcendence of special zeta values.  It identifies the
first non-gauge tournament coordinate—directed triangles—and specifies the
extra information any such arithmetic transfer would have to retain.
