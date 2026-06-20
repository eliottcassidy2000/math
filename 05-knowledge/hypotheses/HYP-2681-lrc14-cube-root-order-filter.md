---
id: HYP-2681
title: LRC(14) cube-root order filter for signed multi-far residuals
status: OPEN; exact algebraic lens, proof value under test
source: codex-2026-06-20-S52
depends_on:
  - HYP-2680
  - THM-548
  - HYP-2679
  - HYP-2677
  - HYP-2676
  - HYP-2639
related:
  - THM-511
  - HYP-2532
  - HYP-2533
  - OPEN-Q-108
---

# HYP-2681 - Cube-Root Order Filter

## Claim Being Tested

The simultaneous-peel residuals from HYP-2680 may carry a useful `C_3`
structure.  For a bounded core `B` and far triple `F={u,v,w}`, write the seven
nonempty residual packets

```text
A = Delta_u  - Phi_1
B = Delta_v  - Phi_1
C = Delta_w  - Phi_1
D = Delta_uv - Phi_2
E = Delta_vw - Phi_2
F = Delta_wu - Phi_2
G = Delta_uvw - Phi_3.
```

The actual simultaneous-peel correction is

```text
H(1) = A+B+C+D+E+F+G.
```

The user's recursion

```text
A+B+C-D-E-F+G
```

is the pair-tax shadow

```text
H(1) - 2*(D+E+F).
```

The cube-root filter keeps more information.  With `omega^2+omega+1=0`, define
Eisenstein cyclic modes

```text
S_omega = A + omega*B + omega^2*C
P_omega = D + omega*E + omega^2*F.
```

The first question is whether dangerous triples are controlled not by scalar
rank alone, but by the scalar pair-tax shadow together with the exact
Eisenstein imbalance `S_omega-P_omega`.

## Challenged Assumption

Do not assume the seven packets should be collapsed immediately to total
residual or absolute mass.  The quotient

```text
(A,B,C ; D,E,F ; G) -> (H(1), A+B+C-D-E-F+G, S_omega-P_omega)
```

preserves the signed order ledger, the pair layer, and the cyclic orientation
imbalance.  It destroys individual packet ownership beyond the chosen cyclic
order, so any use in proof must record the core/far phase address.

## Next Computation

Add an exact Fraction scout that:

1. computes the seven packets for named triples from HYP-2680;
2. computes the pair-tax recursion `A+B+C-D-E-F+G`;
3. computes exact Eisenstein coefficients and norms for `S_omega`, `P_omega`,
   and `S_omega-P_omega`;
4. scans the `(15,16,17)` all-core bank to see whether direct cap risk,
   three-far residual risk, and cube-root imbalance are different orderings.

No proof of LRC(14) is claimed.
