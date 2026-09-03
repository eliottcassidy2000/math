---
title: "JC2 row-ten elliptic sign quotient"
status: >
  INTEGRATED AS THM-4385: PROVED GEOMETRIC COROLLARY OF THM-4380 and
  VERIFIED-EXACT + INDEPENDENTLY AUDITED by a 21-check companion and a
  separate proof/implementation review. The
  fixed weight-at-most-twelve carrier is an elliptic curve minus its
  two-torsion, but higher-weight stability, chart or seam entry, JC(2), and
  DC(2) remain OPEN.
source: root / JC2 continuation session, 2026-09-03
artifacts:
  - 01-canon/theorems/THM-4385-source-normal-row-ten-elliptic-sign-quotient.md
  - 04-computation/jc2_source_normal_row10_elliptic_sign_quotient_thm4385.py
  - 05-knowledge/results/jc2_source_normal_row10_elliptic_sign_quotient_thm4385.out
script_sha256: 780481b1610e281846450edf3b2cbc38a2750c0d7cd36a213b7e90908ed2f4df
output_sha256: 2e68ebaa0c15e44347078adc166a328eb96e64b1664dadf5980dc3df1595aba8
hash_basis: raw LF bytes
related:
  - THM-4380-source-normal-weight-twelve-row-twelve-extinction
  - THM-4377-reciprocal-source-mutation-and-boundary-jet-cokernel
  - THM-4378-bilateral-packet-quotient-reciprocal-eigenlattice
---

# The stratum split was an elliptic curve cut at the wrong coordinate

## Signal

The independent THM-4380 audit caught a false adjective: after setting
`r=eta/Phi`, the row-ten curve is not rational. The correction became useful
when the discarded `Phi=0` stratum was compared with the closure of the
`Phi!=0` equation. The exact identity

```text
D(0,eta)=19 eta F(eta)
```

shows that the two genuine boundary points are already points of the same
smooth cubic. Only the rational origin is extraneous. Projective completion
then reveals the clean intrinsic object:

```text
row-ten source projection = elliptic curve E minus E[2].
```

Simultaneous sign change is elliptic negation. Its quotient is the ratio line,
and the squarefree quartic `A(r)B(r)` is the complete branch sidecar. The
septimic `K` avoids this branch divisor, so seven quotient addresses lift to
fourteen reduced source points for a geometric reason rather than a raw count.

## Inheritance pass and concept board

- **Closest proved mechanism:** THM-4380's exact split into `Phi=0` and
  `Phi!=0`, including the two boundary points and fourteen sign lifts.
- **Canonical hostile:** quotienting by `Phi^2` keeps extinction but loses the
  source-point count; THM-4378 independently warns that a signed two-sheet
  quotient can retain a nontrivial gluing class.
- **Corrected near miss:** “rational carrier” failed because
  `y^2=A(r)B(r)` has a squarefree quartic right side.
- **Least-used sidecar:** the omitted origin and the three points at projective
  infinity, which turn out to be exactly the ramification divisor `E[2]`.

| live object | operation | invariant | loss | decisive test |
|---|---|---|---|---|
| row-ten plane cubic | projective closure | smooth genus one | affine boundary | exact singular ideal and infinity gcd |
| sign pair | quotient by simultaneous negation | ratio `r` | sheet orientation | squarefree `A B` |
| `Phi=0` pair | take `r=infinity` fibre | two unramified points | apparent stratum separation | `D(0,eta)=19eta F` |
| row-eleven carrier | pull back `K=0` | degree fourteen | ramification multiplicity | `gcd(K,A B)=1` |
| row-twelve cut | compare quotient divisors | empty intersection | terminal fibre data | `gcd(K,J)=1` |

The resemblance to THM-4378 is structural only. No map has been proved from
its reciprocal packet parity class to this elliptic double cover.

## Next sharp test

The elliptic carrier may be an artifact of the weight-twelve cap. The first
hostile extension is the omitted residual-weight-13 face

```text
c_51 p^5 y + c_23 p^2 y^3.
```

Its labelled channels can enter before row ten. The cheapest decisive
experiment is therefore to rerun the row-nine/ten compiler with these two
parameters and test the normal response of `D`:

- a transverse response destroys the elliptic target as a seam invariant;
- a tangent or zero response preserves it one step and justifies testing
  weights 14 through 20; and
- only stability through every weight that can reach row ten would justify
  trying to glue chartwise extractions into a global map.

If a relevant seam component has rational normalization and the local
extractions glue to a nonconstant morphism into this elliptic curve, then
Riemann--Hurwitz would obstruct that map. That implication is strictly
conditional: no such global seam map, overlap cocycle, or higher-weight
stability theorem exists yet. `JC(2)` and `DC(2)` remain open.
