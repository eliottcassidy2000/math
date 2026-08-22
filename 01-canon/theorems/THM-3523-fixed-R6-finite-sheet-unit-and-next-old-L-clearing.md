---
id: THM-3523
title: "Fixed-R6 finite-sheet unit and next old-L denominator clearing"
status: >
  PROVED + VERIFIED-EXACT + SPLIT-REPRESENTATION AUDITED.  For the fixed
  sporadic Keller map, R_6=L^1699 N(R_5) is nonzero at the canonical finite
  inverse point q=(2,5/6,-7/8) over L=0.  Hence v_L(N(R_6))=-10663, so
  R_7:=L^10663 N(R_6) is polynomial and coprime to L.  Applying proved
  THM-3522 only after this polynomiality gate gives the complete fixed-chart
  packet A(66907,24255) for R_7.  THM-3528 subsequently proves raw all-level
  polynomial packets.  No image equation, irreducibility, all-level unit law,
  arbitrary-map statement, or general Jacobian-conjecture claim follows.
source: codex/fixed-R6-finite-sheet/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
  - THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate
  - THM-3521-fixed-R5-finite-sheet-unit-and-next-old-L-clearing
  - THM-3522-fixed-keller-five-face-renewal-propagation
related:
  - THM-3504-level-four-sporadic-keller-image-prime-and-four-component-nonproperness
  - THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary
  - THM-3513-fixed-G-hybrid-newton-renewal-faces
  - THM-3527-fixed-R7-finite-sheet-unit-and-next-old-L-clearing
  - MISTAKE-415
scripts:
  - 04-computation/keller_R6_finite_sheet_recursive_norm_probe_20260816.py
  - 04-computation/keller_R6_finite_sheet_split_branch_audit_20260816.py
outputs:
  - 05-knowledge/results/keller_R6_finite_sheet_recursive_norm_probe_20260816.out
  - 05-knowledge/results/keller_R6_finite_sheet_split_branch_audit_20260816.out
script_sha256:
  - 9cacaaa825d9598556f7e8baa35970d30295a139c79ed2dc22c9c5e84ef611e8
  - c6ab159d1ab0ba241306d11e3fbc1273f6579c9e843fd0a00994d8107ca01069
output_sha256:
  - 759057c3f06b77f972f551238891049e261b56d1a2619b786d5230c6d8269d6f
  - cdd902f60147eca39bb708d86e4068f049b80efd08b69db14b2f8bcff16eaa52
semantic_sha256:
  - a2ede01095e73ad727285743b83d6502ff37c0ab772a19e0a03fe9036ba5f7b8
  - 6384038b2f7e386289c7af89ac1c5a020aa18320e6957a59469ff013de396490
hash_basis: raw LF bytes for files; ordered finite-field gate ledgers as printed
---

# THM-3523 -- the next finite sheet is again a unit

**PROVED + VERIFIED-EXACT + SPLIT-REPRESENTATION AUDITED.**

Retain the fixed sporadic Keller map `F:C^3->C^3`, the irreducible target
polynomial `L`, and the cubic function-field norm `N` of THM-2473 and
THM-3495.  Use the proved normalizations

```text
H   =2^6 L N(L),
J   =2^35 L^7 N(H),
G   =L^43 N(J),
R_5 =L^271 N(G),
R_6 =L^1699 N(R_5).                                  (1)
```

THM-3521 proves that `R_6` is polynomial and coprime to `L`.  THM-3522,
applied to that separately proved polynomial, gives its complete five-face
packet

```text
R_6 has A(10663,3867),                                (2)
```

including the complete divergent-sheet face

```text
in_max-lambda(R_6)
 =C_6 x^10663(3xz-2y)^3867,       C_6!=0.             (3)
```

This theorem closes the next finite-sheet and denominator-clearing gate.

## 1. The theorem

At the canonical finite inverse point

```text
p=(2/27,1,1),       L(p)=0,
q=(2,5/6,-7/8),     F(q)=p,                           (4)
```

one has

```text
R_6(q)!=0.                                            (5)
```

Consequently, at the generic divisor `(L)`,

```text
v_L(N(R_6))=-10663.                                  (6)
```

Therefore

```text
R_7:=L^10663 N(R_6) belongs to Q[a,b,c],
gcd(R_7,L)=1.                                         (7)
```

Only now may THM-3522 be applied to `P=R_6` and `Q=R_7`.  It gives

```text
R_7 has the complete packet A(66907,24255),           (8)
```

because

```text
(66907,24255)
 =(7*10663-2*3867, 3*10663-2*3867).                  (9)
```

## 2. Divergent sheets and the exact valuation

At the generic DVR of `(L)`, THM-3498 gives one finite inverse root and two
roots of valuation `-1/2`.  On either divergent sheet, with `u=1/w`,

```text
x=u^-1,
y=D/S+O(u),
z=-3(D/S)u+O(u^2),
3xz-2y=-11D/S+O(u),                                  (10)
```

where `D/S` is a unit and `v_L(u)=1/2`.  Substitution in the complete face
(3) gives

```text
v_L(R_6(q_div,+))=v_L(R_6(q_div,-))=-10663/2.         (11)
```

Completeness is load-bearing: every lower-`lambda` monomial has strictly
higher `u`-order, while the full leading form (3) evaluates to a nonzero
unit.  The two divergent sheets therefore contribute `-10663` to the norm.
The only possible correction is a positive valuation on the finite sheet.

The finite branch is regular over the generic point of `(L)`.  If `R_6`
vanished identically there, every lawful specialization of that branch,
including (4), would vanish.  Equation (5) therefore makes its generic
valuation zero.  Adding all three sheet valuations proves (6).

## 3. The one-rung recursive identity

Let `P_0=L`, `P_1=H`, `P_2=J`, `P_3=G`, `P_4=R_5`, and `P_5=R_6`.
Multiplicativity of the degree-three norm unrolls (1), on the finite-etale
locus, to the exact rational-function identity

```text
R_6
 =2^1431 L^1699 N(L)^271 N^2(L)^43
    *N^3(L)^7 N^4(L) N^5(L).                         (12)
```

The scalar is not extrapolated from a face recurrence.  It is exactly

```text
1431=35*3^3+6*3^4.                                   (13)
```

At a good prime, start with `A_0=F_p`.  If a target point lies in a finite
etale algebra `A_r`, form

```text
A_(r+1)=A_r[w]/(w^3+(T/L)w-2c/L),                    (14)
```

construct the universal inverse point from THM-3495's exact formulas, and
recurse.  Every division by `L` or `S`, and every cubic discriminant, is
certified as a unit by the full regular representation.  Direct substitution
checks `F(source)=target` at every level.

The primary route uses the complete dimensions

```text
1 -> 3 -> 9 -> 27 -> 81 -> 243                       (15)
```

and evaluates the five-term `L` on all `243` terminal sheets.  The six
absolute norm factors in (12) are

```text
p=101: (16,12,72,9,49,97),
p=103: (12,53,22,85,76,94),
p=107: (38,45,28,3,17,17).                           (16)
```

Every entry is nonzero.  Substitution in (12) gives

| prime | `R_6(q)` | frozen-`H` flat norm | omit `64` |
|---:|---:|---:|---:|
| `101` | `26` | `37` | `67` |
| `103` | `70` | `31` | `47` |
| `107` | `69` | `13` | `3` |

Thus any one of the first-column residues proves (5); all three are retained
as independent good-reduction controls.

For completeness, the absolute norm triples
`(N(L),N(S),N(discriminant))` at base dimensions `1,3,9,27,81` are

```text
p=101: (16,78,56), (12,2,41), (72,36,10), (9,9,13), (49,75,58),
p=103: (12,72,91), (53,5,60), (22,1,79), (85,82,68), (76,73,39),
p=107: (38,86,14), (45,105,34), (28,31,100), (3,95,82), (17,72,100).
                                                               (17)
```

All `45` displayed gates are nonzero.

## 4. Independent bottom representation and hostile normalization

A second route stops at `P_1=H`.  It uses only

```text
1 -> 3 -> 9 -> 27 -> 81                              (18)
```

and evaluates the frozen `361`-term polynomial `H` in the bottom algebra.
The transitive cubic norm of this element agrees, at each prime, with the
determinant of its literal `81 by 81` multiplication matrix.  Those
determinants are the third column of the table above.  Thus the bottom
quantity is represented once as `64 L N(L)` and once as the global frozen
polynomial `H`.

Deleting the factor `64` from the first transition changes the final value
by exactly

```text
64^-81,                                               (19)
```

because `H` occurs on all `81` leaves.  The last column records the hostile
values; each differs from the lawful value.  This detects both the nonmonic
normalization and the correct sheet multiplicity.

## 5. Split-outer representation at a fourth prime

The disjoint outer representation does not form one `243`-dimensional
algebra.  Over `F_71`, the inverse cubic above `q` splits into the three
distinct roots

```text
w=10,23,38.                                           (20)
```

The exact inverse points and complete `81`-sheet `R_5` values are

| `w` | `y` | `z` | `R_5(w,y,z)` | bottom `N^4(L)` |
|---:|---:|---:|---:|---:|
| `10` | `0` | `2` | `49` | `6` |
| `23` | `1` | `36` | `22` | `5` |
| `38` | `18` | `60` | `60` | `54` |

Every inverse point is checked directly against `F`.  The twelve inner
unit/discriminant triples, grouped by branch, are

```text
w=10: (26,23,4), (24,26,35), (57,53,66), (62,64,32),
w=23: (57,70,68), (34,34,19), (38,15,26), (40,47,35),
w=38: (20,53,44), (34,56,54), (44,21,29), (36,69,13).
                                                               (21)
```

Since `L(q)=53 mod 71`, direct branch multiplication gives

```text
53^1699 * 49 * 22 * 60 = 9 mod 71.                   (22)
```

The unsplit nested representation also gives `9`.  Omitting the named
`w=38` branch gives `25`, so a proper subuniverse cannot accidentally pass
the norm check.  This split route realizes the same `3*81=243` sheets but
with a scalar outer branch product.

## 6. Denominator clearing and packet propagation

All rational denominators in the recursive charts are units at the four
primes above.  Each residue is therefore a reduction of the same rational
number `R_6(q)`, proving (5).

On

```text
U=Spec(Q[a,b,c,L^-1]),                                (23)
```

THM-2473 makes the inverse cover finite etale.  Hence
`N(R_6)` belongs to `Q[a,b,c,L^-1]`.  The ring `Q[a,b,c]` is a UFD and `L`
is irreducible.  Equation (6) therefore says that the reduced denominator is
exactly `L^10663`; clearing it proves both assertions in (7).  There is no
unexamined `S`, discriminant, or finite-sheet denominator.

THM-3522 is conditional on polynomiality of the cleared norm.  Equation (7)
supplies precisely that missing hypothesis, and only then yields (8)--(9).
No face argument is being used to infer polynomiality.

## 7. Exact boundary

This theorem proves one fixed-map finite-sheet value, one exact old-`L`
valuation, polynomiality and `L`-coprimality of `R_7`, and—through proved
THM-3522—the complete fixed-chart packet of `R_7`.  It does **not** prove:

- irreducibility, squarefreeness, or an integral primitive normalization of
  `R_5`, `R_6`, or `R_7`;
- that any of these polynomials is a new image equation;
- a fifth, sixth, or seventh nonproperness component;
- the still-open degree-`243` eliminant separability/image gate (the
  `243` finite sheets in (15) are an evaluation algebra, not that gate);
- a finite-sheet unit, L-coprimality, or image status at every later rung
  (THM-3528 subsequently proves raw polynomiality only);
- an arbitrary-map norm theorem, a classification of Keller maps, `JC(2)`,
  `DC(2)`, LRC, or any general Jacobian-conjecture conclusion.

Reproduce the two deterministic certificates with

```text
python 04-computation/keller_R6_finite_sheet_recursive_norm_probe_20260816.py
python -O 04-computation/keller_R6_finite_sheet_recursive_norm_probe_20260816.py
python 04-computation/keller_R6_finite_sheet_split_branch_audit_20260816.py
python -O 04-computation/keller_R6_finite_sheet_split_branch_audit_20260816.py
```

Ordinary and optimized transcripts agree exactly with the stored outputs.

**QED.**
