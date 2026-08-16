---
id: THM-3527
title: "Fixed-R7 finite-sheet unit and next old-L denominator clearing"
status: >
  PROVED + VERIFIED-EXACT + SPLIT-REPRESENTATION AUDITED.  For the fixed
  sporadic Keller map, R_7=L^10663 N(R_6) is nonzero at the canonical finite
  inverse point q=(2,5/6,-7/8) over L=0.  Hence v_L(N(R_7))=-66907, so
  R_8:=L^66907 N(R_7) is polynomial and coprime to L.  Applying THM-3522
  only after this polynomiality gate gives the complete fixed-chart packet
  A(419839,152211), exactly the next Pell-57 recurrence row.  No image
  equation, irreducibility, later finite-sheet unit, all-level L-coprime/image
  law, arbitrary-map statement, or general Jacobian-conjecture claim follows.
  THM-3528 subsequently proves the weaker raw all-level polynomial-packet law.
source: codex/fixed-R7-finite-sheet/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
  - THM-3522-fixed-keller-five-face-renewal-propagation
  - THM-3523-fixed-R6-finite-sheet-unit-and-next-old-L-clearing
related:
  - THM-3525-level-five-degree243-separability-and-discriminant-square-class
  - THM-3526-level-six-degree729-separability-and-discriminant-square-class
scripts:
  - 04-computation/keller_R7_finite_sheet_recursive_norm_probe_20260816.py
  - 04-computation/keller_R7_finite_sheet_split_outer_independent_audit_20260816.py
outputs:
  - 05-knowledge/results/keller_R7_finite_sheet_recursive_norm_probe_20260816.out
  - 05-knowledge/results/keller_R7_finite_sheet_split_outer_independent_audit_20260816.out
script_sha256:
  - 80d7a5e105be82d29fb1c3caf28ec0d4d98ea1e51c02bbf080125ef812133463
  - 623a7817ba233d3c3c5e0dfbe218fdb93fb27c3a034d882a25f10e401228a3fd
output_sha256:
  - 100371c34090af4d9e613696b591770e9092e4f199e28e55c14095acd6f1b271
  - 6189cad420d7675eeb0a623ed4e62855d0cc3b14b671fce195eed67c3ac87c8f
semantic_sha256:
  - 82efb24e0c4a6e0df9671f0f5a5009dd0e77d1b0aa8ef2341780dfe23ea28c38
  - d8fd9ddf1f8679b90434d1f6ffa6a717c1725e3dcb5703c040e1ab724081e72b
hash_basis: LF-normalized bytes; ordered finite-field gate ledgers as printed
---

# THM-3527 -- the next finite old-L sheet is again a unit

**PROVED + VERIFIED-EXACT + SPLIT-REPRESENTATION AUDITED.**

Retain the fixed sporadic Keller map `F:C^3->C^3`, the irreducible target
polynomial `L`, and the cubic function-field norm `N` of THM-2473 and
THM-3495.  The proved cleared norms now reach

```text
R_6=L^1699 N(R_5),
R_7=L^10663 N(R_6).                                  (1)
```

THM-3523 proves that `R_7` is polynomial and coprime to `L`; applying
THM-3522 there gives its complete packet

```text
R_7 has A(66907,24255).                               (2)
```

This theorem closes the next finite-sheet and denominator-clearing gate.

## 1. The theorem

At the canonical finite inverse point

```text
p=(2/27,1,1),       L(p)=0,
q=(2,5/6,-7/8),     F(q)=p,                           (3)
```

one has

```text
R_7(q)!=0.                                            (4)
```

Consequently, at the generic divisor `(L)`,

```text
v_L(N(R_7))=-66907.                                  (5)
```

Therefore

```text
R_8:=L^66907 N(R_7) belongs to Q[a,b,c],
gcd(R_8,L)=1.                                         (6)
```

Only now may THM-3522 be applied to `P=R_7` and `Q=R_8`.  It gives

```text
R_8 has the complete packet A(419839,152211),         (7)
```

because

```text
(419839,152211)
 =(7*66907-2*24255,3*66907-2*24255).                 (8)
```

## 2. Divergent sheets and the exact valuation

At the generic DVR of `(L)`, the inverse cubic has one regular finite root and
two roots of valuation `-1/2`.  On either divergent sheet, writing `u=1/w`,

```text
x=u^-1,
y=D/S+O(u),
z=-3(D/S)u+O(u^2),
3xz-2y=-11D/S+O(u),                                  (9)
```

where `D/S` is a unit and `v_L(u)=1/2`.

The complete maximum-`lambda` face supplied by (2) is

```text
in_max-lambda(R_7)
 =C_7 x^66907(3xz-2y)^24255,       C_7!=0.           (10)
```

Substitution in (10) gives

```text
v_L(R_7(q_div,+))=v_L(R_7(q_div,-))=-66907/2.        (11)
```

Completeness is load-bearing: every lower-`lambda` monomial has strictly
higher `u`-order, while the displayed leading form evaluates to a nonzero
unit.  The two divergent sheets contribute `-66907` to the norm.

The finite branch is regular over the generic point of `(L)`.  If `R_7`
vanished identically there, every lawful specialization of that branch,
including (3), would vanish.  Equation (4) therefore gives valuation zero on
the finite sheet, proving (5).

## 3. Complete recursive norm orbit

Norm multiplicativity unrolls the fixed definitions to

```text
R_7
 =2^4293 L^10663 N(L)^1699 N^2(L)^271 N^3(L)^43
          N^4(L)^7 N^5(L) N^6(L).                   (12)
```

The scalar exponent is exact:

```text
4293=35*3^4+6*3^5.                                   (13)
```

The primary route builds the six complete inverse-algebra dimensions

```text
1 -> 3 -> 9 -> 27 -> 81 -> 243 -> 729.               (14)
```

Recursive cubic adjugates certify every division.  At each level, direct
substitution checks the inverse graph, and the regular norms of `L`, the
cubic derivative, the y-chart denominator, and `x^3` are all nonzero.  The
three complete norm orbits are

```text
p=101: (16,12,72, 9,49,97,71),
p=103: (12,53,22,85,76,94,17),
p=107: (38,45,28, 3,17,17,30).                       (15)
```

Substitution in (12) gives

```text
R_7(q) mod (101,103,107)=(72,44,53).                 (16)
```

Every one of the `72` displayed unit gates is nonzero.  Any single residue in
(16) proves (4); all three are retained as independent good-reduction
controls.

Deleting the factor `64` from `H=64LN(L)` changes the three values to

```text
(66,78,19).                                           (17)
```

The ratio is exactly `64^-243`, because the bottom `H` occurs on all 243
leaves.  Equivalently, the scalar exponent in (12) changes from `4293` to
`2835`.  This detects the nonmonic normalization and the complete leaf
multiplicity.

## 4. Independent split-outer representation

The second route never constructs a 729-dimensional algebra.  Over `F_71`,
the outer inverse cubic above `q` splits into the three distinct roots

```text
w=10,23,38.                                           (18)
```

The corresponding inverse points and complete 243-sheet `R_6` values are

| `w` | `y` | `z` | `R_6(w,y,z)` | frozen-`H` flat norm |
|---:|---:|---:|---:|---:|
| `10` | `0` | `2` | `56` | `23` |
| `23` | `1` | `36` | `65` | `64` |
| `38` | `18` | `60` | `10` | `50` |

Each inverse point is checked directly against `F`.  Above each named point,
the first representation descends through five coefficient-algebra layers to
`L`.  The second stops at the frozen 361-term polynomial `H` in dimension 81;
its transitive norm agrees with the determinant of the literal `81 by 81`
multiplication matrix.  Thus every branch value has two bottom
representations.

Since `L(q)=53 mod 71`, the split product is

```text
53^10663 * 56 * 65 * 10 = 56 mod 71.                 (19)
```

It is nonzero, independently proving (4).  Omitting the named `w=38` branch
changes (19) to `34`; deleting the `H` normalization changes it to `46`.
The complete `3*243=729` sheet universe is therefore retained in a geometry
disjoint from the primary nested algebra.

## 5. Denominator clearing and packet propagation

All rational denominators in both recursive charts are units at the four
primes above.  Each residue is a reduction of the same rational number
`R_7(q)`, proving (4).

On

```text
U=Spec(Q[a,b,c,L^-1]),                                (20)
```

the inverse cover is finite etale, so `N(R_7)` belongs to
`Q[a,b,c,L^-1]`.  Since `Q[a,b,c]` is a UFD and `L` is irreducible, (5) says
that the reduced old-`L` denominator is exactly `L^66907`.  Clearing it proves
both assertions in (6).  There is no unexamined `S`, discriminant, or
finite-sheet denominator.

Only after (6) is established does THM-3522 apply, yielding (7)--(8).  The
Pell-57 multiplication law predicted the same pair before this computation,
but did not supply (4)--(6).  This theorem is a prediction/realization match
and not by itself an all-level induction; THM-3528 subsequently supplies the
raw polynomial-packet induction.

## 6. Exact boundary and reproduction

This theorem proves one fixed-map finite-sheet value, one exact old-`L`
valuation, polynomiality and `L`-coprimality of `R_8`, and the complete packet
(7).  It does **not** prove:

- irreducibility, squarefreeness, or image-equation status of any `R_i`;
- a fifth, sixth, seventh, or eighth new nonproperness component;
- a later degree/separability gate or finite-sheet unit;
- an all-level finite-unit, image, or discriminant theorem (THM-3528
  subsequently supplies raw polynomial packets only);
- an arbitrary-map norm theorem, classification of Keller maps, `JC(2)`,
  `DC(2)`, LRC, or any general Jacobian-conjecture conclusion.

Reproduce the two exact certificates with

```text
python -B 04-computation/keller_R7_finite_sheet_recursive_norm_probe_20260816.py
python -B -O 04-computation/keller_R7_finite_sheet_recursive_norm_probe_20260816.py
python -B 04-computation/keller_R7_finite_sheet_split_outer_independent_audit_20260816.py
python -B -O 04-computation/keller_R7_finite_sheet_split_outer_independent_audit_20260816.py
```

**QED.**
