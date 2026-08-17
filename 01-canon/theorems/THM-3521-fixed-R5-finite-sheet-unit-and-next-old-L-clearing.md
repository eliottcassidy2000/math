---
id: THM-3521
title: "Fixed-R5 finite-sheet unit and next old-L denominator clearing"
status: >
  PROVED + VERIFIED-EXACT + REPRESENTATION-INDEPENDENTLY AUDITED.  For the
  fixed sporadic Keller map, R_5=L^271 N(G) is nonzero at the canonical
  finite inverse point q=(2,5/6,-7/8) over L=0.  Hence
  v_L(N(R_5))=-1699, so R_6:=L^1699 N(R_5) is polynomial and coprime to L.
  Its exposed top face is C x^10663(3xz-2y)^3867 with C nonzero.  This theorem
  itself does not prove the renewal faces of R_5; THM-3522 subsequently proves
  fixed-chart renewal and closes the complete packets of R_5 and R_6.
  THM-3523 subsequently closes the next finite-sheet/polynomiality gate and
  gives the complete packet of R_7.  Irreducibility or image status of R_5,
  later images, and every general Jacobian-conjecture claim remain open at
  this theorem state.  THM-3525 later closes the degree-243 gate, THM-3528
  proves the raw all-level polynomial-packet law, and THM-3529 proves all
  complete packets are finite-sheet units.
source: codex/fixed-R5-finite-sheet/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
  - THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate
  - THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary
  - THM-3513-fixed-G-hybrid-newton-renewal-faces
related:
  - THM-3504-level-four-sporadic-keller-image-prime-and-four-component-nonproperness
  - THM-3522-fixed-keller-five-face-renewal-propagation
  - THM-3523-fixed-R6-finite-sheet-unit-and-next-old-L-clearing
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
  - MISTAKE-415
scripts:
  - 04-computation/keller_R5_finite_sheet_recursive_norm_probe_20260816.py
  - 04-computation/keller_R5_finite_sheet_split_global_J_audit_20260816.py
outputs:
  - 05-knowledge/results/keller_R5_finite_sheet_recursive_norm_probe_20260816.out
  - 05-knowledge/results/keller_R5_finite_sheet_split_global_J_audit_20260816.out
script_sha256:
  - a201191410e39d47fbf607191e8bd597453c697f134d3803466694b680d8c60d
  - 1a46b961ab15a61e6e438926c8834fff56c23749f64ed2473e626bebb2fd1d04
output_sha256:
  - 14473ff317e30fb8a90b7d1c0c3879537ce157d7456d31e8ce6901a160d6197f
  - 42ddc36acea979ab859c10787da7a7737362c31a16ca0383cf8e06b8fe27ef2f
semantic_sha256:
  - e4610844c8bd506211662c792f27d5bac2529f9d62d63762aaf39e498c4b8707
  - cd99969bd8949ab971cc9f7ee3fefac8aba835b2a6b4fdd323934ef4294d0589
hash_basis: raw LF bytes for files; ordered finite-field gate ledgers as printed
---

# THM-3521 -- the next finite sheet is a unit

**PROVED + VERIFIED-EXACT + REPRESENTATION-INDEPENDENTLY AUDITED.**

Let `F:C^3->C^3` be the fixed sporadic Keller map of THM-2473.  Retain the
target polynomials and cubic function-field norm `N` used in THM-3495:

```text
L=27a^2c^2-18abc+16a+b^3c-b^2,
T=4-3bc,
S=27ac^2-9bc+8,
E(w)=Lw^3+Tw-2c.
```

Use the proved normalizations

```text
H   =2^6 L N(L),
J   =2^35 L^7 N(H),
G   =L^43 N(J),
R_5 =L^271 N(G).                                      (1)
```

THM-3506 and THM-3513 prove that `G` has the complete packet `A(271,99)`.
They consequently prove that `R_5` is polynomial, is coprime to `L`, and
has the three transported faces at `(e,m)=(1699,615)`, including

```text
in_max-lambda(R_5)=C_5 x^1699(3xz-2y)^615,  C_5!=0,   (2)
```

where `lambda=i-k`.  The theorem below closes the next finite-sheet gate.

## 1. The theorem

At the canonical finite inverse point

```text
p=(2/27,1,1),       L(p)=0,
q=(2,5/6,-7/8),     F(q)=p,                            (3)
```

one has

```text
R_5(q)!=0.                                             (4)
```

Therefore, at the generic divisor `(L)`,

```text
v_L(N(R_5))=-1699.                                    (5)
```

In particular,

```text
R_6:=L^1699 N(R_5) belongs to Q[a,b,c],
gcd(R_6,L)=1.                                          (6)
```

Moreover, the complete opposite-`lambda` and minimum-`beta` faces already
transported to `R_5`, together with (6), give

```text
in_max-lambda(R_6)
 =C_6 x^10663(3xz-2y)^3867,       C_6!=0,              (7)
```

because

```text
(10663,3867)
 =(7*1699-2*615, 3*1699-2*615).                        (8)
```

Equation (7) is only the next exposed top face.  It does not assert that
`R_5` or `R_6` has a complete five-face packet.

## 2. The two divergent sheets

At the generic DVR of `(L)`, THM-3495 gives one finite inverse root and two
roots of valuation `-1/2`.  On either divergent sheet, with `u=1/w`,

```text
x=u^-1,
y=D/S+O(u),
z=-3(D/S)u+O(u^2),
3xz-2y=-11D/S+O(u),                                  (9)
```

where `D/S` is a unit.  Substitution of (9) in the complete face (2) gives

```text
v_L(R_5(q_div))=-1699/2                               (10)
```

on each divergent sheet.  Completeness of (2) rules out cancellation among
equal-weight monomials inside either branch.  Norm valuations add, so these
two branches contribute `-1699` in total.  The only remaining question for
(5) is whether the finite branch contributes a positive valuation.

## 3. A cheap recursive finite-algebra certificate

The point in (3) also satisfies

```text
L(q)=241465/1728!=0.                                   (11)
```

Thus `R_5(q)` can be evaluated entirely inside the finite-etale locus.  At a
good prime, start with `A_0=F_p`.  Given a target point with coordinates in a
finite etale `A_r`, form

```text
A_(r+1)=A_r[w]/(w^3+(T/L)w-2c/L),                     (12)
```

and use the exact inverse formulas of THM-3495 for the universal preimage in
`A_(r+1)^3`.  At every level the companion checks that `L`, `S`, and the
cubic discriminant are units by inverting their full regular-representation
matrices.  It also substitutes the universal point directly into `F`.

Norm multiplicativity also unrolls (1) to the exact localized identity

```text
R_5=2^477 L^271 N(L)^43 N^2(L)^7 N^3(L) N^4(L).       (13)
```

This is an identity of rational functions on the finite-etale locus, not a
Newton-face recurrence.  The companion evaluates its five norm-orbit factors
separately and recovers the recursive value.

The primary route applies (1) recursively through dimensions

```text
1 -> 3 -> 9 -> 27 -> 81                               (14)
```

and stops at the five-term polynomial `L`.  A second route uses dimensions

```text
1 -> 3 -> 9 -> 27                                     (15)
```

and evaluates the frozen `361`-term polynomial `H` directly at the bottom.
The two routes agree at three primes:

| prime | `R_5(q) mod p` | flat norm of bottom `H` | omit the factor `64` |
|---:|---:|---:|---:|
| `101` | `74` | `49` | `60` |
| `103` | `36` | `87` | `91` |
| `107` | `88` | `20` | `96` |

The corresponding rows `(L,N(L),N^2(L),N^3(L),N^4(L))` are

```text
p=101: (16,12,72,9,49),
p=103: (12,53,22,85,76),
p=107: (38,45,28,3,17).                               (16)
```

Every factor in (16) is nonzero, and substitution in (13) gives the second
column of the table.

All entries in the second column are nonzero.  The third column is computed
twice: by three successive cubic norms and by the determinant of the literal
`27 by 27` multiplication matrix.  The fourth column is a hostile control.
Deleting the factor `64` in `H=64LN(L)` changes the answer by exactly
`64^-27`, as required because the normalization occurs on `27` leaves, and
changes the value at all three primes.

For completeness, the absolute norm triples `(N(L),N(S),N(disc))` at base
dimensions `1,3,9,27` are

```text
p=101: (16,78,56), (12,2,41), (72,36,10), (9,9,13),
p=103: (12,72,91), (53,5,60), (22,1,79), (85,82,68),
p=107: (38,86,14), (45,105,34), (28,31,100), (3,95,82). (17)
```

Every entry is nonzero.  These are good-reduction certificates for every
division and every finite-etale extension used in the primary route.

## 4. A representation-disjoint split-`J` audit

The independent audit does not use the recursive descent to `L` or the
explicit `H` route.  It reconstructs THM-3495's full `66,146`-term `J`.
Over `F_71`, the outer inverse cubic above `q` splits with roots

```text
w=10,23,38.                                            (18)
```

The corresponding inverse points and inner ledgers are listed as
`(w,y,z,N(J),L,disc,G)`:

```text
(10, 0,  2, 68, 26,  4, 64),
(23, 1, 36, 14, 57, 68, 66),
(38,18, 60, 40, 20, 44, 19).                           (19)
```

Each point maps directly to `q`; every displayed unit and `G` value is
nonzero.  Branchwise multiplication gives

```text
L(q)^271 product_(F(r)=q) G(r)=43 mod 71.              (20)
```

This agrees with an additional run of the recursive companion at `p=71`,
but it reaches the value through a different global polynomial and a split
branch product.

## 5. Why the reductions prove the valuation

All rational denominators in the displayed inverse charts and
normalizations are units at the four primes used above.  Therefore each
finite-field value is the reduction of the same rational number `R_5(q)`.
One nonzero value would prove (4); equations (16)--(17), (19), and (20) provide
four controlled witnesses by two representations.

The finite branch over the generic point of `(L)` is regular.  If `R_5`
vanished identically on that branch, every lawful specialization such as
(3) would vanish.  Equation (4) therefore proves that the finite branch is
generically a unit.  Combining its valuation zero with the two values in
(10) proves (5).

On `U=Spec(Q[a,b,c,L^-1])`, THM-2473 makes the inverse cover finite etale,
so `N(R_5)` belongs to `Q[a,b,c,L^-1]`.  Since `Q[a,b,c]` is a UFD and `L`
is irreducible, equation (5) says its reduced denominator is exactly
`L^1699`; clearing it gives (6), including coprimality.

Finally, THM-3506's top-face argument uses only the complete
minimum-`lambda` and minimum-`beta` faces of its input plus polynomiality of
the cleared norm.  Those two faces of `R_5` were already transported from
the complete packet of `G`.  Applying that argument with (6) proves (7)--
(8) without assuming either missing renewal face of `R_5`.

## 6. Exact boundary

This theorem closes exactly one gate: the finite old-`L` sheet at the next
rung.  By itself it does **not** prove any of the following:

- the `z`-top or minimum-`gamma` renewal face of `R_5`;
- a complete packet `A(1699,615)` for `R_5`;
- irreducibility, squarefreeness, or a primitive integral normalization of
  `R_5`;
- `closure(F(V(G)))=V(R_5)` or a fifth nonproperness component;
- full degree or squarefreeness of the degree-`243` eliminant;
- a sixth image component attached to `R_6`;
- an all-level L-coprime/image recurrence, arbitrary-map classification, `JC(2)`,
  `DC(2)`, or any general Jacobian-conjecture conclusion.

THM-3522 subsequently closes the first two bullets by proving the complete
packet `A(1699,615)` for `R_5`, and also proves `A(10663,3867)` for `R_6`.
THM-3523 subsequently closes `L^10663N(R_6)`, proves `R_7` polynomial and
`L`-coprime, and gives its complete packet `A(66907,24255)`.  THM-3528 later
closes raw polynomial packets at all levels without changing the geometric,
finite-unit, or image bullets.

Reproduce both exact routes with

```text
python 04-computation/keller_R5_finite_sheet_recursive_norm_probe_20260816.py
python -O 04-computation/keller_R5_finite_sheet_recursive_norm_probe_20260816.py
python 04-computation/keller_R5_finite_sheet_split_global_J_audit_20260816.py
python -O 04-computation/keller_R5_finite_sheet_split_global_J_audit_20260816.py
```

Normal and optimized outputs agree line-for-line with the stored transcripts
after LF normalization.

**QED.**
