---
id: THM-3948
title: "No classified irreducible weight-eight or weight-nine sextic has affine-line normalization"
status: >
  PROVED CANDIDATE + CITED CLASSIFICATION; INDEPENDENT AUDIT PENDING.
  In Degtyarev's weight-eight/nine lists, the five rational simple or
  non-simple families have no one-place line by THM-3943 and THM-3945. Every
  other family has positive normalization genus by the exact delta ledger.
  Hence no line cuts an affine chart whose normalized curve is A1. This does
  not claim that a positive-genus row cannot meet some line at one
  normalization point; its one-point puncture would still not be A1.
source: jc-cohn3709 / THM-3943--3945 plus Degtyarev genus ledger, 2026-08-24
depends_on:
  - THM-3943-rational-weight-eight-four-torus-sextics-have-no-one-place-line
  - THM-3945-nonsimple-weight-eight-j-sextics-have-no-one-place-line
related:
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
  - THM-3882-rational-dual-one-place-wronskian-projection-criterion
---

# THM-3948 -- the classified high-character sextic route has no polynomial boundary

**PROVED CANDIDATE + CITED CLASSIFICATION; INDEPENDENT AUDIT PENDING.**
Work over `C`.  Let `C` be an irreducible plane sextic in one of the
weight-eight or weight-nine families classified by Degtyarev, including the
two non-simple weight-eight `J`-families.  Then there is no projective line
`L` for which

```text
C_tilde minus nu^{-1}(L) ~= A1,                         (1)
```

where `nu:C_tilde->C` is the projective normalization.

Thus none of the classified sextics can be a polynomially parametrized
one-place component in an affine chart.  In particular, increasing from one
torus presentation to the four presentations of weight eight, or the twelve
presentations of weight nine, does not by itself produce the boundary
geometry required of a planar Keller nonproperness component.

The classification input is **CITED** from
[Degtyarev, *Irreducible plane sextics with large fundamental groups*](https://arxiv.org/abs/0712.2290),
Sections 1.1 and 3.2--3.6.  The affine-line conclusion is the following exact
repo synthesis.

## 1. The complete genus ledger

A plane sextic has arithmetic genus

```text
p_a=(6-1)(6-2)/2=10.                                   (2)
```

For the simple singularities appearing here,

```text
delta(A1)=delta(A2)=1,        delta(A5)=delta(E6)=3.    (3)
```

Consequently the seven cited simple weight-eight rows split as follows.

| singularities | total delta | normalization genus | gate |
|---|---:|---:|---|
| `E6+A5+4A2` | 10 | 0 | THM-3943 |
| `2A5+4A2` | 10 | 0 | THM-3943 |
| `A5+6A2+A1` | 10 | 0 | THM-3943 |
| `E6+6A2` | 9 | 1 | genus |
| `A5+6A2` | 9 | 1 | genus |
| `8A2+A1` | 9 | 1 | genus |
| `8A2` | 8 | 2 | genus |

The only simple weight-nine row is `9A2`; it has total delta nine and
normalization genus one.  The two non-simple weight-eight rows are

```text
J_{2,0}+4A2,                    J_{2,3}+3A2,             (4)
```

and THM-3945 treats both uniformly without needing their delta invariants.
This exhausts the cited lists.

## 2. Rational rows fail the line test; positive-genus rows fail A1

For each of the three rational simple rows, THM-3943 proves that every line
pullback has at least two normalization support points.  THM-3945 proves the
same statement for both non-simple rows in `(4)`.  Therefore `(1)` is
impossible in all five rational families.

For every remaining row, `C_tilde` has genus one or two by the table.  Removing
finitely many points does not change its function field or geometric genus.
If `(1)` held, that function field would be `C(t)`, whose smooth projective
model is `P1` and has genus zero.  This contradiction closes all remaining
families.

The distinction in the theorem statement is load-bearing.  Positive genus
alone does **not** prove that `nu^{-1}(L)` cannot be supported at one point: a
positive-genus curve may have a high-contact line.  It proves instead that
even a one-point puncture cannot be `A1`.  The stronger support statement is
used only for the five rational rows, where THM-3943 and THM-3945 establish it
exactly.

## 3. Boundary of the conclusion

This theorem closes the classified weight-eight/nine plane-sextic route to a
polynomial one-place nonproperness component.  It does not classify all torus
sextics of lower weight, arbitrary degree-six curves, non-line affine
boundaries, birational target modifications, or the realization of any curve
as the actual nonproperness set of a Keller map.  It does not prove `JC(2)`.

