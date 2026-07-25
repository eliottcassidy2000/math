---
id: THM-2187
title: "Height-preserving saturation and a universal radix carrier"
status: >
  PROVED. Let m,u be independent integer vectors, with m primitive. The
  saturation of Zm+Zu has a basis (m,n) with
  ||n||_infinity<=max(||m||_infinity,||u||_infinity). Its row map is integrally
  split, hence surjective modulo every integer q>=2, with unrestricted
  fibres of size q^(d-2). Applied to the height-29/43-or-105 LRC relation
  plane, this removes adaptive prime selection without increasing the
  sign-sharp 208875 carry-pair bound: base 14, or any fixed composite bank,
  is available. It also sharpens THM-2178's saturated-basis height from
  286650 to 105. With THM-2185, every zero-safe row moreover has a saturated
  rank-three relation basis of height at most 500, surjective modulo every
  base with unrestricted fibres q^10. The remaining pumping loss is
  phase/current, not radix rank.
source: codex-2026-07-24-LRC-saturated-relation-plane
depends_on:
  - THM-2144
  - THM-2163
  - THM-2164
  - THM-2171
  - THM-2178
  - THM-2185
related:
  - THM-2167
  - THM-2174
  - THM-2184
  - THM-2188
---

# THM-2187 -- height-preserving saturation

The adaptive prime in THM-2167 is sufficient but unnecessary. The missing
operation is to saturate the relation plane while retaining one primitive
short vector.

## 1. A two-dimensional saturation lemma

Let `m,u in Z^d` be rationally independent. Assume

```text
m is primitive,
H_m=||m||_infinity,              H_u=||u||_infinity. (1)
```

Put

```text
L=(Qm+Qu) intersection Z^d.                           (2)
```

Then `L` has a basis `(m,n)` satisfying

```text
||n||_infinity<=max(H_m,H_u).                         (3)
```

Indeed, `m` is primitive in `L`: a nontrivial divisibility in `L` would be
the same divisibility in `Z^d`. Extend it to a basis `(m,n_0)` of `L` and
write, after changing the sign of `n_0` if necessary,

```text
u=a m+D n_0,                         D>=1.             (4)
```

Replacing `n_0` by `n=n_0+k m` replaces `a` by `a-Dk`. Choose `k` so that

```text
|a-Dk|<=D/2.                                           (5)
```

Rename the centered coefficient `a`. If `D=1`, choose `a=0`, so `n=u`.
If `D>=2`, equation (4) gives

```text
||n||_infinity
 <=H_u/D+(|a|/D)H_m
 <=H_u/2+H_m/2
 <=max(H_m,H_u).                                      (6)
```

This proves (3). Notice that the conclusion retains the chosen anchor
`m`; a generic reduced-basis estimate does not.

Both basis vectors are primitive in `Z^d`. They are primitive in `L`
because they belong to a lattice basis, and saturation of `L` promotes
primitivity to the ambient lattice.

## 2. Universal residue fibres

Let `B` be the `2 x d` matrix with rows `m,n`. Since its row lattice `L`
is saturated, Smith normal form has both invariant factors equal to one.
Equivalently, the gcd of the `2 x 2` minors is one, and there is an integral
right inverse

```text
B R=I_2.                                               (7)
```

Reduction of (7) modulo any integer `q>=2` proves that

```text
(Z/qZ)^d -> (Z/qZ)^2,             x |-> Bx            (8)
```

is surjective. Every fibre therefore has exactly

```text
q^(d-2)                                               (9)
```

elements. For prime `q`, the two rows are independent over `F_q`.

There is also an exact live-owner statement. For `O subset {1,...,d}`, let
`B_O` be the restriction to those columns and let `I_O(q)` be its image.
Every compatible owner-supported fibre has size

```text
q^|O|/|I_O(q)|.                                      (10)
```

For `q=14`, Chinese remaindering gives the concrete form

```text
2^(|O|-rho_2(O)) 7^(|O|-rho_7(O)),                   (11)
```

where `rho_p(O)` is the rank of `B_O` modulo `p`. Thus a composite base
does not require a fictitious single field rank.

## 3. Fixed-base LRC relation carrier

Let `v=(v_1,...,v_13)` be a distinct positive zero-safe row. THM-2144 gives
a primitive relation `m` of height at most `29`. THM-2164 gives an
independent companion `u` with

```text
||u||_infinity<=43       if |supp(m)|>=3,
||u||_infinity<=105      if |supp(m)|=2.              (12)
```

Apply Section 1 inside the relation lattice. We obtain a saturated basis
`(m,n)` of relations with exactly the same two height bounds. Therefore,
for **every** chosen base `q>=2`, the simultaneous digit equations in
THM-2163 have unrestricted fibres of size

```text
q^11.                                                 (13)
```

The carry counts do not depend on the base. In the first branch,
primitivity gives

```text
||m||_1<=12*29+28=376,
||n||_1<=12*43+42=558,
```

and hence at most

```text
375*557=208875                                        (14)
```

positive-row carry pairs. In the support-two branch,

```text
||m||_1<=57,             ||n||_1<=1364,
```

so there are at most

```text
56*1363=76328                                             (15)
```

pairs. Thus (14) remains universal. Sorting still gives at most
`14*208875=2924250` carry-owner states, and THM-2171's quotient-tie sidecar
still gives the ordered cap

```text
26*208875=5430750.                                    (16)
```

In particular, base `14` may be fixed from the outset. Its lowest digit is
the complete denominator-`14` residue word, while (11) describes every
later owner fibre exactly. More generally, one may take a composite radix
containing any prescribed finite CRT bank without changing the relation
state counts.

This strictly strengthens the adaptive-prime step of THM-2167. It does not
contradict THM-2161: the radix can retain any fixed finite residue bank, but
no such bank determines the lonely-runner target.

## 4. Sharpened mod-14 rank-harvest fork

THM-2178 began with two height-`105` relations and bounded a saturated basis
by `286650`. Sections 1 and 3 instead provide a saturated relation basis
with

```text
||B||_infinity<=105.                                  (17)
```

Consequently its unconditional fork holds with

```text
H_sat=100*616^22*105^44-1:                            (18)
```

either

```text
dim_Q W_(H_sat)(v)>=3,                                (19)
```

or the saturated plane modulo `14` contains a nonzero word supported on
at most three coordinates. No other part of THM-2178's proof changes.

## 5. Saturated rank three at height 500

THM-2185 now proves unconditionally that a zero-safe distinct positive row
has

```text
dim_Q W_500(v)>=3.                                    (20)
```

Let `(b_1,b_2)` be the saturated basis from Section 3, so

```text
||b_1||_infinity,||b_2||_infinity<=105,
```

and choose a relation `c` outside their rational plane with
`||c||_infinity<=500`. Put

```text
L_3=(Qb_1+Qb_2+Qc) intersection Z^13.                (21)
```

The rank-two lattice `L_2=Zb_1+Zb_2` is primitive in `L_3`: it is already
saturated in `Z^13`, so `L_3/L_2` is a subgroup of the torsion-free group
`Z^13/L_2`. Extend `(b_1,b_2)` to a basis
`(b_1,b_2,n_0)` of `L_3` and write

```text
c=a_1b_1+a_2b_2+D n_0,                 D>=1.          (22)
```

Replace `n_0` by `n_0+k_1b_1+k_2b_2` to center both coefficients:

```text
|a_1|,|a_2|<=D/2.                                    (23)
```

If `D=1`, choose both coefficients zero and take `n=c`. If `D>=2`, then

```text
||n||_infinity
 <=500/D+(|a_1|+|a_2|)105/D
 <=250+105=355.                                      (24)
```

Thus `L_3` has a saturated basis whose first two vectors have height at
most `105` and whose third has height at most `500`. Smith normal form
again supplies an integral right inverse. For every `q>=2`,

```text
(Z/qZ)^13 -> (Z/qZ)^3
```

is surjective, and every unrestricted simultaneous digit fibre has exactly

```text
q^10                                                   (25)
```

words. In base `14`, this is `14^10`; live-owner fibres split exactly into
their mod-`2` and mod-`7` image sizes.

This supersedes the sparse-code terminal in Section 4 as a relation-rank
question. The height-`105` rank-two subcarrier remains useful because its
carry count is much smaller and its first anchor is height `29`.

## 6. Exact boundary

Saturation repairs base selection, unrestricted fibre rank, and the loose
basis height. It does not make a repeated carry state target-preserving.
THM-2171 repairs order, distinctness, and common-gcd normalization.
THM-2174 and THM-2184 show the genuine remaining loss: aligned finite
windows preserve bounded phase data, but higher-order continuation and the
scale-sensitive endpoint term need not repeat.

The faithful carrier is therefore

```text
saturated relation plane
 + optional saturated height-500 third relation
 + exact composite-radix owner image
 + quotient-tie sidecar
 + full phase/current continuation data.             (26)
```

A fixed base is now free; the final sidecar is not. QED.
