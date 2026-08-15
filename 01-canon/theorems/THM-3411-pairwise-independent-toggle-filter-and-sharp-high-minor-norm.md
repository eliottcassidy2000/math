---
id: THM-3411
title: "Pairwise-independent toggle filter and sharp high-minor norm"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF/REPLAY-AUDITED
source: root-2608-crouzeix-puzzle-2026-08-15
audit: combined independent Walsh-sign, demicube-norm, blindness, hostile, normal/-O/stored, hash, dependency, and scope audit clean
depends_on:
  - THM-3396-four-bit-pairwise-independent-fourier-cone
  - THM-3407-hadamard-core-multitoggle-response-plaquette-shells-and-trade-distance
related:
  - THM-3394-twelve-formerly-missing-hadamard-orders-through-2000
verified_by:
  - 04-computation/hadamard_pairwise_independent_toggle_filter_thm3411.py
  - 05-knowledge/results/hadamard_pairwise_independent_toggle_filter_thm3411.out
script: 04-computation/hadamard_pairwise_independent_toggle_filter_thm3411.py
output: 05-knowledge/results/hadamard_pairwise_independent_toggle_filter_thm3411.out
script_sha256: 376c8d5bf5a9f68672bd680ea5df367b5eb752fb30228b2154abdfc37e4d3415
output_sha256: 10fdb3d5cce9fd98e64b8ea78087ec4e66d06097a1c52af0bf518de0388ee93b
semantic_sha256: e702bae75e08cb5bba1e103a5eeb20f21ba540b3300edea9010f9f9c3bad9e02
hash_basis: LF-normalized bytes
---

# THM-3411 -- pairwise-independent toggle filter and sharp high-minor norm

## 1. Statement

Let `H` be a normalized real Hadamard matrix of order `4m`, with sign core
`K` and binary maxdet core `B=(J-K)/2` as in THM-3407.  Choose four distinct
core positions

```text
e_a=(i_a,j_a),                 delta_a=K_(i_a,j_a),       1<=a<=4.   (1)
```

Rows or columns may repeat.  For `S subseteq [4]`, in inherited event order,
put

```text
M_S=(product_(a in S) delta_a)
    det K[(i_a)_(a in S),(j_a)_(a in S)],       M_empty=1. (2)
```

Thus `M_S` is THM-3407's exact Boolean interaction.  For
`x in {+-1}^4`, toggle event `a` precisely when `x_a=-1`, and define

```text
f(x)=det B_(toggle(x))/det B.                                (3)
```

Let `nu` be any unbiased pairwise-independent law on four Rademacher bits.
Using THM-3396's coordinates, write

```text
a_r=E_nu product_(s!=r) X_s,          d=E_nu product_s X_s.   (4)
```

Let `u` denote the uniform law on the four-cube.  Define the five signed
high-minor response coordinates

```text
C_r=(4m M_([4]\{r})-M_[4])/(256m^4),       1<=r<=4,
C_5=M_[4]/(256m^4).                                      (5)
```

Then the exact filtered determinant response is

```text
E_nu f-E_u f=sum_(r=1)^4 a_r C_r+d C_5.                    (6)
```

Moreover the sharp range over all unbiased pairwise-independent laws is

```text
sup_nu |E_nu f-E_u f|
 =max{||C||_infinity,||C||_1/3}                            (7)

 =1/(256m^4) max{
      max_r |4m M_([4]\{r})-M_[4]|,
      [sum_r |4m M_([4]\{r})-M_[4]|+|M_[4]|]/3}.
```

Every extremum in `(7)` is attained by a finite strength-two orthogonal
array: an eight-run half-cube when the first term wins, or a twelve-run array
when the second term wins.  Thus `(7)` is an exact finite compiler, not only
a convex relaxation.

The four-event response is invisible to every pairwise-independent filter
if and only if

```text
M_[4]=0 and M_([4]\{r})=0 for all r.                       (8)
```

Equivalently, all four signed `3 by 3` event minors and the signed `4 by 4`
event minor vanish.  Pairwise plaquettes and Gram data do not decide `(8)`.

## 2. Boolean-to-Walsh proof

Give the four candidate toggles Boolean coordinates `z_a`.  THM-3407 proves
the complete multilinear identity

```text
F(z)=det(B+sum_a z_a delta_a e_(i_a)e_(j_a)^T)/det B
    =sum_(S subseteq [4])(-1/(2m))^|S| M_S product_(a in S)z_a. (9)
```

Set `z_a=(1-x_a)/2`.  Since

```text
product_(a in S) z_a
 =2^(-|S|) sum_(T subseteq S)(-1)^|T| product_(a in T)x_a,   (10)
```

the Walsh coefficient of `f(x)=F((1-x)/2)` at `T` is

```text
fhat(T)=(-1)^|T| sum_(S superseteq T)
                  (-1/(2m))^|S| M_S/2^|S|.                 (11)
```

For `T=[4]\{r}` and `T=[4]`, equation `(11)` gives respectively

```text
fhat([4]\{r})=M_([4]\{r})/(64m^3)-M_[4]/(256m^4)=C_r,
fhat([4])=M_[4]/(256m^4)=C_5.                              (12)
```

The law `nu` and the uniform law have identical Walsh moments in degrees
zero, one, and two.  Their only possible differences are exactly the four
cubic moments and the quartic moment `(4)`.  Fourier pairing with `(12)`
proves `(6)`.

THM-3396 proves that the feasible packet `(a_1,a_2,a_3,a_4,d)` is the polar
of the odd five-demicube and that its exact absolute functional norm is
`max(l-infinity,l-one/3)`.  It also identifies all twenty-six vertices as
ten eight-run and sixteen twelve-run orthogonal-array laws.  Applying that
norm to `(6)` proves `(7)` and the equality statement.  Finally `C_5=0`
forces `M_[4]=0`, after which `C_r=0` is equivalent to
`M_([4]\{r})=0`; this proves `(8)`.  QED.

## 3. Integer OA form

For an `OA(N,4,2,2)` let `(A_1,A_2,A_3,A_4,D)` be its unnormalized cubic
and quartic moments.  Averaging the sixteen toggle responses over its rows
gives the exact integer identity

```text
(1/N)sum_(rows x) f(x)-E_u f
       =[sum_r A_r C_r+D C_5]/N.                            (13)
```

There is no hidden row matching: THM-3396 reconstructs the unique sixteen
cell counts from `(N,A,D)`.  However `(13)` is only an average over four
candidate toggles.  It does not complete the four OA columns to the core
`K`, recover the individual toggle responses, or preserve Hadamard
equivalence.

## 4. Universal size bound and exact order-eight boundary

A `3 by 3` sign determinant belongs to `{0,+-4}`, while a `4 by 4` sign
determinant has absolute value at most `16`.  Hence `(7)` gives the universal
but not asserted sharp bound

```text
|E_nu f-E_u f| <= (4m+5)/(48m^4).                            (14)
```

For the Paley order-eight core (`m=2`), exhaustive enumeration of all
`29,400` nonattacking four-event packets gives `504` blind packets and the
following exact sharp-bias histogram:

| sharp bias | packets |
|---:|---:|
| `0` | `504` |
| `5/1536` | `672` |
| `3/512` | `8064` |
| `5/768` | `2016` |
| `11/1536` | `5376` |
| `1/128` | `8064` |
| `13/1536` | `2016` |
| `5/512` | `2688` |

This census is finite evidence for the response palette, not part of the
universal proof.

Four zero-based order-eight core controls make the sidecar sharp:

```text
events ((0,1),(1,0),(2,3),(3,2)):
  (M_(-1),M_(-2),M_(-3),M_(-4);M_4)=(0,0,0,0;8),

events ((0,0),(1,1),(2,2),(3,4)):
  (0,4,0,4;0),

events ((0,1),(1,2),(2,5),(3,0)):
  (0,4,0,4;8),

events ((0,0),(1,1),(2,3),(5,6)):
  (0,0,0,0;0).                                               (15)
```

Here `M_(-r)=M_([4]\{r})`.  The first two show that neither the four-minor
nor the four three-minors can be deleted.  The middle pair has identical
three-minor data and different four-minor, so the failure is not merely a
support issue.  The last packet is the exact blind control required by `(8)`.

## 5. Connection and loss ledger

| field | content |
|---|---|
| source | four candidate toggles in a normalized Hadamard core |
| target | the expectation of their determinant response under a whitened four-bit law |
| map | Boolean toggle polynomial `(9)` followed by `z=(1-x)/2` and Walsh pairing |
| preserved | normalized determinant response, all cubic/quartic interactions, OA convex mixtures, equality laws |
| destroyed | the individual row word, lower response localization, matrix equivalence, and global OA completion |
| exact sidecar | the four signed `3 by 3` minors and one signed `4 by 4` minor |

The bridge is useful precisely when degree-at-most-two toggle data have been
whitened.  It is counterindicated for biased or pair-correlated columns:
then lower Walsh coefficients survive and the five-coordinate norm is not a
complete response.  No Grothendieck, Crouzeix, LRC(14), JC(2), or Hadamard-
existence conclusion follows.

## 6. Reproduction and status

Run

```text
python3 04-computation/hadamard_pairwise_independent_toggle_filter_thm3411.py
python3 -O 04-computation/hadamard_pairwise_independent_toggle_filter_thm3411.py
```

The standard-library companion performs independent Boolean and Walsh paths,
reconstructs all twenty-six OA vertex laws, checks direct determinants on the
four hostiles, exhausts the order-eight nonattacking universe, and audits an
order-twelve control window.  Normal and optimized runs are byte-identical to
the frozen output.  The LF-normalized script/output hashes are respectively

```text
376c8d5bf5a9f68672bd680ea5df367b5eb752fb30228b2154abdfc37e4d3415
10fdb3d5cce9fd98e64b8ea78087ec4e66d06097a1c52af0bf518de0388ee93b
```

and the semantic digest is
`e702bae75e08cb5bba1e103a5eeb20f21ba540b3300edea9010f9f9c3bad9e02`.
The combined independent immutable-file proof/replay audit reconstructed the
Walsh signs, demicube norm, blindness converse, and sharp high-sidecar hostile,
and found no defect.  Normal, optimized, and stored outputs and every pinned
hash agree.
