---
id: THM-4377
title: "Reciprocal source mutation and boundary-jet cokernel"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4368/4375 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED; JC(2) OPEN.
  Reciprocal reflection of a bilateral diagonal-packet type is induced on
  its source exponent lattice by one explicit Laurent mutation. Its failure
  to preserve polynomial source monomials is exactly the union of three
  named coordinate-wall defects. At boundary scale w>=ceil(ell/2), the
  reflected fibre embeds integrally into the original fibre and the cokernel
  is a finite u-boundary jet; from w>=ceil(ell/2)+d-1 it has the sharp rank
  2d and degrees 0,...,2d-1. The result identifies a missing source sidecar,
  not a bracket, chart-entry, seam, JC(2), or DC(2) obstruction.
source: root + reciprocal_kernel clean-room referee / JC2 continuation session, 2026-09-03
depends_on:
  - THM-4368-diagonal-boundary-valuation-triangular-address-and-simplex-stream-rank
  - THM-4375-bilateral-source-cone-reciprocal-orbits-and-fibre-imbalance
related:
  - THM-4369-source-packet-pascal-circuit-kernel-and-boundary-basis
  - THM-4378-bilateral-packet-quotient-reciprocal-eigenlattice
script: 04-computation/jc2_reciprocal_source_mutation_boundary_jet_thm4377.py
output: 05-knowledge/results/jc2_reciprocal_source_mutation_boundary_jet_thm4377.out
script_sha256: 3e921ddf147dd0b6c82d3e72c4c4188039ab2a0d437e347d73d3863356c9beb4
output_sha256: f6df88dbe69231323aad543d3d794522eb7690f3f75d9b4d61f80f6ec940b872
hash_basis: raw LF bytes
audit: >
  PASS. Direct nonnegative-exponent enumeration is compared with the closed
  source intervals on 332,295 bilateral fibres for ell=2..90,
  ceil(ell/3)<=w<=180, and every 0<=d<ceil(ell/2). The companion performs
  13,370,117 explicit checks of the mutation, inverse image, disjoint
  three-wall decomposition, exact defect counts, matched-core clock,
  boundary-jet stabilization, graded identity, and named hostile controls.
  A separate clean-room derivation recovered the four interval endpoints,
  three disjoint wall counts, sharp full-matching clock, split module
  sequence, and full 2d jet threshold without using the companion.
  Normal, optimized, and nondefault-hash-seed replays reproduce the frozen
  output byte-for-byte.
---

# THM-4377 -- Reciprocal source mutation and boundary-jet cokernel

**PROVED ELEMENTARY RELATIVE TO THM-4368/4375 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. JC(2) REMAINS OPEN.**

## 1. The two source fibres

Fix `ell>=2` and put

```text
s=ceil(ell/2),                     rho=ceil(ell/3).       (1)
```

For a source-realizable THM-4368 boundary type `(u,v)`, its individual
source monomials are indexed by

```text
m_e(u,v)=(a,b,c,e)
        =(2v+e-ell, s+u-v-1-e, v-e, e),                 (2)

max(0,ell-2v)<=e<=min(v,s+u-v-1).                       (3)
```

Write `E_ell(u,v)` for this finite set of exponent tuples. By THM-4375,
every bilateral reciprocal orbit has a unique representative

```text
(u,v)=(w+d,w),          w>=rho,          0<=d<=s-1.      (4)
```

Put

```text
E_+=E_ell(w+d,w),                 E_-=E_ell(w,w+d).      (5)
```

The closest proved mechanism is THM-4375's exact pair of fibre-cardinality
formulas. Its first hostile is already asymmetric: at `ell=10`, the orbit
`(5,4)<->(4,5)` has sizes `3<->4`. The corrected near miss is therefore to
infer a source involution from the ambient trace reflection. The least-used
sidecar is the actual source exponent tuple, rather than its packet type or
fibre size. This theorem retains that tuple and locates the obstruction.

## 2. One Laurent mutation swaps the boundary pair

On the ambient exponent lattice define

```text
T_d(a,b,c,e)=(a+2d,b-2d,c+d,e).                         (6)
```

Direct substitution in `(2)` gives the exact identity

```text
T_d(m_e(w+d,w))=m_e(w,w+d).                             (7)
```

Thus reciprocal reflection is the exponent mutation associated formally to
the Laurent monomial ratio

```text
(x^2 p/u^2)^d.                                          (8)
```

This is an identity in the exponent lattice. It is not a polynomial-algebra
endomorphism: the new `u` exponent may be negative. That distinction is the
whole obstruction below.

Define the matched pieces

```text
E_+^0={m in E_+: b(m)>=2d},
E_-^0={m in E_-: a(m)>=2d and c(m)>=d}.                 (9)
```

Then

```text
T_d:E_+^0 -> E_-^0                                      (10)
```

is a bijection, with inverse `T_(-d)`. Indeed, starting in `E_+`, the
coordinates `a+2d,c+d,e` remain nonnegative, so `T_d(m)` is polynomial
exactly when `b>=2d`. Starting in `E_-`, the inverse can fail polynomiality
only through `a-2d<0` or `c-d<0`, which gives `(9)`.

## 3. The exact three-wall decomposition

The unmatched polynomial tuples are

```text
B_+={m in E_+:b(m)<2d},
A_-={m in E_-:a(m)<2d},
C_-={m in E_-:c(m)<d}.                                 (11)
```

They give disjoint decompositions

```text
E_+=E_+^0 disjoint_union B_+,
E_-=E_-^0 disjoint_union A_- disjoint_union C_-.        (12)
```

Only the disjointness of `A_-` and `C_-` needs comment. If one target tuple
lay in both, its unchanged exponent `e` would satisfy

```text
e<ell-2w                         and e>w.                (13)
```

But `w>=rho=ceil(ell/3)` gives `ell-2w<=w`, a contradiction. The `B_+`
and `C_-` defect species also cannot both be nonempty: `B_+!=empty` requires
`w>=s-d`, whereas `C_-!=empty` requires `w<=s-d-2`.

The decomposition has exact closed counts. To avoid confusing interval
endpoints with the wall sets in `(11)`, put

```text
L_+=max(0,ell-2w),             U_+=min(w,s+d-1),
L_-=max(0,ell-2w-2d),          U_-=min(w+d,s-d-1),       (14)
```

and `[x]_+=max(x,0)`. Then

```text
nu    =|E_+^0|=|E_-^0|
      =[min(w,s-d-1)-L_++1]_+,

beta  =|B_+|=[U_+-max(L_+,s-d)+1]_+,
alpha =|A_-|=[min(U_-,L_+-1)-L_-+1]_+,
gamma =|C_-|=[U_--max(L_-,w+1)+1]_+.                  (15)
```

Consequently

```text
|E_+|=nu+beta,                    |E_-|=nu+alpha+gamma. (16)
```

This refines THM-4375's cardinality imbalance into its actual coordinate-wall
anatomy.

There is also an exact graded version. Let `r` be a formal variable and set

```text
H_+(r)=sum_(m in E_+) r^e(m),
H_-(r)=sum_(m in E_-) r^e(m).                           (17)
```

Since `(7)` preserves `e`, the matched terms cancel and

```text
H_+(r)-H_-(r)
 =sum_(m in B_+)r^e(m)-sum_(m in A_- union C_-)r^e(m).  (18)
```

Thus the three defects are complete for every statistic depending only on
the retained `e` grade, not merely for total fibre size.

## 4. The sharp matching clock

The matched core satisfies

```text
nu<=s-d.                                                (19)
```

For a nonfixed orbit `d>0`, equality occurs exactly when

```text
w>=s.                                                   (20)
```

Indeed, the common `e` interval is

```text
[L_+,min(w,s-d-1)].                                     (21)
```

It contains all `s-d` integers `0,...,s-d-1` exactly when `L_+=0`, which is
equivalent to `w>=ceil(ell/2)=s`; that inequality also settles the upper
endpoint. If `w<s`, then `L_+>0` and `L_-<L_+`; the target tuple with
`e=L_-` lies in the wall set `A_-`. Hence `(20)` is the sharp full-matching
threshold for this Laurent mutation, not just a sufficient large-scale
estimate. This does not exclude an unrelated source operation. For `d=0`,
`T_0` is the identity and `E_+=E_-` at every allowed `w`.

For `d>0` and `w>=s`, formulas `(14)--(15)` simplify to

```text
E_-^0=E_-={m_e(w,w+d):0<=e<=s-d-1},
A_-=C_-=empty,
beta=min(w-s+d+1,2d),
{b(m):m in B_+}={max(0,s+d-1-w),...,2d-1}.             (22)
```

Thus the inverse mutation gives an integral embedding of free source modules

```text
T_(-d): Z[E_-] -> Z[E_+],                               (23)
```

and the monomial bases split the exact sequence

```text
0 -> Z[E_-] -> Z[E_+] -> Z[B_+] -> 0.                  (24)
```

The cokernel is precisely the finite `u`-boundary jet, not an unspecified
failure of reciprocity.

The jet reaches its full size exactly at

```text
w>=s+d-1.                                               (25)
```

At and beyond that scale,

```text
|B_+|=2d,             {b(m):m in B_+}={0,1,...,2d-1},  (26)

H_+(r)-H_-(r)
 =r^(s-d)+r^(s-d+1)+...+r^(s+d-1)
 =r^(s-d)(1-r^(2d))/(1-r).                             (27)
```

Both thresholds are sharp. Equation `(20)` cannot be lowered because an
`a`-wall target defect exists at every `w<s`. In `(26)`, the bound
`|B_+|<=2d` is attained only when the interval endpoint
`U_+=s+d-1`, equivalently
`w>=s+d-1`.

## 5. Exact boundary controls

At `ell=3`, `(s,rho)=(2,1)` and the first nonfixed bilateral orbit has

```text
(u,v)=(2,1)<->(1,2),
E_+={(0,1,0,1)},             E_-={(1,0,2,0)}.           (28)
```

The fibres have equal size but no matched tuple. This is the smallest
hostile to replacing cardinality equality by a polynomial source involution.

At the first THM-4375 asymmetric hostile, `ell=10,w=4,d=1`, the `e` sets are

```text
E_+: {2,3,4},                 E_-:{0,1,2,3},
matched:{2,3},                B_+:{4},
A_-:{0,1},                    C_-:empty.                (29)
```

In particular,

```text
H_+(r)-H_-(r)=r^4-r-1.                                 (30)
```

Both target walls can occur simultaneously. At `ell=15,w=5,d=1`,

```text
E_+:{5},       E_-:{3,4,5,6},       A_-:{3,4}, C_-:{6}. (31)
```

Finally, the sharp-asymmetry orbit of THM-4375 at `ell=10,w=8,d=4` is
already in the stable range:

```text
E_+:{0,1,...,8},              E_-:{0},
matched:{0},                  {b(m):m in B_+}={0,...,7}. (32)
```

Its `8=2d` unmatched monomials are exactly the full boundary jet predicted
by `(26)`.

## 6. Connection contract and scope

```text
source:      individual polynomial source monomials in a bilateral packet
target:      reciprocal exponent-lattice fibre and its boundary cokernel
map:         T_d(a,b,c,e)=(a+2d,b-2d,c+d,e)
preserved:   diagonal intercept ell, reciprocal boundary pair, e grade
destroyed:   nonnegativity after the Laurent mutation
sidecar:     the three coordinate walls a<2d, b<2d, c<d; eventually the
             finite u-boundary jet b=0,...,2d-1
hostiles:    (28) has equal cardinalities and no lift; (29) has mixed signs
```

The mutation is an exact source-semigroup calculation. It does not prove
that a bracket functional factors through the boundary jet, that a
hypothetical counterexample enters this packet chart, or that reciprocal
reflection is a symmetry of an arbitrary source polynomial. In particular,
the jet is a candidate sidecar for the bracket-blind hierarchy, not a proved
bracket obstruction. No seam, `JC(2)`, or `DC(2)` conclusion follows.

## 7. Exact audit

The companion reconstructs each fibre by scanning all allowed nonnegative
`e` exponents and independently compares it with `(3)`. On the declared
universe it checks `(6)--(18)`, the sharp equivalence `(20)`, the split exact
sequence `(24)`, both stabilization clocks, the graded defects, and all four
controls `(28)--(32)`. It uses explicit runtime checks, so optimized Python
does not remove the audit.

```text
python3 -B 04-computation/jc2_reciprocal_source_mutation_boundary_jet_thm4377.py
python3 -B -O 04-computation/jc2_reciprocal_source_mutation_boundary_jet_thm4377.py
PYTHONHASHSEED=271828 python3 -B \
  04-computation/jc2_reciprocal_source_mutation_boundary_jet_thm4377.py
```

All three commands reproduce the frozen output byte-for-byte. **QED.**
