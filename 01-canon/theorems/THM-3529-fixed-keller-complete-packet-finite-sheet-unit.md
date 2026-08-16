---
id: THM-3529
title: "Fixed Keller complete-packet finite-sheet unit obstruction"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT HOSTILE AUDIT.
  The candidate identifies the regular finite inverse divisor with the
  explicit irreducible beta-homogeneous pullback B=F^*L and uses the
  x-free complete minimum-beta packet face to exclude B from every complete
  packet.  No statement in this file is proved canon until audit promotion.
source: codex/finite-branch-beta-obstruction/2026-08-16
depends_on: []
related:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary
  - THM-3528-fixed-keller-all-level-cleared-norm-polynomiality-and-finite-sheet-defect
script: 04-computation/keller_finite_branch_beta_homogeneous_unit_obstruction_20260816.py
output: 05-knowledge/results/keller_finite_branch_beta_homogeneous_unit_obstruction_20260816.out
script_sha256: 0e2e7ba0ee5aa45ca83aa050255cdb2f8ffda8d05c0495b7fc1c513ad1c289dc
output_sha256: 15b026a541afc09d6be0784f00e2c05e3c3b3b148625a1b30dec7c62875d89e6
semantic_sha256: 7c19ce6a3dc08f36fdb19487a7e26411077e1db880afb0201ddf8a67c3df25d1
hash_basis: LF-normalized bytes
---

# THM-3529 -- the finite inverse divisor misses every complete packet

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT HOSTILE
AUDIT.**  The argument below is deliberately not in the proved dependency
graph until its geometric typing and graded-domain step are independently
accepted.

Retain the fixed sporadic Keller map of THM-2473, written in source
coordinates `(x,y,z)`, with `u=1+xy`,

```text
F=(u^3z+y^2u(4+3xy),
   y+3xu^2z+3xy^2(4+3xy),
   2x-3x^2y-x^3z),                    det J_F=-2,       (1)
```

and the irreducible target boundary

```text
L(a,b,c)=27a^2c^2-18abc+16a+b^3c-b^2.                 (2)
```

For a source monomial put

```text
beta(x^i y^j z^k)=i-j-2k.                              (3)
```

## 1. Candidate theorem

Let `P` be any nonzero polynomial with the complete packet `A(e,m)` of
THM-3506.  Let `s_L(P)` be the regular finite-sheet valuation in THM-3528.
Then the candidate conclusion is

```text
s_L(P)=0.                                               (4)
```

Consequently THM-3528's cleared norm satisfies

```text
Q(P)=L^e N(P) in Q[a,b,c],
ord_L(Q(P))=0,
gcd(Q(P),L)=1,                                          (5)
```

and has complete packet `A(7e-2m,3e-2m)`.  For the raw tower

```text
P_0=L,                 P_(n+1)=L^e_n N(P_n),            (6)
```

this would prove

```text
gcd(P_n,L)=1 for every n>=1.                            (7)
```

Notice that `(7)` starts at `n=1`: the seed `P_0=L` is, of course, the old
boundary itself.

## 2. The finite sheet has the explicit source equation `B=F^*L`

Direct exact substitution gives the primitive ten-term polynomial

```text
B:=L(F(x,y,z))
 =-9x^4y^2z^2-54x^3y^3z-18x^3yz^2-81x^2y^4
  -72x^2y^2z-9x^2z^2-54xy^3+6xyz+63y^2+16z.           (8)
```

Every monomial in `(8)` has beta weight `-2`, while

```text
deg_x B=4,                 B(0,y,z)=63y^2+16z.          (9)
```

The polynomial `B` is irreducible over `Q`.  This has a short certificate,
not merely a computer-factorization witness.  Localize at `x` and put

```text
p=xy,                     q=x^2z.                       (10)
```

Then `Q[x^+-1,y,z]=Q[x^+-1,p,q]` and

```text
x^2B=C(p,q)
 =-9(p+1)^2q^2
  -2(27p^3+36p^2-3p-8)q
  -9p^2(9p^2+6p-7).                                    (11)
```

The three coefficients in `(11)` are primitive in `Q[p]`, and

```text
disc_q C=64(3p+4).                                     (12)
```

The linear polynomial `3p+4` has odd valuation at `p=-4/3`, hence is not a
square in `Q(x,p)`.  The quadratic `C` is therefore irreducible over that
field and, by Gauss, in the Laurent polynomial ring.  There `B` is associated
to `C`.  A factor of `B` that became a Laurent unit would be a power of `x`,
but `(9)` shows that `x` does not divide `B`.  Thus `B` is irreducible already
in `Q[x,y,z]`.

Let `C_fin` be THM-3528's closure of the unique regular finite inverse branch
over the generic point of `D=V(L)`.  The branch dominates `D`, so `C_fin` is
an irreducible source divisor.  Equation `(8)` vanishes identically on it,
giving

```text
C_fin subseteq V(B).                                   (13)
```

Both sides of `(13)` are irreducible of dimension two; hence

```text
C_fin=V(B).                                             (14)
```

The exact witness

```text
(2,5/6,-7/8) |-> (2/27,1,1) in V(L),
grad B(2,5/6,-7/8)=(-29/4,33,24)                       (15)
```

checks the regular sheet at the canonical DVR point.  More importantly,
THM-3528 already proves the divisor identity, so `(14)` gives

```text
s_L(P)=ord_(C_fin)(P)=ord_B(P).                         (16)
```

The equation in `(16)` is an order along one prime divisor.  It is not a
count of affine points, reduced components, or inverse branches.

## 3. The minimum-beta face forbids the finite divisor

Suppose for contradiction that `s_L(P)>0`.  By `(16)` and primeness of `B`,
there is a nonzero polynomial `R` with

```text
P=BR.                                                   (17)
```

Write `in_beta,min` for the complete part of least beta weight.  Since all
of `B` has the single weight `-2`, decomposition into beta-homogeneous pieces
in the polynomial domain gives

```text
in_beta,min(P)=B in_beta,min(R).                        (18)
```

There is no cancellation loophole in `(18)`: its right factor is the unique
least-weight part of `R`, and the associated graded ring for the integral
weight `(1,-1,-2)` is again a domain.

Packet completeness is now decisive.  The prescribed nonzero face is, up to
a rational unit,

```text
in_beta,min(P)=
 y^(3e-2m) z^(e-m)
 (y^2+27z)^(2m/3)(y^2+108z)^(m/3).                     (19)
```

It has beta weight `-5e+2m` and belongs to `Q[y,z]`; in particular its
`x`-degree is zero.  But degree in `x` is additive in the domain
`Q[y,z][x]`, so `(9)` and `(18)` imply

```text
deg_x in_beta,min(P)
 =4+deg_x in_beta,min(R) >=4,                           (20)
```

contradicting `(19)`.  Therefore `(17)` is impossible and `(4)` follows.
Equations `(5)`--`(7)` then follow directly from THM-3528.

## 4. What the obstruction says structurally

The complete packet is more than a five-number growth certificate.  Its
minimum-beta face is a **transverse sidecar** that detects the unique affine
return divisor.  The finite branch is beta-homogeneous but carries positive
`x`-degree; the packet face is beta-homogeneous but entirely `x`-free.  These
two facts cannot coexist in a product.  Thus the finite-sheet defect is not
an uncontrolled scalar after all: it is excluded by a coordinate discarded
by the earlier valuation-only recurrence.

The torus quotient `(10)` adds a second interpretation.  The same map factor
`4+3xy=4+3p` occurs as the nonsquare seam in `(12)`.  This is an exact
algebraic recurrence of the fixed map's construction, not evidence for an
unrelated LRC prime or tournament current.

## 5. Equality boundary and hostile controls

Each hypothesis is load-bearing.

1. **Completeness.**  The divisor `B` itself is beta-homogeneous and is
   divisible by `B`, but its beta face has positive `x`-degree.  It is not a
   complete packet of the form `(19)`.
2. **Positive `x`-degree.**  An `x`-free homogeneous divisor such as `y^2`
   can divide an `x`-free face.  Weight homogeneity alone is insufficient.
3. **The fixed source divisor.**  The argument uses the exact identity
   `C_fin=V(F^*L)` and does not transfer automatically to another Keller map.
4. **Minimum-face equality.**  Merely exhibiting one monomial from `(19)`
   would not work; the complete-face assertion excludes additional
   equal-weight `x`-terms.

Even if promoted, this theorem would prove only old-`L` coprimality of the
cleared outputs.  It would not prove that any later `P_n` is irreducible,
that `V(P_n)` is a new image prime, that a later eliminant is separable, that
the discriminant recursion continues, that all nonproperness components are
distinct, or that the construction says anything about the general Jacobian
conjecture.

## Reproduction

```text
python -B 04-computation/keller_finite_branch_beta_homogeneous_unit_obstruction_20260816.py
python -B -O 04-computation/keller_finite_branch_beta_homogeneous_unit_obstruction_20260816.py
```

Normal and optimized transcripts match the stored output exactly.  The
companion verifies `(1)`, `(8)`--`(12)`, `(15)`, all 7,700 admissible packet
exponent rows with `e<=300`, and both hostile boundaries.
