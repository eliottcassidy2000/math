---
id: THM-3876
title: "Monomial marked-root profiles descend exactly at reduced exponent at most two"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For positive exponents
  m,n, reduce by g=gcd(m,n) to coprime M=m/g,N=n/g.  The cubic cusp identity
  forces B=2r^(2N)(3+4r^(M+N)) on the normalization.  This value descends to
  k[r^M,r^N+r^(M+2N)] if and only if M<=2.  For M>=3 an explicit
  primitive-root pair has the same (A,C) and different B.  The already
  audited N=1 specialization is the first boundary of the uniform theorem.
source: jc_sparse_direct_search / post-THM-3873 monomial normalization tower, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit first proved
  the N=1 theorem and then rederived the gcd reduction, finiteness and exact
  normalization, the forced global sign and profile, both constructive
  M=1,2 descents, and the primitive-root collision for every coprime M,N
  with M>=3.  It checked that eta=zeta^N retains exact order M and that the
  branch-ring address condition is used only in the necessary direction.
  The expanded assertion-free companion checks the
  universal cusp and forced-profile identities, both M=1,2 descent
  boundaries for N<=12, and 102 coprime primitive-root hostiles with
  3<=M<=16 and 1<=N<=12.  Normal and optimized runs byte-match the frozen
  361-gate transcript.  The proof for arbitrary coprime M,N uses the exact
  order of eta=zeta^N and does not depend on the finite grid.
related:
  - THM-3866-all-polynomial-graph-branches-force-projective-companion
  - THM-3873-first-nongraph-triangular-parabola-companion
script: 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
output: 05-knowledge/results/jc2_monomial_marked_root_nondescent_thm3876.out
script_sha256: d5d8b4aea098ee0542f1b92db72fefec247fb7635a34fec732a42b2b0973d5d3
output_sha256: f6aa41443f728bddbc4c916adfbd00f04ec8db0365f172d7b3e478670411a894
semantic_sha256: 41b2b480c79ba8182ec5c9531e31eb2d27187cd46eea3dd8e037d8814ae45ffa
hash_basis: raw LF bytes
---

# THM-3876 -- the complete two-exponent monomial descent boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an algebraically closed field `k` of
characteristic zero.  Put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b.                       (1)
```

Start with any positive integers `m,n` and the monomial marked-root data

```text
a=t^m,                         s=t^n,
A=a,                           C=6s(1+as).                     (2)
```

Reduce the exponents by

```text
g=gcd(m,n),             r=t^g,             M=m/g, N=n/g,
gcd(M,N)=1.                                                     (3)
```

The image curve `Gamma_(M,N)` has normalization

```text
nu_(M,N): A1_r -> A2_(A,C),
A=r^M,                  C=6r^N(1+r^(M+N)).                     (4)
```

If `b in k[A,C]` and `Delta_b` vanishes on `Gamma_(M,N)`, then its pullback
to the normalization is uniquely forced:

```text
b(A(r),C(r))=B_(M,N)(r)
             :=2r^(2N)(3+4r^(M+N)).                            (5)
```

The exact descent classification is

```text
B_(M,N) in k[r^M,r^N+r^(M+2N)]       iff       M<=2.           (6)
```

More explicitly:

```text
M=1: A=r,
     B=2A^(2N)(3+4A^(N+1));                            descent,

M=2: N is odd,
     S=C/6-A^(N+1)=r^N,
     B=2S^2(3+4AS);                                    descent,

M>=3: B_(M,N) does not descend.                       nonentry. (7)
```

Consequently, if the reduced `A`-exponent `M` is at least three, there is
**no** polynomial carrier profile `b(A,C)` whose packet `V(Delta_b)`
contains `Gamma_(M,N)`.  This is an obstruction before transverse
factorization or Jacobian-mate analysis: the cusp identity forces a
polynomial function of the normalization parameter which the singular
plane-curve coordinate ring cannot remember.

The theorem covers the complete two-exponent monomial marked-root family
`(2)`.  It does not close arbitrary polynomial pairs `a(r),s(r)`, arbitrary
singular plane curves with affine-line normalization, or non-graph branches
outside this monomial tower.

## 1. Gcd reduction and the actual normalization

If `g>1`, the original `t`-parametrization in `(2)` is not birational: it
factors through `r=t^g`.  The reduced image ring is

```text
R_(M,N)=k[r^M,r^N+r^(M+2N)] subset k[r].                       (8)
```

The factor `6` in `C` is harmless.  The element `r` is integral over
`R_(M,N)`, since it is a root of the monic polynomial
`X^M-A in R_(M,N)[X]`.  Hence `k[r]` is finite over `R_(M,N)`.

We must still prove birationality; an affine-line parametrization need not
be an embedding.  Away from `r=0`, two points with the same `A`-coordinate
have the form

```text
r'=zeta r,                         zeta^M=1.                    (9)
```

Put

```text
eta=zeta^N.                                                        (10)
```

Because `gcd(M,N)=1`, the map `zeta -> zeta^N` permutes the `M`th roots of
unity; in particular `eta=1` if and only if `zeta=1`.  Equality of the two
`C`-coordinates is equivalent to

```text
(eta-1)+(eta^2-1)r^(M+N)=0.                                    (11)
```

For `zeta=1` this is the diagonal.  If `eta=-1`, equation `(11)` is
impossible.  Every remaining nontrivial root gives at most one value of
`r^(M+N)`.  Therefore all off-diagonal coincidences are finite and the
generic fibre of the finite map `(4)` has one point.  Thus

```text
Frac(R_(M,N))=k(r).                                             (12)
```

The ring `k[r]` is exactly the integral closure of `R_(M,N)` in `(12)`.
Indeed it is integral over `R_(M,N)` and integrally closed; conversely, an
element of `k(r)` integral over `R_(M,N)` is also integral over `k[r]`, so
it belongs to `k[r]`.

This step is why the reduced pair `(M,N)`, rather than the unreduced
exponents `(m,n)`, is the invariant in `(6)`.

## 2. The cusp identity forces the marked value

The universal identity behind `(1)` is

```text
A^2 Delta_b=27(P^3-u^2),
P=1+(2/3)AC,                    u=1+AC+A^2b.                    (13)
```

Set

```text
U=r^(M+N).                                                      (14)
```

Along `(4)`, one has

```text
AC=6U(1+U),
P=1+4U(1+U)=(1+2U)^2.                                         (15)
```

If `Delta_b` vanishes on `Gamma_(M,N)`, then in the domain `k[r]`

```text
u^2=(1+2U)^6.
```

Hence `u=(1+2U)^3` or `u=-(1+2U)^3`.  At `r=0`, however, `A=C=0`
and `u(0)=1`, so characteristic zero excludes the negative sign.  Writing
`B(r)=b(A(r),C(r))` and solving `(13)` gives

```text
r^(2M)B(r)
 =(1+2U)^3-1-6U(1+U)
 =2U^2(3+4U),                                                  (16)

B(r)=2r^(2N)(3+4r^(M+N)).                                     (17)
```

Thus `(5)` is necessary for every polynomial carrier.  Conversely, any
polynomial `b(A,C)` having restriction `(17)` makes `(13)` vanish on the
curve, so the existence question is exactly the branch-ring membership in
`(6)`.

## 3. The two and only two descent regimes

### 3.1 Reduced exponent M=1

Here `A=r`, so the entire normalization ring is already generated by `A`.
Both marked coordinates are polynomial in `A`:

```text
C=6A^N(1+A^(N+1)),
B_(1,N)=2A^(2N)(3+4A^(N+1)).                                  (18)
```

Thus `B_(1,N)` descends.  These are forward graph branches.

### 3.2 Reduced exponent M=2

Coprimality forces `N` to be odd.  Define the polynomial sidecar

```text
S=C/6-A^(N+1).                                                 (19)
```

On the normalization,

```text
S=r^N+r^(2+2N)-r^(2N+2)=r^N.                                  (20)
```

Therefore

```text
B_(2,N)=2S^2(3+4AS) in k[A,C].                                (21)
```

This is the only nontrivial descent boundary.  At `N=1`, `(19)` is the
triangular coordinate `S=C/6-A^2=r` and `(21)` is exactly THM-3873.  For
larger odd `N`, the sidecar recovers `r^N`, not necessarily `r`; descent of
the marked value does not assert that the image curve is smoothly embedded.

## 4. Uniform address separation for M>=3

Assume `M>=3`.  Choose a primitive `M`th root of unity `zeta` and put

```text
eta=zeta^N.                                                     (22)
```

Since `gcd(M,N)=1`, `eta` is again primitive of exact order `M`.  Hence

```text
eta !=1,                         eta !=-1,              eta+1 !=0.  (23)
```

Choose `r in k^*` satisfying

```text
r^(M+N)=-1/(eta+1).                                           (24)
```

The algebraic closure and the nonzero right side guarantee such an `r`.
The addresses `r` and `zeta r` are distinct and have the same
`A`-coordinate.  Their `C/6` difference, divided by `r^N`, is

```text
(eta-1)+(eta^2-1)r^(M+N)=0                                    (25)
```

by `(24)`.  Thus

```text
nu_(M,N)(r)=nu_(M,N)(zeta r).                                 (26)
```

The forced marked values are nevertheless different.  Since
`zeta^(M+N)=eta`, direct substitution gives

```text
B_(M,N)(zeta r)-B_(M,N)(r)
 =2r^(2N)[3(eta^2-1)+4(eta^3-1)r^(M+N)]
 =-2r^(2N)(eta-1)^3/(eta+1) !=0.                              (27)
```

Every factor in the last expression is nonzero by `(23)--(24)`.

## 5. Exact non-descent and the iff

Every polynomial `b(A,C)` pulls back to an element of `R_(M,N)`, so it must
take the same value at any two normalization addresses with the same
`(A,C)`.  Equations `(26)--(27)` show that `B_(M,N)` violates this necessary
condition whenever `M>=3`.  Hence

```text
B_(M,N) notin R_(M,N)                         for M>=3.         (28)
```

Sections 3.1--3.2 give explicit branch-ring representatives for every
`M<=2`, while `(28)` excludes every `M>=3`.  This proves the iff `(6)` and,
through the forced-value calculation in Section 2, the polynomial-carrier
nonentry.

No reducedness or squarefreeness of the total discriminant packet is used.
Repeated components and arbitrary transverse terms cannot repair failure
of the restriction itself to descend to the reduced image curve.

## 6. Exact replay and structural boundary

Run

```text
python3 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
python3 -O 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
```

and compare both streams byte-for-byte with the frozen output.  The exact
companion checks the universal identities `(13)--(17)` and `(25)--(27)`,
the `M=1,2` formulas for `N<=12`, and 102 coprime cyclotomic pairs with
`3<=M<=16`, `1<=N<=12`.  The grid is a hostile control, not the logical
scope of the theorem: exact order preservation under
`zeta -> zeta^N` proves `(23)` for every coprime pair.

The sharpened mechanism is **gcd-reduced normalization-address
separation**.  Monomial degree by itself is not the invariant.  The first
reduced exponent determines how much of the normalization survives in the
plane branch ring:

```text
M=1     the base coordinate survives;
M=2     one polynomial sidecar recovers the marked monomial r^N;
M>=3    a primitive-root collision destroys marked-value descent. (29)
```

This ends the whole two-exponent monomial construction grammar.  The next
genuinely open lane requires nonmonomial `a(r),s(r)` or a singular
affine-line normalization whose coincident-address fibres are separated by
the forced marked value in a subtler way.
