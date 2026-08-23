---
id: THM-3884
title: "Cusp residual total-degree and leading-gauge filtration"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.  In
  the complete THM-3881 residual equation, every square survivor with
  nonconstant f satisfies deg T >= deg f+1.  At equality, the highest form is
  243*x^3*f_n^2*((y^2-15x^2)f_n+xT_m)^2; its odd x-valuation forces this form
  to vanish.  Hence x divides f_n and the leading pair lies in the exact
  (K,-a) lift-gauge direction.  This filters the all-degree search but does
  not prove termination or solve the remaining square equation.
source: jc_quartic_c3_construct / post-THM-3881 residual-degree filtration, 2026-08-23
audit: >
  PROVISIONAL EXACT PROOF CANDIDATE.  The companion freezes all fourteen
  nonzero coefficient-degree cells, verifies the all-n strict inequalities,
  factors the equality leading form, checks gcd(x,y^2-15x^2)=1, and verifies
  that the full gauge vector preserves the module address.  Normal and
  optimized runs must byte-match the frozen transcript.  Independent audit
  must recheck the prime-x nonsquare argument, every possible top-degree
  cancellation, and the precise non-invariance warning.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
related:
  - THM-3872-three-cusp-polarization-branches-and-minimal-affine-square-residual-gate
script: 04-computation/jc2_cusp_residual_total_degree_filtration_thm3884.py
output: 05-knowledge/results/jc2_cusp_residual_total_degree_filtration_thm3884.out
script_sha256: bb889185b1c27ba751680d506cb37f14f942e202e8f800a88624864f0e3c8092
output_sha256: 0f4366c1dd234735cafad0722f97c383f084b24f7d7ced65dfcd5e1e7d25249c
semantic_sha256: 0345dc2a0b92899d33c99658ba565f31aaf613cdbca3494b85c1f255171400d1
hash_basis: raw LF bytes
---

# THM-3884 -- total degree exposes the leading lift gauge

**PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.**
Work over an algebraically closed field `k` of characteristic zero.  Retain
the THM-3881 notation

```text
D=k[x,y],        a=x+1,                    L=9x+4,
K=y^2-15x^2-15x-4,                         P=aL^2,
Delta=a^3L^2-K^2,                           r=aT+Kf,
A=KT+aPf.                                                   (1)
```

and let `S(T,f)` be its exact residual

```text
S=L^4+2(3A+3P+r^2)L^2f+(8A+6P+3r^2)(Pf^2-T^2).             (2)
```

The admissible module still has the address condition

```text
T(0,0)=4f(0,0).                                             (3)
```

Suppose `f` is nonconstant and `S(T,f)` is a square in `D`.  Put

```text
n=deg f>=1,                    m=deg T,                      (4)
```

with `deg 0=-infinity`.  Then:

```text
m>=n+1.                                                        (5)
```

If equality holds, write `f_n,T_m` for the highest homogeneous forms and

```text
K_2=y^2-15x^2.                                               (6)
```

Then necessarily

```text
K_2 f_n+xT_m=0.                                              (7)
```

Since `gcd(x,K_2)=1`, there is a homogeneous `q_{n-1}` such that

```text
f_n=xq_{n-1},                    T_m=-K_2q_{n-1}.             (8)
```

Thus the equality seam is precisely the leading homogeneous direction of the
THM-3881 lift gauge `(K,-a)`.

## 1. The complete coefficient-degree packet

Regard `(2)` as a polynomial in `T,f` over `D`.  For every nonzero cell
`T^i f^j`, the total degree of its coefficient is as follows:

```text
(i,j):degree
(4,0):2   (3,1):3   (3,0):2
(2,2):5   (2,1):4   (2,0):3
(1,3):6   (1,2):5   (1,1):4
(0,4):7   (0,3):7   (0,2):6   (0,1):5   (0,0):4.            (9)
```

The four highest-coefficient identities behind the proof are

```text
[f^4]S       =3aK^2L^2,
[Tf^3]S      =6Ka^2L^2,
[T^2f^2]S   =3(a^3L^2-K^2),
[T^4]S       =-3a^2.                                      (10)
```

No omitted cell can compete with the degree comparisons below; `(9)` is the
complete universe.

## 2. The region `m<=n` has odd top degree

Assume `m<=n`.  The `f^4` cell in `(10)` has degree

```text
4n+7.                                                       (11)
```

Every other cell has degree at most `4n+6`.  This follows directly from
`(9)` after replacing `m` by its largest allowed value `n`; the inequality is
already strict at `n=1` and its slope in `n` is nonnegative.  The leading form
is explicitly

```text
243x^3 K_2^2 f_n^4,                                        (12)
```

and is nonzero in the domain `k[x,y]`.  In particular `deg S=4n+7` is odd,
whereas a nonzero polynomial square has even total degree.  This contradiction
proves `(5)`.

## 3. Equality forces the gauge seam

Now take `m=n+1`.  Exactly three cells in `(9)` can reach degree `4n+7`:

```text
f^4,                         Tf^3,                         T^2f^2.            (13)
```

Using `(10)` and the leading forms

```text
a_1=x,                L_1=9x,                K_2=y^2-15x^2,
(a^3L^2-K^2)_5=x^3(9x)^2,                                  (14)
```

their sum factors exactly as

```text
S_{4n+7}
 =3a_1L_1^2 f_n^2(K_2f_n+a_1T_m)^2
 =243x^3f_n^2(K_2f_n+xT_m)^2.                              (15)
```

If the last parenthesis is nonzero, its valuation at the prime `x` is

```text
3+2v_x(f_n)+2v_x(K_2f_n+xT_m),                             (16)
```

which is odd.  Hence `(15)` is not a square in the UFD `k[x,y]`.  But the
highest homogeneous form of a polynomial square must itself be a square.
Therefore `(15)` must vanish.  Since `f_n!=0`, the domain property gives
`(7)`.  Coprimality `gcd(x,K_2)=1` then proves `(8)`.

## 4. Recursive interpretation and exact boundary

The THM-3881 gauge acts by

```text
(T,f) |-> (T+Kq,f-aq).                                     (17)
```

Given `(8)`, choose `q` with highest form `q_{n-1}`.  Then `(17)` cancels both
`T_m` and `f_n`, lowering the displayed leading pair.  The address `(3)` is
preserved because at the origin

```text
K(0,0)=-4,                    a(0,0)=1.                    (18)
```

This supplies an iterative **bookkeeping** operation: every equality-seam
candidate exposes and can peel one leading gauge layer.  It is not a proof
that the process preserves the square equation.  The gauge fixes `r` but
changes `A` by `-Delta q`, so `S` changes by the explicit THM-3881 transport.
No termination, descent of a survivor, or general no-go is claimed.

The theorem leaves two live regions:

1. `deg T>=deg f+2`;
2. equality-seam candidates together with the exact residual debt generated
   when their leading gauge layer is peeled.

Constant `f`, alternate square/cube lifts, a Keller atlas, and JC(2) are also
outside the claim.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_cusp_residual_total_degree_filtration_thm3884.py
python3 -O 04-computation/jc2_cusp_residual_total_degree_filtration_thm3884.py
```

Both streams must byte-match
`05-knowledge/results/jc2_cusp_residual_total_degree_filtration_thm3884.out`.
