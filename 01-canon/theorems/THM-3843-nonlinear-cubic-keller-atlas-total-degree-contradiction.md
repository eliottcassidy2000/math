---
id: THM-3843
title: "The nonlinear cubic Keller atlas degree packet is inconsistent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  No dominant etale
  morphism from the polynomial plane to the THM-3811 nonlinear cubic surface
  exists.  The THM-3836 factor identity and Keller cofactor are inconsistent
  in all three total-degree orders deg(h)>deg(k), deg(k)>deg(h), and equality.
  This closes this candidate surface only; it does not prove JC(2).
source: root + audit_thm3838 / incoming factor-cofactor degree-trichotomy lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / audit_thm3838, 2026-08-23).  The
  audit checked nonzero/nonconstant endpoints, exact product degrees, every
  unequal-degree dominance comparison, the equal-degree cancellation in
  Q=7h^2+3k^2, its noncancelling conjugate factor, the resultant excluding
  simultaneous cancellation in P, the Jacobian degree bound, and the
  dominance-versus-etaleness scope split.  The assertion-free companion
  verifies the polynomial packet, resultants, symbolic degree deficits, a
  216,000-triple hostile census, and the determinant-only/Keller-failing
  constant control.  Normal and optimized replay agree.
depends_on:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
  - THM-3836-cubic-factor-cofactor-darboux-packet
related:
  - THM-3835-polynomial-marked-root-ratio-nonentry
  - THM-3838-root-ratio-numerator-denominator-degree-five-floor
  - THM-3840-forced-cubic-two-arm-jelonek-passport
script: 04-computation/jc2_nonlinear_cubic_keller_atlas_degree_contradiction_thm3843.py
output: 05-knowledge/results/jc2_nonlinear_cubic_keller_atlas_degree_contradiction_thm3843.out
script_sha256: a8bac4a4f391d74b83b47770f8d60c1a5e20aec54efc789ed767dfb1756e2e8b
output_sha256: 406d49f29a932804a7a3f9596dc8fe298f72c3196dcf49efe4eeaed5bdbf9609
semantic_sha256: db743214723ec651d2828740145b62e94d51516c17bc6c94737e2404109a229e
hash_basis: raw LF bytes
---

# THM-3843 -- the factor/cofactor degrees exclude every polynomial atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `K` of characteristic zero.  There is no dominant
etale morphism

```text
psi:A2_(x,y) -> U                                               (1)
```

to the THM-3811 nonlinear cubic surface.  Equivalently, its polynomial Keller
atlas problem is empty.  This theorem excludes this particular candidate
surface, not arbitrary polynomial Keller pairs; `JC(2)` remains open.

## 1. The exact packet and degree coordinates

Suppose `(1)` existed.  Write `h,k,C in K[x,y]` for the pulled-back intrinsic
functions and put

```text
Q=7h^2+3k^2,
B=6h^3+7h^2k-k^3,
P=3h^3+7h^2k+k^3.                                             (2)
```

The intrinsic `h,k,C` are nonzero and nonconstant: the `h`- and `k`-zero
fibres are proper Laurent arms, while THM-3832 gives `K(U)=K(z,C)`.
Dominance makes the coordinate-ring pullback injective, so their images remain
nonzero and nonconstant.  In particular their total degrees

```text
a=deg h,                    b=deg k,                    c=deg C (3)
```

are positive.  Both `psi` and the intrinsic map `U->A2_(A,C)` are etale, so
their composite is etale.  Its polynomial Jacobian is therefore a unit
`lambda in K*`.  THM-3836 supplies

```text
C^2 h^2 Q=CkB+P,                                               (4)
delta(C):=k Jac(h,C)-h Jac(k,C)=lambda P.                       (5)
```

Equation `(4)` is just `P=CS` with `S=Ch^2Q-kB`.  Equation `(5)` is the
Keller input.  The elementary bound used throughout is

```text
deg delta(C) <= a+b+c-2.                                       (6)
```

## 2. Unequal row degrees

Suppose first that `a>b`.  No leading cancellation is possible in `(2)`, so

```text
deg(Q,B,P)=(2a,3a,3a).                                        (7)
```

The left side of `(4)` has degree `2c+4a`.  On the right, `CkB` has degree
`c+b+3a>3a=deg P`, so it uniquely supplies the leading term.  Equality of
degrees gives

```text
2c+4a=c+b+3a,             hence c=b-a<0,                       (8)
```

contrary to `(3)`.

Now suppose `b>a`.  This time

```text
deg(Q,B,P)=(2b,3b,3b),                                        (9)
```

and `(4)` gives

```text
2c+2a+2b=c+4b,            hence c=2(b-a).                     (10)
```

But `(6),(10)` imply

```text
deg delta(C) <= a+b+c-2=3b-a-2<3b=deg P,                      (11)
```

contradicting `(5)`.

## 3. Equal row degrees and the conjugate-factor floor

It remains to take `a=b=d`.  The polynomial `Q` is nonzero: if `Q=0`,
algebraic closedness makes `h` a scalar multiple of `k`, and the intrinsic
unimodular-row law `Ck-mh=1` then makes both scalar, contradicting dominance.
Write `q_0=deg Q`.

Both `B` and `P` have degree at most `3d`, so `(4)` first gives

```text
2c+2d+q_0 <= c+4d,          q_0<=2d-c<2d.                    (12)
```

Thus the degree-`2d` terms of `Q` cancel.  If `H,K` are the leading forms of
`h,k`, there is a `beta in K*` with

```text
H=beta K,                     beta^2=-3/7.                     (13)
```

The cancellation cannot propagate to `P`, because

```text
Res_Z(7Z^2+3,3Z^3+7Z^2+1)=1615!=0.                            (14)
```

Hence `deg P=3d`.  Equation `(5)` and `(6)` now force

```text
3d<=2d+c-2,                  c>=d+2.                           (15)
```

On the other hand, over `K`,

```text
Q=7(h-beta k)(h+beta k).                                      (16)
```

The second factor has degree exactly `d`, while the first is nonzero by the
same unimodular-row argument.  Therefore `q_0>=d`.  Combined with `(12)`,
this gives `c<=d`, contradicting `(15)`.

The three cases exhaust all degree orders, proving the theorem.  Notice the
scope boundary: dominance and the factorization `(4)` alone are not enough;
the cofactor equality `(5)`, hence the Keller/etale hypothesis, is
load-bearing.  Rational source charts and unrelated Keller surfaces are not
excluded.  **QED.**
