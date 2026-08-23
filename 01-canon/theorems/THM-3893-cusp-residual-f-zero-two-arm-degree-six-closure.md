---
id: THM-3893
title: "Cusp-residual f-zero two-arm degree-six closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the THM-3885
  f=0 residual lane, if both exact arm restrictions choose their zero branch,
  then every square survivor of total degree at most six is the base point
  T=0.  Arm divisibility and an odd L-adic response force T=aL^2H; an exact
  cubic univariate square elimination forces x|H; the remaining y-degrees
  zero, one, and two die respectively by missing-y, odd-degree, and odd
  a-adic leading-coefficient parity.  Degree seven, either nonzero arm branch,
  general f, a Keller atlas, and JC(2) remain OPEN.
source: tournament-jc-breakthrough / post-THM-3885 two-arm deletion response, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root plus forgotten_tournament_signal,
  2026-08-23).  Both paths rederived the two-arm factorization, the exact
  L-valuation step, and the complete degree-four/five boundary.  A separate
  root replay built the cubic square root without candidate-coefficient
  division, recomputed the 15-element grevlex Groebner basis from all six
  residual equations, and checked the exhaustive global y-degree split and
  exact a-valuation three at degree six.  The companion's normal and optimized
  runs byte-match the frozen transcript.  No saturation or branch division is
  used.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure
related:
  - THM-3888-f-zero-equianharmonic-jacobian-and-two-section-integrality
  - THM-3886-cusp-residual-equality-seam-second-layer-trichotomy
script: 04-computation/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.py
output: 05-knowledge/results/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.out
script_sha256: bfe09cbdba574546d295fff2cfbd47421ab59b7f6f85e0078fe1db072e9ef9c8
output_sha256: a10844f8243879f02f464eb92fbb4a97ac101c200b3c05b4a91e56e2236ba634
semantic_sha256: 0908602c9e7c1a74d6cf5a82bd7e16904c810e242acd16ea8156d98470d86221
hash_basis: raw LF bytes
---

# THM-3893 -- the simultaneous zero arms are empty through degree six

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero.  Put

```text
D=k[x,y],                  a=x+1,
L=9x+4,                    K=y^2-15x^2-15x-4.             (1)
```

In the `f=0` lane of THM-3881, the residual square equation is

```text
S(T,0)=L^4-6aL^2T^2-8KT^3-3a^2T^4,                      (2)
```

and the module address is

```text
T(0,0)=0.                                                   (3)
```

Suppose `S(T,0)` is a square in `D`, both exact arm restrictions take their
zero branch,

```text
T(-1,y)=0,                 T(-4/9,y)=0,                    (4)
```

and

```text
deg T<=6.                                                       (5)
```

Then

```text
T=0.                                                         (6)
```

## 1. The first arm summary loses a kernel; its L-response restores it

The two evaluations in `(4)` give `a|T` and `L|T`.  Since `a,L` are coprime,
write

```text
T=aLH.                                                       (7)
```

Substitution in `(2)` factors exactly as

```text
S=L^3 B,
B=L(1-6a^3H^2-3a^6H^4)-8a^3KH^3.                         (8)
```

At `L=0`, equivalently `x=-4/9`,

```text
B mod L
 =-8(5/9)^3(y^2-8/27)(H mod L)^3.                         (9)
```

If `L` did not divide `H`, the right side of `(9)` would be nonzero in the
domain `k[y]`.  Then `v_L(S)=3`, impossible for the nonzero square `S`; it is
nonzero because `(3)` gives `S(0,0)=256`.  Therefore

```text
H=LH_1,
T=aL^2H_1,
H_1(0,0)=0.                                                (10)
```

This is the exact deletion-response mechanism.  The two arm values alone
lose the kernel `(aL)`, while the first `L`-adic coefficient of the residual
forces the additional factor.

## 2. Exact cubic specialization at x=0

Since `deg(aL^2)=3`, equations `(5),(10)` give

```text
deg H_1<=3.                                                (11)
```

At `x=0`, write

```text
tau(y)=T(0,y)=py+qy^2+ry^3.                               (12)
```

The constant term vanishes by `(3)`.  Specializing `(2)` gives

```text
F(tau)=256-96tau^2-8(y^2-4)tau^3-3tau^4.                 (13)
```

We need the following exact lemma.

**Cubic-tau lemma.**  If `(13)` is a square in `k[y]`, then

```text
p=q=r=0.                                                   (14)
```

Indeed, `F(0)=256`.  Negate a square root if needed and write its complete
degree-six universe as

```text
G=16+g_2y^2+g_3y^3+g_4y^4+g_5y^5+g_6y^6.                (15)
```

There is no linear term because `[y]F=0`.  For `j=2,...,6`, the coefficient
equation determines `g_j` recursively:

```text
g_j=([y^j]F-[y^j](16+g_2y^2+...+g_(j-1)y^(j-1))^2)/32.  (16)
```

Thus the recursion divides only by a characteristic-zero scalar, never by
`p,q,r`.  Let

```text
E_j=[y^j](G^2-F),                    7<=j<=12.             (17)
```

Exact Groebner reduction in `Q[p,q,r]`, with variable order `(p,q,r)` and
`grevlex`, gives a 15-element basis for `(E_7,...,E_12)`.  It contains

```text
r^4,
q^4-6p^2r^2,
p^5+16p^2r+16pq^2+(448/5)r^3.                             (18)
```

Every field-valued common zero therefore has successively `r=0`, `q=0`, and
`p=0`.  No saturation is involved.  This proves the lemma and hence

```text
H_1(0,y)=0,
H_1=xJ,                         deg J<=2.                  (19)
```

## 3. The three remaining y-profiles

If `J=0`, then `(6)` is immediate.  Otherwise put `d=deg_y J`.

### 3.1. The univariate profile

If `d=0`, then `T` belongs to `k[x]`.  Equation `(2)` has

```text
[y]S=0,                     [y^2]S=-8T^3,                 (20)
```

and `S(0,0)=256`.  If `T!=0` and `S=G_1^2`, then
`G_1=u(x)y+v(x)`.  The missing linear term gives `uv=0`, while `(20)` gives
`u!=0`; hence `v=0`, contradicting the nonzero constant value.  Thus `T=0`.

### 3.2. The linear-y profile

If `d=1`, put `t_1(x)=[y]T`, which is nonzero.  The unique highest y-term in
`(2)` is

```text
[y^5]S=-8t_1^3.                                            (21)
```

Consequently `deg_y S=5`, impossible for a square.

### 3.3. The quadratic-y profile

If `d=2`, the coefficient `b=[y^2]J` is a nonzero scalar because
`deg J<=2`.  Put

```text
t_2=[y^2]T=baL^2x.                                         (22)
```

The complete leading y-coefficient is

```text
[y^8]S=-8t_2^3-3a^2t_2^4
       =-t_2^3(8+3a^2t_2).                                (23)
```

At the prime `a=x+1`, `v_a(t_2)=1`: on `a=0`, both `x=-1` and `L=-5` are
units.  The second factor in `(23)` is `8 mod a`, so

```text
v_a([y^8]S)=3.                                             (24)
```

But the leading coefficient in `y` of a polynomial square must be a square
in `k[x]`, and hence have even valuation at every prime.  This contradiction
closes the last profile and proves `(6)`.

## 4. Exact boundary and elliptic interpretation

THM-3885 already closes the entire `f=0` lane through degree three.  The new
content here is the simultaneous-zero-arm sector in degrees four, five, and
six.  It follows from `(10)` that the first possible nonzero survivor in this
sector has

```text
deg T>=7,                       deg H_1>=4.                (25)
```

At degree seven the `x=0` specialization is quartic and the mixed y-channels
can tie, so neither the cubic-tau elimination nor the unique odd-degree
argument extends formally.

THM-3888 identifies the generic `f=0` equation over `k(x)` with a rank-six
rational elliptic surface and a two-boundary-section integrality problem.
The present theorem supplies a different sidecar: simultaneous evaluation at
the labelled fibres `a=0,L=0` plus the first `L`-adic response enforces global
`k[x,y]` descent through degree six.  It does not enumerate the generic
Mordell--Weil group or prove an all-degree height bound.

Nonzero choices on either arm, degree at least seven, `f!=0`, alternate
square/cube lifts, a polynomial-plane Keller atlas, and `JC(2)` remain
**OPEN**.

Reproduce the exact cubic elimination and complete degree-six profile with

```bash
python3 04-computation/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.py
python3 -O 04-computation/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.py
```

Both streams must byte-match
`05-knowledge/results/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.out`.
