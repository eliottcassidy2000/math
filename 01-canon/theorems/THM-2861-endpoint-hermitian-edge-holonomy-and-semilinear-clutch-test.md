---
id: THM-2861
title: "Endpoint Hermitian edge holonomy and semilinear-clutch test"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT.  The actual THM-2847/2851 carry is
  scalar-trivial on the endpoint field and has no nonzero q7 endpoint in
  the old E3 block.  On the THM-2857 Galois torsor, however, adjacent
  sections have a gauge-invariant bidegree-(1,1) Hermitian edge with fixed
  phase omega^3, thirteen distinct values, and Fourier support {0,3,10}.
  A lawful adjacent-section co-support pair is not supplied.  Independent
  hostile audit is pending.
source: root/lrc-endpoint-hermitian-edge-holonomy-2026-07-28
depends_on:
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2857-endpoint-galois-carry-torsor-and-phase-alignment-sidecar
related:
  - THM-2460-idempotent-semantic-word-copy-and-word-block-cosupport-boundary
  - THM-2744-relative-present-unit-repair-and-root-zero-overlap-clutch
  - THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go
script: 04-computation/lrc14_endpoint_hermitian_edge_holonomy_thm2861.py
output: 05-knowledge/results/lrc14_endpoint_hermitian_edge_holonomy_thm2861.out
script_sha256: 57bad76968ec9c61d2202331e007860f2817d15c606d8ba558ab8b8d3c41f20c
output_sha256: 9bc846b6269b6ca967d32b5b4091ec506b3ede632c58a249c78211e1ecc8b43d
hash_basis: LF-normalized bytes
---

# THM-2861 -- endpoint Hermitian edge holonomy and semilinear-clutch test

**PROVED CANDIDATE + VERIFIED-EXACT.**

The endpoint carry can be detected at bidegree `(1,1)`, rather than by the
tenth power used for q7-to-q11 label alignment.  The detector is invariant
under a common unknown endpoint phase, has only three Fourier channels, and
still separates all thirteen carry sections.  Its input is precisely what
canon lacks: two adjacent semilinear sections on a common physical support.

The theorem also completes THM-2857's cheapest physical test.  The existing
q11-to-q7 ancestry carry does not produce `sigma_1(c)`: algebraically it is
linear over the endpoint field, and physically its q7 endpoint projection
is zero.

## 1. The actual canonical carry has no endpoint image

Write

```text
c=zeta_1183^624-zeta_1183^510.
```

Exact replay of the THM-2806/2829/2847 endpoint builders on all `42`
q3/q11 horn cells gives

```text
origin (0,0):   q3=c, q11=c, q7=0,
origin (12,0):  q3=0, q11=c, q7=0.                       (1)
```

Every nonzero row in `(1)` is one interval with reduced endpoint data

```text
conductor 2366,     exponents (2203,65),
zeta_1183 form      (624,510).                            (2)
```

Let `u` be one allocation-root step and `T_phys=13u` the whole physical
period.  With the exact THM-2847 endpoint frequency and dilation,

```text
Y0*R*u      =83486*N,
Y0*R*T_phys =1085318*N,                                  (3)
```

where `N` is the endpoint conductor before reduction.  Thus both physical
translations add zero to the cyclotomic exponents.  They leave `c` fixed.

On the coefficient/support side, THM-2851's natural lift is

```text
(a,11) --h=9--> (a+1,7).
```

Its proved oriented mask is `h_L(Z)` and the carry is multiplication by
`Z` over `K_0=Q(zeta_1183)`:

```text
T(c h_L)=c Z h_L.                                        (4)
```

This differs from

```text
sigma_1(c)=zeta_1183^624-zeta_1183^783.                   (5)
```

The difference in `(5)` is nonzero already because it is represented by a
nonzero polynomial of degree `783<phi(1183)=936`.

At the fully typed level, there is no alternative nonzero scalar: on `20`
of the `42` cells q7 loses exactly E3, and on the other `22` it loses E3
and at least one ordinary factor.  Equation `(1)` shows that q7 also has
zero endpoint support at both origins.  The physical carried endpoint is
therefore undefined, or zero after projection to the old endpoint/E3 block.

## 2. The adjacent Galois edge

Now use the coefficient-field torsor of THM-2857:

```text
c_r=A-B omega^(3r),       A=zeta^624, B=zeta^510,
omega=zeta^91=zeta_13.                                   (6)
```

Complex conjugation gives the exact identity

```text
bar(c_r)=-c_r/(AB omega^(3r)).                            (7)
```

Consequently the oriented Hermitian edge

```text
E_r=c_(r+1) bar(c_r)
```

satisfies

```text
E_r=omega^3 bar(E_r)
   =omega^3 bar(c_(r+1))c_r.                              (8)
```

In particular,

```text
J_r:=E_r-bar(E_r) !=0,
S_r:=E_r+bar(E_r) !=0,
J_r/S_r=(omega^3-1)/(omega^3+1).                          (9)
```

Thus one complex quadrature detects the orientation of every adjacent
carry edge.  Reversing the edge conjugates `E_r` and changes the phase from
`omega^3` to `omega^-3`.

This observable removes the common-phase debt.  If every endpoint scalar
is multiplied by one unknown `lambda!=0`, then

```text
(lambda c_(r+1)) overline(lambda c_r)
   =|lambda|^2 E_r,                                      (10)
```

so `(8)--(9)` are unchanged.  The same phase law survives arbitrary
positive real amplitude weights on the separate sections.

## 3. Three Fourier channels still separate thirteen sections

Expanding `(8)` gives

```text
E_r =
  1+omega^3
  -(A/B)omega^(-3r)
  -(B/A)omega^(3(r+1)).                                  (11)
```

Hence

```text
supp_Fourier(r |-> E_r)={0,3,10}.                         (12)
```

The thirteen values are nevertheless distinct.  If `E_r=E_s`, put
`u=omega^(3r)` and `v=omega^(3s)`.  Factoring their difference shows that
either `u=v`, hence `r=s`, or

```text
(A/B)^2=omega^3uv.                                       (13)
```

The left side of `(13)` is `zeta^228`, whose exponent is not divisible by
`7`; the right side is a power of `omega` and has exponent divisible by
`91`.  So the second case is impossible.

The antisymmetric current `J_r` in `(9)` has the same support `(12)` and
is nonzero for all thirteen sections.  Therefore the Hermitian edge is a
point-separating, orientation-sensitive carry coordinate at bidegree
`(1,1)`.

## 4. Sharp service boundary

For the existing `K_0`-linear carrier, the scalar sequence is constant:
`d_r=c`.  Its Hermitian edge is the constant real value `|c|^2`, so

```text
J_r=0,                  supp_Fourier(E)={0}.               (14)
```

This is the exact sharp hostile.

There is a second, physical obstruction.  THM-2851 proves
`Lambda` and `Lambda+1` are disjoint (`156` versus `157 mod169`).
Consequently the naive same-address inner product of the two positive horn
masks is zero.  One cannot obtain `(8)` merely by multiplying their
coefficient arrays.

A lawful realization must instead provide an adjacent-section co-support
edge: transport both sections to a common physical atom and common target
phase while retaining the semilinear endpoint action.  Under the separate
support and quadrature hypotheses of THM-2380, polarization of that pair
would recover `E_r`; `(8)` would then detect carry without an absolute
phase reference or a tenfold product.  No such co-support transport is
proved here.

## 5. Exact verification

The companion pins THM-2847, THM-2851, and THM-2857.  It:

1. replays the full `42=20+22` physical horn split and the six endpoint
   rows `(1)`;
2. recomputes `(2)--(5)` and both zero physical phase increments;
3. proves `(8)` for all thirteen sections in the exact group ring of
   `mu_1183`;
4. checks `(8)--(12)`, thirteen-value separation, and common-gauge
   cancellation in an independent exact `F_4733` specialization;
5. verifies the constant-scalar hostile `(14)`.

Normal, optimized, and stored-output replay are byte-identical after LF
normalization.  LF-normalized SHA-256:

```text
script  57bad76968ec9c61d2202331e007860f2817d15c606d8ba558ab8b8d3c41f20c
output  9bc846b6269b6ca967d32b5b4091ec506b3ede632c58a249c78211e1ecc8b43d
```

## 6. Connection contract

```text
source:
  THM-2847 endpoint scalar c and THM-2851 oriented carry;

target:
  a common-support, phase-sensitive adjacent-section carry observable;

algebraic map:
  (c_r,c_(r+1)) |-> E_r=c_(r+1)bar(c_r);

preserved:
  carry orientation, all thirteen section labels, and common phase gauge;

proved failure of the current map:
  T is K0-linear, q7 has zero endpoint support, and the raw positive
  ancestry masks are disjoint;

needed sidecar:
  semilinear endpoint descent plus E3/ordinary-factor co-support transport
  and a lawful Hermitian quadrature;

cheapest decisive next test:
  construct one common physical atom for adjacent sections and check
  whether its polarized edge obeys E=omega^3 bar(E).
```

No row is excluded and the LRC(14) ledger remains `165`.
