---
id: THM-4368
title: "Diagonal boundary-valuation triangular address and simplex-stream rank"
status: >
  PROVED RELATIVE TO THM-4364 + VERIFIED-EXACT + INDEPENDENTLY AUDITED; JC(2)
  OPEN. Complete fixed-order streams are invertible transforms of the raw
  diagonal, consecutive fixed-row order banks have unit minors, and the exact
  hierarchy-bank rank forms a triangular clock. Normalized packet traces are
  classified by two boundary valuations with an exact source-overlap fibre
  and triangular natural address. Reciprocal address reflection is ambient,
  not always source-realizable. No bracket, seam, JC(2), or DC(2) conclusion
  is asserted.
source: root + boundary-address scout + clean-room and hostile referees / next-sharp session, 2026-09-03
depends_on:
  - THM-4364-sharp-binomial-diagonal-annihilator-hierarchy
related:
  - THM-4366-source-normal-u-zero-row-eleven-hierarchy-selected-extinction
  - THM-4367-lrc14-active-first-exit-scale-collision-classification
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
  - THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density
mistake_firewall:
  - MISTAKE-354
primary_script: 04-computation/jc2_diagonal_boundary_address_simplex_stream_rank_thm4368.py
primary_output: 05-knowledge/results/jc2_diagonal_boundary_address_simplex_stream_rank_thm4368.out
primary_script_sha256: b24b707ffd946d6d24fb2b452f55810fffb746a6cc1c7308a0256c39d7777203
primary_output_sha256: 14945ef3c11af0f17dabb599b54351d37a85b73d9560cd1ee1f30018e6459ee8
independent_referee_script: 04-computation/jc2_diagonal_boundary_address_simplex_stream_rank_independent_referee_thm4368.py
independent_referee_output: 05-knowledge/results/jc2_diagonal_boundary_address_simplex_stream_rank_independent_referee_thm4368.out
independent_referee_script_sha256: f80ea54633f7c95c0a67a913c494baf3be5c99dae1883339d3549d8343faed69
independent_referee_output_sha256: a86f08ddba90a66f4c94dcfed63ad7ac48efd653a0c8abc892aef671ac1345d4
hash_basis: raw LF bytes
audit: >
  PASS WITH TRACE-SCOPE, FIXED-BANK, AND AMBIENT-REFLECTION WORDING REPAIRS.
  The 6,790,755-check primary and 14,325,934-check import-free referee rebuild
  the stream transforms, unit minors, rank clock, boundary formulas,
  feasibility and overlap fibres, address/reflection, and hostile examples.
  Normal, optimized, isolated, hash-seeded, and frozen LF streams agree.
---

# THM-4368 -- Diagonal boundary-valuation triangular address and simplex-stream rank

**PROVED RELATIVE TO THM-4364 + VERIFIED-EXACT + INDEPENDENTLY AUDITED. THIS
IS A THEOREM ABOUT THE DIAGONAL CONSUMER HIERARCHY AND ITS MONOMIAL PACKETS.
IT DOES NOT ASSERT A NEW BRACKET OBSTRUCTION, `JC(2)`, OR `DC(2)`.**

## 1. Statement and inheritance

Fix THM-4364's diagonal intercept `ell>=2` and write

```text
s=ceil(ell/2),                     rho=ceil(ell/3).
```

All row/order/depth parameters below are integers with `m>=s` and `q,d>=0`;
all source exponents are nonnegative integers.

For an arbitrary coefficient stream `h`, put

```text
A(n)=(-1)^(n-s) h_(n,2n-ell),                 n>=s,

L_q(m)=sum_(n=s)^m C(m+q-n,q) A(n).                       (1)
```

Thus `L_q(m)=L_(m,ell,q)(h)` in THM-4364. The theorem has four parts.

1. Every complete fixed-`q` row stream is an integral invertible transform
   of the raw diagonal stream. At one fixed row, consecutive `q` banks have
   primitive unit minors.
2. As ambient diagonal-coordinate functionals, the THM-4364 hierarchy
   annihilators of `pi_m(P_d)` have exact independent bank rank

   ```text
   max(0,rho-max(0,ell+d-m)).                              (2)
   ```

3. A normalized nonzero single diagonal packet has trace-equivalence state
   exactly `(N,n0)`, recovered as its two boundary valuations.
4. In positive boundary coordinates `(u,v)`, that state has a triangular
   natural address. Coordinate swap is an ambient reciprocal reflection, but
   it need not preserve the source-monomial cone. Order and primitive-ray
   quotients both have explicit trace hostiles.

The inheritance pass was:

- closest proved mechanism: THM-4364's packet expansion and partial-binomial
  identities;
- canonical hostile: THM-4088's order-blindness pattern, made internal to
  two diagonal packet traces;
- corrected near miss: reflection of the address space need not be realized
  by a source monomial;
- least-used sidecar: the two boundary multiplicities of the whole row trace.

The live board was

```text
row stream | simplex order | Pascal unit minor | boundary valuation
triangular address | reciprocal reflection | primitive ray | source cone. (3)
```

## 2. Full streams and fixed-row Pascal duality

Let `Sigma` denote prefix summation from row `s`, and set

```text
Delta f(m)=f(m)-f(m-1),                     f(s-1)=0.
```

Pascal's identity applied to `(1)` gives

```text
Delta L_0=A,
Delta L_q=L_(q-1),                         q>=1.          (4)
```

Consequently

```text
L_q=Sigma^(q+1)A=Sigma^q L_0,
A=Delta^(q+1)L_q.                                           (5)
```

On any finite row interval `s,...,M`, the map from `A` to a fixed `L_q`
stream is lower unitriangular over the integers. Hence

```text
ker(h -> (A(n))_(n>=s))=ker(h -> (L_q(m))_(m>=s))        (6)
```

for every fixed `q`, and also for the joint family over all `m,q`. Varying
simplex order does not add information when a complete row stream is kept.
This does not identify individual fixed-row functionals or their different
annihilation ranges.

There is an exact fixed-row dual. Fix `m`, put `r=m-s`, reverse columns by
`k=m-n`, and ignore the invertible column signs. For `0<=q,k<=r`, the
consumer matrix is

```text
B_(q,k)=C(q+k,q).                                         (7)
```

If `P_(q,j)=C(q,j)`, Vandermonde convolution gives

```text
B=P P^T.                                                  (8)
```

The Pascal matrix `P` is lower unitriangular, so `det(B)=1`. More generally,
for `q0>=0` and `1<=a<=r+1`, take consecutive orders
`q=q0,...,q0+a-1` and the first `a` reversed columns.
Their determinant is

```text
prod_(i<j)(j-i) / prod_(k=0)^(a-1) k! =1.                (9)
```

Indeed column `k` is the degree-`k` polynomial `C(x+k,k)` with leading
coefficient `1/k!`, evaluated at consecutive integers. Thus every consecutive
bank is primitively independent, not merely independent over the fraction
field.

## 3. The exact hierarchy-annihilator clock

THM-4364 says that at fixed `(ell,m,d)`, the hierarchy orders annihilating
`pi_m(P_d)` are exactly

```text
q0=max(0,ell+d-m), ..., rho-1.                            (10)
```

Equation `(9)` proves their primitive independence as ambient
diagonal-coordinate functionals. The rank of this hierarchy bank is therefore
exactly `(2)`; this is not a claim that the bank is the full annihilator space
of `pi_m(P_d)`.

For fixed `(ell,d)`, put

```text
m0=ell+d-rho+1.                                          (11)
```

The bank is empty before `m0`. At row `m=m0+j`, for `0<=j<=rho-1`, it has
the `j+1` independent orders

```text
q=rho-1-j,...,rho-1.                                     (12)
```

It has stable rank `rho` from row `ell+d` onward. The ramp therefore has
sizes `1,2,...,rho` and exactly `T(rho)` marked row/order positions. The
familiar consumers are its first-entry positions:

```text
(ell,d,m,q)=(8,2,8,2),        triangular row-eight,
              (10,2,9,3),     tetrahedral row-nine A,
              (10,3,10,3),    row-ten C,
              (11,2,10,3),    opposite-order row-ten A. (13)
```

In particular, each displayed THM-4366 fixed-`(ell,m,d)` bank has rank one:
the `q=3` functional is unique within that bank. This is not uniqueness in
the full module nullspace or across other diagonal intercepts. Their `-4/3`
proportionality after source restriction is not implied by `(7)`--`(9)`; it
remains a special fact about that selected source graph.

## 4. A packet is a two-boundary polynomial

For a source monomial

```text
x^a u^b p^c y^e
```

whose packet meets the diagonal, THM-4364 gives

```text
a=2c+3e-ell,
N=c+e,                  n0=b+c+2e,          n1=n0+N.     (14)
```

Here `N>0`. Direct finite summation yields

```text
L_0(m)=(-1)^(m-s) C(N-1,m-n0),       n0<=m<n1,
       0,                             otherwise.          (15)
```

Hence the row generating functions are

```text
F_q(z)=sum_(m>=s)L_q(m)z^(m-s)
      =(-1)^(n0-s) z^(n0-s)(1-z)^(N-q-1).                (16)
```

For `q>=N`, `(16)` is a formal power series. The transitions at
`q<N-1`, `q=N-1`, and `q>=N` are respectively compact support, one impulse,
and an infinite tail. Prefix summation in `(5)` proves all three cases,
including `N=1`.

At `q=0`, the two boundary valuations recover

```text
n0=s+ord_(z=0) F_0,
N =1+ord_(z=1) F_0.                                     (17)
```

Thus `(N,n0)` is exactly the trace-equivalence state for normalized nonzero
single diagonal packets, where trace-equivalence means equality of the
complete `L_0` stream (equivalently any complete fixed-`q` stream). A packet
with arbitrary nonzero scalar coefficient additionally needs that scalar.
Every off-diagonal monomial belongs to one zero-trace class, and no
classification of arbitrary linear combinations is asserted.

The pair `(N,n0)` is source-realizable if and only if

```text
N>=rho,                  n0>=N,             n0+N>=ell.   (18)
```

Necessity follows from `(14)`. For sufficiency choose

```text
e=max(0,ell-2N),       c=N-e,
a=2N+e-ell,            b=n0-N-e.                        (19)
```

The three inequalities in `(18)` make all four exponents nonnegative and
reconstruct `(14)`.

The overlap of source monomials at one valid packet type is also exact. All
realizations, and no others, are indexed by

```text
max(0,ell-2N)<=e<=min(N,n0-N),                            (19a)
c=N-e,              a=2N+e-ell,          b=n0-N-e.
```

Their number is

```text
mu(N,n0)=min(N,n0-N)-max(0,ell-2N)+1.                    (19b)
```

Every realization has the identical expanded packet

```text
x^(2n0-ell)t^n0(1+x^2t)^N.
```

Thus the overlap is genuine trace equivalence, not a collision accidentally
created by the natural-number address. The address does not recover the
individual source exponents.

## 5. Triangular address and reciprocal boundary reflection

Put the two positive boundary coordinates

```text
u=n0-s+1,                         v=N.                   (20)
```

Their antidiagonal address is

```text
Addr(u,v)=T(u+v-2)+u,             T(k)=k(k+1)/2.         (21)
```

This is a bijection from positive ordered pairs to positive natural numbers.
If `tau=u+v-1`, its block is

```text
T(tau-1)+1,...,T(tau),
```

and the decoder is

```text
u=Addr-T(tau-1),                 v=tau+1-u,              (22)
```

where `tau` is the unique integer with `T(tau-1)<Addr<=T(tau)`. In the
integer extension, `T(-tau)=T(tau-1)`.

Swapping the boundary orders gives the ambient involution `R(u,v)=(v,u)`.
It reciprocates `u/v`, preserves `tau`, and reflects the address within its
triangular block:

```text
Addr(u,v)+Addr(v,u)=tau^2+1.                             (23)
```

At trace level, the ambient reflected type satisfies

```text
F_R(z)=(-1)^(v-u) F_0(1-z).                              (24)
```

For odd `tau`, `u=v=(tau+1)/2` is the unique fixed point; for even `tau`
there is none. A fixed point is an honest tie rather than a freely orientable
edge.

Reflection is not automatically a source operation. Starting from an
originally source-realizable pair `(20)`, the reflected type is

```text
N'=u,                         n0'=s+v-1,                 (25)
```

and `(18)` shows that it is source-realizable exactly when

```text
u>=rho,                       u-v<=s-1.                  (26)
```

Equation `(24)` is always an ambient trace identity; it describes a source
trace only under `(26)`.

## 6. Sharp quotient and symmetry hostiles

At `ell=10`, `s=5`, the valid packets

```text
u^3 p^5:       (u,v)=(4,5),   F_0=-z^3(1-z)^4,
x^2u^2p^6:     (u,v)=(4,6),   F_0=-z^3(1-z)^5            (27)
```

have the same strict orientation `u<v` and different traces. Thus the strict
order bit/quotient on the boundary pair does not determine the `q=0` consumer.

The packets

```text
u^3 p^5:        (u,v)=(4,5),
x^10u^2p^10:    (u,v)=(8,10)                             (28)
```

lie on the same primitive rational ray but have different boundary
multiplicities and traces. Unlike THM-4367's metric exit, this consumer reads
the common scale erased by gcd reduction.

The valid packet `x^2 u p^6` has `(u,v)=(3,6)`, but its reflection has
`N'=3<rho=4`; this disproves source closure of the ambient involution. In
contrast, the reflection of the first packet in `(27)` is valid:

```text
u^3p^2y^2:      (u,v)=(5,4),        Addr=33,
u^3p^5:         (u,v)=(4,5),        Addr=32,             (29)
```

and `32+33=8^2+1`. Finally,

```text
u^2p^2y^2:      (u,v)=(4,4),        Addr=25              (30)
```

is a source-realizable fixed-point control.

These are consumer-kernel statements, not an LRC-to-JC or tournament-to-JC
transfer. An injective triangular address retains the exact ordered pair; an
order bit or reduced Stern--Brocot ray does not.

## 7. Audit and scope

The 6,790,755-check primary and 14,325,934-check import-free referee
independently rebuild the Pascal transforms, fixed-row unit minors,
annihilator clock, all formal-series cases, source feasibility and overlap,
address decoder/reflection, and every named hostile, including `ell=2`,
`N=1`, and `q=N-1,N`. Normal, optimized, isolated, hash-seeded, and frozen
LF streams agree.

The theorem is confined to THM-4364's exact diagonal hierarchy and source
monomial packets. It does not produce a new bracket obstruction, prove that
an arbitrary hypothetical counterexample enters this chart, extend a finite
source stratum to all rows, prove polynomial termination, give seam entry or
a Keller pair, or prove `JC(2)` or `DC(2)`.

Reproduce from the repository root:

```text
python3 -B 04-computation/jc2_diagonal_boundary_address_simplex_stream_rank_thm4368.py
python3 -B -O 04-computation/jc2_diagonal_boundary_address_simplex_stream_rank_thm4368.py
python3 -B 04-computation/jc2_diagonal_boundary_address_simplex_stream_rank_independent_referee_thm4368.py
python3 -B -O 04-computation/jc2_diagonal_boundary_address_simplex_stream_rank_independent_referee_thm4368.py
```
