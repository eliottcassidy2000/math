---
id: HYP-8905
title: "Binary symmetric Keller subcase and separate planar descent programs"
status: >
  MIXED / CORRECTED BY MISTAKE-237. The homogeneous binary symmetric
  nilpotent-Hessian subcase is proved and lands in THM-2063. The proposed
  NC2/GMC-to-JC bridge, rank/cycle boundary, and equivalence among VC(4),
  leading-form/Jelonek, and Lame/Newton descent are not proved. JC(2) remains
  open, including higher fiber degree and higher geometric degree strata.
source: codex-2026-07-21-HYP-8905-audit
related:
  - THM-1330
  - THM-1345
  - THM-1830
  - THM-2045
  - THM-2063
  - HYP-8879
  - MISTAKE-235
  - MISTAKE-236
  - MISTAKE-237
---

# HYP-8905 -- what survives and what remains separate

## Exact survivor: the binary homogeneous symmetric subcase

Let `P in C[x,y]` be homogeneous of degree `d>=2`, and suppose its `2x2`
Hessian is nilpotent. Nilpotence gives both

```text
trace Hess(P)=Delta P=0,       det Hess(P)=0.             (1)
```

With `z=x+iy` and `zbar=x-iy`, a homogeneous binary harmonic polynomial is

```text
P=A z^d+B zbar^d.                                      (2)
```

Direct differentiation gives

```text
det Hess(P)
  =-4 d^2(d-1)^2 A B (x^2+y^2)^(d-2).                  (3)
```

Thus `AB=0`: `P` is one-sided. For the symmetric Keller map

```text
F=(x,y)+grad P,                                         (4)
```

one output-pencil member is linear on every fiber in a linear source
direction. For example, when `P=A z^d`,

```text
F_2-iF_1=-iz.                                          (5)
```

The conjugate branch has the analogous identity. THM-2063 therefore gives a
triangular normal form and explicit polynomial inverse. This proves the stated
binary homogeneous symmetric stratum; it does not prove the planar Jacobian
conjecture.

## Why the global bridge fails

The symmetric-Jacobian reduction used for a general planar Keller map changes
the ambient problem: the relevant symmetric/Vanishing-Conjecture target is in
four variables. Solving the binary symmetric target does not solve that image.
Likewise, Laplacian iterates, complex Gaussian moments, and LRC Fourier fibers
are differently typed functionals until an explicit map preserves the claimed
predicate. THM-1830 supplies no chain

```text
NC2 -> GMC(2) -> JC(2).                                 (6)
```

The existence of higher-rank nilpotent matrices in dimension at least three
also does not classify Keller collisions. Unique/coincident cycles and
rank-one/rank-two behavior remain analogy or search scheduling, not a theorem.

## Three live programs, not three equivalent forms

The following may share a descent/termination motif but remain distinct:

1. the four-variable symmetric/Vanishing-Conjecture route;
2. planar leading forms, nonproper values, and inverse-Jelonek geometry;
3. Newton-polygon or continued-fraction descent, including Lame/Fibonacci
   worst-case behavior for the ordinary Euclidean algorithm.

Any future transfer must record source, target, map, preserved predicate,
dimension change, information loss, and a hostile example. The degree-three
covering stratum is the first possible geometric degree after known exclusions,
not the whole remaining problem; higher degrees remain live.

## Reusable mechanism

The valid shared mechanism is **fiber-degree descent**. THM-2063 shows that
fiber degree at most one in any source/output-pencil direction is terminal.
A productive JC(2) search should therefore track the full biprojective table

```text
(source-fiber direction, output-pencil direction) -> exact fiber degree,   (7)
```

not only the displayed components or total degree. For each proposed descent,
record the decreasing quantity and prove that coordinate changes cannot reset
it. This is a lawful connection to return-termination work elsewhere in the
repository; it is not an identification of the underlying objects.
